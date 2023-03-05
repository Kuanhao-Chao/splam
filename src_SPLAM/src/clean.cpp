#include "clean.h"
#include "common.h"
#include "junc.h"
#include "junc_func.h"
#include "util.h"
#include "extract.h"
#include "predict.h"
#include "bundle.h"

// #include "progressbar.h"
#include <progressbar/progressbar.hpp>

#include <unordered_set>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <htslib/htslib/faidx.h>
#include <Python.h>
#include <gclib/GStr.h>

/****************************
* Input : (1) Bed file of junction scores. (2) spliced alignment BAM file.
* Output: (2) unordered_set of reads, (2) hashmap of removed hits.
*****************************/
GStr splamClean() {
    /*********************************************
     * Step 4: SPLAM filtering out reads.
    *********************************************/
    STEP_COUNTER += 1;
    GMessage("\n###########################################\n");
    GMessage("## Step %d: SPLAM filtering out reads\n", STEP_COUNTER);
    GMessage("###########################################\n\n");
    GStr outfname_NH_tag = filterSpurJuncs(infname_scorebed);

    delete outfile_discard;
    delete outfile_cleaned_tmp;

    return outfname_NH_tag;
}


void loadBed(GStr inbedname, robin_hdd_string &rm_juncs) {
    std::ifstream bed_f(inbedname);
    std::string line;
    int bed_counter = 0;
    while (getline(bed_f, line)) {
        bed_counter ++;
        GStr gline = line.c_str();
        GVec<GStr> junc;
        int cnt = 0;
        while (cnt < 7) {
            GStr tmp = gline.split("\t");
            junc.Add(gline);
            gline=tmp;
            cnt++;
        }
        if (junc[6].asDouble() <= threshold) {
            char* chrname =junc[0].detach();
            char* strand =junc[5].detach();
            std::string j = std::to_string(junc[1].asInt()) + "_" + std::to_string(junc[2].asInt()) + "_" + strand + "_" + chrname;
            rm_juncs.insert(j);
        }
    }
}

GStr filterSpurJuncs(GStr outfname_junc_score) {
    robin_hdd_rm_hit rm_hit;
    robin_hdd_rm_algn rm_algn;


    robin_hdd_string rm_juncs;
    loadBed(outfname_junc_score, rm_juncs);
    // infname_scorebed
    if (COMMAND_MODE == CLEAN) {
        int num_samples=in_records.start();
        outfile_discard = new GSamWriter(outfname_discard, in_records.header(), GSamFile_BAM);
        outfile_cleaned_tmp = new GSamWriter(outfname_cleaned_tmp, in_records.header(), GSamFile_BAM);

        BundleData* bundle = new BundleData();
	    GList<CReadAln> readlist;

        uint runoffdist=100;
        GHash<int> hashread; //read_name:pos:hit_index => readlist index
        GStr lastref;
        bool more_alns=true;
        int prev_pos=0;
        int lastref_id=-1; //last seen gseq_id

        bool fr_strand=false;
        bool rf_strand=false;
        int currentstart=0, currentend=0;
        int bundle_counter = 0;

        int max_splice_distance = 20000;

        while (more_alns) {
            bool chr_changed=false;
            int pos=0;
            const char* refseqName=NULL;
            char xstrand=0;
            int nh=1;
            int hi=0;
            int gseq_id=lastref_id;  //current chr id
            bool new_bundle=false;
            //delete brec;
            if ((irec=in_records.next())!=NULL) {
                brec=irec->brec;
                /***********************************
                 * Alignment filtering criteria.
                ************************************/
                // // Unmapped reads.
                // if (brec->isUnmapped()) {
                //     // GMessage("Unmapped read: %s - %d\n", brec->name(), brec->pairOrder());
                //     // removeAlignment(brec, rm_hit);
                //     // continue;
                // }
                // // Invalid mapping reads.
                // if (brec->start<1 || brec->mapped_len<10) {
                //     GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
                //             brec->name(), brec->start, brec->mapped_len);
                //     removeAlignment(brec, rm_hit);
                //     continue;
                // }


                /***********************************
                 * Setting the "chr" "strand" of the current alignment.
                ************************************/
                refseqName=brec->refName();
                // // Invalid mapping reads.
                // if (refseqName==NULL) {
                //     GMessage("Error: cannot retrieve target seq name from BAM record!\n");
                //     removeAlignment(brec, rm_hit);
                // }

                xstrand=brec->spliceStrand(); // tagged strand gets priority
                if(xstrand=='.' && (fr_strand || rf_strand)) { // set strand if stranded library
                    if(brec->isPaired()) { // read is paired
                        if(brec->pairOrder()==1) { // first read in pair
                            if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
                            else xstrand='-';
                        }
                        else {
                            if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='-';
                            else xstrand='+';
                        }
                    }
                    else {
                        if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
                        else xstrand='-';
                    }
                }

                /***********************************
                 * Setting the "chr_changed" and "new_bundle" parameters.
                ************************************/
                pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based

                // GMessage("pos: %d; prev_pos: %d\n", pos, prev_pos);

                chr_changed=(lastref.is_empty() || lastref!=refseqName);
                if (chr_changed) {
                    prev_pos=0;
                }

                if (pos == 0) {
                    // This is an unmapped read
                } else if (pos<prev_pos) {
                    GMessage("[ERROR] %s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
                    brec->name(), brec->start,  pos, refseqName, prev_pos);
                    exit(-1);
                }
                prev_pos=pos;
                nh=brec->tag_int("NH", 0);
                if (nh==0) nh=1;
                hi=brec->tag_int("HI", 0);
                if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
                    new_bundle=true;
                }
            } else { //no more alignments
                more_alns=false;
                new_bundle=true; //fake a new start (end of last bundle)
            }

            /***********************************
             * Process the bundle!
            ************************************/
            if (new_bundle || chr_changed) {
                hashread.Clear();
                if (readlist.Count()>0) {
                    // process reads in previous bundle
                    bundle->getReady(currentstart, currentend);

                    GMessage(">> bundle read count: %d\n", readlist.Count());
                    GMessage(">> bundle start     : %d\n", bundle->start);
                    GMessage(">> bundle end       : %d\n", bundle->end);

                    processBundle(bundle, readlist, rm_juncs, rm_hit, bundle_counter);
                    readlist.Clear();
                } else { 
                    //no read alignments in this bundle?  
                    bundle->Clear();
                    readlist.Clear();
                } //nothing to do with this bundle

                if (chr_changed) {
                    lastref = refseqName;
                    lastref_id = gseq_id;
                    currentend = 0;
                }

                if (!more_alns) {
                    noMoreBundles();
                    break;
                }

                if (brec->start > 0) {
                    currentstart = pos;
                    currentend = brec->end;
                }
                bundle->refseq = lastref;
                bundle->start = currentstart;
                bundle->end = currentend;
            } //<---- new bundle started

            int fragment_end = 0;
            if (brec->refId() == brec->mate_refId()) {

                int insert_size = brec->insertSize();
                int mate_end = brec->mate_start() + insert_size;
                
                if (mate_end > (int)brec->end && insert_size <= max_splice_distance) {
                    fragment_end = mate_end;
                } else {
                    fragment_end = (int)brec->end;
                }
            } else {
                fragment_end = (int)brec->end;
            }

            if (currentend<fragment_end) {
                //current read extends the bundle
                currentend=fragment_end;
            } //adjusted currentend and checked for overlapping reference transcripts

            // GMessage("brec->refName(): %s\n", brec->refName());
            CReadAln* alndata = new CReadAln(brec);
            processRead(currentstart, currentend, readlist, *bundle, hashread, alndata);
        } //for each read alignment




















    } else if (COMMAND_MODE == ALL) {
        GSamReader bam_reader_spliced(outfname_spliced.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

        int counter = 0, prev_tid=-1;
        GStr prev_refname;
        std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
        int b_end=0, b_start=0;
        progressbar bar(ALN_COUNT_SPLICED);
        bar.set_opening_bracket_char("[INFO] SPLAM! Removing junctions with low scores \n\t[");

        while ((brec = bam_reader_spliced.next())!=NULL) {
            bar.update();
            uint32_t dupcount=0;
            int endpos=brec->end;
            bool spur = false;
            if (brec->exons.Count() > 1) {
                for (int e=1; e<2; e++) {
                    char strand = brec->spliceStrand();
                    std::string jnew_sub = std::to_string(brec->exons[e-1].end) + "_" + std::to_string(brec->exons[e].start-1) + "_" + strand + "_" + brec->refName();
                    if (rm_juncs.find(jnew_sub) != rm_juncs.end()) {
                        spur = true;
                        break;
                    }
                }
            }
            if (spur) {
                removeAlignment(brec, rm_hit);
            } else {
                int nh_tag = brec->tag_int("NH", 1);
                if (nh_tag == 1) {
                    outfile_cleaned->write(brec);
                } else {
                    outfile_multimapped->write(brec);
                    ALN_COUNT_NH_UPDATE++;
                }
                ALN_COUNT_GOOD++;
            }
        }
        bam_reader_spliced.bclose();
        delete outfile_multimapped;

        GMessage("\n");
        ALN_COUNT_GOOD += ALN_COUNT_NSPLICED;
    }
    GMessage("[INFO] %d spurious alignments were removed.\n", ALN_COUNT_BAD);
    GStr outfname_NH_tag = writenhHitFile(rm_hit);
    return outfname_NH_tag;
}

void processBundle(BundleData* bundle, GList<CReadAln>& readlist, robin_hdd_string& rm_juncs, robin_hdd_rm_hit& rm_hit, int& bundle_counter) {
    bundle_counter += 1;
    // GMessage("In bundle %d\n", bundle_counter);

    robin_hdd_int hash_good_pair_idx;
    robin_hdd_int hash_bad_pair_idx;

    for (int idx=0; idx<readlist.Count(); idx++) {
        bool remove_algn = false;
        int pair_idx = readlist[idx]->pair_idx;
        GSamRecord& brec_bd = readlist[idx]->brec;


        // GMessage("idx      : %d\n", idx);
        // GMessage("pair_idx : %d\n", pair_idx);
        

        /***********************************
         * Case 1: cannot find its pair.
        ************************************/
        // Check the global hash => only used when its mate is unpaired.
        if (pair_idx == -1) {
            // Check if the read is spuriously spliced.
            bool spur = false;
            spur = alignmentAssessment(&brec_bd, rm_juncs);
            if (spur) {
                removeAlignment(&brec_bd, rm_hit);
                hash_bad_pair_idx.insert(idx);
            } else {
                // make the read single-ended.
                // GMessage("brec_bd      : %d\n", brec_bd.flags());

                if (brec_bd.isUnmapped()) {
                    // Make it single-end but its unmapped => removed
                    removeAlignment(&brec_bd, rm_hit);
                    hash_bad_pair_idx.insert(idx);
                } else {
                    int sngle_end_filter = 3860;
                    brec_bd.set_flags(brec_bd.flags()&sngle_end_filter);
                    // GMessage("brec_bd new  : %d\n", brec_bd.flags());
                    // brec_bd
                    brec_bd.unpair_mate_refName();
                    brec_bd.unpair_mate_start();
                    keepAlignment(&brec_bd);
                    hash_good_pair_idx.insert(idx);
                }
            }
            continue;
        }
        /***********************************
         * Case 2: Find its pair.
        ************************************/
        // Process both the current and paired alginments.
        GSamRecord& brec_bd_p = readlist[readlist[idx]->pair_idx]->brec;

        /***********************************
         * Writing out cleaned & discard alignment
         *  for the *second* seen paired alignments.
        ************************************/
        if (hash_good_pair_idx.find(idx) != hash_good_pair_idx.end()) {
            // It's already been added into the hash => 
            //  It is the second read in a pair being seen.
            // GMessage("Alingment in 'hash_good_pair_idx'\n");
            keepAlignment(&brec_bd);
            continue;
        } 
        if (hash_bad_pair_idx.find(idx) != hash_bad_pair_idx.end()) {
            // It's already been added into the hash => 
            //  It is the second read in a pair being seen.
            // GMessage("Alingment in 'hash_bad_pair_idx'\n");
            removeAlignment(&brec_bd, rm_hit);
            continue;
        } 

        // Invalid mapping reads. => I can skip this 
        if (brec_bd.refName()==NULL || brec_bd_p.refName()==NULL) {
            GMessage("Error: cannot retrieve target seq name from BAM record!\n");

            // GMessage("brec_bd  : %s\n", brec_bd.name());
            // GMessage("brec_bd_p: %s\n", brec_bd_p.name());
            GMessage("Refname null'\n");

            removeAlignment(&brec_bd, rm_hit);
            hash_bad_pair_idx.insert(idx);
            hash_bad_pair_idx.insert(pair_idx);
        }

        /***********************************
         * Checking spurious spliced alignments for 
         *  frist & second alignments.
        ************************************/
        bool spur_m = false;
        bool spur_p = false;
        spur_m = alignmentAssessment(&brec_bd, rm_juncs);
        if (!spur_m) {
            spur_p = alignmentAssessment(&brec_bd_p, rm_juncs);
        }

        // GMessage("Spurious results: %d; %d'\n", spur_m, spur_p);

        /***********************************
         * Writing out cleaned & discard alignment
         *  for the *first* seen paired alignments.
        ************************************/
        if (!spur_m && !spur_p) {
            // Writing out two reads.
            keepAlignment(&brec_bd);
            hash_good_pair_idx.insert(idx);
            hash_good_pair_idx.insert(pair_idx);
        } else {
            // Writing out two alignments into a discarded BAM file.

            // GMessage("brec_bd  : %s\n", brec_bd.name());
            // GMessage("brec_bd_p: %s\n", brec_bd_p.name());

            removeAlignment(&brec_bd, rm_hit);
            hash_bad_pair_idx.insert(idx);
            hash_bad_pair_idx.insert(pair_idx);
        }

        // GMessage("%s_%d_%c\t\t: %d; n: %d;    %d; np: %d\n", readlist[idx]->brec.name(), readlist[idx]->brec.pairOrder(), readlist[idx]->brec.spliceStrand(), readlist[idx]->brec.start, idx, readlist[idx]->brec.mate_start(), pair_idx);
    }
    bundle->Clear();
}

void processRead(int currentstart, int currentend, GList<CReadAln>& readlist, BundleData& bdata, GHash<int>& hashread, CReadAln* alndata) { // some false positives should be eliminated here in order to break the bundle

	GSamRecord& brec=(alndata->brec);			   // bam record
	// GList<CReadAln>& readlist = bdata.readlist;    // list of reads gathered so far
    static GStr _id("", 256); //to prevent repeated reallocation for each parsed read
    static GStr _id_p("", 256); //to prevent repeated reallocation for each parsed read
	
    /*
	{ // DEBUG ONLY
		fprintf(stderr,"Process read %s with exons:", brec->name());
		for (int i=0;i<brec->exons.Count();i++) {
			fprintf(stderr," %d-%d", brec->exons[i].start, brec->exons[i].end);
		}
		fprintf(stderr,"\n");
	}
	*/

	double nm=(double)brec.tag_int("NM"); // read mismatch
	float unitig_cov = unitig_cov=brec.tag_float("YK");

	bool match=false;  // true if current read matches a previous read
	int n = 0;
    // readlist.Count()-1;
	if (bdata.end<currentend) {
		bdata.start=currentstart;
		bdata.end=currentend;
	}
	bdata.numreads++;                         // number of reads gets increased no matter what

    // GMessage("brec->refName(): %s\n", brec->refName());
    // CReadAln* readaln = new CReadAln(&brec);
    
    // // for (int i=0;i<brec->exons.Count();i++) {
    // //     readaln->len+=brec->exons[i].len();
    // //     if(i) {
    // //         int jstrand=strand;
    // //         uint jstart=brec->exons[i-1].end;
    // //         uint jend=brec->exons[i].start;
    // //     }
    // //     readaln->segs.Add(brec->exons[i]);
    // // }
    n=readlist.Add(alndata); // reset n for the case there is no match


	if((int)brec.end>currentend) {
		currentend=brec.end;
	  	bdata.end=currentend;
	}


	// now set up the pairing
	if (brec.refId()==brec.mate_refId()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired

        int self_start = brec.start;
		int pair_start = brec.mate_start();
        int insert_size = brec.insertSize();
        int pair_insert_size = (-1)*insert_size;
        int pair_idx = brec.pairOrder();
        if (brec.isUnmapped()) {
            self_start = pair_start;
        } else if (brec.isMateUnmapped()) {
            pair_start = self_start;
        }



        // GMessage("brecname: %s; self_start: %d;  pair_start: %d;  currentstart: %d; insert_size: %d; pair_insert_size:%d; pair_idx:%d\n", brec.name(), self_start, pair_start, currentstart, insert_size, pair_insert_size, pair_idx);
		if (currentstart<=pair_start) { // if pair_start is in a previous bundle I don't care about it
			//GStr readname();
			//GStr id(brec->name(), 16); // init id with readname
			_id.assign(brec.name()); //assign can be forced to prevent shrinking of the string
            _id_p.assign(brec.name());


            _id+=';';_id+=self_start;
            _id+=';';_id+=pair_start;
            _id+=';';_id+=insert_size;
            _id+=';';_id+=pair_idx;

			_id_p+=';';_id_p+=pair_start;
            _id_p+=';';_id_p+=self_start;
            _id_p+=';';_id_p+=pair_insert_size;
            _id_p+=';';_id_p+=(3-pair_idx);




            int* n_check=hashread[_id.chars()];
            while (n_check) {
                // GMessage("element in readlist \n");
                // GMessage("\tChecking repeat: %s;  %d\n", _id.chars(), *n_check);

                _id+='*';
                _id_p+='*';
                n_check=hashread[_id.chars()];
                // GMessage("\t old n: %d;\n", n);
                // readlist.Remove(alndata);
                // n = n-1;
                // GMessage("\t new n: %d \n", n);
            }




			if(pair_start < self_start) { // if I've seen the pair already <- I might not have seen it yet because the pair starts at the same place
				const int* np=hashread[_id_p.chars()];
                if (np) {
                    // GMessage("\t\tn : %d\n\n", n);
                    // GMessage("\t\tnp: %d\n\n", *np);

                    readlist[*np]->pair_idx = n;
                    readlist[n]->pair_idx = *np;
                } else {
                    // GMessage(">> Pair not in the same bundle\n");
                }
                // hashread.Remove(_id_p.chars());



            } else if (pair_start == self_start) {
				hashread.Add(_id.chars(), n);
				const int* np=hashread[_id_p.chars()];
                // GMessage("\t_equal condition\n");
                // GMessage("\t\tn : %d\n\n", n);
                // GMessage("\t\tnp: %d\n\n", *np);                    

                if (np) {
                    // GMessage("\t\tn : %d\n\n", n);
                    // GMessage("\t\tnp: %d\n\n", *np);
                    readlist[*np]->pair_idx = n;
                    readlist[n]->pair_idx = *np;
                    // hashread.Remove(_id.chars());
                    // hashread.Remove(_id_p.chars());
                // } else {
                //     GMessage("\tHasn't seen the pair yet\n");
                }




            } else { // I might still see the pair in the future

			}


            // GMessage("Adding read to hash: %s;  %d\n", _id.chars(), n);
            hashread.Add(_id.chars(), n);

            const int* n_check_final=hashread[_id.chars()];

            // GMessage("Retrieving read from hash: %s;  %d\n", _id.chars(), *n_check_final);
		}
	} else {

    }
}


void noMoreBundles() {

}

void removeAlignment(GSamRecord* brec, robin_hdd_rm_hit& rm_hit) {
    outfile_discard->write(brec);

    std::string kv = brec->name();
    kv = kv + "_" + std::to_string(brec->pairOrder());
    if (rm_hit.find(kv) == rm_hit.end()) {
        rm_hit[kv] = 1;
    } else {
        rm_hit[kv]++;
    }
    ALN_COUNT_BAD++;
}

void keepAlignment(GSamRecord* brec) {
    // std::string key = get_global_removed_algns_key(brec);
    outfile_cleaned_tmp->write(brec);
    ALN_COUNT_GOOD++;
}

GStr writenhHitFile(robin_hdd_rm_hit& rm_hit) {
    // Writing out auxiliary file 
    std::ofstream NH_tag_f;
    GStr outfname_NH_tag = out_dir + "/TMP/NH_tag_fix.csv";
    NH_tag_f.open(outfname_NH_tag.chars());   
    for (auto ele : rm_hit) {
        NH_tag_f << ele.first << ","  << ele.second  << std::endl;
    }
    NH_tag_f.close();
    return outfname_NH_tag;
}

bool alignmentAssessment(GSamRecord* brec, robin_hdd_string &rm_juncs) {
    bool spur = false;
    if (brec->exons.Count() > 1) {
        for (int e=1; e<2; e++) {
            char strand = brec->spliceStrand();
            std::string jnew_sub = std::to_string(brec->exons[e-1].end) + "_" + std::to_string(brec->exons[e].start-1) + "_" + strand + "_" + brec->refName();
            // GMessage("jnew_sub: %s\n", jnew_sub.c_str());
            if (rm_juncs.find(jnew_sub) != rm_juncs.end()) {
                spur = true;
                return spur;
            }
        }
    }
    return spur;
}

// std::string get_global_removed_algns_key(GSamRecord* brec) {
//     std::string brec_name = brec->name();
//     int brec_start = brec->start;
//     int brec_mate_start = brec->mate_start();
//     int pair = brec->pairOrder();


//     // GMessage("~ brec->refName()     : %s\n", brec->refName());
//     // GMessage("~ brec->mate_refName(): %s\n", brec->mate_refName());
//     // GMessage("~ brec->mapped_len    : %d\n", brec->mapped_len);
//     // brec->mate_refName();
//     // brec->mapped_len;

//     std::string key = brec_name + "_" + std::to_string(brec_start) + "_" + std::to_string(brec_mate_start) + "_" + std::to_string(pair);
//     // bam1_t* brec->b;
//     // GMessage("key: %s\n", key.c_str());
//     return key;
// }

std::string get_global_removed_mate_algns_key(GSamRecord* brec) {
    std::string brec_name = brec->name();
    int brec_start = brec->start;
    int brec_mate_start = brec->mate_start();
    int pair = 3 - brec->pairOrder();
    std::string key = brec_name + "_" + std::to_string(brec_mate_start) + "_" + std::to_string(brec_start) + "_" + std::to_string(pair);

    // GMessage("mate key: %s\n", key.c_str());
    return key;
}