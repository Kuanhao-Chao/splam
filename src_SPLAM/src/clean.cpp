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


void loadBed(GStr inbedname, std::unordered_set<std::string> &rm_juncs) {
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


    std::unordered_set<std::string> rm_juncs;
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
                // Unmapped reads.
                if (brec->isUnmapped()) {
                    // GMessage("Unmapped read: %s - %d\n", brec->name(), brec->pairOrder());
                    // removeAlignment(brec, rm_hit);
                    // continue;
                }
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
                int mate_end = brec->mate_start() + (int)brec->len();
                if (mate_end > (int)brec->end) {
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

void processBundle(BundleData* bundle, GList<CReadAln>& readlist, std::unordered_set<std::string>& rm_juncs, robin_hdd_rm_hit& rm_hit, int& bundle_counter) {
    bundle_counter += 1;
    GMessage("In bundle %d\n", bundle_counter);

    robin_hdd_int hash_good_pair_idx;
    robin_hdd_int hash_bad_pair_idx;

    // GMessage("\treadlist.Count(): %d\n", readlist.Count());
    // GMessage("\treadlist[0]: %p\n", &readlist[0]->brec);
    // GMessage("\treadlist[1]: %p\n", &readlist[1]->brec);
    // GMessage("\treadlist[2]: %p\n", &readlist[2]->brec);
    // GMessage("\treadlist[-1]: %p\n", &(readlist[readlist.Count()-1]->brec));

    // GMessage("\treadlist[0]->brec.start: %u\n", readlist[0]->brec.start);
    // GMessage("\treadlist[0]->brec.end  : %u\n", readlist[0]->brec.end);
    // GMessage("\treadlist[-1]->brec.start: %u\n", readlist[readlist.Count()-1]->brec.start);
    // GMessage("\treadlist[-1]->brec.end  : %u\n", readlist[readlist.Count()-1]->brec.end);

    for (int idx=0; idx<readlist.Count(); idx++) {
        int pair_idx = readlist[idx]->pair_idx;
        GMessage(">> pair_idx  : %d\n", pair_idx);

        // Just write out the read to cleaned BAM
        //    if it is unpaired.
        if (pair_idx == -1) {

        }

        GSamRecord& brec_bd = readlist[idx]->brec;
        GSamRecord& brec_bd_p = readlist[readlist[idx]->pair_idx]->brec;

        GMessage(">> brec_bd  : %s (%d)\n", brec_bd.name(), idx);
        GMessage(">> brec_bd_p: %s (%d)\n", brec_bd_p.name(), pair_idx);

        /***********************************
         * Writing out cleaned & discard alignment
         *  for the *second* seen paired alignments.
        ************************************/
        if (hash_good_pair_idx.find(idx) != hash_good_pair_idx.end()) {
            // It's already been added into the hash => 
            //  It is the second read in a pair being seen.
            keepAlignment(&brec_bd);
            continue;
        } 
        if (hash_bad_pair_idx.find(idx) != hash_bad_pair_idx.end()) {
            // It's already been added into the hash => 
            //  It is the second read in a pair being seen.
            continue;
        } 

        if (brec_bd.isUnmapped() || brec_bd_p.isUnmapped()) {
            // Writing out alignments into the file in a pair.
            // GMessage("Unmapped read: %s - %d\n", brec_bd.name(), brec_bd.pairOrder());

            GMessage("brec_bd  : %s\n", brec_bd.name());
            GMessage("brec_bd_p: %s\n", brec_bd_p.name());

            removeAlignment(&brec_bd, rm_hit);
            removeAlignment(&brec_bd_p, rm_hit);
            hash_bad_pair_idx.insert(idx);
            hash_bad_pair_idx.insert(pair_idx);
            // continue;
        }


        // Invalid mapping reads.
        if (brec_bd.refName()==NULL || brec_bd_p.refName()==NULL) {
            GMessage("Error: cannot retrieve target seq name from BAM record!\n");

            // GMessage("brec_bd  : %s\n", brec_bd.name());
            // GMessage("brec_bd_p: %s\n", brec_bd_p.name());

            removeAlignment(&brec_bd, rm_hit);
            removeAlignment(&brec_bd_p, rm_hit);
            hash_bad_pair_idx.insert(idx);
            hash_bad_pair_idx.insert(pair_idx);
        }

        /***********************************
         * Checking spurious spliced alignments for 
         *  frist & second alignments.
        ************************************/
        bool spur_m = false;
        bool spur_p = false;
        if (brec_bd.exons.Count() > 1) {
            for (int e=1; e<2; e++) {
                char strand = brec_bd.spliceStrand();
                std::string jnew_sub = std::to_string(brec_bd.exons[e-1].end) + "_" + std::to_string(brec_bd.exons[e].start-1) + "_" + strand + "_" + brec_bd.refName();
                // GMessage("jnew_sub: %s\n", jnew_sub.c_str());
                if (rm_juncs.find(jnew_sub) != rm_juncs.end()) {
                    spur_m = true;
                    break;
                }
            }
        }

        if (!spur_m) {
            // There are no spurious junctions inside the first alignment in a pair. => check it's mate.
            // GMessage("brec_bd_p.exons.Count(): %d\n", brec_bd_p.exons.Count());
            if (brec_bd_p.exons.Count() > 1) {
                for (int e=1; e<2; e++) {
                    char strand = brec_bd_p.spliceStrand();
                    // GMessage("strand: %c\n", strand);

                    // GMessage("brec_bd_p.exons[e-1].end: %u\n", brec_bd_p.exons[e-1].end);
                    // GMessage("brec_bd_p.exons[e].start-1: %u\n", brec_bd_p.exons[e].start-1);
                    // GMessage("brec_bd_p.refName(): %s\n", brec_bd_p.refName());

                    std::string jnew_sub_p = std::to_string(brec_bd_p.exons[e-1].end) + "_" + std::to_string(brec_bd_p.exons[e].start-1) + "_" + strand + "_" + brec_bd_p.refName();
                    // GMessage("jnew_sub_p: %s\n", jnew_sub_p.c_str());
                    if (rm_juncs.find(jnew_sub_p) != rm_juncs.end()) {
                        spur_p = true;
                        break;
                    }
                }
            }
        }

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
            removeAlignment(&brec_bd_p, rm_hit);
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
	int readstart=brec.start;

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

    // GMessage("refName: %s\n", bdata.readlist[n]->brec.refName());

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Add read %s with  and exons:",brec->name());
		for (int i=0;i<brec->exons.Count();i++) {
			fprintf(stderr," %d-%d", brec->exons[i].start, brec->exons[i].end);
		}
		fprintf(stderr,"\n");
		//fprintf(stderr,"Read %s is at n=%d with unitig_cov=%f and strand=%d\n",brec->name(),n,unitig_cov,strand);
	}
	*/

	if((int)brec.end>currentend) {
		currentend=brec.end;
	  	bdata.end=currentend;
	}

	// now set up the pairing
	if (brec.refId()==brec.mate_refId()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired
        int selfstart = brec.start;
		int pairstart = brec.mate_start();
		if (currentstart<=pairstart) { // if pairstart is in a previous bundle I don't care about it
			//GStr readname();
			//GStr id(brec->name(), 16); // init id with readname
			_id.assign(brec.name()); //assign can be forced to prevent shrinking of the string
            _id_p.assign(brec.name());
            _id+='-';_id+=selfstart;_id+='-';_id+=pairstart;
			_id_p+='-';_id_p+=pairstart;_id_p+='-';_id_p+=selfstart;

            // GMessage("_id.chars()  : %s\n", _id.chars());
            // GMessage("_id_p.chars(): %s\n", _id_p.chars());
            // GMessage("n: %d\n", n);
			if(pairstart<=readstart) { // if I've seen the pair already <- I might not have seen it yet because the pair starts at the same place
				const int* np=hashread[_id_p.chars()];
                if (np) {
                    // GMessage("np: %d\n\n", *np);

                    readlist[*np]->pair_idx = n;
                    readlist[n]->pair_idx = *np;
                }
                hashread.Remove(_id_p.chars());
            }
			else { // I might still see the pair in the future
                // GMessage("Bigger!\n");
				hashread.Add(_id.chars(), n);
			}
		}
	} //<-- if mate is mapped on the same chromosome
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