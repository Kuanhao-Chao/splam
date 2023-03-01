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

bool exonmatch(GVec<GSeg> &prevexons, GVec<GSeg> &exons) {
	if(prevexons.Count() != exons.Count()) return false;
	for(int i=0;i<exons.Count();i++) {
		if(prevexons[i].end!=exons[i].end || prevexons[i].start!=exons[i].start) return false;
	}
	return true;
}

GStr filterSpurJuncs(GStr outfname_junc_score) {
    robin_hdd_rm_hit rm_hit;
    robin_hdd_rm_algn rm_algn;

    // infname_scorebed
    if (COMMAND_MODE == CLEAN) {

        // GSamReader bam_reader_spliced(ofn .chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
        int num_samples=in_records.start();


        outfile_discard = new GSamWriter(outfname_discard, in_records.header(), GSamFile_BAM);




        BundleData* bundle = new BundleData();

        uint runoffdist=200;
        GHash<int> hashread;      //read_name:pos:hit_index => readlist index
        GStr lastref;
        bool more_alns=true;
        int prev_pos=0;
        int lastref_id=-1; //last seen gseq_id

        bool fr_strand=false;
        bool rf_strand=false;
        int currentstart=0, currentend=0;


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
                // GMessage("Read: %s; %s\n", brec->refName(), brec->name());
                if (brec->isUnmapped()) {
                    outfile_discard->write(brec);
                    continue;
                }
                if (brec->start<1 || brec->mapped_len<10) {
                    if (verbose) GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
                            brec->name(), brec->start, brec->mapped_len);
                    outfile_discard->write(brec);
                    continue;
                }


                /***********************************
                 * Setting the "chr" "strand" of the current alignment.
                ************************************/
                refseqName=brec->refName();
                if (refseqName==NULL) {
                    GMessage("Error: cannot retrieve target seq name from BAM record!\n");
                    outfile_discard->write(brec);
                }

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
                chr_changed=(lastref.is_empty() || lastref!=refseqName);
                if (chr_changed) {
                    prev_pos=0;
                }
                if (pos<prev_pos) {
                    GMessage("[ERROR] %s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
                    brec->name(), brec->start,  pos, refseqName, prev_pos);
                    exit(-1);
                }
                prev_pos=pos;
                nh=brec->tag_int("NH");
                if (nh==0) nh=1;
                hi=brec->tag_int("HI");
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
                if (bundle->readlist.Count()>0) {
                    // process reads in previous bundle
                    bundle->getReady(currentstart, currentend);
                    processBundle(bundle);
                } else { 
                    //no read alignments in this bundle?  
                    bundle->Clear();
                } //nothing to do with this bundle

                if (chr_changed) {
                    lastref=refseqName;
                    lastref_id=gseq_id;
                    currentend=0;
                }

                if (!more_alns) {
                    noMoreBundles();
                    break;
                }
                currentstart=pos;
                currentend=brec->end;
                bundle->refseq=lastref;
                bundle->start=currentstart;
                bundle->end=currentend;
            } //<---- new bundle started

            if (currentend<(int)brec->end) {
                //current read extends the bundle
                //this might not happen if a longer guide had already been added to the bundle
                currentend=brec->end;
            } //adjusted currentend and checked for overlapping reference transcripts

            
            
            
            // GReadAlnData alndata(brec, 0, nh, hi);
            CReadAln alndata(brec);

            // if (xstrand=='+') alndata.strand=1;
            // else if (xstrand=='-') alndata.strand=-1;
            //const char* bname=brec->name();
            //GMessage("%s\t%c\t%d\thi=%d\n",bname, xstrand, alndata.strand,hi);
            //countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
            //fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);

            processRead(currentstart, currentend, *bundle, hashread, alndata);


        } //for each read alignment

        //cleaning up
        // delete brec;





































        // int prev_tid=-1;
        // GStr prev_refname;
        // int b_end=0, b_start=0;
        // while ((irec=in_records.next())!=NULL) {
        //     brec=irec->brec;
        //     int endpos=brec->end;
        //     GMessage("brec: %s\n", brec->refName());
        //     if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
        //         flushJuncs(joutf);
        //         b_start=brec->start;
        //         b_end=endpos;
        //         prev_tid=brec->refId();
        //         prev_refname=(char*)brec->refName();
        //     } else { //extending current bundle
        //         if (b_end<endpos) {
        //             b_end=endpos;
        //         }
        //     }
        //     int accYC = 0;
        //     accYC = brec->tag_int("YC", 1);
        //     if (brec->exons.Count()>1) {
        //         // Spliced reads



        //         char strand = brec->spliceStrand();
        //         std::string jnew_sub = std::to_string(brec->exons[e-1].end) + "_" + std::to_string(brec->exons[e].start-1) + "_" + strand + "_" + brec->refName();
        //         // std::string jnew_sub = std::to_string(brec->exons[e-1].end) + "_" + std::to_string(brec->exons[e].start-1) + "_" + brec->refName();
        //         if (rm_juncs.find(jnew_sub) != rm_juncs.end()) {
        //             spur = true;
        //             break;
        //         }



        //         outfile_spliced->write(brec);
        //         ALN_COUNT_SPLICED++;


        //     } else {
        //         // Non-spliced reads.
        //         // Not spliced => check their NH tags!
        //         if (brec->isUnmapped()) continue;


        //         int new_nh = brec->tag_int("NH", 0);
        //         if (new_nh == 1) {

        //         } else if (new_nh == 0){
        //             GMessage("\t\t brec->name(): %s !!!\n", brec->name());
        //             GMessage("\t\t NH tag is zero !!!: %d\n", new_nh);
        //         } else {

        //         }

        //     }
        //     ALN_COUNT++;
        //     if (ALN_COUNT % 1000000 == 0) {
        //         GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
        //     }
        // }









    } else if (COMMAND_MODE == ALL) {
        infname_scorebed;

        GSamReader bam_reader_spliced(outfname_spliced.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

        std::unordered_set<std::string> rm_juncs;
        // GMessage("Before rm_juncs.size()  %d\n", rm_juncs.size());
        loadBed(outfname_junc_score, rm_juncs);
        // GMessage("After rm_juncs.size()  %d\n", rm_juncs.size());

        int counter = 0, prev_tid=-1;
        GStr prev_refname;
        std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
        int b_end=0, b_start=0;
        progressbar bar(ALN_COUNT_SPLICED);
        bar.set_opening_bracket_char("[INFO] SPLAM! Removing junctions with low scores \n\t[");

        while ((brec=bam_reader_spliced.next())!=NULL) {
            bar.update();

            uint32_t dupcount=0;
            int endpos=brec->end;
            bool spur = false;
            if (brec->exons.Count() > 1) {
                for (int e=1; e<2; e++) {
                    char strand = brec->spliceStrand();
                    std::string jnew_sub = std::to_string(brec->exons[e-1].end) + "_" + std::to_string(brec->exons[e].start-1) + "_" + strand + "_" + brec->refName();
                    // std::string jnew_sub = std::to_string(brec->exons[e-1].end) + "_" + std::to_string(brec->exons[e].start-1) + "_" + brec->refName();
                    if (rm_juncs.find(jnew_sub) != rm_juncs.end()) {
                        spur = true;
                        break;
                    }
                }
            }
            if (spur) {
                std::string kv = brec->name();
                kv = kv + "_" + std::to_string(brec->pairOrder());

                GMessage("brec->mate_refId(): %s\n", brec->mate_refName());
                GMessage("brec->mate_start(): %d\n", brec->mate_start());
                GMessage("brec->mate_start(): %d\n", brec);

                GMessage("kv: %s\n", kv.c_str());

                if (rm_hit.find(kv) == rm_hit.end()) {
                    rm_hit[kv] = 1;
                } else {
                    rm_hit[kv]++;
                }
                // char* seq = brec->sequence();
                // char* cigar_seq = brec->cigar();
                // kv = kv + "_" + seq + "_" + cigar_seq + "_" + std::to_string(brec->flags()) + "_" + std::to_string(brec->start);
                // rm_rd_set.insert(kv);
                outfile_discard->write(brec);
                ALN_COUNT_BAD++;
                // free(seq);
                // free(cigar_seq);
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
        GMessage("[INFO] %d spurious alignments were removed.\n", ALN_COUNT_BAD);
        // Writing out auxiliary file 
        std::ofstream NH_tag_f;
        GStr outfname_NH_tag = out_dir + "/NH_tag_fix.csv";
        NH_tag_f.open(outfname_NH_tag.chars());   
        for (auto ele : rm_hit) {
            NH_tag_f << ele.first << ","  << ele.second  << std::endl;
        }
        NH_tag_f.close();
        return outfname_NH_tag;
    }
    return GStr("");
}

void processBundle(BundleData* bundle) {
    for (int i=0; i<bundle->readlist.Count(); i++) {
        int pair_idx = bundle->readlist[i]->pair_idx;
        GMessage("%s_%d_%c\t\t: %d; n: %d;    %d; np: %d\n", bundle->readlist[i]->brec.name(), bundle->readlist[i]->brec.pairOrder(), bundle->readlist[i]->brec.spliceStrand(), bundle->readlist[i]->brec.start, i, bundle->readlist[i]->brec.mate_start(), pair_idx);
    }
}


void processRead(int currentstart, int currentend, BundleData& bdata, GHash<int>& hashread, CReadAln& alndata) { // some false positives should be eliminated here in order to break the bundle
	GSamRecord& brec=(alndata.brec);			   // bam record
	GList<CReadAln>& readlist = bdata.readlist;    // list of reads gathered so far
	int readstart=brec.start;

    static GStr _id("", 256); //to prevent repeated reallocation for each parsed read
    static GStr _id_p("", 256); //to prevent repeated reallocation for each parsed read
	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Process read %s with exons:", brec.name());
		for (int i=0;i<brec.exons.Count();i++) {
			fprintf(stderr," %d-%d", brec.exons[i].start, brec.exons[i].end);
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

    CReadAln* readaln = new CReadAln(&brec);
    
    // for (int i=0;i<brec.exons.Count();i++) {
    //     readaln->len+=brec.exons[i].len();
    //     if(i) {
    //         int jstrand=strand;
    //         uint jstart=brec.exons[i-1].end;
    //         uint jend=brec.exons[i].start;
    //     }
    //     readaln->segs.Add(brec.exons[i]);
    // }
    n=readlist.Add(readaln);  // reset n for the case there is no match

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Add read %s with strand=%d and exons:",brec.name(),strand);
		for (int i=0;i<brec.exons.Count();i++) {
			fprintf(stderr," %d-%d", brec.exons[i].start, brec.exons[i].end);
		}
		fprintf(stderr,"\n");
		//fprintf(stderr,"Read %s is at n=%d with unitig_cov=%f and strand=%d\n",brec.name(),n,unitig_cov,strand);
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
			//GStr id(brec.name(), 16); // init id with readname
			_id.assign(brec.name()); //assign can be forced to prevent shrinking of the string
            _id_p.assign(brec.name());
            _id+='-';_id+=selfstart;_id+='-';_id+=pairstart;
			_id_p+='-';_id_p+=pairstart;_id_p+='-';_id_p+=selfstart;

            GMessage("_id.chars()  : %s\n", _id.chars());
            GMessage("_id_p.chars(): %s\n", _id_p.chars());
            // GMessage("n: %d\n", n);
			if(pairstart<=readstart) { // if I've seen the pair already <- I might not have seen it yet because the pair starts at the same place
				const int* np=hashread[_id_p.chars()];
                if (np) {
                    GMessage("np: %d\n\n", *np);

                    readlist[*np]->pair_idx = n;
                    readlist[n]->pair_idx = *np;
                }
                hashread.Remove(_id_p.chars());
				// if(np) { // the pair was stored --> why wouldn't it be? : only in the case that the pair starts at the same position
				// // 	// if(readlist[*np]->nh>nh && !nomulti) rdcount=float(1)/readlist[*np]->nh;
                // //     readlist[*np]->pair_idx



				// 	bool notfound=true;

                //     if(readlist[*np]->pair_idx[i]==n) {
                //     }


				// 	for(int i=0;i<readlist[*np]->pair_idx.Count();i++)
				// 		if(readlist[*np]->pair_idx[i]==n) {
				// 			readlist[*np]->pair_count[i]+=rdcount;
				// 			notfound=false;
				// 			break;
				// 		}
				// 	if(notfound) { // I didn't see the pairing before
				// 		readlist[*np]->pair_idx.Add(n);
				// 		readlist[*np]->pair_count.Add(rdcount);
				// 	}

				// // 	notfound=true;
				// // 	for(int i=0;i<readlist[n]->pair_idx.Count();i++)
				// // 		if(readlist[n]->pair_idx[i]==*np) {
				// // 			readlist[n]->pair_count[i]+=rdcount;
				// // 			notfound=false;
				// // 			break;
				// // 		}
				// // 	if(notfound) { // I didn't see the pairing before
				// // 		int i=*np;
				// // 		readlist[n]->pair_idx.Add(i);
				// // 		readlist[n]->pair_count.Add(rdcount);
				// // 	}
				// // 	hashread.Remove(_id.chars());
				// }
            }
			else { // I might still see the pair in the future
                GMessage("Bigger!\n");
				hashread.Add(_id.chars(), n);
			}
		}
	} //<-- if mate is mapped on the same chromosome
}


void noMoreBundles() {

}