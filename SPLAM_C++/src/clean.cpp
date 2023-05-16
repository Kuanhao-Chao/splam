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

    if (verbose) {
        GMessage("\n###########################################\n");
        GMessage("## Step %d: SPLAM filtering out spurious alignments\n", STEP_COUNTER);
        GMessage("###########################################\n");
    }
    GStr outfname_NH_tag = filterSpurJuncs(infname_scorebed);
    return outfname_NH_tag;
}

void loadBed(GStr inbedname, robin_hdd_string &rm_juncs) {
    GMessage("threshold: %f\n", threshold);
    std::ifstream bed_f(inbedname);
    std::string line;
    while (getline(bed_f, line)) {
        GStr gline = line.c_str();
        GVec<GStr> junc;
        int cnt = 0;
        while (cnt < 8) {
            GStr tmp = gline.split("\t");
            junc.Add(gline);
            gline=tmp;
            cnt++;
        }
        double donor_score = junc[6].asDouble();
        double acceptor_score = junc[7].asDouble();
        if (donor_score <= threshold && acceptor_score <= threshold) {
            char* chrname =junc[0].detach();
            char* strand =junc[5].detach();
            std::string j = std::to_string(junc[1].asInt()) + "_" + std::to_string(junc[2].asInt()) + "_" + strand + "_" + chrname;
            rm_juncs.insert(j);
            JUNC_COUNT_BAD ++;
        } else {
            JUNC_COUNT_GOOD ++;
        }
    }
    // for (auto ele : rm_juncs) {
    //     GMessage("rm_juncs : %s\n", ele.c_str());
    //     // GMessage("ele.second : %d\n", ele.second);
    // }
}

GStr filterSpurJuncs(GStr outfname_junc_score) {
    robin_hdd_rm_hit rm_hit;
    robin_hdd_rm_algn rm_algn;


    robin_hdd_string rm_juncs;
    loadBed(outfname_junc_score, rm_juncs);
    // infname_scorebed
    if (g_paired_removal) {
        /*********************************************
         * Cleaning up alignments by pairs.
        *********************************************/
        
        /*********************************
         * Processonmg uniquely mapped paired spliced reads
        *********************************/
        GSamReader reader_s_uniq_map(outfname_s_uniq_map.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

        bool uniq_next_main_aln = true;
        GSamRecord* uniq_brec_1st;

        progressbar bar_uniq(ALN_COUNT_SPLICED_UNIQ);
        bar_uniq.set_opening_bracket_char("[INFO] SPLAM! Filtering unique spliced alignments (paired) \n\t[");
        while ( (brec = reader_s_uniq_map.next())!=NULL ) {
            if (verbose) {
                bar_uniq.update();
            }
            bool spur = alignmentAssessment(brec, rm_juncs);
            uniq_brec_1st = new GSamRecord(*brec);
            brec = reader_s_uniq_map.next();
            bool spur_pair = alignmentAssessment(brec, rm_juncs);
            if (!g_2_stage_run) {
                if (spur && spur_pair) {
                    // Make sure both reads are unmapped.
                    update_flag_paired_remove_both(uniq_brec_1st, brec);
                    removeAlignment(outfile_discard_s_uniq_map, uniq_brec_1st, rm_hit);
                    removeAlignment(outfile_discard_s_uniq_map, brec, rm_hit);
                    ALN_COUNT_SPLICED_UNIQ_DISCARD += 2;
                } else if (spur && !spur_pair) {
                    update_flag_paired_remove_one(uniq_brec_1st, brec);
                    removeAlignment(outfile_discard_s_uniq_map, uniq_brec_1st, rm_hit);
                    keepAlignment(outfile_cleaned, brec);
                    ALN_COUNT_SPLICED_UNIQ_DISCARD += 1;
                } else if (!spur && spur_pair) {
                    update_flag_paired_remove_one(brec, uniq_brec_1st);
                    keepAlignment(outfile_cleaned, uniq_brec_1st);
                    removeAlignment(outfile_discard_s_uniq_map, brec, rm_hit);
                    ALN_COUNT_SPLICED_UNIQ_DISCARD += 1;
                } else if (!spur && !spur_pair) {
                    keepAlignment(outfile_cleaned, uniq_brec_1st);
                    keepAlignment(outfile_cleaned, brec);
                }
            } else {
                if (spur && spur_pair) {
                    // Make sure both reads are unmapped.
                    update_flag_paired_remove_both(uniq_brec_1st, brec);
                    removeAlignment(outfile_discard_s_uniq_map, uniq_brec_1st, rm_hit);
                    removeAlignment(outfile_discard_s_uniq_map, brec, rm_hit);
                    ALN_COUNT_SPLICED_UNIQ_DISCARD += 2;
                } else if (spur && !spur_pair) {
                    update_flag_paired_remove_one(uniq_brec_1st, brec);
                    removeAlignment(outfile_discard_s_uniq_map, uniq_brec_1st, rm_hit);
                    keepAlignment(outfile_cleaned_2stage, brec);
                    ALN_COUNT_SPLICED_UNIQ_DISCARD += 1;
                } else if (!spur && spur_pair) {
                    update_flag_paired_remove_one(brec, uniq_brec_1st);
                    keepAlignment(outfile_cleaned_2stage, uniq_brec_1st);
                    removeAlignment(outfile_discard_s_uniq_map, brec, rm_hit);
                    ALN_COUNT_SPLICED_UNIQ_DISCARD += 1;
                } else if (!spur && !spur_pair) {
                    keepAlignment(outfile_cleaned_2stage, uniq_brec_1st);
                    keepAlignment(outfile_cleaned_2stage, brec);
                }
            }
            delete uniq_brec_1st;
        }
        reader_s_uniq_map.bclose();
        if (verbose) GMessage("\n");

        /*********************************
         * Processonmg multi-mapped paired spliced reads
        *********************************/
        GSamReader reader_s_multi_map(outfname_s_multi_map.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

        bool multi_next_main_aln = true;
        GSamRecord* multi_brec_1st;

        progressbar bar_multi(ALN_COUNT_SPLICED_MULTI);
        bar_multi.set_opening_bracket_char("[INFO] SPLAM! Filtering multi-mapped spliced alignments (paired) \n\t[");
        while ( (brec = reader_s_multi_map.next())!=NULL ) {
            if (verbose) {
                bar_multi.update();
            }

            bool spur = alignmentAssessment(brec, rm_juncs);
            multi_brec_1st = new GSamRecord(*brec);
            brec = reader_s_multi_map.next();
            bool spur_pair = alignmentAssessment(brec, rm_juncs);
            if (spur && spur_pair) {
                update_flag_paired_remove_both(multi_brec_1st, brec);
                removeAlignment(outfile_discard_s_multi_map, multi_brec_1st, rm_hit);
                removeAlignment(outfile_discard_s_multi_map, brec, rm_hit);
                ALN_COUNT_SPLICED_MULTI_DISCARD += 2;
            } else if (spur && !spur_pair) {
                update_flag_paired_remove_one(multi_brec_1st, brec);
                removeAlignment(outfile_discard_s_multi_map, multi_brec_1st, rm_hit);
                outfile_s_multi_map_tmp->write(brec);
                ALN_COUNT_SPLICED_MULTI_DISCARD += 1;
            } else if (!spur && spur_pair) {
                update_flag_paired_remove_one(brec, multi_brec_1st);
                outfile_s_multi_map_tmp->write(multi_brec_1st);
                removeAlignment(outfile_discard_s_multi_map, brec, rm_hit);
                ALN_COUNT_SPLICED_MULTI_DISCARD += 1;
            } else if (!spur && !spur_pair) {
                outfile_s_multi_map_tmp->write(multi_brec_1st);
                outfile_s_multi_map_tmp->write(brec);
            }
            delete multi_brec_1st;
        }
        reader_s_multi_map.bclose();
        if (verbose) GMessage("\n");


        /*********************************
         * Processing uniq-mapped unpaired alignments
        *********************************/
        GSamReader reader_s_uniq_unpair(outfname_s_uniq_unpair.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
        progressbar bar_uniq_unpair(ALN_COUNT_SPLICED_UNIQ_UNPAIR);
        GMessage("ALN_COUNT_SPLICED_UNIQ_UNPAIR: %d\n", ALN_COUNT_SPLICED_UNIQ_UNPAIR);

        bar_uniq_unpair.set_opening_bracket_char("[INFO] SPLAM! Filtering unique spliced alignments (unpaired) \n\t[");
        while ( (brec = reader_s_uniq_unpair.next())!=NULL ) {
            if (verbose) {
                bar_uniq_unpair.update();
            }
            bool spur = alignmentAssessment(brec, rm_juncs);
            if (spur) {
                update_flag_unpair_remove(brec);
                removeAlignment(outfile_discard_s_uniq_map, brec, rm_hit);
                ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD += 1;
            } else {
                update_flag_unpair_kept(brec);
                if (!g_2_stage_run) {
                    keepAlignment(outfile_cleaned, brec);
                } else {
                    keepAlignment(outfile_cleaned_2stage, brec);
                }
            }
        }
        reader_s_uniq_unpair.bclose();
        if (verbose) GMessage("\n");


        /*********************************
         * Processing multi-mapped unpaired alignments
        *********************************/
        GSamReader reader_s_multi_unpair(outfname_s_multi_unpair.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
        progressbar bar_multi_unpair(ALN_COUNT_SPLICED_MULTI_UNPAIR);

        bar_multi_unpair.set_opening_bracket_char("[INFO] SPLAM! Filtering multi-mapped spliced alignments (unpaired) \n\t[");
        while ( (brec = reader_s_multi_unpair.next())!=NULL ) {
            if (verbose) {
                bar_multi_unpair.update();
            }
            bool spur = alignmentAssessment(brec, rm_juncs);
            if (spur) {
                update_flag_unpair_remove(brec);
                removeAlignment(outfile_discard_s_multi_map, brec, rm_hit);
                ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD += 1;
            } else {
                outfile_s_multi_unpair_tmp->write(brec);
            }
        }
        reader_s_multi_unpair.bclose();
        if (verbose) GMessage("\n");
        delete outfile_s_multi_unpair_tmp;





    } else {

        /*********************************************
         * Cleaning up alignments by individuals (not-paired).
        *********************************************/

        /*********************************
         * Processing unique spliced alignments
        *********************************/
        GSamReader reader_s_uniq_map(outfname_s_uniq_map.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
        progressbar bar_uniq(ALN_COUNT_SPLICED_UNIQ);
        bar_uniq.set_opening_bracket_char("[INFO] SPLAM! Filtering unique spliced alignments \n\t[");
        while ( (brec = reader_s_uniq_map.next())!=NULL ) {
            if (verbose) {
                bar_uniq.update();
            }
            bool spur = alignmentAssessment(brec, rm_juncs);
            if (spur) {
                update_flag_unpair_remove(brec);
                removeAlignment(outfile_discard_s_uniq_map, brec, rm_hit);
                ALN_COUNT_SPLICED_UNIQ_DISCARD += 1;
            } else {
                if (!g_2_stage_run) {
                    keepAlignment(outfile_cleaned, brec);
                } else {
                    keepAlignment(outfile_cleaned_2stage, brec);
                }
            }
        }
        reader_s_uniq_map.bclose();
        if (verbose) GMessage("\n");

        /*********************************
         * Processing multi-mapped spliced alignments
        *********************************/
        GSamReader reader_s_multi_map(outfname_s_multi_map.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
        progressbar bar_multi(ALN_COUNT_SPLICED_MULTI);
        bar_multi.set_opening_bracket_char("[INFO] SPLAM! Filtering multi-mapped spliced alignments \n\t[");
        while ( (brec = reader_s_multi_map.next())!=NULL ) {
            if (verbose) {
                bar_multi.update();
            }
            bool spur = alignmentAssessment(brec, rm_juncs);
            if (spur) {
                update_flag_unpair_remove(brec);
                removeAlignment(outfile_discard_s_multi_map, brec, rm_hit);
                ALN_COUNT_SPLICED_MULTI_DISCARD += 1;
            } else {
                outfile_s_multi_map_tmp->write(brec);
            }
        }
        reader_s_multi_map.bclose();
        if (verbose) GMessage("\n");
    }

    for (auto ele : rm_hit) {
        GMessage("ele.first : %s\n", ele.first.c_str());
        GMessage("ele.second: %d\n\n", ele.second);
    }

    delete outfile_s_multi_map_tmp;
    delete outfile_discard_s_uniq_map;
    delete outfile_discard_s_multi_map;

    ALN_COUNT_BAD = ALN_COUNT_SPLICED_UNIQ_DISCARD + ALN_COUNT_SPLICED_MULTI_DISCARD + ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD + ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD;
    if (verbose) GMessage("[INFO] %d spurious alignments were removed.\n", ALN_COUNT_BAD);
    GStr outfname_NH_tag = writenhHitFile(rm_hit);
    return outfname_NH_tag;
}

void processRead(int currentstart, int currentend, GList<CReadAln>& readlist, BundleData& bdata, GHash<int>& hashread, CReadAln* alndata) { // some false positives should be eliminated here in order to break the bundle

	GSamRecord& brec=(alndata->brec);			   // bam record
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

void removeAlignment(GSamWriter* outfile_target, GSamRecord* brec, robin_hdd_rm_hit& rm_hit) {
    outfile_target->write(brec);
    std::string kv = brec->name();
    kv = kv + "_" + std::to_string(brec->pairOrder());
    if (rm_hit.find(kv) == rm_hit.end()) {
        rm_hit[kv] = 1;
    } else {
        rm_hit[kv]++;
    }
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
    // GMessage("%d - %d;  brec->hasIntrons(): %d\n", brec->start, brec->end, brec->hasIntrons());
    if (brec->hasIntrons()) {
        for (int e=1; e<brec->exons.Count(); e++) {
            char strand = brec->spliceStrand();
            // GMessage("\tIntron %d - %d\n", brec->exons[e-1].end, brec->exons[e].start-1);
            std::string jnew_sub = std::to_string(brec->exons[e-1].end) + "_" + std::to_string(brec->exons[e].start-1) + "_" + strand + "_" + brec->refName();
            // GMessage("\t jnew_sub: %s\n", jnew_sub.c_str());
            if (rm_juncs.find(jnew_sub) != rm_juncs.end()) {
                spur = true;
                return spur;
            }
        }
    }
    return spur;
}

std::string get_global_removed_mate_algns_key(GSamRecord* brec) {
    std::string brec_name = brec->name();
    int brec_start = brec->start;
    int brec_mate_start = brec->mate_start();
    int pair = 3 - brec->pairOrder();
    std::string key = brec_name + "_" + std::to_string(brec_mate_start) + "_" + std::to_string(brec_start) + "_" + std::to_string(pair);

    // GMessage("mate key: %s\n", key.c_str());
    return key;
}

void update_flag_paired_remove_both(GSamRecord* brec_1, GSamRecord* brec_2) {
    
    // GMessage("Before brec_1->flags(): %d (%s)\n", brec_1->flags(), brec_1->name());
    // First read in a pair
    int brec_1_flag_update = 0;
    if (brec_1->isPaired()) {
        brec_1_flag_update -= 1;
    }
    if (brec_1->isProperPaired()) {
        brec_1_flag_update -= 2;
    }
    if (brec_1->isMapped()) {
        brec_1_flag_update += 4;
    }
    if (brec_1->isMateMapped()) {
        brec_1_flag_update += 8;
    }
    brec_1->set_flags(brec_1->flags() + brec_1_flag_update);
    // GMessage("After brec_1->flags(): %d\n", brec_1->flags());


    // GMessage("Before brec_2->flags(): %d (%s)\n", brec_2->flags(), brec_2->name());
    // Second read in a pair
    int brec_2_flag_update = 0;
    if (brec_2->isPaired()) {
        brec_2_flag_update -= 1;
    }
    if (brec_2->isProperPaired()) {
        brec_2_flag_update -= 2;
    }
    if (brec_2->isMapped()) {
        brec_2_flag_update += 4;
    }
    if (brec_2->isMateMapped()) {
        brec_2_flag_update += 8;
    }
    brec_2->set_flags(brec_2->flags() + brec_2_flag_update);
    // GMessage("After brec_2->flags(): %d\n", brec_2->flags());
}

void update_flag_paired_remove_one(GSamRecord* removed, GSamRecord* kept) {
    // GMessage("Before removed->flags(): %d (%s)\n", removed->flags(), removed->name());
    // Removed alignment in a pair
    int removed_flag_update = 0;
    if (removed->isPaired()) {
        removed_flag_update -= 1;
    }
    if (removed->isProperPaired()) {
        removed_flag_update -= 2;
    }
    if (removed->isMapped()) {
        removed_flag_update += 4;
    }
    removed->set_flags(removed->flags() + removed_flag_update);
    // GMessage("After removed->flags(): %d\n", removed->flags());

    // Kept alignment in a pair
    int kept_flag_update = 0;
    if (kept->isPaired()) {
        kept_flag_update -= 1;
    }
    if (kept->isProperPaired()) {
        kept_flag_update -= 2;
    }
    if (kept->isMateMapped()) {
        kept_flag_update += 8;
    }
    kept->set_flags(kept->flags() + kept_flag_update);
}

void update_flag_unpair_remove(GSamRecord* removed) {

    // GMessage("Before removed->flags(): %d (%s)\n", removed->flags(), removed->name());
    int removed_flag_update = 0;
    if (removed->isPaired()) {
        removed_flag_update -= 1;
    }
    if (removed->isProperPaired()) {
        removed_flag_update -= 2;
    }
    // Unmap the current read
    if (removed->isMapped()) {
        removed_flag_update += 4;
    }
    removed->set_flags(removed->flags() + removed_flag_update);
    // GMessage("After removed->flags(): %d\n", removed->flags());
}

void update_flag_unpair_kept(GSamRecord* kept) {
    int kept_flag_update = 0;
    if (kept->isPaired()) {
        kept_flag_update -= 1;
    }
    if (kept->isProperPaired()) {
        kept_flag_update -= 2;
    }
    kept->set_flags(kept->flags() + kept_flag_update);
}
