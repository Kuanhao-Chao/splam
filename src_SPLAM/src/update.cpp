#include "update.h"
#include <fstream>
#include <iostream>
#include <cctype>
#include <sstream>
#include <filesystem>
#include <progressbar/progressbar.hpp>

// #include "predict.h"
// #include "extract.h"
#include "common.h"
#include "util.h"
#include "clean.h"

/****************************
* Input : (1) unordered_set of reads, (2) hashmap of removed hits.
* Output: (1) cleamed BAM file.
*****************************/
GStr splamNHUpdate() {
    STEP_COUNTER += 1;
    if (verbose) {
        GMessage("\n###########################################\n");
        GMessage("** Step %d: Updating NH tags in final clean BAM file\n", STEP_COUNTER);
        GMessage("###########################################\n");
    }
    robin_hdd_rm_hit rm_hit;
    readnhHitFile(rm_hit);

        // for (auto const& x : rm_hit)
        // {
        //     std::cout << x.first  // string (key)
        //             << ':' 
        //             << x.second // string's value 
        //             << std::endl;
        // }

    if (COMMAND_MODE == CLEAN) {
        // outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);
        // GSamReader bam_reader_nh_updater(outfname_cleaned_tmp.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
        // int bam_clean_counter=0;
        // GMessage("[INFO] Processing BAM file ...\n");

        // progressbar bar(ALN_COUNT_GOOD);
        // bar.set_opening_bracket_char("[INFO] SPLAM! Output the final clean BAM file \n\t[");
        // // while ((irec=final_bam_records.next())!=NULL) {
        // while ((brec=bam_reader_nh_updater.next())!=NULL) {

        //     bar.update();
        //     std::string kv = brec->name();
            
        //     kv = kv + "_" + std::to_string(brec->pairOrder());

        //     // GMessage("kv: %s\n", kv.c_str());

        //     if (rm_hit.find(kv) != rm_hit.end()) {
        //         // GMessage("Original brec NH tag: %d\n", brec->tag_int("NH", 0));
        //         int new_nh = brec->tag_int("NH", 0) - rm_hit[kv];
        //         brec->add_int_tag("NH", new_nh);
        //         // GMessage("After    brec NH tag: %d\n", brec->tag_int("NH", 0));
        //     }
        //     outfile_cleaned->write(brec);
        // }
        // delete outfile_cleaned;
    } else if (COMMAND_MODE == ALL) {

        // TInputFiles final_bam_records;
        // final_bam_records.setup(VERSION, argc, argv);
        // for (int i=0; i<in_records.freaders.Count(); i++) {
        //     GStr fname = in_records.freaders[i]->fname.chars();
        //     GMessage(">> fname: %s\n", fname.chars());
        //     final_bam_records.addFile(fname.chars());
        // }
        // int num_samples=final_bam_records.start();
        // outfile_cleaned = new GSamWriter(outfname_cleaned, final_bam_records.header(), GSamFile_BAM);
        // outfile_discard = new GSamWriter(outfname_discard, final_bam_records.header(), GSamFile_BAM);

        int bam_clean_counter=0;
        GSamReader reader_s_multi_map_tmp(outfname_s_multi_map_tmp.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

        progressbar bar_s(ALN_COUNT_SPLICED_MULTI - ALN_COUNT_SPLICED_MULTI_DISCARD);
        bar_s.set_opening_bracket_char("[INFO] SPLAM! Processing multi-mapped spliced alignments \n\t[");
        // while ((irec=final_bam_records.next())!=NULL) {
        while ((brec=reader_s_multi_map_tmp.next())!=NULL) {
            if (verbose) {
                bar_s.update();
            }
            std::string kv = brec->name();
            // GMessage("kv: %s\n", kv.c_str());
            if (rm_hit.find(kv) != rm_hit.end()) {
                // GMessage("rm_hit[kv]: %d\n", rm_hit[kv]);
                // GMessage("Before update NH tag: %d\n", brec->tag_int("NH", 0));
                int new_nh = brec->tag_int("NH", 0) - rm_hit[kv];
                brec->add_int_tag("NH", new_nh);
                // GMessage("After update NH tag: %d\n", brec->tag_int("NH", 0));
            }
            outfile_cleaned->write(brec);   
            ALN_COUNT_GOOD++;
        }
        GMessage("\n");


        GSamReader reader_ns_multi_map(outfname_ns_multi_map .chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

        progressbar bar_ns(ALN_COUNT_NSPLICED_MULTI);
        bar_ns.set_opening_bracket_char("[INFO] SPLAM! Processing multi-mapped nonspliced alignments \n\t[");
        // while ((irec=final_bam_records.next())!=NULL) {
        while ((brec=reader_ns_multi_map.next())!=NULL) {
            if (verbose) {
                bar_ns.update();
            }
            std::string kv = brec->name();
            // GMessage("kv: %s\n", kv.c_str());
            if (rm_hit.find(kv) != rm_hit.end()) {
                int new_nh = brec->tag_int("NH", 0) - rm_hit[kv];

                // GMessage("Before update NH tag: %d\n", new_nh);
                brec->add_int_tag("NH", new_nh);
                // GMessage("After update NH tag: %d\n", brec->tag_int("NH", 0));
            }
            outfile_cleaned->write(brec); 
            ALN_COUNT_GOOD++;  
        }
        GMessage("\n");
        delete outfile_cleaned;

        if (verbose) {
            GMessage("[INFO] %d alignments processed.\n", ALN_COUNT_SPLICED_MULTI - ALN_COUNT_SPLICED_MULTI_DISCARD + ALN_COUNT_NSPLICED_MULTI);
        }
    }
    return outfname_cleaned;
}

void readnhHitFile(robin_hdd_rm_hit& rm_hit) {
    GStr NH_tag_fname = out_dir + "/TMP/NH_tag_fix.csv";
    std::ifstream ref_f(NH_tag_fname.chars());
    std::string line;
    while(getline(ref_f, line)){
        std::string read_id;
        int hits;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::stringstream ss(line);
        ss >> read_id;
        ss >> hits;
        rm_hit[read_id] = hits;
    }   
}