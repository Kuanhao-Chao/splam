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

    int processed_aln = 0;
    /*********************************
     * Processing multip-mapped spliced temporary alignments (paired)
    *********************************/
    GSamReader reader_s_multi_map_tmp(outfname_s_multi_map_tmp.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    progressbar bar_s(ALN_COUNT_SPLICED_MULTI - ALN_COUNT_SPLICED_MULTI_DISCARD);
    bar_s.set_opening_bracket_char("[INFO] SPLAM! Processing multi-mapped spliced alignments\n\t[");

    while ((brec=reader_s_multi_map_tmp.next())!=NULL) {
        processed_aln += 1;
        if (verbose) {
            bar_s.update();
        }
        std::string kv = brec->name();
        kv = kv + "_" + std::to_string(brec->pairOrder());
        // GMessage("kv: %s\n", kv.c_str());
        if (rm_hit.find(kv) != rm_hit.end()) {
            // GMessage("\nkv: %s\n", kv.c_str());
            // GMessage("rm_hit[kv]: %d\n", rm_hit[kv]);
            // GMessage("Before update NH tag: %d\n", brec->tag_int("NH", 0));
            int new_nh = brec->tag_int("NH", 0) - rm_hit[kv];
            brec->add_int_tag("NH", new_nh);
            // GMessage("After update NH tag: %d\n", brec->tag_int("NH", 0));
        }
        keepAlignment(outfile_cleaned, brec);
    }
    reader_s_multi_map_tmp.bclose();
    GMessage("\n");

    /*********************************
     * Processing multip-mapped non-spliced alignments (paired)
    *********************************/
    GSamReader reader_ns_multi_map(outfname_ns_multi_map.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    progressbar bar_ns(ALN_COUNT_NSPLICED_MULTI);
    bar_ns.set_opening_bracket_char("[INFO] SPLAM! Processing multi-mapped nonspliced alignments\n\t[");

    while ((brec=reader_ns_multi_map.next())!=NULL) {
        processed_aln += 1;
        if (verbose) {
            bar_ns.update();
        }
        std::string kv = brec->name();
        kv = kv + "_" + std::to_string(brec->pairOrder());
        // GMessage("kv: %s\n", kv.c_str());
        if (rm_hit.find(kv) != rm_hit.end()) {
            // GMessage("\nkv: %s\n", kv.c_str());
            // GMessage("rm_hit[kv]: %d\n", rm_hit[kv]);
            // GMessage("Before update NH tag: %d\n", brec->tag_int("NH", 0));
            int new_nh = brec->tag_int("NH", 0) - rm_hit[kv];
            brec->add_int_tag("NH", new_nh);
            // GMessage("After update NH tag: %d\n", brec->tag_int("NH", 0));
        }
        keepAlignment(outfile_cleaned, brec);
    }
    reader_ns_multi_map.bclose();
    GMessage("\n");

    if (g_paired_removal) {
        /*********************************
         * Processing multip-mapped spliced temporary alignments (unpaired)
        *********************************/
        GSamReader reader_s_multi_unpair_tmp(outfname_s_multi_unpair_tmp.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

        progressbar bar_s_unpair(ALN_COUNT_SPLICED_MULTI_UNPAIR - ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD);
        bar_s_unpair.set_opening_bracket_char("[INFO] SPLAM! Processing multi-mapped spliced alignments (unpaired) \n\t[");

        while ((brec=reader_s_multi_unpair_tmp.next())!=NULL) {
            processed_aln += 1;
            if (verbose) {
                bar_s_unpair.update();
            }
            std::string kv = brec->name();
            kv = kv + "_" + std::to_string(brec->pairOrder());
            // GMessage("kv: %s\n", kv.c_str());
            if (rm_hit.find(kv) != rm_hit.end()) {
                // GMessage("\nkv: %s\n", kv.c_str());
                // GMessage("rm_hit[kv]: %d\n", rm_hit[kv]);
                // GMessage("Before update NH tag: %d\n", brec->tag_int("NH", 0));
                int new_nh = brec->tag_int("NH", 0) - rm_hit[kv];
                brec->add_int_tag("NH", new_nh);
                // GMessage("After update NH tag: %d\n", brec->tag_int("NH", 0));
            }
            keepAlignment(outfile_cleaned, brec);
        }
        reader_s_multi_unpair_tmp.bclose();
        GMessage("\n");

        /*********************************
         * Processing multip-mapped non-spliced alignments (unpaired)
        *********************************/
        GSamReader reader_ns_multi_unpair(outfname_ns_multi_unpair.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

        progressbar bar_ns_unpair(ALN_COUNT_NSPLICED_MULTI_UNPAIR);
        bar_ns_unpair.set_opening_bracket_char("[INFO] SPLAM! Processing multi-mapped nonspliced alignments (unpaired) \n\t[");

        while ((brec=reader_ns_multi_unpair.next())!=NULL) {
            processed_aln += 1;
            if (verbose) {
                bar_ns_unpair.update();
            }
            std::string kv = brec->name();
            kv = kv + "_" + std::to_string(brec->pairOrder());
            // GMessage("kv: %s\n", kv.c_str());
            if (rm_hit.find(kv) != rm_hit.end()) {
                // GMessage("\nkv: %s\n", kv.c_str());
                // GMessage("rm_hit[kv]: %d\n", rm_hit[kv]);
                // GMessage("Before update NH tag: %d\n", brec->tag_int("NH", 0));
                int new_nh = brec->tag_int("NH", 0) - rm_hit[kv];
                brec->add_int_tag("NH", new_nh);
                // GMessage("After update NH tag: %d\n", brec->tag_int("NH", 0));
            }
            keepAlignment(outfile_cleaned, brec);
        }
        reader_ns_multi_unpair.bclose();
        GMessage("\n");
    }

    delete outfile_cleaned;
    if (verbose) {
        GMessage("[INFO] %d alignments processed.\n", processed_aln);
        // ALN_COUNT_SPLICED_MULTI - ALN_COUNT_SPLICED_MULTI_DISCARD + 
        // ALN_COUNT_NSPLICED_MULTI +
        // ALN_COUNT_SPLICED_MULTI_UNPAIR - ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD + 
        // ALN_COUNT_NSPLICED_MULTI_UNPAIR);
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