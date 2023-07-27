/*  update.cpp -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#include "update.h"
#include <fstream>
#include <iostream>
#include <cctype>
#include <sstream>
#include <filesystem>
#include <progressbar/progressbar.hpp>

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

    if (g_paired_removal) {
        /*********************************
         * Processing multip-mapped spliced temporary alignments (paired)
        *********************************/
        update_NH_tag_write_alignment(outfname_s_multi_map_cleaned, outfile_s_multi_map_cleaned_nh_updated, processed_aln, rm_hit, "multi-mapped spliced alignments (paired)");

        /*********************************
         * Processing multip-mapped non-spliced alignments (paired)
        *********************************/
        update_NH_tag_write_alignment(outfname_ns_multi_map, outfile_ns_multi_map_nh_updated, processed_aln, rm_hit, "multi-mapped non-spliced alignments (paired)");

        // /*********************************
        //  * Processing uniq-mapped non-spliced alignments (paired)
        // *********************************/
        // update_NH_tag_write_alignment(outfname_ns_uniq_map, processed_aln, rm_hit, ALN_COUNT_NSPLICED_UNIQ);
        // ALN_COUNT_NSPLICED_UNIQ += 1;

        /*********************************
         * Processing multip-mapped spliced temporary alignments (unpaired)
        *********************************/
        update_NH_tag_write_alignment(outfname_s_multi_unpair_cleaned, outfile_s_multi_unpair_cleaned_nh_updated, processed_aln, rm_hit, "multi-mapped spliced alignments (unpaired)");

        /*********************************
         * Processing multip-mapped non-spliced alignments (unpaired)
        *********************************/
        update_NH_tag_write_alignment(outfname_ns_multi_unpair, outfile_ns_multi_unpair_nh_updated, processed_aln, rm_hit, "multi-mapped non-spliced alignments (unpaired)");

        // /*********************************
        //  * Processing uniq-mapped non-spliced alignments (unpaired)
        // *********************************/
        // update_NH_tag_write_alignment(outfname_ns_uniq_unpair, processed_aln, rm_hit, ALN_COUNT_NSPLICED_UNIQ);
        // ALN_COUNT_NSPLICED_UNIQ_UNPAIR += 1;

    } else {
        /*********************************
         * Processing multip-mapped spliced temporary alignments (unpaired)
        *********************************/
        update_NH_tag_write_alignment(outfname_s_multi_map_cleaned, outfile_s_multi_map_cleaned_nh_updated, processed_aln, rm_hit, "multi-mapped spliced alignments");

        /*********************************
         * Processing multip-mapped non-spliced alignments (unpaired)
        *********************************/
        update_NH_tag_write_alignment(outfname_ns_multi_map, outfile_ns_multi_map_nh_updated, processed_aln, rm_hit, "multi-mapped non-spliced alignments");

        /*********************************
         * Processing uniq-mapped non-spliced alignments (unpaired)
        *********************************/
        // update_NH_tag_write_alignment(outfname_ns_uniq_map, processed_aln, rm_hit, ALN_COUNT_NSPLICED_UNIQ);
        // ALN_COUNT_NSPLICED_UNIQ += 1;
    }

    if (verbose) {
        GMessage("[INFO] %d alignments processed.\n", 
        ALN_COUNT_SPLICED_MULTI - ALN_COUNT_SPLICED_MULTI_DISCARD + 
        ALN_COUNT_NSPLICED_MULTI +
        ALN_COUNT_SPLICED_MULTI_UNPAIR - ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD + 
        ALN_COUNT_NSPLICED_MULTI_UNPAIR);
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

void update_NH_tag_write_alignment(GStr infname, GSamWriter *outfile, int& processed_aln, robin_hdd_rm_hit rm_hit, GStr log) {
    GSamReader reader(infname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    GStr msg = "[INFO] SPLAM! Processing " + log + "\n";
    fprintf(stderr, msg);
    while ((brec=reader.next())!=NULL) {
        processed_aln += 1;
        std::string kv = brec->name();
        kv = kv + "_" + std::to_string(brec->pairOrder());
        // GMessage("kv: %s\n", kv.c_str());
        if (rm_hit.find(kv) != rm_hit.end()) {
            // GMessage("\nkv: %s\n", kv.c_str());
            // GMessage("rm_hit[kv]: %d\n", rm_hit[kv]);
            // GMessage("Before update NH tag: %d\n", brec->tag_int("NH", 0));
            int new_nh = brec->tag_int("NH", 0) - rm_hit[kv];
            brec->add_int_tag("NH", new_nh);
        }
        int tmp = 0;
        keepAlignment(outfile, brec, tmp);
    }
    reader.bclose();
    delete outfile;
}