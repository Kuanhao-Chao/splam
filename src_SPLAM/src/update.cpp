#include "update.h"
#include <fstream>
#include <iostream>
#include <cctype>
#include <sstream>
#include <filesystem>

// #include "predict.h"
// #include "extract.h"
#include "common.h"
#include "util.h"
#include "clean.h"

GStr splamNHUpdate(int argc, char* argv[]) {

    splamClean(argc, argv);

    GMessage("\n###########################################\n");
    GMessage("** Step 5: Updating NH tags in final clean BAM file\n");
    GMessage("###########################################\n");
    std::unordered_map<std::string, int> nh_hm;

    GStr NH_tag_fname = out_dir + "/NH_tag_fix.csv";
    std::ifstream ref_f(NH_tag_fname.chars());
    std::string line;
    while(getline(ref_f, line)){
        std::string read_id;
        int hits;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::stringstream ss(line);
        ss >> read_id;
        ss >> hits;
        nh_hm[read_id] = hits;
    }   
    // for (auto i : nh_hm) {
    //   std::cout << i.first << " ---- " << i.second << std::endl;
    // }

    GStr outfname_fix = out_dir + "/cleaned.fix.bam";
    GSamWriter* outfile_fix = new GSamWriter(outfname_fix, in_records.header(), GSamFile_BAM);
    GSamReader bam_reader_cleaned(outfname_cleaned.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    // Reading BAM file.
    // int prev_tid=-1;
    // GStr prev_refname;
    // int b_end=0, b_start=0;
    
    int bam_clean_counter=0;
    GMessage("[INFO] Processing BAM file ...\n");
    while ((brec=bam_reader_cleaned.next())!=NULL) {
        int endpos=brec->end;

        std::string kv = brec->name();
        kv = kv + "_" + std::to_string(brec->pairOrder());
        if (nh_hm.find(kv) != nh_hm.end()) {
            GMessage("\tBefore update NH: %d\n", brec->tag_int("NH", 0));
            int new_nh = brec->tag_int("NH", 0) - nh_hm[kv];
            brec->add_int_tag("NH", new_nh);
            GMessage("\tAfter update NH: %d\n\n", brec->tag_int("NH", 0));
        }
        // if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
        //     b_start=brec->start;
        //     b_end=endpos;
        //     prev_tid=brec->refId();
        //     prev_refname=(char*)brec->refName();
        // } else { //extending current bundle
        //     if (b_end<endpos) {
        //         b_end=endpos;
        //     }
        // }
        outfile_fix->write(brec);
        bam_clean_counter++;
        if (bam_clean_counter % 1000000 == 0) {
            GMessage("\t\t%d alignments processed.\n", bam_clean_counter);
        }
    }
    bam_reader_cleaned.bclose();
    GMessage("\t\t%d alignments processed.\n", bam_clean_counter);
    delete outfile_fix;
    return outfname_fix;
}