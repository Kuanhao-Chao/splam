#include "clean.h"
#include "common.h"
#include "junc.h"
#include "junc_func.h"
#include "util.h"
#include "extract.h"
#include "predict.h"
// #include "progressbar.h"
#include <progressbar/progressbar.hpp>

#include <unordered_set>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <htslib/htslib/faidx.h>
#include <Python.h>
#include <gclib/GStr.h>

void splamClean(int argc, char* argv[]) {
    GStr outfname_junc_score = splamPredict();

    /*********************************************
     * Step 4: SPLAM filtering out reads.
    *********************************************/
    // GMessage("> outfname_discard: %s\n", outfname_discard.chars());
    // GMessage("> outfname_cleaned: %s\n", outfname_cleaned.chars());
    // GMessage("> outfname_spliced: %s\n", outfname_spliced.chars());
    // GMessage("> outfname_nspliced: %s\n", outfname_nspliced.chars());

    GMessage("\n###########################################\n");
    GMessage("## Step 4: SPLAM filtering out reads\n");
    GMessage("###########################################\n");
    robin_hdd_hm rm_rd_hm;
    std::unordered_set<std::string> rm_rd_set;

    GMessage(">> rm_rd_set size %d\n", rm_rd_set.size());
    GStr outfname_spliced_good = filterSpurJuncs(outfname_junc_score, rm_rd_hm, rm_rd_set);
    GMessage(">> rm_rd_set size %d\n", rm_rd_set.size());
    // for (auto ele : rm_rd_set) {
    //     GMessage("ele: %s \n", ele.c_str());
    // }


    GMessage("\n###########################################\n");
    GMessage("** Step 5: Updating NH tags in final clean BAM file\n");
    GMessage("###########################################\n");


    TInputFiles final_bam_records;
    final_bam_records.setup(VERSION, argc, argv);

    for (int i=0; i<in_records.freaders.Count(); i++) {
        GStr fname = in_records.freaders[i]->fname.chars();
        GMessage(">> fname: %s\n", fname.chars());
        final_bam_records.addFile(fname.chars());
    }
    // final_bam_records.addFile(get_full_path(outfname_nspliced.chars()).c_str());
    // final_bam_records.addFile(get_full_path(outfname_spliced_good.chars()).c_str());
    int num_samples=final_bam_records.start();

    outfile_cleaned = new GSamWriter(outfname_cleaned, final_bam_records.header(), GSamFile_BAM);
    outfile_discard = new GSamWriter(outfname_discard, final_bam_records.header(), GSamFile_BAM);

    // Reading BAM file.
    int counter = 0, prev_tid=-1;
    GStr prev_refname;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    int b_end=0, b_start=0;

    // int final_b = 0;
    // double percentage = 0;
    // double aln_count_good_counter = 0;
    // progressbar *progress = progressbar_new("Loading", ALN_COUNT_GOOD);
    progressbar bar(ALN_COUNT);
    bar.set_opening_bracket_char("[INFO] SPLAM! Updating NH tags \n\t[");

    while ((irec=final_bam_records.next())!=NULL) {
        bar.update();

        brec=irec->brec;
        int endpos=brec->end;

        std::string kv = brec->name();
        kv = kv + "_" + std::to_string(brec->pairOrder());



        char* seq = brec->sequence();
        char* cigar_seq = brec->cigar();

        std::string rm_rd_key = kv + "_" + seq + "_" + cigar_seq + "_" + std::to_string(brec->flags()) + "_" + std::to_string(brec->start);

        free(seq);
        free(cigar_seq);

        if (rm_rd_set.find(rm_rd_key) != rm_rd_set.end()) {
            // The alignment is found in the removed set.
            outfile_discard->write(brec);
        } else {
            if (rm_rd_hm.find(kv) != rm_rd_hm.end()) {
                // Update NH tag.
                // GMessage("Before updating NH tag: %d\n", brec->tag_int("NH", 0));   
                int new_nh = brec->tag_int("NH", 0) - rm_rd_hm[kv];
                // GMessage("New NH tag            : %d\n", new_nh);   

                brec->add_int_tag("NH", new_nh);
                // GMessage("After updating NH tag: %d\n\n\n", brec->tag_int("NH", 0));
            }
            outfile_cleaned->write(brec);
        }
    }
    GMessage("\n");

    final_bam_records.stop();

    delete outfile_discard;
    delete outfile_cleaned;
    // std::cout << "Done delete outfile_cleaned!" << std::endl;





    // extern int ALN_COUNT;
    // extern int ALN_COUNT_SPLICED;
    // extern int ALN_COUNT_NSPLICED;
    // extern int ALN_COUNT_GOOD;


    GMessage("\n\n[INFO] Total number of alignments\t:\t%d\n", ALN_COUNT);
    GMessage("[INFO]     spliced alignments\t\t:\t%d\n", ALN_COUNT_SPLICED);
    GMessage("[INFO]     non-spliced alignments\t:\t%d\n", ALN_COUNT_NSPLICED);
    GMessage("[INFO] Number of removed alignments\t:\t%d\n", ALN_COUNT_BAD);
    GMessage("[INFO] Number of kept alignments\t:\t%d\n", ALN_COUNT_GOOD);
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
        // GStr chr_str(chrname);
        if (junc[6].asDouble() <= threshold) {

            char* chrname =junc[0].detach();
            char* strand =junc[5].detach();
            // CJunc j(junc[1].asInt()+1, junc[2].asInt(), *junc[5].detach(), chr_str);
            std::string j = std::to_string(junc[1].asInt()+1) + "_" + std::to_string(junc[2].asInt()) + "_" + strand + "_" + chrname;

            // GMessage(">> j: %s\n", j.c_str());
            rm_juncs.insert(j);
            // std::cout << "rm_juncs.size: " << rm_juncs.Count() << std::endl;
        }
    }
    // std::cout << "bed_counter: " << bed_counter << std::endl;
}


GStr filterSpurJuncs(GStr outfname_junc_score, robin_hdd_hm &rm_rd_hm, std::unordered_set<std::string> &rm_rd_set) {
    GStr outfname_spliced_good;
    // GSamWriter* outfile_spliced_good = NULL;
    outfname_spliced_good = out_dir + "/TMP/spliced_good.bam";

    // outfile_spliced_good = new GSamWriter(outfname_spliced_good, in_records.header(), GSamFile_BAM);

    GSamReader bam_reader_spliced(outfname_spliced.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    int spur_cnt = 0;

    std::unordered_set<std::string> rm_juncs;
    GMessage("Before rm_juncs.size()  %d\n", rm_juncs.size());
    loadBed(outfname_junc_score, rm_juncs);
    GMessage("After rm_juncs.size()  %d\n", rm_juncs.size());

    int counter = 0, prev_tid=-1;
    GStr prev_refname;
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    int b_end=0, b_start=0;



    // extern int ALN_COUNT;
    // extern int ALN_COUNT_SPLICED;
    // extern int ALN_COUNT_NSPLICED;
    // extern int ALN_COUNT_GOOD;



    // int final_b = 0;
    double percentage = 0;
    // std::cout << "ALN_COUNT_SPLICED:" << ALN_COUNT_SPLICED << std::endl;

    progressbar bar(ALN_COUNT_SPLICED);
    bar.set_opening_bracket_char("[INFO] SPLAM! Removing junctions with low scores \n\t[");

    while ((brec=bam_reader_spliced.next())!=NULL) {
        bar.update();

        uint32_t dupcount=0;
        int endpos=brec->end;
        int r_exon_count = brec->exons.Count();
        bool spur = false;
        if (r_exon_count > 1) {
            for (int e=1; e<r_exon_count; e++) {
	            char strand = brec->spliceStrand();
                // CJunc jnew_sub(brec->exons[e-1].end+1, brec->exons[e].start-1, strand, GStr(brec->refName()));

                std::string jnew_sub = std::to_string(brec->exons[e-1].end+1) + "_" + std::to_string(brec->exons[e].start-1) + "_" + strand + "_" + brec->refName();

                // GMessage(">> jnew_sub: %s\n", jnew_sub.c_str());

                if (rm_juncs.find(jnew_sub) != rm_juncs.end()) {
                    spur = true;
                    break;
                }
            }
        }
        if (spur) {
            spur_cnt++;
            // std::cout << "~~ SPLAM!" << std::endl;
            std::string kv = brec->name();
            kv = kv + "_" + std::to_string(brec->pairOrder());
            if (rm_rd_hm.find(kv) == rm_rd_hm.end()) {
                rm_rd_hm[kv] = 1;
            } else {
                rm_rd_hm[kv]++;
            }
            char* seq = brec->sequence();
            char* cigar_seq = brec->cigar();

            kv = kv + "_" + seq + "_" + cigar_seq + "_" + std::to_string(brec->flags()) + "_" + std::to_string(brec->start);
            rm_rd_set.insert(kv);
            // GMessage(">> rm_rd_set size: %d;  ALN_COUNT_BAD: %d\n", rm_rd_set.size(), ALN_COUNT_BAD);

            ALN_COUNT_BAD++;
            free(seq);
            free(cigar_seq);
        } else {
            ALN_COUNT_GOOD++;
        }
    }
    GMessage("\n");
    ALN_COUNT_GOOD += ALN_COUNT_NSPLICED;
    GMessage("[INFO] %d spurious alignments were removed.\n", ALN_COUNT_BAD);
    return outfname_spliced_good;
}