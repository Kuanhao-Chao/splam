#include "clean.h"
#include "common.h"
#include "junc.h"
#include "junc_func.h"
#include "util.h"
#include "extract.h"
#include "predict.h"
#include "progressbar.h"

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

    GMessage("\n********************************************\n");
    GMessage("** Step 4: SPLAM filtering out reads\n");
    GMessage("********************************************\n");
    robin_hdd_hm rm_rd_hm;
    GStr outfname_spliced_good = filterSpurJuncs(outfname_junc_score, rm_rd_hm);

    GMessage("\n********************************************\n");
    GMessage("** Step 5: Updating NH tags in final clean BAM file\n");
    GMessage("********************************************\n");
    TInputFiles final_bam_records;
    final_bam_records.setup(VERSION, argc, argv);
    final_bam_records.addFile(get_full_path(outfname_nspliced.chars()).c_str());
    final_bam_records.addFile(get_full_path(outfname_spliced_good.chars()).c_str());
    int num_samples=final_bam_records.start();
    outfile_cleaned = new GSamWriter(outfname_cleaned, final_bam_records.header(), GSamFile_BAM);

    // Reading BAM file.
    int counter = 0, prev_tid=-1;
    GStr prev_refname;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    int b_end=0, b_start=0;

    // int final_b = 0;
    double percentage = 0;
    double aln_count_good_counter = 0;
    while ((irec=final_bam_records.next())!=NULL) {
        aln_count_good_counter ++;
        percentage = aln_count_good_counter/ALN_COUNT_GOOD;
        printProgress(percentage);
        brec=irec->brec;
        int endpos=brec->end;

        std::string kv = brec->name();
        kv = kv + ";" + std::to_string(brec->pairOrder());
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
    GMessage("\n");

    final_bam_records.stop();
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


void loadBed(GStr inbedname, GArray<CJunc> &spur_juncs) {

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
        char* chrname =junc[0].detach();
        GStr chr_str(chrname);
        if (junc[6].asDouble() <= threshold) {
            CJunc j(junc[1].asInt()+1, junc[2].asInt(), *junc[5].detach(), chr_str);
            spur_juncs.Add(j);
            // std::cout << "spur_juncs.size: " << spur_juncs.Count() << std::endl;
        }
    }
    // std::cout << "bed_counter: " << bed_counter << std::endl;
}


GStr filterSpurJuncs(GStr outfname_junc_score, robin_hdd_hm &rm_rd_hm) {
    GStr outfname_spliced_good;
    GSamWriter* outfile_spliced_good = NULL;
    outfname_spliced_good = out_dir + "/TMP/spliced_good.bam";

    outfile_spliced_good = new GSamWriter(outfname_spliced_good, in_records.header(), GSamFile_BAM);
    outfile_discard = new GSamWriter(outfname_discard, in_records.header(), GSamFile_BAM);

    GSamReader bam_reader_spliced(outfname_spliced.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    auto start=std::chrono::high_resolution_clock::now();
    int spur_cnt = 0;

    GArray<CJunc> rm_juncs;
    loadBed(outfname_junc_score, rm_juncs);

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
    int aln_spliced_counter = 0;

    // std::cout << "ALN_COUNT_SPLICED:" << ALN_COUNT_SPLICED << std::endl;
    while ((brec=bam_reader_spliced.next())!=NULL) {
        aln_spliced_counter++;
        GMessage("aln_spliced_counter: %d;  ALN_COUNT_SPLICED: %d\n", aln_spliced_counter, ALN_COUNT_SPLICED);
        percentage = aln_spliced_counter/ALN_COUNT_SPLICED;
        printProgress(percentage);

        uint32_t dupcount=0;
        int endpos=brec->end;
        int r_exon_count = brec->exons.Count();
        bool spur = false;
        if (r_exon_count > 1) {
            for (int e=1; e<r_exon_count; e++) {
	            char strand = brec->spliceStrand();
                CJunc jnew_sub(brec->exons[e-1].end+1, brec->exons[e].start-1, strand, GStr(brec->refName()));
                if (rm_juncs.Exists(jnew_sub)) {
                    spur = true;
                    break;
                }
            }
        }
        if (spur) {
            spur_cnt++;
            // std::cout << "~~ SPLAM!" << std::endl;
            std::string kv = brec->name();
            kv = kv + ";" + std::to_string(brec->pairOrder());
            if (rm_rd_hm.find(kv) == rm_rd_hm.end()) {
                rm_rd_hm[kv] = 1;
            } else {
                rm_rd_hm[kv]++;
            }
            ALN_COUNT_BAD++;
            outfile_discard->write(brec);
        } else {
            ALN_COUNT_GOOD++;
            outfile_spliced_good->write(brec);
        }
    }
    GMessage("\n");

    ALN_COUNT_GOOD += ALN_COUNT_NSPLICED;
    
    delete outfile_discard;
    // std::cout << "Done delete outfile_discard!" << std::endl;
    delete outfile_spliced_good;
    // std::cout << "Done delete outfile_spliced_good!" << std::endl;


    auto end =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    GMessage("[INFO] %d spurious alignments were removed.\n", spur_cnt);
    GMessage("[INFO] Completed in %f seconds.\n", duration.count());

    // std::cout << spur_cnt << " spurious alignments were removed." << std::endl;
    // std::cout << "Completed in " << duration.count() << " seconds" << std::endl;
    
    return outfname_spliced_good;
}