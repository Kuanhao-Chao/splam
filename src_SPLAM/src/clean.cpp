#include "clean.h"
#include "common.h"
#include "junc.h"
#include "junc_func.h"
#include "util.h"
#include "extract.h"
#include "predict.h"

#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <htslib/htslib/faidx.h>
#include <Python.h>
#include <gclib/GStr.h>
#include <robin_hood/robin_hood.h>

typedef robin_hood::unordered_map<std::string, int> robin_hdd_hm;

void splamClean(int argc, char* argv[]) {
    GStr outfname_junc_score = splamPredict();

    /*********************************************
     * Step 4: SPLAM filtering out reads.
    *********************************************/
    // GMessage("> outfname_discard: %s\n", outfname_discard.chars());
    // GMessage("> outfname_cleaned: %s\n", outfname_cleaned.chars());
    // GMessage("> outfname_spliced: %s\n", outfname_spliced.chars());
    // GMessage("> outfname_nspliced: %s\n", outfname_nspliced.chars());

    GStr outfname_spliced_good;
    GStr outfname_spliced_bad;

    GSamWriter* outfile_spliced_good = NULL;
    GSamWriter* outfile_spliced_bad = NULL;

    outfname_spliced_good = out_dir + "/TMP/spliced_good.bam";
    outfname_spliced_bad = out_dir + "/TMP/spliced_bad.bam";

    outfile_spliced_good = new GSamWriter(outfname_spliced_good, in_records.header(), GSamFile_BAM);
    outfile_spliced_bad = new GSamWriter(outfname_spliced_bad, in_records.header(), GSamFile_BAM);

    GMessage("\n********************************************\n");
    GMessage("** Step 4: SPLAM filtering out reads\n");
    GMessage("********************************************\n");
    GSamReader bam_reader_spliced(outfname_spliced.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    auto start=std::chrono::high_resolution_clock::now();
    int spur_cnt = 0;

    GArray<CJunc> rm_juncs;
    loadBed(outfname_junc_score, rm_juncs);

    int counter = 0, prev_tid=-1;
    GStr prev_refname;
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    int b_end=0, b_start=0;

    robin_hdd_hm rm_rd_hm;

    int bam_counter = 0;
    while ((brec=bam_reader_spliced.next())!=NULL) {
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
            // std::cout << "rm_juncs.Exists! " << std::endl;
            spur_cnt++;
            // std::cout << "~~ SPLAM!" << std::endl;
            std::string kv = brec->name();
            kv = kv + ";" + std::to_string(brec->pairOrder());
            if (rm_rd_hm.find(kv) == rm_rd_hm.end()) {
                rm_rd_hm[kv] = 1;
            } else {
                rm_rd_hm[kv]++;
            }
            outfile_spliced_bad->write(brec);
        } else {
            outfile_spliced_good->write(brec);
        }
        bam_counter++;
    }

    for (auto i : rm_rd_hm)
        std::cout << i.first << "   " << i.second
             << std::endl;
    GMessage(">> bam_counter: %d\n", bam_counter);

    
    delete outfile_spliced_bad;
    std::cout << "Done delete outfile_spliced_bad!" << std::endl;
    delete outfile_spliced_good;
    std::cout << "Done delete outfile_spliced_good!" << std::endl;






    auto end =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << spur_cnt << " spurious alignments were removed." << std::endl;
    std::cout << "Completed in " << duration.count() << " seconds" << std::endl;
    
    
    delete outfile_discard;
    std::cout << "Done delete outfile_discard!" << std::endl;
    delete outfile_cleaned;
    std::cout << "Done delete outfile_cleaned!" << std::endl;
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
            // std::cout << "junc[6].asDouble(): " << junc[6].asDouble() << std::endl;

        	// CJunc(int vs=0, int ve=0, char vstrand='+', std::string vref=".", uint64_t dcount=1):
            CJunc j(junc[1].asInt()+1, junc[2].asInt(), *junc[5].detach(), chr_str);
            spur_juncs.Add(j);
            // std::cout << "spur_juncs.size: " << spur_juncs.Count() << std::endl;
        }
    }
    // std::cout << "bed_counter: " << bed_counter << std::endl;
}

