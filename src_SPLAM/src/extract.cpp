#include "extract.h"
#include "common.h"
#include "extract.h"
#include "junc_func.h"
#include "util.h"
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <gclib/GBase.h>

GStr splamJExtract() {
    int num_samples=in_records.start();
    GStr outfname_junc_bed = out_dir + "/junction.bed";
    
    outfile_spliced = new GSamWriter(outfname_spliced, in_records.header(), GSamFile_BAM);
    outfile_nspliced = new GSamWriter(outfname_nspliced, in_records.header(), GSamFile_BAM);

    GMessage("[INFO] Extracting junctions ...\n");
    GMessage("[INFO] Number of samples\t: %d\n", num_samples);
    GMessage("[INFO] Output directory\t\t: %s\n", out_dir.chars());
    GMessage("[INFO] Output Junction file\t: %s\n", outfname_junc_bed.chars());

    // Creating the output junction bed file
    if (!outfname_junc_bed.is_empty()) {
        if (strcmp(outfname_junc_bed.substr(outfname_junc_bed.length()-4, 4).chars(), ".bed")!=0) {
            outfname_junc_bed.append(".bed");
        }
        joutf = fopen(outfname_junc_bed.chars(), "w");
        if (joutf==NULL) GError("Error creating file %s\n", outfname_junc_bed.chars());
        // fprintf(joutf, "track name=junctions\n");
    }

    // Reading BAM file.
    int prev_tid=-1;
    GStr prev_refname;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    int b_end=0, b_start=0;

    GMessage("[INFO] Processing BAM file ...\n");
    while ((irec=in_records.next())!=NULL) {
        brec=irec->brec;
        uint32_t dupcount=0;
        std::vector<int> cur_samples;
        int endpos=brec->end;
        if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
            if (joutf) {
                flushJuncs(joutf);
            } // TODO: write the last column to 3 dec places
            b_start=brec->start;
            b_end=endpos;
            prev_tid=brec->refId();
            prev_refname=(char*)brec->refName();
        } else { //extending current bundle
            if (b_end<endpos) {
                b_end=endpos;
                bcov.setCount(b_end-b_start+1, (int)0);
            }
        }
        int accYC = 0;
        accYC = brec->tag_int("YC", 1);
        if (joutf && brec->exons.Count()>1) {
            addJunction(*brec, accYC, prev_refname);
            outfile_spliced->write(brec);
            ALN_COUNT_SPLICED++;
        } else {
            outfile_nspliced->write(brec);
            ALN_COUNT_NSPLICED++;
        }
        ALN_COUNT++;
        if (ALN_COUNT % 1000000 == 0) {
            GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
        }
    }
    GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
    in_records.stop();
    flushJuncs(joutf);
    fclose(joutf);

    delete outfile_spliced;
    // std::cout << "Done delete outfile_spliced!" << std::endl;
    delete outfile_nspliced;
    // std::cout << "Done delete outfile_nspliced!" << std::endl;

    GMessage("[INFO] SPLAM! Total number of junctions: %d\n", JUNC_COUNT);	
    return outfname_junc_bed;
// extern GSamWriter* outfile_discard;
// extern GSamWriter* outfile_spliced;
// extern GSamWriter* outfile_cleaned;
}
