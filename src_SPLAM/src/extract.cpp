#include "extract.h"
#include "common.h"
#include "extract.h"
#include "junc_func.h"
#include "util.h"
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <gclib/GBase.h>

/****************************
* Input : BAM file (s)
* Output: Junction bed file.
*****************************/
GStr splamJExtract() {
    STEP_COUNTER += 1;
    GMessage("###########################################\n");
    GMessage("## Step %d: generating spliced junctions in BED\n", STEP_COUNTER);
    GMessage("###########################################\n");

    int num_samples=in_records.start();

    // This is normal workflow for writing out all junctions.
    outfile_multimapped = new GSamWriter(outfname_multimapped, in_records.header(), GSamFile_BAM);
    outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);
    
    GStr outfname_junc_bed = out_dir + "/junction.bed";
    outfile_spliced = new GSamWriter(outfname_spliced, in_records.header(), GSamFile_BAM);

    GMessage("[INFO] Extracting junctions ...\n");
    GMessage("[INFO] Output directory\t\t: %s\n", out_dir.chars());
    GMessage("[INFO] Output Junction file\t: %s\n", outfname_junc_bed.chars());

    /****************************
    * Creating junction bed files.
    *****************************/
    // Creating the output junction bed file
    if (!outfname_junc_bed.is_empty()) {
        if (strcmp(outfname_junc_bed.substr(outfname_junc_bed.length()-4, 4).chars(), ".bed")!=0) {
            outfname_junc_bed.append(".bed");
        }
        joutf = fopen(outfname_junc_bed.chars(), "w");
        if (joutf==NULL) GError("Error creating file %s\n", outfname_junc_bed.chars());
    }

    // This is additional workflow to write out junctions above & below the thresholds.
    if (j_extract_threshold > 0) {
        GStr outfname_junc_above_bed = out_dir + "/junction_above.bed";
        GStr outfname_junc_below_bed = out_dir + "/junction_below.bed";

        outfile_above_spliced = new GSamWriter(outfname_junc_above_bed, in_records.header(), GSamFile_BAM);
        outfile_below_spliced = new GSamWriter(outfname_junc_below_bed, in_records.header(), GSamFile_BAM);

        GMessage("[INFO] Extracting junctions ...\n");
        GMessage("[INFO] Output directory\t\t: %s\n", out_dir.chars());
        GMessage("[INFO] Output Junction file\t: %s; %s\n", outfname_junc_above_bed.chars(), outfname_junc_below_bed.chars());

        // Creating the output junction bed file
        if (!outfname_junc_above_bed.is_empty()) {
            if (strcmp(outfname_junc_above_bed.substr(outfname_junc_above_bed.length()-4, 4).chars(), ".bed")!=0) {
                outfname_junc_above_bed.append(".bed");
            }
            joutf_above = fopen(outfname_junc_above_bed.chars(), "w");
            if (joutf_above==NULL) GError("Error creating file %s\n", outfname_junc_above_bed.chars());
        }
        if (!outfname_junc_below_bed.is_empty()) {
            if (strcmp(outfname_junc_below_bed.substr(outfname_junc_below_bed.length()-4, 4).chars(), ".bed")!=0) {
                outfname_junc_below_bed.append(".bed");
            }
            joutf_below = fopen(outfname_junc_below_bed.chars(), "w");
            if (joutf_below==NULL) GError("Error creating file %s\n", outfname_junc_below_bed.chars());
        }
    }

    /****************************
    * Iterating BAM file(s) and write out junctions.
    *****************************/
    // Reading BAM file.
    int prev_tid=-1;
    GStr prev_refname;
    int b_end=0, b_start=0;

    GMessage("[INFO] Processing BAM file ...\n");
    GMessage("\t\tBefore Hash map size: %d\n", read_hashmap.size());
    while ((irec=in_records.next())!=NULL) {
        brec=irec->brec;
        int endpos=brec->end;
        if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
            flushJuncs(joutf);
            if (j_extract_threshold > 0) {
                flushJuncs(joutf_above, joutf_below);
            }
            junctions.Clear();
            junctions.setCapacity(128);
            b_start=brec->start;
            b_end=endpos;
            prev_tid=brec->refId();
            prev_refname=(char*)brec->refName();
        } else { //extending current bundle
            if (b_end<endpos) {
                b_end=endpos;
            }
        }
        int accYC = 0;
        accYC = brec->tag_int("YC", 1);
        if (joutf && brec->exons.Count()>1) {
            // Spliced reads
            addJunction(*brec, accYC, prev_refname);
            outfile_spliced->write(brec);
            ALN_COUNT_SPLICED++;
        } else {
            // Non-spliced reads.
            // Not spliced => check their NH tags!
            if (brec->isUnmapped()) continue;
            int new_nh = brec->tag_int("NH", 0);
            if (new_nh == 1) {
                outfile_cleaned->write(brec);
            } else if (new_nh == 0){
                GMessage("\t\t brec->name(): %s !!!\n", brec->name());
                GMessage("\t\t NH tag is zero !!!: %d\n", new_nh);
            } else {
                outfile_multimapped->write(brec);
                ALN_COUNT_NH_UPDATE++;
            }
            ALN_COUNT_NSPLICED++;
        }
        ALN_COUNT++;
        if (ALN_COUNT % 1000000 == 0) {
            GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
        }
    }
    GMessage("\t\tAfter Hash map size: %d\n", read_hashmap.size());
    // for (auto it : read_hashmap) {
    //     std::cout << " " << it.first << ":" << "(NH tag) " << it.second.NH_tag_bound << std::endl;
    //     for (int i=0; i<it.second.sam_list.Count(); i++) {
    //         std::cout << it.second.sam_list[i].name() << std::endl;
    //     }
    // }
    GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
    in_records.stop();
    flushJuncs(joutf);
    if (j_extract_threshold > 0) {
        flushJuncs(joutf_above, joutf_below);
    }
    fclose(joutf);
    if (j_extract_threshold > 0) {
        fclose(joutf_above);
        fclose(joutf_below);
    }
    junctions.Clear();
    junctions.setCapacity(128);
    delete outfile_spliced;
    GMessage("[INFO] SPLAM! Total number of junctions: %d\n", JUNC_COUNT);	
    return outfname_junc_bed;
}
