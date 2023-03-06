#include "extract.h"
#include "common.h"
#include "extract.h"
#include "junc_func.h"
#include "util.h"
#include "bundle.h"
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <gclib/GBase.h>
#include <gclib/GHashMap.hh>

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
    // outfile_multimapped = new GSamWriter(outfname_multimapped, in_records.header(), GSamFile_BAM);
    outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);
    
    GStr outfname_junc_bed = out_dir + "/junction.bed";
    // outfile_spliced = new GSamWriter(outfname_spliced, in_records.header(), GSamFile_BAM);

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




    outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);
    outfile_ns_multi_map = new GSamWriter(outfname_ns_multi_map, in_records.header(), GSamFile_BAM);
    outfile_s_uniq_map = new GSamWriter(outfname_s_uniq_map, in_records.header(), GSamFile_BAM);
    outfile_s_multi_map = new GSamWriter(outfname_s_multi_map, in_records.header(), GSamFile_BAM);
    outfile_discard_unpair = new GSamWriter(outfname_discard_unpair, in_records.header(), GSamFile_BAM);


    
    BundleData* bundle = new BundleData();
	GList<CReadAln> readlist;

    uint runoffdist=100;
    GHash<int> hashread; //read_name:pos:hit_index => readlist index
    GStr lastref;
    bool more_alns = true;
    int prev_pos = 0;
    int lastref_id = -1; //last seen gseq_id

    bool fr_strand = false;
    bool rf_strand = false;
    int currentstart = 0, currentend = 0;
    int bundle_counter = 0;

    int max_s_distance = 20000;
        // while (more_alns) {
        //     bool chr_changed=false;
        //     int pos=0;
        //     const char* refseqName=NULL;
        //     char xstrand=0;
        //     int nh=1;
        //     int hi=0;
        //     int gseq_id=lastref_id;  //current chr id
        //     bool new_bundle=false;
        //     //delete brec;
        //     if ((irec=in_records.next())!=NULL) {
        //         brec=irec->brec;

        //         /***********************************
        //          * Setting the "chr" "strand" of the current alignment.
        //         ************************************/
        //         refseqName=brec->refName();
        //         xstrand=brec->spliceStrand(); // tagged strand gets priority
        //         if(xstrand=='.' && (fr_strand || rf_strand)) { // set strand if stranded library
        //             if(brec->isPaired()) { // read is paired
        //                 if(brec->pairOrder()==1) { // first read in pair
        //                     if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
        //                     else xstrand='-';
        //                 }
        //                 else {
        //                     if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='-';
        //                     else xstrand='+';
        //                 }
        //             }
        //             else {
        //                 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
        //                 else xstrand='-';
        //             }
        //         }

        //         /***********************************
        //          * Setting the "chr_changed" and "new_bundle" parameters.
        //         ************************************/
        //         pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
        //         chr_changed=(lastref.is_empty() || lastref!=refseqName);
        //         if (chr_changed) {
        //             prev_pos=0;
        //         }

        //         if (pos == 0) {
        //             // This is an unmapped read
        //         } else if (pos<prev_pos) {
        //             GMessage("[ERROR] %s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
        //             brec->name(), brec->start,  pos, refseqName, prev_pos);
        //             exit(-1);
        //         }
        //         prev_pos=pos;
        //         nh=brec->tag_int("NH", 0);
        //         if (nh==0) nh=1;
        //         hi=brec->tag_int("HI", 0);
        //         if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
        //             new_bundle=true;
        //         }
        //     } else { //no more alignments
        //         more_alns=false;
        //         new_bundle=true; //fake a new start (end of last bundle)
        //     }

        //     /***********************************
        //      * Process the bundle!
        //     ************************************/
        //     if (new_bundle || chr_changed) {
        //         hashread.Clear();
        //         if (readlist.Count()>0) {
        //             // process reads in previous bundle
        //             bundle->getReady(currentstart, currentend);

        //             GMessage(">> bundle read count: %d\n", readlist.Count());
        //             GMessage(">> bundle start     : %d\n", bundle->start);
        //             GMessage(">> bundle end       : %d\n", bundle->end);

        //             processBundle(bundle, readlist, rm_juncs, rm_hit, bundle_counter);
        //             readlist.Clear();
        //         } else { 
        //             //no read alignments in this bundle?  
        //             bundle->Clear();
        //             readlist.Clear();
        //         } //nothing to do with this bundle

        //         if (chr_changed) {
        //             lastref = refseqName;
        //             lastref_id = gseq_id;
        //             currentend = 0;
        //         }

        //         if (!more_alns) {
        //             noMoreBundles();
        //             break;
        //         }

        //         if (brec->start > 0) {
        //             currentstart = pos;
        //             currentend = brec->end;
        //         }
        //         bundle->refseq = lastref;
        //         bundle->start = currentstart;
        //         bundle->end = currentend;
        //     } //<---- new bundle started

        //     int fragment_end = 0;
        //     if (brec->refId() == brec->mate_refId()) {

        //         int insert_size = brec->insertSize();
        //         int mate_end = brec->mate_start() + insert_size;
                
        //         if (mate_end > (int)brec->end && insert_size <= max_s_distance) {
        //             fragment_end = mate_end;
        //         } else {
        //             fragment_end = (int)brec->end;
        //         }
        //     } else {
        //         fragment_end = (int)brec->end;
        //     }

        //     if (currentend<fragment_end) {
        //         //current read extends the bundle
        //         currentend=fragment_end;
        //     } //adjusted currentend and checked for overlapping reference transcripts

        //     // GMessage("brec->refName(): %s\n", brec->refName());
        //     CReadAln* alndata = new CReadAln(brec);
        //     processRead(currentstart, currentend, readlist, *bundle, hashread, alndata);
        // } //for each read alignment






































    // This is additional workflow to write out junctions above & below the thresholds.
    if (g_j_extract_threshold > 0) {
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
    // GMessage("\t\tBefore Hash map size: %d\n", read_hashmap.size());
    while ((irec=in_records.next())!=NULL) {
        brec=irec->brec;
        int endpos=brec->end;
        if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
            flushJuncs(joutf);
            if (g_j_extract_threshold > 0) {
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
            // outfile_spliced->write(brec);
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
                // outfile_multimapped->write(brec);
                ALN_COUNT_NH_UPDATE++;
            }
            ALN_COUNT_NSPLICED++;
        }
        ALN_COUNT++;
        if (ALN_COUNT % 1000000 == 0) {
            GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
        }
    }
    // GMessage("\t\tAfter Hash map size: %d\n", read_hashmap.size());
    // for (auto it : read_hashmap) {
    //     std::cout << " " << it.first << ":" << "(NH tag) " << it.second.NH_tag_bound << std::endl;
    //     for (int i=0; i<it.second.sam_list.Count(); i++) {
    //         std::cout << it.second.sam_list[i].name() << std::endl;
    //     }
    // }
    GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
    in_records.stop();
    flushJuncs(joutf);
    if (g_j_extract_threshold > 0) {
        flushJuncs(joutf_above, joutf_below);
    }
    fclose(joutf);
    if (g_j_extract_threshold > 0) {
        fclose(joutf_above);
        fclose(joutf_below);
    }
    junctions.Clear();
    junctions.setCapacity(128);
    // delete outfile_spliced;
    GMessage("[INFO] SPLAM! Total number of junctions: %d\n", JUNC_COUNT);	
    return outfname_junc_bed;
}
