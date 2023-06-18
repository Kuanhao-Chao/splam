// #define DEBUG
#include <iostream>
#include <vector>
#include <memory>
#include <filesystem>

#include "extract.h"
#include "predict.h"
#include "clean.h"

#include "common.h"
#include "tmerge.h"
#include "util.h"
#include "junc.h"
#include "update.h"

#include <gclib/GArgs.h>
#include <gclib/GBase.h>
#include <gclib/GStr.h>
#include <robin_hood/robin_hood.h>
#include <Python.h>

void processOptions(int argc, char* argv[]);
void processOptionsClean(GArgs& args);

void optionsOutput(GArgs& args);
void optionsWriteTMP(GArgs& args);
// void options2StageRun(GArgs& args);

CommandMode COMMAND_MODE = UNSET;
GStr command_str;

// input file names 
GStr infname_model_name("");
GStr infname_reffa("");
GStr infname_bam("");
GStr infname_juncbed("");
GStr infname_scorebed("");
GStr infname_NH_tag("");
GStr out_dir;

bool verbose = false;
TInputFiles in_records;
TInputRecord* irec=NULL;

float threshold = 0.1;
int aln_num_thr = 4;
GSamRecord* brec=NULL;

// output file names 
GStr outfname_cleaned;

// Paired
GStr outfname_ns_multi_map;
GStr outfname_ns_uniq_map;
GStr outfname_s_multi_map;
GStr outfname_s_uniq_map;
GStr outfname_s_multi_map_tmp;
// Unpaired
GStr outfname_ns_multi_unpair;
GStr outfname_ns_uniq_unpair;
GStr outfname_s_multi_unpair;
GStr outfname_s_uniq_unpair;
GStr outfname_s_multi_unpair_tmp;

GStr outfname_discard_s_uniq_map;
GStr outfname_discard_s_multi_map;

// GSamWriter 
GSamWriter* outfile_cleaned = NULL;
GSamWriter* outfile_cleaned_2stage = NULL;
// Paired
GSamWriter* outfile_ns_multi_map = NULL;
GSamWriter* outfile_ns_uniq_map = NULL;
GSamWriter* outfile_s_multi_map = NULL;
GSamWriter* outfile_s_uniq_map = NULL;
GSamWriter* outfile_s_multi_map_tmp = NULL;
// Unpaired
GSamWriter* outfile_ns_multi_unpair = NULL;
GSamWriter* outfile_ns_uniq_unpair = NULL;
GSamWriter* outfile_s_multi_unpair = NULL;
GSamWriter* outfile_s_uniq_unpair = NULL;
GSamWriter* outfile_s_multi_unpair_tmp = NULL;

GSamWriter* outfile_discard_s_uniq_map = NULL;
GSamWriter* outfile_discard_s_multi_map = NULL;

FILE* joutf=NULL;

// ALN summary
int ALN_COUNT = 0;
int ALN_COUNT_BAD = 0;
int ALN_COUNT_GOOD = 0;
int ALN_COUNT_GOOD_CAL = 0;

// JUNC summary
int JUNC_COUNT = 0;
int JUNC_COUNT_GOOD = 0;
int JUNC_COUNT_BAD = 0;

// Paired
int ALN_COUNT_SPLICED = 0;
int ALN_COUNT_NSPLICED = 0;
int ALN_COUNT_PAIRED = 0;
int ALN_COUNT_SPLICED_UNIQ = 0;
int ALN_COUNT_SPLICED_MULTI = 0;
int ALN_COUNT_SPLICED_UNIQ_DISCARD = 0;
int ALN_COUNT_SPLICED_MULTI_DISCARD = 0;
int ALN_COUNT_NSPLICED_UNIQ = 0;
int ALN_COUNT_NSPLICED_MULTI = 0;


// Unpaired
int ALN_COUNT_SPLICED_UNPAIR = 0;
int ALN_COUNT_NSPLICED_UNPAIR = 0;
int ALN_COUNT_UNPAIR = 0;
int ALN_COUNT_SPLICED_UNIQ_UNPAIR = 0;
int ALN_COUNT_SPLICED_MULTI_UNPAIR = 0;
int ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD = 0;
int ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD = 0;
int ALN_COUNT_NSPLICED_UNIQ_UNPAIR = 0;
int ALN_COUNT_NSPLICED_MULTI_UNPAIR = 0;

robin_hood::unordered_map<std::string, int>  CHRS;
// robin_hood::unordered_map<std::string, GSamRecordList> read_hashmap;

int STEP_COUNTER = 0;

// j-extract parameters:
int g_max_splice = 100000;
int g_bundle_gap = 100000;
GSamWriter* outfile_above_spliced = NULL;
GSamWriter* outfile_below_spliced = NULL;
FILE* joutf_above=NULL;
FILE* joutf_below=NULL;

// predict parameters:
bool write_bam = true;

// clean parameters:
bool g_paired_removal = false;
bool g_2_stage_run = false;

int main(int argc, char* argv[]) {
    GMessage(
            "==========================================================================================\n"
            "An accurate spliced alignment pruner and spliced junction predictor.\n"
            "==========================================================================================\n");
    const char *banner = R"""(
  ███████╗██████╗ ██╗      █████╗ ███╗   ███╗
  ██╔════╝██╔══██╗██║     ██╔══██╗████╗ ████║
  ███████╗██████╔╝██║     ███████║██╔████╔██║
  ╚════██║██╔═══╝ ██║     ██╔══██║██║╚██╔╝██║
  ███████║██║     ███████╗██║  ██║██║ ╚═╝ ██║
  ╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝
    )""";
    std::cout << banner << std::endl;
    
    in_records.setup(VERSION, argc, argv);
    processOptions(argc, argv);

    outfname_cleaned = out_dir + "/cleaned.bam";
    /*********************
     * For paired uniq- / multi- mapped alignments
    *********************/
    outfname_ns_multi_map = out_dir + "/tmp/nonsplice_multi_map.bam";
    outfname_ns_uniq_map = out_dir + "/tmp/nonsplice_uniq_map.bam";
    outfname_s_multi_map = out_dir + "/tmp/splice_multi_map.bam";
    outfname_s_uniq_map = out_dir + "/tmp/splice_uniq_map.bam";
    outfname_s_multi_map_tmp = out_dir + "/tmp/splice_multi_map_tmp.bam";
        
    /*********************
     * For unpaired uniq- / multi- mapped alignments
    *********************/
    outfname_ns_multi_unpair = out_dir + "/tmp/nonsplice_multi_unpair.bam";
    outfname_ns_multi_unpair = out_dir + "/tmp/nonsplice_uniq_unpair.bam";
    outfname_s_multi_unpair = out_dir + "/tmp/splice_multi_unpair.bam";
    outfname_s_uniq_unpair = out_dir + "/tmp/splice_uniq_unpair.bam";
    outfname_s_multi_unpair_tmp = out_dir + "/tmp/splice_multi_unpair_tmp.bam";

    outfname_discard_s_uniq_map = out_dir + "/discard/discard_splice_uniq_map.bam";
    outfname_discard_s_multi_map = out_dir + "/discard/discard_splice_multi_map.bam";

    // The junction score file
    infname_juncbed = out_dir + "/junction_score.bed";

    GStr discard_dir(out_dir + "/discard");
    std::filesystem::create_directories(out_dir.chars());
    create_CHRS();
    
    // Start reading the files
    int num_samples=in_records.start();

    /*********************
     * Directory creating + Sam writer creation
    *********************/
    // Creating the directory
    std::filesystem::create_directories(discard_dir.chars());
    // The tmp files only for CLEANING
    outfile_s_multi_map_tmp = new GSamWriter(outfname_s_multi_map_tmp, in_records.header(), GSamFile_BAM);
    if (g_paired_removal) {
        outfile_s_multi_unpair_tmp = new GSamWriter(outfname_s_multi_unpair_tmp, in_records.header(), GSamFile_BAM);    
    }

    // cleaned BAM file
    outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);

    // discarded files
    outfile_discard_s_uniq_map = new GSamWriter(outfname_discard_s_uniq_map, in_records.header(), GSamFile_BAM);
    outfile_discard_s_multi_map = new GSamWriter(outfname_discard_s_multi_map, in_records.header(), GSamFile_BAM);   

    GMessage(">> Finish creating GSamWriter\n");

    /*********************
     * Main algorithms
    *********************/
    infname_scorebed = out_dir + "/junction_score.bed";
    infname_NH_tag = splamClean();
    splamNHUpdate();

    /*********************
     * Final statistics printing
    *********************/
    if (COMMAND_MODE == CLEAN) {
        ALN_COUNT_SPLICED = ALN_COUNT_SPLICED_UNIQ + ALN_COUNT_SPLICED_MULTI;
        ALN_COUNT_NSPLICED = ALN_COUNT_NSPLICED_UNIQ + ALN_COUNT_NSPLICED_MULTI;
        ALN_COUNT_SPLICED_UNPAIR = ALN_COUNT_SPLICED_UNIQ_UNPAIR + ALN_COUNT_SPLICED_MULTI_UNPAIR;
        ALN_COUNT_NSPLICED_UNPAIR = ALN_COUNT_NSPLICED_UNIQ_UNPAIR + ALN_COUNT_NSPLICED_MULTI_UNPAIR;

        ALN_COUNT_UNPAIR = ALN_COUNT_SPLICED_UNPAIR + ALN_COUNT_NSPLICED_UNPAIR;

        ALN_COUNT_BAD = ALN_COUNT_SPLICED_UNIQ_DISCARD + ALN_COUNT_SPLICED_MULTI_DISCARD + ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD + ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD;

        ALN_COUNT_GOOD_CAL = ALN_COUNT - ALN_COUNT_BAD;
        if (g_paired_removal) {
            GMessage("\n[INFO] Total number of alignments\t:%10d \n", ALN_COUNT);

            // Printing for paired alignment
            GMessage("           paired alignments\t\t:%10d \n", ALN_COUNT - ALN_COUNT_UNPAIR);
            GMessage("               spliced alignments\t:%10d \n", ALN_COUNT_SPLICED);
            GMessage("                   - uniquely mapped\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_UNIQ, ALN_COUNT_SPLICED_UNIQ-ALN_COUNT_SPLICED_UNIQ_DISCARD, ALN_COUNT_SPLICED_UNIQ_DISCARD);
            GMessage("                   - multi-mapped\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_MULTI, ALN_COUNT_SPLICED_MULTI-ALN_COUNT_SPLICED_MULTI_DISCARD, ALN_COUNT_SPLICED_MULTI_DISCARD);
            GMessage("               non-spliced alignments\t:%10d \n", ALN_COUNT_NSPLICED);
            GMessage("                   - uniquely mapped\t:%10d\n", ALN_COUNT_NSPLICED_UNIQ);
            GMessage("                   - multi-mapped\t:%10d\n\n", ALN_COUNT_NSPLICED_MULTI);

            // Printing for unpaired alignment
            GMessage("           unpaired alignments\t\t:%10d \n", ALN_COUNT_UNPAIR);
            GMessage("               spliced alignments\t:%10d \n", ALN_COUNT_SPLICED_UNPAIR);
            GMessage("                   - uniquely mapped\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_UNIQ_UNPAIR, ALN_COUNT_SPLICED_UNIQ_UNPAIR-ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD, ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD);
            GMessage("                   - multi-mapped\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_MULTI_UNPAIR, ALN_COUNT_SPLICED_MULTI_UNPAIR-ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD, ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD);
            GMessage("               non-spliced alignments\t:%10d \n", ALN_COUNT_NSPLICED_UNPAIR);
            GMessage("                   - uniquely mapped\t:%10d\n", ALN_COUNT_NSPLICED_UNIQ_UNPAIR);
            GMessage("                   - multi-mapped\t:%10d\n\n", ALN_COUNT_NSPLICED_MULTI_UNPAIR);



            GMessage("\n[INFO] Number of junctions\t\t:%10d   (good: %d / bad: %d / unstranded: %d)\n", JUNC_COUNT, JUNC_COUNT_GOOD, JUNC_COUNT_BAD, JUNC_COUNT-JUNC_COUNT_GOOD-JUNC_COUNT_BAD);
            GMessage("\n[INFO] Number of removed alignments\t:%10d \n", ALN_COUNT_BAD);
            GMessage("[INFO] Number of kept alignments\t:%10d \n", ALN_COUNT_GOOD);

            if (ALN_COUNT_GOOD_CAL != ALN_COUNT_GOOD) GMessage("Num of cleaned alignments do not agree with each other. Calculated: %d; iter: %d\n", ALN_COUNT_GOOD_CAL, ALN_COUNT_GOOD);
        } else {
            GMessage("\n[INFO] Total number of alignments\t:%10d \n", ALN_COUNT);
            GMessage("           spliced alignments\t\t:%10d \n", ALN_COUNT_SPLICED);
            GMessage("               - uniquely mapped\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_UNIQ, ALN_COUNT_SPLICED_UNIQ-ALN_COUNT_SPLICED_UNIQ_DISCARD, ALN_COUNT_SPLICED_UNIQ_DISCARD);
            GMessage("               - multi-mapped\t\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_MULTI, ALN_COUNT_SPLICED_MULTI-ALN_COUNT_SPLICED_MULTI_DISCARD, ALN_COUNT_SPLICED_MULTI_DISCARD);
            
            GMessage("           non-spliced alignments\t:%10d \n", ALN_COUNT_NSPLICED);
            GMessage("               - uniquely mapped\t:%10d\n", ALN_COUNT_NSPLICED_UNIQ);
            GMessage("               - multi-mapped\t\t:%10d\n", ALN_COUNT_NSPLICED_MULTI);

            GMessage("\n[INFO] Number of junctions\t\t:%10d   (good: %d / bad: %d / unstranded: %d)\n", JUNC_COUNT, JUNC_COUNT_GOOD, JUNC_COUNT_BAD, JUNC_COUNT-JUNC_COUNT_GOOD-JUNC_COUNT_BAD);
            GMessage("\n[INFO] Number of removed alignments\t:%10d \n", ALN_COUNT_BAD);
            GMessage("[INFO] Number of kept alignments\t:%10d \n", ALN_COUNT_GOOD);

            if (ALN_COUNT_GOOD_CAL != ALN_COUNT_GOOD) GMessage("Num of cleaned alignments do not agree with each other. Calculated: %d; iter: %d\n", ALN_COUNT_GOOD_CAL, ALN_COUNT_GOOD);
        }
    }
    return 0;
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;cite;verbose;version;paired-removal;junction;no-write-bam;2-stage-run;model=;output=;score=;max-splice=;bundle-gap=;hvcVSPJo:N:Q:m:r:s:M:g:");
    // args.printError(usage_clean, true);
    command_str=args.nextNonOpt();
    if (argc == 0) {
        usage_clean();
        GERROR("\n[ERROR] No command provide. The subcommand must be 'j-extract', 'predict', or 'clean'.\n");
        exit(1);   
    }

    if (args.getOpt('h') || args.getOpt("help")) {
        usage_clean();
        exit(0);
    }

    if (args.getOpt('v') || args.getOpt("version")) {
        fprintf(stdout,"SPLAM v.%s\n", VERSION);
        exit(0);
    }

    if (args.getOpt('c') || args.getOpt("cite")) {
        fprintf(stdout,"%s\n", "Kuan-Hao Chao, Mihaela Pertea, and Steven Salzberg, \033[1m\x1B[3mSPLAM: accurate deep-learning-based splice site predictor to clean up spurious spliced alignments\x1B[0m\033[0m, (2023), GitHub repository, https://github.com/Kuanhao-Chao/SPLAM");
        exit(0);
    }

    if (args.getOpt('P') || args.getOpt("paired-removal") ) {
        g_paired_removal = true;
        GMessage("g_paired_removal: %d\n", g_paired_removal);
    }


    if (verbose) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
    }
    
    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (verbose) {
        // fprintf(stderr, "Running SPLAM " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }

    /********************************
     * Process the arguments
    *********************************/
    processOptionsClean(args);

// #ifdef DEBUG
    GMessage(">>  command_str      : %s\n", command_str.chars());
    GMessage(">> infname_model_name: %s\n", infname_model_name.chars());
    GMessage(">> infname_juncbed   : %s\n", infname_juncbed.chars());
    GMessage(">> infname_scorebed  : %s\n", infname_scorebed.chars());
    GMessage(">> infname_reffa     : %s\n", infname_reffa.chars());
    GMessage(">> infname_bam       : %s\n", infname_bam.chars());
    GMessage(">> out_dir           : %s\n", out_dir.chars());
    GMessage(">> g_paired_removal  : %d\n", g_paired_removal);
    GMessage(">> g_max_splice      : %d\n", g_max_splice);
    GMessage(">> g_bundle_gap      : %d\n", g_bundle_gap);
    GMessage(">> verbose           : %d\n", verbose);
// #endif

    args.startNonOpt();
}

void processOptionsClean(GArgs& args) {
    optionsOutput(args);
}


void optionsOutput(GArgs& args) {
    // -o / --output
    out_dir=args.getOpt('o');    
    if (out_dir.is_empty()) {
        out_dir=args.getOpt("output");
        if (out_dir.is_empty()) {
            usage_clean();
            GMessage("\n[ERROR] output directory must be provided (-o / --output)!\n");
            exit(1);
        }
    }
}


void optionsWriteTMP(GArgs& args) {
    //--no-write-bam
    write_bam = (args.getOpt("no-write-bam")==NULL);
    GMessage("write_bam: %d\n", write_bam);
    if (write_bam) {
        GMessage(">>  Running write_bam mode\n");
    } else {
        GMessage(">>  Running predict mode\n");
    }
}


// if (fileExists(infname_juncbed.chars())>1) {
//     // guided=true;
// } else {
//     GError("[ERROR] junction bed file (%s) not found.\n",
//         infname_juncbed.chars());
// }