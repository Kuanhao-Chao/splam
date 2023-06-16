// #define DEBUG
#include <iostream>
#include <vector>
#include <memory>
#include <filesystem>

#include "extract.h"

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
void processOptionsJExtract(GArgs& args);

void optionsMaxSplice(GArgs& args);
void optionsBundleGap(GArgs& args);
void optionsOutput(GArgs& args);
void checkJunction(GArgs& args);
void optionsJunction(GArgs& args);
void optionsWriteTMP(GArgs& args);
// void options2StageRun(GArgs& args);

CommandMode COMMAND_MODE = UNSET;
GStr command_str;

// input file names 
GStr infname_bam("");
GStr infname_juncbed("");
GStr infname_scorebed("");
GStr infname_NH_tag("");
GStr out_dir;

bool verbose = false;
bool predict_junc_mode = false;
TInputFiles in_records;
TInputRecord* irec=NULL;

float threshold = 0.1;
int aln_num_thr = 4;
GSamRecord* brec=NULL;

// output file names 
GStr outfname_cleaned;
GStr outfname_cleaned_2stage;

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

robin_hood::unordered_map<std::string, int>  CHRS;

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
  ███████╗██████╗ ██╗      █████╗ ███╗   ███╗    ██╗
  ██╔════╝██╔══██╗██║     ██╔══██╗████╗ ████║    ██║
  ███████╗██████╔╝██║     ███████║██╔████╔██║    ██║
  ╚════██║██╔═══╝ ██║     ██╔══██║██║╚██╔╝██║    ╚═╝
  ███████║██║     ███████╗██║  ██║██║ ╚═╝ ██║    ██╗
  ╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝    ╚═╝
    )""";
    std::cout << banner << std::endl;
    
    in_records.setup(VERSION, argc, argv);
    processOptions(argc, argv);

    outfname_cleaned = out_dir + "/cleaned.bam";
    outfname_cleaned_2stage = out_dir + "/cleaned_2stage.bam";
    /*********************
     * For paired uniq- / multi- mapped alignments
    *********************/
    outfname_ns_multi_map = out_dir + "/TMP/nonsplice_multi_map.bam";
    outfname_ns_uniq_map = out_dir + "/TMP/nonsplice_uniq_map.bam";
    outfname_s_multi_map = out_dir + "/TMP/splice_multi_map.bam";
    outfname_s_uniq_map = out_dir + "/TMP/splice_uniq_map.bam";
    outfname_s_multi_map_tmp = out_dir + "/TMP/splice_multi_map_tmp.bam";
        
    /*********************
     * For unpaired uniq- / multi- mapped alignments
    *********************/
    outfname_ns_multi_unpair = out_dir + "/TMP/nonsplice_multi_unpair.bam";
    outfname_ns_multi_unpair = out_dir + "/TMP/nonsplice_uniq_unpair.bam";
    outfname_s_multi_unpair = out_dir + "/TMP/splice_multi_unpair.bam";
    outfname_s_uniq_unpair = out_dir + "/TMP/splice_uniq_unpair.bam";
    outfname_s_multi_unpair_tmp = out_dir + "/TMP/splice_multi_unpair_tmp.bam";

    outfname_discard_s_uniq_map = out_dir + "/discard/discard_splice_uniq_map.bam";;
    outfname_discard_s_multi_map = out_dir + "/discard/discard_splice_multi_map.bam";;

    GStr tmp_dir(out_dir + "/TMP");
    GStr discard_dir(out_dir + "/discard");
    std::filesystem::create_directories(out_dir.chars());
    create_CHRS();
    
    // Start reading the files
    int num_samples=in_records.start();

    /*********************
     * Directory creating + Sam writer creation
    *********************/
   
    if (write_bam) {
        // Creating the directory
        std::filesystem::create_directories(tmp_dir.chars());

        // Paired
        outfile_ns_multi_map = new GSamWriter(outfname_ns_multi_map, in_records.header(), GSamFile_BAM);
        outfile_ns_uniq_map = new GSamWriter(outfname_ns_uniq_map, in_records.header(), GSamFile_BAM);
        outfile_s_uniq_map = new GSamWriter(outfname_s_uniq_map, in_records.header(), GSamFile_BAM);
        outfile_s_multi_map = new GSamWriter(outfname_s_multi_map, in_records.header(), GSamFile_BAM);

        if (g_paired_removal) {
            // Unpaired
            outfile_ns_multi_unpair = new GSamWriter(outfname_ns_multi_unpair, in_records.header(), GSamFile_BAM);
            outfile_ns_uniq_unpair = new GSamWriter(outfname_ns_uniq_unpair, in_records.header(), GSamFile_BAM);
            outfile_s_multi_unpair = new GSamWriter(outfname_s_multi_unpair, in_records.header(), GSamFile_BAM);
            outfile_s_uniq_unpair = new GSamWriter(outfname_s_uniq_unpair, in_records.header(), GSamFile_BAM);
        }
    }

    GMessage(">> Finish creating GSamWriter\n");

    /*********************
     * Main algorithms
    *********************/
    // infname_juncbed = splamJExtract();
    GMessage("\n[INFO] Total number of junctions\t:%d\n", JUNC_COUNT);
    return 0;
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;cite;verbose;version;paired-removal;junction;no-write-bam;model=;output=;score=;max-splice=;bundle-gap=;hvcVSPJo:N:Q:m:r:s:M:g:");

    command_str=args.nextNonOpt();
    if (argc == 0) {
        usage_extract();
        GERROR("\n[ERROR] No command provide. The subcommand must be 'j-extract', 'predict', or 'clean'.\n");
        exit(1);   
    }

    if (args.getOpt('h') || args.getOpt("help")) {
        usage_extract();
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

    // if (verbose) {
    //     GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
    // }
    
    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (verbose) {
        // fprintf(stderr, "Running SPLAM " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }

    /********************************
     * Process arguments by COMMAND_MODE
    *********************************/
    processOptionsJExtract(args);    

// #ifdef DEBUG
    GMessage(">>  command_str      : %s\n", command_str.chars());
    GMessage(">> infname_juncbed   : %s\n", infname_juncbed.chars());
    GMessage(">> infname_scorebed  : %s\n", infname_scorebed.chars());
    GMessage(">> infname_bam       : %s\n", infname_bam.chars());
    GMessage(">> out_dir           : %s\n", out_dir.chars());
    GMessage(">> g_paired_removal  : %d\n", g_paired_removal);
    GMessage(">> g_max_splice      : %d\n", g_max_splice);
    GMessage(">> g_bundle_gap      : %d\n", g_bundle_gap);
    GMessage(">> verbose           : %d\n", verbose);
// #endif

    args.startNonOpt();

    if (args.getNonOptCount()==1) {
        usage_extract();
        GMessage("\n[ERROR] no input provided!\n");
        exit(1);
    }
    args.nextNonOpt(); 

    const char* ifn=NULL;
    while ( (ifn=args.nextNonOpt())!=NULL) {
        //input alignment files
        std::string absolute_ifn = get_full_path(ifn);
        in_records.addFile(absolute_ifn.c_str());
    }
    checkJunction(args);
}


void processOptionsJExtract(GArgs& args) {
    optionsOutput(args);
    optionsMaxSplice(args);
    optionsBundleGap(args);
    optionsWriteTMP(args);
}


void processOptionsPredict(GArgs& args) {
    optionsOutput(args);
    optionsJunction(args);
    optionsWriteTMP(args);
}


void processOptionsClean(GArgs& args) {
    optionsOutput(args);
    // options2StageRun(args);
}


void optionsMaxSplice(GArgs& args) {
    // -M / --max-splice
    GStr s;
    s = args.getOpt('M');
    if (!s.is_empty()) {
        g_max_splice = s.asInt();
    } else {
        s=args.getOpt("max-splice");
        if (!s.is_empty()) {
            // Use the default max-splice
            g_max_splice = s.asInt();
        }
    }
}


void optionsBundleGap(GArgs& args) {
    // -g / --bundle-gap
    GStr s;
    s = args.getOpt('g');
    if (!s.is_empty()) {
        g_bundle_gap = s.asInt();
    } else {
        s=args.getOpt("bundle-gap");
        if (!s.is_empty()) {
            // Use the default bundle-gap
            g_bundle_gap = s.asInt();
        }
    }
}


void optionsOutput(GArgs& args) {
    // -o / --output
    out_dir=args.getOpt('o');    
    if (out_dir.is_empty()) {
        out_dir=args.getOpt("output");
        if (out_dir.is_empty()) {
            usage_extract();
            GMessage("\n[ERROR] output directory must be provided (-o / --output)!\n");
            exit(1);
        }
    }
}


void optionsJunction(GArgs& args) {
    // -J / --junction
    predict_junc_mode = (args.getOpt("junction")!=NULL || args.getOpt('J')!=NULL);
    // if (predict_junc_mode) {
    //     GMessage(">>  SPLAM predict [Junction mode]\n");
    // } else {
    //     GMessage(">>  SPLAM predict [Default mode]\n");
    // }
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

void checkJunction(GArgs& args) {
    infname_juncbed = args.nextNonOpt(); 
    if (infname_juncbed.is_empty()) {
        usage_extract();
        GMessage("\n[ERROR] junction input bed file must be provided!\n");
        exit(1);
    } else {
        if (fileExists(infname_juncbed.chars())>1) {
            // guided=true;
        } else {
            GError("[ERROR] junction bed file (%s) not found.\n",
                infname_juncbed.chars());
        }
    }
}