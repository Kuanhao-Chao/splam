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
// #include "filter.h"


// #include "splam_stream.h"


#include <gclib/GArgs.h>
#include <gclib/GBase.h>
#include <gclib/GStr.h>
#include <robin_hood/robin_hood.h>

#include <Python.h>


void processOptions(int argc, char* argv[]);
void processOptionsJExtract(GArgs& args);
void processOptionsPredict(GArgs& args);
void processOptionsClean(GArgs& args);
void processOptionsNHUpdate(GArgs& args);

void optionsMaxSplice(GArgs& args);
void optionsBundleGap(GArgs& args);
void optionsModel(GArgs& args);
void optionsRef(GArgs& args);
void optionsOutput(GArgs& args);
void checkJunction(GArgs& args);
void optionsJunction(GArgs& args);

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
bool predict_junc_mode = false;
TInputFiles in_records;
TInputRecord* irec=NULL;

float threshold = 0.3;
int aln_num_thr = 4;

GSamRecord* brec=NULL;

// output file names 
GStr outfname_cleaned;
// Paired
GStr outfname_ns_multi_map;
GStr outfname_s_uniq_map;
GStr outfname_s_multi_map;
GStr outfname_s_multi_map_tmp;
// Unpaired
GStr outfname_ns_multi_unpair;
GStr outfname_s_uniq_unpair;
GStr outfname_s_multi_unpair;
GStr outfname_s_multi_unpair_tmp;

GStr outfname_discard_s_uniq_map;
GStr outfname_discard_s_multi_map;

// GSamWriter 
GSamWriter* outfile_cleaned = NULL;
// Paired
GSamWriter* outfile_ns_multi_map = NULL;
GSamWriter* outfile_s_uniq_map = NULL;
GSamWriter* outfile_s_multi_map = NULL;
GSamWriter* outfile_s_multi_map_tmp = NULL;
// Unpaired
GSamWriter* outfile_ns_multi_unpair = NULL;
GSamWriter* outfile_s_uniq_unpair = NULL;
GSamWriter* outfile_s_multi_unpair = NULL;
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

// clean parameters:
bool g_paired_removal = false;

bool g_is_single_end = false;

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
    /*********************
     * For paired uniq- / multi- mapped alignments
    *********************/
    outfname_ns_multi_map = out_dir + "/TMP/nonsplice_multi_map.bam";
    outfname_s_uniq_map = out_dir + "/TMP/splice_uniq_map.bam";
    outfname_s_multi_map = out_dir + "/TMP/splice_multi_map.bam";
    outfname_s_multi_map_tmp = out_dir + "/TMP/splice_multi_map_tmp.bam";
        
    // outfname_unpair = out_dir + "/TMP/unpair.bam";
    /*********************
     * For unpaired uniq- / multi- mapped alignments
    *********************/
    outfname_ns_multi_unpair = out_dir + "/TMP/nonsplice_multi_unpair.bam";
    outfname_s_uniq_unpair = out_dir + "/TMP/splice_uniq_unpair.bam";
    outfname_s_multi_unpair = out_dir + "/TMP/splice_multi_unpair.bam";
    outfname_s_multi_unpair_tmp = out_dir + "/TMP/splice_multi_unpair_tmp.bam";

    outfname_discard_s_uniq_map = out_dir + "/discard/discard_splice_uniq_map.bam";;
    outfname_discard_s_multi_map = out_dir + "/discard/discard_splice_multi_map.bam";;

    GStr tmp_dir(out_dir + "/TMP");
    GStr discard_dir(out_dir + "/discard");
    std::filesystem::create_directories(out_dir.chars());
    create_CHRS();
    
    // Start reading the files
    int num_samples=in_records.start();
    if (COMMAND_MODE == CLEAN) {
        std::filesystem::create_directories(tmp_dir.chars());
        std::filesystem::create_directories(discard_dir.chars());
        // The only two modes that need to write out files.
        outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);

        // Paired
        outfile_ns_multi_map = new GSamWriter(outfname_ns_multi_map, in_records.header(), GSamFile_BAM);
        outfile_s_uniq_map = new GSamWriter(outfname_s_uniq_map, in_records.header(), GSamFile_BAM);
        outfile_s_multi_map = new GSamWriter(outfname_s_multi_map, in_records.header(), GSamFile_BAM);
        outfile_s_multi_map_tmp = new GSamWriter(outfname_s_multi_map_tmp, in_records.header(), GSamFile_BAM);
        // outfile_unpair = new GSamWriter(outfname_unpair, in_records.header(), GSamFile_BAM);

        // Unpaired
        outfile_ns_multi_unpair = new GSamWriter(outfname_ns_multi_unpair, in_records.header(), GSamFile_BAM);
        outfile_s_uniq_unpair = new GSamWriter(outfname_s_uniq_unpair, in_records.header(), GSamFile_BAM);
        outfile_s_multi_unpair = new GSamWriter(outfname_s_multi_unpair, in_records.header(), GSamFile_BAM);
        outfile_s_multi_unpair_tmp = new GSamWriter(outfname_s_multi_unpair_tmp, in_records.header(), GSamFile_BAM);
        
        outfile_discard_s_uniq_map= new GSamWriter(outfname_discard_s_uniq_map, in_records.header(), GSamFile_BAM);
        outfile_discard_s_multi_map= new GSamWriter(outfname_discard_s_multi_map, in_records.header(), GSamFile_BAM);
    }
    

    if (COMMAND_MODE == J_EXTRACT) {
        infname_juncbed = splamJExtract();
        GMessage("\n[INFO] Total number of junctions\t:%d\n", JUNC_COUNT);
    } else if (COMMAND_MODE == PREDICT) {
        if (!predict_junc_mode) {
            infname_juncbed = splamJExtract();
        }
        infname_scorebed = splamPredict();
    } else if (COMMAND_MODE == CLEAN) {
        infname_juncbed = splamJExtract();
        infname_scorebed = splamPredict();
        infname_NH_tag = splamClean();
        splamNHUpdate();
    }

    if (COMMAND_MODE == CLEAN) {

        ALN_COUNT_SPLICED = ALN_COUNT_SPLICED_UNIQ + ALN_COUNT_SPLICED_MULTI;
        ALN_COUNT_NSPLICED = ALN_COUNT_NSPLICED_UNIQ + ALN_COUNT_NSPLICED_MULTI;
        if (g_paired_removal) {

            ALN_COUNT_SPLICED_UNPAIR = ALN_COUNT_SPLICED_UNIQ_UNPAIR + ALN_COUNT_SPLICED_MULTI_UNPAIR;
            ALN_COUNT_NSPLICED_UNPAIR = ALN_COUNT_NSPLICED_UNIQ_UNPAIR + ALN_COUNT_NSPLICED_MULTI_UNPAIR;
            ALN_COUNT_UNPAIR = ALN_COUNT_SPLICED_UNIQ_UNPAIR + ALN_COUNT_SPLICED_MULTI_UNPAIR + ALN_COUNT_NSPLICED_UNIQ_UNPAIR + ALN_COUNT_NSPLICED_MULTI_UNPAIR;

            ALN_COUNT_BAD = ALN_COUNT_SPLICED_UNIQ_DISCARD + ALN_COUNT_SPLICED_MULTI_DISCARD;
            ALN_COUNT_GOOD_CAL = ALN_COUNT - ALN_COUNT_BAD;

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
            ALN_COUNT_BAD = ALN_COUNT_SPLICED_UNIQ_DISCARD + ALN_COUNT_SPLICED_MULTI_DISCARD;
            ALN_COUNT_GOOD_CAL = ALN_COUNT - ALN_COUNT_BAD;
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
    GArgs args(argc, argv, "help;cite;verbose;version;single-end;paired-removal;junction;model=;output=;score=;max-splice=;bundle-gap=;hvcVSPJo:N:Q:m:r:s:M:g:");
    // args.printError(USAGE, true);
    command_str=args.nextNonOpt();
    if (argc == 0) {
        usage();
        GERROR("\n[ERROR] No command provide. The subcommand must be 'j-extract', 'predict', or 'clean'.\n");
        exit(1);   
    }

    if (args.getOpt('h') || args.getOpt("help")) {
        usage();
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

    if (args.getOpt('S') || args.getOpt("single-end") ) {
        g_is_single_end = true;
        GMessage("g_is_single_end: %d\n", g_is_single_end);
    }

    if (args.getOpt('P') || args.getOpt("paired-removal") ) {
        g_paired_removal = true;
        GMessage("g_paired_removal: %d\n", g_paired_removal);
    }

    if (strcmp(command_str.chars(), "j-extract") == 0) {
        COMMAND_MODE = J_EXTRACT;
    } else if (strcmp(command_str.chars(), "predict") == 0) {
        COMMAND_MODE = PREDICT;
    } else if (strcmp(command_str.chars(), "clean") == 0) {
        COMMAND_MODE = CLEAN;
    } else {
        usage();
        GERROR("\n[ERROR] The subcommand must be 'j-extract', 'predict', or 'clean'.\n");
        exit(1);   
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
     * Process arguments by COMMAND_MODE
    *********************************/
    if (COMMAND_MODE == J_EXTRACT) {
        processOptionsJExtract(args);    
    } else if (COMMAND_MODE == PREDICT) {
        GMessage(">> Inside PREDICT:\n");
        processOptionsPredict(args);
    } else if (COMMAND_MODE == CLEAN) {
        processOptionsClean(args);
    }

    GMessage(">>  command_str      : %s\n", command_str.chars());
    GMessage(">> infname_model_name: %s\n", infname_model_name.chars());
    GMessage(">> infname_juncbed   : %s\n", infname_juncbed.chars());
    GMessage(">> infname_scorebed  : %s\n", infname_scorebed.chars());
    GMessage(">> infname_reffa     : %s\n", infname_reffa.chars());
    GMessage(">> infname_bam       : %s\n", infname_bam.chars());
    GMessage(">> out_dir           : %s\n", out_dir.chars());
    GMessage(">> g_is_single_end   : %d\n", g_is_single_end);
    GMessage(">> g_paired_removal  : %d\n", g_paired_removal);
    GMessage(">> g_max_splice      : %d\n", g_max_splice);
    GMessage(">> g_bundle_gap      : %d\n", g_bundle_gap);
    GMessage(">> verbose           : %d\n", verbose);

    args.startNonOpt();

    if (args.getNonOptCount()==1) {
        usage();
        GMessage("\n[ERROR] no input provided!\n");
        exit(1);
    }
    args.nextNonOpt(); 

    if (!predict_junc_mode) {
        const char* ifn=NULL;
        while ( (ifn=args.nextNonOpt())!=NULL) {
            //input alignment files
            std::string absolute_ifn = get_full_path(ifn);
            in_records.addFile(absolute_ifn.c_str());
        }
    } else if (predict_junc_mode && COMMAND_MODE == PREDICT) {
        checkJunction(args);
    }

    // } 
    // else if (COMMAND_MODE == PREDICT) {
    //     const char* ifn=NULL;
    //     while ( (ifn=args.nextNonOpt())!=NULL) {
    //         //input alignment files
    //         std::string absolute_ifn = get_full_path(ifn);
    //         std::cout << "absolute_ifn: " << absolute_ifn << std::endl;
    //         in_records.addFile(absolute_ifn.c_str());
    //     }
    // } else if (COMMAND_MODE == CLEAN) {
    // }

    // int num_samples=in_records.start();
    // outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);
    // outfile_ns_multi_map = new GSamWriter(outfname_ns_multi_map, in_records.header(), GSamFile_BAM);
    // outfile_s_uniq_map = new GSamWriter(outfname_s_uniq_map, in_records.header(), GSamFile_BAM);
    // outfile_s_multi_map = new GSamWriter(outfname_s_multi_map, in_records.header(), GSamFile_BAM);
    // outfile_unpair = new GSamWriter(outfname_unpair, in_records.header(), GSamFile_BAM);
}


void processOptionsJExtract(GArgs& args) {
    optionsOutput(args);
    optionsMaxSplice(args);
    optionsBundleGap(args);
}


void processOptionsPredict(GArgs& args) {
    optionsModel(args);
    optionsRef(args);
    optionsOutput(args);
    optionsJunction(args);
}


void processOptionsClean(GArgs& args) {
    optionsModel(args);
    optionsRef(args);
    optionsOutput(args);
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


void optionsModel(GArgs& args) {
    // -m / --model
    infname_model_name=args.getOpt('m');
    if (infname_model_name.is_empty()) {
        infname_model_name=args.getOpt("model");
        if (infname_model_name.is_empty()) {
            usage();
            GMessage("\n[ERROR] model file must be provided (-m)!\n");
            exit(1);
        } else {
            if (fileExists(infname_model_name.chars())>1) {
                // guided=true;
            } else {
                GError("[ERROR] model file (%s) not found.\n",
                    infname_model_name.chars());
            }
        }
    }
} 


void optionsRef(GArgs& args) {
    // -r / --ref
    infname_reffa=args.getOpt('r');        
    if (infname_reffa.is_empty()) {
        infname_reffa=args.getOpt("ref");
        if (infname_reffa.is_empty()) {
            usage();
            GMessage("\n[ERROR] reference fasta file must be provided (-r)!\n");
            exit(1);
        } else {
            if (fileExists(infname_reffa.chars())>1) {
                // guided=true;
            } else {
                GError("[ERROR] reference fasta file (%s) not found.\n",
                    infname_reffa.chars());
            }
        }
    }
}


void optionsOutput(GArgs& args) {
    // -o / --output
    out_dir=args.getOpt('o');    
    if (out_dir.is_empty()) {
        out_dir=args.getOpt("output");
        if (out_dir.is_empty()) {
            usage();
            GMessage("\n[ERROR] output directory must be provided (-o / --output)!\n");
            exit(1);
        }
    }
}


void optionsJunction(GArgs& args) {
    // -J / --junction
    predict_junc_mode = (args.getOpt("junction")!=NULL || args.getOpt('J')!=NULL);
    if (predict_junc_mode) {
        GMessage(">>  SPLAM predict [Junction mode]\n");
    } else {
        GMessage(">>  SPLAM predict [Default mode]\n");
    }
}


void checkJunction(GArgs& args) {
    infname_juncbed = args.nextNonOpt(); 
    if (infname_juncbed.is_empty()) {
        usage();
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


// void optionsScore(GArgs& args) {
//     // -s / --score
//     infname_scorebed = args.getOpt('s');
//     if (infname_scorebed.is_empty()) {
//         infname_scorebed=args.getOpt("score");
//         if (infname_scorebed.is_empty()) {
//             usage();
//             GMessage("\n[ERROR] junction score bed file must be provided (-s / --score)!\n");
//             exit(1);
//         } else {
//             if (fileExists(infname_scorebed.chars())>1) {
//                 // guided=true;
//             } else {
//                 GError("[ERROR] junction score bed file (%s) not found.\n",
//                     infname_scorebed.chars());
//             }
//         }
//     } else {
//         if (fileExists(infname_scorebed.chars())>1) {
//             // guided=true;
//         } else {
//             GError("[ERROR] junction score bed file (%s) not found.\n",
//                 infname_scorebed.chars());
//         }
//     }
// }
