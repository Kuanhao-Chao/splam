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
void processOptionsAll(GArgs& args);
void processOptionsNHUpdate(GArgs& args);


void optionsModel(GArgs& args);
void optionsRef(GArgs& args);
void optionsOutput(GArgs& args);
void optionsJunction(GArgs& args);
void optionsScore(GArgs& args);

CommandMode COMMAND_MODE = UNSET;
GStr command_str;

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

float threshold = 0.3;
int aln_num_thr = 4;

GStr outfname_multimapped;
GStr outfname_spliced;
GStr outfname_discard;
GStr outfname_cleaned;
GStr outfname_cleaned_tmp;

GSamRecord* brec=NULL;
GSamWriter* outfile_multimapped = NULL;
GSamWriter* outfile_spliced = NULL;
GSamWriter* outfile_discard = NULL;
GSamWriter* outfile_cleaned = NULL;
GSamWriter* outfile_cleaned_tmp = NULL;
FILE* joutf=NULL;

int JUNC_COUNT = 0;
int ALN_COUNT = 0;
int ALN_COUNT_SPLICED = 0;
int ALN_COUNT_NSPLICED = 0;
int ALN_COUNT_BAD = 0;
int ALN_COUNT_GOOD = 0;
int ALN_COUNT_NH_UPDATE = 0;

robin_hood::unordered_map<std::string, int>  CHRS;
// robin_hood::unordered_map<std::string, GSamRecordList> read_hashmap;

int STEP_COUNTER = 0;
// j-extract parameters:
int j_extract_threshold = 0;
GSamWriter* outfile_above_spliced = NULL;
GSamWriter* outfile_below_spliced = NULL;
FILE* joutf_above=NULL;
FILE* joutf_below=NULL;

// robin_hood::unordered_set<std::string>* rm_rd_set;

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


    outfname_cleaned_tmp = out_dir + "/TMP/cleaned_tmp.bam";
    outfname_multimapped = out_dir + "/TMP/multimapped.bam";
    outfname_spliced = out_dir + "/TMP/spliced.bam";
    outfname_cleaned = out_dir + "/cleaned.bam";
    outfname_discard = out_dir + "/discard.bam";
    
    GStr tmp_dir(out_dir + "/TMP");
    std::filesystem::create_directories(out_dir.chars());
    std::filesystem::create_directories(tmp_dir.chars());

    create_CHRS();

    if (COMMAND_MODE == J_EXTRACT) {
        infname_juncbed = splamJExtract();
    } else if (COMMAND_MODE == PREDICT) {
        infname_scorebed = splamPredict();
    } else if (COMMAND_MODE == CLEAN) {
        infname_NH_tag = splamClean();
        splamNHUpdate();
        
        GMessage("\n\n[INFO] Total number of alignments\t:\t%d\n", ALN_COUNT);
        GMessage("[INFO]     spliced alignments\t\t:\t%d\n", ALN_COUNT_SPLICED);
        GMessage("[INFO]     non-spliced alignments\t:\t%d\n", ALN_COUNT_NSPLICED);
        GMessage("[INFO] Number of removed alignments\t:\t%d\n", ALN_COUNT_BAD);
        GMessage("[INFO] Number of kept alignments\t:\t%d\n", ALN_COUNT_GOOD);
        
    } else if (COMMAND_MODE == ALL) {
        infname_juncbed = splamJExtract();
        infname_scorebed = splamPredict();
        infname_NH_tag = splamClean();
        splamNHUpdate();

        // GStr outfname_cleaned = splamNHUpdate();
        
        GMessage("\n\n[INFO] Total number of alignments\t:\t%d\n", ALN_COUNT);
        GMessage("[INFO]     spliced alignments\t\t:\t%d\n", ALN_COUNT_SPLICED);
        GMessage("[INFO]     non-spliced alignments\t:\t%d\n", ALN_COUNT_NSPLICED);
        GMessage("[INFO] Number of removed alignments\t:\t%d\n", ALN_COUNT_BAD);
        GMessage("[INFO] Number of kept alignments\t:\t%d\n", ALN_COUNT_GOOD);
    }

    return 0;
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;cite;verbose;version;model=;junction=;threshold=;output=;score=;hvcVo:t:N:Q:m:j:r:s:");
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
        fprintf(stdout,"%s\n", "Paper to SPLAM");
        exit(0);
    }

    if (strcmp(command_str.chars(), "j-extract") == 0) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
        COMMAND_MODE = J_EXTRACT;
    } else if (strcmp(command_str.chars(), "predict") == 0) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
        COMMAND_MODE = PREDICT;
    } else if (strcmp(command_str.chars(), "clean") == 0) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
        COMMAND_MODE = CLEAN;
    } else if (strcmp(command_str.chars(), "all") == 0) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
        COMMAND_MODE = ALL;
    } else {
        usage();
        GERROR("\n[ERROR] The subcommand must be 'j-extract', 'predict', or 'clean'.\n");
        exit(1);   
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
    } else if (COMMAND_MODE == ALL) {
        processOptionsAll(args);
    }

    GMessage(">>  command_str      : %s\n", command_str.chars());
    GMessage(">> infname_model_name: %s\n", infname_model_name.chars());
    GMessage(">> infname_juncbed   : %s\n", infname_juncbed.chars());
    GMessage(">> infname_scorebed  : %s\n", infname_scorebed.chars());
    GMessage(">> infname_reffa     : %s\n", infname_reffa.chars());
    GMessage(">> infname_bam       : %s\n", infname_bam.chars());
    GMessage(">> out_dir           : %s\n", out_dir.chars());

    args.startNonOpt();

    if (args.getNonOptCount()==1) {
        usage();
        GMessage("\n[ERROR] no input provided!\n");
        exit(1);
    }
    infname_bam=args.nextNonOpt(); 


    // if (COMMAND_MODE == J_EXTRACT) {
    const char* ifn=NULL;
    while ( (ifn=args.nextNonOpt())!=NULL) {
        //input alignment files
        std::string absolute_ifn = get_full_path(ifn);
        std::cout << "absolute_ifn: " << absolute_ifn << std::endl;
        in_records.addFile(absolute_ifn.c_str());
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
}


void processOptionsJExtract(GArgs& args) {
    optionsOutput(args);

    // -t / --threshold
    GStr s;
    s = args.getOpt('t');
    if (s.is_empty()) {
        s = args.getOpt("threshold");
        if (s.is_empty()) {
            j_extract_threshold = 0;
        } else {
            j_extract_threshold = s.asInt();
        }
    } else {
        j_extract_threshold = s.asInt();       
    }
}


void processOptionsPredict(GArgs& args) {
    optionsModel(args);
    optionsJunction(args);
    optionsRef(args);
    optionsOutput(args);
}


void processOptionsClean(GArgs& args) {
    optionsModel(args);
    optionsRef(args);
    optionsScore(args);
    optionsOutput(args);
}

void processOptionsAll(GArgs& args) {
    optionsModel(args);
    optionsRef(args);
    optionsOutput(args);
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
    // -j / --junction
    infname_juncbed = args.getOpt('j');
    if (infname_juncbed.is_empty()) {
        infname_juncbed=args.getOpt("junction");
        if (infname_juncbed.is_empty()) {
            usage();
            GMessage("\n[ERROR] junction bed file must be provided (-j / --junction)!\n");
            exit(1);
        } else {
            if (fileExists(infname_juncbed.chars())>1) {
                // guided=true;
            } else {
                GError("[ERROR] junction bed file (%s) not found.\n",
                    infname_juncbed.chars());
            }
        }
    } else {
        if (fileExists(infname_juncbed.chars())>1) {
            // guided=true;
        } else {
            GError("[ERROR] junction bed file (%s) not found.\n",
                infname_juncbed.chars());
        }
    }
}


void optionsScore(GArgs& args) {
    // -s / --score
    infname_scorebed = args.getOpt('s');
    if (infname_scorebed.is_empty()) {
        infname_scorebed=args.getOpt("score");
        if (infname_scorebed.is_empty()) {
            usage();
            GMessage("\n[ERROR] junction score bed file must be provided (-s / --score)!\n");
            exit(1);
        } else {
            if (fileExists(infname_scorebed.chars())>1) {
                // guided=true;
            } else {
                GError("[ERROR] junction score bed file (%s) not found.\n",
                    infname_scorebed.chars());
            }
        }
    } else {
        if (fileExists(infname_scorebed.chars())>1) {
            // guided=true;
        } else {
            GError("[ERROR] junction score bed file (%s) not found.\n",
                infname_scorebed.chars());
        }
    }
}
