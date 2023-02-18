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
#include "filter.h"


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

CommandMode COMMAND_MODE = UNSET;
GStr command_str;

GStr infname_model_name;
GStr infname_reffa;
GStr infname_bam;

GStr out_dir;

bool verbose = false;
TInputFiles in_records;
TInputRecord* irec=NULL;

float threshold = 0.2;
// int juncCount = 0;
// GArray<CJunc> junctions(64, true);

GSamRecord* brec=NULL;
GSamWriter* outfile_discard = NULL;
GSamWriter* outfile_spliced = NULL;
GSamWriter* outfile_cleaned = NULL;
FILE* joutf=NULL;


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
    std::filesystem::create_directories(out_dir.chars());

    if (COMMAND_MODE == J_EXTRACT) {
        splamJExtract();
    } else if (COMMAND_MODE == PREDICT) {
        splamPredict();
    } else if (COMMAND_MODE == CLEAN) {
        splamClean(argc, argv);
    }
    return 0;
}

void processOptions(int argc, char* argv[]) {

    GArgs args(argc, argv, "help;cite;verbose;version;SLPEDVvhco:N:Q:m:r:");
    // args.printError(USAGE, true);
    command_str=args.nextNonOpt();
    GMessage(">> command_str       : %s\n", command_str.chars());
    // command_str=args.nextNonOpt();
    // GMessage(">> command_str       : %s\n", command_str.chars());
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
    } else {
        usage();
        GERROR("\n[ERROR] The subcommand must be 'j-extract', 'predict', or 'clean'.\n");
        exit(1);   
    }


    printf("COMMAND_MODE: %d\n", COMMAND_MODE);


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

        GMessage(">>  command_str      : %s\n", command_str.chars());
        GMessage(">> infname_model_name: %s\n", infname_model_name.chars());
        GMessage(">> infname_reffa     : %s\n", infname_reffa.chars());
        GMessage(">> infname_bam       : %s\n", infname_bam.chars());
        GMessage(">> out_dir           : %s\n", out_dir.chars());
    } else if (COMMAND_MODE == PREDICT) {
        processOptionsPredict(args);

    } else if (COMMAND_MODE == CLEAN) {
        processOptionsClean(args);
    }

        GMessage(">> args.startNonOpt()       : %d\n", args.startNonOpt());

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




void processOptionsPredict(GArgs& args) {
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


void processOptionsClean(GArgs& args) {
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