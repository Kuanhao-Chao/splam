#include <iostream>
#include <vector>
#include <memory>
#include <unordered_map>

#include "clean.h"
#include "common.h"
#include "tmerge.h"
#include "util.h"
#include "junc.h"
#include "filter.h"

#include <gclib/GArgs.h>
#include <gclib/GBase.h>
#include <gclib/GStr.h>

#define VERSION "0.0.1"

void processOptions(int argc, char* argv[]);

CommandMode COMMAND_MODE;
GStr out_dir;

GStr model_name;
GStr infname_reffa;
GStr outfname_junction;

bool verbose = false;
TInputFiles in_records;
TInputRecord* irec=NULL;

float threshold = 0.2;
// int juncCount = 0;
// GArray<CJunc> junctions(64, true);

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
    
    if (COMMAND_MODE == PREDICT) {

    } else if (COMMAND_MODE == CLEAN) {
        splamClean();
    } else if (COMMAND_MODE == JUNC_EXTRACT) {

    }

    return 0;
}

void processOptions(int argc, char* argv[]) {

    GArgs args(argc, argv, "help;cite;verbose;version;SLPEDVvhco:N:Q:m:r:");
    // args.printError(USAGE, true);

    if (strcmp(argv[1], "junc-extract") == 0) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
        COMMAND_MODE = JUNC_EXTRACT;
    } else if (strcmp(argv[1], "predict") == 0) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
        COMMAND_MODE = PREDICT;
    } else if (strcmp(argv[1], "clean") == 0) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
        COMMAND_MODE = CLEAN;
    } else {
        usage();
        GERROR("\n[ERROR] The subcommand must be 'junc-extract', 'predict', or 'clean'.\n");
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

    printf("COMMAND_MODE: %d\n", COMMAND_MODE);


    if (args.startNonOpt()==0) {
        usage();
        GMessage("\n[ERROR] no input provided!\n");
        exit(1);
    }

    /********************************
     * Process arguments by COMMAND_MODE
    *********************************/
    if (COMMAND_MODE == PREDICT) {

    } else if (COMMAND_MODE == CLEAN) {
        
    } else if (COMMAND_MODE == JUNC_EXTRACT) {

    }

    // -m / --model
    model_name=args.getOpt('m');
    if (model_name.is_empty()) {
        model_name=args.getOpt('model');
        if (model_name.is_empty()) {
            usage();
            GMessage("\n[ERROR] model file must be provided (-m)!\n");
            exit(1);
        } else {
            if (fileExists(model_name.chars())>1) {
                // guided=true;
            } else {
                GError("[ERROR] model file (%s) not found.\n",
                    model_name.chars());
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
    outfname_junction = out_dir + "/bed/junction.bed";




    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running SPLAM " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    const char* ifn=NULL;
    while ( (ifn=args.nextNonOpt())!=NULL) {
        //input alignment files
        std::string absolute_ifn = get_full_path(ifn);
        std::cout << "absolute_ifn: " << absolute_ifn << std::endl;
        in_records.addFile(absolute_ifn.c_str());
    }
}