/*  splam_extract.cpp -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

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

#include <string>
#include <vector>
#include <gclib/GArgs.h>
#include <gclib/GBase.h>
#include <gclib/GStr.h>
#include <robin_hood/robin_hood.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

void processOptions(int argc, char* argv[]);
void processOptionsJExtract(GArgs& args);

void optionsMaxSplice(GArgs& args);
void optionsBundleGap(GArgs& args);
void optionsOutput(GArgs& args);
void optionsWriteJuncOnly(GArgs& args);

CommandMode COMMAND_MODE = UNSET;

// input file names 
GStr infname_bam("");
GStr infname_juncbed("");
GStr out_dir;

bool verbose = false;
TInputFiles in_records;
TInputRecord* irec=NULL;

float threshold = 0.1;
GSamRecord* brec=NULL;


/******************************
 * Parameters for output files
******************************/
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
// Discarded
GStr outfname_discard_s_uniq_map;
GStr outfname_discard_s_multi_map;

// GSamWriter 
GSamWriter* outfile_cleaned = NULL;
// Paired
GSamWriter* outfile_ns_multi_map = NULL;
GSamWriter* outfile_ns_uniq_map = NULL;
GSamWriter* outfile_s_multi_map = NULL;
GSamWriter* outfile_s_uniq_map = NULL;
GSamWriter* outfile_s_multi_map_cleaned = NULL;
// Unpaired
GSamWriter* outfile_ns_multi_unpair = NULL;
GSamWriter* outfile_ns_uniq_unpair = NULL;
GSamWriter* outfile_s_multi_unpair = NULL;
GSamWriter* outfile_s_uniq_unpair = NULL;
GSamWriter* outfile_s_multi_unpair_tmp = NULL;
FILE* joutf=NULL;


/******************************
 * ALN counters 
******************************/
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
int STEP_COUNTER = 0;


/******************************
 * j-extract parameters:
******************************/
int g_max_splice = 100000;
int g_bundle_gap = 1000;
bool write_bam = true;
bool g_paired_removal = false;

robin_hood::unordered_map<std::string, int>  CHRS;

namespace py = pybind11;

int splam_extract(py::list args_pyls) {
    std::vector<std::string> args;
    // 1. Convert each element of the py::list to std::string
    for (const auto& item : args_pyls) {
        std::string str = py::cast<std::string>(item);
        args.push_back(str);
    }
    int argc = args.size();
    char** argv = new char*[args.size() + 1]; // +1 for the null terminator
    // 2. Convert each string to a C-style string and copy it to argv
    for (size_t i = 0; i < args.size(); ++i) {
        argv[i] = new char[args[i].size() + 1]; // +1 for the null terminator
        std::strcpy(argv[i], args[i].c_str());
    }
    // 3. Add a null terminator at the end
    argv[args.size()] = nullptr;
    in_records.setup(VERSION, argc, argv);
    processOptions(argc, argv);

    /*********************
     * For paired uniq- / multi- mapped alignments
    *********************/
    outfname_ns_uniq_map = out_dir + "/tmp/ns_uniq.bam";
    outfname_ns_multi_map = out_dir + "/tmp/ns_multi.bam";
    outfname_s_uniq_map = out_dir + "/tmp/s_uniq.bam";
    outfname_s_multi_map = out_dir + "/tmp/s_multi.bam";

    /*********************
     * For unpaired uniq- / multi- mapped alignments
    *********************/
    outfname_ns_uniq_unpair = out_dir + "/tmp/ns_uniq_unpair.bam";
    outfname_ns_multi_unpair = out_dir + "/tmp/ns_multi_unpair.bam";
    outfname_s_uniq_unpair = out_dir + "/tmp/s_uniq_unpair.bam";
    outfname_s_multi_unpair = out_dir + "/tmp/s_multi_unpair.bam";

    GStr tmp_dir(out_dir + "/tmp");
    std::filesystem::create_directories(out_dir.chars());
    create_CHRS();
    
    // Start reading the files
    int num_samples=in_records.start();

    /*********************
     * Sam writer creation
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

    /*********************
     * Main algorithms
    *********************/
    infname_juncbed = splamJExtract();

    /*********************
     * Printing results
    *********************/
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
        GMessage("                   - uniquely mapped\t:%10d\n", ALN_COUNT_SPLICED_UNIQ);
        GMessage("                   - multi-mapped\t:%10d\n", ALN_COUNT_SPLICED_MULTI);
        GMessage("               non-spliced alignments\t:%10d \n", ALN_COUNT_NSPLICED);
        GMessage("                   - uniquely mapped\t:%10d\n", ALN_COUNT_NSPLICED_UNIQ);
        GMessage("                   - multi-mapped\t:%10d\n\n", ALN_COUNT_NSPLICED_MULTI);

        // Printing for unpaired alignment
        GMessage("           unpaired alignments\t\t:%10d \n", ALN_COUNT_UNPAIR);
        GMessage("               spliced alignments\t:%10d \n", ALN_COUNT_SPLICED_UNPAIR);
        GMessage("                   - uniquely mapped\t:%10d\n", ALN_COUNT_SPLICED_UNIQ_UNPAIR);
        GMessage("                   - multi-mapped\t:%10d\n", ALN_COUNT_SPLICED_MULTI_UNPAIR);
        GMessage("               non-spliced alignments\t:%10d \n", ALN_COUNT_NSPLICED_UNPAIR);
        GMessage("                   - uniquely mapped\t:%10d\n", ALN_COUNT_NSPLICED_UNIQ_UNPAIR);
        GMessage("                   - multi-mapped\t:%10d\n\n", ALN_COUNT_NSPLICED_MULTI_UNPAIR);
    } else {
        GMessage("\n[INFO] Total number of alignments\t:%10d \n", ALN_COUNT);
        GMessage("           spliced alignments\t\t:%10d \n", ALN_COUNT_SPLICED);
        GMessage("               - uniquely mapped\t:%10d\n", ALN_COUNT_SPLICED_UNIQ);
        GMessage("               - multi-mapped\t\t:%10d\n", ALN_COUNT_SPLICED_MULTI);
        
        GMessage("           non-spliced alignments\t:%10d \n", ALN_COUNT_NSPLICED);
        GMessage("               - uniquely mapped\t:%10d\n", ALN_COUNT_NSPLICED_UNIQ);
        GMessage("               - multi-mapped\t\t:%10d\n", ALN_COUNT_NSPLICED_MULTI);
    }
    GMessage("\n[INFO] Total number of junctions\t:%10d\n", JUNC_COUNT);
    return 0;
}


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


PYBIND11_MODULE(splam_extract, m) {
    m.def("splam_extract", &splam_extract, R"pbdoc(
        Extracting splice junctions
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}



void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;cite;verbose;paired;write-junctions-only;output=;max-splice=;bundle-gap=;hcPVno:M:g:");

    if (args.getOpt('h') || args.getOpt("help")) {
        usage_extract();
        exit(0);
    }

    if (args.getOpt('c') || args.getOpt("cite")) {
        fprintf(stdout,"%s\n", "Kuan-Hao Chao, Mihaela Pertea, and Steven Salzberg, \033[1m\x1B[3mSPLAM: accurate deep-learning-based splice site predictor to clean up spurious spliced alignments\x1B[0m\033[0m, (2023), GitHub repository, https://github.com/Kuanhao-Chao/SPLAM");
        exit(0);
    }

    if (args.getOpt('P') || args.getOpt("paired") ) {
        g_paired_removal = true;
    }
    
    verbose=(args.getOpt('V')!=NULL || args.getOpt("verbose")!=NULL);
    if (verbose) {
        // fprintf(stderr, "Running SPLAM " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }

    /********************************
     * Process arguments by COMMAND_MODE
    *********************************/
    processOptionsJExtract(args);    

#ifdef DEBUG
    GMessage(">> infname_juncbed   : %s\n", infname_juncbed.chars());
    GMessage(">> infname_bam       : %s\n", infname_bam.chars());
    GMessage(">> out_dir           : %s\n", out_dir.chars());
    GMessage(">> g_paired_removal  : %d\n", g_paired_removal);
    GMessage(">> g_max_splice      : %d\n", g_max_splice);
    GMessage(">> g_bundle_gap      : %d\n", g_bundle_gap);
    GMessage(">> verbose           : %d\n", verbose);
#endif

    args.startNonOpt();
    if (args.getNonOptCount()==0) {
        usage_extract();
        GMessage("\n[ERROR] no input provided!\n");
        exit(1);
    }

    const char* ifn=NULL;
    while ( (ifn=args.nextNonOpt())!=NULL) {
        //input alignment files
        std::string absolute_ifn = get_full_path(ifn);
        in_records.addFile(absolute_ifn.c_str());
    }
}


void processOptionsJExtract(GArgs& args) {
    optionsOutput(args);
    optionsMaxSplice(args);
    optionsBundleGap(args);
    optionsWriteJuncOnly(args);
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
            out_dir = "tmp_out";
        }
    }
    if (out_dir[out_dir.length()-1] == '/' ) {
		out_dir.cut(out_dir.length()-1, 1);
    }
}

void optionsWriteJuncOnly(GArgs& args) {
    // -n / --write-junctions-only
	if (args.getOpt('n') || args.getOpt("write-junctions-only")) {
        write_bam = false;
	}
}