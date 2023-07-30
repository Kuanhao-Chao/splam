/*  splam_clean.cpp -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#include <iostream>
#include <vector>
#include <memory>
#include <filesystem>

#include "extract.h"
#include "clean.h"

#include "common.h"
#include "tmerge.h"
#include "util.h"
#include "junc.h"
#include "update.h"
#include "bam_sort.h"

#include <gclib/GArgs.h>
#include <gclib/GBase.h>
#include <gclib/GStr.h>
#include <robin_hood/robin_hood.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Include htslib
#include <htslib/htslib/sam.h>
#include <htslib/htslib/hts.h>

void splam_clean_valid_precheck ();
void processOptions(int argc, char* argv[]);
void processOptionsClean(GArgs& args);

void optionsOutput(GArgs& args);
void optionsThreads(GArgs& args);
void optionsThreshold(GArgs& args);

GStr thread_num = 1;

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
GSamRecord* brec=NULL;

/*********************
 * File name
*********************/
GStr outfname_cleaned;
// Paired
GStr outfname_ns_uniq_map;
GStr outfname_ns_multi_map;
GStr outfname_ns_multi_map_nh_updated;
GStr outfname_s_uniq_map;
GStr outfname_s_uniq_map_cleaned;
GStr outfname_s_multi_map;
GStr outfname_s_multi_map_cleaned;
GStr outfname_s_multi_map_cleaned_nh_updated;
// Unpaired
GStr outfname_ns_uniq_unpair;
GStr outfname_ns_multi_unpair;
GStr outfname_ns_multi_unpair_nh_updated;
GStr outfname_s_uniq_unpair;
GStr outfname_s_uniq_unpair_cleaned;
GStr outfname_s_multi_unpair;
GStr outfname_s_multi_unpair_cleaned;
GStr outfname_s_multi_unpair_cleaned_nh_updated;
// Discarded
GStr outfname_discard_s_uniq_map;
GStr outfname_discard_s_multi_map;

/*********************
 * GSamWriter
*********************/
GSamWriter* outfile_cleaned = NULL;
// Paired
GSamWriter* outfile_ns_uniq_map = NULL;
GSamWriter* outfile_ns_multi_map = NULL;
GSamWriter* outfile_ns_multi_map_nh_updated = NULL;
GSamWriter* outfile_s_uniq_map = NULL;
GSamWriter* outfile_s_uniq_map_cleaned = NULL;
GSamWriter* outfile_s_multi_map = NULL;
GSamWriter* outfile_s_multi_map_cleaned = NULL;
GSamWriter* outfile_s_multi_map_cleaned_nh_updated = NULL;
// Unpaired
GSamWriter* outfile_ns_uniq_unpair = NULL;
GSamWriter* outfile_ns_multi_unpair = NULL;
GSamWriter* outfile_ns_multi_unpair_nh_updated = NULL;
GSamWriter* outfile_s_uniq_unpair;
GSamWriter* outfile_s_uniq_unpair_cleaned = NULL;
GSamWriter* outfile_s_multi_unpair;
GSamWriter* outfile_s_multi_unpair_cleaned = NULL;
GSamWriter* outfile_s_multi_unpair_cleaned_nh_updated = NULL;
// Discarded
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
int ALN_COUNT_SPLICED_UNIQ_KEEP = 0;
int ALN_COUNT_SPLICED_MULTI_KEEP = 0;
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
int ALN_COUNT_SPLICED_UNIQ_UNPAIR_KEEP = 0;
int ALN_COUNT_SPLICED_MULTI_UNPAIR_KEEP = 0;
int ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD = 0;
int ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD = 0;
int ALN_COUNT_NSPLICED_UNIQ_UNPAIR = 0;
int ALN_COUNT_NSPLICED_MULTI_UNPAIR = 0;

robin_hood::unordered_map<std::string, int>  CHRS;

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
// predict parameters:
bool write_bam = true;

int array_size(char *test[]) {
    int size = 0;
    while (test[size] != NULL) {
        size++;
    }
    return size;
}

namespace py = pybind11;

int splam_clean(py::list args_pyls) {
    /*********************
     * Creating argv
    *********************/
    std::vector<std::string> args;
    // Convert each element of the py::list to std::string
    for (const auto& item : args_pyls) {
        std::string str = py::cast<std::string>(item);
        args.push_back(str);
    }
    int argc = args.size();
    char** argv = new char*[args.size() + 1]; // +1 for the null terminator
    // Convert each string to a C-style string and copy it to argv
    for (size_t i = 0; i < args.size(); ++i) {
        argv[i] = new char[args[i].size() + 1]; // +1 for the null terminator
        std::strcpy(argv[i], args[i].c_str());
    }
    // Add a null terminator at the end
    argv[args.size()] = nullptr;

    /*********************
     * Process options.
    *********************/
    processOptions(argc, argv);

    outfname_cleaned = out_dir + "/cleaned.bam";
    /*********************
     * For paired uniq- / multi- mapped alignments
    *********************/
    /*********************
     * Paired File name
    *********************/
    // None
    outfname_ns_uniq_map = out_dir + "/tmp/ns_uniq.bam";
    // (1) NH updated
    outfname_ns_multi_map = out_dir + "/tmp/ns_multi.bam";
    outfname_ns_multi_map_nh_updated = out_dir + "/tmp/ns_multi_nh_updated.bam";
    // (1) cleaned
    outfname_s_uniq_map = out_dir + "/tmp/s_uniq.bam";
    outfname_s_uniq_map_cleaned = out_dir + "/tmp/s_uniq_cleaned.bam";
    // (1) cleaned, (2) NH updated
    outfname_s_multi_map = out_dir + "/tmp/s_multi.bam";
    outfname_s_multi_map_cleaned = out_dir + "/tmp/s_multi_cleaned.bam";
    outfname_s_multi_map_cleaned_nh_updated = out_dir + "/tmp/s_multi_cleaned_nh_updated.bam";
    /*********************
     * Unpaired File name
    *********************/
    // None
    outfname_ns_uniq_unpair = out_dir + "/tmp/ns_uniq_unpair.bam";
    // (1) NH updated
    outfname_ns_multi_unpair = out_dir + "/tmp/ns_multi_unpair.bam";
    outfname_ns_multi_unpair_nh_updated = out_dir + "/tmp/ns_multi_unpair_nh_updated.bam";
    // (1) cleaned
    outfname_s_uniq_unpair = out_dir + "/tmp/s_uniq_unpair.bam";
    outfname_s_uniq_unpair_cleaned = out_dir + "/tmp/s_uniq_unpair_cleaned.bam";
    // (1) cleaned, (2) NH updated
    outfname_s_multi_unpair = out_dir + "/tmp/s_multi_unpair.bam";
    outfname_s_multi_unpair_cleaned = out_dir + "/tmp/s_multi_unpair_cleaned.bam";
    outfname_s_multi_unpair_cleaned_nh_updated = out_dir + "/tmp/s_multi_unpair_cleaned_nh_updated.bam";
    // Discard
    outfname_discard_s_uniq_map = out_dir + "/discard/discard_splice_uniq_map.bam";
    outfname_discard_s_multi_map = out_dir + "/discard/discard_splice_multi_map.bam";

    // The junction score file
    infname_juncbed = out_dir + "/junction_score.bed";

    /*********************
     * Precheck
    *********************/
    splam_clean_valid_precheck();

    GStr discard_dir(out_dir + "/discard");
    std::filesystem::create_directories(out_dir.chars());
    create_CHRS();
    
    /*********************
     * Directory creating + Sam writer creation
    *********************/
    // Creating the directory
    std::filesystem::create_directories(discard_dir.chars());
    GSamReader bam_reader = GSamReader(outfname_ns_uniq_map);

    /*********************
     * Creating GSamWriter
    *********************/
    outfile_ns_multi_map_nh_updated = new GSamWriter(outfname_ns_multi_map_nh_updated, bam_reader.header(), GSamFile_BAM);
    outfile_s_uniq_map_cleaned = new GSamWriter(outfname_s_uniq_map_cleaned, bam_reader.header(), GSamFile_BAM);
    outfile_s_multi_map_cleaned = new GSamWriter(outfname_s_multi_map_cleaned, bam_reader.header(), GSamFile_BAM);
    outfile_s_multi_map_cleaned_nh_updated = new GSamWriter(outfname_s_multi_map_cleaned_nh_updated, bam_reader.header(), GSamFile_BAM);

    if (g_paired_removal) {
        outfile_ns_multi_unpair_nh_updated = new GSamWriter(outfname_ns_multi_unpair_nh_updated, bam_reader.header(), GSamFile_BAM);
        outfile_s_uniq_unpair_cleaned = new GSamWriter(outfname_s_uniq_unpair_cleaned, bam_reader.header(), GSamFile_BAM);
        outfile_s_multi_unpair_cleaned = new GSamWriter(outfname_s_multi_unpair_cleaned, bam_reader.header(), GSamFile_BAM);
        outfile_s_multi_unpair_cleaned_nh_updated = new GSamWriter(outfname_s_multi_unpair_cleaned_nh_updated, bam_reader.header(), GSamFile_BAM);
    }

    // discarded files
    outfile_discard_s_uniq_map = new GSamWriter(outfname_discard_s_uniq_map, bam_reader.header(), GSamFile_BAM);
    outfile_discard_s_multi_map = new GSamWriter(outfname_discard_s_multi_map, bam_reader.header(), GSamFile_BAM);   

    /*********************
     * Main algorithms
    *********************/
    infname_juncbed = out_dir + "/junction.bed";
    infname_scorebed = out_dir + "/junction_score.bed";
    infname_NH_tag = splamClean();
    splamNHUpdate();

    /*********************
     * Merging all files into a clean BAM file
    *********************/

    fprintf(stderr, "Outside of merging all files into a clean BAM file!\n");
    if (g_paired_removal) {
        // Step 1: Sort BAM file
        int argc_sort = 0;
        int res = -1;

        // I do not need to sort the nonsplice alignment
        // /************************
        // * ns_uniq.bam (paired)
        // ************************/        
        // outfname_sort = out_dir + "/tmp/ns_uniq.sort.bam";
        // fprintf(stderr, "outfname_sort: %s\n", outfname_sort.chars());
        // char *ns_uniq_argv[] = {"samtools", "sort", "-@", (char*)thread_num.chars(), (char*)outfname_ns_uniq_map.chars(), "-o", (char*)outfname_sort.chars()};
        // res =bam_sort(argc_sort-1, ns_uniq_argv+1);
        // fprintf(stderr, ">> bam_sort res: %d\n", res);

        // /************************
        // * ns_multi_nh_updated.bam (paired)
        // ************************/ 
        // outfname_sort = out_dir + "/tmp/ns_multi_nh_updated.sort.bam";
        // fprintf(stderr, "outfname_sort: %s\n", outfname_sort.chars());
        // char *ns_multi_argv[] = {"samtools", "sort", "-@", (char*)thread_num.chars(), (char*)outfname_ns_multi_map_nh_updated.chars(), "-o", (char*)outfname_sort.chars()};
        // res = bam_sort(argc_sort-1, ns_multi_argv+1);
        // fprintf(stderr, ">> bam_sort res: %d\n", res);


        /************************
        * Sorting s_uniq_cleaned.bam (paired)
        ************************/ 
        GStr outfname_s_uniq_map_cleaned_sort = out_dir + "/tmp/s_uniq_cleaned.sort.bam";
        char *s_uniq_argv[] = {"samtools", "sort", "-@", (char*)thread_num.chars(), (char*)outfname_s_uniq_map_cleaned.chars(), "-o", (char*)outfname_s_uniq_map_cleaned_sort.chars()};
        argc_sort = sizeof(s_uniq_argv) / sizeof(s_uniq_argv[0]);
        res = bam_sort(argc_sort-1, s_uniq_argv+1);
        // [Important!] Reset argument
        optind = 1;
        
        /************************
        * Sorting s_multi_cleaned_nh_updated.bam (paired)
        ************************/ 
        GStr outfname_s_multi_map_cleaned_nh_updated_sort = out_dir + "/tmp/s_multi_cleaned_nh_updated.sort.bam";
        char *s_multi_argv[] = {"samtools", "sort", "-@", (char*)thread_num.chars(), (char*)outfname_s_multi_map_cleaned_nh_updated.chars(), "-o", (char*)outfname_s_multi_map_cleaned_nh_updated_sort.chars()};
        argc_sort = sizeof(s_multi_argv) / sizeof(s_multi_argv[0]);
        res = bam_sort(argc_sort-1, s_multi_argv+1);
        // [Important!] Reset argument
        optind = 1;

        // Step 2: Merge BAM file
        char *merge_argv[] = {"samtools", "merge", "-f", "-@", (char*)thread_num.chars(), "-o", (char*)outfname_cleaned.chars(), 
        (char*)outfname_ns_uniq_map.chars(), 
        (char*)outfname_ns_multi_map_nh_updated.chars(), 

        (char*)outfname_s_uniq_map_cleaned_sort.chars(), 
        (char*)outfname_s_multi_map_cleaned_nh_updated_sort.chars(), 

        (char*)outfname_ns_uniq_unpair.chars(), 
        (char*)outfname_ns_multi_unpair_nh_updated.chars(), 
        (char*)outfname_s_uniq_unpair_cleaned.chars(), 
        (char*)outfname_s_multi_unpair_cleaned_nh_updated.chars()};

        int argc_merge = sizeof(merge_argv) / sizeof(merge_argv[0]);
        fprintf(stderr, "Numver of input argument!: %d\n", argc_merge);
        bam_merge(argc_merge-1, merge_argv+1);

    } else {
        int argc_merge = 9;
        char *test[] = {"samtools", "merge", "-f", "-@", (char*)thread_num.chars(), "-o", (char*)outfname_cleaned.chars(), (char*)outfname_ns_uniq_map.chars(), (char*)outfname_ns_multi_map_nh_updated.chars(), (char*)outfname_s_uniq_map_cleaned.chars(), (char*)outfname_s_multi_map_cleaned_nh_updated.chars()};
        bam_merge(argc_merge-1, test+1);
    }
    /*********************
     * End of merging all files into a clean BAM file
    *********************/


    /*********************
     * Final statistics printing
    *********************/
    ALN_COUNT_SPLICED = ALN_COUNT_SPLICED_UNIQ + ALN_COUNT_SPLICED_MULTI;
    ALN_COUNT_SPLICED_UNPAIR = ALN_COUNT_SPLICED_UNIQ_UNPAIR + ALN_COUNT_SPLICED_MULTI_UNPAIR;

    ALN_COUNT_BAD = ALN_COUNT_SPLICED_UNIQ_DISCARD + ALN_COUNT_SPLICED_MULTI_DISCARD + ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD + ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD;
    
    ALN_COUNT_GOOD = ALN_COUNT_SPLICED_UNIQ_KEEP + ALN_COUNT_SPLICED_MULTI_KEEP + ALN_COUNT_SPLICED_UNIQ_UNPAIR_KEEP + ALN_COUNT_SPLICED_MULTI_UNPAIR_KEEP;
    ALN_COUNT = ALN_COUNT_GOOD + ALN_COUNT_BAD;
    if (g_paired_removal) {
        GMessage("\n[INFO] Total number of spliced alignments\t:%10d \n", ALN_COUNT);

        // Printing for paired alignment
        GMessage("           paired spliced alignments\t\t:%10d \n", ALN_COUNT_SPLICED);
        GMessage("               - uniquely mapped\t\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_UNIQ, ALN_COUNT_SPLICED_UNIQ_KEEP, ALN_COUNT_SPLICED_UNIQ_DISCARD);
        GMessage("               - multi-mapped\t\t\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_MULTI, ALN_COUNT_SPLICED_MULTI_KEEP, ALN_COUNT_SPLICED_MULTI_DISCARD);

        // Printing for unpaired alignment
        GMessage("           unpaired spliced alignments\t\t:%10d \n", ALN_COUNT_SPLICED_UNPAIR);
        GMessage("               - uniquely mapped\t\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_UNIQ_UNPAIR, ALN_COUNT_SPLICED_UNIQ_UNPAIR_KEEP, ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD);
        GMessage("               - multi-mapped\t\t\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_MULTI_UNPAIR, ALN_COUNT_SPLICED_MULTI_UNPAIR_KEEP, ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD);

        GMessage("\n[INFO] Number of junctions\t\t\t:%10d   (good: %d / bad: %d / unstranded: %d)\n", JUNC_COUNT, JUNC_COUNT_GOOD, JUNC_COUNT_BAD, JUNC_COUNT-JUNC_COUNT_GOOD-JUNC_COUNT_BAD);
        GMessage("\n[INFO] Number of removed spliced alignments\t:%10d \n", ALN_COUNT_BAD);
        GMessage("[INFO] Number of kept spliced alignments\t:%10d \n", ALN_COUNT_GOOD);
    } else {
        GMessage("\n[INFO] Total number of spliced alignments\t:%10d \n", ALN_COUNT);
        GMessage("           spliced alignments\t\t\t:%10d \n", ALN_COUNT_SPLICED);
        GMessage("               - uniquely mapped\t\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_UNIQ, ALN_COUNT_SPLICED_UNIQ_KEEP, ALN_COUNT_SPLICED_UNIQ_DISCARD);
        GMessage("               - multi-mapped\t\t\t:%10d   (kept: %d / removed: %d )\n", ALN_COUNT_SPLICED_MULTI, ALN_COUNT_SPLICED_MULTI_KEEP, ALN_COUNT_SPLICED_MULTI_DISCARD);


        GMessage("\n[INFO] Number of junctions\t\t\t:%10d   (good: %d / bad: %d / unstranded: %d)\n", JUNC_COUNT, JUNC_COUNT_GOOD, JUNC_COUNT_BAD, JUNC_COUNT-JUNC_COUNT_GOOD-JUNC_COUNT_BAD);
        GMessage("\n[INFO] Number of removed spliced alignments\t:%10d \n", ALN_COUNT_BAD);
        GMessage("[INFO] Number of kept spliced alignments\t:%10d \n", ALN_COUNT_GOOD);
    }
    return 0;
}

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


PYBIND11_MODULE(splam_clean, m) {
    // m.doc() = R"pbdoc(
    //     Pybind11 example plugin
    //     -----------------------

    //     .. currentmodule:: bind_test

    //     .. autosummary::
    //        :toctree: _generate

    //        add
    //        subtract
    // )pbdoc";

    m.def("splam_clean", &splam_clean, R"pbdoc(
        Extracting splice junctions
    )pbdoc");


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}


void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;cite;verbose;paired;no-write-bam;output=;threads=;threshold=;hcVSPJo:N:Q:m:r:s:M:g:@:t:");
    // args.printError(usage_clean, true);
    if (argc == 0) {
        usage_clean();
        GERROR("\n[ERROR] No command provide. The subcommand must be 'j-extract', 'predict', or 'clean'.\n");
        exit(1);   
    }

    if (args.getOpt('h') || args.getOpt("help")) {
        usage_clean();
        exit(0);
    }

    if (args.getOpt('c') || args.getOpt("cite")) {
        fprintf(stdout,"%s\n", "Kuan-Hao Chao, Mihaela Pertea, and Steven Salzberg, \033[1m\x1B[3mSPLAM: accurate deep-learning-based splice site predictor to clean up spurious spliced alignments\x1B[0m\033[0m, (2023), GitHub repository, https://github.com/Kuanhao-Chao/SPLAM");
        exit(0);
    }

    if (args.getOpt('P') || args.getOpt("paired") ) {
        g_paired_removal = true;
    }


    if (verbose) {
        GMessage("[INFO] Running in '%s' mode\n\n", argv[1]);
    }
    
    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (verbose) {
        args.printCmdLine(stderr);
    }

    /********************************
     * Process the arguments
    *********************************/
    processOptionsClean(args);

#ifdef DEBUG
    GMessage(">> infname_model_name: %s\n", infname_model_name.chars());
    GMessage(">> infname_juncbed   : %s\n", infname_juncbed.chars());
    GMessage(">> infname_scorebed  : %s\n", infname_scorebed.chars());
    GMessage(">> infname_reffa     : %s\n", infname_reffa.chars());
    GMessage(">> infname_bam       : %s\n", infname_bam.chars());
    GMessage(">> out_dir           : %s\n", out_dir.chars());
    GMessage(">> thread_num        : %s\n", thread_num.chars());
    GMessage(">> g_paired_removal  : %d\n", g_paired_removal);
    GMessage(">> g_max_splice      : %d\n", g_max_splice);
    GMessage(">> g_bundle_gap      : %d\n", g_bundle_gap);
    GMessage(">> verbose           : %d\n", verbose);
#endif

    args.startNonOpt();
}

void processOptionsClean(GArgs& args) {
    optionsOutput(args);
    optionsThreads(args);
    optionsThreshold(args);
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


void optionsThreads(GArgs& args) {
    // -@ / --threads
    GStr s;
    s = args.getOpt('@');
    if (!s.is_empty()) {
        thread_num = s;
    } else {
        s=args.getOpt("threads");
        if (!s.is_empty()) {
            // Use the default bundle-gap
            thread_num = s;
        }
    }
}


void optionsThreshold(GArgs& args) {
    // -t / --threshold
    GStr s;
    s = args.getOpt('t');
    if (!s.is_empty()) {
        threshold = std::atof(s.chars());
        s.asInt();
    } else {
        s=args.getOpt("threshold");
        if (!s.is_empty()) {
            threshold = std::atof(s.chars());
        }
    }
}


void splam_clean_valid_precheck () {
    bool ns_multi_map_bool = fileExists(outfname_ns_multi_map.chars())>1;
    bool ns_uniq_map_bool = fileExists(outfname_ns_uniq_map.chars())>1;
    bool s_multi_map_bool = fileExists(outfname_s_multi_map.chars())>1;
    bool s_uniq_map_bool = fileExists(outfname_s_uniq_map.chars())>1;

    bool ns_multi_unpair_bool = fileExists(outfname_ns_multi_unpair.chars())>1;
    bool ns_uniq_unpair_bool = fileExists(outfname_ns_uniq_unpair.chars())>1;
    bool s_multi_unpair_bool = fileExists(outfname_s_multi_unpair.chars())>1;
    bool s_uniq_unpair_bool = fileExists(outfname_s_uniq_unpair.chars())>1;
    bool infname_juncbed_bool = fileExists(infname_juncbed.chars())>1;

    if (ns_multi_map_bool && ns_uniq_map_bool && s_multi_map_bool && s_uniq_map_bool && 
        ns_multi_unpair_bool && ns_uniq_unpair_bool && s_multi_unpair_bool && s_uniq_unpair_bool && 
        infname_juncbed_bool) {
        g_paired_removal = true;
    } else if (ns_multi_map_bool && ns_uniq_map_bool && s_multi_map_bool && s_uniq_map_bool && 
        !ns_multi_unpair_bool && !ns_uniq_unpair_bool && !s_multi_unpair_bool && !s_uniq_unpair_bool && 
        infname_juncbed_bool) {
        g_paired_removal = false;
    }

    if (g_paired_removal) {
        if (!ns_multi_map_bool || !ns_uniq_map_bool || !s_multi_map_bool || !s_uniq_map_bool || 
            !ns_multi_unpair_bool || !ns_uniq_unpair_bool || !s_multi_unpair_bool || !s_uniq_unpair_bool || 
            !infname_juncbed_bool) {
            
            if (!ns_multi_map_bool) GMessage("[Error] %s' does not exist.\n", outfname_ns_multi_map.chars());
            if (!ns_uniq_map_bool) GMessage("[Error] %s' does not exist.\n", outfname_ns_uniq_map.chars());
            if (!s_multi_map_bool) GMessage("[Error] %s' does not exist.\n", outfname_s_multi_map.chars());
            if (!s_uniq_map_bool) GMessage("[Error] %s' does not exist.\n", outfname_s_uniq_map.chars());
            if (!ns_multi_unpair_bool) GMessage("[Error] %s' does not exist.\n", outfname_ns_multi_unpair.chars());
            if (!ns_uniq_unpair_bool) GMessage("[Error] %s' does not exist.\n", outfname_ns_uniq_unpair.chars());
            if (!s_multi_unpair_bool) GMessage("[Error] %s' does not exist.\n", outfname_s_multi_unpair.chars());
            if (!s_uniq_unpair_bool) GMessage("[Error] %s' does not exist.\n", outfname_s_uniq_unpair.chars());

            if (!infname_juncbed_bool) GMessage("'%s' does not exist.\n", infname_juncbed.chars());
            GError("[Error] file missing.\n");            
        }
    } else {
        if (!ns_multi_map_bool || !ns_uniq_map_bool || !s_multi_map_bool || !s_uniq_map_bool || 
            !infname_juncbed_bool) {
            
            if (!ns_multi_map_bool) GMessage("[Error] %s' does not exist.\n", outfname_ns_multi_map.chars());
            if (!ns_uniq_map_bool) GMessage("[Error] %s' does not exist.\n", outfname_ns_uniq_map.chars());
            if (!s_multi_map_bool) GMessage("[Error] %s' does not exist.\n", outfname_s_multi_map.chars());
            if (!s_uniq_map_bool) GMessage("[Error] %s' does not exist.\n", outfname_s_uniq_map.chars());
            if (!infname_juncbed_bool) GMessage("'%s' does not exist.\n", infname_juncbed.chars());
            GError("[Error] file missing.\n");            
        }
    }
}
