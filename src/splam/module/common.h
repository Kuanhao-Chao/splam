/*  common.h -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#ifndef _COMMON_H_
#define _COMMON_H_

#include "var.h"
#include "tmerge.h"
#include <gclib/GStr.h>
#include <robin_hood/robin_hood.h>
#define VERSION "0.0.1"

typedef robin_hood::unordered_map<std::string, int> robin_hdd_rm_hit;
typedef robin_hood::unordered_set<std::string> robin_hdd_rm_algn;
typedef robin_hood::unordered_set<int> robin_hdd_int;
typedef robin_hood::unordered_set<std::string> robin_hdd_string;

extern robin_hood::unordered_map<std::string, int>  CHRS;

extern TInputFiles in_records;
extern TInputRecord* irec;

extern GStr thread_num;

extern GStr infname_model_name;
extern GStr infname_reffa;
extern GStr infname_bam;
extern GStr infname_juncbed;
extern GStr infname_scorebed;
extern GStr infname_NH_tag;

extern GStr out_dir;

extern bool verbose;
extern float threshold;
extern GSamRecord* brec;

/*********************
 * File name
*********************/
extern GStr outfname_cleaned;
// Paired
extern GStr outfname_ns_uniq_map;
extern GStr outfname_ns_multi_map;
extern GStr outfname_ns_multi_map_nh_updated;
extern GStr outfname_s_uniq_map;
extern GStr outfname_s_uniq_map_cleaned;
extern GStr outfname_s_multi_map;
extern GStr outfname_s_multi_map_cleaned;
extern GStr outfname_s_multi_map_cleaned_nh_updated;
// Unpaired
extern GStr outfname_ns_uniq_unpair;
extern GStr outfname_ns_multi_unpair;
extern GStr outfname_ns_multi_unpair_nh_updated;
extern GStr outfname_s_uniq_unpair;
extern GStr outfname_s_uniq_unpair_cleaned;
extern GStr outfname_s_multi_unpair;
extern GStr outfname_s_multi_unpair_cleaned;
extern GStr outfname_s_multi_unpair_cleaned_nh_updated;
// Discarded
extern GStr outfname_discard_s_uniq_map;
extern GStr outfname_discard_s_multi_map;

/*********************
 * GSamWriter
*********************/
extern GSamWriter* outfile_cleaned;
// Paired
extern GSamWriter* outfile_ns_uniq_map;
extern GSamWriter* outfile_ns_multi_map;
extern GSamWriter* outfile_ns_multi_map_nh_updated;
extern GSamWriter* outfile_s_uniq_map;
extern GSamWriter* outfile_s_uniq_map_cleaned;
extern GSamWriter* outfile_s_multi_map;
extern GSamWriter* outfile_s_multi_map_cleaned;
extern GSamWriter* outfile_s_multi_map_cleaned_nh_updated;
// Unpaired
extern GSamWriter* outfile_ns_uniq_unpair;
extern GSamWriter* outfile_ns_multi_unpair;
extern GSamWriter* outfile_ns_multi_unpair_nh_updated;
extern GSamWriter* outfile_s_uniq_unpair;
extern GSamWriter* outfile_s_uniq_unpair_cleaned;
extern GSamWriter* outfile_s_multi_unpair;
extern GSamWriter* outfile_s_multi_unpair_cleaned;
extern GSamWriter* outfile_s_multi_unpair_cleaned_nh_updated;
// Discarded
extern GSamWriter* outfile_discard_s_uniq_map;
extern GSamWriter* outfile_discard_s_multi_map;

extern FILE* joutf;

// ALN summary
extern int ALN_COUNT;
extern int ALN_COUNT_BAD;
extern int ALN_COUNT_GOOD;
extern int ALN_COUNT_GOOD_CAL;


// JUNC summary
extern int JUNC_COUNT;
extern int JUNC_COUNT_GOOD;
extern int JUNC_COUNT_BAD;

// Paired
extern int ALN_COUNT_SPLICED;
extern int ALN_COUNT_NSPLICED;
extern int ALN_COUNT_PAIRED;
extern int ALN_COUNT_SPLICED_UNIQ;
extern int ALN_COUNT_SPLICED_MULTI;
extern int ALN_COUNT_SPLICED_UNIQ_KEEP;
extern int ALN_COUNT_SPLICED_MULTI_KEEP;
extern int ALN_COUNT_SPLICED_UNIQ_DISCARD;
extern int ALN_COUNT_SPLICED_MULTI_DISCARD;
extern int ALN_COUNT_NSPLICED_UNIQ;
extern int ALN_COUNT_NSPLICED_MULTI;

// Unpaired
extern int ALN_COUNT_SPLICED_UNPAIR;
extern int ALN_COUNT_NSPLICED_UNPAIR;
extern int ALN_COUNT_UNPAIR;
extern int ALN_COUNT_SPLICED_UNIQ_UNPAIR;
extern int ALN_COUNT_SPLICED_MULTI_UNPAIR;
extern int ALN_COUNT_SPLICED_UNIQ_UNPAIR_KEEP;
extern int ALN_COUNT_SPLICED_MULTI_UNPAIR_KEEP;
extern int ALN_COUNT_SPLICED_UNIQ_UNPAIR_DISCARD;
extern int ALN_COUNT_SPLICED_MULTI_UNPAIR_DISCARD;
extern int ALN_COUNT_NSPLICED_UNIQ_UNPAIR;
extern int ALN_COUNT_NSPLICED_MULTI_UNPAIR;


extern int STEP_COUNTER;

// j-extract parameters.
extern int g_max_splice;
extern int g_bundle_gap;
extern GSamWriter* outfile_above_spliced;
extern GSamWriter* outfile_below_spliced;
extern FILE* joutf_above;
extern FILE* joutf_below;

// predict parameters
extern bool write_bam;

// clean parameters
extern bool g_paired_removal;
extern bool g_2_stage_run;

#endif /* TIEBRUSH_TMERGE_H_ */