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

extern CommandMode COMMAND_MODE;

extern GStr infname_model_name;
extern GStr infname_reffa;
extern GStr infname_bam;
extern GStr infname_juncbed;
extern GStr infname_scorebed;
extern GStr infname_NH_tag;

extern GStr out_dir;

extern bool verbose;
extern int aln_num_thr;
extern float threshold;

extern GSamRecord* brec;


// output file names 
extern GStr outfname_cleaned;

extern GStr outfname_ns_multi_map;
extern GStr outfname_s_uniq_map;
extern GStr outfname_s_multi_map;
extern GStr outfname_s_multi_map_tmp;
extern GStr outfname_discard_unpair;
extern GStr outfname_discard_s_uniq_map;
extern GStr outfname_discard_s_multi_map;

// GSamWriter 
extern GSamWriter* outfile_cleaned;

extern GSamWriter* outfile_ns_multi_map;
extern GSamWriter* outfile_s_uniq_map;
extern GSamWriter* outfile_s_multi_map;
extern GSamWriter* outfile_s_multi_map_tmp;
extern GSamWriter* outfile_discard_unpair;
extern GSamWriter* outfile_discard_s_uniq_map;
extern GSamWriter* outfile_discard_s_multi_map;


extern FILE* joutf;

extern int ALN_COUNT;
extern int JUNC_COUNT;
extern int JUNC_COUNT_GOOD;
extern int JUNC_COUNT_BAD;
extern int ALN_COUNT_UNPAIRED;
extern int ALN_COUNT_SPLICED;
extern int ALN_COUNT_NSPLICED;
extern int ALN_COUNT_SPLICED_UNIQ;
extern int ALN_COUNT_SPLICED_MULTI;
extern int ALN_COUNT_SPLICED_UNIQ_DISCARD;
extern int ALN_COUNT_SPLICED_MULTI_DISCARD;
extern int ALN_COUNT_NSPLICED_UNIQ;
extern int ALN_COUNT_NSPLICED_MULTI;
extern int ALN_COUNT_NSPLICED_UNIQ_DISCARD;
extern int ALN_COUNT_NSPLICED_MULTI_DISCARD;
extern int ALN_COUNT_BAD;
extern int ALN_COUNT_GOOD;
extern int ALN_COUNT_GOOD_CAL;

extern int STEP_COUNTER;

// j-extract parameters.
extern int g_j_extract_threshold;
extern int g_max_splice;
extern int g_bundle_gap;
extern GSamWriter* outfile_above_spliced;
extern GSamWriter* outfile_below_spliced;
extern FILE* joutf_above;
extern FILE* joutf_below;

// clean parameters
// extern robin_hood::unordered_map<std::string, GSamRecordList> read_hashmap;
// extern robin_hood::unordered_set<std::string>* rm_rd_set;

extern bool g_is_single_end;

#endif /* TIEBRUSH_TMERGE_H_ */