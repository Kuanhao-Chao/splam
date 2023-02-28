#ifndef _COMMON_H_
#define _COMMON_H_

#include "var.h"
#include "tmerge.h"
#include <gclib/GStr.h>
#include <unordered_map>
#define VERSION "0.0.1"

extern std::unordered_map<std::string, int>  CHRS;

extern TInputFiles in_records;
extern TInputRecord* irec;

extern CommandMode COMMAND_MODE;

extern GStr infname_model_name;
extern GStr infname_reffa;
extern GStr infname_bam;
extern GStr infname_juncbed;
extern GStr infname_scorebed;

extern GStr out_dir;

extern bool verbose;
extern int aln_num_thr;
extern float threshold;

extern GSamRecord* brec;

extern GStr outfname_multimapped;
extern GStr outfname_spliced;
extern GStr outfname_discard;
extern GStr outfname_cleaned;

extern GSamWriter* outfile_multimapped;
extern GSamWriter* outfile_spliced;
extern GSamWriter* outfile_discard;
extern GSamWriter* outfile_cleaned;
extern FILE* joutf;

extern int JUNC_COUNT;
extern int ALN_COUNT;
extern int ALN_COUNT_SPLICED;
extern int ALN_COUNT_NSPLICED;
extern int ALN_COUNT_BAD;
extern int ALN_COUNT_GOOD;
extern int ALN_COUNT_NH_UPDATE;

extern std::unordered_map<std::string, GSamRecordList> read_hashmap;

extern int STEP_COUNTER;

// j-extract parameters.
extern int j_extract_threshold;
extern GSamWriter* outfile_above_spliced;
extern GSamWriter* outfile_below_spliced;
extern FILE* joutf_above;
extern FILE* joutf_below;

#endif /* TIEBRUSH_TMERGE_H_ */