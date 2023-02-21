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

extern GStr out_dir;

extern bool verbose;
extern int aln_num_thr;
extern float threshold;

extern GSamRecord* brec;

extern GStr outfname_spliced;
extern GStr outfname_discard;
extern GStr outfname_cleaned;

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

extern bool skip_extact;

#endif /* TIEBRUSH_TMERGE_H_ */