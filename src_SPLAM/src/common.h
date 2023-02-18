#ifndef _COMMON_H_
#define _COMMON_H_

#include "var.h"
#include "tmerge.h"
#include <gclib/GStr.h>
#define VERSION "0.0.1"

extern TInputFiles in_records;
extern TInputRecord* irec;

extern CommandMode COMMAND_MODE;

extern GStr infname_model_name;
extern GStr infname_reffa;
extern GStr infname_bam;

extern GStr out_dir;

extern bool verbose;
extern float threshold;

extern GSamRecord* brec;

extern GStr outfname_spliced;
extern GStr outfname_nspliced;
extern GStr outfname_discard;
extern GStr outfname_cleaned;

extern GSamWriter* outfile_spliced;
extern GSamWriter* outfile_nspliced;
extern GSamWriter* outfile_discard;
extern GSamWriter* outfile_cleaned;
extern FILE* joutf;

#endif /* TIEBRUSH_TMERGE_H_ */