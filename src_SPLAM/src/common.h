#ifndef _COMMON_H_
#define _COMMON_H_

#include "var.h"
#include "tmerge.h"
#include <gclib/GStr.h>

extern TInputFiles in_records;
extern TInputRecord* irec;

extern CommandMode COMMAND_MODE;
extern GStr out_dir;

extern GStr model_name;
extern GStr infname_reffa;
extern GStr outfname_junction;

extern bool verbose;
extern float threshold;

#endif /* TIEBRUSH_TMERGE_H_ */