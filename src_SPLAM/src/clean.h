#ifndef _CLEAN_H_
#define _CLEAN_H_

#include <gclib/GStr.h>
#include "junc.h"
#include <robin_hood/robin_hood.h>

typedef robin_hood::unordered_map<std::string, int> robin_hdd_hm;

void splamClean(int argc, char* argv[]);
GStr filterSpurJuncs(GStr outfname_junc_score, robin_hdd_hm &rm_rd_hm);
void loadBed(GStr inbedname, GArray<CJunc> &spur_juncs);
#endif