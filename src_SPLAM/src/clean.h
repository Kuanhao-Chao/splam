#ifndef _CLEAN_H_
#define _CLEAN_H_

#include <gclib/GStr.h>
#include "junc.h"
#include <robin_hood/robin_hood.h>
#include <string>
#include <unordered_set>

typedef robin_hood::unordered_map<std::string, int> robin_hdd_hm;

std::unordered_set<std::string>* splamClean(int argc, char* argv[]);
GStr filterSpurJuncs(GStr outfname_junc_score, robin_hdd_hm &rm_rd_hm, std::unordered_set<std::string> &rm_rd_set);
void loadBed(GStr inbedname, std::unordered_set<std::string> &spur_juncs);
#endif