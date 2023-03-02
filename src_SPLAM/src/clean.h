#ifndef _CLEAN_H_
#define _CLEAN_H_

#include <gclib/GStr.h>
#include <gclib/GHashMap.hh>
#include "junc.h"
#include "bundle.h"
#include <robin_hood/robin_hood.h>
#include <string>
#include <unordered_set>

typedef robin_hood::unordered_map<std::string, int> robin_hdd_rm_hit;
typedef robin_hood::unordered_set<std::string> robin_hdd_rm_algn;
typedef robin_hood::unordered_set<int> robin_hdd_int;

GStr splamClean();
GStr filterSpurJuncs(GStr outfname_junc_score);
void loadBed(GStr inbedname, std::unordered_set<std::string> &spur_juncs);
void processBundle(BundleData* bundle, GList<CReadAln>& readlist, std::unordered_set<std::string>& rm_juncs, int& bundle_counter);
void processRead(int currentstart, int currentend, GList<CReadAln>& readlist, BundleData& bdata, GHash<int>& hashread, CReadAln* alndata);
void noMoreBundles();
#endif