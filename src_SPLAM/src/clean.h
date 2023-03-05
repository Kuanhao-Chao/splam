#ifndef _CLEAN_H_
#define _CLEAN_H_

#include <gclib/GStr.h>
#include <gclib/GHashMap.hh>
#include "junc.h"
#include "bundle.h"
#include "common.h"
#include <robin_hood/robin_hood.h>
#include <string>
#include <unordered_set>

GStr splamClean();
GStr filterSpurJuncs(GStr outfname_junc_score);
void loadBed(GStr inbedname, robin_hdd_string &spur_juncs);
void processBundle(BundleData* bundle, GList<CReadAln>& readlist, robin_hdd_string& rm_juncs, robin_hdd_rm_hit& rm_hit, int& bundle_counter);
void processRead(int currentstart, int currentend, GList<CReadAln>& readlist, BundleData& bdata, GHash<int>& hashread, CReadAln* alndata);
void noMoreBundles();
void removeAlignment(GSamRecord* brec, robin_hdd_rm_hit& rm_hit);
void keepAlignment(GSamRecord* brec);
GStr writenhHitFile(robin_hdd_rm_hit& rm_hit);
bool alignmentAssessment(GSamRecord* brec, robin_hdd_string &rm_juncs);
// std::string get_global_removed_algns_key(GSamRecord* brec);
std::string get_global_removed_mate_algns_key(GSamRecord* brec);
#endif