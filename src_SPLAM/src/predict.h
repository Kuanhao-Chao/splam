#ifndef _PREDICT_H_
#define _PREDICT_H_
#include <gclib/GStr.h>
#include <htslib/htslib/faidx.h>
#include <robin_hood/robin_hood.h>

typedef robin_hood::unordered_map<std::string, int> robin_hdd_rm_hit;

GStr splamPredict();
faidx_t *fastaIndex();
GStr splamCreateFasta(GStr outfname_junction, robin_hdd_rm_hit &doner_dimers, robin_hdd_rm_hit &acceptor_dimers, faidx_t *ref_faidx);
#endif