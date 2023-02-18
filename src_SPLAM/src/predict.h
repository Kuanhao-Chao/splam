#ifndef _PREDICT_H_
#define _PREDICT_H_
#include <gclib/GStr.h>
#include <htslib/htslib/faidx.h>
#include <robin_hood/robin_hood.h>

typedef robin_hood::unordered_map<std::string, int> dimer_hm;

GStr splamPredict();
faidx_t *fastaIndex();
GStr splamCreateFasta(GStr outfname_junction, dimer_hm &doner_dimers, dimer_hm &acceptor_dimers, faidx_t *ref_faidx);
#endif