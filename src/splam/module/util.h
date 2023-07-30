/*  util.h -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#ifndef _UTIL_H_
#define _UTIL_H_

#include "junc.h"
// #include "common.h"

#include <string>
#include <unordered_map>

std::string get_full_path(std::string fname);

void keepAlignment(GSamWriter* outfile_target, GSamRecord* brec, int& counter);

void reverse_complement(char *str, const hts_pos_t len);

wchar_t *GetWC(const char *c);

int usage_extract();

int usage_clean();

std::unordered_map<std::string, int> get_hg38_chrom_size(std::string target);


void flushBrec(GVec<GSamRecord*> &pbrecs, std::unordered_map<std::string, int> &hits, GSamWriter* outfile_cleaned);

void loadBed(GStr inbedname, GArray<CJunc> &spur_juncs);
void create_CHRS();

#endif /* TIEBRUSH_TMERGE_H_ */
