/*  extract.h -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#ifndef _EXTRACT_H_
#define _EXTRACT_H_
#include "bundle.h"
#include <gclib/GStr.h>
#include <gclib/GHashMap.hh>

GStr splamJExtract();
void processBundle_jext(BundleData* bundle, GList<CReadAln>& readlist, int& bundle_counter);
void processRead_jext(int currentstart, int currentend, GList<CReadAln>& readlist, BundleData& bdata, GHash<int>& hashread, CReadAln* alndata);
#endif