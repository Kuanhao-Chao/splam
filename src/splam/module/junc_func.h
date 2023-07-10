/*  junc_func.h -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#ifndef _JUNC_FUNC_H
#define _JUNC_FUNC_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <gclib/GStr.h>
#include <string>

#include "common.h"
#include "GSam.h"
#include "junc.h"

extern int juncCount;
extern GArray<CJunc> junctions;
// extern Gset<CJunc> junction;

void addJunction(GSamRecord& r, int dupcount, GStr ref);

void flushJuncs(FILE* f);

#endif