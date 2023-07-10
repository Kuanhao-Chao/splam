/*  update.h -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#ifndef _UPDATE_H_
#define _UPDATE_H_
#include <gclib/GStr.h>
#include "common.h"

GStr splamNHUpdate();

void readnhHitFile(robin_hdd_rm_hit& rm_hit);

void update_NH_tag_write_alignment(GStr infname, GSamWriter *outfile, int& processed_aln, robin_hdd_rm_hit rm_hit, GStr log);
#endif