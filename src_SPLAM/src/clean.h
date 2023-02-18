#ifndef _CLEAN_H_
#define _CLEAN_H_

#include <gclib/GStr.h>
#include "junc.h"

void splamClean(int argc, char* argv[]);

void loadBed(GStr inbedname, GArray<CJunc> &spur_juncs);
#endif