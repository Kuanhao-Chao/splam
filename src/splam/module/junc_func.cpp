/*  junc_func.cpp -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#include "junc.h"
#include "junc_func.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <gclib/GStr.h>
#include <string>

#include "GSam.h"

GArray<CJunc> junctions;

void addJunction(GSamRecord& r, int dupcount, GStr ref) {
	char strand = r.spliceStrand();    
	GSamRecord *r_copy = new GSamRecord(r);
	for (int i=1; i<r.exons.Count(); i++) {
		if (r.exons[i].start-1 - r.exons[i-1].end+1 + 1 <= g_max_splice) {
			CJunc j(r.exons[i-1].end+1, r.exons[i].start-1, strand, ref,
					dupcount);
			int ei;
			int res=junctions.AddIfNew(j, &ei);
			if (res==-1) {
				//existing junction => update junc count
				junctions[ei].add(j);
			} else {
				// new junction 
			}
		}
	}
}

void flushJuncs(FILE* f) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f);
    }
    junctions.Clear();
	// for (int i = 0; i < junctions.Count(); i++) {
	// 	std::cout << i <<  " Junction name: " << junctions[i].start << " - " << junctions[i].end << std::endl;
	// 	std::cout << ">> Read count: " << junctions[i].read_ls.size() << std::endl;
	// 	for (int r = 0; r < junctions[i].read_ls.size(); r++) {
	// 		std::cout << "\tRead " <<  r << " : " << junctions[i].read_ls[r]->cigar() << std::endl;
	// 	}
	// 	std::cout << std::endl;
	// }
	// if (g_j_extract_threshold == 0) {
	// 	junctions.Clear();
	// 	junctions.setCapacity(128);
	// }
}
