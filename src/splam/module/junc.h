/*  junc.h -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#ifndef _JUNC_H
#define _JUNC_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <gclib/GStr.h>
#include <string>

#include "common.h"
#include "GSam.h"

class CJunc {

	public:
		GStr ref;
		int start, end;
		char strand;
		uint64_t dupcount;
		
		CJunc(int vs=0, int ve=0, char vstrand='+', GStr vref=".", uint64_t dcount=1):
		start(vs), end(ve), strand(vstrand), ref(vref), dupcount(dcount) { }

		bool operator==(const CJunc& a) {
			/************************
			 * Printing check
			************************/
			// if (start==a.start && end==a.end) {
			// 	std::cout << " >> 1 == : " << ref.chars() << ";  strand: " << strand << ";  start: " << start << ";  end: " << end << std::endl;
			// 	std::cout << " >> 2 == : " << a.ref.chars() << ";  strand: " << a.strand << ";  start: " << a.start << ";  end: " << a.end << std::endl;
			// }

			// std::cout << " >> == : " << strcmp(ref.c_str(), a.ref.c_str())==0 && strand==a.strand && start==a.start && end==a.end << std::endl;
			// std::cout << "(strcmp(ref.c_str(), a.ref.c_str())==0 : " << (strcmp(ref.c_str(), a.ref.c_str())==0 ) << std::endl;
			// std::cout << "strand==a.strand : " << (strand==a.strand) << std::endl;
			// std::cout << "end==a.end       : " << (end==a.end) << std::endl;

			// std::cout << "(strcmp(ref.c_str(), a.ref.c_str())==0 && strand==a.strand && start==a.start && end==a.end): " << (strcmp(ref.c_str(), a.ref.c_str())==0 && strand==a.strand && start==a.start && end==a.end) << std::endl;
			return (ref==a.ref && strand==a.strand && start==a.start && end==a.end);
		}

		bool operator<(const CJunc& a) { // sort by strand last
			if (start==a.start){
				if(end==a.end){
					if ((strand=='+') &&  (a.strand=='+')) {
						return false;
					} else if ((strand=='-') &&  (a.strand=='-')) {
						return false;
					} else if ((strand=='+') &&  (a.strand=='-')) {
						return true;
					} else if ((strand=='-') &&  (a.strand=='+')) {
						return false;
					}
					return false;
				}
				else{
					return (end<a.end);
				}
			}
			else{
				return (start<a.start);
			}
		}

		void add(CJunc& j) {
			dupcount+=j.dupcount;
		}

		void write(FILE* f) {
			JUNC_COUNT++;
			fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%ld\t%c\n",
					ref.chars(), start-1, end, JUNC_COUNT, (long)dupcount, strand);
		}
};

#endif