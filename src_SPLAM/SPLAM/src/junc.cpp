#include <iostream>
#include <cstdlib>

#include "GSam.h"

// bool verbose=false;
bool bigwig=false;
int juncCount=0;

struct CJunc {
	int start, end;
	char strand;
	uint64_t dupcount;
	CJunc(int vs=0, int ve=0, char vstrand='+', uint64_t dcount=1):
	  start(vs), end(ve), strand(vstrand), dupcount(dcount) { }

	bool operator==(const CJunc& a) {
		return (strand==a.strand && start==a.start && end==a.end);
	}

    bool operator<(const CJunc& a) { // sort by strand last
        if (start==a.start){
            if(end==a.end){
                return strand<a.strand;
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

	void write(FILE* f, const char* chr) {
		juncCount++;
		fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%ld\t%c\n",
				chr, start-1, end, juncCount, (long)dupcount, strand);
	}
};

GArray<CJunc> junctions(64, true);

void addJunction(GSamRecord& r, int dupcount) {
	char strand = r.spliceStrand();
	for (int i=1;i<r.exons.Count();i++) {
		CJunc j(r.exons[i-1].end+1, r.exons[i].start-1, strand,
				dupcount);
		int ei;
		int r=junctions.AddIfNew(j, &ei);
		if (r==-1) {//existing junction, update
			junctions[ei].add(j);
		}
	}
}

void flushJuncs(FILE* f, const char* chr) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f,chr);
    }
    junctions.Clear();
    junctions.setCapacity(128);
}



