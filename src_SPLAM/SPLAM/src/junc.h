#include <iostream>
#include <cstdlib>

#include "GSam.h"
#include <gclib/GStr.h>

#include <string>

// using namespace std;

int juncCount=0;


class CJunc {

public:
	GStr ref;
	int start, end;
	char strand;
	uint64_t dupcount;
	
	std::vector<GSamRecord*> read_ls;
	// GArray<CJunc> junctionss;


	CJunc(int vs=0, int ve=0, char vstrand='+', GStr vref=".", uint64_t dcount=1):
	  start(vs), end(ve), strand(vstrand), ref(vref), dupcount(dcount) { }

	bool operator==(const CJunc& a) {
		// std::cout << " >> 1 == : " << ref.c_str() << ";  strand: " << strand << ";  start: " << start << ";  end: " << end << std::endl;
		// std::cout << " >> 2 == : " << a.ref.c_str() << ";  strand: " << a.strand << ";  start: " << a.start << ";  end: " << a.end << std::endl;
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

	void add_read(CJunc& j, GSamRecord* r) {
		read_ls.push_back(r);
	}


	void write(FILE* f) {
		juncCount++;
		fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%ld\t%c\n",
				ref.chars(), start-1, end, juncCount, (long)dupcount, strand);
	}
};

// The reference is in UCSC chr name.
GArray<CJunc> junctions(64, true);

void addJunction(GSamRecord& r, int dupcount, GStr ref) {
	char strand = r.spliceStrand();
//	if (strand!='+' && strand!='-') return; // TODO: should we output .?
    
	// GSamRecord *r_copy = new GSamRecord(r);

	GSamRecord *r_copy = new GSamRecord(r);
	for (int i=1; i<r.exons.Count(); i++) {
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
		// Adding read into the read_ls
		junctions[ei].add_read(j, r_copy);
	}
}


void writeJuncs(FILE* f) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f);
    }

	std::cout << ">> Junction count: " << junctions.Count() << std::endl;
	for (int i = 0; i < junctions.Count(); i++) {
		std::cout << i <<  " Junction name: " << junctions[i].start << " - " << junctions[i].end << std::endl;
		std::cout << ">> Read count: " << junctions[i].read_ls.size() << std::endl;
		for (int r = 0; r < junctions[i].read_ls.size(); r++) {
			std::cout << "\tRead " <<  r << " : " << junctions[i].read_ls[r]->cigar() << std::endl;
		}
		std::cout << std::endl;
	}
}

void flushJuncs(FILE* f) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f);
    }

	std::cout << ">> Junction count: " << junctions.Count() << std::endl;
	for (int i = 0; i < junctions.Count(); i++) {
		std::cout << i <<  " Junction name: " << junctions[i].start << " - " << junctions[i].end << std::endl;
		std::cout << ">> Read count: " << junctions[i].read_ls.size() << std::endl;
		for (int r = 0; r < junctions[i].read_ls.size(); r++) {
			std::cout << "\tRead " <<  r << " : " << junctions[i].read_ls[r]->cigar() << std::endl;
		}
		std::cout << std::endl;
	}

    junctions.Clear();
    junctions.setCapacity(128);
}

