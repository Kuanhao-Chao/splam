#include <iostream>
#include <cstdlib>

#include "GSam.h"
#include <string>

using namespace std;

int juncCount=0;


class CJunc {

public:
	string ref;
	int start, end;
	char strand;
	uint64_t dupcount;
	
	vector<GSamRecord*> read_ls;
	// GArray<CJunc> junctionss;


	CJunc(int vs=0, int ve=0, char vstrand='+', string vref=".", uint64_t dcount=1):
	  start(vs), end(ve), strand(vstrand), ref(vref), dupcount(dcount) { }

	bool operator==(const CJunc& a) {
		// cout << " >> 1 == : " << ref.c_str() << ";  strand: " << strand << ";  start: " << start << ";  end: " << end << endl;
		// cout << " >> 2 == : " << a.ref.c_str() << ";  strand: " << a.strand << ";  start: " << a.start << ";  end: " << a.end << endl;
		// cout << " >> == : " << strcmp(ref.c_str(), a.ref.c_str())==0 && strand==a.strand && start==a.start && end==a.end << endl;
		// cout << "(strcmp(ref.c_str(), a.ref.c_str())==0 : " << (strcmp(ref.c_str(), a.ref.c_str())==0 ) << endl;
		// cout << "strand==a.strand : " << (strand==a.strand) << endl;
		// cout << "end==a.end       : " << (end==a.end) << endl;

		// cout << "(strcmp(ref.c_str(), a.ref.c_str())==0 && strand==a.strand && start==a.start && end==a.end): " << (strcmp(ref.c_str(), a.ref.c_str())==0 && strand==a.strand && start==a.start && end==a.end) << endl;
		return (strcmp(ref.c_str(), a.ref.c_str())==0 && strand==a.strand && start==a.start && end==a.end);
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
				ref.c_str(), start-1, end, juncCount, (long)dupcount, strand);
	}
};

// The reference is in UCSC chr name.
GArray<CJunc> junctions(64, true);

void addJunction(GSamRecord& r, int dupcount, string ref) {
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

	cout << ">> Junction count: " << junctions.Count() << endl;
	for (int i = 0; i < junctions.Count(); i++) {
		cout << i <<  " Junction name: " << junctions[i].start << " - " << junctions[i].end << endl;
		cout << ">> Read count: " << junctions[i].read_ls.size() << endl;
		for (int r = 0; r < junctions[i].read_ls.size(); r++) {
			cout << "\tRead " <<  r << " : " << junctions[i].read_ls[r]->cigar() << endl;
		}
		cout << endl;
	}
}

void flushJuncs(FILE* f) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f);
    }

	cout << ">> Junction count: " << junctions.Count() << endl;
	for (int i = 0; i < junctions.Count(); i++) {
		cout << i <<  " Junction name: " << junctions[i].start << " - " << junctions[i].end << endl;
		cout << ">> Read count: " << junctions[i].read_ls.size() << endl;
		for (int r = 0; r < junctions[i].read_ls.size(); r++) {
			cout << "\tRead " <<  r << " : " << junctions[i].read_ls[r]->cigar() << endl;
		}
		cout << endl;
	}

    junctions.Clear();
    junctions.setCapacity(128);
}

