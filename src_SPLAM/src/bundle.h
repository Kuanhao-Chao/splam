#ifndef _BUNDLE_H_
#define _BUNDLE_H_

#include "gclib/GStr.h"
#include "gclib/GList.hh"
#include "GSam.h"

// struct CReadAln:public GSeg {
// 	//DEBUG ONLY:
// 	// GStr name;
// 	char strand; // 1, 0 (unkown), -1 (reverse)
// 	short int nh;
// 	uint len;
// 	float read_count;       // keeps count for all reads (including paired and unpaired)
// 	bool unitig:1;			// set if read come from an unitig
// 	bool longread:1;	    // set if read comes from long read data
// 	GVec<float> pair_count;   // keeps count for all paired reads
// 	GVec<int> pair_idx;     // keeps indeces for all pairs in assembly mode, or all reads that were collapsed in merge mode
// 	GVec<GSeg> segs; //"exons"
// 	// GPVec<CJunction> juncs;

// 	CReadAln(char _strand=0, short int _nh=0,
// 			int rstart=0, int rend=0): GSeg(rstart, rend), //name(rname),
// 					strand(_strand),nh(_nh), len(0), read_count(0), unitig(false),longread(false),pair_count(),pair_idx(),
// 					segs() {}
//                     // , juncs(false) { }
// 	CReadAln(CReadAln &rd):GSeg(rd.start,rd.end) { // copy contructor
// 		strand=rd.strand;
// 		nh=rd.nh;
// 		len=rd.len;
// 		read_count=rd.read_count;
// 		unitig=rd.unitig;
// 		longread=rd.longread;
// 		pair_count=rd.pair_count;
// 		pair_idx=rd.pair_idx;
// 	}
// 	int overlapSegLen(CReadAln* r) {

// 		if (r->start>end || start>r->end) return 0;

// 		int i=0;
// 		int j=0;
// 		int len=0;
// 		while(i<segs.Count()) {
// 			if(segs[i].end<r->segs[j].start) i++;
// 			else if(r->segs[j].end<segs[i].start) j++;
// 			else { // there is overlap
// 				len+=segs[i].overlapLen(r->segs[j].start,r->segs[j].end);
// 				if(segs[i].end<r->segs[j].end) i++;
// 				else j++;
// 			}
// 			if(j==r->segs.Count()) break;
// 		}
// 		return len;
// 	}
// 	~CReadAln() { }
// };

struct GReadAlnData {
	GSamRecord* brec;
	char strand; //-1, 0, 1
	int nh;
	int hi;
	// GPVec<CJunction> juncs;
	//GPVec< GVec<RC_ExonOvl> > g_exonovls; //>5bp overlaps with guide exons, for each read "exon"
	GReadAlnData(GSamRecord* bamrec=NULL, char nstrand=0, int num_hits=0,
			int hit_idx=0):brec(bamrec), strand(nstrand),
					nh(num_hits), hi(hit_idx) { } //, g_exonovls(true)
	~GReadAlnData() { }
};

// bundle data structure, holds all data needed for
// infering transcripts from a bundle
enum BundleStatus {
	BUNDLE_STATUS_CLEAR=0, //available for loading/prepping
	BUNDLE_STATUS_LOADING, //being prepared by the main thread (there can be only one)
	BUNDLE_STATUS_READY //ready to be processed, or being processed
};

// bundle data structure, holds all input data parsed from BAM file
struct BundleData {
	BundleStatus status;
	//int64_t bamStart; //start of bundle in BAM file
	int idx; //index in the main bundles array
	int start;
	int end;
	unsigned long numreads; // number of reads in this bundle
	/*
	float wnumreads; // NEW: weighted numreads; a multi-mapped read mapped in 2 places will contribute only 0.5
	double sumreads; // sum of all reads' lengths in bundle
	double sumfrag; // sum of all fragment lengths (this includes the insertion so it is an estimate)
	float num_reads; // number of all reads in bundle that we considered (weighted)
	float num_cov; // how many coverages we added (weighted) to obtain sumcov
	float num_frag; // how many fragments we added to obtain sumfrag
	double num_fragments3;
	double sum_fragments3;
	*/
	double num_fragments; //aligned read/pairs
	double frag_len;
	double sum_cov; // sum of all transcripts coverages --> needed to compute TPMs
	char covflags;

	GStr refseq; //reference sequence name
	char* gseq; //actual genomic sequence for the bundle
	GList<GSamRecord> readlist;
	//  GList<CJunction> junction;
	BundleData():status(BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0),
			numreads(0),
			num_fragments(0), frag_len(0),sum_cov(0),covflags(0),
			refseq(), gseq(NULL), readlist(false,true) {
				// , //bpcov(1024), junction(true, true, true) {
	}

	void getReady(int currentstart, int currentend) {
		GMessage("getReady bundle!\n");
		//this is only called when the bundle is valid and ready to be processed
		start=currentstart;
		end=currentend;
		//refseq=ref;
		//tag all these guides
		status=BUNDLE_STATUS_READY;
		GMessage("Done getReady bundle!\n");
	}

	//bool evalReadAln(GSamRecord& brec, char& strand, int nh); //, int hi);
	//  bool evalReadAln(GReadAlnData& alndata, char& strand);

	void Clear() {
		GMessage("Clearing bundle!\n");
		readlist.Clear();
		// junction.Clear();
		start=0;
		end=0;
		status=BUNDLE_STATUS_CLEAR;
		numreads=0;
		num_fragments=0;
		frag_len=0;
		sum_cov=0;
		covflags=0;
		GFREE(gseq);
		GMessage("~ Done Clearing bundle!\n");
	}

	~BundleData() {
		Clear();
	}
};


#endif