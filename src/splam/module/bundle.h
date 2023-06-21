#ifndef _BUNDLE_H_
#define _BUNDLE_H_

#include "gclib/GStr.h"
#include "gclib/GList.hh"
#include "GSam.h"

struct CReadAln:public GSeg {
	//DEBUG ONLY:
	GSamRecord brec;
	int pair_idx;     // keeps index for the pair of the alignment.
	// GVec<GSeg> segs; //"exons"

	CReadAln(GSamRecord* bamrec): 
			GSeg(bamrec->start, bamrec->end), brec(NULL), pair_idx(-1) {

		// GMessage("Inside CReadAln constructore!\n");
		// GMessage("CReadAln constructor called\n");

		// GMessage("bamrec: %s\n", bamrec->refName());
		// GMessage("==> Before      : %p\n", &bamrec);
		// GMessage("==> Before brec  : %s\n", bamrec->refName());
		// GSamRecord brec_tmp = *bamrec;
		this->brec = (*bamrec);

		// this->brec = bamrec;

		// this->brec = &brec_tmp;

		// * brec = brec_tmp;

		// brec = &brec_tmp;
		// brec = bamrec;
		// GMessage("After bamrec: %s\n", brec_tmp.refName());
		// GMessage("After brec  : %p\n", &brec_tmp);
		// GMessage("==> After brec  : %p\n", &brec);
		// GMessage("==> After brecrefname: %s\n", this->brec.refName());
		// GMessage("==> After materefname: %s\n", this->brec.mate_refName());
		// GMessage("==> After brec start : %d\n", brec.start);
		// GMessage("==> After brec start : %p\n", &(brec.start));
		// GMessage("==> After brec end : %d\n", brec.end);
		// GMessage("==> After brec end : %p\n", &(brec.end));
	}
	
	// CReadAln(CReadAln &rd):GSeg(rd.start,rd.end) { // copy contructor
	// 	pair_idx=rd.pair_idx;
	// 	GSamRecord brec_tmp = GSamRecord(rd.brec);
	// 	brec = &brec_tmp;
	// }

	~CReadAln() { 
		// GMessage("CReadAln destructor called\n");
		// delete brec;
	}
};

// struct GReadAlnData {
// 	GSamRecord* brec;
// 	char strand; //-1, 0, 1
// 	int nh;
// 	int hi;
// 	// GPVec<CJunction> juncs;
// 	//GPVec< GVec<RC_ExonOvl> > g_exonovls; //>5bp overlaps with guide exons, for each read "exon"
// 	GReadAlnData(GSamRecord* bamrec=NULL, char nstrand=0, int num_hits=0,
// 			int hit_idx=0):brec(bamrec), strand(nstrand),
// 					nh(num_hits), hi(hit_idx) { } //, g_exonovls(true)
// 	~GReadAlnData() { }
// };

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
	// GList<CReadAln> readlist;
	//  GList<CJunction> junction;
	BundleData():status(BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0),
			numreads(0),
			num_fragments(0), frag_len(0),sum_cov(0),covflags(0),
			refseq() {
				// , //bpcov(1024), junction(true, true, true) {
	}

	void getReady(int currentstart, int currentend) {
		//this is only called when the bundle is valid and ready to be processed
		start=currentstart;
		end=currentend;
		//refseq=ref;
		//tag all these guides
		status=BUNDLE_STATUS_READY;
	}

	//bool evalReadAln(GSamRecord& brec, char& strand, int nh); //, int hi);
	//  bool evalReadAln(GReadAlnData& alndata, char& strand);

	void Clear() {
		start=0;
		end=0;
		status=BUNDLE_STATUS_CLEAR;
		numreads=0;
		num_fragments=0;
		frag_len=0;
		sum_cov=0;
		covflags=0;
	}

	~BundleData() {
		Clear();
	}
};


#endif