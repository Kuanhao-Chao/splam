/*  bundle.h -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */
#ifndef _BUNDLE_H_
#define _BUNDLE_H_

#include "gclib/GStr.h"
#include "gclib/GList.hh"
#include "GSam.h"

struct CReadAln:public GSeg {
	GSamRecord brec;
	int pair_idx;     // keeps index for the pair of the alignment.
	CReadAln(GSamRecord* bamrec): 
			GSeg(bamrec->start, bamrec->end), brec(NULL), pair_idx(-1) {
		this->brec = (*bamrec);
	}

	~CReadAln() { 
	}
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
	int idx; //index in the main bundles array
	int start;
	int end;
	int numreads;
	GStr refseq; //reference sequence name

	BundleData():status(BUNDLE_STATUS_CLEAR), idx(0), start(0), end(0), numreads(0), refseq() {
	}

	void getReady(int currentstart, int currentend) {
		//this is only called when the bundle is valid and ready to be processed
		start=currentstart;
		end=currentend;
		//tag all these guides
		status=BUNDLE_STATUS_READY;
	}

	void Clear() {
		start=0;
		end=0;
		status=BUNDLE_STATUS_CLEAR;
		numreads = 0;
	}

	~BundleData() {
		Clear();
	}
};


#endif