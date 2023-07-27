/*  tmerge.h

    Copyright (C) 2016 Mihaela Pertea & Geo Pertea

    Author: Mihaela Pertea <mpertea@jhu.edu>
    Author: Geo Pertea

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef _TMERGE_H_
#define _TMERGE_H_

#include <map>
#include <gclib/GStr.h>
#include <gclib/GVec.hh>
#include <gclib/GList.hh>
#include "GSam.h"
#include <htslib/htslib/khash.h>

struct TSamReader {
	GStr fname;
	GSamReader* samreader;
	bool tbMerged; //based on the header, is the file a product of TieBrush?
	TSamReader(const char* fn=NULL, GSamReader* samr=NULL):
		fname(fn), samreader(samr), tbMerged(false) {}
	~TSamReader() {
		delete samreader;
	}
};

struct TInputRecord {
	GSamRecord* brec;
	int fidx; //file index in files and readers
	bool tbMerged; //is it from a TieBrush generated file?
	bool operator<(TInputRecord& o) {
		 //decreasing location sort
		 GSamRecord& r1=*brec;
		 GSamRecord& r2=*(o.brec);
		 int r1_tid=r1.refId();
		 int r2_tid=r2.refId();
		 if (r1_tid==r2_tid) {
		 //higher coords first
			if (r1.start!=r2.start)
				 return (r1.start>r2.start);
			else {
				if (r1.end!=r2.end)
				   return (r1.end>r2.end);
				else if (fidx==o.fidx)
						return strcmp(r1.name(), r2.name())>0;
					else return fidx>o.fidx;
			}
		 }
		 else {
			 return (r1_tid>r2_tid);
		 }
	}
	bool operator==(TInputRecord& o) {
		 GSamRecord& r1=*brec;
		 GSamRecord& r2=*(o.brec);
		 return ( r1.refId()==r2.refId() && r1.start==r2.start && r1.end==r2.end
				 && fidx==o.fidx && r1.get_b()->l_data==r2.get_b()->l_data &&
				 memcmp(r1.get_b()->data, r1.get_b()->data, r1.get_b()->l_data)==0
				 );
	}
    void disown() {
    	brec=NULL;
    }
	TInputRecord(GSamRecord* b=NULL, int i=0, bool tb_merged=false):brec(b),
			fidx(i),tbMerged(tb_merged) {}
	~TInputRecord() {
		delete brec;
	}
};

struct TInputFiles {
 protected:
	TInputRecord* crec;
	// use that to check if each input SAM file has the refseqs sorted by coordinate
	sam_hdr_t* mHdr; //merged output header data
	char* pg_ver;
	GStr pg_args;
 public:
	GPVec<TSamReader> freaders;
	void addFile(const char* fn);
	bool addSam(GSamReader* r, int fidx); //update mHdr data
	GList<TInputRecord> recs; //next record for each
	TInputFiles():crec(NULL), mHdr(NULL), pg_ver(NULL), pg_args(),
			freaders(true), recs(true, true, true) { }

	sam_hdr_t* header() { return mHdr; }

	void setup(const char* ver, int argc, char** argv) {
		if (ver) pg_ver=Gstrdup(ver);
		if (argc>0 && argv!=NULL) {
			 for (int i=0;i<argc;i++) {
			   pg_args.append(argv[i]);
			   if (i<argc-1) pg_args.append(' ');
			 }
		}
	}

	~TInputFiles() {
		GFREE(pg_ver);
		sam_hdr_destroy(mHdr);
	}

	int count() { return freaders.Count(); }
	int start(); //open all files, load 1 record from each
	TInputRecord* next();
	void stop(); //

	// index declarations
    bool add_tb_tag_if_not_exists(sam_hdr_t *bh); // adds a line to the header which tells whether the file has been processed with tiebrush before
    void delete_all_hdr_with_tag(sam_hdr_t *hdr,std::string tag1, std::string tag2);
    std::string get_full_path(std::string fname);
    void load_hdr_samples(sam_hdr_t* hdr,std::string filename,bool tbMerged,bool donor); // returns true if ID:SAMPLE present
    bool get_sample_from_line(std::string& line);
	std::string headerfilename; // filename of the file which was used to construct the header
	bool headerfiletbMerged; // whether the input file from which header was borrowed was processed by tiebrush
	int max_sample_id = 0; // current line number of the last sample in the merged header
	std::map<std::string,std::tuple<int,int,std::string,bool>> sample2lineno; // value: first int is the line number; second int is the linenumber in the input file to which sample correspods; third string is the filename of the corresponding index; fourth bool is true if the sample is the one which donated the header
	std::pair<std::map<std::string,std::tuple<int,int,std::string,bool>>::iterator,bool> s2l_it; // check for no duplicate samples
	std::map<int,std::tuple<std::string,int,std::string,bool>> lineno2sample; // value: first string is the sample name; second int is the linenumber in the input file to which sample correspods; third string is the filename of the corresponding index; fourth bool is true if the sample is the one which donated the header
};


#endif /* TIEBRUSH_TMERGE_H_ */
