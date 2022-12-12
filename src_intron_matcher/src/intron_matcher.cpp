#include <iostream>
#include <gclib/GArgs.h>
#include <gclib/gff.h>
#include "tmerge.h"
#include "commons.h"
#include <map>
#include <set>
#define VERSION "0.0.1"

using namespace std;

TInputFiles inRecords;

GStr outfname;
GStr guidegff; // "G" tag
GSamWriter* outfile=NULL;
bool verbose=false;

GStrSet<> excludeGseqs; // "x" tag. hash of chromosomes/contigs to exclude (e.g. chrM)
bool skipGseq=false;

void processOptions(int argc, char* argv[]);

struct intron_key {
	char* seqname;
	int strand;
	uint start;
	uint end;

	// intron_key();
	inline bool operator==(const intron_key& a) const {
		if (strcmp(a.seqname, seqname)==0 && a.strand==strand && a.start==start && a.end==end)
			return true;
		else
			return false;
	}

	bool operator<(const intron_key& a) const
    {
        return (this->start < a.start);
    }
};

const char* USAGE = "Intron_matcher v" VERSION "\n"
                              "==========================================================================================\n"
                              "This is a simple script to compare the intron similarity between BAM and GFF annotation.\n"
                              "==========================================================================================\n";

int main(int argc, char* argv[]) {
	int total_intron_num = 0;
	int matched_intron_num = 0;
	int unmatched_intron_num = 0;

	const char *banner = R"""(
	██╗███╗   ██╗████████╗██████╗  ██████╗ ███╗   ██╗        ███╗   ███╗ █████╗ ████████╗ ██████╗██╗  ██╗███████╗██████╗ 
	██║████╗  ██║╚══██╔══╝██╔══██╗██╔═══██╗████╗  ██║        ████╗ ████║██╔══██╗╚══██╔══╝██╔════╝██║  ██║██╔════╝██╔══██╗
	██║██╔██╗ ██║   ██║   ██████╔╝██║   ██║██╔██╗ ██║        ██╔████╔██║███████║   ██║   ██║     ███████║█████╗  ██████╔╝
	██║██║╚██╗██║   ██║   ██╔══██╗██║   ██║██║╚██╗██║        ██║╚██╔╝██║██╔══██║   ██║   ██║     ██╔══██║██╔══╝  ██╔══██╗
	██║██║ ╚████║   ██║   ██║  ██║╚██████╔╝██║ ╚████║███████╗██║ ╚═╝ ██║██║  ██║   ██║   ╚██████╗██║  ██║███████╗██║  ██║
	╚═╝╚═╝  ╚═══╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝╚═╝     ╚═╝╚═╝  ╚═╝   ╚═╝    ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝
	)""";
	cout << banner << endl;

	inRecords.setup(VERSION, argc, argv);
	processOptions(argc, argv);
	int numSamples=inRecords.start();
	cout << "numSamples: " << numSamples << endl;
	cout << "outfname  : " << outfname << endl;
	outfile = new GSamWriter(outfname, inRecords.header(), GSamFile_BAM);
	TInputRecord* irec=NULL;
	GSamRecord* brec=NULL;

	/***************************
	 * Reading GFF file.
	 ***************************/
	FILE* f=fopen(guidegff.chars(),"r");
	if (f==NULL) GError("Error: could not open reference annotation file (%s)!\n",
		guidegff.chars());
	//                transcripts_only    sort by location?
	GffReader gffr(f,       true,             true); //loading only recognizable transcript features
	gffr.setRefAlphaSorted(); //alphabetical sorting of refseq IDs
	gffr.showWarnings(verbose);
	//        keepAttrs    mergeCloseExons   noExonAttrs
	gffr.readAll(false,          true,        true);

	set<intron_key> GFF_set;
	set<intron_key> BAM_set;
	set<intron_key> BAM_unmapped_set;

	char* last_ref_seq_name = "";
	for (int i=0;i<gffr.gflst.Count();i++) {
		GffObj* m=gffr.gflst[i];
		char* gffseq_name = strdup(m->getGSeqName());

		int strand = 0;
		if (m->strand == '+') strand = 1;
		if (m->strand == '-') strand = -1;

		for (int j=1; j < m->exons.Count(); j++) {

		   intron_key* intronkey = new intron_key;
		   intronkey->seqname = gffseq_name;
		   intronkey->strand = strand;
		   intronkey->start = m->exons[j-1]->end+1;
		   intronkey->end = m->exons[j]->start-1;

			// intron_key intronkey {
			// 	gffseq_name,
			// 	strand,
			// 	m->exons[j-1]->end+1,
			// 	m->exons[j]->start-1,
			// };

			// cout << "\tGFF Intron: " << m->getGSeqName() << "; " << m->strand << "; " << m->exons[j-1]->end+1 << " - " << m->exons[j]->start-1 << endl;	
			// cout << "\tGFF Exon: " << m->exons[j]->start << " - " << m->exons[j]->end << endl;
			// int k = GFF_hash.Add(*intronkey);
			GFF_set.insert(*intronkey);

			// cout << "\tintronkey->seqname: " << intronkey->seqname << endl;
			// cout << "\tintronkey->strand : " << intronkey->strand << endl;
			// cout << "\tintronkey->start  : " << intronkey->start << endl;
			// cout << "\tintronkey->end    : " << intronkey->end << endl;
			// cout << "k: " << k << endl;
			// cout << "GFF_hash.size(): " << GFF_hash.size() << endl;
		}
	}

	cout << "GFF_set size: " << GFF_set.size() << endl;
	// for(auto& key: GFF_set)
	// {
	// 	std::cout << "\t" << key.seqname << endl;
	// 	std::cout << "\t" << key.strand << endl;
	// 	std::cout << "\t" << key.start << endl;
	// 	std::cout << "\t" << key.end << endl;
	// }

	/***************************
	 * Reading BAM file.
	 ***************************/
	int counter = 0;
	while ((irec=inRecords.next())!=NULL) {
		brec=irec->brec;
		// cout << irec->fidx << endl;
		if (brec->hasIntrons()) {
			// This is a spliced read => start processing it!
			// char strand = brec->spliceStrand();
			// cout << "strand       : " << strand << endl;
			// cout << "brec->cigar(): " << brec->cigar() << endl;
			for (int i=1;i<brec->exons.Count();i++) {
				int strand = 0;
				if (brec->spliceStrand() == '+') strand = 1;
				if (brec->spliceStrand() == '-') strand = -1;
				// cout << "brec->refName(): " << brec->refName()<< endl;
				char* bamseq_name = strdup(brec->refName());
				// cout << "bamseq_name: " << bamseq_name << endl;

				intron_key* bamkey = new intron_key;
				bamkey->seqname = bamseq_name;
				bamkey->strand = strand;
				bamkey->start = brec->exons[i-1].end+1;
				bamkey->end = brec->exons[i].start-1;

				// intron_key bamkey {
				// 	bamseq_name,
				// 	strand,
				// 	brec->exons[i-1].end+1,
				// 	brec->exons[i].start-1,
				// };

				if (GFF_set.count(*bamkey)) {
				// if (GFF_hash.hasKey(bamkey)) {
					// BAM_hash.Add(bamkey);
					BAM_set.insert(*bamkey);
					counter += 1;
					// cout << "Splam!  " << counter << endl;
					// cout << "brec->spliceStrand()   : " << brec->spliceStrand() << endl;
					// cout << "brec->refName(),       : " << brec->refName() << endl;
					// cout << "brec->exons[i-1].end+1 : " << brec->exons[i-1].end+1 << endl;
					// cout << "brec->exons[i].start-1 : " << brec->exons[i].start-1 << endl;

				} else {
					BAM_unmapped_set.insert(*bamkey);
				}
				// cout << "\tIntron: " << brec->refName() << "; " << brec->spliceStrand() << "; " << brec->exons[i-1].end+1 << " - " << brec->exons[i].start-1 << endl;	
				// CJunc j(brec->exons[i-1].end+1, brec->exons[i].start-1, strand,
				// 		dupcount);
			}
		}
	}

	total_intron_num = GFF_set.size();
	matched_intron_num = BAM_set.size();
	unmatched_intron_num = BAM_unmapped_set.size();
	
	cout << "GFF_hash.size()             : " << total_intron_num << endl;
	cout << "BAM_hash.size()             : " << matched_intron_num << endl;
	cout << "BAM_hash_unmapped.size()    : " << unmatched_intron_num << endl;

	cout << "\n\nIntron matching precision   : " << (float)matched_intron_num / (float)(matched_intron_num + unmatched_intron_num) << endl;
	cout << "Intron matching sensitivity : " << (float)matched_intron_num / (float)(total_intron_num) << endl;

	inRecords.stop();

  return 1;
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;full;clip;exon;keep-supp;keep-unmap;SMLPEDVho:N:Q:F:G:");
    args.printError(USAGE, true);

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s",USAGE);
        exit(0);
    }

    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }

    if (args.startNonOpt()==0) {
        GMessage(USAGE);
        GMessage("\nError: no input provided!\n");
        exit(1);
    }
    outfname=args.getOpt('o');
    if (outfname.is_empty()) {
        GMessage(USAGE);
        GMessage("\nError: output filename must be provided (-o)!\n");
        exit(1);
    }

	if (args.getOpt('G')) {
		guidegff=args.getOpt('G');
		if (fileExists(guidegff.chars())>1) {
			// guided=true;
		} else {
			GError("Error: reference annotation file (%s) not found.\n",
					guidegff.chars());
		}
	} else {
        GMessage(USAGE);
        GMessage("\nError: gff reference file must be provided (-G)!\n");
        exit(1);
	}

    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running intron_matcher " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    const char* ifn=NULL;
    while ( (ifn=args.nextNonOpt())!=NULL) {
        //input alignment files
        std::string absolute_ifn = get_full_path(ifn);
        cout << "absolute_ifn: " << absolute_ifn << endl;
        inRecords.addFile(absolute_ifn.c_str());
    }
}
