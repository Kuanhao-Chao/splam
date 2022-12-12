#include <torch/torch.h>
#include <torch/script.h>
#include <iostream>
#include <memory>
#include <gclib/GArgs.h>
#include <gclib/GBase.h>
#include <gclib/GStr.h>
#include "tmerge.h"
#include "util.h"

#define VERSION "0.0.1"

using namespace std;

const char* USAGE = "SPLAM v" VERSION "\n"
                              "==========================================================================================\n"
                              "This is a file to clean up BAM file.\n"
                              "==========================================================================================\n";

GStr outfdir;
GStr outfname_discard;
GStr outfname_cleaned;

GStr guidegff; // "G" tag
GSamWriter* outfile_discard=NULL;
GSamWriter* outfile_cleaned=NULL;

GStr modelname;
GStr junctionfname;
bool verbose = false;
TInputFiles inRecords;
TInputRecord* irec=NULL;
GSamRecord* brec=NULL;
float threshold = 0.3;

void processOptions(int argc, char* argv[]);

struct intron_key {
	string seqname;
	char strand;
	uint start;
	uint end;

	// intron_key();
	inline bool operator==(const intron_key& a) const {
		if (a.seqname.compare(seqname)==0 && a.strand==strand && a.start==start && a.end==end)
			return true;
		else
			return false;
	}

	bool operator<(const intron_key& a) const
    {
        return (this->start < a.start);
    }
};


int main(int argc, char* argv[]) {
  inRecords.setup(VERSION, argc, argv);
	processOptions(argc, argv);
  int numSamples=inRecords.start();
	cout << "numSamples: " << numSamples << endl;
	cout << "outfdir  : " << outfdir << endl;

	outfname_discard = outfdir + "/discard.bam";
	outfname_cleaned = outfdir + "/cleaned.bam";

	outfile_discard = new GSamWriter(outfname_discard, inRecords.header(), GSamFile_BAM);
	outfile_cleaned = new GSamWriter(outfname_cleaned, inRecords.header(), GSamFile_BAM);

  // torch::Tensor tensor = torch::rand({2, 3});
  // std::cout << tensor << std::endl;
  const char *banner = R"""(
  ███████╗██████╗ ██╗      █████╗ ███╗   ███╗██╗
  ██╔════╝██╔══██╗██║     ██╔══██╗████╗ ████║██║
  ███████╗██████╔╝██║     ███████║██╔████╔██║██║
  ╚════██║██╔═══╝ ██║     ██╔══██║██║╚██╔╝██║╚═╝
  ███████║██║     ███████╗██║  ██║██║ ╚═╝ ██║██╗
  ╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝
  )""";
  std::cout << banner << std::endl;
  
  /************************
   * Loading the model
   ************************/
  torch::jit::script::Module module;
  module = torch::jit::load(modelname.chars());

  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
    std::cout << "Loading "<< modelname.chars() << std::endl;
    module = torch::jit::load(modelname.chars());
  }
  catch (const c10::Error& e) {
    std::cerr << "error loading the model\n";
    return -1;
  }
  std::cout << "Model "<< modelname.chars() <<" loaded fine\n";

  /************************
   * Processing JUNCTION file.
   ************************/
  set<intron_key> bad_intron_set;
  set<intron_key> good_intron_set;
  string line;
  ifstream myfile (junctionfname); // this is equivalent to the above method
  if ( myfile.is_open() ) { // always check whether the file is open

      // chr18	21682058	21596619	JUNC_144610	0	+	8.2311524e-10	4.0248174e-12
      while ( std::getline(myfile, line) ) {

        int parser_counter = 0;

        string chr = "";
        int start = 0;
        int end = 0;
        string junc_name = "";
        int tmp = 0;
        char strand = ' ';
        float d_score = .0;
        float a_score = .0;

        string token;
        istringstream ss(line);

        // then read each element by delimiter
        while (std::getline(ss, token, '\t')) {
            // std::cout << token << std::endl;
            if (parser_counter == 0) chr = token;
            else if (parser_counter == 1) start = std::stoi(token);
            else if (parser_counter == 2) end = std::stoi(token);
            else if (parser_counter == 3) junc_name = token;
            else if (parser_counter == 4) tmp = std::stoi(token);
            else if (parser_counter == 5) strand = token[0];
            else if (parser_counter == 6) d_score = stof(token);
            else if (parser_counter == 7) a_score = stof(token);
            parser_counter += 1;
        }

        // cout << "chr: " << chr << endl;
        // cout << "start: " << start << endl;
        // cout << "end: " << end << endl;
        // cout << "junc_name: " << junc_name << endl;
        // cout << "tmp: " << tmp << endl;
        // cout << "strand: " << strand << endl;
        // cout << "d_score: " << d_score << endl;
        // cout << "a_score: " << a_score << endl;

        intron_key* intronkey = new intron_key;
        intronkey->seqname = chr;
        intronkey->strand = strand;
        intronkey->start = start;
        intronkey->end = end;

        if (d_score > threshold && a_score > threshold) {
          good_intron_set.insert(*intronkey);
        } else {
  			  bad_intron_set.insert(*intronkey);
        }
      }
      myfile.close();
  }
  cout << "bad_intron_set.size(): " << bad_intron_set.size() << endl;
  cout << "good_intron_set.size(): " << good_intron_set.size() << endl;




  /************************
   * Processing BAM file.
   ************************/
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
				// int strand = 0;
				// if (brec->spliceStrand() == '+') strand = 1;
				// if (brec->spliceStrand() == '-') strand = -1;
				// cout << "brec->refName(): " << brec->refName()<< endl;
				string bamseq_name(brec->refName());
				// cout << "bamseq_name: " << bamseq_name << endl;

				intron_key* bamkey = new intron_key;
				bamkey->seqname = bamseq_name;
				bamkey->strand = brec->spliceStrand();
				bamkey->start = brec->exons[i-1].end+1;
				bamkey->end = brec->exons[i].start-1;

				// // intron_key bamkey {
				// // 	bamseq_name,
				// // 	strand,
				// // 	brec->exons[i-1].end+1,
				// // 	brec->exons[i].start-1,
				// // };

				if (bad_intron_set.count(*bamkey)) {
				// if (GFF_hash.hasKey(bamkey)) {
					// BAM_hash.Add(bamkey);
					// BAM_set.insert(*bamkey);
					counter += 1;
					cout << "Splam!  " << counter << endl;
	        outfile_discard->write(brec);
					// cout << "brec->spliceStrand()   : " << brec->spliceStrand() << endl;
					// cout << "brec->refName(),       : " << brec->refName() << endl;
					// cout << "brec->exons[i-1].end+1 : " << brec->exons[i-1].end+1 << endl;
					// cout << "brec->exons[i].start-1 : " << brec->exons[i].start-1 << endl;

				} else {
					// BAM_unmapped_set.insert(*bamkey);
          outfile_cleaned->write(brec);
				}
				// cout << "\tIntron: " << brec->refName() << "; " << brec->spliceStrand() << "; " << brec->exons[i-1].end+1 << " - " << brec->exons[i].start-1 << endl;	
				// CJunc j(brec->exons[i-1].end+1, brec->exons[i].start-1, strand,
				// 		dupcount);
			}
		} else {
      // outfile_cleaned->write(brec);
    }
	}
  delete outfile_discard;
  delete outfile_cleaned;

  return 0;
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;SLPEDVho:N:Q:F:M:J:");
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
    outfdir=args.getOpt('o');
    if (outfdir.is_empty()) {
        GMessage(USAGE);
        GMessage("\nError: output filename must be provided (-o)!\n");
        exit(1);
    }

    if (args.getOpt('J')) {
      junctionfname=args.getOpt('J');
      if (fileExists(junctionfname.chars())>1) {
        // guided=true;
      } else {
        GError("Error: model file (%s) not found.\n",
            junctionfname.chars());
      }
    } else {
          GMessage(USAGE);
          GMessage("\nError: junction fa file must be provided (-J)!\n");
          exit(1);
    }

    if (args.getOpt('M')) {
      modelname=args.getOpt('M');
      if (fileExists(modelname.chars())>1) {
        // guided=true;
      } else {
        GError("Error: model file (%s) not found.\n",
            modelname.chars());
      }
    } else {
          GMessage(USAGE);
          GMessage("\nError: model file must be provided (-M)!\n");
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
