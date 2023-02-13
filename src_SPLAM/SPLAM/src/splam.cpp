// #include <torch/torch.h>
// #include <torch/script.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <unordered_map>
// #include <sys/stat.h>

#include <gclib/GArgs.h>
#include <gclib/GBase.h>
#include <gclib/GStr.h>
#include "tmerge.h"
#include "util.h"
#include "junc.h"
#include "filter.h"

#include <htslib/htslib/faidx.h>
#include <Python.h>


#define VERSION "0.0.1"

using namespace std;

const char* USAGE = "SPLAM v" VERSION "\n"
                              "==========================================================================================\n"
                              "An accurate spliced alignment pruner and spliced junction predictor.\n"
                              "==========================================================================================\n";

GStr outfdir;
GStr outfname_discard;
GStr outfname_spliced;
GStr outfname_cleaned;

GStr guidegff; // "G" tag
GSamWriter* outfile_discard = NULL;
GSamWriter* outfile_spliced = NULL;
GSamWriter* outfile_cleaned = NULL;

GStr modelname;
GStr junctionfname;
FILE* joutf=NULL;

bool verbose = false;
TInputFiles inRecords;
TInputRecord* irec=NULL;
GSamRecord* brec=NULL;
float threshold = 0.3;

void processOptions(int argc, char* argv[]);
void loadBed(GStr inbedname, GArray<CJunc> &spur_juncs, unordered_map<string, string> &chrs_refseq_2_ucsc);
void flushBrec(GVec<GSamRecord*> &pbrecs, unordered_map<string, int> &hits, GSamWriter* outfile_cleaned);

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

unordered_map<string, int> get_hg38_chrom_size(string target) {
    unordered_map<string, int> chrs;
    string ref_file;
    if (target == "STAR") {
        ref_file = "../../../src/hg38_chrom_size_refseq.tsv";
        // f_chrs = open("../hg38_chrom_size_refseq.tsv", "r")
    } else {
        ref_file = "../../../src/hg38_chrom_size.tsv";
    }
    cout << "ref_file: " << ref_file << endl;
    std::ifstream ref_f(ref_file);


    string line;
    while(getline(ref_f, line)){
        // cout << line << endl;
        string chromosome;
        int len;

        std::replace(line.begin(), line.end(), '\t', ' ');

        stringstream ss(line);

        ss >> chromosome;
        ss >> len;

        // cout << "chromosome: " << chromosome << " ";
        // cout << "len: " << len << " " << endl;
        chrs[chromosome] = len;
    }   

    // for (auto i : chrs) {
    //   cout << i.first << " ---- " << i.second << endl;
    // }
    return chrs;
}

void get_Refseq_2_UCSC_chr_names(unordered_map<string, string> &chrs_refseq_2_ucsc, unordered_map<string, string> &chrs_ucsc_2_refseq) {

    string ref_file;
    ref_file = "../../../src/Refseq_2_UCSU_chromosome_names.tsv";
    cout << "ref_file: " << ref_file << endl;
    std::ifstream ref_f(ref_file);


    string line;
    while(getline(ref_f, line)){
        // cout << line << endl;
        string refseq;
        string ucsc;

        std::replace(line.begin(), line.end(), '\t', ' ');

        stringstream ss(line);

        ss >> refseq;
        ss >> ucsc;

        // cout << "chromosome: " << chromosome << " ";
        // cout << "len: " << len << " " << endl;
        chrs_ucsc_2_refseq[ucsc] = refseq;
        chrs_refseq_2_ucsc[refseq] = ucsc;
    }   

    for (auto i : chrs_refseq_2_ucsc) {
      cout << i.first << " ---- " << i.second << endl;
    }

    for (auto i : chrs_ucsc_2_refseq) {
      cout << i.first << " ---- " << i.second << endl;
    }
    // return chrs_refseq_2_ucsc;
}

int main(int argc, char* argv[]) {
    const char *banner = R"""(
  ███████╗██████╗ ██╗      █████╗ ███╗   ███╗    ██╗
  ██╔════╝██╔══██╗██║     ██╔══██╗████╗ ████║    ██║
  ███████╗██████╔╝██║     ███████║██╔████╔██║    ██║
  ╚════██║██╔═══╝ ██║     ██╔══██║██║╚██╔╝██║    ╚═╝
  ███████║██║     ███████╗██║  ██║██║ ╚═╝ ██║    ██╗
  ╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝    ╚═╝
    )""";
    std::cout << banner << std::endl;
    
    
    inRecords.setup(VERSION, argc, argv);
    processOptions(argc, argv);
    int numSamples=inRecords.start();
    cout << "numSamples: " << numSamples << endl;
    cout << "outfdir        : " << outfdir << endl;
    cout << "junctionfname  : " << junctionfname << endl;


    mkdir(outfdir, 0777);
    mkdir(outfdir+"/bed", 0777);
    mkdir(outfdir+"/fasta", 0777);


    outfname_discard = outfdir + "/discard.bam";
    outfname_spliced = outfdir + "/spliced.bam";
    outfname_cleaned = outfdir + "/cleaned.bam";

    outfile_discard = new GSamWriter(outfname_discard, inRecords.header(), GSamFile_BAM);
    outfile_spliced = new GSamWriter(outfname_spliced, inRecords.header(), GSamFile_BAM);
    outfile_cleaned = new GSamWriter(outfname_cleaned, inRecords.header(), GSamFile_BAM);


    unordered_map<string, int> chrs = get_hg38_chrom_size("HISAT2");

    unordered_map<string, string> chrs_refseq_2_ucsc;
    unordered_map<string, string> chrs_ucsc_2_refseq;
    get_Refseq_2_UCSC_chr_names(chrs_refseq_2_ucsc, chrs_ucsc_2_refseq);

    string fa_name = "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic.fna";

    FILE *ref_fa_f;
    if (ref_fa_f = fopen(fa_name.c_str(), "r")) {
        fclose(ref_fa_f);
        cout << fa_name << "i has been created!" << endl;
    } else {
        int res = fai_build(fa_name.c_str());
        cout << "Creating " << fa_name << "i" << endl;
        cout << ">> res: " << res << endl;
    }

    faidx_t * ref_faidx = fai_load(fa_name.c_str());
    cout << ">> ref_faidx: " << ref_faidx << endl;



    /*********************************************
     * Step 1: generating spliced junctions in BED
    *********************************************/
    cout << "********************************************" << endl;
    cout << "** Step 1: generating spliced junctions in BED" << endl;
    cout << "********************************************" << endl;

    /***************************
     * Creating the output junction bed file
    ***************************/    
    if (!junctionfname.is_empty()) {
        if (strcmp(junctionfname.substr(junctionfname.length()-4, 4).chars(), ".bed")!=0) {
            junctionfname.append(".bed");
        }
        joutf = fopen(junctionfname.chars(), "w");
        if (joutf==NULL) GError("Error creating file %s\n", junctionfname.chars());
        fprintf(joutf, "track name=junctions\n");
    }

    /***************************
     * Reading BAM file.
     ***************************/
    int counter = 0;
    int prev_tid=-1;
    // char* prev_refname;
    string prev_refname;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    //    std::vector<std::set<int>> bsam_idx(2048*1024,std::set<int>{}); // for indexed runs
    int b_end=0; //bundle start, end (1-based)
    int b_start=0; //1 based

    while ((irec=inRecords.next())!=NULL) {
        brec=irec->brec;
        // cout << brec->refId() << endl;
        uint32_t dupcount=0;
        std::vector<int> cur_samples;
        int endpos=brec->end;

        if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
            // if (joutf) {
            //     flushJuncs(joutf, prev_refname);
            // } // TODO: write the last column to 3 dec places
            b_start=brec->start;
            b_end=endpos;
            prev_tid=brec->refId();

            prev_refname=(char*)brec->refName();
        } else { //extending current bundle
            if (b_end<endpos) {
                b_end=endpos;
                bcov.setCount(b_end-b_start+1, (int)0);
            }
        }
        int accYC = 0;
        accYC = brec->tag_int("YC", 1);
        // cout << "accYC: " << accYC << endl;
        if (joutf && brec->exons.Count()>1) {
            // cout << "1 prev_refname: " << prev_refname << endl;
            // prev_refname = chrs_convert[prev_refname];
            // cout << "2 prev_refname: " << prev_refname << endl;
            addJunction(*brec, accYC, prev_refname);
            outfile_spliced->write(brec);
        } else {
            outfile_cleaned->write(brec);
        }
    }
    writeJuncs(joutf);
    fclose(joutf);

    // The reference is in refseq chr name.
    /*********************************************
     * Step 2: (1) getting coordinates of donors and acceptors
     *         (2) Writing FASTA file of donors and acceptors
     *         (3) checking the junctions. (GT-AG ratio)
    *********************************************/
    cout << "********************************************" << endl;
    cout << "** Step 2: getting coordinates of donors and acceptors" << endl;
    cout << "********************************************" << endl;

    string threshold = "100";
    string SEQ_LEN="800";
    int QUOTER_SEQ_LEN = stoi(SEQ_LEN)/4; // 4

    /*********************************************
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    *******************************************/
    GStr donor_bed = outfdir + "/bed/donor.bed";
    GStr acceptor_bed = outfdir + "/bed/acceptor.bed";

    ofstream fw_donor(donor_bed);
    ofstream fw_acceptor(acceptor_bed);
    ofstream fw_da(outfdir + "/bed/d_a.bed");

    GStr junc_fasta = outfdir + "/fasta/junction.fa";
    GStr donor_fasta = outfdir + "/fasta/donor.fa";
    GStr accceptor_fasta = outfdir + "/fasta/acceptor.fa";

    ofstream fw_fa_junc(junc_fasta);
    ofstream fw_fa_donor(donor_fasta);
    ofstream fw_fa_acceptor(accceptor_fasta);


    ifstream fr_junc(junctionfname);
    string line;
    while(getline(fr_junc, line)){
        // cout << line << endl;

        string chromosome;
        int start;
        int end;
        string junc_name;
        int num_alignment;
        string strand;

        std::replace(line.begin(), line.end(), '\t', ' ');

        stringstream ss(line);

        ss >> chromosome;
        ss >> start;
        ss >> end;
        ss >> junc_name;
        ss >> num_alignment;
        ss >> strand;

        // cout << "** chromosome: " << chromosome << " ";
        // cout << "start: " << start << " ";
        // cout << "end: " << end << " ";
        // cout << "junc_name: " << junc_name << " ";
        // cout << "num_alignment: " << num_alignment << " ";
        // cout << "strand: " << strand << " " << endl;

        int splice_junc_len = 0;
        int flanking_size = QUOTER_SEQ_LEN;

        int donor = 0;
        int acceptor = 0;

        int donor_s = 0;
        int donor_e = 0;

        int acceptor_s = 0;
        int acceptor_e = 0;

        splice_junc_len = end - start;
        if (splice_junc_len < QUOTER_SEQ_LEN) {
            flanking_size = splice_junc_len;
        }

        // cout << "QUOTER_SEQ_LEN: " << QUOTER_SEQ_LEN << endl;
        // cout << "flanking_size : " << flanking_size << endl;
        // cout << endl << endl;

        if (strand == "+") {
            donor = start;
            acceptor = end;

            donor_s = donor - QUOTER_SEQ_LEN;
            donor_e = donor + flanking_size;
            acceptor_s = acceptor - flanking_size;
            acceptor_e = acceptor + QUOTER_SEQ_LEN;
        } else if (strand == "-") {
            donor = end;
            acceptor = start;

            donor_s = donor - flanking_size;
            donor_e = donor + QUOTER_SEQ_LEN;
            acceptor_s = acceptor - QUOTER_SEQ_LEN;
            acceptor_e = acceptor + flanking_size;

        } else if (strand == ".") {
            continue;
        }

        if (donor_e >= chrs[chromosome] or acceptor_e >= chrs[chromosome]) {
            cout << "Skip!!" << endl;
            continue;
        }
        if (donor_s < 0 or acceptor_s < 0) {
            cout << "Skip!!" << endl;
            continue;
        }


        fw_donor << chromosome << "\t" + to_string(donor_s) + "\t" + to_string(donor_e) + "\t" + junc_name+"_donor" + "\t" + to_string(num_alignment) + "\t" + strand + "\n";

        fw_acceptor << chromosome << "\t" + to_string(acceptor_s) + "\t" + to_string(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + to_string(num_alignment) + "\t" + strand + "\n";


        cout << ">>> Before chromosome: " << chromosome << endl;
        chromosome = chrs_ucsc_2_refseq[chromosome];
        cout << ">>> After chromosome: " << chromosome << endl;
        cout << "strand: " << strand << endl;

        int donor_len = donor_e - donor_s;
        char* donor_seq = faidx_fetch_seq(ref_faidx, chromosome.c_str(), donor_s, donor_e-1, &donor_len);

        int acceptor_len = acceptor_e - acceptor_s;
        char* acceptor_seq = faidx_fetch_seq(ref_faidx, chromosome.c_str(), acceptor_s, acceptor_e-1, &acceptor_len);

        if (strand == "+") {

        } else if (strand == "-") {
            hts_pos_t donor_len_hts = donor_e - donor_s;
            hts_pos_t acceptor_len_hts = acceptor_e - acceptor_s;

            reverse_complement(donor_seq, donor_len_hts);
            reverse_complement(acceptor_seq, acceptor_len_hts);
        }

        if ((donor_seq == NULL) ){
            printf("c is empty\n");
            continue;
        }

        if ((acceptor_seq == NULL) ) {
            printf("c is empty\n");
            continue;
        }

        fw_fa_donor << ">" << chromosome << endl;
        fw_fa_donor << donor_seq << endl;
        fw_fa_acceptor << ">" << chromosome << endl;
        fw_fa_acceptor << acceptor_seq << endl;

        cout << "strlen(donor_seq)   : " << strlen(donor_seq) << endl;
        cout << "strlen(acceptor_seq): " << strlen(acceptor_seq) << endl;


        fw_fa_junc << ">" << chromosome << ";" << to_string(donor) << ";" << to_string(acceptor) << ";" << strand << endl;
        if (strlen(donor_seq) >= 400) {
            fw_fa_junc << donor_seq << acceptor_seq << endl;
        } else {
            fw_fa_junc << donor_seq << string(2*(400-(int)strlen(donor_seq)), 'N') << acceptor_seq << endl;
        }

        cout << "donor   : " << donor_seq[200] << donor_seq[201] << endl;
        cout << "acceptor: " << acceptor_seq[198] << acceptor_seq[199] << endl;

        if (strand == "+") {
            fw_da << chromosome + "\t" + to_string(donor) + "\t" + to_string(acceptor+1) + "\tJUNC\t" + to_string(num_alignment) + "\t" + strand + "\n";
        } else if (strand == "-") {
            fw_da << chromosome + "\t" + to_string(acceptor) + "\t" + to_string(donor+1) + "\tJUNC\t" + to_string(num_alignment) + "\t" + strand + "\n";
        }
        cout << endl;
    }  






    // /*********************************************
    //  * Step 3: SPLAM model prediction
    // *********************************************/
    // Py_Initialize();
    // PyRun_SimpleString("import sys");

    // string python_f = "../script.py";
    // PyObject *obj = Py_BuildValue("s", python_f.c_str());
    // FILE *file = _Py_fopen_obj(obj, "r+");
    // if(file != NULL) {
    //     PyRun_SimpleFile(file, python_f.c_str());
    // }

    // FILE *fd = fopen(python_f.c_str(), "r");
    // PyRun_SimpleFile(fd, python_f.c_str()); // last parameter == 1 means to close the
    //                                     // file before returning.


    // /*********************************************
    //  * Step 4: SPLAM filtering out reads.
    // *********************************************/
    // // GSamReader bamreader(inbamname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    // // outfile=new GSamWriter(outfname, bamreader.header(), GSamFile_BAM);
    // // GSamRecord brec;


    cout << "********************************************" << endl;
    cout << "** Step 4: SPLAM filtering out reads." << endl;
    cout << "********************************************" << endl;

    auto start=std::chrono::high_resolution_clock::now();
    int spur_cnt = 0;
    std::cout << "brrrm! identifying alignment records with spurious splice junctions" << std::endl;

    GArray<CJunc> spur_juncs;
    GVec<GSamRecord*> kept_brecs;
    std::unordered_map<std::string, int> hits;
    GStr inbedname(outfdir + "/output/junc_scores.bed");
    loadBed(inbedname, spur_juncs, chrs_refseq_2_ucsc);

    cout << ">> (4) Junction count: " << junctions.Count() << endl;
	for (int i = 0; i < junctions.Count(); i++) {
		cout << i <<  " (4) Junction name: " << junctions[i].start << " - " << junctions[i].end << endl;
		cout << ">> (4) Read count: " << junctions[i].read_ls.size() << endl;

        cout << "junctions[i].ref: " << junctions[i].ref << endl;
        cout << "junctions[i].start: " << junctions[i].start << endl;
        cout << "junctions[i].end: " << junctions[i].end << endl;
        cout << "junctions[i].strand: " << junctions[i].strand << endl;
        cout << endl;

        CJunc jnew(junctions[i].start, junctions[i].end, junctions[i].strand, junctions[i].ref);
        
        cout << "spur_juncs.Exists(jnew):  " << spur_juncs.Exists(jnew) << endl;
        if (spur_juncs.Exists(jnew)) {
            // spur = true;
            cout << "spur_juncs.Exists! " << endl;
            for (int j=0; j<junctions[i].read_ls.size(); j++) {
                cout << "~~ SPLAM!" << endl;
                outfile_discard->write(junctions[i].read_ls.at(j));
                // delete junctions[i].read_ls.at(j);
            }
            cout << "spur_juncs.Exists Done! " << endl;
        } else {
            for (int j=0; j<junctions[i].read_ls.size(); j++) {
                bool spur = false;
                int r_exon_count = junctions[i].read_ls.at(j)->exons.Count();
                if (r_exon_count > 1) {
                    for (int e=1; e<r_exon_count; e++) {
                        CJunc jnew_sub(junctions[i].read_ls.at(j)->exons[e-1].end, junctions[i].read_ls.at(j)->exons[e-1].start-1, junctions[i].strand, junctions[i].ref);
                        if (spur_juncs.Exists(jnew_sub)) {
                            spur = true;
                            break;
                        }
                    }
                }



    //         if (!spur) {
    //             // cout << "Not spurious!" << endl;
    //             GSamRecord *rec = new GSamRecord(*brec);
    //             PBRec *newpbr = new PBRec(rec);
    //             kept_brecs.Add(newpbr);
    //         } else {
    //             cout << "Spurious!" << endl;
    //             spur_cnt++;
    //             std::string kv = brec->name();
    //             std::string tmp = std::to_string(brec->pairOrder());
    //             kv += ";";
    //             kv += tmp;
    //             // key not present
    //             if (hits.find(kv) == hits.end()) {
    //                 hits[kv] = 1;
    //             } else {
    //                 int val = hits[kv];
    //                 val++;
    //                 hits[kv] = val;
    //             }
    //         }

                if (spur) {
                    cout << "spur_juncs.Exists! " << endl;
                    cout << "~~ SPLAM!" << endl;
                    std::string kv = brec->name();
                    std::string tmp = std::to_string(brec->pairOrder());
                    kv += ";";
                    kv += tmp;

                    if (hits.find(kv) == hits.end()) {
                        hits[kv] = 1;
                    } else {
                        int val = hits[kv];
                        val++;
                        hits[kv] = val;
                    }

                    outfile_discard->write(junctions[i].read_ls.at(j));
                    // delete junctions[i].read_ls.at(j);
                } else {
                    cout << "spur_juncs not Exists! " << endl;
                    kept_brecs.Add(junctions[i].read_ls.at(j));

                    // cout << "~~ Clean!" << endl;
                    // outfile_cleaned->write(junctions[i].read_ls.at(j));
                    // delete junctions[i].read_ls.at(j);
                }
            }
            cout << "~~~ Done! " << endl;
        }

        cout << "Done!" << endl;

		
        // for (int r = 0; r < junctions[i].read_ls.size(); r++) {
		// 	cout << "\tRead " <<  r << " : " << junctions[i].read_ls[r]->cigar() << endl;
		// }
		// cout << endl;
	}

    flushBrec(kept_brecs, hits, outfile_cleaned);
    auto end =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << spur_cnt << " spurious alignments were removed." << std::endl;
    std::cout << "Completed in " << duration.count() << " seconds" << std::endl;
    


	// while ((irec=inRecords.next())!=NULL) {
	// 	brec=irec->brec;
	// 	// cout << irec->fidx << endl;
	// 	if (brec->hasIntrons()) {
    //         const char* chrname=brec->refName();
    //         char chr = chrname[strlen(chrname) - 1];
    //         char strand = brec->spliceStrand();
    //         bool spur = false;

    //         // cout << "chrname: " << chrname << endl;
    //         for (int i = 1; i < brec->exons.Count(); i++) {
    //             CJunc j(brec->exons[i-1].end, brec->exons[i].start-1, strand, chr);
    //             if (spur_juncs.Exists(j)) {
    //                 spur = true;
    //                 break;
    //             }
    //         }
    //         if (!spur) {
    //             // cout << "Not spurious!" << endl;
    //             GSamRecord *rec = new GSamRecord(*brec);
    //             PBRec *newpbr = new PBRec(rec);
    //             kept_brecs.Add(newpbr);
    //         } else {
    //             cout << "Spurious!" << endl;
    //             spur_cnt++;
    //             std::string kv = brec->name();
    //             std::string tmp = std::to_string(brec->pairOrder());
    //             kv += ";";
    //             kv += tmp;
    //             // key not present
    //             if (hits.find(kv) == hits.end()) {
    //                 hits[kv] = 1;
    //             } else {
    //                 int val = hits[kv];
    //                 val++;
    //                 hits[kv] = val;
    //             }
    //         }
	// 	} else {
    //         GSamRecord *rec = new GSamRecord(*brec);
    //         PBRec* newpbr = new PBRec(rec);
    //         kept_brecs.Add(newpbr);
    //     }
	// }
  
    // std::cout << "vacuuming completed. writing only clean bam records to the output file." << std::endl;
    // flushBrec(kept_brecs, hits, outfile_cleaned);
    // auto end =std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    // std::cout << spur_cnt << " spurious alignment records were removed." << std::endl;
    // std::cout << "Vacuuming completed in " << duration.count() << " seconds" << std::endl;
    
        


    /************************
     * Loading the model
     ************************/
    // cout << ">> Loading the model: " << endl;
    // torch::jit::script::Module module;
    // // module = torch::jit::load(modelname.chars());

    // try {
    //     // Deserialize the ScriptModule from a file using torch::jit::load().
    //     std::cout << "Loading "<< modelname.chars() << std::endl;
    //     module = torch::jit::load(modelname.chars());
    // }
    // catch (const c10::Error& e) {
    //     std::cerr << "error loading the model\n";
    //     return -1;
    // }
    // std::cout << "Model "<< modelname.chars() <<" loaded fine\n";

  // /************************
  //  * Processing JUNCTION file.
  //  ************************/
  // set<intron_key> bad_intron_set;
  // set<intron_key> good_intron_set;
  // string line;
  // ifstream myfile (junctionfname); // this is equivalent to the above method
  // if ( myfile.is_open() ) { // always check whether the file is open

  //     // chr18	21682058	21596619	JUNC_144610	0	+	8.2311524e-10	4.0248174e-12
  //     while ( std::getline(myfile, line) ) {

  //       int parser_counter = 0;

  //       string chr = "";
  //       int start = 0;
  //       int end = 0;
  //       string junc_name = "";
  //       int tmp = 0;
  //       char strand = ' ';
  //       float d_score = .0;
  //       float a_score = .0;

  //       string token;
  //       istringstream ss(line);

  //       // then read each element by delimiter
  //       while (std::getline(ss, token, '\t')) {
  //           // std::cout << token << std::endl;
  //           if (parser_counter == 0) chr = token;
  //           else if (parser_counter == 1) start = std::stoi(token);
  //           else if (parser_counter == 2) end = std::stoi(token);
  //           else if (parser_counter == 3) junc_name = token;
  //           else if (parser_counter == 4) tmp = std::stoi(token);
  //           else if (parser_counter == 5) strand = token[0];
  //           else if (parser_counter == 6) d_score = stof(token);
  //           else if (parser_counter == 7) a_score = stof(token);
  //           parser_counter += 1;
  //       }

  //       // cout << "chr: " << chr << endl;
  //       // cout << "start: " << start << endl;
  //       // cout << "end: " << end << endl;
  //       // cout << "junc_name: " << junc_name << endl;
  //       // cout << "tmp: " << tmp << endl;
  //       // cout << "strand: " << strand << endl;
  //       // cout << "d_score: " << d_score << endl;
  //       // cout << "a_score: " << a_score << endl;

  //       intron_key* intronkey = new intron_key;
  //       intronkey->seqname = chr;
  //       intronkey->strand = strand;
  //       intronkey->start = start;
  //       intronkey->end = end;

  //       if (d_score > threshold && a_score > threshold) {
  //         good_intron_set.insert(*intronkey);
  //       } else {
  // 			  bad_intron_set.insert(*intronkey);
  //       }
  //     }
  //     myfile.close();
  // }
  // cout << "bad_intron_set.size(): " << bad_intron_set.size() << endl;
  // cout << "good_intron_set.size(): " << good_intron_set.size() << endl;




  // /************************
  //  * Processing BAM file.
  //  ************************/
	// int counter = 0;
	// while ((irec=inRecords.next())!=NULL) {
	// 	brec=irec->brec;
	// 	// cout << irec->fidx << endl;
	// 	if (brec->hasIntrons()) {
	// 		// This is a spliced read => start processing it!
	// 		// char strand = brec->spliceStrand();
	// 		// cout << "strand       : " << strand << endl;
	// 		// cout << "brec->cigar(): " << brec->cigar() << endl;
	// 		for (int i=1;i<brec->exons.Count();i++) {
	// 			// int strand = 0;
	// 			// if (brec->spliceStrand() == '+') strand = 1;
	// 			// if (brec->spliceStrand() == '-') strand = -1;
	// 			// cout << "brec->refName(): " << brec->refName()<< endl;
	// 			string bamseq_name(brec->refName());
	// 			// cout << "bamseq_name: " << bamseq_name << endl;

	// 			intron_key* bamkey = new intron_key;
	// 			bamkey->seqname = bamseq_name;
	// 			bamkey->strand = brec->spliceStrand();
	// 			bamkey->start = brec->exons[i-1].end+1;
	// 			bamkey->end = brec->exons[i].start-1;

	// 			// // intron_key bamkey {
	// 			// // 	bamseq_name,
	// 			// // 	strand,
	// 			// // 	brec->exons[i-1].end+1,
	// 			// // 	brec->exons[i].start-1,
	// 			// // };

	// 			if (bad_intron_set.count(*bamkey)) {
	// 			// if (GFF_hash.hasKey(bamkey)) {
	// 				// BAM_hash.Add(bamkey);
	// 				// BAM_set.insert(*bamkey);
	// 				counter += 1;
	// 				cout << "Splam!  " << counter << endl;
	//         outfile_discard->write(brec);
	// 				// cout << "brec->spliceStrand()   : " << brec->spliceStrand() << endl;
	// 				// cout << "brec->refName(),       : " << brec->refName() << endl;
	// 				// cout << "brec->exons[i-1].end+1 : " << brec->exons[i-1].end+1 << endl;
	// 				// cout << "brec->exons[i].start-1 : " << brec->exons[i].start-1 << endl;

	// 			} else {
	// 				// BAM_unmapped_set.insert(*bamkey);
  //         outfile_cleaned->write(brec);
	// 			}
	// 			// cout << "\tIntron: " << brec->refName() << "; " << brec->spliceStrand() << "; " << brec->exons[i-1].end+1 << " - " << brec->exons[i].start-1 << endl;	
	// 			// CJunc j(brec->exons[i-1].end+1, brec->exons[i].start-1, strand,
	// 			// 		dupcount);
	// 		}
	// 	} else {
  //     // outfile_cleaned->write(brec);
  //   }
	// }

    
    delete outfile_discard;
    cout << "Done delete outfile_discard!" << endl;
    delete outfile_spliced;
    cout << "Done delete outfile_spliced!" << endl;
    delete outfile_cleaned;
    cout << "Done delete outfile_cleaned!" << endl;

    return 0;
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;SLPEDVho:N:Q:F:M:J:");
    args.printError(USAGE, true);

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s",USAGE);
        exit(0);
    }

    if (args.getOpt('v') || args.getOpt("version")) {
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

    junctionfname = outfdir + "/bed/junction.bed";
    // if (args.getOpt('J')) {
    //   junctionfname=args.getOpt('J');
    //   if (fileExists(junctionfname.chars())>1) {
    //       GMessage("\nJunction bed file: ", junctionfname.chars());
    //   } else {
    //     GMessage("Error: bed file (%s) not found.\n",
    //         junctionfname.chars());
    //   }
    // } else {
    //       GMessage(USAGE);
    //       GMessage("\nError: junction fa file must be provided (-J)!\n");
    //       exit(1);
    // }








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
        fprintf(stderr, "Running SPLAM " VERSION ". Command line:\n");
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

void flushBrec(GVec<GSamRecord*> &pbrecs, unordered_map<string, int> &hits, GSamWriter* outfile_cleaned) {
    if (pbrecs.Count()==0) return;
    for (int i=0; i < pbrecs.Count(); i++) {
        std::string kv = pbrecs[i]->name();

        cout << "kv: " << kv << endl;
        std::string tmp = std::to_string(pbrecs[i]->pairOrder());
        kv += ";";
        kv += tmp;
        if (hits.find(kv) != hits.end()) {
            cout << "Update NH tage!!" << endl;
            int new_nh = pbrecs[i]->tag_int("NH", 0) - hits[kv];
            pbrecs[i]->add_int_tag("NH", new_nh);
        }
// using GArray
//        CRead tmp(pbrecs[i]->pairOrder(), pbrecs[i]->name());
//        int idx;
//        if (spur_reads.Found(tmp, idx)) {
//            int new_nh = pbrecs[i]->tag_int("NH", 0) - spur_reads[idx].spurCnt;
//            pbrecs[i]->add_int_tag("NH", new_nh);
//        }
// for testing
       if (!strcmp(pbrecs[i]->name(), "ERR188044.24337229")) {
           std::cout << pbrecs[i]->tag_int("NH", 0) << std::endl;
       }
        outfile_cleaned->write(pbrecs[i]);
    }
}

void loadBed(GStr inbedname, GArray<CJunc> &spur_juncs, unordered_map<string, string> &chrs_refseq_2_ucsc) {

    // Bed score is in ncbi chr name.
    std::ifstream bed_f(inbedname);
    std::string line;
    int bed_counter = 0;
    while (getline(bed_f, line)) {
        bed_counter ++;
        // cout << "line: " << line << endl;
        GStr gline = line.c_str();
        GVec<GStr> junc;
        int cnt = 0;
        while (cnt < 7) {
            GStr tmp = gline.split("\t");
            // cout << "tmp: " << tmp << endl;
            junc.Add(gline);
            gline=tmp;
            cnt++;
        }
        char* chrname =junc[0].detach();
        // char chr = chrname[strlen(chrname) - 1];
        string chr_str(chrname);
        cout << "1 chr_str: " << chr_str << endl;
        chr_str = chrs_refseq_2_ucsc[chr_str];
        cout << "2 chr_str: " << chr_str << endl;


        // cout << ">> chrs_convert: " << endl;
        // for (auto i : chrs_convert) {
        //     cout << i.first << " ---- " << i.second << endl;
        // }

        // cout << "1 chr_str: " << chr_str << endl;
        // chr_str = chrs_convert[chr_str];
        // cout << "2 chr_str: " << chr_str << endl;
        
        // if (true) {
        //     cout << "junc[1].asInt(): " << junc[1].asInt() << endl;
        //     cout << "junc[2].asInt(): " << junc[2].asInt() << endl;
        //     cout << "junc[5].detach(): " << *junc[5].detach() << endl;
        //     cout << "junc[6].: " << junc[6].asDouble() << endl;
        //     cout << "chr: " << chr << endl;
        // }
        if (junc[6].asDouble() <= 0.2) {
            cout << "junc[6].asDouble(): " << junc[6].asDouble() << endl;

        	// CJunc(int vs=0, int ve=0, char vstrand='+', string vref=".", uint64_t dcount=1):
            CJunc j(junc[1].asInt()+1, junc[2].asInt(), *junc[5].detach(), chr_str);

            spur_juncs.Add(j);
            cout << "spur_juncs.size: " << spur_juncs.Count() << endl;
        }
    }

    cout << "bed_counter: " << bed_counter << endl;
}
