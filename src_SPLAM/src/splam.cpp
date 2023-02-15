#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <unordered_map>
#include <filesystem>

#include "tmerge.h"
#include "util.h"
#include "junc.h"
#include "filter.h"

#include <gclib/GArgs.h>
#include <gclib/GBase.h>
#include <gclib/GStr.h>

#include <htslib/htslib/faidx.h>
#include <Python.h>

#define VERSION "0.0.1"

void processOptions(int argc, char* argv[]);
void loadBed(GStr inbedname, GArray<CJunc> &spur_juncs);
void flushBrec(GVec<GSamRecord*> &pbrecs, std::unordered_map<std::string, int> &hits, GSamWriter* outfile_cleaned);


wchar_t *GetWC(const char *c)
{
    const size_t cSize = strlen(c)+1;
    wchar_t* wc = new wchar_t[cSize];
    mbstowcs (wc, c, cSize);

    return wc;
}


const char* USAGE = "SPLAM v" VERSION "\n"
                              "==========================================================================================\n"
                              "An accurate spliced alignment pruner and spliced junction predictor.\n"
                              "==========================================================================================\n";


GStr out_dir;
GSamWriter* outfile_discard = NULL;
GSamWriter* outfile_spliced = NULL;
GSamWriter* outfile_cleaned = NULL;

GStr model_name;
GStr infname_reffa;
GStr outfname_junction;
FILE* joutf=NULL;

bool verbose = false;
TInputFiles in_records;
TInputRecord* irec=NULL;
GSamRecord* brec=NULL;
float threshold = 0.2;

std::unordered_map<std::string, int> get_hg38_chrom_size(std::string target) {
    std::unordered_map<std::string, int> chrs;
    std::string ref_file;
    if (target == "STAR") {
        ref_file = "./hg38_chrom_size_refseq.tsv";
    } else {
        ref_file = "./hg38_chrom_size.tsv";
    }
    std::ifstream ref_f(ref_file);
    std::string line;
    while(getline(ref_f, line)){
        std::string chromosome;
        int len;
        std::replace(line.begin(), line.end(), '\t', ' ');
        std::stringstream ss(line);
        ss >> chromosome;
        ss >> len;
        chrs[chromosome] = len;
    }   
    // for (auto i : chrs) {
    //   std::cout << i.first << " ---- " << i.second << std::endl;
    // }
    return chrs;
}

int main(int argc, char* argv[]) {

    for(int i=0;i<argc-1;i++)
      std::cout << argv[i] << " "; 
    std::cout << std::endl;


    std::cout << "************************" << std::endl;
    std::cout << "* Start of the program: " << std::endl;
    std::cout << "************************" << std::endl;


    const char *banner = R"""(
  ███████╗██████╗ ██╗      █████╗ ███╗   ███╗    ██╗
  ██╔════╝██╔══██╗██║     ██╔══██╗████╗ ████║    ██║
  ███████╗██████╔╝██║     ███████║██╔████╔██║    ██║
  ╚════██║██╔═══╝ ██║     ██╔══██║██║╚██╔╝██║    ╚═╝
  ███████║██║     ███████╗██║  ██║██║ ╚═╝ ██║    ██╗
  ╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝    ╚═╝
    )""";
    std::cout << banner << std::endl;
    
    in_records.setup(VERSION, argc, argv);
    processOptions(argc, argv);
    
    int num_samples=in_records.start();
    std::cout << "* Number of samples  : " << num_samples << std::endl;
    std::cout << "* Output directory   : " << out_dir << std::endl;
    std::cout << "* Junction file name : " << outfname_junction << std::endl;
    
    //  Creating directories for bed & fasta files.
    GStr bed_dir(out_dir+"/bed");
    GStr fasta_dir(out_dir+"/fasta");
    GStr score_dir(out_dir+"/scores");
    std::filesystem::create_directories(bed_dir.chars());
    std::filesystem::create_directories(fasta_dir.chars());
    std::filesystem::create_directories(score_dir.chars());

    // Initializing the BAM files.
    GStr outfname_discard(out_dir + "/discard.bam");
    GStr outfname_spliced(out_dir + "/spliced.bam");
    GStr outfname_cleaned(out_dir + "/cleaned.bam");
    outfile_discard = new GSamWriter(outfname_discard, in_records.header(), GSamFile_BAM);
    outfile_spliced = new GSamWriter(outfname_spliced, in_records.header(), GSamFile_BAM);
    outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);

    std::unordered_map<std::string, int> chrs = get_hg38_chrom_size("HISAT2");

    FILE *ref_fa_f;
    if ((ref_fa_f = fopen(infname_reffa.chars(), "r"))) {
        fclose(ref_fa_f);
        std::cout << "> FASTA index \"" << infname_reffa << "i\" has been created!" << std::endl;
    } else {
        int res = fai_build(infname_reffa.chars());
        std::cout << "> Creating FASTA index \"" << infname_reffa << "i\"" << std::endl;
    }

    faidx_t * ref_faidx = fai_load(infname_reffa.chars());
    std::cout << "> Loading FASTA file" << std::endl;

    /*********************************************
     * Step 1: generating spliced junctions in BED
    *********************************************/
    std::cout << "********************************************" << std::endl;
    std::cout << "** Step 1: generating spliced junctions in BED" << std::endl;
    std::cout << "********************************************" << std::endl;

    // Creating the output junction bed file
    if (!outfname_junction.is_empty()) {
        if (strcmp(outfname_junction.substr(outfname_junction.length()-4, 4).chars(), ".bed")!=0) {
            outfname_junction.append(".bed");
        }
        joutf = fopen(outfname_junction.chars(), "w");
        if (joutf==NULL) GError("Error creating file %s\n", outfname_junction.chars());
        // fprintf(joutf, "track name=junctions\n");
    }

    // Reading BAM file.
    int counter = 0;
    int prev_tid=-1;
    GStr prev_refname;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    int b_end=0;
    int b_start=0;

    while ((irec=in_records.next())!=NULL) {
        brec=irec->brec;
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
        if (joutf && brec->exons.Count()>1) {
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
    std::cout << "********************************************" << std::endl;
    std::cout << "** Step 2: getting coordinates of donors and acceptors" << std::endl;
    std::cout << "********************************************" << std::endl;

    int SEQ_LEN = 800;
    int QUOTER_SEQ_LEN = SEQ_LEN/4;

    /*********************************************
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    *******************************************/
    GStr donor_bed(out_dir + "/bed/donor.bed");
    GStr acceptor_bed(out_dir + "/bed/acceptor.bed");
    GStr da_bed(out_dir + "/bed/d_a.bed");

    std::ofstream outfile_bed_donor(donor_bed);
    std::ofstream outfile_bed_acceptor(acceptor_bed);
    std::ofstream outfile_bed_da(da_bed);

    GStr junc_fasta(out_dir + "/fasta/junction.fa");
    GStr donor_fasta(out_dir + "/fasta/donor.fa");
    GStr accceptor_fasta(out_dir + "/fasta/acceptor.fa");

    std::ofstream outfile_fa_junc(junc_fasta);
    std::ofstream outfile_fa_donor(donor_fasta);
    std::ofstream outfile_fa_acceptor(accceptor_fasta);

    std::ifstream fr_junc(outfname_junction);
    std::string line;
    while(getline(fr_junc, line)){
        std::cout << "fr_junc: " << line << std::endl;

        std::string chromosome;
        int start = 0, end = 0;
        std::string junc_name;
        int num_alignment;
        std::string strand;
        std::replace(line.begin(), line.end(), '\t', ' ');
        std::stringstream ss(line);

        ss >> chromosome >> start >> end >> junc_name >> num_alignment >> strand;


        std::cout << "** chromosome: " << chromosome << " ";
        std::cout << "start: " << start << " ";
        std::cout << "end: " << end << " ";
        std::cout << "junc_name: " << junc_name << " ";
        std::cout << "num_alignment: " << num_alignment << " ";
        std::cout << "strand: " << strand << " " << std::endl;

        int splice_junc_len = 0;
        int flanking_size = QUOTER_SEQ_LEN;

        int donor = 0, acceptor = 0;
        int donor_s = 0, donor_e = 0;
        int acceptor_s = 0, acceptor_e = 0;
        splice_junc_len = end - start;
        if (splice_junc_len < QUOTER_SEQ_LEN) {
            flanking_size = splice_junc_len;
        }

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

        // std::cout << "chrs[chromosome]: " << chrs[chromosome] << std::endl;
        // std::cout << "donor_e         : " << donor_e << std::endl;

        // std::cout << "chrs[chromosome]: " << chrs[chromosome] << std::endl;
        // std::cout << "acceptor_e      : " << acceptor_e << std::endl;

        if (donor_e >= chrs[chromosome] or acceptor_e >= chrs[chromosome]) {
            std::cout << "Skip!!" << std::endl;
            continue;
        }
        if (donor_s < 0 or acceptor_s < 0) {
            std::cout << "Skip!!" << std::endl;
            continue;
        }


        outfile_bed_donor << chromosome << "\t" + std::to_string(donor_s) + "\t" + std::to_string(donor_e) + "\t" + junc_name+"_donor" + "\t" + std::to_string(num_alignment) + "\t" + strand + "\n";
        outfile_bed_acceptor << chromosome << "\t" + std::to_string(acceptor_s) + "\t" + std::to_string(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + std::to_string(num_alignment) + "\t" + strand + "\n";

        int donor_len = donor_e - donor_s;
        char* donor_seq = faidx_fetch_seq(ref_faidx, chromosome.c_str(), donor_s, donor_e-1, &donor_len);

        int acceptor_len = acceptor_e - acceptor_s;
        char* acceptor_seq = faidx_fetch_seq(ref_faidx, chromosome.c_str(), acceptor_s, acceptor_e-1, &acceptor_len);

        if (strand == "-") {
            hts_pos_t donor_len_hts = donor_e - donor_s;
            hts_pos_t acceptor_len_hts = acceptor_e - acceptor_s;
            reverse_complement(donor_seq, donor_len_hts);
            reverse_complement(acceptor_seq, acceptor_len_hts);
        }

        if (donor_seq == NULL){
            printf("c is empty\n");
            continue;
        }

        if (acceptor_seq == NULL) {
            printf("c is empty\n");
            continue;
        }

        outfile_fa_donor << ">" << chromosome << std::endl;
        outfile_fa_donor << donor_seq << std::endl;
        outfile_fa_acceptor << ">" << chromosome << std::endl;
        outfile_fa_acceptor << acceptor_seq << std::endl;

        outfile_fa_junc << ">" << chromosome << ";" << std::to_string(donor) << ";" << std::to_string(acceptor) << ";" << strand << std::endl;
        if (strlen(donor_seq) >= 400) {
            outfile_fa_junc << donor_seq << acceptor_seq << std::endl;
        } else {
            outfile_fa_junc << donor_seq << std::string(2*(400-(int)strlen(donor_seq)), 'N') << acceptor_seq << std::endl;
        }

        std::cout << "donor   : " << donor_seq[200] << donor_seq[201] << std::endl;
        std::cout << "acceptor: " << acceptor_seq[198] << acceptor_seq[199] << std::endl;

        if (strand == "+") {
            outfile_bed_da << chromosome + "\t" + std::to_string(donor) + "\t" + std::to_string(acceptor+1) + "\tJUNC\t" + std::to_string(num_alignment) + "\t" + strand + "\n";
        } else if (strand == "-") {
            outfile_bed_da << chromosome + "\t" + std::to_string(acceptor) + "\t" + std::to_string(donor+1) + "\tJUNC\t" + std::to_string(num_alignment) + "\t" + strand + "\n";
        }
        std::cout << std::endl;
    }  










    /*********************************************
     * Step 3: SPLAM model prediction
    *********************************************/

    Py_Initialize();
    GStr python_f = "./script/prediction.py";
    GStr outfile_junction_score(score_dir + "/junction_score.bed");

    std::cout << "outfile_junction_score: " << outfile_junction_score.chars() << std::endl;
    wchar_t *argvv[] = {GetWC(python_f.chars()), L"-f", GetWC(junc_fasta.chars()), L"-o", GetWC(outfile_junction_score.chars()), L"-m", GetWC(model_name.chars())};
    PySys_SetArgv(7, argvv);
    
    FILE *file = _Py_fopen(python_f.chars(), "r");

    if(file != NULL) {

        int re = PyRun_SimpleFileExFlags(file, python_f.chars(), false, {0});

        std::cout << "************************" << std::endl;
        std::cout << "Results: " << re << std::endl;
        std::cout << "************************" << std::endl;
    }


    // Py_Initialize();
    // PyRun_SimpleString("from time import time,ctime\n"
    //                     "print('Today is',ctime(time()))\n");
    Py_Finalize();












    // /*********************************************
    //  * Step 4: SPLAM filtering out reads.
    // *********************************************/
    // // GSamReader bamreader(inbamname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    // // outfile=new GSamWriter(outfname, bamreader.header(), GSamFile_BAM);
    // // GSamRecord brec;

    // std::cout << "********************************************" << std::endl;
    // std::cout << "** Step 4: SPLAM filtering out reads." << std::endl;
    // std::cout << "********************************************" << std::endl;

    // auto start=std::chrono::high_resolution_clock::now();
    // int spur_cnt = 0;
    // std::cout << "brrrm! identifying alignment records with spurious splice junctions" << std::endl;

    // GArray<CJunc> spur_juncs;
    // GVec<GSamRecord*> kept_brecs;
    // std::unordered_map<std::string, int> hits;
    // GStr inbedname(out_dir + "/output/junc_scores.bed");
    // loadBed(inbedname, spur_juncs);

    // std::cout << ">> (4) Junction count: " << junctions.Count() << std::endl;
	// for (int i = 0; i < junctions.Count(); i++) {
	// 	std::cout << i <<  " (4) Junction name: " << junctions[i].start << " - " << junctions[i].end << std::endl;
	// 	std::cout << ">> (4) Read count: " << junctions[i].read_ls.size() << std::endl;

    //     std::cout << "junctions[i].ref: " << junctions[i].ref << std::endl;
    //     std::cout << "junctions[i].start: " << junctions[i].start << std::endl;
    //     std::cout << "junctions[i].end: " << junctions[i].end << std::endl;
    //     std::cout << "junctions[i].strand: " << junctions[i].strand << std::endl;
    //     std::cout << std::endl;

    //     CJunc jnew(junctions[i].start, junctions[i].end, junctions[i].strand, junctions[i].ref);
        
    //     std::cout << "spur_juncs.Exists(jnew):  " << spur_juncs.Exists(jnew) << std::endl;
    //     if (spur_juncs.Exists(jnew)) {
    //         // spur = true;
    //         std::cout << "spur_juncs.Exists! " << std::endl;
    //         for (int j=0; j<junctions[i].read_ls.size(); j++) {
    //             std::cout << "~~ SPLAM!" << std::endl;
    //             outfile_discard->write(junctions[i].read_ls.at(j));
    //             // delete junctions[i].read_ls.at(j);
    //         }
    //         std::cout << "spur_juncs.Exists Done! " << std::endl;
    //     } else {
    //         for (int j=0; j<junctions[i].read_ls.size(); j++) {
    //             bool spur = false;
    //             int r_exon_count = junctions[i].read_ls.at(j)->exons.Count();
    //             if (r_exon_count > 1) {
    //                 for (int e=1; e<r_exon_count; e++) {
    //                     CJunc jnew_sub(junctions[i].read_ls.at(j)->exons[e-1].end, junctions[i].read_ls.at(j)->exons[e-1].start-1, junctions[i].strand, junctions[i].ref);
    //                     if (spur_juncs.Exists(jnew_sub)) {
    //                         spur = true;
    //                         break;
    //                     }
    //                 }
    //             }



    // //         if (!spur) {
    // //             // std::cout << "Not spurious!" << std::endl;
    // //             GSamRecord *rec = new GSamRecord(*brec);
    // //             PBRec *newpbr = new PBRec(rec);
    // //             kept_brecs.Add(newpbr);
    // //         } else {
    // //             std::cout << "Spurious!" << std::endl;
    // //             spur_cnt++;
    // //             std::string kv = brec->name();
    // //             std::string tmp = std::to_string(brec->pairOrder());
    // //             kv += ";";
    // //             kv += tmp;
    // //             // key not present
    // //             if (hits.find(kv) == hits.end()) {
    // //                 hits[kv] = 1;
    // //             } else {
    // //                 int val = hits[kv];
    // //                 val++;
    // //                 hits[kv] = val;
    // //             }
    // //         }

    //             if (spur) {
    //                 std::cout << "spur_juncs.Exists! " << std::endl;
    //                 std::cout << "~~ SPLAM!" << std::endl;
    //                 std::string kv = brec->name();
    //                 std::string tmp = std::to_string(brec->pairOrder());
    //                 kv += ";";
    //                 kv += tmp;

    //                 if (hits.find(kv) == hits.end()) {
    //                     hits[kv] = 1;
    //                 } else {
    //                     int val = hits[kv];
    //                     val++;
    //                     hits[kv] = val;
    //                 }

    //                 outfile_discard->write(junctions[i].read_ls.at(j));
    //                 // delete junctions[i].read_ls.at(j);
    //             } else {
    //                 std::cout << "spur_juncs not Exists! " << std::endl;
    //                 kept_brecs.Add(junctions[i].read_ls.at(j));

    //                 // std::cout << "~~ Clean!" << std::endl;
    //                 // outfile_cleaned->write(junctions[i].read_ls.at(j));
    //                 // delete junctions[i].read_ls.at(j);
    //             }
    //         }
    //         std::cout << "~~~ Done! " << std::endl;
    //     }

    //     std::cout << "Done!" << std::endl;

		
    //     // for (int r = 0; r < junctions[i].read_ls.size(); r++) {
	// 	// 	std::cout << "\tRead " <<  r << " : " << junctions[i].read_ls[r]->cigar() << std::endl;
	// 	// }
	// 	// std::cout << std::endl;
	// }

    // flushBrec(kept_brecs, hits, outfile_cleaned);
    // auto end =std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    // std::cout << spur_cnt << " spurious alignments were removed." << std::endl;
    // std::cout << "Completed in " << duration.count() << " seconds" << std::endl;
    


	// while ((irec=in_records.next())!=NULL) {
	// 	brec=irec->brec;
	// 	// std::cout << irec->fidx << std::endl;
	// 	if (brec->hasIntrons()) {
    //         const char* chrname=brec->refName();
    //         char chr = chrname[strlen(chrname) - 1];
    //         char strand = brec->spliceStrand();
    //         bool spur = false;

    //         // std::cout << "chrname: " << chrname << std::endl;
    //         for (int i = 1; i < brec->exons.Count(); i++) {
    //             CJunc j(brec->exons[i-1].end, brec->exons[i].start-1, strand, chr);
    //             if (spur_juncs.Exists(j)) {
    //                 spur = true;
    //                 break;
    //             }
    //         }
    //         if (!spur) {
    //             // std::cout << "Not spurious!" << std::endl;
    //             GSamRecord *rec = new GSamRecord(*brec);
    //             PBRec *newpbr = new PBRec(rec);
    //             kept_brecs.Add(newpbr);
    //         } else {
    //             std::cout << "Spurious!" << std::endl;
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
    // std::cout << ">> Loading the model: " << std::endl;
    // torch::jit::script::Module module;
    // // module = torch::jit::load(model_name.chars());

    // try {
    //     // Deserialize the ScriptModule from a file using torch::jit::load().
    //     std::cout << "Loading "<< model_name.chars() << std::endl;
    //     module = torch::jit::load(model_name.chars());
    // }
    // catch (const c10::Error& e) {
    //     std::cerr << "error loading the model\n";
    //     return -1;
    // }
    // std::cout << "Model "<< model_name.chars() <<" loaded fine\n";

  // /************************
  //  * Processing JUNCTION file.
  //  ************************/
  // set<intron_key> bad_intron_set;
  // set<intron_key> good_intron_set;
  // std::string line;
  // std::ifstream myfile (outfname_junction); // this is equivalent to the above method
  // if ( myfile.is_open() ) { // always check whether the file is open

  //     // chr18	21682058	21596619	JUNC_144610	0	+	8.2311524e-10	4.0248174e-12
  //     while ( std::getline(myfile, line) ) {

  //       int parser_counter = 0;

  //       std::string chr = "";
  //       int start = 0;
  //       int end = 0;
  //       std::string junc_name = "";
  //       int tmp = 0;
  //       char strand = ' ';
  //       float d_score = .0;
  //       float a_score = .0;

  //       std::string token;
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

  //       // std::cout << "chr: " << chr << std::endl;
  //       // std::cout << "start: " << start << std::endl;
  //       // std::cout << "end: " << end << std::endl;
  //       // std::cout << "junc_name: " << junc_name << std::endl;
  //       // std::cout << "tmp: " << tmp << std::endl;
  //       // std::cout << "strand: " << strand << std::endl;
  //       // std::cout << "d_score: " << d_score << std::endl;
  //       // std::cout << "a_score: " << a_score << std::endl;

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
  // std::cout << "bad_intron_set.size(): " << bad_intron_set.size() << std::endl;
  // std::cout << "good_intron_set.size(): " << good_intron_set.size() << std::endl;




  // /************************
  //  * Processing BAM file.
  //  ************************/
	// int counter = 0;
	// while ((irec=in_records.next())!=NULL) {
	// 	brec=irec->brec;
	// 	// std::cout << irec->fidx << std::endl;
	// 	if (brec->hasIntrons()) {
	// 		// This is a spliced read => start processing it!
	// 		// char strand = brec->spliceStrand();
	// 		// std::cout << "strand       : " << strand << std::endl;
	// 		// std::cout << "brec->cigar(): " << brec->cigar() << std::endl;
	// 		for (int i=1;i<brec->exons.Count();i++) {
	// 			// int strand = 0;
	// 			// if (brec->spliceStrand() == '+') strand = 1;
	// 			// if (brec->spliceStrand() == '-') strand = -1;
	// 			// std::cout << "brec->refName(): " << brec->refName()<< std::endl;
	// 			std::string bamseq_name(brec->refName());
	// 			// std::cout << "bamseq_name: " << bamseq_name << std::endl;

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
	// 				std::cout << "Splam!  " << counter << std::endl;
	//         outfile_discard->write(brec);
	// 				// std::cout << "brec->spliceStrand()   : " << brec->spliceStrand() << std::endl;
	// 				// std::cout << "brec->refName(),       : " << brec->refName() << std::endl;
	// 				// std::cout << "brec->exons[i-1].end+1 : " << brec->exons[i-1].end+1 << std::endl;
	// 				// std::cout << "brec->exons[i].start-1 : " << brec->exons[i].start-1 << std::endl;

	// 			} else {
	// 				// BAM_unmapped_set.insert(*bamkey);
  //         outfile_cleaned->write(brec);
	// 			}
	// 			// std::cout << "\tIntron: " << brec->refName() << "; " << brec->spliceStrand() << "; " << brec->exons[i-1].end+1 << " - " << brec->exons[i].start-1 << std::endl;	
	// 			// CJunc j(brec->exons[i-1].end+1, brec->exons[i].start-1, strand,
	// 			// 		dupcount);
	// 		}
	// 	} else {
  //     // outfile_cleaned->write(brec);
  //   }
	// }

    
    delete outfile_discard;
    std::cout << "Done delete outfile_discard!" << std::endl;
    delete outfile_spliced;
    std::cout << "Done delete outfile_spliced!" << std::endl;
    delete outfile_cleaned;
    std::cout << "Done delete outfile_cleaned!" << std::endl;

    return 0;
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;SLPEDVho:N:Q:F:M:J:R:");
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
    out_dir=args.getOpt('o');
    if (out_dir.is_empty()) {
        GMessage(USAGE);
        GMessage("\nError: output filename must be provided (-o)!\n");
        exit(1);
    }
    outfname_junction = out_dir + "/bed/junction.bed";
    // if (args.getOpt('J')) {
    //   outfname_junction=args.getOpt('J');
    //   if (fileExists(outfname_junction.chars())>1) {
    //       GMessage("\nJunction bed file: ", outfname_junction.chars());
    //   } else {
    //     GMessage("Error: bed file (%s) not found.\n",
    //         outfname_junction.chars());
    //   }
    // } else {
    //       GMessage(USAGE);
    //       GMessage("\nError: junction fa file must be provided (-J)!\n");
    //       exit(1);
    // }

    if (args.getOpt('R')) {
        infname_reffa=args.getOpt('R');
        if (fileExists(infname_reffa.chars())>1) {
            // guided=true;
        } else {
            GError("Error: reference fasta file (%s) not found.\n",
                infname_reffa.chars());
        }
    } else {
          GMessage(USAGE);
          GMessage("\nError: reference fasta file must be provided (-R)!\n");
          exit(1);
    }

    if (args.getOpt('M')) {
        model_name=args.getOpt('M');
        if (fileExists(model_name.chars())>1) {
            // guided=true;
        } else {
            GError("Error: model file (%s) not found.\n",
                model_name.chars());
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
        std::cout << "absolute_ifn: " << absolute_ifn << std::endl;
        in_records.addFile(absolute_ifn.c_str());
    }
}

void flushBrec(GVec<GSamRecord*> &pbrecs, std::unordered_map<std::string, int> &hits, GSamWriter* outfile_cleaned) {
    if (pbrecs.Count()==0) return;
    for (int i=0; i < pbrecs.Count(); i++) {
        std::string kv = pbrecs[i]->name();

        std::cout << "kv: " << kv << std::endl;
        std::string tmp = std::to_string(pbrecs[i]->pairOrder());
        kv += ";";
        kv += tmp;
        if (hits.find(kv) != hits.end()) {
            std::cout << "Update NH tage!!" << std::endl;
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

void loadBed(GStr inbedname, GArray<CJunc> &spur_juncs) {

    // Bed score is in ncbi chr name.
    std::ifstream bed_f(inbedname);
    std::string line;
    int bed_counter = 0;
    while (getline(bed_f, line)) {
        bed_counter ++;
        // std::cout << "line: " << line << std::endl;
        GStr gline = line.c_str();
        GVec<GStr> junc;
        int cnt = 0;
        while (cnt < 7) {
            GStr tmp = gline.split("\t");
            // std::cout << "tmp: " << tmp << std::endl;
            junc.Add(gline);
            gline=tmp;
            cnt++;
        }
        char* chrname =junc[0].detach();
        // char chr = chrname[strlen(chrname) - 1];
        GStr chr_str(chrname);
        std::cout << "1 chr_str: " << chr_str.chars() << std::endl;


        // std::cout << ">> chrs_convert: " << std::endl;
        // for (auto i : chrs_convert) {
        //     std::cout << i.first << " ---- " << i.second << std::endl;
        // }

        // std::cout << "1 chr_str: " << chr_str << std::endl;
        // chr_str = chrs_convert[chr_str];
        // std::cout << "2 chr_str: " << chr_str << std::endl;
        
        // if (true) {
        //     std::cout << "junc[1].asInt(): " << junc[1].asInt() << std::endl;
        //     std::cout << "junc[2].asInt(): " << junc[2].asInt() << std::endl;
        //     std::cout << "junc[5].detach(): " << *junc[5].detach() << std::endl;
        //     std::cout << "junc[6].: " << junc[6].asDouble() << std::endl;
        //     std::cout << "chr: " << chr << std::endl;
        // }
        if (junc[6].asDouble() <= threshold) {
            std::cout << "junc[6].asDouble(): " << junc[6].asDouble() << std::endl;

        	// CJunc(int vs=0, int ve=0, char vstrand='+', std::string vref=".", uint64_t dcount=1):
            CJunc j(junc[1].asInt()+1, junc[2].asInt(), *junc[5].detach(), chr_str);

            spur_juncs.Add(j);
            std::cout << "spur_juncs.size: " << spur_juncs.Count() << std::endl;
        }
    }

    std::cout << "bed_counter: " << bed_counter << std::endl;
}


// void get_Refseq_2_UCSC_chr_names(std::unordered_map<std::string, std::string> &chrs_refseq_2_ucsc, std::unordered_map<std::string, std::string> &chrs_ucsc_2_refseq) {

//     std::string ref_file;
//     ref_file = "../../../src/Refseq_2_UCSU_chromosome_names.tsv";
//     std::cout << "ref_file: " << ref_file << std::endl;
//     std::ifstream ref_f(ref_file);
//     std::string line;
//     while(getline(ref_f, line)){
//         std::string refseq;
//         std::string ucsc;
//         std::replace(line.begin(), line.end(), '\t', ' ');
//         std::stringstream ss(line);
//         ss >> refseq;
//         ss >> ucsc;
//         chrs_ucsc_2_refseq[ucsc] = refseq;
//         chrs_refseq_2_ucsc[refseq] = ucsc;
//     }   

//     for (auto i : chrs_refseq_2_ucsc) {
//       std::cout << i.first << " ---- " << i.second << std::endl;
//     }

//     for (auto i : chrs_ucsc_2_refseq) {
//       std::cout << i.first << " ---- " << i.second << std::endl;
//     }
// }