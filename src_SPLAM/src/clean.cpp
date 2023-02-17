#include "clean.h"
#include "common.h"
#include "junc.h"
#include "junc_func.h"
#include "util.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <htslib/htslib/faidx.h>
#include <Python.h>

GSamRecord* brec=NULL;
GSamWriter* outfile_discard = NULL;
GSamWriter* outfile_spliced = NULL;
GSamWriter* outfile_cleaned = NULL;
FILE* joutf=NULL;
// GArray<CJunc> junctions;

void splamClean() {
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
    std::cout << "> Loading FASTA file\n" << std::endl;


    /*********************************************
     * Step 1: generating spliced junctions in BED
    *********************************************/
    std::cout << "********************************************" << std::endl;
    std::cout << "** Step 1: generating spliced junctions in BED" << std::endl;
    std::cout << "********************************************\n" << std::endl;

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
    std::unordered_map<std::string, int> doner_dimers;
    std::unordered_map<std::string, int> acceptor_dimers;


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
        std::string chromosome;
        int start = 0, end = 0;
        std::string junc_name;
        int num_alignment;
        std::string strand;
        std::replace(line.begin(), line.end(), '\t', ' ');
        std::stringstream ss(line);

        ss >> chromosome >> start >> end >> junc_name >> num_alignment >> strand;
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

        // std::cout << "donor   : " << donor_seq[200] << donor_seq[201] << std::endl;
        // std::cout << "acceptor: " << acceptor_seq[198] << acceptor_seq[199] << std::endl;

        
        // char *donor_dim = (char*)malloc(2);
        // std::memcpy(donor_dim, donor_seq[200], 2);

        char donor_dim[3];
        memcpy(donor_dim, &donor_seq[200], 2);
        donor_dim[2] = '\0';
        std::string donor_str(donor_dim);

        char acceptor_dim[3];
        memcpy(acceptor_dim, &acceptor_seq[198], 2);
        acceptor_dim[2] = '\0';
        std::string acceptor_str(acceptor_dim);

        if (doner_dimers.find(donor_str) == doner_dimers.end()) {
            doner_dimers[donor_str] = 1;
        } else {
            doner_dimers[donor_str] += 1;
        }

        if (acceptor_dimers.find(acceptor_str) == acceptor_dimers.end()) {
            acceptor_dimers[acceptor_str] = 1;
        } else {
            acceptor_dimers[acceptor_str] += 1;
        }
        if (strand == "+") {
            outfile_bed_da << chromosome + "\t" + std::to_string(donor) + "\t" + std::to_string(acceptor+1) + "\tJUNC\t" + std::to_string(num_alignment) + "\t" + strand + "\n";
        } else if (strand == "-") {
            outfile_bed_da << chromosome + "\t" + std::to_string(acceptor) + "\t" + std::to_string(donor+1) + "\tJUNC\t" + std::to_string(num_alignment) + "\t" + strand + "\n";
        }
    }  

    std::cout << ">> Donor dimers: " << std::endl;
    for (auto i : doner_dimers) {
        std::cout << "\t" << i.first << ": " << i.second << std::endl;
    }

    std::cout << ">> Acceptor dimers: " << std::endl;
    for (auto i : acceptor_dimers) {
        std::cout << "\t" << i.first << ": " << i.second << std::endl;
    }
        

    /*********************************************
     * Step 3: SPLAM model prediction
    *********************************************/
    std::cout << "\n********************************************" << std::endl;
    std::cout << "** Step 3: SPLAM! model prediction" << std::endl;
    std::cout << "********************************************\n" << std::endl;
    Py_Initialize();
    GStr python_f = "./script/splam.py";
    GStr outfile_junction_score(score_dir + "/junction_score.bed");
    wchar_t *argvv[] = {GetWC(python_f.chars()), L"-f", GetWC(junc_fasta.chars()), L"-o", GetWC(outfile_junction_score.chars()), L"-m", GetWC(infname_model_name.chars())};
    PySys_SetArgv(7, argvv);
    FILE *file = _Py_fopen(python_f.chars(), "r");
    if(file != NULL) {
        int re = PyRun_SimpleFileExFlags(file, python_f.chars(), false, {0});
    }
    Py_Finalize();


    /*********************************************
     * Step 4: SPLAM filtering out reads.
    *********************************************/
    std::cout << "********************************************" << std::endl;
    std::cout << "** Step 4: SPLAM filtering out reads." << std::endl;
    std::cout << "********************************************\n" << std::endl;
    // GSamReader bamreader(inbamname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    // outfile=new GSamWriter(outfname, bamreader.header(), GSamFile_BAM);
    // GSamRecord brec;

    std::cout << "********************************************" << std::endl;
    std::cout << "** Step 4: SPLAM filtering out reads." << std::endl;
    std::cout << "********************************************" << std::endl;

    auto start=std::chrono::high_resolution_clock::now();
    int spur_cnt = 0;
    std::cout << "brrrm! identifying alignment records with spurious splice junctions" << std::endl;

    GArray<CJunc> spur_juncs;
    GVec<GSamRecord*> kept_brecs;
    std::unordered_map<std::string, int> hits;
    GStr inbedname(out_dir + "/output/junc_scores.bed");
    loadBed(inbedname, spur_juncs);

    std::cout << ">> (4) Junction count: " << junctions.Count() << std::endl;
	for (int i = 0; i < junctions.Count(); i++) {
		std::cout << i <<  " (4) Junction name: " << junctions[i].start << " - " << junctions[i].end << std::endl;
		std::cout << ">> (4) Read count: " << junctions[i].read_ls.size() << std::endl;

        std::cout << "junctions[i].ref: " << junctions[i].ref << std::endl;
        std::cout << "junctions[i].start: " << junctions[i].start << std::endl;
        std::cout << "junctions[i].end: " << junctions[i].end << std::endl;
        std::cout << "junctions[i].strand: " << junctions[i].strand << std::endl;
        std::cout << std::endl;

        CJunc jnew(junctions[i].start, junctions[i].end, junctions[i].strand, junctions[i].ref);
        
        std::cout << "spur_juncs.Exists(jnew):  " << spur_juncs.Exists(jnew) << std::endl;
        if (spur_juncs.Exists(jnew)) {
            // spur = true;
            std::cout << "spur_juncs.Exists! " << std::endl;
            for (int j=0; j<junctions[i].read_ls.size(); j++) {
                std::cout << "~~ SPLAM!" << std::endl;
                outfile_discard->write(junctions[i].read_ls.at(j));
                // delete junctions[i].read_ls.at(j);
            }
            std::cout << "spur_juncs.Exists Done! " << std::endl;
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
                if (spur) {
                    std::cout << "spur_juncs.Exists! " << std::endl;
                    std::cout << "~~ SPLAM!" << std::endl;
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
                    std::cout << "spur_juncs not Exists! " << std::endl;
                    kept_brecs.Add(junctions[i].read_ls.at(j));
                    // std::cout << "~~ Clean!" << std::endl;
                }
            }
            std::cout << "~~~ Done! " << std::endl;
        }

        std::cout << "Done!" << std::endl;
	}

    flushBrec(kept_brecs, hits, outfile_cleaned);
    auto end =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << spur_cnt << " spurious alignments were removed." << std::endl;
    std::cout << "Completed in " << duration.count() << " seconds" << std::endl;
    
    
    delete outfile_discard;
    std::cout << "Done delete outfile_discard!" << std::endl;
    delete outfile_spliced;
    std::cout << "Done delete outfile_spliced!" << std::endl;
    delete outfile_cleaned;
    std::cout << "Done delete outfile_cleaned!" << std::endl;
}