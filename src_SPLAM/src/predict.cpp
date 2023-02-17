#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>

#include "predict.h"
#include "extract.h"
#include "common.h"
#include "util.h"
#include "splam_stream.h"

#include <gclib/GBase.h>
#include <htslib/htslib/faidx.h>
#include <robin_hood/robin_hood.h>

#include <Python.h>

typedef robin_hood::unordered_map<std::string, int> dimer_hm;

void splamPredict() {
    splamJExtract();

    std::unordered_map<std::string, int> chrs = get_hg38_chrom_size("HISAT2");
    faidx_t * ref_faidx = fastaIndex();
    GMessage("[INFO] Predicting ...\n");


    /*********************************************
     * Step 2: (1) getting coordinates of donors and acceptors
     *         (2) Writing FASTA file of donors and acceptors
     *         (3) checking the junctions. (GT-AG ratio)
    *********************************************/
    GMessage("********************************************\n");
    GMessage("** Step 2: getting coordinates of donors and acceptors\n");
    GMessage("********************************************\n\n");

    int SEQ_LEN = 800;
    int QUOTER_SEQ_LEN = SEQ_LEN/4;

    /*********************************************
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    *******************************************/
    dimer_hm doner_dimers;
    dimer_hm acceptor_dimers;

    GStr bed_dir(out_dir + "/bed");
    GStr fa_dir(out_dir + "/fasta");

    std::filesystem::create_directories(bed_dir.chars());
    std::filesystem::create_directories(fa_dir.chars());

    GStr donor_bed(bed_dir + "/donor.bed");
    GStr acceptor_bed(bed_dir + "/acceptor.bed");
    GStr da_bed(bed_dir + "/d_a.bed");


    std::ofstream outfile_bed_donor(donor_bed);
    std::ofstream outfile_bed_acceptor(acceptor_bed);
    std::ofstream outfile_bed_da(da_bed);

    GStr junc_fasta(fa_dir + "/junction.fa");
    GStr donor_fasta(fa_dir + "/donor.fa");
    GStr accceptor_fasta(fa_dir + "/acceptor.fa");

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

        // if (donor_e >= chrs[chromosome] or acceptor_e >= chrs[chromosome]) {
        //     std::cout << "Skip!!" << std::endl;
        //     continue;
        // }
        // if (donor_s < 0 or acceptor_s < 0) {
        //     std::cout << "Skip!!" << std::endl;
        //     continue;
        // }


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
    GMessage("********************************************\n");
    GMessage("** Step 3: SPLAM model prediction\n");
    GMessage("********************************************\n");
    Py_Initialize();
    // GStr python_f = "./script/splam.py";
    GStr outfile_junction_score(out_dir + "/junction_score.bed");
    
    // wchar_t *argvv[] = {GetWC(python_f.chars()), GetWC("-f"), GetWC(junc_fasta.chars()), GetWC("-o"), GetWC(outfile_junction_score.chars()), GetWC("-m"), GetWC(infname_model_name.chars())};
    // PySys_SetArgv(7, argvv);


    wchar_t *argvv[] = {GetWC("./script/splam.py"), GetWC("-f"), GetWC(junc_fasta.chars()), GetWC("-o"), GetWC(outfile_junction_score.chars()), GetWC("-m"), GetWC(infname_model_name.chars())};
    PySys_SetArgv(7, argvv);

    GMessage("infname_model_name: %s\n", infname_model_name.chars());
    PyRun_SimpleString(python_script);
    // FILE *file = _Py_fopen(python_f.chars(), "r");
    // if(file != NULL) {
    //     int re = PyRun_SimpleFileExFlags(file, python_f.chars(), false, NULL);
    // }
    Py_Finalize();

}


faidx_t *fastaIndex() {
    FILE *ref_fa_f;
    if ((ref_fa_f = fopen(infname_reffa.chars(), "r"))) {
        fclose(ref_fa_f);
        GMessage("[INFO] FASTA index \"%si\" has been created!\n", infname_reffa.chars());
    } else {
        int res = fai_build(infname_reffa.chars());
        GMessage("[INFO] Creating FASTA index \"%si\"\n", infname_reffa.chars());
    }

    faidx_t *ref_faidx = fai_load(infname_reffa.chars());
    GMessage("[INFO] Loading FASTA file\n");
    return ref_faidx;
}