#include <fstream>
#include <iostream>
#include <cctype>
#include <sstream>
#include <filesystem>

#include "predict.h"
#include "extract.h"
#include "common.h"
#include "util.h"
#include "splam_stream.h"

#include <progressbar/progressbar.hpp>
#include <gclib/GBase.h>
#include <htslib/htslib/faidx.h>

#include <Python.h>

GStr splamPredict() {
    /*********************************************
     * Step 1: (1) Iterating through a BAM file / BAM files
     *         (2) Accumulate junctions in a BED file.
    *********************************************/
    GMessage("###########################################\n");
    GMessage("## Step 1: generating spliced junctions in BED\n");
    GMessage("###########################################\n");

    GStr outfname_junc_bed;
    if (!skip_extact) {
        outfname_junc_bed = splamJExtract();
    } else {
        outfname_junc_bed = out_dir + "/junction.bed";
    }
    // std::unordered_map<std::string, int> chrs = get_hg38_chrom_size("HISAT2");
    faidx_t * ref_faidx = fastaIndex();
    // GMessage("[INFO] Predicting ...\n");

    /*********************************************
     * Step 2: (1) getting coordinates of donors and acceptors
     *         (2) Writing FASTA file of donors and acceptors
     *         (3) checking the junctions. (GT-AG ratio)
    *********************************************/
    GMessage("\n###########################################\n");
    GMessage("## Step 2: getting coordinates of donors and acceptors\n");
    GMessage("###########################################\n");

    robin_hdd_hm doner_dimers;
    robin_hdd_hm acceptor_dimers;
    GStr outfname_junc_fa = splamCreateFasta(outfname_junc_bed, doner_dimers, acceptor_dimers, ref_faidx);


    // std::cout << ">> Donor dimers: " << std::endl;
    // for (auto i : doner_dimers) {
    //     std::cout << "\t" << i.first << ": " << i.second << std::endl;
    // }

    // std::cout << ">> Acceptor dimers: " << std::endl;
    // for (auto i : acceptor_dimers) {
    //     std::cout << "\t" << i.first << ": " << i.second << std::endl;
    // }

    std::cout << "[Info] Donor dimers: " << std::endl;
    int cnl_donor_ct=0, ncnl_donor_ct=0, cnl_acceptor_ct=0, ncnl_acceptor_ct = 0;
    for (auto i : doner_dimers) {
        if (std::strcmp(i.first.c_str(), "GT") == 0) {
            cnl_donor_ct += i.second;
        } else {
            ncnl_donor_ct += i.second;
        }
        // std::cout << "\t" << i.first << ": " << i.second << std::endl;
    }
    GMessage("\tCanonical donor\t\t: %d\n", cnl_donor_ct);
    GMessage("\tNoncanonical donor\t: %d\n", ncnl_donor_ct);

    std::cout << "[Info] Acceptor dimers: " << std::endl;
    for (auto i : acceptor_dimers) {
        if (std::strcmp(i.first.c_str(), "AG") == 0) {
            cnl_acceptor_ct += i.second;
        } else {
            ncnl_acceptor_ct += i.second;
        }
        // std::cout << "\t" << i.first << ": " << i.second << std::endl;
    }

    GMessage("\tCanonical acceptor\t: %d\n", cnl_acceptor_ct);
    GMessage("\tNoncanonical donor\t: %d\n", ncnl_acceptor_ct);


    /*********************************************
     * Step 3: SPLAM model prediction
    *********************************************/
    GMessage("\n###########################################\n");
    GMessage("## Step 3: SPLAM model prediction\n");
    GMessage("###########################################\n");
    Py_Initialize();
    // GStr python_f = "./script/splam.py";
    GStr outfname_junc_score(out_dir + "/junction_score.bed");

    wchar_t *argvv[] = {GetWC("."), GetWC("-f"), GetWC(outfname_junc_fa.chars()), GetWC("-o"), GetWC(outfname_junc_score.chars()), GetWC("-m"), GetWC(infname_model_name.chars())};
    PySys_SetArgv(7, argvv);

    PyRun_SimpleString(python_script.c_str());
    Py_Finalize();

    return outfname_junc_score;
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

GStr splamCreateFasta(GStr outfname_junc_bed, robin_hdd_hm &doner_dimers, robin_hdd_hm &acceptor_dimers, faidx_t *ref_faidx) {
    int SEQ_LEN = 800;
    int QUOTER_SEQ_LEN = SEQ_LEN/4;

    /*********************************************
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    *******************************************/
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

    GStr outfname_junc_fa(fa_dir + "/junction.fa");
    GStr outfname_donor_fa(fa_dir + "/donor.fa");
    GStr outfname_accceptor_fa(fa_dir + "/acceptor.fa");

    std::ofstream outfile_fa_junc(outfname_junc_fa);
    std::ofstream outfile_fa_donor(outfname_donor_fa);
    std::ofstream outfile_fa_acceptor(outfname_accceptor_fa);

    std::ifstream fr_junc(outfname_junc_bed);
    std::string line;

    progressbar bar(JUNC_COUNT);
    bar.set_opening_bracket_char("[INFO] SPLAM! Writing junction BED file \n\t[");

    while(getline(fr_junc, line)){
        bar.update();

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

        if (donor_e >= CHRS[chromosome] or acceptor_e >= CHRS[chromosome]) {
            continue;
        }
        if (donor_s < 0 or acceptor_s < 0) {
            continue;
        }

        outfile_bed_donor << chromosome << "\t" + std::to_string(donor_s) + "\t" + std::to_string(donor_e) + "\t" + junc_name+"_donor" + "\t" + std::to_string(num_alignment) + "\t" + strand + "\n";
        outfile_bed_acceptor << chromosome << "\t" + std::to_string(acceptor_s) + "\t" + std::to_string(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + std::to_string(num_alignment) + "\t" + strand + "\n";

        char* donor_seq = NULL;
        char* acceptor_seq = NULL;

        try {
            int donor_len = donor_e - donor_s;
            donor_seq = faidx_fetch_seq(ref_faidx, chromosome.c_str(), donor_s, donor_e-1, &donor_len);
        }
        catch (int a) {
            continue;
        }

        try {
            int acceptor_len = acceptor_e - acceptor_s;
            acceptor_seq = faidx_fetch_seq(ref_faidx, chromosome.c_str(), acceptor_s, acceptor_e-1, &acceptor_len);
        }
        catch (int a) {
            continue;
        }

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


        if (strlen(donor_seq) > 400 || strlen(acceptor_seq) > 400) {
            continue;
        } else {
            outfile_fa_junc << ">" << chromosome << ";" << std::to_string(donor) << ";" << std::to_string(acceptor) << ";" << strand << std::endl;
            
            if (strlen(donor_seq) == 400 &&  strlen(acceptor_seq) == 400) {
                outfile_fa_junc << donor_seq << acceptor_seq << std::endl;
            } else if (strlen(donor_seq) < 400 && strlen(donor_seq) == strlen(acceptor_seq)){
                outfile_fa_junc << donor_seq << std::string(2*(400-(int)strlen(donor_seq)), 'N') << acceptor_seq << std::endl;
            }
        }



        char donor_dim[3];
        memcpy(donor_dim, &donor_seq[200], 2);
        donor_dim[2] = '\0';
        for(int i=0;i<strlen(donor_dim);i++)
            donor_dim[i] = toupper(donor_dim[i]);
        std::string donor_str(donor_dim);

        char acceptor_dim[3];
        memcpy(acceptor_dim, &acceptor_seq[198], 2);
        acceptor_dim[2] = '\0';
        for(int i=0;i<strlen(acceptor_dim);i++)
            acceptor_dim[i] = toupper(acceptor_dim[i]);
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
    GMessage("\n");
    return outfname_junc_fa;
}