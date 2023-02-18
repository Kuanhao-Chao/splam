#include "clean.h"
#include "common.h"
#include "junc.h"
#include "junc_func.h"
#include "util.h"
#include "extract.h"
#include "predict.h"

#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <htslib/htslib/faidx.h>
#include <Python.h>
#include <gclib/GStr.h>
#include <robin_hood/robin_hood.h>

void splamClean(int argc, char* argv[]) {
    GStr outfname_junc_score = splamPredict();

    /*********************************************
     * Step 4: SPLAM filtering out reads.
    *********************************************/
    GMessage("\n********************************************\n");
    GMessage("** Step 4: SPLAM filtering out reads\n");
    GMessage("********************************************\n\n");
    // GSamReader bamreader(inbamname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    // outfile=new GSamWriter(outfname, bamreader.header(), GSamFile_BAM);
    // GSamRecord brec;

    auto start=std::chrono::high_resolution_clock::now();
    int spur_cnt = 0;

    GArray<CJunc> spur_juncs;
    GVec<GSamRecord*> kept_brecs;
    std::unordered_map<std::string, int> hits;
    // GStr inbedname(out_dir + "/output/junc_scores.bed");
    loadBed(outfname_junc_score, spur_juncs);

    TInputFiles in_records_clean;

    for (int i=0;i<in_records.freaders.Count();++i) {
        in_records_clean.addFile(in_records.freaders[i]->fname.chars());
    }

    // in_records_clean.setup(VERSION, argc, argv);
    for (int i=0;i<in_records_clean.freaders.Count();++i) {
        GMessage("%s\n", in_records_clean.freaders[i]->fname.chars());
    }

    int num_samples=in_records_clean.start();
    GMessage("num_samples: %d\n", num_samples);



    int counter = 0;
    int prev_tid=-1;
    GStr prev_refname;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    int b_end=0;
    int b_start=0;





    dimer_hm removed_rd_ls;

    while ((irec=in_records_clean.next())!=NULL) {
        brec=irec->brec;
        uint32_t dupcount=0;
        std::vector<int> cur_samples;
        int endpos=brec->end;

        // std::cout << brec->refId() << std::endl;
        // printf(">> %d \n", brec->refId());

        // Message(">> brec->flags(): %d \n", brec->flags());
        // if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
        //     // if (joutf) {



        //     // } // TODO: write the last column to 3 dec places
        //     b_start=brec->start;
        //     b_end=endpos;
        //     prev_tid=brec->refId();
        //     prev_refname=(char*)brec->refName();
        // } else { //extending current bundle
        //     if (b_end<endpos) {
        //         b_end=endpos;


        //     }
        // }

        // spur_juncs

        int r_exon_count = brec->exons.Count();


        bool spur = false;
        if (r_exon_count > 1) {
            for (int e=1; e<r_exon_count; e++) {
                char strand = '+';
                if (brec->revStrand()) strand = '-';
		        // CJunc j(r.exons[i-1].end+1, r.exons[i].start-1, strand, ref,

                CJunc jnew_sub(brec->exons[e-1].end+1, brec->exons[e].start-1, strand, GStr(brec->refName()));

			    // std::cout << " >> @@ == : " << brec->refName() << ";  strand: " << strand << ";  start: " << brec->exons[e-1].end+1 << ";  end: " << brec->exons[e].start-1 << std::endl;

                if (spur_juncs.Exists(jnew_sub)) {
                    spur = true;
                    break;
                }
            }
        }
        if (spur) {
            // std::cout << "spur_juncs.Exists! " << std::endl;
            spur_cnt++;
            std::cout << "~~ SPLAM!" << std::endl;
            std::string kv = brec->name();
            std::string tmp = std::to_string(brec->pairOrder());
            kv += ";";
            kv += tmp;

            if (removed_rd_ls.find(kv) == removed_rd_ls.end()) {
                removed_rd_ls[kv] = 1;
            } else {
                int val = removed_rd_ls[kv];
                val++;
                removed_rd_ls[kv] = val;
            }

            // outfile_discard->write(junctions[i].read_ls.at(j));
            // delete junctions[i].read_ls.at(j);
        } else {
            // std::cout << "spur_juncs not Exists! " << std::endl;
            // kept_brecs.Add(junctions[i].read_ls.at(j));
            // std::cout << "~~ Clean!" << std::endl;
        }
    }


    // for (auto const& x : removed_rd_ls) {
    //     std::cout << x.first  // string (key)
    //             << ':' 
    //             << x.second // string's value 
    //             << std::endl;
    // }

    in_records_clean.stop();
    // flushJuncs(joutf);
    // fclose(joutf);

    
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
    //             }
    //         }
    //         std::cout << "~~~ Done! " << std::endl;
    //     }

    //     std::cout << "Done!" << std::endl;
	// }

    // flushBrec(kept_brecs, hits, outfile_cleaned);
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