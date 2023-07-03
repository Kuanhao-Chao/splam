/*  extract.cpp -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */


#include "extract.h"
#include "common.h"
#include "extract.h"
#include "junc_func.h"
#include "util.h"
#include "bundle.h"
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <gclib/GBase.h>

/****************************
* Input : BAM file (s)
* Output: Junction bed file.
*****************************/
GStr splamJExtract() {
    STEP_COUNTER += 1;

    if (verbose) {
        GMessage("###########################################\n");
        GMessage("## Step %d: generating spliced junctions in BED\n", STEP_COUNTER);
        GMessage("###########################################\n");
    }
    
    GStr outfname_junc_bed = out_dir + "/junction.bed";

    if (verbose) {
        GMessage("[INFO] Extracting junctions ...\n");
        GMessage("[INFO] Processing BAM file ...\n");
    }

    /****************************
    * Creating junction bed files.
    *****************************/
    if (!outfname_junc_bed.is_empty()) {
        if (strcmp(outfname_junc_bed.substr(outfname_junc_bed.length()-4, 4).chars(), ".bed")!=0) {
            outfname_junc_bed.append(".bed");
        }
        joutf = fopen(outfname_junc_bed.chars(), "w");
        if (joutf==NULL) GError("Error creating file %s\n", outfname_junc_bed.chars());
    }

    BundleData* bundle = new BundleData();
	GList<CReadAln> readlist;

    GHash<int> hashread; //read_name:pos:hit_index => readlist index
    GStr lastref;
    bool more_alns = true;
    int prev_pos = 0;
    int lastref_id = -1; //last seen gseq_id

    bool fr_strand = false;
    bool rf_strand = false;
    int currentstart = 0, currentend = 0;
    int bundle_counter = 0;

    /****************************
    * Workflow to extract junctions
    *****************************/
    if (g_paired_removal) {
        while (more_alns) {
            bool chr_changed=false;
            int pos=0;
            const char* refseqName=NULL;
            char xstrand=0;
            int nh=1;
            int gseq_id=lastref_id;  //current chr id
            bool new_bundle=false;

            /***********************************
             * (1) Assessing the current alignment. See if it's good
            ************************************/
            if ((irec=in_records.next())!=NULL) {
                ALN_COUNT ++;
                brec=irec->brec;

                /***********************************
                 * Setting the "chr" "strand" of the current alignment.
                ************************************/
                refseqName=brec->refName();
                xstrand=brec->spliceStrand(); // tagged strand gets priority
                if(xstrand=='.' && (fr_strand || rf_strand)) { // set strand if stranded library
                    if(brec->isPaired()) { // read is paired
                        if(brec->pairOrder()==1) { // first read in pair
                            if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
                            else xstrand='-';
                        }
                        else {
                            if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='-';
                            else xstrand='+';
                        }
                    }
                    else {
                        if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
                        else xstrand='-';
                    }
                }

                /***********************************
                 * Setting the "chr_changed" and "new_bundle" parameters.
                ************************************/
                pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
                chr_changed=(lastref.is_empty() || lastref!=refseqName);
                if (chr_changed) {
                    prev_pos=0;
                }

                if (pos == 0) {
                    // This is an unmapped read => do nothing for now
                } else if (pos<prev_pos) {
                    GMessage("[ERROR] %s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
                    brec->name(), brec->start,  pos, refseqName, prev_pos);
                    exit(-1);
                } else {
                    prev_pos=pos;
                }
                nh=brec->tag_int("NH", 1);
                if (!chr_changed && currentend>0 && pos>currentend+g_bundle_gap) {
                    GMessage("New bundle: currentend(%d - %d), %d, g_bundle_gap: %d\n ", currentend, currentend, pos, g_bundle_gap);
                    new_bundle=true;
                }
            } else { //no more alignments
                more_alns=false;
                new_bundle=true; //fake a new start (end of last bundle)
            }


            /***********************************
             * (2) Process the bundle!
            ************************************/
            if (new_bundle || chr_changed) {
                if (currentstart < currentend) {
                    GMessage("current boundaries: %d - %d\n", currentstart, currentend);
                }
                hashread.Clear();
                if (readlist.Count()>0) {
                    // process reads in previous bundle
                    bundle->getReady(currentstart, currentend);
                    processBundle_jext(bundle, readlist, bundle_counter);
                    readlist.Clear();
                } else { 
                    //no read alignments in this bundle?  
                    bundle->Clear();
                    readlist.Clear();
                } //nothing to do with this bundle

                if (chr_changed) {
                    lastref = refseqName;
                    lastref_id = gseq_id;
                    currentend = 0;
                }

                if (!more_alns) {
                    // noMoreBundles();
                    break;
                }

                if (brec->start > 0) {
                    currentstart = pos;
                    currentend = brec->end;
                }
                bundle->refseq = lastref;
                bundle->start = currentstart;
                bundle->end = currentend;
            } //<---- new bundle started

            int fragment_end = 0;
            if (brec->refId() == brec->mate_refId()) {

                int insert_size = brec->insertSize();
                int mate_end = brec->mate_start() + insert_size;
                
                if (mate_end > (int)brec->end) {
                    fragment_end = mate_end;
                } else {
                    fragment_end = (int)brec->end;
                }
            } else {
                fragment_end = (int)brec->end;
            }

            if (currentend<fragment_end) {
                //current read extends the bundle => adjusted currentend
                currentend=fragment_end;
            }

            /***********************************
             * (3) Process the alignment!
            ************************************/
            CReadAln* alndata = new CReadAln(brec);
            processRead_jext(currentstart, currentend, readlist, *bundle, hashread, alndata);
        } //for each read alignment
    } else {
        int prev_tid=-1;
        GStr prev_refname;
        int b_end=0;
        int Read_counter = 0;
        while ((irec=in_records.next())!=NULL) {
            Read_counter += 1;
            ALN_COUNT += 1;
            brec=irec->brec;
            int endpos=brec->end;
            if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
                flushJuncs(joutf);
                junctions.setCapacity(128);
                b_end=endpos;
                prev_tid=brec->refId();
                prev_refname=(char*)brec->refName();
            } else { //extending current bundle
                if (b_end<endpos) {
                    b_end=endpos;
                }
            }

            if (joutf && brec->exons.Count()>1) {
                /*************************
                 * Spliced reads
                *************************/
                int accYC = 0;
                accYC = brec->tag_int("YC", 1);
                addJunction(*brec, accYC, prev_refname);
                if (write_bam) {
                    int new_nh = brec->tag_int("NH", 0);
                    if (new_nh <= 1) {
                        outfile_s_uniq_map->write(brec);
                        ALN_COUNT_SPLICED_UNIQ += 1;
                    } else if (new_nh > 1) {
                        outfile_s_multi_map->write(brec);
                        ALN_COUNT_SPLICED_MULTI += 1;
                    }
                }
            } else {
                /*************************
                 * Non-spliced reads
                *************************/
                // Non-spliced reads.
                if (write_bam) {
                    int new_nh = brec->tag_int("NH", 0);
                    if (new_nh <= 1) {
                        outfile_ns_uniq_map->write(brec);
                        ALN_COUNT_NSPLICED_UNIQ += 1;
                        ALN_COUNT_GOOD += 1;
                    } else if (new_nh > 1){
                        outfile_ns_multi_map->write(brec);
                        ALN_COUNT_NSPLICED_MULTI += 1;
                    }
                }
            }
            // if (ALN_COUNT % 1000000 == 0) {
            //     GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
            // }
        }

        GMessage("[Info] Total number of alignments processed: %d!\n", Read_counter);
    }

    in_records.stop();
    flushJuncs(joutf);
    fclose(joutf);
    junctions.setCapacity(128);
    if (write_bam) {
        delete outfile_ns_multi_map;
        delete outfile_ns_uniq_map;
        delete outfile_s_multi_map;
        delete outfile_s_uniq_map;
        if (g_paired_removal) {
            delete outfile_ns_multi_unpair;
            delete outfile_ns_uniq_unpair;
            delete outfile_s_multi_unpair;
            delete outfile_s_uniq_unpair;
        }
    }
    return outfname_junc_bed;
}




void processBundle_jext(BundleData* bundle, GList<CReadAln>& readlist, int& bundle_counter) {
    bundle_counter += 1;

    if (verbose) {
        GMessage("\t* In bundle %d (%s)\n", bundle_counter, readlist[0]->brec.refName());
        GMessage("\t\t>> bundle read count: %d\n", readlist.Count());
        GMessage("\t\t>> bundle start     : %d\n", bundle->start);
        GMessage("\t\t>> bundle end       : %d\n", bundle->end);
    }

    robin_hdd_int hash_processed;

    // Writing out junctions.
    flushJuncs(joutf);
    junctions.Clear();
    junctions.setCapacity(128);

    // Writing out alignments
    int Read_counter = 0;
    for (int idx=0; idx<readlist.Count(); idx++) {
        Read_counter += 1;
        GSamRecord& brec_bd = readlist[idx]->brec;

        if (hash_processed.find(idx) != hash_processed.end()) {
            // The read has already been processed.
            continue;
        }

        bool remove_algn = false;
        int pair_idx = readlist[idx]->pair_idx;
    
        /***********************************
         * Processing unpaired reads.
         *   => Write it out into the unpaired BAM.
        ************************************/
        // Check the global hash => only used when its mate is unpaired.
        if (pair_idx == -1) {
            int brec_tag = brec_bd.tag_int("NH", 0);
            if (brec_bd.hasIntrons()) {
                if (brec_tag > 1)  {
                    outfile_s_multi_unpair->write(&brec_bd);
                    ALN_COUNT_SPLICED_MULTI_UNPAIR += 1;
                } else {
                    outfile_s_uniq_unpair->write(&brec_bd);
                    ALN_COUNT_SPLICED_UNIQ_UNPAIR += 1;
                }
            } else {
                if (brec_tag > 1)  {
                    outfile_ns_multi_unpair->write(&brec_bd);
                    ALN_COUNT_NSPLICED_MULTI_UNPAIR += 1;
                } else {
                    outfile_ns_uniq_unpair->write(&brec_bd);
                    ALN_COUNT_NSPLICED_UNIQ_UNPAIR += 1;
                    ALN_COUNT_GOOD += 1;
                }
            }
            continue;
        }
        // fprintf(stderr, ">> idx: %d;  pair_idx: %d\n", idx, pair_idx);

        /***********************************
         * Processing paired reads.
         *   => Write it out into the paired BAM.
        ************************************/
        // Process both the current and paired alginments.
        GSamRecord& brec_bd_p = readlist[readlist[idx]->pair_idx]->brec;

        int brec_bd_tag = brec_bd.tag_int("NH", 0);
        int brec_bd_p_tag = brec_bd_p.tag_int("NH", 0);

        if (write_bam) {
            if ((!brec_bd.hasIntrons() && !brec_bd_p.hasIntrons())) {
                /***********************************
                 * Both alignment does not have introns.
                ************************************/
                if ((brec_bd_tag==1 || brec_bd_p_tag==1)) {
                    // a, b nonspliced, NH == 1
                    outfile_ns_uniq_map->write(&brec_bd);
                    ALN_COUNT_NSPLICED_UNIQ +=1;
                } else if ((brec_bd_tag>1 || brec_bd_p_tag>1)) {
                    // a, b nonspliced, NH > 1
                    outfile_ns_multi_map->write(&brec_bd);
                    ALN_COUNT_NSPLICED_MULTI+=1;
                }
                if (pair_idx < idx) {
                    brec_bd.clear();
                    brec_bd_p.clear();
                } 
            } else if ((brec_bd.hasIntrons() || brec_bd_p.hasIntrons()))  {
                /***********************************
                 * Any of the alignment has introns.
                ************************************/
                if ((brec_bd_tag==1 || brec_bd_p_tag==1)) {
                    // a, b spliced, NH = 1
                    outfile_s_uniq_map->write(&brec_bd);
                    outfile_s_uniq_map->write(&brec_bd_p);
                    ALN_COUNT_SPLICED_UNIQ+=2;
                } else if ((brec_bd_tag>1 || brec_bd_p_tag>1)) {
                    // a, b spliced, NH > 1
                    outfile_s_multi_map->write(&brec_bd);
                    outfile_s_multi_map->write(&brec_bd_p);
                    ALN_COUNT_SPLICED_MULTI+=2;
                }
                hash_processed.insert(pair_idx);
                brec_bd.clear();
                brec_bd_p.clear();
            }
            // else {
            //     // Execption case
            //     GMessage("There are wierd things happening! Check this alignment.\n");
            //     outfile_ns_uniq_map->write(&brec_bd);
            //     outfile_ns_uniq_map->write(&brec_bd_p);
            //     ALN_COUNT_GOOD += 2;
            // }
        }
    }

    GMessage("\t* # alignments processed: %d!\n", Read_counter);
    bundle->Clear();
}

void processRead_jext(int currentstart, int currentend, GList<CReadAln>& readlist, BundleData& bdata, GHash<int>& hashread, CReadAln* alndata) { 
	
    GSamRecord& brec=(alndata->brec);
    static GStr _id("", 256); //to prevent repeated reallocation for each parsed read
    static GStr _id_p("", 256); //to prevent repeated reallocation for each parsed read
	
    /*
	{ // DEBUG ONLY
		fprintf(stderr,"Process read %s with exons:", brec.name());
		for (int i=0;i<brec.exons.Count();i++) {
			fprintf(stderr," %d-%d", brec.exons[i].start, brec.exons[i].end);
		}
		fprintf(stderr,"\n");
	}
	*/

	bool match=false;  // true if current read matches a previous read
	int n = 0;
	if (bdata.end<currentend) {
		bdata.start=currentstart;
		bdata.end=currentend;
	}
	bdata.numreads++;  // number of reads gets increased no matter what

    /************************
    * Adding splice sites
    ************************/
    if (joutf && brec.exons.Count()>1) {
        // Spliced reads
        int accYC = 0;
        accYC = brec.tag_int("YC", 1);
        addJunction(brec, accYC, brec.refName());
        ALN_COUNT_SPLICED++;
    }

    n = readlist.Add(alndata); // reset n for the case there is no match

	if((int)brec.end>currentend) {
		currentend=brec.end;
	  	bdata.end=currentend;
	}

	// now set up the pairing
	if (brec.refId()==brec.mate_refId()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired
        int self_start = brec.start;
		int pair_start = brec.mate_start();
        int insert_size = brec.insertSize();
        int pair_insert_size = (-1)*insert_size;
        int pair_idx = brec.pairOrder();
        if (brec.isUnmapped()) {
            self_start = pair_start;
        } else if (brec.isMateUnmapped()) {
            pair_start = self_start;
        }

        // GMessage("brecname: %s; self_start: %d;  pair_start: %d;  currentstart: %d; insert_size: %d; pair_insert_size:%d; pair_idx:%d\n", brec.name(), self_start, pair_start, currentstart, insert_size, pair_insert_size, pair_idx);
		if (currentstart<=pair_start) { // if pair_start is in a previous bundle I don't care about it
			//GStr readname();
			//GStr id(brec->name(), 16); // init id with readname
			_id.assign(brec.name()); //assign can be forced to prevent shrinking of the string
            _id_p.assign(brec.name());
            _id+=';';_id+=self_start;
            _id+=';';_id+=pair_start;
            _id+=';';_id+=insert_size;
            _id+=';';_id+=pair_idx;

			_id_p+=';';_id_p+=pair_start;
            _id_p+=';';_id_p+=self_start;
            _id_p+=';';_id_p+=pair_insert_size;
            _id_p+=';';_id_p+=(3-pair_idx);

            int* n_check=hashread[_id.chars()];
            while (n_check) {
                // GMessage("\tChecking repeat: %s;  %d\n", _id.chars(), *n_check);
                _id+='*';
                _id_p+='*';
                n_check=hashread[_id.chars()];
            }

			if(pair_start < self_start) { // if I've seen the pair already <- I might not have seen it yet because the pair starts at the same place
				const int* np=hashread[_id_p.chars()];
                if (np) {
                    // GMessage("\t\tn : %d\n", n);
                    // GMessage("\t\tnp: %d\n", *np);
                    readlist[*np]->pair_idx = n;
                    readlist[n]->pair_idx = *np;
                } else {
                    // GMessage(">> Pair not in the same bundle\n");
                }

            } else if (pair_start == self_start) {
				hashread.Add(_id.chars(), n);
				const int* np=hashread[_id_p.chars()];
                // GMessage("\t_equal condition\n");
                // GMessage("\t\tn : %d\n\n", n);
                // GMessage("\t\tnp: %d\n\n", *np);                    

                if (np) {
                    // GMessage("\t\tn : %d\n", n);
                    // GMessage("\t\tnp: %d\n", *np);
                    readlist[*np]->pair_idx = n;
                    readlist[n]->pair_idx = *np;
                    // hashread.Remove(_id.chars());
                    // hashread.Remove(_id_p.chars());
                }
            } // Do not process the pair that is larger.
            hashread.Add(_id.chars(), n);
            const int* n_check_final=hashread[_id.chars()];
            // GMessage("Retrieving read from hash: %s;  %d\n", _id.chars(), *n_check_final);
		}
	}
}