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
    
    // This is normal workflow for writing out all junctions.
    // outfile_multimapped = new GSamWriter(outfname_multimapped, in_records.header(), GSamFile_BAM);
    
    GStr outfname_junc_bed = out_dir + "/junction.bed";
    // outfile_spliced = new GSamWriter(outfname_spliced, in_records.header(), GSamFile_BAM);

    if (verbose) {
        GMessage("[INFO] Extracting junctions ...\n");
    }
    // GMessage("[INFO] Output directory\t\t: %s\n", out_dir.chars());
    // GMessage("[INFO] Output Junction file\t: %s\n", outfname_junc_bed.chars());

    /****************************
    * Creating junction bed files.
    *****************************/
    // Creating the output junction bed file
    if (!outfname_junc_bed.is_empty()) {
        if (strcmp(outfname_junc_bed.substr(outfname_junc_bed.length()-4, 4).chars(), ".bed")!=0) {
            outfname_junc_bed.append(".bed");
        }
        joutf = fopen(outfname_junc_bed.chars(), "w");
        if (joutf==NULL) GError("Error creating file %s\n", outfname_junc_bed.chars());
    }



    // outfile_cleaned = new GSamWriter(outfname_cleaned, in_records.header(), GSamFile_BAM);
    // outfile_ns_multi_map = new GSamWriter(outfname_ns_multi_map, in_records.header(), GSamFile_BAM);
    // outfile_s_uniq_map = new GSamWriter(outfname_s_uniq_map, in_records.header(), GSamFile_BAM);
    // outfile_s_multi_map = new GSamWriter(outfname_s_multi_map, in_records.header(), GSamFile_BAM);
    // outfile_s_multi_map_tmp = new GSamWriter(outfname_s_multi_map_tmp, in_records.header(), GSamFile_BAM);
    // outfile_discard_unpair = new GSamWriter(outfname_discard_unpair, in_records.header(), GSamFile_BAM);
    // outfile_discard_s_uniq_map= new GSamWriter(outfname_discard_s_uniq_map, in_records.header(), GSamFile_BAM);
    // outfile_discard_s_multi_map= new GSamWriter(outfname_discard_s_multi_map, in_records.header(), GSamFile_BAM);

    
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

    while (more_alns) {
        bool chr_changed=false;
        int pos=0;
        const char* refseqName=NULL;
        char xstrand=0;
        int nh=1;
        int hi=0;
        int gseq_id=lastref_id;  //current chr id
        bool new_bundle=false;
        //delete brec;
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
                // This is an unmapped read
            } else if (pos<prev_pos) {
                GMessage("[ERROR] %s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
                brec->name(), brec->start,  pos, refseqName, prev_pos);
                exit(-1);
            }
            prev_pos=pos;
            nh=brec->tag_int("NH", 0);
            if (nh==0) nh=1;
            hi=brec->tag_int("HI", 0);
            if (!chr_changed && currentend>0 && pos>currentend+g_max_splice) {
                new_bundle=true;
            }
        } else { //no more alignments
            more_alns=false;
            new_bundle=true; //fake a new start (end of last bundle)
        }

        /***********************************
         * Process the bundle!
        ************************************/
        if (new_bundle || chr_changed) {
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
            
            if (mate_end > (int)brec->end && insert_size <= g_max_splice) {
                fragment_end = mate_end;
            } else {
                fragment_end = (int)brec->end;
            }
        } else {
            fragment_end = (int)brec->end;
        }

        if (currentend<fragment_end) {
            //current read extends the bundle
            currentend=fragment_end;
        } //adjusted currentend and checked for overlapping reference transcripts

        // GMessage("brec->refName(): %s\n", brec->refName());
        CReadAln* alndata = new CReadAln(brec);
        processRead_jext(currentstart, currentend, readlist, *bundle, hashread, alndata);
    } //for each read alignment






































    // // This is additional workflow to write out junctions above & below the thresholds.
    // if (g_j_extract_threshold > 0) {
    //     GStr outfname_junc_above_bed = out_dir + "/junction_above.bed";
    //     GStr outfname_junc_below_bed = out_dir + "/junction_below.bed";

    //     outfile_above_spliced = new GSamWriter(outfname_junc_above_bed, in_records.header(), GSamFile_BAM);
    //     outfile_below_spliced = new GSamWriter(outfname_junc_below_bed, in_records.header(), GSamFile_BAM);

    //     GMessage("[INFO] Extracting junctions ...\n");
    //     GMessage("[INFO] Output directory\t\t: %s\n", out_dir.chars());
    //     GMessage("[INFO] Output Junction file\t: %s; %s\n", outfname_junc_above_bed.chars(), outfname_junc_below_bed.chars());

    //     // Creating the output junction bed file
    //     if (!outfname_junc_above_bed.is_empty()) {
    //         if (strcmp(outfname_junc_above_bed.substr(outfname_junc_above_bed.length()-4, 4).chars(), ".bed")!=0) {
    //             outfname_junc_above_bed.append(".bed");
    //         }
    //         joutf_above = fopen(outfname_junc_above_bed.chars(), "w");
    //         if (joutf_above==NULL) GError("Error creating file %s\n", outfname_junc_above_bed.chars());
    //     }
    //     if (!outfname_junc_below_bed.is_empty()) {
    //         if (strcmp(outfname_junc_below_bed.substr(outfname_junc_below_bed.length()-4, 4).chars(), ".bed")!=0) {
    //             outfname_junc_below_bed.append(".bed");
    //         }
    //         joutf_below = fopen(outfname_junc_below_bed.chars(), "w");
    //         if (joutf_below==NULL) GError("Error creating file %s\n", outfname_junc_below_bed.chars());
    //     }
    // }

    // /****************************
    // * Iterating BAM file(s) and write out junctions.
    // *****************************/
    // // Reading BAM file.
    // int prev_tid=-1;
    // GStr prev_refname;
    // int b_end=0, b_start=0;

    // GMessage("[INFO] Processing BAM file ...\n");
    // // GMessage("\t\tBefore Hash map size: %d\n", read_hashmap.size());
    // while ((irec=in_records.next())!=NULL) {
    //     brec=irec->brec;
    //     int endpos=brec->end;
    //     if (brec->refId()!=prev_tid || (int)brec->start>b_end) {
    //         flushJuncs(joutf);
    //         if (g_j_extract_threshold > 0) {
    //             flushJuncs(joutf_above, joutf_below);
    //         }
    //         junctions.Clear();
    //         junctions.setCapacity(128);
    //         b_start=brec->start;
    //         b_end=endpos;
    //         prev_tid=brec->refId();
    //         prev_refname=(char*)brec->refName();
    //     } else { //extending current bundle
    //         if (b_end<endpos) {
    //             b_end=endpos;
    //         }
    //     }
    //     int accYC = 0;
    //     accYC = brec->tag_int("YC", 1);
    //     if (joutf && brec->exons.Count()>1) {
    //         // Spliced reads
    //         addJunction(*brec, accYC, prev_refname);
    //         // outfile_spliced->write(brec);
    //         ALN_COUNT_SPLICED++;
    //     } else {
    //         // Non-spliced reads.
    //         // Not spliced => check their NH tags!
    //         if (brec->isUnmapped()) continue;
    //         int new_nh = brec->tag_int("NH", 0);
    //         if (new_nh == 1) {
    //             outfile_cleaned->write(brec);
    //         } else if (new_nh == 0){
    //             GMessage("\t\t brec->name(): %s !!!\n", brec->name());
    //             GMessage("\t\t NH tag is zero !!!: %d\n", new_nh);
    //         } else {
    //             // outfile_multimapped->write(brec);
    //             ALN_COUNT_NH_UPDATE++;
    //         }
    //         ALN_COUNT_NSPLICED++;
    //     }
    //     ALN_COUNT++;
    //     if (ALN_COUNT % 1000000 == 0) {
    //         GMessage("\t\t%d alignments processed.\n", ALN_COUNT);
    //     }
    // }
    // // GMessage("\t\tAfter Hash map size: %d\n", read_hashmap.size());
    // // for (auto it : read_hashmap) {
    // //     std::cout << " " << it.first << ":" << "(NH tag) " << it.second.NH_tag_bound << std::endl;
    // //     for (int i=0; i<it.second.sam_list.Count(); i++) {
    // //         std::cout << it.second.sam_list[i].name() << std::endl;
    // //     }
    // // }

    if (verbose) {
        GMessage("[INFO] SPLAM! %d alignments processed.\n", ALN_COUNT);
    }
    in_records.stop();
    flushJuncs(joutf);
    // if (g_j_extract_threshold > 0) {
    //     flushJuncs(joutf_above, joutf_below);
    // }
    fclose(joutf);
    if (g_j_extract_threshold > 0) {
        fclose(joutf_above);
        fclose(joutf_below);
    }
    junctions.Clear();
    junctions.setCapacity(128);
    // delete outfile_spliced;
    if (verbose) {
        GMessage("[INFO] SPLAM! Total number of junctions: %d\n", JUNC_COUNT);	
    }
    


    // delete outfile_cleaned;
    
    delete outfile_ns_multi_map;
    delete outfile_s_uniq_map;
    delete outfile_s_multi_map;
    delete outfile_discard_unpair;

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
    for (int idx=0; idx<readlist.Count(); idx++) {
        // GMessage("\tidx %d\n", idx);
        GSamRecord& brec_bd = readlist[idx]->brec;
    
        if (joutf && brec_bd.hasIntrons()) {
            // Spliced reads
            int accYC = 0;
            accYC = brec_bd.tag_int("YC", 1);
            // GMessage("Add junction!\n");
            addJunction(brec_bd, accYC, brec_bd.refName());
            // outfile_spliced->write(brec);
            ALN_COUNT_SPLICED++;
        }


        if (hash_processed.find(idx) != hash_processed.end()) {
            // The read has already been processed.
            continue;
        }

        bool remove_algn = false;
        int pair_idx = readlist[idx]->pair_idx;
        // GMessage("\t\tpair_idx %d\n", pair_idx);

        // GMessage("idx      : %d\n", idx);
        // GMessage("pair_idx : %d\n", pair_idx);
    

        /***********************************
         * Cannot find its pair.
         *  => Write it out into the unpaired BAM.
        ************************************/
        // Check the global hash => only used when its mate is unpaired.
        if (pair_idx == -1) {
            // This read is unpaired.
            // Check the NH tag hit.
            outfile_discard_unpair->write(&brec_bd);
            ALN_COUNT_UNPAIRED++;
            continue;
        }
        hash_processed.insert(pair_idx);

        /***********************************
         * Find its pair.
        ************************************/
        // Process both the current and paired alginments.
        GSamRecord& brec_bd_p = readlist[readlist[idx]->pair_idx]->brec;

        /***********************************
         * Find its pair.
        ************************************/
        int brec_bd_tag = brec_bd.tag_int("NH", 0);
        int brec_bd_p_tag = brec_bd_p.tag_int("NH", 0);

        // GMessage("\tbrec_bd_tag   %d\n", brec_bd_tag);
        // GMessage("\tbrec_bd_p_tag %d\n", brec_bd_p_tag);
        // GMessage("\tbrec_bd_refName   %s\n", brec_bd.name());
        // GMessage("\tbrec_bd_p_refName %s\n", brec_bd_p.name());

        if ( (!brec_bd.hasIntrons() && !brec_bd_p.hasIntrons()) && (brec_bd_tag==1 || brec_bd_p_tag==1)) {
            // a, b nonspliced, NH == 1
            outfile_cleaned->write(&brec_bd);
            outfile_cleaned->write(&brec_bd_p);
            ALN_COUNT_NSPLICED_UNIQ +=2;
            ALN_COUNT_GOOD += 2;
        } else if ( (!brec_bd.hasIntrons() && !brec_bd_p.hasIntrons()) && (brec_bd_tag>1 || brec_bd_p_tag>1)) {
            // a, b nonspliced, NH > 1
            outfile_ns_multi_map->write(&brec_bd);
            outfile_ns_multi_map->write(&brec_bd_p);
            ALN_COUNT_NSPLICED_MULTI+=2;
        } else if ( (brec_bd.hasIntrons() || brec_bd_p.hasIntrons()) && (brec_bd_tag==1 || brec_bd_p_tag==1)) {
            // a, b spliced, NH = 1
            outfile_s_uniq_map->write(&brec_bd);
            outfile_s_uniq_map->write(&brec_bd_p);
            ALN_COUNT_SPLICED_UNIQ+=2;
        } else if ( (brec_bd.hasIntrons() || brec_bd_p.hasIntrons()) && (brec_bd_tag>1 || brec_bd_p_tag>1)) {
            // a, b spliced, NH > 1
            outfile_s_multi_map->write(&brec_bd);
            outfile_s_multi_map->write(&brec_bd_p);
            ALN_COUNT_SPLICED_MULTI+=2;
        } else {
            // Execption case
            outfile_discard_unpair->write(&brec_bd);
            outfile_discard_unpair->write(&brec_bd_p);
            ALN_COUNT_UNPAIRED+=2;
        }
        brec_bd.clear();
        brec_bd_p.clear();
    }
    bundle->Clear();
}

void processRead_jext(int currentstart, int currentend, GList<CReadAln>& readlist, BundleData& bdata, GHash<int>& hashread, CReadAln* alndata) { // some false positives should be eliminated here in order to break the bundle

	GSamRecord& brec=(alndata->brec);			   // bam record
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

	double nm=(double)brec.tag_int("NM"); // read mismatch
	float unitig_cov = unitig_cov=brec.tag_float("YK");
	bool match=false;  // true if current read matches a previous read
	int n = 0;
	if (bdata.end<currentend) {
		bdata.start=currentstart;
		bdata.end=currentend;
	}
	bdata.numreads++;  // number of reads gets increased no matter what

    // GMessage("brec->refName(): %s\n", brec->refName());
    // CReadAln* readaln = new CReadAln(&brec);
    
    // // for (int i=0;i<brec->exons.Count();i++) {
    // //     readaln->len+=brec->exons[i].len();
    // //     if(i) {
    // //         int jstrand=strand;
    // //         uint jstart=brec->exons[i-1].end;
    // //         uint jend=brec->exons[i].start;
    // //     }
    // //     readaln->segs.Add(brec->exons[i]);
    // // }
    n=readlist.Add(alndata); // reset n for the case there is no match

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
                // GMessage("\t old n: %d;\n", n);
                // readlist.Remove(alndata);
                // n = n-1;
                // GMessage("\t new n: %d \n", n);
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
                // hashread.Remove(_id_p.chars());

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
            // GMessage("Adding read to hash: %s;  %d\n", _id.chars(), n);
            hashread.Add(_id.chars(), n);
            const int* n_check_final=hashread[_id.chars()];
            // GMessage("Retrieving read from hash: %s;  %d\n", _id.chars(), *n_check_final);
		}
	} else {

    }
}