/*  tmerge.cpp

    Copyright (C) 2016 Mihaela Pertea & Geo Pertea

    Author: Mihaela Pertea <mpertea@jhu.edu>
    Author: Geo Pertea

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include "tmerge.h"
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cstring>

bool check_id(std::string& line, std::string id_tag){
    std::stringstream *line_stream = new std::stringstream(line);
    std::string col;

    // first tag should have already been checked
    std::getline(*line_stream,col,'\t');

    // check if ID == SAMPLE
    std::getline(*line_stream,col,'\t');
    if(strcmp(col.c_str(),id_tag.c_str())!=0){
        delete line_stream;
        return false;
    }
    delete line_stream;
    return true;
}


bool check_id_full(sam_hdr_t *hdr,std::string tag1,std::string tag2){
    bool found_line = false;
    int line_pos = 0;
    std::string line;
    while(true){
        kstring_t str = KS_INITIALIZE;
        if (sam_hdr_find_line_pos(hdr, tag1.c_str(), line_pos, &str)!=0) {
            break;
        }
        else{
            // parse line to check if tb tag is present
            line = std::string(str.s);
            bool ret = check_id(line,tag2.c_str());
            if(ret){
                found_line=true;
            }
            line_pos++;
        }
        ks_free(&str);
    }
    return found_line;
}


void TInputFiles::addFile(const char* fn) {
    GStr sfn(fn);
    if (sfn!="-" && fileExists(fn)<2) {
        GError("Error: input file %s cannot be found!\n",
               fn);
    }

    freaders.Add(new TSamReader(fn));
}


bool TInputFiles::addSam(GSamReader* r, int fidx) {
    //requirement: all files must have the same number of SQ entries in the same order!
    kstring_t hd_line = KS_INITIALIZE;
    int res = sam_hdr_find_hd(r->header(), &hd_line);
    if (res < 0) GError("Error: failed to get @HD line from header!\n");
    //check for SO:coordinate
    kstring_t str = KS_INITIALIZE;
    if (sam_hdr_find_tag_hd(r->header(), "SO", &str)
        || !str.s
        || strcmp(str.s, "coordinate") )
        GError("Error: %s file not coordinate-sorted!\n", r->fileName());
    ks_free(&hd_line);
    bool tb_file=false; //was this file a product of TieBrush? (already merged)
    res=sam_hdr_find_tag_id(r->header(), "PG", "PN", "TieBrush", "VN", &str);
    if (res<-1) GError("Error: failed to lookup @PG line in header!\n");
    if (res==0) {
        tb_file=true;
#ifdef _DEBUG
        GMessage("DEBUG info: %s is a TieBrush (merged) file.\n", freaders[fidx]->fname.chars());
#endif
    }
    if (mHdr==NULL) { //first file
        headerfilename = r->fileName();
        headerfiletbMerged = tb_file;
        mHdr=sam_hdr_dup(r->header());
    }
    else { //check if this file has the same SQ entries in the same order
        //if it has more seqs, make it the main header
        int r_numrefs=sam_hdr_nref(r->header());
        int m_nrefs=sam_hdr_nref(mHdr);
        bool swapHdr=(r_numrefs>m_nrefs);
        sam_hdr_t *loHdr = swapHdr ? mHdr : r->header();
        sam_hdr_t *hiHdr = swapHdr ? r->header() : mHdr;
        int loNum = swapHdr ? m_nrefs : r_numrefs;
        //if (r_numrefs!=m_nrefs)
        // GError("Error: file %s has different number of reference sequences (%d)!\n", r->fileName(), r_numrefs);
        for (int i = 0; i < loNum; ++i) {
            str.l = 0;
            res = sam_hdr_find_tag_pos(loHdr, "SQ", i, "SN", &str);
            if (res < 0)
                GError("Error: failed to get @SQ SN #%d from header\n", i + 1);
            int m_tid = sam_hdr_name2tid(hiHdr, str.s);
            if (m_tid < -1)
                GError("Error: unexpected ref lookup failure (%s)!\n", str.s);
            if (m_tid < 0)
                GError("Error: ref %s not seen before!\n", str.s);
            int r_tid = sam_hdr_name2tid(loHdr, str.s);
            if (r_tid != m_tid)
                GError("Error: ref %s from file %s does not have the expected id#!", str.s, r->fileName());
        }
        if (swapHdr) {
            sam_hdr_destroy(mHdr);
            headerfilename = r->fileName();
            headerfiletbMerged = tb_file;
            mHdr=sam_hdr_dup(r->header());
        }
    }

    ks_free(&str);
    freaders[fidx]->samreader=r;
    freaders[fidx]->tbMerged=tb_file;

    if (fidx==freaders.Count()-1) { // last samreader entry
        // add any currently existing ID:SAMPLE lines to the maps
        load_hdr_samples(mHdr,this->headerfilename,this->headerfiletbMerged,true);

        for(int fi=0;fi<this->freaders.Count();fi++){ // add metadata about the files being collapsed and the index of each of them
            if(std::strcmp(this->freaders[fi]->fname.chars(),this->headerfilename.c_str())==0){ // this is the file which contributed header to the merged result - can safely
                continue;
            }
            else{
                load_hdr_samples(this->freaders[fi]->samreader->header(),this->freaders[fi]->fname.chars(),this->freaders[fi]->tbMerged,false);
            }
        }
        // now that we have a full list of samples - we can add them to the header
        for(auto& ls : this->lineno2sample){ // sorted order of the map by line number
            if(this->headerfiletbMerged && std::get<3>(ls.second)){ // if donor - can skip since already in the header
                continue;
            }
            int res_rg = sam_hdr_add_line(mHdr, "CO", ("SAMPLE:"+std::get<0>(ls.second)).c_str(),NULL);
                if(res_rg==-1){
                    std::cerr<<"unable to complete adding CO tags for file names"<<std::endl;
                    exit(-1);
                }
        }
        sam_hdr_add_pg(mHdr, "TieBrush",
                       "VN", pg_ver, "CL", pg_args.chars(), NULL);
        // sam_hdr_rebuild(mHdr); -- is this really needed?
    }
    return tb_file;
}


void TInputFiles::load_hdr_samples(sam_hdr_t* hdr,std::string filename,bool tbMerged,bool donor){
    int sample_line_pos = 0;
    if(tbMerged){
        bool found_line = false;
        int line_pos = 0;
        std::string line;
        while(true){
            kstring_t str = KS_INITIALIZE;
            if (sam_hdr_find_line_pos(hdr,"CO", line_pos, &str)!=0) {
                break;
            }
            else{
                line = std::string(str.s);
                bool ret = get_sample_from_line(line);

                if(ret){
                    found_line=true;
                    this->s2l_it = this->sample2lineno.insert(std::make_pair(line,std::make_tuple(this->max_sample_id,sample_line_pos,filename,donor)));
                    
                    // if(!this->s2l_it.second){ // not inserted
                    //     std::cerr<<"duplicate entries detected"<<std::endl;
                    //     exit(-1);
                    // }
                    
                    this->lineno2sample.insert(std::make_pair(this->max_sample_id,std::make_tuple(line,sample_line_pos,filename,donor)));
                    sample_line_pos++;
                    this->max_sample_id++;
                }
                line_pos++;
            }
            ks_free(&str);
        }
        if(!found_line){
            std::cerr<<"Collapsed file does not have any CO: lines in the header"<<std::endl;
            exit(-1);
        }
    }
    else{ // no sample was found - need to add current name to the header
        this->s2l_it = this->sample2lineno.insert(std::make_pair(get_full_path(filename),std::make_tuple(this->max_sample_id,sample_line_pos,"",donor)));
        if(!this->s2l_it.second){ // not inserted
            GMessage("this->s2l_it.second: %d\n", this->s2l_it.second);
            std::cerr<<"duplicate entries detected"<<std::endl;
            exit(-1);
        }
        this->lineno2sample.insert(std::make_pair(this->max_sample_id,std::make_tuple(get_full_path(filename),sample_line_pos,"",donor)));
        sample_line_pos++;
        this->max_sample_id++;
    }
}


std::string TInputFiles::get_full_path(std::string fname){
    const char *cur_path = fname.c_str();
    char *actualpath;


    actualpath = realpath(cur_path, NULL);
    if (actualpath != NULL){
        std::string ret = actualpath;
        free(actualpath);
        return ret;
    }
    else{
        std::cerr<<"could not resolve path: "<<fname<<std::endl;
        exit(-1);
    }
}


bool TInputFiles::get_sample_from_line(std::string& line){ // returns true if is sample pg line
    std::stringstream *line_stream = new std::stringstream(line);
    std::string col;

    // make sure it's CO - comment line
    std::getline(*line_stream,col,'\t');
    if(std::strcmp(col.c_str(),"@CO")!=0){
        delete line_stream;
        return false;
    }

    // check if ID == SAMPLE
    std::getline(*line_stream,col,':');
    if(std::strcmp(col.c_str(),"SAMPLE")!=0){
        delete line_stream;
        return false;
    }

    std::getline(*line_stream,col,'\t');
    line = col;
    delete line_stream;
    return true;
}


// adds a line to the header which tells whether the file has been processed with tiebrush before
bool TInputFiles::add_tb_tag_if_not_exists(sam_hdr_t *hdr){ // returns true if the tb_tag already existed and false if a new one was created
    // check if already exists
    bool found_tb_tag_line = false;
    int line_pos = 0;
    std::string line;
    while(true){
        kstring_t str = KS_INITIALIZE;
        if (sam_hdr_find_line_pos(hdr, "PG", line_pos, &str)!=0) {
            break;
        }
        else{
            // parse line to check if tb tag is present
            line = std::string(str.s);
            bool ret = check_id(line,"ID:TB_TAG");
            if(ret){
                if(found_tb_tag_line){
                    std::cerr<<"multiple tb tag lines found"<<std::endl;
                    exit(-1);
                }
                else{
                    found_tb_tag_line=true;
                }
            }
            line_pos++;
        }
        ks_free(&str);
    }
    if(!found_tb_tag_line){ // no tb_tag line - can now append the tag line
        int res_rg = sam_hdr_add_line(mHdr, "PG", "ID","TB_TAG",
                                      "VL","1",
                                      "XN","1",NULL);
        if(res_rg==-1){
            std::cerr<<"unable to complete adding pg tags for file names"<<std::endl;
            exit(-1);
        }
        return false;
    }
    return true;
}


// removed all lines with matching tags
void TInputFiles::delete_all_hdr_with_tag(sam_hdr_t *hdr,std::string tag1, std::string tag2){
    // check if already exists
    bool found_tb_tag_line = false;
    int line_pos = 0;
    std::string line;
    while(true){
        kstring_t str = KS_INITIALIZE;
        if (sam_hdr_find_line_pos(hdr, tag1.c_str(), line_pos, &str)!=0) {
            break;
        }
        else{
            // parse line to check if tb tag is present
            line = std::string(str.s);
            bool ret = check_id(line,tag2.c_str());
            if(ret){ // found line - remove
                int ret_rm = sam_hdr_remove_line_pos(hdr,"PG",line_pos);
                if(ret_rm!=0){
                    std::cerr<<"could not find requested header line"<<std::endl;
                    exit(-1);
                }
            }
            line_pos++;
        }
        ks_free(&str);
    }
}


// todo: merge header PG tags
int TInputFiles::start(){
    if (this->freaders.Count()==1) {
        //special case, if it's only one file it might be a list of file paths
        GStr& fname= this->freaders.First()->fname;
        //try to open it as a SAM/BAM/CRAM
        htsFile* hf=hts_open(fname.chars(), "r");
        if (hf==NULL || hf->format.category!=sequence_data) {
            //must be a list of file paths
            if (hf) hts_close(hf);
            FILE* flst=fopen(fname.chars(),"r");
            if (flst==NULL) GError("Error: could not open input file %s!\n",
                                   fname.chars());
            char* line=NULL;
            int lcap=5000;
            GMALLOC(line, lcap);
            freaders.Clear();
            while (fgetline(line,lcap,flst)) {
                GStr s(line);
                s.trim();
                if (s.length()<2 || s[0]=='#') continue; //skip comments/header in the list file, if any
                if (!fileExists(s.chars()))
                    GError("Error: cannot find alignment file %s !\n",s.chars());
                freaders.Add(new TSamReader(s.chars()));
            } //for each line in the list file
            GFREE(line);
            fclose(flst);
        }
        else { // single alignment file
            hts_close(hf);
        }
    }

    for (int i=0;i<freaders.Count();++i) {
        GSamReader* samrd=new GSamReader(freaders[i]->fname.chars(),
                                         SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
        bool tb_merged=addSam(samrd, i); //merge SAM headers etc.

        GSamRecord* brec=samrd->next();
        if (brec)
            recs.Add(new TInputRecord(brec, i, tb_merged));
    }
    return freaders.Count();
}


TInputRecord* TInputFiles::next() {
    //must free old current record first
    delete crec;
    crec=NULL;
    if (recs.Count()>0) {
        crec=recs.Pop();//lowest coordinate
        GSamRecord* rnext=freaders[crec->fidx]->samreader->next();
        if (rnext)
            recs.Add(new TInputRecord(rnext,crec->fidx, crec->tbMerged));
        return crec;
    }
    else return NULL;
}


void TInputFiles::stop() {
    for (int i=0;i<freaders.Count();++i) {
        freaders[i]->samreader->bclose();
    }
}