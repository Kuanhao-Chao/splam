//
// Created by sparrow on 11/9/20.
//

#include <vector>
#include <math.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <stdlib.h>
#include <fstream>

#include <gclib/GStr.h>
#include "GSam.h"

#ifndef TIEBRUSH_COMMONS_H
#define TIEBRUSH_COMMONS_H

bool parse_pg_sample_line(std::string& line){ // returns true if is sample pg line
    std::stringstream *line_stream = new std::stringstream(line);
    std::string col;

    // make sure it's CO
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

    std::getline(*line_stream,col,'\n');
    line = col;
    delete line_stream;
    return true;
}

void load_sample_info(sam_hdr_t* hdr,std::vector<std::string>& info){
    bool found_sample_line = false;
    int line_pos = 0;
    std::string line;
    while(true){
        kstring_t str = KS_INITIALIZE;
        if (sam_hdr_find_line_pos(hdr, "CO", line_pos, &str)!=0) {
            if(!found_sample_line){
                GError("Error: no sample lines found in header");
            }
            break;
        }
        else{
            // parse line to check if indeed SAMPLE CO
            line = std::string(str.s);
            bool ret = parse_pg_sample_line(line);
            if(ret){
                found_sample_line=true;
                info.push_back(line);
            }
            line_pos++;
        }
        ks_free(&str);
    }
}

std::string get_full_path(std::string fname){
    const char *cur_path = fname.c_str();
    char *actualpath;


    actualpath = realpath(cur_path, NULL);
    if (actualpath != NULL){
        return std::string(actualpath);
        free(actualpath);
    }
    else{
        std::cerr<<"could not resolve path: "<<fname<<std::endl;
        exit(-1);
    }
}

struct Index_Builder{
public:
    Index_Builder() = default;
    ~Index_Builder() = default;

    void add(uint64_t dupcount,std::fstream* out_fp){
        buffer[nr*4]   = (dupcount >> 24) & 0xFF;
        buffer[(nr*4)+1] = (dupcount >> 16) & 0xFF;
        buffer[(nr*4)+2] = (dupcount >> 8) & 0xFF;
        buffer[(nr*4)+3] = dupcount & 0xFF;

        nr += 1;

        if(nr*4 == 1024*4096){
            write(out_fp);
        }
    }
    void clear(std::fstream* out_fp){
        write(out_fp);
    }
private:
    char buffer[1024*4096];
    int nr = 0;

    void write(std::fstream* out_fp){
        out_fp->write(buffer,nr*4);
        nr=0;
    }
};

class Index_Loader{
public:
    Index_Loader() = default;
    Index_Loader(GStr& ifname):index_fname(ifname.chars()){
        this->index_ss.open(ifname.chars(),std::ios::in | std::ios::binary);
        this->index_ss.unsetf(std::ios::skipws);
    };
    ~Index_Loader(){
        if(index_ss.is_open()){
            index_ss.close();
        }
        for(int i=0;i<this->tbd_streams.size();i++){
            delete this->tbd_streams[i];
        }
    };

    void load(GStr& ifname){
        index_fname = ifname.chars();
        this->index_ss.open(ifname.chars(),std::ios::in | std::ios::binary);
        // read header of the index to determine byte-offsets of each sample line
        std::string header;
        std::getline(this->index_ss,header);
        std::stringstream hss(header);
        std::string num_bytes;
        this->header_end_pos = this->index_ss.tellg();
        while(std::getline(hss,num_bytes,'|')){
            this->byte_offsets.push_back(std::atoi(num_bytes.c_str()));
        }

        this->index_ss.unsetf(std::ios::skipws);
    }

    void init(std::vector<int>& lst){ // initializes file streams to the begining of each sample
        // initialize a vector of ifstreams and seek to the starting byte of each sample
        for(int i=0;i<lst.size();i++){ // for each requested sample in the list
            std::fstream* nss;
            this->tbd_streams.push_back(nss);

            this->tbd_streams.back() = new std::fstream();
            this->tbd_streams.back()->open(this->index_fname,std::ios::in | std::ios::binary);
            this->tbd_streams.back()->unsetf(std::ios::skipws);

            // skip to the start of the current sample
            this->tbd_streams.back()->seekg(this->header_end_pos+this->byte_offsets[lst[i]],std::ios::beg);
        }
    }

    void next(uint32_t &dup_val, std::vector<int>& samples) { // lst is the list of samples for which to extract the values
        for(int i=0;i<this->tbd_streams.size();i++){
            if (!this->tbd_streams[i]->is_open() || !this->tbd_streams[i]->good())
                GError("Warning: Index::next() called with no open file.\n");
            char buffer[4];
            this->tbd_streams[i]->read(buffer,4);
            if(!(*this->tbd_streams[i])) {
                GError("Error: only %d bytes could be loaded\n",this->tbd_streams[i]->gcount());
            }
            int tmp = ((uint8_t)buffer[0] << 24) | ((uint8_t)buffer[1] << 16) | ((uint8_t)buffer[2] << 8) | (uint8_t)buffer[3];
            dup_val+=tmp;
            if(tmp>0){
                samples.push_back(i);
            }
        }
    }

private:
    std::string index_fname = "";
    std::fstream index_ss;

    int header_end_pos = 0;

    std::vector<int> byte_offsets; // offsets of samples
    std::vector<std::fstream*> tbd_streams; // streams for lines within the index file corresponding to desired samples
};

#endif //TIEBRUSH_COMMONS_H
