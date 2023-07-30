/*  GSam.h

    Copyright (C) 2019 Geo Pertea
    
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

#ifndef _G_SAM_H
#define _G_SAM_H

#include <assert.h>
#include <string>
#include <iostream>
#include <gclib/GBase.h>
#include <gclib/GList.hh>
#include <htslib/htslib/kstring.h>
#include <htslib/htslib/sam.h>
#include <htslib/htslib/cram.h>

class GSamReader;
class GSamWriter;

enum GSamFileType {
   GSamFile_SAM=1,
   GSamFile_UBAM,
   GSamFile_BAM,
   GSamFile_CRAM
};

class GSamRecord: public GSeg {
   friend class GSamReader;
   friend class GSamWriter;
   bam1_t* b=NULL;
   // b->data has the following strings concatenated:
   //  qname (including the terminal \0)
   //  +cigar (each event encoded on 32 bits)
   //   +seq  (4bit-encoded)
   //    +qual
   //     +aux

   union {
	  uint16_t iflags=0;
	struct {
    	  bool novel         :1; //if set, the destructor must free b
    	  bool hard_Clipped  :1;
    	  bool soft_Clipped  :1;
    	  bool has_Introns   :1;
      };
   };
   sam_hdr_t* b_hdr=NULL;
 public:
   GVec<GSeg> exons; //coordinates will be 1-based
   int clipL=0; //soft clipping data, as seen in the CIGAR string
   int clipR=0;
   int mapped_len=0; //sum of exon lengths
   //FIXME: DEBUG only fields
#ifdef _DEBUG
   char* _cigar=NULL;
   const char* _read=NULL;
#endif
   // -- DEBUG only fields
   bool isHardClipped() { return hard_Clipped; }
   bool isSoftClipped() { return soft_Clipped; }
   bool hasIntrons() { return has_Introns; }
   //created from a reader:
   void bfree_on_delete(bool b_free=true) { novel=b_free; }
   GSamRecord() { }
   GSamRecord(bam1_t* from_b, sam_hdr_t* b_header=NULL, bool b_free=true):b_hdr(b_header) {
      if (from_b==NULL) {
           b=bam_init1();
           novel=true;
      }
      else {
           b=from_b; //it'll take over from_b
           novel=b_free;
#ifdef _DEBUG
           _cigar=cigar();
           _read=name();
#endif
           setupCoordinates();//set 1-based coordinates (start, end and exons array)
      }
   }

   void init(bam1_t* from_b, sam_hdr_t* b_header=NULL, bool adopt_b=false) {
	   clear();
	   novel=adopt_b;
	   b=from_b;
#ifdef _DEBUG
       _cigar=cigar();
       _read=name();
#endif
	   b_hdr=b_header;
	   setupCoordinates();
   }

   //deep copy constructor:
   GSamRecord(GSamRecord& r):GSeg(r.start, r.end), iflags(r.iflags), b_hdr(r.b_hdr),
		   exons(r.exons), clipL(r.clipL), clipR(r.clipR), mapped_len(r.mapped_len)
		   {
            // GMessage("@@ Inside 'GSamRecord' deep copy constructor!\n");
	      //makes a new copy of the bam1_t record etc.
	      b=bam_dup1(r.b);
	      novel=true; //will also free b when destroyed
#ifdef _DEBUG
	      _cigar=Gstrdup(r._cigar);
	      _read=r._read;
#endif
   }

   const GSamRecord& operator=(GSamRecord& r) {
      //copy operator
      //makes a new copy of the bam1_t struct etc.
      // GMessage("@@ Inside 'GSamRecord' copy constructor!\n");

      clear();
      // GMessage("\tr.refName()  : %s\n", r.refName());
      b=bam_dup1(r.b);
      b_hdr = r.b_hdr;
      // GMessage("\tnew.refName(): %s\n", this->refName());

      iflags=r.iflags;
      novel=true; //will also free b when destroyed
      start=r.start;
      end=r.end;
      exons = r.exons;
      clipL = r.clipL;
      clipR = r.clipR;
      mapped_len=r.mapped_len;
#ifdef _DEBUG
      _cigar=Gstrdup(r._cigar);
      _read=r._read;
#endif
      return *this;
   }
   

   bool operator == (const GSamRecord& brec) {
      // if (brec.start == start && brec.end == end && brec.b->data == b->data) {
      if((bam_get_qual(brec.b)==bam_get_qual(b)) && (bam_get_cigar(brec.b)==bam_get_cigar(b)) && (bam_get_seq(b)==bam_get_seq(b)) && (bam_get_aux(b)==bam_get_aux(b)) ) {
         return true;
      } else {
         return false;
      }
   }

   bool operator<(const GSamRecord& brec) { // sort by strand last
      if (start==brec.start){
         if(end==brec.end) {
				if(end==brec.end) {
					return false;
				}
				else{
					return (end < brec.end);
				}
			}
      }
      return (start < brec.start);
   }

     void setupCoordinates();

     void clear() {
        if (novel) {
           bam_destroy1(b);
           //novel=false;
        }
        b=NULL;
        exons.Clear();
        mapped_len=0;
        b_hdr=NULL;
        iflags=0;
#ifdef _DEBUG
        GFREE(_cigar);
        _read=NULL;
#endif
    }

    ~GSamRecord() {
       clear();
    }

    void print_cigar(bam1_t *al){
        for (uint8_t c=0;c<al->core.n_cigar;++c){
            uint32_t *cigar_full=bam_get_cigar(al);
            int opcode=bam_cigar_op(cigar_full[c]);
            int length=bam_cigar_oplen(cigar_full[c]);
            std::cout<<length<<bam_cigar_opchr(opcode);
        }
        std::cout<<std::endl;
    }

    void print_seq(bam1_t *new_rec){
        int32_t qlen = new_rec->core.l_qseq;
        int8_t *buf = NULL;
        buf = static_cast<int8_t *>(realloc(buf, qlen+1));
        buf[qlen] = '\0';
        uint8_t* seq = bam_get_seq(new_rec);
        for (int i = 0; i < qlen; ++i)
            buf[i] = bam_seqi(seq, i);
        for (int i = 0; i < qlen; ++i) {
            buf[i] = seq_nt16_str[buf[i]];
        }
        std::string str_seq((char*)(char*)buf);
        std::cout<<str_seq<<std::endl;
    }

    // taken from samtools/bam_import.c
    static inline uint8_t * alloc_data(bam1_t *b, size_t size)
    {
        if (b->m_data < size)
        {
            b->m_data = size;
            kroundup32(b->m_data);
            b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
        return b->data;
    }

    bam1_t * bam_update(bam1_t * b,
                        const size_t nbytes_old,
                        const size_t nbytes_new,
                        uint8_t * field_start){ // from pysam
        int d = nbytes_new - nbytes_old;
        int new_size;
        size_t nbytes_before;
        uint8_t * retval = NULL;

        // no change
        if (d == 0)
            return b;

        // new size of total data
        new_size = d + b->l_data;

        // fields before field in data
        nbytes_before = field_start - b->data;

        if (b->l_data != 0)
        {
            assert(nbytes_before >= 0);
            assert(nbytes_before <= b->l_data);
        }

        // increase memory if required
        if (d > 0)
        {
            retval = alloc_data(b, new_size);
            if (retval == NULL)
                return NULL;
            field_start = b->data + nbytes_before;
        }

        // move data after field to new location
        memmove(field_start + nbytes_new,
                field_start + nbytes_old,
                b->l_data - (nbytes_before + nbytes_old));

        // adjust l_data
        b->l_data = new_size;

        return b;
    }

    void replace_qname(int id){ // replace the name with an ID
        char * p = bam_get_qname(b);

        std::string qname = std::to_string(id);
        int l = qname.size()+1;
        int l_extranul = 0;
        if (l % 4 != 0){
            l_extranul = 4 - l % 4;
        }

        bam1_t * retval = bam_update(b,b->core.l_qname,l + l_extranul,(uint8_t*)p);
        if (retval == NULL){
            std::cerr<<"could not allocate memory"<<std::endl;
            exit(-1);
        }

        b->core.l_extranul = l_extranul;
        b->core.l_qname = l + l_extranul;

        p = bam_get_qname(b);

        strcpy(p,qname.c_str());
        uint16_t x = 0;

        for (int x=l;x<l+l_extranul;x++){
            p[x] = '\0';
        }
    }

    void parse_error(const char* s) {
      GError("SAM parsing error: %s\n", s);
    }

    bam1_t* get_b() { return b; }

    void set_mdata(int32_t mtid, int32_t m0pos, //0-based coordinate, -1 if not available
                     int32_t isize=0) { //mate info for current record
      b->core.mtid=mtid;
      b->core.mpos=m0pos; // should be -1 if '*'
      b->core.isize=isize; //should be 0 if not available
    }

    void set_flags(uint16_t samflags) {
      b->core.flag=samflags;
    }

    void add_aux(const char* str); //adds one aux field in plain SAM text format (e.g. "NM:i:1")
    int  add_aux(const char tag[2], char atype, int len, uint8_t *data) {
      //IMPORTANT:  strings (Z,H) should include the terminal \0
     return bam_aux_append(b, tag, atype, len, data);
    }

    int add_tag(const char tag[2], char atype, int len, uint8_t *data) {
      //same with add_aux()
      //IMPORTANT:  strings type (Z,H) should include the terminal \0
      return bam_aux_append(b, tag, atype, len, data);
    }

    int add_int_tag(const char tag[2], int64_t val) { //add or update int tag
    	return bam_aux_update_int(b, tag, val);
    }
    int remove_tag(const char tag[2]);
    inline int delete_tag(const char tag[2]) { return remove_tag(tag); }

 //--query methods:
 uint32_t flags() { return b->core.flag; } //return SAM flags
 bool isUnmapped() { return ((b->core.flag & BAM_FUNMAP) != 0); }
 bool isMapped() { return ((b->core.flag & BAM_FUNMAP) == 0); }

 bool isMateUnmapped() { return ((b->core.flag & BAM_FMUNMAP) != 0); }
 bool isMateMapped() { return ((b->core.flag & BAM_FMUNMAP) == 0); }
 bool isPaired() { return ((b->core.flag & BAM_FPAIRED) != 0); }
 bool isProperPaired() { return ((b->core.flag & BAM_FPROPER_PAIR) != 0); }

 const char* name() { return bam_get_qname(b); }
 int pairOrder() {
    //which read in the pair: 0 = unpaired, 1=first read, 2=second read
    int r=0;
    if ((b->core.flag & BAM_FREAD1) != 0) r=1;
    else if ((b->core.flag & BAM_FREAD2) != 0) r=2;
    return r;
    }
 bool revStrand() {
   //this is the raw alignment strand, NOT the transcription/splice strand
   return ((b->core.flag & BAM_FREVERSE) != 0);
 }
 const char* refName() {
   return (b_hdr!=NULL) ?
         ((b->core.tid<0) ? "*" : b_hdr->target_name[b->core.tid]) : NULL;
   }
 inline int32_t refId() { return b->core.tid; }
 inline int32_t mate_refId() { return b->core.mtid; }
 const char* mate_refName() {
    return (b_hdr!=NULL) ?
       ((b->core.mtid<0) ? "*" : b_hdr->target_name[b->core.mtid]) : NULL;
    }
 inline int32_t insertSize() { return b->core.isize; }
 inline int32_t mate_start() { return b->core.mpos<0? 0 : b->core.mpos+1; }
 inline uint8_t mapq() { return b->core.qual; }
 //int find_tag(const char tag[2], uint8_t* & s, char& tag_type);
 uint8_t* find_tag(const char tag[2]);

 char* tag_str(const char tag[2]); //return tag value for tag type 'Z'
 int64_t tag_int(const char tag[2], int nfval=0); //return numeric value of tag (for numeric types)
 double tag_float(const char tag[2]); //return float value of tag (for float types)
 char tag_char(const char tag[2]); //return char value of tag (for type 'A')
 char tag_char1(const char tag[2]);
 char spliceStrand(); // '+', '-' from the XS tag, or '.' if no XS tag
 char* sequence(); //user should free after use
 char* qualities();//user should free after use
 char* cigar(); //returns text version of the CIGAR string; user must deallocate

 void unpair_mate_refName() {
   b->core.mtid = -1;
 }
 void unpair_mate_start() {
   b->core.mpos = -1;
 }

};

// from sam.c:
#define FTYPE_BAM  1
#define FTYPE_READ 2

class GSamReader {
   htsFile* hts_file;
   char* fname;
   sam_hdr_t* hdr;
   bam1_t* b_next; //for light next(GBamRecord& b)
 public:
   void bopen(const char* filename, int32_t required_fields,
		   const char* cram_refseq=NULL) {
	      hts_file=hts_open(filename, "r");
	      if (hts_file==NULL)
	         GError("Error: could not open alignment file %s \n",filename);
	      if (hts_file->is_cram && cram_refseq!=NULL) {
	              hts_set_opt(hts_file, CRAM_OPT_REFERENCE, cram_refseq);
    	  }

          hts_set_opt(hts_file, CRAM_OPT_REQUIRED_FIELDS,
	    			  required_fields);
	      fname=Gstrdup(filename);
	      hdr=sam_hdr_read(hts_file);
   }

   void bopen(const char* filename, const char* cram_refseq=NULL) {
      hts_file=hts_open(filename, "r");
      if (hts_file==NULL)
         GError("Error: could not open alignment file %s \n",filename);
      if (hts_file->is_cram) {
    	  if (cram_refseq!=NULL) {
              hts_set_opt(hts_file, CRAM_OPT_REFERENCE, cram_refseq);
    	  }
    	  else hts_set_opt(hts_file, CRAM_OPT_REQUIRED_FIELDS,
    		SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_MAPQ|SAM_CIGAR|
			SAM_RNEXT|SAM_PNEXT|SAM_TLEN|SAM_AUX );

      }
      fname=Gstrdup(filename);
      hdr=sam_hdr_read(hts_file);
   }

   GSamReader(const char* fn, int32_t required_fields,
		   const char* cram_ref=NULL):hts_file(NULL),fname(NULL), hdr(NULL), b_next(NULL) {
      bopen(fn, required_fields, cram_ref);
   }

   GSamReader(const char* fn, const char* cram_ref=NULL):hts_file(NULL),fname(NULL),
		   hdr(NULL), b_next(NULL) {
      bopen(fn, cram_ref);
   }

   sam_hdr_t* header() {
      return hts_file ? hdr : NULL;
   }
   const char* fileName() {
      return fname;
   }

   const char* refName(int tid) {
	   if (!hts_file) return NULL;
	   return hdr->target_name[tid];
   }

   void bclose() {
      if (hts_file) {
   	    if (hdr!=NULL) sam_hdr_destroy(hdr);
   	    hdr=NULL;
        hts_close(hts_file);
        hts_file=NULL;
        }
    }

   ~GSamReader() {
      if (b_next) bam_destroy1(b_next);
      bclose();
      GFREE(fname);
   }

   void rewind() {
     if (fname==NULL) {
       GMessage("Warning: GSamReader::rewind() called without a file name.\n");
       return;
     }
     bclose();
     char* ifname=fname;
     bopen(ifname);
     GFREE(ifname);
  }

   //the caller has to FREE the created GSamRecord
   GSamRecord* next() {
      if (hts_file==NULL)
        GError("Warning: GSamReader::next() called with no open file.\n");
      bam1_t* b = bam_init1();
      if (sam_read1(hts_file, hdr, b) >= 0) {
        GSamRecord* bamrec=new GSamRecord(b, hdr, true);
        return bamrec;
      }
      bam_destroy1(b);
      return NULL;
   }

   bool next(GSamRecord& rec) {
       if (hts_file==NULL)
	        GError("Warning: GSamReader::next() called with no open file.\n");
	   if (b_next==NULL) b_next=bam_init1();
       if (sam_read1(hts_file, hdr, b_next) >= 0) {
	        rec.init(b_next, hdr, false);
	        return true;
	   }
       return false;
   }
};

//basic BAM/SAM/CRAM writer class
class GSamWriter {
   htsFile* bam_file;
   sam_hdr_t* hdr;
 public:
   void create(const char* fname, sam_hdr_t* bh, GSamFileType ftype=GSamFile_BAM) {
     hdr=sam_hdr_dup(bh);
     create(fname, ftype);
   }

   GSamWriter(const char* fname, sam_hdr_t* bh, GSamFileType ftype=GSamFile_BAM):
	                                    bam_file(NULL),hdr(NULL) {
      create(fname, bh, ftype);
   }

   void create(const char* fname, GSamFileType ftype=GSamFile_BAM) {
      if (hdr==NULL)
         GError("Error: no header data provided for GSamWriter::create()!\n");
	  kstring_t mode=KS_INITIALIZE;
      kputc('w', &mode);
      switch (ftype) {
         case GSamFile_BAM:
        	kputc('b', &mode);
        	break;
         case GSamFile_UBAM:
        	kputs("bu", &mode);
        	break;
         case GSamFile_CRAM:
        	kputc('c', &mode);
        	break;
         case GSamFile_SAM:
         	break;
         default:
      	   GError("Error: unrecognized output file type!\n");
      }
      bam_file = hts_open(fname, mode.s);
      if (bam_file==NULL)
         GError("Error: could not create output file %s\n", fname);
      if (sam_hdr_write(bam_file, hdr)<0)
    	  GError("Error writing header data to file %s\n", fname);
      ks_free(&mode);
   }

   sam_hdr_t* header() { return hdr; }
   GSamWriter(const char* fname, const char* hdr_file, GSamFileType ftype=GSamFile_BAM):
	                                             bam_file(NULL),hdr(NULL) {
	  //create an output file fname with the SAM header copied from hdr_file
      htsFile* samf=hts_open(hdr_file, "r");
      if (samf==NULL)
    	  GError("Error: could not open SAM file %s\n", hdr_file);
      hdr=sam_hdr_read(samf);
      if (hdr==NULL)
    	  GError("Error: could not read header data from %s\n", hdr_file);
      hts_close(samf);
      create(fname, ftype);
   }

   ~GSamWriter() {
      hts_close(bam_file);
      sam_hdr_destroy(hdr);
   }

   sam_hdr_t* get_header() { return hdr; }

   int32_t get_tid(const char *seq_name) {
      if (hdr==NULL)
         GError("Error: missing SAM header (get_tid())\n");
      return sam_hdr_name2tid(hdr, seq_name);
      }

   void write(GSamRecord* brec) {
      if (brec!=NULL) {
          if (sam_write1(this->bam_file,this->hdr, brec->b)<0)
        	  GError("Error writing SAM record!\n");
      }
   }

   void write(bam1_t* xb) {
     if (sam_write1(this->bam_file, this->hdr, xb)<0)
    	 GError("Error writing SAM record!\n");
   }
};
#endif
