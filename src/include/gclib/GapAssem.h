#ifndef G_GAP_ASSEM_DEFINED
#define G_GAP_ASSEM_DEFINED
#include "GFastaFile.h"
#include "gdna.h"
#include "GList.hh"
#include "GHash.hh"
#include <ctype.h>

class GSeqAlign;
class MSAColumns;

extern const unsigned char GA_flag_IS_REF;
extern const unsigned char GA_flag_HAS_PARENT;
extern const unsigned char GA_flag_BAD_ALIGN;
extern const unsigned char GA_flag_PREPPED;

struct SeqDelOp {
  int pos;
  bool revcompl;
  SeqDelOp(int p, bool r) { pos=p; revcompl=r;}
  //
  bool operator==(SeqDelOp& d){
     return (this==&d);
     }
  bool operator<(SeqDelOp& d){
     return (this<&d);
     }

  };

class GASeq : public FastaSeq {
protected:
   int numgaps; //total number of accumulated gaps in this sequence
   short *ofs; //array of gaps at each position;
              //a negative value (-1) means DELETION of the nucleotide
              //at that position!

  #ifdef ALIGN_COVERAGE_DATA
   int* cov; //coverage of every nucleotide of this seq
             // it starts with 0 by itself
             // it'll be decreased by -1 for mismatching ends!
  #endif
  GList<SeqDelOp> delops; //delete operations
  void prepSeq(); //reverse complement if needed, and apply deletions (delops)
  //should only be called once when the MSA is complete (by GSeqAlign::finalize())
public:
  unsigned char flags; //8 general purpose boolean flags (bits)
  // bad_align flag is the last bit -- i.e. bit 7
  // all the others (0..6) are free for custom use
  GSeqAlign* msa;
  int msaidx; //actual index at which this sequence is to be found in GASeqAlign;
  int seqlen; // exactly the size of ofs[]
  int offset; //offset in the layout
  int ng_ofs; //non-gapped offset in the layout
              //(approx, for clipping constraints only)
  char revcompl; //0 = forward, 1=reverse
  int ext5; // layout-positional extension at 5' end
  int ext3; // layout-positional extension at 3' end
  //--
  int clp5; //ever "disproved" region at 5' end
  int clp3; //ever "disproved" region at 3' end
  //------- comparison operators (for GList) :
  //sorting by offset in cluster
  void allupper() {
     for (int i=0;i<len;i++) {
       seq[i]=toupper(seq[i]);
       }
     }
  void lowercase() {
     for (int i=0;i<len;i++) {
       seq[i]=tolower(seq[i]);
       }
     }

  void reverseComplement() {
    if (len==0) return;
    //ntCompTableInit();
    reverseChars(seq,len);
    for (int i=0;i<len;i++) seq[i]=ntComplement(seq[i]);
    }
  bool operator==(GASeq& d){
     return (offset==d.offset && strcmp(id,d.id)==0);
     }
  bool operator>(GASeq& d){
     return (offset>d.offset);
     }
  bool operator<(GASeq& d){
     return (offset<d.offset);
     }
  friend class GSeqAlign;
  //-------------------------------
  GASeq(FastaSeq& faseq, bool takeover=false);
  GASeq(GASeq& aseq); //copy constructor
  GASeq(const char* sname=NULL,int slen=0):FastaSeq(sname),
		  numgaps(0), ofs(NULL), delops(false, true, false),flags(0),msa(NULL),
		  msaidx(-1), seqlen(slen), offset(0),
		  ng_ofs(0),revcompl(0), ext5(0), ext3(0),
		  clp5(0), clp3(0) {
   if (seqlen>0) {
	GCALLOC(ofs, seqlen * sizeof(short));
	#ifdef ALIGN_COVERAGE_DATA
		GCALLOC(cov,seqlen*sizeof(int));
	#endif
   }
  };
  GASeq(const char* sname, const char* sdescr=NULL, const char* sseq=NULL, int slen=0, int soffset=0);
  GASeq(const char* sname, int soffset, int slen, int sclipL=0, int sclipR=0, char rev=0);
  ~GASeq();
  void refineClipping(GDynArray<char>& cons, int cpos, bool skipDels=false);
  void setGap(int pos, short gaplen=1); // set the gap in this pos
  void addGap(int pos, short gapadd); //extend the gap in this pos
  //bitno is 0 based here, for simplicity:
  inline void setFlag(unsigned char bitno) { flags |= ((unsigned char)1 << bitno); }
  inline void clearFlag(unsigned char bitno) { flags ^= ((unsigned char)1 << bitno); }
  inline bool hasFlag(unsigned char bitno) { return ( (((unsigned char)1 << bitno) & flags) !=0 ); }
  int getNumGaps() { return numgaps;  }
  int gap(int pos) { return ofs[pos];  }
  void removeBase(int pos); //remove the nucleotide at that position
  int endOffset() { return offset+seqlen+numgaps; }
  int endNgOffset() { return ng_ofs+seqlen; }
  int removeClipGaps(); //remove gaps within clipped regions
                        //offset should be corrected appropriately!
  void printGappedSeq(FILE* f, int baseoffs=0);
  void printGappedSeq(int baseoffs=0) { printGappedSeq(stdout, baseoffs); }
  void printGappedFasta(FILE* f);
  void printMFasta(FILE* f, int llen=60); //offset padded
  //void loadProcessing(); //to be called immediately after loading the sequence
                   // it will revCompl if needed and apply delops
  void finalize(); //delete inserts and reverse complement sequence if needed
  #ifdef ALIGN_COVERAGE_DATA
  void addCoverage(GASeq* s);
  #endif
  void reverseGaps(); //don't update offset and flags
                      //useful after reading mgblast gap info
  void revComplement(int alignlen=0);
  void toMSA(MSAColumns& msa, int nucValue=1);
};

// -- nucleotide origin -- for every nucleotide in a MSA column
//  this info is needed by SNP reporting
class NucOri {
public:
 GASeq* seq;
 int pos; //0-based position of nucleotide letter
 NucOri() { seq=NULL; pos=0; }
 NucOri(GASeq* s, int p) { seq=s; pos=p; }
 bool operator==(NucOri& d){
    return (strcmp(seq->id,d.seq->id)==0 && pos==d.pos);
    }
 bool operator>(NucOri& d){
    int cmp=strcmp(seq->id,d.seq->id);
    if (cmp==0) return pos>d.pos;
      else return cmp>0;
    }
 bool operator<(NucOri& d){
     int cmp=strcmp(seq->id,d.seq->id);
     if (cmp==0) return pos<d.pos;
       else return cmp<0;
    }
};

struct SeqClipOp {
  GASeq* seq;
  int clp[2];
  SeqClipOp(GASeq* s, int newclp5, int newclp3=-1) {
     seq=s;
     clp[0]=newclp5;
     clp[1]=newclp3;
     }
  //--
  bool operator==(SeqClipOp& d){
     return (this==&d);
     }
  bool operator>(SeqClipOp& d){
     return (this>&d);
     }
  bool operator<(SeqClipOp& d){
     return (this<&d);
     }
  };

class AlnClipOps :public GList<SeqClipOp> {
 public:
  char q_rev;
  int d5;
  int d3;
  //---
  int total;
  AlnClipOps():GList<SeqClipOp>(false,true,false) {
   total=0;
   d5=0;
   d3=0;
   q_rev=false;
   }
  bool add5(GASeq* s, int clp, float clipmax) {
   if (s->clp5<clp) {
     if (clipmax>0) {
        int maxovh = clipmax>1 ? (int)clipmax : iround(clipmax * (float)s->seqlen);
        if (clp>maxovh) return false;
        }
     //----- base test: the read should be left with no less than 25% of its length
     if (s->seqlen-s->clp3-clp < (s->seqlen >> 2)) return false;
     total+=10000+clp-s->clp5;
     Add(new SeqClipOp(s,clp,-1));
     }
   return true;
   }
  bool add3(GASeq* s, int clp, float clipmax) {
   if (s->clp3<clp) {
     if (clipmax>0) {
        int maxovh = clipmax>1 ? (int)clipmax : iround(clipmax * (float)s->seqlen);
        if (clp>maxovh) return false;
        }
     //----- base test: if the read is left with less than 25% of its length
     if (s->seqlen-s->clp5-clp < (s->seqlen >> 2)) return false;
     total+= 10000+clp-s->clp3;
     Add(new SeqClipOp(s,-1,clp));
     }
    return true;
   }
  bool add(GASeq* s, int clp5, int clp3, float clipmax) {
   int newclp5=-1;
   int newclp3=-1;
   int add=0;
   if (s->clp5<clp5) {
     if (clipmax>0) {
        int maxovh = clipmax>1 ? (int)clipmax : iround(clipmax * (float)s->seqlen);
        if (clp5>maxovh) return false;
        }
      //----- base test: if the read is left with less than 25% of its length!
      if (s->seqlen-s->clp3-clp5 < (s->seqlen >> 2)) return false;
      add+= 10000+clp5-s->clp5;
      newclp5=clp5;
      }
     else clp5=s->clp5;
   if (s->clp3<clp3) {
     if (clipmax>0) {
        int maxovh = clipmax>1 ? (int)clipmax : iround(clipmax * (float)s->seqlen);
        if (clp3>maxovh) return false;
        }
      //----- base test: if the read is left with less than 25% of its length!
      if (s->seqlen-clp5-clp3 < (s->seqlen >> 2)) return false;
      add+=10000+clp3-s->clp3;
      newclp3=clp3;
      }
   if (add>0) {
     total+=add;
     Add(new SeqClipOp(s,newclp5,newclp3));
     }
   return true;
   }
};

class GAlnColumn {
 protected:
   struct NucCount {
     char nt; // A, C, G, T, N or -
             // precisely in this order (except after qsort)
     int count;
     void set(char l, int num=0) { nt=l;count=num; }
     };
  enum { ncA=0, ncC, ncG, ncT, ncN, ncGap };
  NucCount counts[6];
  bool countsSorted;
 public:
  bool hasClip;
  char consensus;
  int layers; //total "thickness"
  NucOri* clipnuc;
  GList<NucOri>* nucs;
  friend int qsortnuc(const void* p1, const void* p2);
  //int total() { return numgaps+numN+numA()+numC()+numG()+numT(); }
  GAlnColumn():countsSorted(false),hasClip(false), consensus(0),
		  layers(0), clipnuc(NULL), nucs(NULL) {        //sorted?, free?, unique?
   nucs=new GList<NucOri>(false,true,false);
   /*lstC=new GList<NucOri>(false,true,false);
   lstG=new GList<NucOri>(false,true,false);
   lstT=new GList<NucOri>(false,true,false);*/
   counts[ncA].set('A');
   counts[ncC].set('C');
   counts[ncG].set('G');
   counts[ncT].set('T');
   counts[ncN].set('N');
   counts[ncGap].set('-');
  }
  ~GAlnColumn() {
   delete nucs;
   if (clipnuc!=NULL) delete clipnuc;
  }
  void addGap(int nucVal=1) { counts[ncGap].count+=nucVal;
                  layers++; //-- Not a "layer", actually
                  //numgaps++;
                  }
  void addNuc(GASeq* seq, int pos, bool clipped=false, short nucVal=1) {
   //assumes the seq is already loaded and reverse complemented if necessary
   //position is precisely where it should be
   if (clipped) {
      if (hasClip==false) {
            hasClip=true;
            clipnuc = new NucOri(seq,pos);
            }
       return;
       }
   char c=(char)toupper(seq->seq[pos]);
   switch (c) {
       case 'A':nucs->Add(new NucOri(seq,pos));
                counts[ncA].count+=nucVal;
                layers++;
                break;
       case 'C':nucs->Add(new NucOri(seq,pos));
                counts[ncC].count+=nucVal;
                layers++;
                break;
       case 'G':nucs->Add(new NucOri(seq,pos));
                counts[ncG].count+=nucVal;
                layers++;
                break;
       case 'T':nucs->Add(new NucOri(seq,pos));
                counts[ncT].count+=nucVal;
                layers++;
                break;
       case '-': //this shouldn't be the case!
       case '*':counts[ncGap].count+=nucVal;
                layers++;
                //numgaps++;
                break;
       default: nucs->Add(new NucOri(seq,pos));
                counts[ncN].count+=nucVal;
                layers++;
                //numN++;
       }//switch
   }

  char bestChar(int16_t *qscore=NULL);
  void remove(); //removes a nucleotide from all involved sequences
                 //adjust all affected offsets in the alignment
};

// A MSA columns container
class MSAColumns {
   int size;
 public:
   static bool removeConsGaps;
   static bool refineClipping;
   GAlnColumn* columns;
   int baseoffset;
   int mincol;
   int maxcol;
   MSAColumns(int len, int baseofs=0) {
     columns=new GAlnColumn[len];
     size=len;
     baseoffset=baseofs;
     mincol=INT_MAX;
     maxcol=0;
     }
   ~MSAColumns() {
     size=0;
     baseoffset=0;
     delete[] columns;
     }
   GAlnColumn& operator[](int idx) {
    if (idx<0 || idx>=size)
         GError("MSAColumns op[]: bad index %d (size=%d)\n", idx,size);
    return columns[idx];
    }
   int len() { return maxcol-mincol+1; }
   void updateMinMax(int minc, int maxc) {
    if (minc<mincol) mincol=minc;
    if (maxc>maxcol) maxcol=maxc;
    }
};


//-----------------------------------------------
// a sequence alignment: could be pairwise or MSA
class GSeqAlign :public GList<GASeq> {
   //static unsigned int counter;
   int length;
   int minoffset;
   //int consensus_cap;
   void buildMSA(bool refWeighDown=false);
   void ErrZeroCov(int col);
 public:
    bool refinedMSA; //if refineMSA() was applied
    MSAColumns* msacolumns;
    unsigned int ordnum; //order number -- when it was created
              // the lower the better (earlier=higher score)
   int ng_len;     //ungapped length and minoffset (approximative,
   int ng_minofs;  //  for clipping constraints only)
   int badseqs;
   //char* consensus; //consensus sequence (built by refineMSA())
   GDynArray<char> consensus;
   GDynArray<int16_t> consensus_bq;
   //int consensus_len;
   friend class GASeq;
   bool operator==(GSeqAlign& d){
     return (this==&d);
     }
  bool operator>(GSeqAlign& d){
     return (this>&d);
     }
  bool operator<(GSeqAlign& d){
     return (this<&d);
     }
  //--
  //GSeqAlign():GList<GASeq>(true,true,false), length(0), minoffset(0),
  GSeqAlign():GList<GASeq>(false,true,false), length(0), minoffset(0),
  		refinedMSA(false), msacolumns(NULL), ordnum(0),
  		ng_len(0),ng_minofs(0), badseqs(0), consensus(512), consensus_bq(512) {
    //default is: sorted by GASeq offset, free nodes, non-unique
    }
  GSeqAlign(bool sorted, bool free_elements=true, bool beUnique=false)
     :GList<GASeq>(sorted,free_elements,beUnique), length(0), minoffset(0),
  		refinedMSA(false), msacolumns(NULL), ordnum(0),
  		ng_len(0),ng_minofs(0), badseqs(0), consensus(512), consensus_bq(512) {
    }
  //void incOrd() { ordnum = ++counter; }
  void incOrd() { ordnum++; }
  //first time creation from a pairwise alignment:
  #ifdef ALIGN_COVERAGE_DATA
  GSeqAlign(GASeq* s1, int l1, int r1, GASeq* s2, int l2, int r2);
  #else
  GSeqAlign(GASeq* s1, GASeq* s2);
  #endif
  ~GSeqAlign() {
    if (msacolumns!=NULL) delete msacolumns;
    }
  int len() { return length; }

  void revComplement();
  void addSeq(GASeq* s, int soffs, int ngofs);
  void injectGap(GASeq* seq, int pos, int xgap);
  void removeBase(GASeq* seq, int pos);
  void extendConsensus(char c, int16_t bq=SHRT_MIN);
  //try to propagate the planned trimming of a read
  //to the whole MSA containing it
  // returns false if too much is trimmed of any component read
  void applyClipping(AlnClipOps& clipops);

  bool evalClipping(GASeq* seq, int c5, int c3, float clipmax,
                             AlnClipOps& clipops);

  //merge other alignment into this msa
  //seq->id MUST be the same with oseq->id
  // *if OK, gaps AND coverage values are propagated
  //  and omsa is deallocated
  // *if not OK <=> the layout doesn't accept the merge
  //  due to clipmax constraint, then nothing happens
  bool addAlign(GASeq* seq, GSeqAlign* omsa, GASeq* oseq);
  void finalize(); //delete inserts and reverse complement sequences as needed
  void print(FILE* f, char c=0); //debug printing one-line alignments
  void print() { print(stdout); }
  void removeColumn(int column);
  void freeMSA();
  void refineMSA(bool refWeighDown=false, bool redo_ends=false);
      // find consensus, refine clipping, remove gap-columns
  void writeACE(FILE* f, const char* name, bool refWeighDown=false);
  void writeMSA(FILE* f, int linelen=60); //write as multi-FASTA (MAF?) file
  void writeInfo(FILE* f, const char* name, bool refWeighDown=false);
};

int compareOrdnum(void* p1, void* p2);
int compareCounts(void* p1, void* p2);

#endif
