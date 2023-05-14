#ifndef _GALIGNEXTEND_H

//greedy gapped alignment extension
//(mostly lifted from NCBI's blast gapped extension code)

#include "GBase.h"
#include "GList.hh"
#include "gdna.h"

enum {
    gxEDIT_OP_MASK = 0x3,
    gxEDIT_OP_ERR  = 0x0,
    gxEDIT_OP_INS  = 0x1,
    gxEDIT_OP_DEL  = 0x2,
    gxEDIT_OP_REP  = 0x3
};

#define GX_EDITOP_VAL(op) ((op) >> 2)
#define GX_EDITOP_GET(op) ((op) & gxEDIT_OP_MASK)
#define GX_EDITOP_CONS(op, val) (((val) << 2) | ((op) & gxEDIT_OP_MASK))

enum {c_black=0,
      c_red, c_green,c_brown,c_blue,c_magenta,c_cyan,c_white
      };

void color_fg(int c, FILE* f=stderr);
void color_bg(int c, FILE* f=stderr);
void color_resetfg(FILE* f=stderr);
void color_resetbg(FILE* f=stderr);
void color_reset(FILE* f=stderr);
void color_normal(FILE* f=stderr);

struct GXEditScript{
    uint32 *ops;               // array of edit operations
    uint32 opsize, opnum;         // size of allocation, number in use
    uint32 oplast;              // most recent operation added
    //methods

    GXEditScript() {
       init();
       }
    ~GXEditScript() {
       GFREE(ops);
       }
    void init() {
      ops = NULL;
      opsize = 0;
      opnum = 0;
      oplast = 0;
      getReady(8);
      }

    int getReady(uint32 n) {
        uint32 m = n + n/2;
        if (opsize <= n) {
            GREALLOC(ops, m*sizeof(uint32));
            opsize = m;
            }
        return 1;
       }

    int getReady2(uint32 n) {
        if (opsize - opnum <= n)
            return getReady(n + opnum);
        return 1;
       }

    int Put(uint32 op, uint32 n) {
        if (!getReady2(2))
            return 0;
        oplast = op;
        ops[opnum] = GX_EDITOP_CONS(op, n);
        opnum += 1;
        ops[opnum] = 0; // sentinel
        return 1;
        }
    uint32* First() {
        return opnum > 0 ? & ops[0] : NULL;
        }

    uint32* Next(uint32 *op) {
        //  assumes flat address space !
        if (&ops[0] <= op && op < &ops[opnum-1])
            return op+1;
        else
            return 0;
    }

   int More(uint32 op, uint32 k) {
        if (op == gxEDIT_OP_ERR) {
            GError("GXEditScript::opMore: bad opcode %d:%d", op, k);
            return -1;
            }

        if (GX_EDITOP_GET(oplast) == op) {
            uint32 l=ops[opnum-1];
            ops[opnum-1]=GX_EDITOP_CONS((GX_EDITOP_GET(l)),
               (GX_EDITOP_VAL(l) + k));
            }
        else {
            Put(op, k);
            }

        return 0;
    }

    GXEditScript* Append(GXEditScript *et) {
      uint32 *op;
      for (op = et->First(); op; op = et->Next(op))
        More(GX_EDITOP_GET(*op), GX_EDITOP_VAL(*op));
      return this;
    }

    int opDel(uint32 k) {
        return More(gxEDIT_OP_DEL, k);
        }
    int opIns(uint32 k) {
        return More(gxEDIT_OP_INS, k);
        }
    int opRep(uint32 k) {
        return More(gxEDIT_OP_REP, k);
        }

    GXEditScript *reverse() {
      const uint32 mid = opnum/2;
      const uint32 end = opnum-1;
      for (uint32 i = 0; i < mid; ++i) {
          const uint32 t = ops[i];
          ops[i] = ops[end-i];
          ops[end-i] = t;
          }
      return this;
      }
};


/** Bookkeeping structure for greedy alignment. When aligning
    two sequences, the members of this structure store the
    largest offset into the second sequence that leads to a
    high-scoring alignment for a given start point */
struct SGreedyOffset {
    int insert_off;    // Best offset for a path ending in an insertion
    int match_off;     // Best offset for a path ending in a match
    int delete_off;    // Best offset for a path ending in a deletion
};

// ----- pool allocator -----
// works as a linked list of allocated memory blocks
struct GXMemPool {
    SGreedyOffset* memblock;
    int used, size;
    GXMemPool *next;
    static int kMinSpace;
    //methods
    GXMemPool(int num_offsp=0) { //by default allocate a large block here (10M)
    	num_offsp=GMAX(kMinSpace, num_offsp);
        GMALLOC(memblock, num_offsp*sizeof(SGreedyOffset));
        if (memblock == NULL) {
           GError("Failed to allocated GXMemPool(%d) for greedy extension!\n",num_offsp);
           return;
           }
        used = 0;
        size = num_offsp;
        next = NULL;
        }

    void refresh() {
       GXMemPool* sp=this;
       while (sp) {
          sp->used = 0;
          sp = sp->next;
          }
       }
    ~GXMemPool() {
        GXMemPool* next_sp;
        GXMemPool* sp=this->next;
        while (sp) {
           next_sp = sp->next;
           GFREE(sp->memblock);
           delete sp;
           sp = next_sp;
           }
      GFREE(memblock);
      }

  SGreedyOffset* getSpace(int num_alloc) { // SGreedyOffset[num_alloc] array
     //can use the first found memory block with enough room,
    // or allocate a new large block
    SGreedyOffset* v;
    if (num_alloc < 0) return NULL;
    GXMemPool* S=this;
    while (used+num_alloc > S->size) {
        //no room in current block, get a new mem block
        if (next == NULL) {
           next=new GXMemPool(num_alloc); //allocates a large contiguous memory block
           }
        S = S->next;
        }
    v = S->memblock+S->used;
    S->used += num_alloc;
    //align to first 8-byte boundary
    int m8 = S->used & 7; //modulo 8
    if (m8)
         S->used += 8 - m8;
    return v;
    }

  void* getByteSpace(int byte_size) { //amount to use or allocate memory, in bytes
    return (void*)getSpace(byte_size/sizeof(SGreedyOffset));
    }

};

#define GREEDY_MAX_COST_FRACTION 8
/* (was 2) sequence_length / (this number) is a measure of how hard the
  alignment code will work to find the optimal alignment; in fact
  this gives a worst case bound on the number of loop iterations */

#define GREEDY_MAX_COST 1000
// The largest diff distance (max indels+mismatches) to be examined for an optimal alignment
// (should be increased for large sequences)

#define GX_GALLOC_ERROR "Error: failed to allocate memory for greedy alignment!\n"

// all auxiliary memory needed for the greedy extension algorithm
class CGreedyAlignData {
   int d_diff;
   int max_d;
public:
   int** last_seq2_off;              // 2-D array of distances
   int* max_score;                   // array of maximum scores
   GXMemPool* space;                // local memory pool for SGreedyOffset structs
   //
   bool scaled; //scores are all x2
   int match_reward;
   int mismatch_penalty;
   int x_drop;
   int xdrop_ofs;
   // Allocate memory for the greedy gapped alignment algorithm
   CGreedyAlignData(int reward, int penalty, int xdrop) {
	  scaled=false;
	  xdrop_ofs = 0;
      //int max_d, diff_d;
      if (penalty<0) penalty=-penalty;
      if (reward % 2) {
         //scale params up
    	 scaled=true;
         match_reward = reward << 1;
         mismatch_penalty = (penalty << 1);
         x_drop = xdrop<<1;
         }
      else {
         match_reward=reward;
         mismatch_penalty = penalty;
         x_drop=xdrop;
         }
      xdrop_ofs=(x_drop + (match_reward>>1)) /
          (match_reward + mismatch_penalty) + 1;
      //if (gap_open == 0 && gap_extend == 0)
      //   gap_extend = (reward >> 1) + penalty;
      const int max_dbseq_length=255; //adjust this accordingly
      max_d = GMIN(GREEDY_MAX_COST,
                  (max_dbseq_length/GREEDY_MAX_COST_FRACTION + 1));

      last_seq2_off=NULL;   // 2-D array of distances
      max_score=NULL;       // array of maximum scores
      space=NULL;           // local memory pool for SGreedyOffset structs
      //if (score_params.gap_open==0 && score_params.gap_extend==0) {
         //non-affine, simpler Greedy algorithm
       d_diff = (x_drop+match_reward/2)/(mismatch_penalty+match_reward)+1;
       GMALLOC(last_seq2_off, ((max_d + 2) * sizeof(int*)));
       if (!last_seq2_off)
         GError(GX_GALLOC_ERROR);
       GCALLOC(last_seq2_off[0], ((max_d + max_d + 6) * sizeof(int) * 2));
       //allocates contiguous memory for 2 rows here
       if (!last_seq2_off[0])
               GError(GX_GALLOC_ERROR);
       last_seq2_off[1] = last_seq2_off[0] + max_d + max_d + 6; //memory allocated already for this row

      GCALLOC(max_score, (sizeof(int) * (max_d + 1 + d_diff)));
      space = new GXMemPool();
      if (!max_score || !space)
         GError(GX_GALLOC_ERROR);
    } //consructor

   void reset() {
     space->refresh();
     if (last_seq2_off) {
         GFREE((last_seq2_off[0]));
         }
     GFREE(max_score);
     GCALLOC(last_seq2_off[0], ((max_d + max_d + 6) * sizeof(int) * 2));
     if (!last_seq2_off[0]) GError(GX_GALLOC_ERROR);
     //allocates contiguous memory for 2 rows here
     last_seq2_off[1] = last_seq2_off[0] + max_d + max_d + 6;
     GCALLOC(max_score, (sizeof(int) * (max_d + 1 + d_diff)));
     if (!max_score) GError(GX_GALLOC_ERROR);
     }

   ~CGreedyAlignData() {
     if (last_seq2_off) {
         GFREE(last_seq2_off[0]);
         GFREE(last_seq2_off);
         }
     GFREE(max_score);
     delete space;
     }

};


#define GAPALIGN_SUB ((unsigned char)0)  /*op types within the edit script*/
#define GAPALIGN_INS ((unsigned char)1)
#define GAPALIGN_DEL ((unsigned char)2)
#define GAPALIGN_DECLINE ((unsigned char)3)

struct GapXEditScript {
    unsigned char op_type;  // GAPALIGN_SUB, GAPALIGN_INS, or GAPALIGN_DEL
    int num;                // Number of operations
    GapXEditScript* next;
    GapXEditScript() {
        op_type=0;
        num=0;
        next=NULL;
        }
    void print();
};

class CSeqGap { //
 public:
   int offset;
   int len;
   CSeqGap(int gofs=0,int glen=1) {
     offset=gofs;
     len=glen;
     }
};

class CAlnGapInfo {
    int a_ofs; //alignment start on seq a (0 based)
    int b_ofs; //alignment start on seq b (0 based)
    int a_len; //length of alignment on seq a
    int b_len; //length of alignment on seq b
  public:
    GVec<CSeqGap> a_gaps;
    GVec<CSeqGap> b_gaps;
    CAlnGapInfo(GXEditScript* ed_script, int astart=0, int bstart=0):a_gaps(),b_gaps() {
       a_ofs=astart;
       b_ofs=bstart;
       a_len=0;
	   b_len=0;
	   if (ed_script==NULL) return;
	   for (uint32 i=0; i<ed_script->opnum; i++) {
		  int num=((ed_script->ops[i]) >> 2);
		  char op_type = 3 - ( ed_script->ops[i] & gxEDIT_OP_MASK );
		  if (op_type == 3 || op_type < 0 )
			 GError("Error: encountered op_type %d in ed_script?!\n", (int)op_type);
		  CSeqGap gap;
		  switch (op_type) {
			 case GAPALIGN_SUB: a_len+=num;
								b_len+=num;
								break;
			 case GAPALIGN_INS: a_len+=num;
								gap.offset=b_ofs+b_len;
								gap.len=num;
								b_gaps.Add(gap);
								break;
			 case GAPALIGN_DEL: b_len+=num;
								gap.offset=a_ofs+a_len;
								gap.len=num;
								a_gaps.Add(gap);
								break;
			 }
		  }
	   }

#ifdef TRIMDEBUG
	void printAlignment(FILE* f, const char* sa, int sa_len,
		             const char* sb, int sb_len) {
		//print seq A
	   char al[1024]; //display buffer for seq A
	   int ap=0; //index in al[] for current character printed
	   int g=0;
	   int aend=a_ofs+a_len;
	   if (a_ofs<b_ofs) {
		   for (int i=0;i<b_ofs-a_ofs;i++) {
    		   fprintf(f, " ");
    		   al[++ap]=' ';
    	       }
           }
       //int curbg=c_blue;
       for (int i=0;i<aend;i++) {
         if (g<a_gaps.Count() && a_gaps[g].offset==i) {
            color_bg(c_red, f);
            for (int j=0;j<a_gaps[g].len;j++) {
            	 fprintf(f, "-");
            	 al[++ap]='-';
                 }
            color_bg(c_blue, f);
            g++;
            }
         if (i==a_ofs) color_bg(c_blue,f);
         fprintf(f, "%c", sa[i]);
         al[++ap]=sa[i];
         }
       color_reset(f);
       if (aend<sa_len)
         fprintf(f, &sa[aend]);
       fprintf(f, "\n");
       //print seq B
       ap=0;
       g=0;
       int bend=b_ofs+b_len;
       if (a_ofs>b_ofs) {
    	   for (int i=0;i<a_ofs-b_ofs;i++) {
    		   fprintf(f, " ");
    		   ap++;
    	       }
           }
       for (int i=0;i<b_ofs;i++) {
    	   fprintf(f, "%c", sb[i]);
    	   ap++;
           }
       for (int i=b_ofs;i<bend;i++) {
         if (g<b_gaps.Count() && b_gaps[g].offset==i) {
            color_bg(c_red,f);
            for (int j=0;j<b_gaps[g].len;j++) {
            	fprintf(f, "-");
            	ap++;
                }
            color_bg(c_blue,f);
            //curbg=c_blue;
            g++;
            }
         if (i==b_ofs) color_bg(c_blue,f);
         ap++;
         bool mismatch=(sb[i]!=al[ap] && al[ap]!='-');
         if (mismatch) color_bg(c_red,f);
         fprintf(f, "%c", sb[i]);
         if (mismatch) color_bg(c_blue,f);
         }
       color_reset(f);
       if (bend<sb_len)
         fprintf(f, &sb[bend]);
       fprintf(f, "\n");
       }
#endif
  };

struct GXAlnInfo {
 const char *qseq;
 int ql,qr;
 const char *sseq;
 int sl,sr;
 int score;
 double pid;
 bool strong;
 int udata;
 GXEditScript* editscript;
 CAlnGapInfo* gapinfo;
 GXAlnInfo(const char* q, int q_l, int q_r, const char* s, int s_l, int s_r,
      int sc=0, double percid=0) {
    qseq=q;
    sseq=s;
    ql=q_l;
    qr=q_r;
    sl=s_l;
    sr=s_r;
    score=sc;
    pid=percid;
    strong=false;
    udata=0;
    editscript=NULL;
    gapinfo=NULL;
    }
  ~GXAlnInfo() {
    delete editscript;
    delete gapinfo;
    }
  bool operator<(GXAlnInfo& d) {
    return ((score==d.score)? pid>d.pid : score>d.score);
    }
  bool operator==(GXAlnInfo& d) {
    return (score==d.score && pid==d.pid);
    }

};



struct GXSeed {
   int b_ofs; //0-based coordinate on seq b (x coordinate)
   int a_ofs; //0-based coordinate on seq a (y coordinate)
   int len; //length of exact match after extension
   bool operator<(GXSeed& d){
      return ((b_ofs==d.b_ofs) ? a_ofs<d.a_ofs : b_ofs<d.b_ofs);
      }
   bool operator==(GXSeed& d){
      return (b_ofs==d.b_ofs && a_ofs==d.a_ofs); //should never be the case, seeds are uniquely constructed
      }
   GXSeed(int aofs=0, int bofs=0, int l=4) {
     a_ofs=aofs;
     b_ofs=bofs;
     len=l;
     }
};

int cmpSeedDiag(const pointer p1, const pointer p2);
  //seeds are "equal" if they're on the same diagonal (for selection purposes only)

int cmpSeedScore(const pointer p1, const pointer p2); //also takes position into account
   //among seeds with same length, prefer those closer to the left end of the read (seq_b)

struct GXBand {
  //bundle of seed matches on 3 adjacent diagonals
  int diag; //first anti-diagonal (b_ofs-a_ofs) in this group of 3
   //seeds for this, and diag+1 and diag+2 are stored here
  int min_a, max_a; //maximal coordinates of the bundle
  int min_b, max_b;
  int w_min_b; //weighted average of left start coordinate
  int avg_len;
  GList<GXSeed> seeds; //sorted by x coordinate (b_ofs)
  int score; //sum of seed scores (- overlapping_bases/2 - gaps)
  bool tested;
  GXBand(int start_diag=-1, GXSeed* seed=NULL):seeds(true, false, false) {
	  diag=start_diag;
	  min_a=MAX_INT;
	  min_b=MAX_INT;
	  max_a=0;
	  max_b=0;
	  score=0;
	  avg_len=0;
	  w_min_b=0;
	  tested=false;
    if (seed!=NULL) addSeed(seed);
    }
  void addSeed(GXSeed* seed) {
     seeds.Add(seed);
     score+=seed->len;
     avg_len+=seed->len;
     w_min_b+=seed->b_ofs * seed->len;
     //if (diag<0) diag=seed->diag; //should NOT be done like this
     if (seed->a_ofs < min_a) min_a=seed->a_ofs;
     if (seed->a_ofs+ seed->len > max_a) max_a=seed->a_ofs+seed->len;
     if (seed->b_ofs < min_b) min_b=seed->b_ofs;
     if (seed->b_ofs+seed->len > max_b) max_b=seed->b_ofs+seed->len;
     }

  void finalize() {
	  //!! to be called only AFTER all seeds have been added
	  // seeds are sorted by b_ofs
	  //penalize seed gaps and overlaps on b sequence
    if (avg_len==0) return;
    w_min_b/=avg_len;
    avg_len>>=1;
	  for (int i=1;i<seeds.Count();i++) {
        GXSeed& sprev=*seeds[i-1];
        GXSeed& scur=*seeds[i];
        if (scur==sprev) GError("Error: duplicate seeds found (%d-%d:%d-%d)!\n",
        		      scur.a_ofs+1, scur.a_ofs+scur.len, scur.b_ofs+1, scur.b_ofs+scur.len);
        int b_gap=scur.b_ofs-sprev.b_ofs-sprev.len;
        int a_gap=scur.a_ofs-sprev.a_ofs-sprev.len;
        int max_gap=b_gap;
        int min_gap=a_gap;
        if (min_gap>max_gap) Gswap(max_gap, min_gap);
        int _penalty=0;
        if (min_gap<0) { //overlap
               if (max_gap>0) { _penalty=GMAX((-min_gap), max_gap); }
                  else _penalty=-min_gap;
               }
            else { //gap
               _penalty=max_gap;
               }
        score-=(_penalty>>1);
        //score-=_penalty;
 	      }//for each seed
     }

  //bands will be sorted by decreasing score eventually, after all seeds are added
  //more seeds better than one longer seed?
  bool operator<(GXBand& d){
     //return ((score==d.score) ? seeds.Count()>d.seeds.Count() : score>d.score);
     return ((score==d.score) ? w_min_b<d.w_min_b : score>d.score);
     }
  bool operator==(GXBand& d){
    //return (score==d.score && seeds.Count()==d.seeds.Count());
     return (score==d.score && w_min_b==d.w_min_b);
     }

};

class GXBandSet:public GList<GXBand> {
  public:
   GXSeed* qmatch; //long match (mismatches allowed) if a very good match was extended well
   GXSeed* tmatch_r; //terminal match to be used if there is no better alignment
   GXSeed* tmatch_l; //terminal match to be used if there is no better alignment
   int idxoffset; //global anti-diagonal->index offset (a_len-1)
   //used to convert a diagonal to an index
   //diagonal is always b_ofs-a_ofs, so the minimum value is -a_len+1
   //hence offset is a_len-1
   GXBand* band(int diag) { //retrieve the band for given anti-diagonal (b_ofs-a_ofs)
      return Get(diag+idxoffset);
      }
   GXBand* band(int a_ofs, int b_ofs) { //retrieve the band for given anti-diagonal (b_ofs-a_ofs)
      return Get(b_ofs-a_ofs+idxoffset);
      }
   GXBandSet(int a_len, int b_len):GList<GXBand>(a_len+b_len-1, false, true, false) {
      idxoffset=a_len-1;
      qmatch=NULL;
      tmatch_l=NULL; //terminal match to be used if everything else fails
      tmatch_r=NULL;
	  //diag will range from -a_len+1 to b_len-1, so after adjustment
	  //by idxoffset we get a max of a_len+b_len-2
      int bcount=a_len+b_len-1;
      for (int i=0;i<bcount;i++)
	      this->Add(new GXBand(i-idxoffset));
           //unsorted, this should set fList[i]
      }
   ~GXBandSet() {
      delete qmatch;
      }
   void addSeed(GXSeed* seed) {
	 //MUST be unsorted !!!
	 int idx=(seed->b_ofs-seed->a_ofs)+idxoffset;
     fList[idx]->addSeed(seed);
     if (idx>0) fList[idx-1]->addSeed(seed);
     if (idx<fCount-1) fList[idx+1]->addSeed(seed);
     }
};

inline int calc_safelen(int alen) {
	 int r=iround(double(alen*0.6));
	 return (r<22)? 22 : r;
  }

struct GXSeqData {
  const char* aseq;
  int alen;
  const char* bseq;
  int blen;
  GVec<uint16>** amers;
  int amlen; //minimum alignment length that's sufficient to
             //trigger the quick extension heuristics
  GXSeqData(const char* sa=NULL, int la=0, const char* sb=NULL, int lb=0,
  GVec<uint16>* mers[]=NULL):aseq(sa), alen(la),
     bseq(sb),  blen(lb), amers(mers), amlen(0) {
   calc_amlen();
   calc_bmlen();
   }

  void calc_amlen() {
    if (alen) {
       int ah=calc_safelen(alen);
       if (amlen>ah) amlen=ah;
       }
    }
  void calc_bmlen() {
    if (blen) {
      int bh = iround(double(blen)*0.6);
      if (bh<22) bh=22;
      if (amlen>bh) amlen=bh;
      }
    }
  void update(const char* sa, int la, GVec<uint16>** mers,
	  const char* sb, int lb, int mlen=0) {
     aseq=sa;
     alen=la;
     amers=mers;
     if (mlen) {
       amlen=mlen;
       }
       else calc_amlen();
     if (sb==bseq && blen==lb) return;
     bseq=sb;
     blen=lb;
     calc_bmlen();
     }
  /*
  void update_b(const char* sb, int lb) {
     if (sb==bseq && blen==lb) return;
     bseq=sb;
     blen=lb;
     calc_bmlen();
     }*/
};

uint16 get6mer(char* p);
void table6mers(const char* s, int slen, GVec<uint16>* amers[]);

void printEditScript(GXEditScript* ed_script);


int GXGreedyExtend(const char* seq1, int len1,
                       const char* seq2, int len2,
                       bool reverse, //int xdrop_threshold, int match_cost, int mismatch_cost,
                       int& seq1_align_len, int& seq2_align_len,
                       CGreedyAlignData& aux_data,
                       GXEditScript *edit_block);


enum GAlnTrimType {
  //Describes trimming intent
  galn_None=0, //no trimming, just alignment
  galn_TrimLeft,
  galn_TrimRight,
  galn_TrimEither //adapter should be trimmed from either end
  };

struct CAlnTrim {
  GAlnTrimType type;
  int minMatch; //minimum terminal exact match that will be removed from ends
  int l_boundary; //base index (either left or right) excluding terminal poly-A stretches
  int r_boundary; //base index (either left or right) excluding terminal poly-A stretches
  int alen; //query/adapter seq length (for validate())
  int safelen; //alignment length > amlen should be automatically validated
  int seedlen;
  void prepare(const char* s, int s_len) {
    //type=trim_type;
    //amlen=smlen;
    l_boundary=0;
    r_boundary=0;
    //if (type==galn_TrimLeft) {
        int s_lbound=0;
        if (s[0]=='A' && s[1]=='A' && s[2]=='A') {
           s_lbound=3;
           while (s_lbound<s_len-1 && s[s_lbound]=='A') s_lbound++;
           }
        else if (s[1]=='A' && s[2]=='A' && s[3]=='A') {
           s_lbound=4;
           while (s_lbound<s_len-1 && s[s_lbound]=='A') s_lbound++;
           }
        l_boundary=s_lbound+3;
    //    return;
    //    }
    //if (type==galn_TrimRight) {
       int r=s_len-1;
       if (s[r]=='A' && s[r-1]=='A' && s[r-2]=='A') {
          r-=3;
          while (r>0 && s[r]=='A') r--;
          }
       else if (s[r-1]=='A' && s[r-2]=='A' && s[r-3]=='A') {
          r-=4;
          while (r>0 && s[r]=='A') r--;
          }
       r_boundary=r-3;
    //   }
    }

  CAlnTrim(GAlnTrimType trim_type, const char* s, int s_len, int a_len, int minEndTrim, int smlen):
	           type(trim_type), minMatch(minEndTrim), l_boundary(0), r_boundary(0),
	           alen(a_len), safelen(smlen) {
    prepare(s, s_len);
    }

  bool validate_R(int sr, int admax, int badj, int adist) {
	if (adist>admax) return false;
	return (sr>=r_boundary+badj);
   }

  bool validate_L(int sl, int alnlen, int admax, int badj, int alnpid, int adist) {
	if (adist>admax) return false;
    //left match should be more stringent (5')
    if (alnpid<93) {
      if (alnlen<13 || alnlen<minMatch) return false;
      admax=0;
      badj++;
      }
    return (sl<=l_boundary-badj);
  }

  bool validate(GXAlnInfo* alninfo) {
   int alnlen=alninfo->sr - alninfo->sl + 1;
/*   #ifdef TRIMDEBUG
     GMessage("\t::: alnlen=%d, safelen=%d, pid=%4.2f\n",
          alnlen, safelen, alninfo->pid);
   #endif
*/
   if (alninfo->pid>90.0 && alnlen>=safelen) {
	   //special case: heavy match, could be in the middle
	   if (alninfo->pid>94)
		    alninfo->strong=true;
	   return true;
       }
   int sl=alninfo->sl;
   int sr=alninfo->sr;
   sl--;sr--; //boundary is 0-based
   int badj=0; //default boundary is 3 bases distance to end
   int admax=1;
   if (alnlen>20) ++admax;
   if (alnlen<13) {
      //stricter boundary check
      if (alninfo->pid<90) return false;
      badj=2;
      if (alnlen<=7) { badj++; admax=0; }
      }
   if (type==galn_TrimLeft) {
	 return validate_L(sl, alnlen, admax, badj, alninfo->pid, alen-alninfo->qr);
     }
   else if (type==galn_TrimRight) {
	 return validate_R(sr, admax, badj, alninfo->ql-1);
     }
   else if (type==galn_TrimEither) {
     return (validate_R(sr, admax, badj, alninfo->ql-1) ||
    	   validate_L(sl, alnlen, admax, badj, alninfo->pid, alen-alninfo->qr));
     }
   return true;
   /*
   if (type==galn_TrimRight) {
      return (sr>=boundary+badj);
      }
   else {
      //left match should be more stringent (5')
      if (alnpid<93) {
        if (alnlen<13) return false;
        admax=0;
        badj++;
        }
      return (sl<=boundary-badj);
      }
    */
   }
};




//GXBandSet* collectSeeds_R(GList<GXSeed>& seeds, GXSeqData& sd); //for overlap at 3' end of seqb

GXBandSet* collectSeeds(GList<GXSeed>& seeds, GXSeqData& sd, GAlnTrimType trim_intent); //for overlap at 5' end of seqb
                                                                //=galn_None
// reward MUST be >1 for this function
GXAlnInfo* GreedyAlignRegion(const char* q_seq, int q_alnstart, int q_max,
                  const char* s_seq, int s_alnstart, int s_max,
                  int reward, int penalty, int xdrop, CGreedyAlignData* gxmem=NULL,
                  CAlnTrim* trim=NULL, bool editscript=false);
GXAlnInfo* GreedyAlignRegion(const char* q_seq, int q_alnstart, int q_max,
                       const char* s_seq, int s_alnstart, int s_max, CGreedyAlignData* gxmem,
                       CAlnTrim* trim=NULL, bool editscript=false);

GXAlnInfo* GreedyAlign(const char* q_seq,  int q_alnstart, const char* s_seq, int s_alnstart,
        bool editscript=false, int reward=2, int penalty=10, int xdrop=32);

GXAlnInfo* match_adapter(GXSeqData& sd, GAlnTrimType trim_type, int minMatch,
	                        CGreedyAlignData* gxmem=NULL, double min_pid=90);
//GXAlnInfo* match_RightEnd(GXSeqData& sd, CGreedyAlignData* gxmem=NULL, int min_pid=90);
#endif
