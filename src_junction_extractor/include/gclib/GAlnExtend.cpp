#include "GAlnExtend.h"

//greedy gapped alignment extension
//(mostly lifted from NCBI's megablast gapped extension code)

int GXMemPool::kMinSpace = 1000000;

// ifdef TRIMDEBUG
 char COLOR_buf[6]={0x1B,'[', 'n','m','m','\0'};

void color_fg(int c,FILE* f) {
 if (f!=stderr && f!=stdout) return;
 sprintf((char *)(&COLOR_buf[2]),"%dm",c+30);
 fwrite(COLOR_buf,1,strlen(COLOR_buf), f);
}

void color_bg(int c, FILE* f) {
 if (f!=stderr && f!=stdout) return;
 sprintf((char *)(&COLOR_buf[2]),"%dm",c+40);
 fwrite(COLOR_buf,1,strlen(COLOR_buf),f);
};

void color_resetfg(FILE* f) {
 if (f!=stderr && f!=stdout) return;
 sprintf((char *)(&COLOR_buf[2]),"39m");
 fwrite(COLOR_buf,1,strlen(COLOR_buf), f);
};

void color_resetbg(FILE* f) {
 if (f!=stderr && f!=stdout) return;
 sprintf((char *)(&COLOR_buf[2]),"49m");
 fwrite(COLOR_buf,1,strlen(COLOR_buf), f);
}

void color_reset(FILE* f) {
 if (f!=stderr && f!=stdout) return;
 sprintf((char *)(&COLOR_buf[2]),"0m");
 fwrite(COLOR_buf,1,strlen(COLOR_buf), f);
};

void color_normal(FILE* f) {
 if (f!=stderr && f!=stdout) return;
 sprintf((char *)(&COLOR_buf[2]),"22m");
 fwrite(COLOR_buf,1,strlen(COLOR_buf), f);
};

// endif


char xgapcodes[4]={'S','I', 'D', 'X'};

int get_last(int **flast_d, int d, int diag, int *row1) {
  if (flast_d[d-1][diag-1] > GMAX(flast_d[d-1][diag], flast_d[d-1][diag+1])) {
      *row1 = flast_d[d-1][diag-1];
      return diag-1;
      }
  if (flast_d[d-1][diag] > flast_d[d-1][diag+1]) {
      *row1 = flast_d[d-1][diag];
      return diag;
      }
  *row1 = flast_d[d-1][diag+1];
  return diag+1;
}

void GapXEditScript::print() { //debug
      GapXEditScript* p=this;
      do {
        GMessage("%d%c ",p->num, xgapcodes[p->op_type]);
        } while ((p=p->next)!=NULL);
      GMessage("\n");
    }


int BLAST_Gcd(int a, int b) {
   int   c;

   b = abs(b);
   if (b > a)
      c=a, a=b, b=c;

   while (b != 0) {
      c = a%b;
      a = b;
      b = c;
   }
   return a;
}

int BLAST_Gdb3(int* a, int* b, int* c) {
    int g;
    if (*b == 0)
        g = BLAST_Gcd(*a, *c);
    else
        g = BLAST_Gcd(*a, BLAST_Gcd(*b, *c));
    if (g > 1) {
		*a /= g;
        *b /= g;
        *c /= g;
    }
    return g;
}


uint16 get6mer(char* p) {
  uint16 r=gdna2bit(p,3);
  r <<= 6;
  r |= gdna2bit(p,3);
  return r;
}


void table6mers(const char* s, int slen, GVec<uint16>* amers[]) {
 for (uint16 i=0; i <= slen-6; i++) {
   char* p = (char*)(s+i);
   uint16 v=get6mer(p);
   if (amers[v]==NULL) {
      amers[v]=new GVec<uint16>(1);
      }
   amers[v]->Add(i);
 }
}

GVec<uint16>* match6mer(char* start, GVec<uint16>* amers[]) {
  //careful: this is broken if start+5 falls beyond the end of the string!
  uint16 r=get6mer(start);
  return amers[r];
}

//signal that a diagonal is invalid
static const int kInvalidOffset = -2;

int s_FindFirstMismatch(const char *seq1, int len1,
                        const char *seq2, int len2,
                           int seq1_index, int seq2_index,
                              //bool &fence_hit,
                               bool reverse) {
    int start_index = seq1_index;
    /* Sentry detection here should be relatively inexpensive: The
       sentry value cannot appear in the query, so detection only
       needs to be done at exit from the subject-query matching loop.
       For uncompressed sequences, ambiguities in the query (i.e. seq1)
       always count as mismatches */
    if (reverse) {
              while (seq1_index < len1 && seq2_index < len2 &&
                     //seq1[len1-1 - seq1_index] < 4 &&
                       seq1[len1-1 - seq1_index] == seq2[len2-1 - seq2_index]) {
                    ++seq1_index;
                    ++seq2_index;
                    }
            //if (seq2_index < len2 && seq2[len2-1-seq2_index] == FENCE_SENTRY) {
              //if len2-1-seq2_index<=0) {
              //   fence_hit = true;
             //    }
      }
    else { //forward lookup
             while (seq1_index < len1 && seq2_index < len2 &&
                   //seq1[seq1_index] < 4 &&
                   seq1[seq1_index] == seq2[seq2_index]) {
                ++seq1_index;
                ++seq2_index;
                }
            //if (seq2_index < len2 && seq2[seq2_index] == FENCE_SENTRY) {
             //if (seq2_index==len2) {
             //   fence_hit = true;
            //}
    }
    return seq1_index - start_index;
}



/** During the traceback for a non-affine greedy alignment,
    compute the diagonal that will result from the next
    traceback operation

    @param last_seq2_off Array of offsets into the second sequence;
                        last_seq2_off[d][k] gives the largest offset into
                        the second sequence that lies on diagonal k and
                        has distance d [in]
    @param d Starting distance [in]
    @param diag Index of diagonal that produced the starting distance [in]
    @param seq2_index The offset into the second sequence after the traceback
                operation has completed [out]
    @return The diagonal resulting from the next traceback operation
                being applied
*/
int s_GetNextNonAffineTback(int **last_seq2_off, int d,
                        int diag, int *seq2_index) {
    // choose the traceback operation that results in the
    //   largest seq2 offset at this point, then compute the
    //   new diagonal that is implied by the operation
    if (last_seq2_off[d-1][diag-1] >
                GMAX(last_seq2_off[d-1][diag], last_seq2_off[d-1][diag+1])) {
        *seq2_index = last_seq2_off[d-1][diag-1];
        return diag - 1;    // gap in seq2
    }
    if (last_seq2_off[d-1][diag] > last_seq2_off[d-1][diag+1]) {
        *seq2_index = last_seq2_off[d-1][diag];
        return diag;        // match
    }
    *seq2_index = last_seq2_off[d-1][diag+1];
    return diag + 1;        // gap in seq1
}


int GXGreedyExtend(const char* seq1, int len1,
                       const char* seq2, int len2,
                       bool reverse, //int xdrop_threshold, int match_cost, int mismatch_cost,
                       int& seq1_align_len, int& seq2_align_len,
                       CGreedyAlignData& aux_data,
                       GXEditScript *edit_block) {
                       //GapPrelimEditBlock *edit_block,
                       //bool& fence_hit, SGreedySeed *seed) {
  int seq1_index;
    int seq2_index;
    int index;
    int d;
    int k;
    int diag_lower, diag_upper;
    int max_dist;
    int diag_origin;
    int best_dist;
    int best_diag;
    int** last_seq2_off;
    int* max_score;
    int xdrop_offset;
    int longest_match_run;
    bool end1_reached, end2_reached;
    GXMemPool* mem_pool;

    /* ordinary dynamic programming alignment, for each offset
       in seq1, walks through offsets in seq2 until an X-dropoff
       test fails, saving the best score encountered along
       the way. Instead of score, this code tracks the 'distance'
       (number of mismatches plus number of gaps) between seq1
       and seq2. Instead of walking through sequence offsets, it
       walks through diagonals that can achieve a given distance.

       Note that in what follows, the numbering of diagonals implies
       a dot matrix where increasing seq1 offsets go to the right on
       the x axis, and increasing seq2 offsets go up the y axis.
       The gapped alignment thus proceeds up and to the right in
       the graph, and diagonals are numbered increasing to the right */

    best_dist = 0;
    best_diag = 0;

    /* set the number of distinct distances the algorithm will
       examine in the search for an optimal alignment. The
       heuristic worst-case running time of the algorithm is
       O(max_dist**2 + (len1+len2)); for sequences which are
       very similar, the average running time will be sig-
       nificantly better than this */

    max_dist = GMIN(GREEDY_MAX_COST,
                   (len2/GREEDY_MAX_COST_FRACTION + 1));

    /* the main loop assumes that the index of all diagonals is
       biased to lie in the middle of allocated bookkeeping
       structures */

    diag_origin = max_dist + 2;

    // last_seq2_off[d][k] is the largest offset into seq2 that
    //   lies on diagonal k and has distance d

    last_seq2_off = aux_data.last_seq2_off;

    /* Instead of tracking the best alignment score and using
       xdrop_theshold directly, track the best score for each
       unique distance and use the best score for some previously
       computed distance to implement the X-dropoff test.

       xdrop_offset gives the distance backwards in the score
       array to look */

    xdrop_offset = aux_data.xdrop_ofs;

    // find the offset of the first mismatch between seq1 and seq2

    index = s_FindFirstMismatch(seq1, len1,  seq2, len2, 0, 0, reverse);
                               // fence_hit, reverse, rem);

    // update the extents of the alignment, and bail out
    //   early if no further work is needed

    seq1_align_len = index;
    seq2_align_len = index;
    seq1_index = index;
    /*
    seed->start_q = 0;
    seed->start_s = 0;
    seed->match_length = index;
    */
    longest_match_run = index;

    if (index == len1 || index == len2) {
        /* Return the number of differences, which is zero here */
        if (edit_block != NULL)
            //GapPrelimEditBlockAdd(edit_block, eGapAlignSub, index);
            edit_block->opRep(index);
        return 0;
    }

    // set up the memory pool
    mem_pool = aux_data.space;
    if (edit_block == NULL) {
       mem_pool = NULL;
      }
    else if (mem_pool == NULL) {
       aux_data.space = mem_pool = new GXMemPool();
    }
    else {
        mem_pool->refresh();
    }

    /* set up the array of per-distance maximum scores. There
       are max_diags + xdrop_offset distances to track, the first
       xdrop_offset of which are 0 */

    max_score = aux_data.max_score + xdrop_offset;
    for (index = 0; index < xdrop_offset; index++)
        aux_data.max_score[index] = 0;

    // fill in the initial offsets of the distance matrix

    last_seq2_off[0][diag_origin] = seq1_index;
    max_score[0] = seq1_index * aux_data.match_reward;
    diag_lower = diag_origin - 1;
    diag_upper = diag_origin + 1;
    end1_reached = end2_reached = false;

    // for each distance
    for (d = 1; d <= max_dist; d++) {
        int xdrop_score;
        int curr_score;
        int curr_extent = 0;
        int curr_seq2_index = 0;
        int curr_diag = 0;
        int tmp_diag_lower = diag_lower;
        int tmp_diag_upper = diag_upper;

        // Assign impossible seq2 offsets to any diagonals that
        // are not in the range (diag_lower,diag_upper).
        // These will serve as sentinel values for the inner loop
        last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
        last_seq2_off[d - 1][diag_lower] = kInvalidOffset;
        last_seq2_off[d - 1][diag_upper] = kInvalidOffset;
        last_seq2_off[d - 1][diag_upper+1] = kInvalidOffset;

        // compute the score for distance d corresponding to the X-dropoff criterion

        xdrop_score = max_score[d - xdrop_offset] +
                      (aux_data.match_reward + aux_data.mismatch_penalty) * d - aux_data.x_drop;
        xdrop_score = (int)ceil((double)xdrop_score / (aux_data.match_reward>>1));

        // for each diagonal of interest
        for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
            /* find the largest offset into seq2 that increases
               the distance from d-1 to d (i.e. keeps the alignment
               from getting worse for as long as possible), then
               choose the offset into seq1 that will keep the
               resulting diagonal fixed at k

               Note that this requires kInvalidOffset+1 to be smaller
               than any valid offset into seq2, i.e. to be negative */

            seq2_index = GMAX(last_seq2_off[d - 1][k + 1],
                              last_seq2_off[d - 1][k    ]) + 1;
            seq2_index = GMAX(seq2_index, last_seq2_off[d - 1][k - 1]);
            seq1_index = seq2_index + k - diag_origin;

            if (seq2_index < 0 || seq1_index + seq2_index < xdrop_score) {

                // if no valid diagonal can reach distance d, or the
                //   X-dropoff test fails, narrow the range of diagonals
                //   to test and skip to the next diagonal
                if (k == diag_lower)
                    diag_lower++;
                else
                    last_seq2_off[d][k] = kInvalidOffset;
                continue;
            }
            diag_upper = k;

            /* slide down diagonal k until a mismatch
               occurs. As long as only matches are encountered,
               the current distance d will not change */

            index = s_FindFirstMismatch(seq1, len1, seq2, len2,
                                        seq1_index, seq2_index, reverse);
                                        //fence_hit, reverse, rem);
            if (index > longest_match_run) {
                //seed->start_q = seq1_index;
                //seed->start_s = seq2_index;
                //seed->match_length = index;
                longest_match_run = index;
                }
            seq1_index += index;
            seq2_index += index;

            // set the new largest seq2 offset that achieves
            //   distance d on diagonal k

            last_seq2_off[d][k] = seq2_index;

            // since all values of k are constrained to have the
            // same distance d, the value of k which maximizes the
            // alignment score is the one that covers the most of seq1 and seq2
            if (seq1_index + seq2_index > curr_extent) {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
               }

            /* clamp the bounds on diagonals to avoid walking off
               either sequence. Because the bounds increase by at
               most one for each distance, diag_lower and diag_upper
               can each be of size at most max_diags+2 */

            if (seq2_index == len2) {
                diag_lower = k + 1;
                end2_reached = true;
            }
            if (seq1_index == len1) {
                diag_upper = k - 1;
                end1_reached = true;
            }
        }  // end loop over diagonals

        // compute the maximum score possible for distance d
        curr_score = curr_extent * (aux_data.match_reward / 2) -
                        d * (aux_data.match_reward + aux_data.mismatch_penalty);
        // if this is the best score seen so far, update the
        // statistics of the best alignment
        if (curr_score > max_score[d - 1]) {
            max_score[d] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            seq2_align_len = curr_seq2_index;
            seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        }
        else {
            max_score[d] = max_score[d - 1];
        }

        // alignment has finished if the lower and upper bounds
        //   on diagonals to check have converged to each other

        if (diag_lower > diag_upper)
            break;

        /* set up for the next distance to examine. Because the
           bounds increase by at most one for each distance,
           diag_lower and diag_upper can each be of size at
           most max_diags+2 */

        if (!end2_reached)
            diag_lower--;
        if (!end1_reached)
            diag_upper++;

        if (edit_block == NULL) {
           // if no traceback is specified, the next row of
           //   last_seq2_off can reuse previously allocated memory
           //WARNING The following assumes two arrays of
           //  at least max_dist+4 int's have already been allocated
            last_seq2_off[d + 1] = last_seq2_off[d - 1];
            }
        else {
            // traceback requires all rows of last_seq2_off to be saved,
            // so a new row must be allocated
            last_seq2_off[d + 1] = (int*)mem_pool->getByteSpace((diag_upper - diag_lower + 7)*sizeof(int));
            // move the origin for this row backwards
            // dubious pointer arithmetic ?!
            //last_seq2_off[d + 1] = last_seq2_off[d + 1] - diag_lower + 2;
            }
    }   // end loop over distinct distances


    if (edit_block == NULL)
        return best_dist;

    //----  perform traceback
    d = best_dist;
    seq1_index = seq1_align_len;
    seq2_index = seq2_align_len;
    // for all positive distances

    //if (fence_hit && *fence_hit)
    //    goto done;
    if (index==len1 || index==len2) d=0;
    while (d > 0) {
        int new_diag;
        int new_seq2_index;

        /* retrieve the value of the diagonal after the next
           traceback operation. best_diag starts off with the
           value computed during the alignment process */

        new_diag = s_GetNextNonAffineTback(last_seq2_off, d,
                                           best_diag, &new_seq2_index);

        if (new_diag == best_diag) {
            // same diagonal: issue a group of substitutions
            if (seq2_index - new_seq2_index > 0) {
                  edit_block->opRep(seq2_index - new_seq2_index);
            }
        }
        else if (new_diag < best_diag) {
            // smaller diagonal: issue a group of substitutions
            //   and then a gap in seq2 */
            if (seq2_index - new_seq2_index > 0) {
              edit_block->opRep(seq2_index - new_seq2_index);
            }
            //GapPrelimEditBlockAdd(edit_block, eGapAlignIns, 1);
            edit_block->opIns(1);
        }
        else {
            // larger diagonal: issue a group of substitutions
            //   and then a gap in seq1
            if (seq2_index - new_seq2_index - 1 > 0) {
                edit_block->opRep(seq2_index - new_seq2_index - 1);
            }
            edit_block->opDel(1);
        }
        d--;
        best_diag = new_diag;
        seq2_index = new_seq2_index;
    }
//done:
    // handle the final group of substitutions back to distance zero,
    //   i.e. back to offset zero of seq1 and seq2
    //GapPrelimEditBlockAdd(edit_block, eGapAlignSub,
    //                      last_seq2_off[0][diag_origin]);
    edit_block->opRep(last_seq2_off[0][diag_origin]);
    if (!reverse)
      edit_block->reverse();
    return best_dist;
}

void printEditScript(GXEditScript* ed_script) {
   uint i;
   if (ed_script==NULL || ed_script->opnum == 0)
      return;
   for (i=0; i<ed_script->opnum; i++) {
      int num=((ed_script->ops[i]) >> 2);
      unsigned char op_type = 3 - ( ed_script->ops[i] & gxEDIT_OP_MASK );
      if (op_type == 3)
         GError("Error: printEditScript encountered op_type 3 ?!\n");
      GMessage("%d%c ", num, xgapcodes[op_type]);
      }
    GMessage("\n");
  }

GXAlnInfo* GreedyAlign(const char* q_seq,  int q_alnstart, const char* s_seq, int s_alnstart,
        bool editscript, int reward, int penalty, int xdrop) {
 int q_max=strlen(q_seq); //query
 int s_max=strlen(s_seq); //subj
 return GreedyAlignRegion(q_seq, q_alnstart, q_max,
          s_seq,  s_alnstart, s_max, reward, penalty, xdrop, NULL, NULL, editscript);
}

struct GXSeedTable {
  int a_num, b_num;
  int a_cap, b_cap;
  char* xc;
  GXSeedTable(int a=12, int b=255) {
    a_cap=0;
    b_cap=0;
    a_num=0;
    b_num=0;
    xc=NULL;
    init(a,b);
    }
  ~GXSeedTable() {
     GFREE(xc);
     }
  void init(int a, int b) {
    a_num=a;
    b_num=b;
    bool resize=false;
    if (b_num>b_cap) { resize=true; b_cap=b_num;}
    if (a_num>a_cap) { resize=true; a_cap=a_num;}
    if (resize) {
      GFREE(xc);
      GCALLOC(xc, (a_num*b_num));
      }
     else {
      //just clear up to a_max, b_max
      memset((void*)xc, 0, (a_num*b_num));
      }
    }
   char& x(int ax, int by) {
       return xc[by*a_num+ax];
       }

};

const int a_m_score=2; //match score
const int a_mis_score=-3; //mismatch
const int a_dropoff_score=7;
const int a_min_score=12; //at least 6 bases full match

// ------------------ adapter matching - simple k-mer seed & extend, no indels for now
//when a k-mer match is found, simply try to extend the alignment using a drop-off scheme
//check minimum score and
//for 3' adapter trimming:
//     require that the right end of the alignment for either the adapter OR the read must be
//     < 3 distance from its right end
// for 5' adapter trimming:
//     require that the left end of the alignment for either the adapter OR the read must
//     be at coordinate < 3 from start

bool extendUngapped(const char* a, int alen, int ai,
                 const char* b, int blen, int bi, int& mlen, int& l5, int& l3, bool end5=false) {
 //so the alignment starts at ai in a, bi in b, with a perfect match of length mlen
 //if (debug) {
 //  GMessage(">> in %s\n\textending hit: %s at position %d\n", a, (dbg.substr(bi, mlen)).chars(), ai);
 //  }
 int a_l=ai; //alignment coordinates on a
 int a_r=ai+mlen-1;
 int b_l=bi; //alignment coordinates on b
 int b_r=bi+mlen-1;
 int ai_maxscore=ai;
 int bi_maxscore=bi;
 int score=mlen*a_m_score;
 int maxscore=score;
 int mism5score=a_mis_score;
 if (end5 && ai<(alen>>1)) mism5score-=2; // increase penalty for mismatches at 5' end
 //try to extend to the left first, if possible
 while (ai>0 && bi>0) {
   ai--;
   bi--;
   score+= (a[ai]==b[bi])? a_m_score : mism5score;
   if (score>maxscore) {
       ai_maxscore=ai;
       bi_maxscore=bi;
       maxscore=score;
       }
     else if (maxscore-score>a_dropoff_score) break;
   }
 a_l=ai_maxscore;
 b_l=bi_maxscore;
 //now extend to the right
 ai_maxscore=a_r;
 bi_maxscore=b_r;
 ai=a_r;
 bi=b_r;
 score=maxscore;
 //sometimes there are extra As at the end of the read, ignore those
 if (a[alen-2]=='A' && a[alen-1]=='A') {
    alen-=2;
    while (a[alen-1]=='A' && alen>ai) alen--;
    }
 while (ai<alen-1 && bi<blen-1) {
   ai++;
   bi++;
   //score+= (a[ai]==b[bi])? a_m_score : a_mis_score;
   if (a[ai]==b[bi]) { //match
      score+=a_m_score;
      if (ai>=alen-2) {
           score+=a_m_score-(alen-ai-1);
           }
      }
    else { //mismatch
      score+=a_mis_score;
      }
   if (score>maxscore) {
       ai_maxscore=ai;
       bi_maxscore=bi;
       maxscore=score;
       }
     else if (maxscore-score>a_dropoff_score) break;
   }
  a_r=ai_maxscore;
  b_r=bi_maxscore;
  int a_ovh3=alen-a_r-1;
  int b_ovh3=blen-b_r-1;
  int mmovh3=(a_ovh3<b_ovh3)? a_ovh3 : b_ovh3;
  int mmovh5=(a_l<b_l)? a_l : b_l;
  if (maxscore>=a_min_score && mmovh3<2 && mmovh5<2) {
     if (a_l<a_ovh3) {
        //adapter closer to the left end (typical for 5' adapter)
        l5=a_r+1;
        l3=alen-1;
        }
      else {
        //adapter matching at the right end (typical for 3' adapter)
        l5=0;
        l3=a_l-1;
        }
     return true;
     }
  //do not trim:
  l5=0;
  l3=alen-1;
  return false;
 }

GXBandSet* collectSeeds(GList<GXSeed>& seeds, GXSeqData& sd, GAlnTrimType trim_intent) {
 int bimin=GMAX(0,(sd.blen-sd.alen-6)); //from collectSeeds_R
 int bimax=GMIN((sd.alen+2), (sd.blen-6));
 int b_start = (trim_intent==galn_TrimRight) ? bimin : 0;
 int b_end = (trim_intent==galn_TrimLeft) ?  bimax : sd.blen-6;
 //gx.init(a_maxlen, b_maxlen);
 GXSeedTable gx(sd.alen, sd.blen);
 GXBandSet* diagstrips=new GXBandSet(sd.alen, sd.blen); //set of overlapping 3-diagonal strips
 for (int bi=b_start;bi<=b_end;bi++) {
   //for each 6-mer of seqb
   uint16 bv = get6mer((char*) & (sd.bseq[bi]));
   GVec<uint16>* alocs = sd.amers[bv];
   if (alocs==NULL) continue;
   //extend each hit
   for (int h=0;h<alocs->Count();h++) {
	 int ai=alocs->Get(h);
	 //word match
	 if (gx.x(ai,bi))
	   //already have a previous seed covering this region of this diagonal
	   continue;
	 if (trim_intent==galn_TrimLeft && sd.blen>sd.alen+6 && bi>ai+6)
	       continue; //improper positioning for 5' trimming
	 if (trim_intent==galn_TrimRight && sd.blen>sd.alen+6 && bi<ai-6)
	       continue; //improper positioning for 5' trimming

	 //TODO: there could be Ns in this seed, should we count/adjust score?
	 for (int i=0;i<6;i++)
	    gx.x(ai+i,bi+i)=1;
	 //see if we can extend to the right
	 int aix=ai+6;
	 int bix=bi+6;
	 int len=6;
	 while (bix<sd.blen && aix<sd.alen && sd.aseq[aix]==sd.bseq[bix]) {
		 gx.x(aix,bix)=1;
		 aix++;bix++;
		 len++;
		 }
	 if (len>sd.amlen) {
		//heuristics: very likely the best we can get
		//quick match shortcut
		diagstrips->qmatch=new GXSeed(ai,bi,len);
		return diagstrips;
		}
	 if (bi>bimax && bi<bimin && len<9)
		 //ignore mid-sequence seeds that are not high scoring
	     continue;

	 GXSeed* newseed=new GXSeed(ai,bi,len);
	 seeds.Add(newseed);
	 diagstrips->addSeed(newseed);//add it to all 3 adjacent diagonals
     //keep last resort terminal match to be used if no better alignment is there
     if (bi<2 && ai+len>=sd.alen-1 &&
		 (!diagstrips->tmatch_l || diagstrips->tmatch_l->len<len))
		      diagstrips->tmatch_l=newseed;
     //collectSeeds_R:
	 if (ai<2 && bi+len>sd.blen-2 &&
		 (!diagstrips->tmatch_r || diagstrips->tmatch_r->len<len))
	          diagstrips->tmatch_r=newseed;
     }
 } //for each 6-mer of the read
 for (int i=0;i<diagstrips->Count();i++) {
    diagstrips->Get(i)->finalize(); //adjust scores due to overlaps or gaps between seeds
    }
 diagstrips->setSorted(true); //sort by score
 return diagstrips;
}

int cmpSeedScore(const pointer p1, const pointer p2) {
  //return (((GXSeed*)s2)->len-((GXSeed*)s1)->len);
  GXSeed* s1=(GXSeed*)p1;
  GXSeed* s2=(GXSeed*)p2;
  if (s1->len==s2->len) {
     return (s1->b_ofs-s2->b_ofs);
     }
  else return (s2->len-s1->len);
}

int cmpSeedScore_R(const pointer p1, const pointer p2) {
  //return (((GXSeed*)s2)->len-((GXSeed*)s1)->len);
  GXSeed* s1=(GXSeed*)p1;
  GXSeed* s2=(GXSeed*)p2;
  if (s1->len==s2->len) {
     return (s2->b_ofs-s1->b_ofs);
     }
  else return (s2->len-s1->len);
}


int cmpSeedDiag(const pointer p1, const pointer p2) {
  GXSeed* s1=(GXSeed*)p1;
  GXSeed* s2=(GXSeed*)p2;
  return ((s1->b_ofs-s1->a_ofs)-(s2->b_ofs-s2->a_ofs));
}


int cmpDiagBands_R(const pointer p1, const pointer p2) {
  //return (((GXSeed*)s2)->len-((GXSeed*)s1)->len);
  GXBand* b1=(GXBand*)p1;
  GXBand* b2=(GXBand*)p2;
  if (b1->score==b2->score) {
     return (b2->w_min_b-b1->w_min_b);
     }
  else return (b2->score-b1->score);
}



GXAlnInfo* GreedyAlignRegion(const char* q_seq, int q_alnstart, int q_max,
    const char* s_seq, int s_alnstart, int s_max,
    int reward, int penalty, int xdrop, CGreedyAlignData* gxmem,
    CAlnTrim* trim, bool editscript) {
  GXEditScript* ed_script_fwd = NULL;
  GXEditScript* ed_script_rev = NULL;
  if ( q_alnstart>q_max || q_alnstart<1 || s_alnstart>s_max || s_alnstart<1 )
     GError("GreedyAlign() Error: invalid anchor coordinate.\n");
  q_alnstart--;
  s_alnstart--;
  if (q_seq==NULL || q_seq[0]==0 || s_seq==NULL || s_seq[0]==0)
    GError("GreedyAlign() Error: attempt to use an empty sequence string!\n");
  /*if (q_seq[q_alnstart]!=s_seq[s_alnstart])
    GError("GreedyAlign() Error: improper anchor (mismatch):\n%s (start %d len %d)\n%s (start %d len %d)\n",
    	   q_seq, q_alnstart, q_max, s_seq, s_alnstart, s_max);
    	   */
  int q_ext_l=0, q_ext_r=0, s_ext_l=0, s_ext_r=0;
  const char* q=q_seq+q_alnstart;
  int q_avail=q_max-q_alnstart;
  const char* s=s_seq+s_alnstart;
  int s_avail=s_max-s_alnstart;
  if (penalty<0) penalty=-penalty;
  GXAlnInfo* alninfo=NULL;
  bool freeAlnMem=(gxmem==NULL);
  if (freeAlnMem) {
     gxmem=new CGreedyAlignData(reward, penalty, xdrop);
     reward=gxmem->match_reward;
     penalty=gxmem->mismatch_penalty;
     xdrop=gxmem->x_drop;
     }
   else
	 gxmem->reset();
  int minMatch= trim ? trim->minMatch : 6;
  int MIN_GREEDY_SCORE=minMatch*reward; //minimum score for an alignment to be reported for 0 diffs
  int retscore = 0;
  int numdiffs = 0;
  if (trim!=NULL && trim->type==galn_TrimLeft) {
    //intent: trimming the left side
    if (editscript)
       ed_script_rev=new GXEditScript();

    int numdiffs_l = GXGreedyExtend(s_seq, s_alnstart, q_seq, q_alnstart, true, // xdrop, reward, penalty,
    	            s_ext_l, q_ext_l, *gxmem, ed_script_rev);
    //check this extension here and bail out if it's not a good extension
    if (s_ext_l+(trim->seedlen>>1) < trim->safelen &&
    	q_alnstart+1-q_ext_l>1 &&
    	s_alnstart+1-s_ext_l>trim->l_boundary) {
      #ifdef TRIMDEBUG 
       GMessage(".. 5' greedy alignment rejected (trim_len=%d, trim_safelen=%d, "
          "q_delta=%d, s_delta=%d, trim_l_boundary=%d)\n",
          s_ext_l+(trim->seedlen>>1), trim->safelen, 
          q_alnstart+1-q_ext_l, s_alnstart+1-s_ext_l, trim->l_boundary);
      #endif

      delete ed_script_rev;
      if (freeAlnMem) delete gxmem;
      return NULL;
      }

    if (editscript)
       ed_script_fwd=new GXEditScript();
    int numdiffs_r = GXGreedyExtend(s, s_avail, q, q_avail, false, //xdrop, reward, penalty,
    	                             s_ext_r, q_ext_r, *gxmem, ed_script_fwd);
    numdiffs=numdiffs_r+numdiffs_l;
    //convert num diffs to actual score
    retscore = (q_ext_r + s_ext_r + q_ext_l + s_ext_l)*(reward>>1) - numdiffs*(reward+penalty);
    if (editscript)
       ed_script_rev->Append(ed_script_fwd); //combine the two extensions
    }
  else {
    if (editscript) {
       ed_script_fwd=new GXEditScript();
       }
    int numdiffs_r = GXGreedyExtend(s, s_avail, q, q_avail, false, // xdrop, reward, penalty,
    	                             s_ext_r, q_ext_r, *gxmem, ed_script_fwd);
    //check extension here and bail out if not a good right extension
    //assuming s_max is really at the right end of s_seq
    if (trim!=NULL && trim->type==galn_TrimRight &&
        s_ext_r+(trim->seedlen>>1) < trim->safelen &&
            q_alnstart+q_ext_r<q_max-2 &&
            s_alnstart+s_ext_r<trim->r_boundary) {
      delete ed_script_fwd;
      if (freeAlnMem) delete gxmem;
      return NULL;
      }
    if (editscript)
       ed_script_rev=new GXEditScript();
    int numdiffs_l =  GXGreedyExtend(s_seq, s_alnstart, q_seq, q_alnstart, true, // xdrop, reward, penalty,
    	                          s_ext_l, q_ext_l, *gxmem, ed_script_rev);
    //convert num diffs to actual score
    numdiffs=numdiffs_r+numdiffs_l;
    retscore = (q_ext_r + s_ext_r + q_ext_l + s_ext_l)*(reward>>1) - numdiffs*(reward+penalty);
    if (editscript)
       ed_script_rev->Append(ed_script_fwd); //combine the two extensions
  }
  if (retscore>=MIN_GREEDY_SCORE) {
    alninfo=new GXAlnInfo(q_seq, q_alnstart+1-q_ext_l, q_alnstart+q_ext_r, s_seq, s_alnstart+1-s_ext_l, s_alnstart+s_ext_r);
    int hsp_length = GMIN(q_ext_l+q_ext_r, s_ext_l+s_ext_r);
    alninfo->score=retscore;
    if (gxmem->scaled) alninfo->score >>= 1;
    alninfo->pid = 100 * (1 - ((double) numdiffs) / hsp_length);
#ifdef TRIMDEBUG
    //if (ed_script_rev) {
    //   GMessage("Final Edit script ::: ");
    //   printEditScript(ed_script_rev);
    //   }
#endif
    alninfo->editscript=ed_script_rev;
    alninfo->gapinfo = new CAlnGapInfo(ed_script_rev, alninfo->ql-1, alninfo->sl-1);
    }
  else {
    #ifdef TRIMDEBUG 
      GMessage(".. greedy extension rejected (score=%d)\n", retscore);
    #endif
    //if (freeAlnMem) delete gxmem;
    delete ed_script_rev;
    delete alninfo;
    alninfo=NULL;
    }
  if (freeAlnMem) delete gxmem;
  delete ed_script_fwd;
  return alninfo;
 }

GXAlnInfo* GreedyAlignRegion(const char* q_seq, int q_alnstart, int q_max,
                       const char* s_seq, int s_alnstart, int s_max, CGreedyAlignData* gxmem,
                       CAlnTrim* trim, bool editscript) {
int reward=2;
int penalty=10;
int xdrop=32;
if (gxmem) {
   reward=gxmem->match_reward;
   penalty=gxmem->mismatch_penalty;
   xdrop=gxmem->x_drop;
   }
 return GreedyAlignRegion(q_seq, q_alnstart, q_max, s_seq, s_alnstart, s_max,
     reward, penalty, xdrop, gxmem, trim, editscript);
}

GXAlnInfo* match_adapter(GXSeqData& sd, GAlnTrimType trim_type, int minMatch,
	                             CGreedyAlignData* gxmem, double min_pid) {
  bool editscript=false;
  #ifdef TRIMDEBUG
   editscript=true;
   if (trim_type==galn_TrimLeft) {
	 GMessage("=======> searching left (5') end : %s\n", sd.aseq);
     }
   else if (trim_type==galn_TrimRight) {
     GMessage("=======> searching right(3') end : %s\n", sd.aseq);
     }
   else if (trim_type==galn_TrimEither) {
     GMessage("==========> searching  both ends : %s\n", sd.aseq);
     }
  #endif
  CAlnTrim trimInfo(trim_type, sd.bseq, sd.blen, sd.alen, minMatch, sd.amlen);
  GList<GXSeed> rseeds(true,true,false);
  GXBandSet* alnbands=collectSeeds(rseeds, sd, trim_type);
  GList<GXSeed> anchor_seeds(cmpSeedDiag, NULL, true); //stores unique seeds per diagonal
  //did we find a shortcut?
  if (alnbands->qmatch) {
    #ifdef TRIMDEBUG
     GMessage("::: Found a quick long match at %d, len %d\n",
          alnbands->qmatch->b_ofs, alnbands->qmatch->len);
    #endif
    anchor_seeds.Add(alnbands->qmatch);
    }
  else {
    int max_top_bands=5;
    int top_band_count=0;
    for (int b=0;b<alnbands->Count();b++) {
       if (alnbands->Get(b)->score<6) break;
       //#ifdef TRIMDEBUG
       //GMessage("\tBand %d score: %d\n", b, alnbands->Get(b)->score);
       //#endif
       top_band_count++;
       GXBand& band=*(alnbands->Get(b));
       band.seeds.setSorted(cmpSeedScore);
       anchor_seeds.Add(band.seeds.First());
       //band.tested=true;
       if (anchor_seeds.Count()>2 || top_band_count>max_top_bands) break;
       }
    //#ifdef TRIMDEBUG
    //GMessage("::: Collected %d anchor seeds.\n",anchor_seeds.Count());
    //#endif
    }
  GList<GXAlnInfo> galns(true,true,false);
  for (int i=0;i<anchor_seeds.Count();i++) {
    GXSeed& aseed=*anchor_seeds[i];
    int a1=aseed.a_ofs+(aseed.len>>1)+1;
    int a2=aseed.b_ofs+(aseed.len>>1)+1;
    trimInfo.seedlen=aseed.len;
#ifdef TRIMDEBUG
    GMessage("\t::: align from seed (%d, %d) of len %d.\n",aseed.a_ofs, aseed.b_ofs,
    	   aseed.len);
#endif
    GXAlnInfo* alninfo=GreedyAlignRegion(sd.aseq, a1, sd.alen,
                            sd.bseq, a2, sd.blen, gxmem, &trimInfo, editscript);

#ifdef TRIMDEBUG
     if (alninfo) {
        GMessage("\t::: aln pid=%4.2f (vs. min_pid=%.2f)\n", alninfo->pid, min_pid);
        alninfo->gapinfo->printAlignment(stderr, sd.aseq, sd.alen, sd.bseq, sd.blen);
        }
     else
       GMessage("\t::: GreedyAlignRegion failed.\n");
#endif

    if (alninfo && alninfo->pid>=min_pid && trimInfo.validate(alninfo))
             galns.AddIfNew(alninfo, true);
        else delete alninfo;
    }

  if (galns.Count()==0) {
	 //last resort: look for weaker terminal seeds
	  GPVec<GXSeed> tmatches(2,false);
	  if (trim_type!=galn_TrimRight) {
		 if (alnbands->tmatch_l)
		    tmatches.Add(alnbands->tmatch_l);
	     }
	  if (trim_type!=galn_TrimLeft) {
		 if (alnbands->tmatch_r)
		    tmatches.Add(alnbands->tmatch_r);
	     }
	  for (int i=0;i<tmatches.Count();i++) {
		GXSeed& aseed=*tmatches[i];
		int halfseed=aseed.len>>1;
		int a1=aseed.a_ofs+halfseed+1;
		int a2=aseed.b_ofs+halfseed+1;
		trimInfo.seedlen=aseed.len;
#ifdef TRIMDEBUG
    GMessage("\t::: align from terminal seed (%d, %d)of len %d.\n",aseed.a_ofs, aseed.b_ofs,
    	   aseed.len);
#endif
        GXAlnInfo* alninfo=GreedyAlignRegion(sd.aseq, a1, sd.alen,
                                sd.bseq, a2, sd.blen, gxmem, &trimInfo, editscript);
        if (alninfo && alninfo->pid>=min_pid && trimInfo.validate(alninfo))
                 galns.AddIfNew(alninfo, true);
             else delete alninfo;
        }//for each terminal seed
      }
  //---- found all alignments
  delete alnbands;
  
  #ifdef TRIMDEBUG
  //print all valid alignments found
  for (int i=0;i<galns.Count();i++) {
    GXAlnInfo* alninfo=galns[i];
    GMessage("a(%d..%d) align to b(%d..%d), score=%d, pid=%4.2f\n", alninfo->ql, alninfo->qr,
                         alninfo->sl, alninfo->sr, alninfo->score, alninfo->pid);
    if (alninfo->gapinfo!=NULL) {
      GMessage("Alignment:\n");
      alninfo->gapinfo->printAlignment(stderr, sd.aseq, sd.alen, sd.bseq, sd.blen);
      }
    }
  #endif
  
  if (galns.Count()) {
    GXAlnInfo* bestaln=galns.Shift();
    #ifdef TRIMDEBUG
      GMessage("Best alignment: a(%d..%d) align to b(%d..%d), score=%d, pid=%4.2f\n", bestaln->ql, bestaln->qr,
          bestaln->sl, bestaln->sr, bestaln->score, bestaln->pid);
      if (bestaln->gapinfo!=NULL) {
        bestaln->gapinfo->printAlignment(stderr, sd.aseq, sd.alen, sd.bseq, sd.blen);
        }
    #endif
    return bestaln;
    }
  else return NULL;
}
