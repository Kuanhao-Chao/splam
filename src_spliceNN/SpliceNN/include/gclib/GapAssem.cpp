#include "GapAssem.h"

const unsigned char GA_flag_IS_REF=0;
const unsigned char GA_flag_HAS_PARENT=1;
const unsigned char GA_flag_BAD_ALIGN=7;
const unsigned char GA_flag_PREPPED=5;

//bool GASeq::debug=false;
bool MSAColumns::removeConsGaps = true;
bool MSAColumns::refineClipping = true;
//unsigned int GSeqAlign::counter = 0;

int qsortnuc(const void* p1, const void* p2) {
	GAlnColumn::NucCount* c1 = (GAlnColumn::NucCount*) p1;
	GAlnColumn::NucCount* c2 = (GAlnColumn::NucCount*) p2;
	return (c1->count > c2->count) ? -1 : ((c1->count < c2->count) ? 1 : 0);
}

int compareOrdnum(void* p1, void* p2) {
	int v1 = ((GSeqAlign*) p1)->ordnum;
	int v2 = ((GSeqAlign*) p2)->ordnum;
	return (v1 < v2) ? -1 : ((v1 > v2) ? 1 : 0);
}

int compareCounts(void* p1, void* p2) {
	int v1 = ((GSeqAlign*) p1)->Count();
	int v2 = ((GSeqAlign*) p2)->Count();
	return (v1 > v2) ? -1 : ((v1 < v2) ? 1 : 0);
}

GASeq::GASeq(const char* sname, const char* sdescr, const char* sseq, int slen, int soffset) :
		FastaSeq(sname, sdescr, sseq, slen),
		  numgaps(0), ofs(NULL), delops(false, true, false),flags(0),msa(NULL),
		  msaidx(-1), seqlen(slen), offset(soffset),
		  ng_ofs(soffset),revcompl(0), ext5(0), ext3(0),
		  clp5(0), clp3(0) {
	seqlen = len; //FastaSeq constructor settles it
	if (seqlen>0) {
      GCALLOC(ofs, seqlen * sizeof(short));
#ifdef ALIGN_COVERAGE_DATA
	GCALLOC(cov,seqlen*sizeof(int));
#endif
	}
}

GASeq::GASeq(const char* sname, int soffset, int slen, int sclipL, int sclipR,
    char rev) : FastaSeq(sname),
  		  numgaps(0), ofs(NULL), delops(false, true, false),flags(0),msa(NULL),
  		  msaidx(-1), seqlen(slen), offset(soffset),
  		  ng_ofs(soffset),revcompl(rev), ext5(0), ext3(0),
  		  clp5(sclipL), clp3(sclipR) {
	if (seqlen>0) {
	  GCALLOC(ofs, seqlen * sizeof(short));
#ifdef ALIGN_COVERAGE_DATA
	  GCALLOC(cov,seqlen*sizeof(int));
#endif
	}
}

GASeq::GASeq(GASeq& aseq): FastaSeq(aseq.id, aseq.descr, aseq.seq, aseq.len),
		  numgaps(0), ofs(NULL), delops(false, true, false),flags(0),msa(NULL),
		  msaidx(-1), seqlen(aseq.len), offset(0),
		  ng_ofs(0),revcompl(0), ext5(0), ext3(0),
		  clp5(0), clp3(0) {
	if (seqlen>0) {
      GCALLOC(ofs, seqlen * sizeof(short));
#ifdef ALIGN_COVERAGE_DATA
	  GCALLOC(cov,seqlen*sizeof(int));
#endif
	}
}

GASeq::GASeq(FastaSeq& faseq, bool takeover):FastaSeq(faseq, takeover),
		  numgaps(0), ofs(NULL), delops(false, true, false),flags(0),msa(NULL),
		  msaidx(-1), seqlen(len), offset(0),
		  ng_ofs(0),revcompl(0), ext5(0), ext3(0),
		  clp5(0), clp3(0) {
	if (seqlen>0) {
      GCALLOC(ofs, seqlen * sizeof(short));
#ifdef ALIGN_COVERAGE_DATA
	  GCALLOC(cov,seqlen*sizeof(int));
#endif
	}
}

GASeq::~GASeq() {
	GFREE(ofs);
#ifdef ALIGN_COVERAGE_DATA
	GFREE(cov);
#endif
}
/*
void GASeq::loadProcessing() {
	//process all delops
	for (int i = 0; i < delops->Count(); i++) {
		SeqDelOp& delop = *delops->Get(i);
		int pos = delop.revcompl ? len - delop.pos - 1 : delop.pos;
		removeBase(pos);
	}
	if (revcompl == 1)
		reverseComplement();
}
*/

void GASeq::finalize() {
	 if (this->len==0)  GError("Error: sequence for %s not loaded!\n",this->getId());
	 if (!this->hasFlag(GA_flag_PREPPED)) this->prepSeq();
}

void GASeq::prepSeq() {
	//should only be called once (use hasFlag() before calling)
	//apply all deletions to the sequence
	for (int i = 0; i < delops.Count(); i++) {
		SeqDelOp& delop = *delops.Get(i);
		int pos = delop.revcompl ? len - delop.pos - 1 : delop.pos;
		removeBase(pos);
	}
	if (revcompl == 1)
		reverseComplement();

	setFlag(GA_flag_PREPPED);
}

//set the gap length in this position
void GASeq::setGap(int pos, short gaplen) {
	if (pos < 0 || pos >= seqlen)
		GError("Error: invalid gap position (%d) given for sequence %s\n", pos + 1,
		    id);
	numgaps -= ofs[pos];
	ofs[pos] = gaplen;
	numgaps += gaplen;
}

//add to the existing gap length in this position
void GASeq::addGap(int pos, short gapadd) {
	if (pos < 0 || pos >= seqlen)
		GError("Error: invalid gap position (%d) given for sequence %s\n", pos + 1,
		    id);
	numgaps += gapadd;
	ofs[pos] += gapadd;
}

void GASeq::removeBase(int pos) {
	if (pos < 0 || pos >= seqlen)
		GError("Error: invalid gap position (%d) given for sequence %s\n", pos + 1,
		    id);
//if there is a gap at that position, remove the gap
//otherwise, remove the actual nucleotide!
//if (ofs[pos]>0) {
	ofs[pos]--;
	numgaps--;
//  return;
//  }
	/* if it's end base or within clipping --
	 //   don't remove, just adjust clipping
	 if (revcompl!=0) { //reversed in this alignment
	 if (pos<=clp3) {
	 if (pos==clp3) clp3++;
	 offset--;
	 return;
	 }
	 if (pos>=seqlen-clp5-1) {
	 if (pos==seqlen-clp5-1) clp5++;
	 return;
	 }
	 }
	 else {//forward
	 if (pos<=clp5) {
	 if (pos==clp5) clp5++;
	 offset--;
	 return;
	 }
	 if (pos>=seqlen-clp3-1) {
	 if (pos==seqlen-clp3-1) clp3++;
	 return;
	 }
	 }
	 */
//-- couldn't just this do it ?
	/*
	 ofs[pos]--;
	 numgaps--;
	 return;
	 */
	/*
	 //===========================================
	 // worst case: modify/rebuild the whole sequence..
	 //short* newofs;
	 seqlen--;
	 memmove(ofs+pos, ofs+pos+1, seqlen-pos);
	 //do the same for the actual sequence, if loaded
	 if (len>0) {//sequence loaded
	 len--;
	 memmove(seq+pos, seq+pos+1, len-pos);
	 endSeq();
	 }
	 else { //push this sequence editing information on a stack for later processing
	 delops->Add(new SeqDelOp(pos, revcompl!=0));
	 }
	 */
}

void GASeq::refineClipping(GDynArray<char>& cons, int cpos, bool skipDels) {
	//check if endings match consensus..
	//adjust clipping as appropriate
	//int clipL, clipR;
	if (clp3 == 0 && clp5 == 0)
		return;
	int& clipL = (revcompl != 0) ? clp3 : clp5;
	int& clipR = (revcompl != 0) ? clp5 : clp3;
	//build the gapped sequence string in memory
	char* gseq;
	int glen = seqlen + numgaps;
	int allocsize = glen;
	int gclipR = clipR;
	int gclipL = clipL;
	if (skipDels) {
		for (int i = 1; i <= clipR; i++) {
			if (ofs[seqlen - i] < 0)
				allocsize++;
			else
				gclipR += ofs[seqlen - i];
		}
		for (int i = 0; i < clipL; i++) {
			if (ofs[i] < 0)
				allocsize++;
			else
				gclipL += ofs[i];
		}
	} else {
		for (int i = 1; i <= clipR; i++)
			gclipR += ofs[seqlen - i];
		for (int i = 0; i < clipL; i++)
			gclipL += ofs[i];
	}
	int* gxpos; //mapping of positions from gseq to seq
	GMALLOC(gxpos, allocsize * sizeof(int));
	GMALLOC(gseq, allocsize + 1);
	gseq[allocsize] = 0;
	int gseqpos = 0;
	for (int i = 0; i < seqlen; i++) {
		//bool notClip=(i>=clipL && i<seqlen-clipR);
		if (ofs[i] < 0) {
			if (!skipDels)
				continue; //always skip gaps
			if (i >= clipL && i < seqlen - clipR) //in non-clipped region
				continue; //skip gaps
			else
				glen++;
		}
		for (int j = 0; j < ofs[i]; j++) {
			gseq[gseqpos] = '*';
			gseqpos++;
		}
		gseq[gseqpos] = seq[i];
		gxpos[gseqpos] = i;
		gseqpos++;
	}
	gseq[allocsize] = 0;
	if (glen != allocsize)
		GError(
		    "Length mismatch (allocsize %d vs. glen %d) while refineClipping for seq %s !\n",
		    allocsize, glen, id);
	//adjust end clipping, using a simple X-drop algorithm
	// match_reward=1, mismatch_penalty=-3
#define XDROP -16
#define MATCH_SC 1
#define MISMATCH_SC -3
	if (clipR > 0) {
		//-------------- clipR -------------------------
		// actual span of clipped regions in the alignment
		// could be larger then clipR/clipL due to gaps propagated
		// WITHIN the clipped regions

		//******** right end adjust
		// cp = last "matching" position on consensus
		int cp = cpos + glen - gclipR - 1;
		// sp = corresponding last match position on read
		int sp = glen - gclipR - 1;
		//we could be on a mismatch or a gap
		// so, first go backward to find the first match
		while (gseq[sp] != cons[cp] || gseq[sp] == '*') {
			if (gseq[sp] != '*')
				clipR++;
			sp--;
			cp--;
			if (sp < gclipL) {
				GMessage(
				    "Warning: reached clipL trying to find an initial match on %s!\n",
				    id);
				GFREE(gseq);
				return;     //break
			}
		}

		//now go forward for as much as we can, using the dropoff test

		int score = MATCH_SC;
		int maxscore = MATCH_SC;
		int startpos = sp;
		int bestpos = sp; //new Right clipping position for maxscore
		while (score > XDROP && ++cp < (int)cons.Count() && ++sp < glen) {
			if (gseq[sp] == cons[cp]) {
				if (gseq[sp] != '*') { //real match
					score += MATCH_SC;
					//clpinc++;
					if (score > maxscore) {
						bestpos = sp;
						maxscore = score;
					} //better score than before
				} //real match
			} //match
			else { //mismatch
				if (gseq[sp] != '*') {
					score += MISMATCH_SC; //clipinc++;
				}
			}
		} //while XDROP
		if (bestpos > startpos)
			clipR = seqlen - gxpos[bestpos] - 1;
	} //<------- was clipR

	// ******** left end adjust
	if (clipL > 0) {
		// cp = last "matching" position on consensus
		int cp = cpos + gclipL;
		// sp = corresponding last match position on read
		int sp = gclipL;
		// we could be on a mismatch or a gap
		// so, first go backward to find the first match
		while (gseq[sp] != cons[cp] || gseq[sp] == '*') {
			if (gseq[sp] != '*')
				clipL++;
			sp++;
			cp++;
			if (sp >= glen - gclipR) {
				GMessage(
				    "Warning: reached clipR trying to find an initial match on %s!\n",
				    id);
				GFREE(gseq);
				return;     //break
			}
		}
		//-- now go backward for as much as we can, using the dropoff test
		int score = MATCH_SC;
		int maxscore = MATCH_SC;
		int startpos = sp;
		int bestpos = sp;
		while (score > XDROP && --cp >= 0 && --sp >= 0) {
			if (gseq[sp] == cons[cp]) {
				if (gseq[sp] != '*') {     //real match
					score += MATCH_SC;
					if (score > maxscore) {
						bestpos = sp;
						maxscore = score;
					}     //better score than before
				}     //real match
			}     //match
			else {     //mismatch
				if (gseq[sp] != '*') {
					score += MISMATCH_SC;
				}
			}
		}     //while XDROP
		if (bestpos < startpos)
			clipL = gxpos[bestpos];
	}     //is clipL
	GFREE(gseq);
	GFREE(gxpos);
}

void GASeq::reverseGaps() {
	//--when reading mgblast alignments and gap info
	//the gap positions are reversed starting and shifted by 1
	//because the first ofs is always 0
	int l = 1;
	int r = seqlen - 1;
	while (l < r) {
		short c = ofs[l];
		ofs[l] = ofs[r];
		ofs[r] = c;
		l++;
		r--;
	}
}

void GASeq::revComplement(int alignlen) {
	if (alignlen > 0) {  //sequence is in an alignment
		offset = alignlen - endOffset();
		if (msa != NULL) {
			ng_ofs = msa->ng_len - endNgOffset();
			if (msa->minoffset > offset)
				msa->minoffset = offset;
			if (msa->ng_minofs > ng_ofs)
				msa->ng_minofs = ng_ofs;

		}
	}
	revcompl = !revcompl;
	if (len == seqlen)  //sequence is loaded, complement it
		reverseComplement();
	reverseGaps();
	//-- also reverse the coverage array:
#ifdef ALIGN_COVERAGE_DATA
	int l=0;int r=seqlen-1;
	while (l<r) {
		int c=cov[l];
		cov[l]=cov[r];
		cov[r]=c;
		l++;r--;
	}
#endif
}

#ifdef ALIGN_COVERAGE_DATA
void GASeq::addCoverage(GASeq* s) {
	// add coverage information from the same sequence involved
	// in another pairwise alignment
	if (seqlen!=s->seqlen)
	GError("GSeqAlign Error: invalid addCoverage %s(len %d) vs %s(len %d)\n",
			name(),seqlen, s->name(), s->seqlen);
	if (s->revcompl!=revcompl) {
		for (int i=0;i<seqlen;i++)
		cov[i]+=s->cov[seqlen-i-1];
	}
	else
	for (int i=0;i<seqlen;i++)
	cov[i]+=s->cov[i];

}
#endif

void GASeq::printGappedSeq(FILE* f, int baseoffs) {
	// for now, simple console printing of short sequences -- for testing
	if (len == 0 || len != seqlen)
		GError(
		    "GASeq print Error: invalid sequence data '%s' (len=%d, seqlen=%d)\n",
		    id, len, seqlen);
	int i;
	int clipL, clipR;
	if (revcompl != 0) {
		clipL = clp3;
		clipR = clp5;
	} else {
		clipL = clp5;
		clipR = clp3;
	}
	for (i = 0; i < (offset - baseoffs); i++)
		fprintf(f, " ");
	for (i = 0; i < seqlen; i++) {
		if (ofs[i] < 0)
			continue; //deleted base
		for (int j = 0; j < ofs[i]; j++)
			fprintf(f, "-");
		char c = seq[i];
		if (i < clipL || i >= seqlen - clipR)
			c = (char) tolower(c);
		fprintf(f, "%c", c);
	} //for each base
	fprintf(f, "\n");
}

void GASeq::printGappedFasta(FILE* f) {
	if (len == 0 || len != seqlen)
		GError(
		    "GASeq print Error: invalid sequence data '%s' (len=%d, seqlen=%d)\n",
		    id, len, seqlen);
	int i;
	/*  //FIXME TESTME - original mblaor had this uncommented!
	int clipL, clipR;
	if (revcompl != 0) {
		clipL = clp3;
		clipR = clp5;
	} else {
		clipL = clp5;
		clipR = clp3;
	}
	*/
	int printed = 0;
	for (i = 0; i < seqlen; i++) {
		if (ofs[i] < 0)
			continue; //deleted base
		for (int j = 0; j < ofs[i]; j++) {
			fprintf(f, "*");
			printed++;
			if (printed == 60) {
				fprintf(f, "\n");
				printed = 0;
			}
		}
		char c = seq[i];
		printed++;
		if (printed == 60) {
			fprintf(f, "%c\n", c);
			printed = 0;
		} else
			fprintf(f, "%c", c);
	} //for each base
	if (printed < 60)
		fprintf(f, "\n");
}

void GASeq::printMFasta(FILE* f, int llen) {
	if (len == 0 || len != seqlen)
		GError("GASeq print Error: invalid sequence data '%s' (len=%d, seqlen=%d)\n",
		    id, len, seqlen);
	if (this->descrlen>0) fprintf(f, ">%s %s\n", id, descr);
	else
	  fprintf(f, ">%s\n", id);
	int i;
	int printed = 0;
	for (i=0;i<offset;i++) {
		fprintf(f, "-");
		printed++;
		if (printed == llen) {
			fprintf(f, "\n");
			printed = 0;
		}
	}
	for (i = 0; i < seqlen; i++) {
		if (ofs[i] < 0)
			continue; //deleted base
		for (int j = 0; j < ofs[i]; j++) {
			fprintf(f, "-");
			printed++;
			if (printed == llen) {
				fprintf(f, "\n");
				printed = 0;
			}
		}
		char c = seq[i];
		printed++;
		if (printed == llen) {
			fprintf(f, "%c\n", c);
			printed = 0;
		} else
			fprintf(f, "%c", c);
	} //for each base
	if (printed < llen)
		fprintf(f, "\n");
}

int GASeq::removeClipGaps() { //remove gaps within clipped regions
	//offset is also corrected appropriately!
	int clipL;
	int clipR;
	if (revcompl != 0) {
		clipL = clp3;
		clipR = clp5;
	} else {
		clipL = clp5;
		clipR = clp3;
	}
	int delgapsL = 0;
	int delgapsR = 0;
	for (int i = 0; i < seqlen; i++) {
		if (i <= clipL) { // within left clipping
			delgapsL += ofs[i];
			ofs[i] = 0;
			continue;
		} //within left clipping
		if (i >= seqlen - clipR) { //with right clipping
			delgapsR += ofs[i];
			ofs[i] = 0;
		}
	} //for
	offset += delgapsL;
	numgaps -= (delgapsL + delgapsR);
	return delgapsL + delgapsR;
}

void GASeq::toMSA(MSAColumns& msacols, int nucValue) {
	if (len == 0 || len != seqlen)
		GError(
		    "GASeq::toMSA Error: invalid sequence data '%s' (len=%d, seqlen=%d)\n",
		    id, len, seqlen);
	int i;
	int clipL, clipR;
	if (revcompl != 0) {
		clipL = clp3;
		clipR = clp5;
	} else {
		clipL = clp5;
		clipR = clp3;
	}
	int mincol = INT_MAX;
	int maxcol = 0;
	//for (i=0;i<(offset-msa.baseoffset);i++) col++;
	int col = offset - msa->minoffset;
	for (i = 0; i < seqlen; i++) {
		bool clipped;
		if (i < clipL || i >= seqlen - clipR)
			clipped = true;
		else {
			clipped = false;
			if (mincol == INT_MAX)
				mincol = col;
		}
		for (int j = 0; j < ofs[i]; j++) {
			//storing gap
			if (!clipped)
				msacols[col].addGap(nucValue);
			col++;
		}
		msacols[col].addNuc(this, i, clipped, nucValue);
		if (!clipped)
			maxcol = col;
		col++;
	} //for each base
	//update msa's min-max
	msacols.updateMinMax(mincol, maxcol);
}

//=================================== GSeqAlign ===============================

// -- creation from a pairwise alignment
// l1-r1 = coordinates of alignment region on s1
// l2-r2 = coordinates of alignment region on s2
// coordinates MUST be 0-based
#ifdef ALIGN_COVERAGE_DATA
GSeqAlign::GSeqAlign(GASeq* s1, int l1, int r1,
		GASeq* s2, int l2, int r2)
//:GList<GASeq>(true,true,false) {
:GList<GASeq>(false,true,false) {
#else
GSeqAlign::GSeqAlign(GASeq* s1, GASeq* s2) :
		GList<GASeq>(false, true, false), length(0), minoffset(0),
  		refinedMSA(false), msacolumns(NULL), ordnum(0),
  		ng_len(0),ng_minofs(0), badseqs(0), consensus(512), consensus_bq(512) {
#endif
	s1->msa = this;
	s2->msa = this;
	//the offset for at least one sequence is 0
	this->Add(s1);
	this->Add(s2);
	minoffset = GMIN(s1->offset, s2->offset);
	ng_minofs = minoffset; //no gaps in the clipped regions for now
	length = GMAX(s1->endOffset(), s2->endOffset());
	length -= minoffset;
	ng_len = GMAX(s1->endNgOffset(), s2->endNgOffset());
	ng_len -= ng_minofs;
	//-- according to the alignment, update the coverage for each sequence
	//-- overlaps are granted +1 bonus
#ifdef ALIGN_COVERAGE_DATA
	for (int i=l1;i<r1;i++) s1->cov[i]++;
	for (int i=l2;i<r2;i++) s2->cov[i]++;
	//-- mismatch regions at the left end
	int msml=(l2>l1) ? l1:l2;
	for (int i=1;i<=msml;i++) {
		s1->cov[l1-msml]--;
		s2->cov[l2-msml]--;
	}
	//-- mismatch regions at the right end
	int cr1=s1->seqlen-r1-1;
	int cr2=s2->seqlen-r2-1;
	int msmr=(cr2>cr1) ? cr1:cr2;
	for (int i=1;i<=msmr;i++) {
		s1->cov[r1+msmr]--;
		s2->cov[r2+msmr]--;
	}
#endif
}

//merge other alignment omsa into this msa
//seq->id MUST be the same with oseq->id
bool GSeqAlign::addAlign(GASeq* seq, GSeqAlign* omsa, GASeq* oseq) {
	//error checking -- could be disabled to speed it up a bit
	if (seq->seqlen != oseq->seqlen)
		GError("GSeqAlign Error: invalid merge %s(len %d) vs %s(len %d)\n",
		    seq->getName(), seq->seqlen, oseq->getName(), oseq->seqlen);
	// for this merge to work, the shared sequence MUST have
	// the same orientation in both MSAs
	if (seq->revcompl != oseq->revcompl)
		omsa->revComplement(); //reverse-complement all sequences in omsa
#ifdef ALIGN_COVERAGE_DATA
	//add coverage values:
	seq->addCoverage(oseq);
#endif
	//--- now propagate gaps as appropriate
	for (int i = 0; i < seq->seqlen; i++) {
		int d = seq->gap(i) - oseq->gap(i);
		if (d > 0) { //extra gap in seq
		             //propagate into all omsa
			omsa->injectGap(oseq, i, d);
			continue;
		} //extra gap in seq
		if (d < 0) { //extra gap in oseq
		             //propagate into all msa
			injectGap(seq, i, -d);
			continue;
		}              //extra gap in oseq
	}              //--for each base position
	//--now add the sequences from omsa to this MSA
	for (int i = 0; i < omsa->Count(); i++) {
		GASeq* s = omsa->Get(i);
		if (s == oseq)
			continue;
		//adjust offset -- which can be extended by gaps in seq BEFORE s->offset
		//the offsets had been adjusted already (by injectGap() method)
		// to account for propagated gaps in both MSAs!
		//--add this sequence
		//cluster minoffset and length will be updated too!
		addSeq(s, seq->offset + s->offset - oseq->offset,
		    seq->ng_ofs + s->ng_ofs - oseq->ng_ofs);

	}
	omsa->setFreeItem(false);
	//delete omsa; //we no longer need this alignment
	delete oseq; //also deletes oseq
	return true;
}

//just to automatically set the offset, msa,
//and to update the MSA length if needed
void GSeqAlign::addSeq(GASeq* s, int soffs, int ngofs) {
	s->offset = soffs;
	s->ng_ofs = ngofs;
	s->msa = this;
	this->Add(s);
	//keep track of minimum offset
	//this also adjusts length!
	if (soffs < minoffset) {
		length += minoffset - soffs;
		minoffset = soffs;
	}
	if (ngofs < ng_minofs) {
		ng_len += ng_minofs - ngofs;
		ng_minofs = ngofs;
	}

	//adjust length of alignment if very long sequences were added
	if (s->endOffset() - minoffset > length)
		length = s->endOffset() - minoffset;
	if (s->endNgOffset() - ng_minofs > ng_len)
		ng_len = s->endNgOffset() - ng_minofs;

}

//propagate a gap in a sequence into the whole alignment containing it
// offsets of all seqs after the gap MUST be adjusted too!
void GSeqAlign::injectGap(GASeq* seq, int pos, int xgap) {
	//find the actual alignment position of this pos in the layout
	int alpos = seq->offset + pos;
	for (int i = 0; i <= pos; i++)
		alpos += seq->gap(i);
	//now alpos = the exact offset of seq[pos] in this MSA
	for (int i = 0; i < Count(); i++) {
		GASeq* s = Get(i);
		int spos = 0; // finding out position of gap in seq s
		if (s == seq)
			spos = pos;
		else {
			//walk to lpos on sequence s
			int salpos = s->offset;
			if (salpos >= alpos) {
				//s->offset is AFTER this gap, so only the offset is affected
				s->offset += xgap;
				continue;
			}
			while (spos < s->seqlen) {
				salpos += 1 + s->gap(spos);
				if (salpos > alpos)
					break;
				spos++;
			}
			if (spos >= s->seqlen) //spos is AFTER the end of sequence s
				continue; // s not affected
			//--it is a valid position for this sequence
			//--TO DO: clipping? first/last positions?
		}
		s->addGap(spos, xgap);
	}      //for each sequense in MSA
	length += xgap;
}

void GSeqAlign::removeColumn(int column) {
	int alpos = column + minoffset;
	for (int i = 0; i < Count(); i++) {
		GASeq* s = Get(i);
		int spos = 0; // finding out position of this base in seq s
		int salpos = s->offset;
		if (salpos >= alpos) {
			//s->offset is AFTER this gap, so only the offset is affected
			s->offset--; //deletion of 1
			continue;
		}
		while (spos < s->seqlen) {
			salpos += 1 + s->gap(spos);
			if (salpos > alpos)
				break;
			spos++;
		}
		if (spos >= s->seqlen) //spos is AFTER the end of sequence s
			continue; // s not affected
		//--now spos is a valid position for this sequence
		//    }
		s->removeBase(spos);
	}  //for each sequence in MSA
	length--;
}

void GSeqAlign::removeBase(GASeq* seq, int pos) {
	//find the actual alignment position of this pos in the layout
	int alpos = seq->offset + pos;
	for (int i = 0; i <= pos; i++)
		alpos += seq->gap(i);
	//now alpos = the exact offset of seq[pos] in this MSA
	for (int i = 0; i < Count(); i++) {
		GASeq* s = Get(i);
		int spos = 0; // finding out position of this base in seq s
		if (s == seq)
			spos = pos;
		else { //walk to lpos on sequence s
			int salpos = s->offset;
			if (salpos >= alpos) {
				//s->offset is AFTER this gap, so only the offset is affected
				s->offset--; //deletion of 1
				continue;
			}
			while (spos < s->seqlen) {
				salpos += 1 + s->gap(spos);
				if (salpos > alpos)
					break;
				spos++;
			}
			if (spos >= s->seqlen) //spos is AFTER the end of sequence s
				continue; // s not affected
			//--now spos is a valid position for this sequence
		}
		s->removeBase(spos);
	}      //for each sequence in MSA
	length--;
}

void GSeqAlign::applyClipping(AlnClipOps& clipops) {
	for (int i = 0; i < clipops.Count(); i++) {
		SeqClipOp& cop = *clipops.Get(i);
		if (cop.clp[0] >= 0)
			cop.seq->clp5 = cop.clp[0];
		if (cop.clp[1] >= 0)
			cop.seq->clp3 = cop.clp[1];
	}
}
bool GSeqAlign::evalClipping(GASeq* seq, int c5, int c3, float clipmax,
    AlnClipOps& clipops) {
	//propagate trimming of a read to the rest of this container MSA
	//-- returns false if any of the reads in this MSA are clipped too much!
	//GList<SeqClipOp> clipops(false,true,false);
	if (c5 >= 0) {
		//the position of the first/last non-clipped letter
		int pos = (seq->revcompl != 0) ? seq->seqlen - c5 - 1 : c5;
		//find the actual alignment position of this pos in the layout
		int alpos = seq->offset + pos;
		for (int i = 0; i <= pos; i++)
			alpos += seq->gap(i);
		//alpos = the position of seq[pos] in this MSA
		for (int i = 0; i < Count(); i++) {
			GASeq* s = Get(i);
			if (s == seq) {
				if (!clipops.add5(s, c5, clipmax))
					return false;
				continue;
			}
			int spos = 0; // finding out position in seq s
			//walk to lpos on sequence s
			//salpos is going to be the position of seq[pos] in s
			int salpos = s->offset;
			if (salpos >= alpos) {
				//-- s starts AFTER this alpos position
				if (seq->revcompl != 0) { //clipping right side
					// which means ALL of s is to the right => clipped entirely!
					// !!! TODO:
					return false;
					//clipops.Add(new SeqClipOp(s, s->seqlen));
				}
				continue;
			}
			while (spos < s->seqlen) {
				salpos += 1 + s->gap(spos);
				if (salpos > alpos)
					break;
				spos++;
			}
			if (spos >= s->seqlen) {
				//s ends BEFORE this alpos position
				if (seq->revcompl == 0) { //clipping left side
					// which means ALL of s is to the left => clipped entirely!
					return false;
				}
				continue; // s not affected
			}
			//--it is a valid position for this sequence
			//now spos is in the corresponding position of pos
			//trim s here
			if (seq->revcompl != 0) { //trim right side in this msa
				if (s->revcompl != 0) {
					if (!clipops.add5(s, s->seqlen - spos - 1, clipmax))
						return false;
					/*if (s->clp5<newclp) {
					 if (badClipping(s,newclp,s->clp3,clipmax)) return false;
					 clipops.Add(new SeqClipOp(s,newclp));
					 }*/
				} else {
					if (!clipops.add3(s, s->seqlen - spos - 1, clipmax))
						return false;
					/*if (s->clp3<newclp) {
					 if (badClipping(s,s->clp5,newclp,clipmax)) return false;
					 clipops.Add(new SeqClipOp(s,-1,newclp));
					 }*/
				}
			} else { //trim left side in this msa
				if (s->revcompl != 0) {
					if (!clipops.add3(s, spos, clipmax))
						return false;
					/*if (s->clp3<spos) {
					 if (badClipping(s,s->clp5,spos,clipmax)) return false;
					 clipops.Add(new SeqClipOp(s,-1,spos));
					 }*/
				} else {
					if (!clipops.add5(s, spos, clipmax))
						return false;
					/*if (s->clp5<spos) {
					 if (badClipping(s,spos,s->clp3,clipmax)) return false;
					 clipops.Add(new SeqClipOp(s,spos));
					 }*/
				}
			}
		} //for each sequense in MSA
	} // 5' clipping case
//---------------
	if (c3 >= 0) {
		//the position of the first/last non-clipped letter
		int pos = (seq->revcompl != 0) ? c3 : seq->seqlen - c3 - 1;
		//find the actual alignment position of this pos in the layout
		int alpos = seq->offset + pos;
		for (int i = 0; i <= pos; i++)
			alpos += seq->gap(i);
		//now alpos = the exact offset of seq[pos] in this MSA
		for (int i = 0; i < Count(); i++) {
			GASeq* s = Get(i);
			if (s == seq) {
				if (!clipops.add3(s, c3, clipmax))
					return false;
				/*if (s->clp3<c3) {
				 if (badClipping(s,s->clp5,c3,clipmax)) return false;
				 clipops.Add(new SeqClipOp(s,-1,c3));
				 }*/
				continue;
			}
			int spos = 0; // finding out position in seq s
			//walk to lpos on sequence s
			int salpos = s->offset;
			if (salpos >= alpos) {
				//-- s starts AFTER this alpos position
				if (seq->revcompl == 0) { //clipping right side
					// which means ALL of s is to the right => clipped entirely!
					return false;
					//clipops.Add(new SeqClipOp(s, s->seqlen));
				}
				continue;
			}
			while (spos < s->seqlen) {
				salpos += 1 + s->gap(spos);
				if (salpos > alpos)
					break;
				spos++;
			}
			if (spos >= s->seqlen) {
				//s ends BEFORE this alpos position
				if (seq->revcompl != 0) { //clipping left side
					// which means ALL of s is to the left => clipped entirely!
					return false;
					//clipops.Add(new SeqClipOp(s, s->seqlen));
				}
				continue; // s not affected
			}
			//--it is a valid position for this sequence
			//now spos is in the corresponding position of pos
			//trim s here
			if (seq->revcompl != 0) { //trim left side in this msa
				if (s->revcompl != 0) {
					if (!clipops.add3(s, spos, clipmax))
						return false;
					/*if (s->clp3<spos) {
					 if (badClipping(s,s->clp5,spos,clipmax)) return false;
					 clipops.Add(new SeqClipOp(s,-1, spos));
					 }*/
				} else {
					if (!clipops.add5(s, spos, clipmax))
						return false;
					/*if (s->clp5<spos) {
					 if (badClipping(s,spos,s->clp3,clipmax)) return false;
					 clipops.Add(new SeqClipOp(s,spos));
					 }*/
				}
			} else { //trim right side in this msa
				//int newclp=s->seqlen-spos-1;
				if (s->revcompl != 0) {
					if (!clipops.add5(s, s->seqlen - spos - 1, clipmax))
						return false;
					/*if (s->clp5<newclp) {
					 if (badClipping(s,newclp,s->clp3,clipmax)) return false;
					 clipops.Add(new SeqClipOp(s,newclp));
					 }*/
				} else {
					if (!clipops.add3(s, s->seqlen - spos - 1, clipmax))
						return false;
					/*if (s->clp3<newclp) {
					 if (badClipping(s,s->clp5,newclp,clipmax)) return false;
					 clipops.Add(new SeqClipOp(s,-1, newclp));
					 }*/
				}
			}
		}         //for each sequense in MSA
	}         // 3' clipping
	return true;
}

void GSeqAlign::revComplement() {
	for (int i = 0; i < Count(); i++) {
		GASeq* s = Get(i);
		s->revComplement(length);
	}
	Sort();
}

void GSeqAlign::finalize() { //prepare for printing
  for (int i=0;i<Count();i++) {
	 GASeq* s=Get(i);
	 if (s->len==0)  GError("Error: sequence for %s not loaded!\n",s->getId());
	 if (!s->hasFlag(GA_flag_PREPPED)) s->prepSeq();
  }
}
void GSeqAlign::print(FILE* f, char c) {
	finalize(); //this calls prepSeq as needed to reverse complement sequence etc.
	int max = 0;
	for (int i = 0; i < Count(); i++) {
		int n = Get(i)->getNameLen();
		if (n > max)
			max = n;
	}
	char fmtstr[128];
	fmtstr[0] = '%';
	sprintf(&fmtstr[1], "%d", max);
	strcat(fmtstr, "s %c ");
	if (c != 0) {         // draw a separator line built from c
		fprintf(f, fmtstr, " ", ' ');
		for (int k = 0; k < length; k++)
			fprintf(f, "%c", c);
		fprintf(f, "\n");
	}
	for (int i = 0; i < Count(); i++) {
		GASeq* s = Get(i);
		char orientation = s->revcompl == 1 ? '-' : '+';
		fprintf(f, fmtstr, s->name(), orientation);
		s->printGappedSeq(f, minoffset);
	}
}

void GSeqAlign::writeMSA(FILE* f, int linelen) {
	finalize();
	for (int i = 0; i < Count(); i++) {
		GASeq& s = *Get(i);
		//
		s.printMFasta(f, linelen);
	}
}

char GAlnColumn::bestChar(int16_t *qscore) {  //returns most frequent char -- could be a gap!
	if (layers == 0)
		return 0;
	if (consensus != 0)
		return consensus;
	if (!countsSorted) {
		qsort(counts, 6, sizeof(NucCount), qsortnuc);
		countsSorted = true;
	}
	int r = 0;
	char best = counts[0].nt;
	int bq = counts[0].count;
	for (;r<5;) {
		//if gap or N have the same freq as a real base, pick the base instead
		if ((best == '-' || best == 'N') &&
		     counts[r].count == counts[r + 1].count) {
			r++;
			best = counts[r].nt;
			bq = counts[r].count;
		} else
			break;
	}
	if (qscore!=NULL) {
	   for (int q=0;q<6;++q) {
		   if (q!=r) bq-=counts[q].count;
	   }
	   if (bq<=SHRT_MIN) bq=SHRT_MIN+1;
	   if (bq>=SHRT_MAX) bq=SHRT_MAX-1;
	   *qscore=(short)bq;
	}
	consensus = best;
	return best;
}

void GAlnColumn::remove() {
	if (hasClip) {
		clipnuc->seq->msa->removeBase(clipnuc->seq, clipnuc->pos);
		return;
	}
	if (nucs->Count() > 0) {
		NucOri* n = nucs->Get(0);
		n->seq->msa->removeBase(n->seq, n->pos);
		//this should also be enough to propagate the deletion
		// to all involved sequences!
		// (all affected ofs[] and offsets)
		return;
	}
	GMessage(
	    "Warning: column remove() couldn't find a sequence at that position!\n");
}

void GSeqAlign::buildMSA(bool refWeighDown) {
	if (msacolumns != NULL)
		GError("Error: cannot call buildMSA() twice!\n");
	msacolumns = new MSAColumns(length, minoffset);
	for (int i = 0; i < Count(); i++) {
		GASeq* seq = Get(i);
		seq->msaidx = i; // if GSeqAlign is sorted by offset
		                 // this could speed up some later adjustments
		if (seq->seqlen - seq->clp3 - seq->clp5 < 1) {
			GMessage("Warning: sequence %s (length %d) was trimmed too badly (%d,%d)"
					 " -- should be removed from MSA w/ %s!\n", seq->id, seq->seqlen,
			    seq->clp5, seq->clp3, Get(0)->id);
			seq->setFlag(GA_flag_BAD_ALIGN); //bad-align flag!
			badseqs++;
		}
		int incVal=1;
		if (refWeighDown && !seq->hasFlag(GA_flag_IS_REF)) {
			incVal = 10;
		}
		seq->toMSA(*msacolumns, incVal);
	}
	//this->Pack();
}

void GSeqAlign::freeMSA() {
	if (msacolumns != NULL) {
		delete msacolumns;
		msacolumns = NULL;
	}
	//free sequence data too!
	for (int i = 0; i < Count(); i++) {
		GASeq* seq = Get(i);
		char* p = seq->detachSeqPtr();
		GFREE(p);
	}
}

void GSeqAlign::ErrZeroCov(int col) {
	int cnt = Count();
	fprintf(stderr,
	    "WARNING: 0 coverage column %d (mincol=%d) found within alignment of %d seqs!\n",
	    col, msacolumns->mincol, cnt);
	for (int i = 0; i < cnt; i++) {
		GASeq* seq = Get(i);
		fprintf(stderr, "%s\n", seq->id);
	}
	exit(5);
}

void GSeqAlign::refineMSA(bool refWeighDown, bool redo_ends) {
	if (redo_ends) {
		//TODO:
		//recompute consensus at the ends of MSA INCLUDING trimmed sequence
		//and recompute the trimming accordingly
	} else { //freeze end trimming
		//if (msacolumns==NULL)
		buildMSA(refWeighDown); //populate MSAColumns only based on existing trimming
	}
	//==> remove columns and build consensus
	int cols_removed = 0;
	for (int col = msacolumns->mincol; col <= msacolumns->maxcol; col++) {
		int16_t qscore=0;
		char c = msacolumns->columns[col].bestChar(&qscore);
		if (c == 0) { //should never be the case!
			ErrZeroCov(col);
			c = '*';
		}
		if (c == '-' || c == '*') {
			c = '*';
			if (MSAColumns::removeConsGaps) {
				removeColumn(col - cols_removed);
				cols_removed++;
				//this will delete the corresponding nucleotides
				//from every involved read, also updating the offsets of
				//every read AFTER this column
				continue;//don't add this gap to the consensus
			}
		}
		extendConsensus(c, qscore);
	}
	//make sure consensus is 0 terminated:
	char e=0;consensus.Add(e);consensus.Pop();
	//-- refine clipping and remove gaps propagated in the clipping regions
	for (int i = 0; i < Count(); i++) {
		GASeq* seq = Get(i);
		//if (seq->hasFlag(7)) continue; -- checking the badalign flag..
		//refine clipping -- first pass:
		if (MSAColumns::refineClipping)
			seq->refineClipping(consensus, //consensus_len,
			    seq->offset - minoffset - msacolumns->mincol);

		//..remove any "gaps" in the non-aligned (trimmed) regions
		int grem = 0;
		if (MSAColumns::removeConsGaps)
			grem = seq->removeClipGaps();
		//if any gaps were removed, take one more shot at
		//refining the clipping -- we may get lucky and realign better..
		if (grem != 0 && MSAColumns::refineClipping)
			seq->refineClipping(consensus, //consensus_len,
			    seq->offset - minoffset - msacolumns->mincol, true);
	}
	refinedMSA = true;
}

void GSeqAlign::extendConsensus(char c, int16_t bq) {
	/*
	int newlen = consensus_len + 1;
	if (newlen >= consensus_cap) {
		consensus_cap += 128;
		if (consensus_len == 0) {
			GMALLOC(consensus, consensus_cap);
		} else {
			GREALLOC(consensus, consensus_cap);
		}
	}
	consensus[consensus_len] = c;
	consensus[newlen] = 0;
	consensus_len++;
	*/
	consensus.Add(c);
	if (bq!=SHRT_MIN && consensus_bq.Count()==consensus.Count()-1) {
		consensus_bq.Add(bq);
	}
}

void GSeqAlign::writeACE(FILE* f, const char* name, bool refWeighDown) {
	//--build a consensus sequence
	if (!refinedMSA)
		refineMSA(refWeighDown);

	//FastaSeq conseq((char*)name);
	//conseq.setSeqPtr(consensus, consensus_len, consensus_cap);
	int fwd = 0; //number of reversed reads
	int rvs = 0; //number of forward reads
	for (int i = 0; i < Count(); i++) {
		GASeq* seq = Get(i);
		if (seq->revcompl != 0)
			rvs++;
		else
			fwd++;
	}
	char consDir = (rvs > fwd) ? 'C' : 'U';

	fprintf(f, "CO %s %d %d 0 %c\n", name, consensus.Count(), Count(), consDir);
	//   conseq.len, Count()-badseqs, consDir);
	//conseq.fprint(f,60);
	FastaSeq::write(f, NULL, NULL, consensus(), 60, consensus.Count());
	//fprintf(f, "\nBQ \n\n"); //TODO: print consensus_bq array values here!
	fprintf(f, "\nBQ\n");
	int bl=0;
	for (uint i=0;i<consensus_bq.Count();++i) {
		if (bl) fprintf(f, " %d", consensus_bq[i]);
		   else fprintf(f,  "%d", consensus_bq[i]);
		++bl;
		if (bl==60) {
			fprintf(f,"\n");
			bl=0;
		}
	}
	if (bl) fprintf(f, "\n\n");
	   else fprintf(f, "\n");
	for (int i = 0; i < Count(); i++) {
		GASeq* seq = Get(i);
		//if (seq->hasFlag(7)) continue; -- checking the badalign flag..
		char sc = (seq->revcompl == 0) ? 'U' : 'C';
		fprintf(f, "AF %s %c %d\n", seq->id, sc,
		    seq->offset - minoffset - msacolumns->mincol + 1);
	}
	fprintf(f, "\n");
	// a second pass to write the actual read entries..
	for (int i = 0; i < Count(); i++) {
		GASeq* seq = Get(i);
		//if (seq->hasFlag(7)) continue; //badalign flag set
		int gapped_len = seq->seqlen + seq->numgaps;
		fprintf(f, "RD %s %d 0 0\n", seq->id, gapped_len);
		seq->printGappedFasta(f);
		int clpl, clpr;
		if (seq->revcompl == 0) { //forward
			clpl = seq->clp5;
			clpr = seq->clp3;
		} else { //reverse complement
			clpl = seq->clp3;
			clpr = seq->clp5;
		}
		int l = clpl;
		int r = clpr;
		for (int j = 1; j <= r; j++)
			clpr += seq->ofs[seq->seqlen - j];
		for (int j = 0; j <= l; j++)
			clpl += seq->ofs[j];
		int seql = clpl + 1;
		int seqr = gapped_len - clpr;
		if (seqr < seql) {
			fprintf(stderr, "Bad trimming for %s of gapped len %d (%d, %d)\n",
			    seq->id, gapped_len, seql, seqr);
			seqr = seql + 1;
		}
		fprintf(f, "\nQA %d %d %d %d\nDS \n\n", seql, seqr, seql, seqr);
	}

}

void GSeqAlign::writeInfo(FILE* f, const char* name, bool refWeighDown) {
	/*
	 File format should match assembly & asmbl_link tables in our db:

	 >contig_name seq_count contig_sequence
	 seqname seqlen offset asm_lend asm_rend seq_lend seq_rend pid alndata

	 Notes:
	 seq_lend>seq_rend if the sequence is reverse complemented

	 */
//--build the actual MSA and a consensus sequence, if not done yet:
// this will also remove the consensus gaps as appropriate (unless disabled)
	if (!refinedMSA)
		refineMSA(refWeighDown);
//-- also compute this, just in case:
// redundancy = sum(asm_rend-asm_lend+1)/contig_len
//(and also the pid for each reads vs. consensus)
	fprintf(f, ">%s %d %s\n", name, Count(), consensus());
	float redundancy = 0; // = sum(asm_rend-asm_lend+1)/contig_len
	for (int i = 0; i < Count(); i++) {
		GASeq* seq = Get(i);
		//GStr alndata;
		//if (seq->hasFlag(7)) continue; //badalign flag set
		int gapped_len = seq->seqlen + seq->numgaps;
		//fprintf(f, "RD %s %d 0 0\n", seq->id, gapped_len);
		/*seq->printGappedFasta(f);*/
		int seqoffset = seq->offset - minoffset - msacolumns->mincol + 1;
		int clpl, clpr;
		int asml = seqoffset + 1;
		int asmr = asml - 1;
		float pid = 0;
		if (seq->revcompl == 0) { //forward
			clpl = seq->clp5;
			clpr = seq->clp3;
		} else { //reverse complement
			clpl = seq->clp3;
			clpr = seq->clp5;
		}
		int aligned_len = 0;
		//int indel_ofs = 0; //distance to last indel position
		for (int j = seq->clp5; j < seq->seqlen - seq->clp3; j++) {
			int indel = seq->ofs[j];
			//char indel_type = 0;
			asmr += indel + 1;
			if (indel < 0) { //deletion
				//indel_type = 'd';
				indel = -indel;
			} else { //  indel>=0
				//actually aligned nucleotide here
				//if (indel > 0)
				//	indel_type = 'g';
				//else
				//	// indel==0, no indel at all
				//	indel_ofs++;
				if (toupper(seq->seq[j]) == toupper(consensus[asmr - 1]))
					pid++;
				aligned_len++;
			}
			/*
			if (indel_type) {
				if (indel > 2)
					alndata.appendfmt("%d%c%d-", indel_ofs, indel_type, indel);
				else
					for (int r = 0; r < indel; r++)
						alndata += indel_type;
				indel_ofs = 0;
			}
			*/
		}
		pid = (pid * 100.0) / (float) aligned_len;
		redundancy += aligned_len;
		/*int l=clpl;
		 int r=clpr;
		 for (int j=1;j<=r;j++) clpr+=seq->ofs[seq->seqlen-j];
		 for (int j=0;j<=l;j++) clpl+=seq->ofs[j];*/
		int seql = clpl + 1;
		int seqr = seq->len - clpr;
		if (seqr < seql) {
			fprintf(stderr,
			    "WARNING: Bad trimming for %s of gapped len %d (%d, %d)\n", seq->id,
			    gapped_len, seql, seqr);
			seqr = seql + 1;
		}
		if (seq->revcompl)
			Gswap(seqr, seql);
		//         id ln of al ar sl sr pi an
		fprintf(f, "%s %d %d %d %d %d %d %4.2f\n", seq->id, seq->len, seqoffset,
		    asml, asmr, seql, seqr, pid); //alndata.chars());
	}
	redundancy /= (float) consensus.Count();
}
