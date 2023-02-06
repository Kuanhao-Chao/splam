#include "GBase.h"
#include "GArgs.h"
#include "GStr.h"
#include "GBitVec.h"
#include "GList.hh"
#include "GHash.hh"

#define USAGE "Usage:\n\
gtest [--bit-test|-g|--genomic-fasta <genomic_seqs_fasta>] [-c|COV=<cov%>] \n\
 [-s|--seq <seq_info.fsize>] [-o|--out <outfile.gff>] [--disable-flag] [-t|--test <string>]\n\
 [-p|PID=<pid%>] file1 [file2 file3 ..]\n\
 "
enum {
 OPT_HELP=1,
 OPT_GENOMIC,
 OPT_COV,
 OPT_SEQ,
 OPT_OUTFILE,
 OPT_DISABLE_FLAG,
 OPT_TEST,
 OPT_PID,
 OPT_BITVEC,
 OPT_NUM
};

GArgsDef opts[] = {
{"help",          'h', 0, OPT_HELP},
{"genomic-fasta", 'g', 1, OPT_GENOMIC},
{"COV",           'c', 1, OPT_COV},
{"seq",           's', 1, OPT_SEQ},
{"out",           'o', 1, OPT_OUTFILE},
{"disable-flag",   0,  0, OPT_DISABLE_FLAG},
{"test",          't', 1, OPT_TEST},
{"PID",           'p', 1, OPT_PID},
{"bit-test",      'B', 0, OPT_BITVEC},
{"bignum",        'n', 1, OPT_NUM},
{0,0,0,0}
};

void bitError(int idx) {
 GError("Error bit checking (index %d)!\n", idx);
}

struct Gint {
  int v;
  Gint(int vv=0):v(vv) {}
  int val() { return v; }
  ~Gint() {
  GMessage("Gint with val %d getting destroyed\n", v);
  }
};

int cmpGint(pointer p1, pointer p2) {
 int v1=((Gint*)p1)->v;
 int v2=((Gint*)p2)->v;
 if (v1<v2) return -1;
 return (v1>v2)? 1 : 0;
}

void testGPVec() {
 GPVec<Gint> vecs[3];
 vecs[1].Add(new Gint(2));
 vecs[2].Add(new Gint(3));
 GMessage("Added to vecs[1]:%d\n", vecs[1][0]->val());
 GMessage("Added to vecs[2]:%d\n", vecs[2][0]->val());
}

int main(int argc, char* argv[]) {
 //GArgs args(argc, argv, "hg:c:s:t:o:p:help;genomic-fasta=COV=PID=seq=out=disable-flag;test=");
 GArgs args(argc, argv, opts);
 fprintf(stderr, "Command line was:\n");
 args.printCmdLine(stderr);
 args.printError(USAGE, true);
 //if (args.getOpt('h') || args.getOpt("help"))
 GVec<int> transcripts(true);
 transcripts.cAdd(0);
 fprintf(stderr,"after add transcript counts=%d\n",transcripts.Count());
 exit(0);
 
 if (args.getOpt(OPT_HELP))
     {
     GMessage("%s\n", USAGE);
     exit(1);
     }

 if (args.getOpt(OPT_NUM)) {
      GStr snum(args.getOpt(OPT_NUM));
      int num=snum.asInt();
      char* numstr=commaprintnum(num);
      GMessage("Number %d written with commas: %s\n", num, numstr);
      GFREE(numstr);
 }
 //---
 GHash<GVec<int> > ends;
 
 /*
 testGPVec();
 //exit(0);

 //uint pos=3;
 //GStr spos((int)pos);
 //GVec<int> *ev=ends[spos.chars()];
 
 GPVec<Gint> v;
 int r(5);
 int rr=v.Add(new Gint(3));
 //if (rr<0) {
 // GMessage("Error adding 0! (code %d)\n",rr);
 // }
 v.Add(new Gint(r));
 v.Add(new Gint(2));
 v.Add(new Gint(1));
 v.Add(new Gint(4));
 rr=v.Add(new Gint(0));
 v[rr]->v=-1;
 v.Sort(cmpGint);
 GMessage("collection has %d elements:\n",v.Count());
 for (int i=0;i<v.Count();i++) {
   GMessage("v[%d]=%d;\n",i,v[i]->v);
 }
 exit(0);
 */
 //---
 int numopts=args.startOpt();
 if (numopts)
   GMessage("#### Recognized %d option arguments:\n", numopts);
 int optcode=0;
 while ((optcode=args.nextCode())) {
   char* r=args.getOpt(optcode);
   GMessage("%14s\t= %s\n", args.getOptName(optcode), (r[0]==0)?"True":r);
   }
 int numargs=args.startNonOpt();
 if (numargs>0) {
   GMessage("\n#### Found %d non-option arguments given:\n", numargs);
   char* a=NULL;
   while ((a=args.nextNonOpt())) {
     GMessage("%s\n",a);
     }
   }
 GStr s=args.getOpt('t');
 if (!s.is_empty()) {
    GStr token;
    GMessage("Tokens in \"%s\" :\n",s.chars());
    s.startTokenize(";,: \t");
    int c=1;
    while (s.nextToken(token)) {
      GMessage("token %2d : \"%s\"\n",c,token.chars());
      c++;
      }
    }
 if (args.getOpt(OPT_BITVEC)) {
    uint numbits=4156888234;
    GBitVec bits(numbits);
    GMessage(">>> -- BitVec(%u) created (size=%u, mem=%lu) -- \n", numbits, bits.size(),
    		bits.getMemorySize());
    bits[405523342]=true;
    GMessage("      memory size: %lu , size()=%u, count()=%d \n", bits.getMemorySize(), bits.size(), bits.count());
    /*
    //GMessage(">>> -- Start BitVec Test -- \n");
    if (bits[1092]) bitError(1092);
    bits.resize(2049);
    if (bits[2048]) bitError(2048);
    bits[2048]=true;
    if (!bits[2048]) bitError(2048);
    bits.resize(4097);
    if (!bits[2048]) bitError(2048);
    if (bits[4096]) bitError(4096);
    bits[4096]=true;
    if (!bits[4096]) bitError(4096);
    GBitVec bits2(64);
    Gswap(bits, bits2);
    if (!bits2[2048]) bitError(2048);
    if (!bits2[4096]) bitError(4096);
    */
    //GMessage("<<< -- End BitVec Test (size: %d, count: %d, bits2 size=%d, count=%d) --\n",
    ///       bits.size(), bits.count(), bits2.size(), bits2.count());
    }
}
