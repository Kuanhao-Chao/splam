#include "GBase.h"
#include "GArgs.h"
#include "GStr.h"
#include "GVec.hh"
namespace old {
 #include "GHash.hh"
}

#include "GResUsage.h"
#include <iostream>
//#include "tsl/hopscotch_map.h"
//#include "tsl/robin_map.h"
#include <unordered_map>
//#include "ska/bytell_hash_map.hpp"
#include "GHashMap.hh"

#define USAGE "Usage:\n\
  htest [-Q] [-C] [-n num_clusters] textfile.. \n\
  E.g. quick query test: ./htest -Q qtest_str.dta\n\
  \n\
 "
//quick query test: ./htest -Q qtest_str.dta

bool qryMode=false;
bool checkRM=false;
int numClusters=500;

struct HStrData {
	int cmd; // 0=add, 1=remove, 2=clear
	GStr str;
	HStrData(char* s=NULL, int c=0):cmd(c), str(s) { }
};

int loadStrings(FILE* f, GPVec<HStrData>& strgsuf, GPVec<HStrData>& strgs, int toLoad) {
  int num=0;
  GLineReader lr(f);
  char* line=NULL;
  int numcl=0;
  while ((line=lr.nextLine())!=NULL) {
	  int len=strlen(line);
	  if (len<3) continue;
	  if (line[0]=='>') {
		  numcl++;
		  if (toLoad && numcl>toLoad) {
			  break;
		  }
		  continue;
	  }
	  if (strcmp(line, "HCLR")==0) {
		  strgs.Add(new HStrData(NULL, 2));
		  strgsuf.Add(new HStrData(NULL, 2));
		  continue;
	  }
	  if (startsWith(line, "RM ")) {
	     strgsuf.Add(new HStrData(line+3,1) );
	     line[len-3]=0;
	     strgs.Add(new HStrData(line+3,1));
	     continue;
	  }
      strgsuf.Add(new HStrData(line));
      line[len-3]=0;
      strgs.Add(new HStrData(line));
	  num++;
  } //while line
  return num;
}

void showTimings(GResUsage swatch) {
 char *wtime=commaprintnum((uint64_t)swatch.elapsed());
 char *utime=commaprintnum((uint64_t)swatch.u_elapsed());
 char *stime=commaprintnum((uint64_t)swatch.s_elapsed());
 char *smem=commaprintnum((uint64_t)swatch.memoryUsed());
 GMessage("Elapsed time (microseconds): %12s us\n", wtime);
 GMessage("                  user time: %12s us\n", utime);
 GMessage("                system time: %12s us\n", stime);
 GMessage("                  mem usage: %12s KB\n", smem);

 GFREE(wtime);GFREE(utime);GFREE(stime); GFREE(smem);
}

// default values recommended by http://isthe.com/chongo/tech/comp/fnv/
const uint32_t Prime = 0x01000193; // 16777619
const uint32_t Seed  = 0x811C9DC5; // 2166136261
/// hash a single byte
inline uint32_t fnv1a(unsigned char b, uint32_t h = Seed) {   return (b ^ h) * Prime; }

/// hash a C-style string
uint32_t fnv1a(const char* text, uint32_t hash = Seed) {
	while (*text)
		hash = fnv1a((unsigned char)*text++, hash);
	return hash;
}

struct cstr_eq {
    inline bool operator()(const char* x, const char* y) const {
        return (strcmp(x, y) == 0);
    }
};

struct cstr_hash {
	inline uint32_t operator()(const char* s) const {
		return XXH32(s, std::strlen(s),0);
		//return fnv1a(s);
	}
};

void run_GHash(GResUsage& swatch, GPVec<HStrData> & hstrs, const char* label) {
	old::GHash<int> ghash; // @suppress("Type cannot be resolved")
	int num_add=0, num_rm=0, num_clr=0;
	GMessage("----------------- %s ----------------\n", label);
	ghash.Clear();
	swatch.start();
	int cl_i=0;
	int prevcmd=2;
	for (int i=0;i<hstrs.Count();i++) {
			  if (hstrs[i]->cmd==prevcmd) {
				  if (prevcmd==2) continue;
			  } else prevcmd=hstrs[i]->cmd;
			  switch (hstrs[i]->cmd) {
				case 0:
					if (cl_i==0) cl_i=i;
					ghash.fAdd(hstrs[i]->str.chars(), new int(i));
					num_add++;
					break;
				case 1:
					if (qryMode) break;
					ghash.Remove(hstrs[i]->str.chars());
					num_rm++;
					break;
				case 2:
					//run tests here
					if (qryMode) {
						//run some query tests here
						for(int j=cl_i;j<i;j+=3) {
							if (hstrs[j]->cmd) continue;
							int* v=ghash[hstrs[j]->str.chars()];
							if (v==NULL)
								GError("Error at <%s>, key %s not found (count:%d, cl_i=%d, i=%d)!\n",label, hstrs[j]->str.chars(),
										ghash.Count(), cl_i, i );
							if (*v!=j)
								GError("Error at <%s>, invalid value for key %s!\n",label, hstrs[j]->str.chars() );
						}
					}
					cl_i=0;
					ghash.Clear();
					num_clr++;
					break;
			  }
	}
	swatch.stop();
	ghash.Clear();
	GMessage("  (%d inserts, %d deletions, %d clears)\n", num_add, num_rm, num_clr);
}
/*
void run_Hopscotch(GResUsage& swatch, GPVec<HStrData> & hstrs, const char* label) {
  int num_add=0, num_rm=0, num_clr=0;
  //tsl::hopscotch_map<const char*, int, cstr_hash, cstr_eq> hsmap;
  tsl::hopscotch_map<const char*, int, cstr_hash, cstr_eq,
      std::allocator<std::pair<const char*, int>>,  30, true> hsmap;
  GMessage("----------------- %s ----------------\n", label);
  swatch.start();
  int cl_i=0;
  int prevcmd=2;
  for (int i=0;i<hstrs.Count();i++) {
	  if (hstrs[i]->cmd==prevcmd) {
		  if (prevcmd==2) continue;
	  } else prevcmd=hstrs[i]->cmd;
	  switch (hstrs[i]->cmd) {
		case 0:
			if (cl_i==0) cl_i=i;
			hsmap.insert({hstrs[i]->str.chars(), i});
			num_add++;
			break;
		case 1:
			if (qryMode) break;
			hsmap.erase(hstrs[i]->str.chars());
			num_rm++; break;
		case 2:
			if (qryMode) {
				//run some query tests here
				//with strings from hstrs[cl_i .. i-1] range
				for(int j=cl_i;j<i;j+=3) {
					if (hstrs[j]->cmd) continue;
					int v=hsmap[hstrs[j]->str.chars()];
					if (v!=j)
						GError("Error at <%s>, invalid value for key %s! (got %d, expected %d)\n",label,
								hstrs[j]->str.chars(), v, j );
				}
			}
			cl_i=0;
			hsmap.clear();
			num_clr++;
			break;
	  }
  }
  swatch.stop();
  hsmap.clear();
  GMessage("  (%d inserts, %d deletions, %d clears)\n", num_add, num_rm, num_clr);
}

void run_Robin(GResUsage& swatch, GPVec<HStrData> & hstrs, const char* label) {
  int num_add=0, num_rm=0, num_clr=0;
  //tsl::hopscotch_map<const char*, int, cstr_hash, cstr_eq> hsmap;
  tsl::robin_map<const char*, int, cstr_hash, cstr_eq,
      std::allocator<std::pair<const char*, int>>,  true> rmap;
  GMessage("----------------- %s ----------------\n", label);
  swatch.start();
  int cl_i=0;
  int prevcmd=2;
  for (int i=0;i<hstrs.Count();i++) {
	  if (hstrs[i]->cmd==prevcmd) {
		  if (prevcmd==2) continue;
	  } else prevcmd=hstrs[i]->cmd;
	  switch (hstrs[i]->cmd) {
		case 0:
			if (cl_i==0) cl_i=i;
			rmap.insert({hstrs[i]->str.chars(), i});
			num_add++;
			break;
		case 1: if (qryMode) break;
			rmap.erase(hstrs[i]->str.chars()); num_rm++; break;
		case 2:
			if (qryMode) {
				//run some query tests here
				//with strings from hstrs[cl_i .. i-1] range
				for(int j=cl_i;j<i;j+=3) {
					if (hstrs[j]->cmd) continue;
					int v=rmap[hstrs[j]->str.chars()];
					if (v!=j)
						GError("Error at <%s>, invalid value for key %s!\n",label, hstrs[j]->str.chars() );
				}
			}
			cl_i=0;
			rmap.clear(); num_clr++; break;
	  }
  }
  swatch.stop();
  rmap.clear();
  GMessage("  (%d inserts, %d deletions, %d clears)\n", num_add, num_rm, num_clr);
}

void run_Bytell(GResUsage& swatch, GPVec<HStrData> & hstrs, const char* label) {
  int num_add=0, num_rm=0, num_clr=0;
  ska::bytell_hash_map<const char*, int, cstr_hash, cstr_eq> bmap;
  GMessage("----------------- %s ----------------\n", label);
  swatch.start();
  for (int i=0;i<hstrs.Count();i++) {
	  switch (hstrs[i]->cmd) {
		case 0:bmap.insert({hstrs[i]->str.chars(), 1}); num_add++; break;
		case 1:bmap.erase(hstrs[i]->str.chars()); num_rm++; break;
		case 2:bmap.clear(); num_clr++; break;
	  }
  }
  swatch.stop();
  bmap.clear();
  GMessage("  (%d inserts, %d deletions, %d clears)\n", num_add, num_rm, num_clr);
}
*/
void run_Khashl(GResUsage& swatch, GPVec<HStrData> & hstrs, const char* label) {
  int num_add=0, num_rm=0, num_clr=0;
  klib::KHashMapCached<const char*, int, cstr_hash, cstr_eq > khmap;
  GMessage("----------------- %s ----------------\n", label);
  swatch.start();
  int cl_i=0;
  int prevcmd=2;
  for (int i=0;i<hstrs.Count();i++) {
	  if (hstrs[i]->cmd==prevcmd) {
		  if (prevcmd==2) continue;
	  } else prevcmd=hstrs[i]->cmd;
	  switch (hstrs[i]->cmd) {
		case 0:if (cl_i==0) cl_i=i;
			khmap[hstrs[i]->str.chars()]=i; num_add++; break;
		case 1:if (qryMode) break;
			khmap.del(khmap.get(hstrs[i]->str.chars())); num_rm++; break;
		case 2:
			if (qryMode) {
				//run some query tests here
				for(int j=cl_i;j<i;j+=3) {
					if (hstrs[j]->cmd) continue;
					int v=khmap[hstrs[j]->str.chars()];
					if (v!=j)
						GError("Error at <%s>, invalid value for key %s!\n",label, hstrs[j]->str.chars() );
				}
			}
			cl_i=0;
			khmap.clear(); num_clr++; break;
	  }
  }
  swatch.stop();
  khmap.clear();
  GMessage("  (%d inserts, %d deletions, %d clears)\n", num_add, num_rm, num_clr);
}

void run_GHashMap(GResUsage& swatch, GPVec<HStrData> & hstrs, const char* label) {
  int num_add=0, num_rm=0, num_clr=0;
  //GKHashSet<const char*> khset;
  //GHashSet<> khset;
  //GHash<int, cstr_hash, GHashKey_Eq<const char*>, uint32_t> khset;
  GHashMap<const char*, int> khset;
  GMessage("----------------- %s ----------------\n", label);
  int cl_i=0;
  swatch.start();
  int prevcmd=2;
  for (int i=0;i<hstrs.Count();i++) {
	  if (hstrs[i]->cmd==prevcmd) {
		  if (prevcmd==2) continue;
	  } else prevcmd=hstrs[i]->cmd;
	  switch (hstrs[i]->cmd) {
		case 0: if (cl_i==0) cl_i=i;
			khset.Add(hstrs[i]->str.chars(), i); num_add++; break;
		case 1:if (qryMode) break;
			if (khset.Remove(hstrs[i]->str.chars())<0)
				if (checkRM) GMessage("Warning: key %s could not be removed!\n", hstrs[i]->str.chars());
			num_rm++;
			break;
		case 2:
			if (qryMode) {
				//run some query tests here
				//with strings from hstrs[cl_i .. i-1] range
				for(int j=cl_i;j<i;j+=3) {
					if (hstrs[j]->cmd) continue;
					int* v=khset[hstrs[j]->str.chars()];
					if (*v!=j)
						GError("Error at <%s>, invalid value for key %s!\n",label, hstrs[j]->str.chars() );
				}
			}
			cl_i=0;
			khset.Clear(); num_clr++; break;
	  }
  }
  swatch.stop();
  khset.Clear();
  GMessage("  (%d inserts, %d deletions, %d clears)\n", num_add, num_rm, num_clr);
}

void run_GxxHashMap(GResUsage& swatch, GPVec<HStrData> & hstrs, const char* label) {
  int num_add=0, num_rm=0, num_clr=0;
  //GHash<int> khset;
  GHashMap<const char*, int, GHashKey_xxHash32<const char*>,
          GHashKey_Eq<const char*>, uint32_t > khset;
  GMessage("----------------- %s ----------------\n", label);
  int cl_i=0;
  swatch.start();
  int prevcmd=2;
  for (int i=0;i<hstrs.Count();i++) {
	  if (hstrs[i]->cmd==prevcmd) {
		  if (prevcmd==2) continue;
	  } else prevcmd=hstrs[i]->cmd;
	  switch (hstrs[i]->cmd) {
		case 0: if (cl_i==0) cl_i=i;
			khset.Add(hstrs[i]->str.chars(), i); num_add++; break;
		case 1:if (qryMode) break;
			if (khset.Remove(hstrs[i]->str.chars())<0)
				if (checkRM) GMessage("Warning: key %s could not be removed!\n", hstrs[i]->str.chars());
			num_rm++;
			break;
		case 2:
			if (qryMode) {
				//run some query tests here
				//with strings from hstrs[cl_i .. i-1] range
				for(int j=cl_i;j<i;j+=3) {
					if (hstrs[j]->cmd) continue;
					int* v=khset[hstrs[j]->str.chars()];
					if (*v!=j)
						GError("Error at <%s>, invalid value for key %s!\n",label, hstrs[j]->str.chars() );
				}
			}
			cl_i=0;
			khset.Clear(); num_clr++; break;
	  }
  }
  swatch.stop();
  khset.Clear();
  GMessage("  (%d inserts, %d deletions, %d clears)\n", num_add, num_rm, num_clr);
}

void run_GHashMapShk(GResUsage& swatch, GPVec<HStrData> & hstrs, const char* label) {
  int num_add=0, num_rm=0, num_clr=0;
  GHashMap<const char*, int> khset;
  GMessage("----------------- %s ----------------\n", label);
  int cl_i=0;
  swatch.start();
  int prevcmd=2;
  for (int i=0;i<hstrs.Count();i++) {
	  if (hstrs[i]->cmd==prevcmd) {
		  if (prevcmd==2) continue;
	  } else prevcmd=hstrs[i]->cmd;
	  switch (hstrs[i]->cmd) {
		case 0: if (cl_i==0) cl_i=i;
			khset.Add(hstrs[i]->str.chars(),i); num_add++; break;
		case 1:if (qryMode) break;
			if (khset.Remove(hstrs[i]->str.chars())<0)
				if (checkRM) GMessage("Warning: key %s could not be removed!\n", hstrs[i]->str.chars());
			num_rm++;
			break;
		case 2:
			if (qryMode) {
				//run some query tests here
				//with strings from hstrs[cl_i .. i-1] range
				for(int j=cl_i;j<i;j+=3) {
					if (hstrs[j]->cmd) continue;
					int* v=khset[hstrs[j]->str.chars()];
					if (*v!=j)
						GError("Error at <%s>, invalid value for key %s!\n",label, hstrs[j]->str.chars() );
				}
			}
			cl_i=0;
			khset.Clear(); num_clr++; break;
	  }
  }
  swatch.stop();
  khset.Clear();
  GMessage("  (%d inserts, %d deletions, %d clears)\n", num_add, num_rm, num_clr);
}

struct SObj {
  GStr atr;
  int val;
  SObj(const char* a=NULL, const int v=0):atr(a),val(v) { }
  bool operator<(const SObj& o) const { return val<o.val; }
  bool operator==(const SObj& o) const {
	  return (atr==o.atr && val==o.val);
  }
};

int main(int argc, char* argv[]) {
 GPVec<HStrData> strs;
 GPVec<HStrData> sufstrs;
 //GArgs args(argc, argv, "hg:c:s:t:o:p:help;genomic-fasta=COV=PID=seq=out=disable-flag;test=");
 GArgs args(argc, argv, "hQCn:");
 //fprintf(stderr, "Command line was:\n");
 //args.printCmdLine(stderr);
 args.printError(USAGE, true);
 if (args.getOpt('h') || args.getOpt("help")) GMessage(USAGE);
 GStr s=args.getOpt('n');
 if (!s.is_empty()) {
	 numClusters=s.asInt();
	 if (numClusters<=0)
		 GError("%s\nError: invalid value for -n !\n", USAGE);
 }
 qryMode=(args.getOpt('Q'));
 checkRM=(args.getOpt('C'));
 int numargs=args.startNonOpt();
 const char* a=NULL;
 FILE* f=NULL;
 int total=0;
//==== quick test area
 /*
 std::unordered_map<SObj*, int > umap;
 GHash<int> gh;
 GPVec<SObj> ptrs(false);
 GQHash<int, SObj*> ihash(false);
 GQHash<SObj*, int> phash;
 GQStrHash<SObj*> shash;
 const char* tstrs[6] = {"twelve", "five", "nine", "eleven", "three", "nope"};
 int  vals[6]     = { 12,  5, 9, 11, 3, 777 };
 char buf[20];
 for (int i=0;i<5;i++) {
   SObj* o=new SObj(tstrs[i], vals[i]*10);
   ptrs.Add(o);
   sprintf(buf, "%lx", o);
   GMessage("SObj (%s, %d) pointer added: %s\n",tstrs[i], o->val, buf);
   gh.Add(buf, new int(vals[i]));
   shash.Add(tstrs[i], o);
   ihash.Add(vals[i], o);
   phash.Add(o, vals[i]);
   umap[o]=vals[i];
 }
 ptrs.Sort();
 GMessage("shash has now %d entries.\n", shash.Count());
 //enumerate shash entries:
 {
   shash.startIterate();
   SObj* iv=NULL;
   while (const char* k=shash.Next(iv)) {
	   GMessage("Enumerating shash entry: (%s => %lx)\n",
			    k, iv);
   }
 }
 //qry:
 for (int i=0;i<ptrs.Count();i++) {
   SObj* o=ptrs[i];
   //test tset
   SObj* v=shash.Find(o->atr.chars());
   if (v==NULL)
	   GMessage("key <%s> not found in shash!\n", o->atr.chars());
   int* iv=phash.Find(o);
   if (iv==NULL)
	   GMessage("key <%lx> not found in phash!\n", o);
   //if (!oset[*o])
//   GMessage("struct {%s, %d} not found in oset!\n", o->atr.chars(), o->val);

   //sprintf(buf, "%lx", o);
   //int* hv=gh[buf];
   //GMessage("Item {%s, %d} : GHash retrieved flag = %d, umap retrieved flag = %d\n",
   //  o->atr.chars(), o->val, *hv, umap[o]);
 }
//SObj* n=new SObj("test", 10);
//if (!pset[n])
//	   GMessage("key <%lx> not found in pset!\n", n);

 for (int i=0;i<6;i++) {
	SObj* o=shash[tstrs[i]];
	if (o==NULL) GMessage("key <%s> not found in shash!\n", tstrs[i]);
	if (o && i<5) {
		if (o->atr!=tstrs[i])
			GMessage("shash value does not match key <%s!\n", tstrs[i]);
	}
 }

 //delete n;
 //int v=umap[n];
 //GMessage("Non-existing test entry returned value %d\n", v);
  */
 /*
 auto found=umap.find(n);
 if (found!=umap.end()) {
   GMessage("Found flags %d for entry {\"%s\", %d}\n", found->second,
         n->atr.chars(), found->first->val );
 } else  GMessage("New test obj not found !\n");

 return(0);
*/
//==== quick test area end
 if (numargs==0) {
	 //a="htest_data.lst";
	 a="htest_over500.lst";
	 f=fopen(a, "r");
	 if (f==NULL) GError("Error: could not open file %s !\n", a);
	 GMessage("loading %d clusters from file..\n", numClusters);
	 int num=loadStrings(f, sufstrs, strs, numClusters);
	 total+=num;
	 fclose(f);
 }
 else {
	   while ((a=args.nextNonOpt())) {
		   f=fopen(a, "r");
		   if (f==NULL) GError("Error: could not open file %s !\n", a);
		   int num=loadStrings(f, sufstrs, strs, numClusters);
		   total+=num;
		   fclose(f);
	   }
  }
   GResUsage swatch;


   run_GHash(swatch, sufstrs, "GHash w/ suffix");
   showTimings(swatch);
   //run_GHash(swatch, strs, "GHash no suffix");
   //showTimings(swatch);

/*
   run_Hopscotch(swatch, sufstrs, "hopscotch w/ suffix");
   showTimings(swatch);
   run_Hopscotch(swatch, strs, "hopscotch no suffix");
   showTimings(swatch);
*/
/*
   run_Robin(swatch, sufstrs, "robin w/ suffix");
   showTimings(swatch);
   run_Robin(swatch, strs, "robin no suffix");
   showTimings(swatch);
*/

   run_Khashl(swatch, sufstrs, "khashl w/ suffix");
   showTimings(swatch);
/*
   run_Khashl(swatch, strs, "khashl no suffix");
   showTimings(swatch);
*/
   run_GHashMap(swatch, sufstrs, "GHashMap default w/ suffix");
   showTimings(swatch);

   run_GxxHashMap(swatch, sufstrs, "GHashMap xxHash32 w/ suffix");
   showTimings(swatch);

   //run_GHashMap(swatch, strs, "GHashMap no suffix");
   //showTimings(swatch);

   //run_GHashMapShk(swatch, sufstrs, "GHashSetShk w/ suffix");
   //showTimings(swatch);
   //run_GHashMapShk(swatch, strs, "GHashSetShk no suffix");
   //showTimings(swatch);

/*
   run_Bytell(swatch, sufstrs, "bytell w/ suffix");
   showTimings(swatch);
   run_Bytell(swatch, strs, "bytell no suffix");
   showTimings(swatch);
*/

}
