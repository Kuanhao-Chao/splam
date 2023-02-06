#include "GSam.h"

// struct CJunc {
// 	int start, end;
// 	char strand;
// 	uint64_t dupcount;
// 	CJunc(int vs=0, int ve=0, char vstrand='+', uint64_t dcount=1):
// 	  start(vs), end(ve), strand(vstrand), dupcount(dcount) { }

// 	bool operator==(const CJunc& a) {
// 		return (strand==a.strand && start==a.start && end==a.end);
// 	}

//     bool operator<(const CJunc& a) { // sort by strand last
//         if (start==a.start){
//             if(end==a.end){
//                 return strand<a.strand;
//             }
//             else{
//                 return (end<a.end);
//             }
//         }
//         else{
//             return (start<a.start);
//         }
//     }

// 	void add(CJunc& j) {
//        dupcount+=j.dupcount;
// 	}

// 	void write(FILE* f, const char* chr) {
// 		juncCount++;
// 		fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%ld\t%c\n",
// 				chr, start-1, end, juncCount, (long)dupcount, strand);
// 	}
// };

void processOptions(int argc, char* argv[]);

int main(int argc, char*argv[]) {
    processOptions(argc, argv);

    GSamReader samreader();
}

void processOptions(int argc, char* argv[]) {

}