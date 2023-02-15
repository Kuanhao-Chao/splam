#include <iostream>
#include <fstream>
// #include <vector>
// #include <string>
#include "tmerge.h"

struct CRead {
    char* read;
    int pair, spurCnt;
    CRead(int vp=0, const char* vr= (char*) "ERR", int vsc=1):
            pair(vp), spurCnt(vsc){
        read = (char *) malloc(strlen(vr) + 1);
        strcpy(read, vr);
    }

    // overload operators
    bool operator==(const CRead& a) {
        return (pair == a.pair && !strcmp(read, a.read));
    }

    bool operator<(const CRead& a) {
        if (!strcmp(read, a.read)) {
            return (pair < a.pair);
        } else {
            return (strcmp(read, a.read) < 0);
        }
    }

    void incSpurCnt() {
        spurCnt++;
    }

    void clear() {
        free(read);
    }
};

struct PBRec {
    GSamRecord* r;
    PBRec(GSamRecord *rec=NULL):
    r(rec){ }
};
