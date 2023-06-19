#include <thread>
#include <vector>
#include <htslib/htslib/sam.h>
#include <htslib/htslib/hts.h>

#include "bam_sort.h"

int main() {

    char *test[] = {"Hello", "World", "Example"};


    bam_sort(10, test);
    // bam_merge(10, test);

    return 0;
}
