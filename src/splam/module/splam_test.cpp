/*  splam_test.cpp -- 

    Copyright (C) 2023 Kuan-Hao Chao

    Author: Kuan-Hao Chao <kuanhao.chao@gmail.com> */

#include <thread>
#include <vector>
#include <htslib/htslib/sam.h>
#include <htslib/htslib/hts.h>

#include "bam_sort.h"


int array_size(char *test[]) {
  int size = 0;
  while (test[size] != NULL) {
    fprintf(stderr, "s_multi_argv (%d): %s\n", size, test[size]);
    size++;
  }
  return size;
  // 'size' now contains the size of the array
}

int main(int argc, char *argv[]) {

    for (int i=0; i<argc; i++) {
        fprintf(stderr, "%s ", argv[i]);
    }

    char *s_multi_argv[] = {"samtools", "sort", "-@", "1", "SRR1352129_chr9_sub_p/tmp/s_multi_cleaned_nh_updated.bam", "-o", "SRR1352129_chr9_sub_p/tmp/s_multi_cleaned_nh_updated.sort.bam"};
    int argc_sort = sizeof(s_multi_argv) / sizeof(s_multi_argv[0]);

    // int argc_sort = array_size(s_multi_argv);

    fprintf(stderr, "argc_sort count: %d\n", argc_sort);

    // res = bam_sort(argc_sort-1, s_multi_argv+1);

    // fprintf(stderr, ">> bam_sort res: %d\n\n\n", res);
    // char *test[] = {"Hello", "World", "Example"};


    // bam_sort(10, test);
    // bam_merge(10, test);

    return 0;
}
