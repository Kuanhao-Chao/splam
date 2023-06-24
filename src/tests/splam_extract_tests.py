import pandas as pd
import os 
import argparse
import sys
from splam import prediction, config, parse, chr_size
import splam_extract
def main(argv=None):
    print("Hello world!")
    splam_extract.splam_extract(["splam-extract", "-o", "SRR1352129_chr9_sub", "--paired", "../../Dataset/SRR1352129_chr9_sub/SRR1352129_chr9_sub.bam"])
    
if __name__ == '__main__':
    main()