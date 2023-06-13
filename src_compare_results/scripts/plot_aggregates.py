import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys
from progress.bar import Bar

def run_plotter(db, type):

    # obtain the aggregated data from all 5 versions of spliceai
    aggregate_data(db, type)

def aggregate_data(db, type):

    aggregate_score = pd.DataFrame(columns=['seqid', 'start', 'end', 'strand', ''])
    for version_num in range(1, 6):
        ifp = f'./output/comparison/full_data.v{version_num}.{type}.{db}.csv'
        ofp = f'./output/aggregate/full_data.v{version_num}.{type}.{db}.csv'
        os.makedirs(os.path.dirname(ofp), exist_ok=True)

        full_data = pd.read_csv(ifp)


        # open the


if __name__ == '__main__':

    type = 'noN'
    databases = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    nums = [3]

    for num in nums:
        run_plotter(databases[num], type)