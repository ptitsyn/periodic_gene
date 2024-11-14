import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

import sys
import os
sys.path.append(os.path.abspath("C:\\Andrey\\perpy"))
from per import *

#import scipy as sp

df=pd.read_csv('liver_count.csv', index_col='gene_id', nrows=1000)
num_genes=df.shape[0]
data_cols=df.columns

#for consistency, both period and framelength are numbers of timepoints, although the framelength must be a multiple of
#periods and presumably a multiple of single gene profiles, such as 3 genes at a time, 30 points per gene: framelength=90
tested_df=phase_continuum_test(df, period=6, framelength=60, cutoff=0.05, method=1).drop('signal/noise', axis=1)
#let's use the same dataframe to output all test results. For this, make new test data keepint the same order of raws as in phase continuum results
new_df=tested_df[data_cols]
#allocate arrays
fg=np.zeros(num_genes)
pt=np.zeros(num_genes)
ac=np.zeros(num_genes)
i=0
for gene in tqdm(df.index):
    fg[i] = Fg_test(df.loc[gene],6)
    pt[i] = Pt_test(df.loc[gene],6)
    phase, ac[i] = autocorr(df.loc[gene],6)
#add test results columns to the dataframe
tested_df['Fishers g-test']=fg
tested_df['Pt-test']=pt
tested_df['Autocorr test']=ac
tested_df.to_csv('all_test_results.csv')
print(tested_df)
exit(0)

