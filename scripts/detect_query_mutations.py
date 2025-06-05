import pandas as pd
import os
import numpy as np
from tqdm import tqdm
import datetime
# get metadata


times_df = pd.read_csv('../sample_metadata.csv', skipinitialspace=True).dropna().reset_index(drop=True)
times_df = times_df[[tc for tc in times_df.columns if "Unnamed" not in tc]]

times_df['Sequence_ID'] = times_df['Sequence_ID'].apply(lambda x:x.replace('_','-').replace('ENV-','').split('.')[0])

times_df['LabNumber'] = times_df['LabNumber'].apply(lambda x:x.replace('ENV-',''))
times_df = times_df.drop_duplicates()
dupMeta = times_df.loc[times_df['LabNumber'].duplicated(keep='first'),'LabNumber'].to_list()
if len(dupMeta)>0:
    print('lab numbers are duplicated.')
    print(dupMeta)

times_df = times_df.set_index('LabNumber')
times_df['SampleCollectionDate'] = pd.to_datetime(times_df['SampleCollectionDate'])
times_df = times_df.loc[~times_df.duplicated(keep='first')]

daysIncluded=180
cut_date = pd.to_datetime(times_df['SampleCollectionDate'].max()-datetime.timedelta(days=daysIncluded))
times_df = times_df[times_df['SampleCollectionDate']>=cut_date]

def probOK(s):
    if '+' in s or '-' in s:
        if ((len(s)-1)%3) != 0:
            return False
        else:
            return True
    else:
        return True

dir0='../variants/'
depthDir='../depths/'
inLoop = False

query = ['T27395A', 'G27382C', 'C635T','T4579A']
querySites = [int(t[1:-1]) for t in query]
for j,fn in tqdm(enumerate(os.listdir(dir0))):
    sname = fn.replace('_','-').split('-S')[0].split('ENV-')[-1]
    if sname not in times_df.index:
        continue
    df = pd.read_csv(dir0+fn, sep='\t')

    df = df[df['POS'].isin(querySites)]
    # #restrict to spike
    # df = df[(df['POS']>=21563) & (df['POS']<=25384)]
    #require at least ten reads and greater that 2% prevalence. 
    df = df[(df['ALT_DP']>10) & (df['ALT_FREQ']>0.01)]
    # ignore likely seq errors (indels of length not divisible by 3)
    df = df[df['ALT'].apply(lambda x: probOK(x))]
    # only nonsynonymous muts for now. 
    if df.shape[0]==0:
        continue
    df = df[ df['REF_AA']!=df['ALT_AA']]
    df = df[['REF','POS','ALT','ALT_FREQ','REF_AA','POS_AA','ALT_AA']]
    # add metadata
    if sname not in times_df.index:
        print(sname,' not in')
        continue
    df['Province'] = times_df.loc[sname,'SiteProvince']
    df['District'] = times_df.loc[sname,'DistrictName']
    df['Date'] = times_df.loc[sname,'SampleCollectionDate']
    df['sample'] = sname

    df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']

    if not inLoop:
        inLoop=True
        all_df = df.copy()
    else:
        all_df = pd.concat([all_df,df],axis=0)

def sortFun(x):
    # sort based on nuc position, ignoring nuc identities
    if '+' in x:
        return int(x[1:(x.index('+'))])
    elif '-' in x:
        return int(x[1:(x.index('-'))])
    else:
        return int(x[1:(len(x)-1)])

# all_df = all_df[all_df['District']=='Mangaung MM']
import plotly.express as px
df = pd.pivot_table(all_df, values='ALT_FREQ', index='sample', columns=['mutName'], sort=False)

for j,fn in enumerate(os.listdir(depthDir)):
    sname = fn.replace('_','-').split('-S')[0].split('ENV-')[-1]
    if sname in df.index:
        df_depth = pd.read_csv(depthDir+fn, sep='\t', header=None, index_col=1)
        pos = [sortFun(dfi) for dfi in df.columns]
        for k, val in enumerate(df.loc[sname]):
            if np.isnan(val):
                if df_depth.loc[pos[k],3] >10:
                    df.loc[sname,df.columns[k]] = 0.

import freyja.read_analysis_utils as rutils
gff = rutils.parse_gff('NC_045512_Hu-1.gff')
AA_muts = rutils.translate_snvs(df.columns,'NC_045512_Hu-1.fasta',gff)

AA_muts_r = {v:k for k,v in zip(AA_muts.keys(),AA_muts.values())}

df = df[[c for c in query]]

df2 = df[(df>0.01).mean(axis=1)>=0.75]