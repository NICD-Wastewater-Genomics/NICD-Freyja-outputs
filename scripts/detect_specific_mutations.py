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
for j,fn in tqdm(enumerate(os.listdir(dir0))):
    sname = fn.replace('_','-').split('-S')[0].split('ENV-')[-1]
    if sname not in times_df.index:
        continue
    df = pd.read_csv(dir0+fn, sep='\t')

    #restrict to spike
    df = df[(df['POS']>=21563) & (df['POS']<=25384)]
    #require at least ten reads and greater that 2% prevalence. 
    df = df[(df['ALT_DP']>10) & (df['ALT_FREQ']>0.02)]
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

updateBarcodes = False
if updateBarcodes or ('usher_barcodes.csv' not in os.listdir('.')):
    os.system('freyja update --outdir .')

df_barcodes = pd.read_csv('usher_barcodes.csv',index_col=0)
lineages = ['BA.2.86','JN.1','BA.2.87.1','XBB.1.5','BA.3','BA.2.15']
for lin in lineages:
    target = df_barcodes.loc[lin]
    target = target[target>0]
    targetMuts = list(target.index)
    targetMuts.sort(key=sortFun)
    df = pd.concat([df,pd.DataFrame({t: (1 if (t in targetMuts) else np.nan if '-' in t or '+' in t else 0) for t in df.columns },index=[lin])],axis=0)

#.fillna(0)
muts= list(df.columns)
muts.sort(key=sortFun)
df = df[muts]
df = df.loc[:,df.iloc[0:df.shape[0]-len(lineages)].sum(axis=0)>0.1]#.fillna('No Sequencing Coverage')


import freyja.read_analysis_utils as rutils
gff = rutils.parse_gff('NC_045512_Hu-1.gff')
AA_muts = rutils.translate_snvs(df.columns,'NC_045512_Hu-1.fasta',gff)

AA_muts_r = {v:k for k,v in zip(AA_muts.keys(),AA_muts.values())}


# check for mutations of interest.
tMuts = ['S:P26L','S:N164K','S:S172F','T21967-TGTAATGATCCATTTTTGGGTGTTTATTACCACAAA(S:DEL136/147)','S:I326V','S:K795T','S:D1184E']

tMutsNUC = [AA_muts_r[t] for t in tMuts]
#require 30% of mutations to be present. 
df_target = df[(df[tMutsNUC]>0.001).mean(axis=1)>=0.3]#df[(df[tMutsNUC]>0.001).all(axis=1)]

df_target = df_target[df_target.columns[df_target.sum(axis=0)>0]]
df_target = df_target[df_target.columns[df_target.mean(axis=0)<0.95]]

df_target2 = df_target.copy()
df_target2.index = [dfi+'/'+times_df.loc[dfi,'SiteProvince']+'/'+times_df.loc[dfi,'DistrictName'] for dfi in df_target.index]
df_target_sub = df_target[tMutsNUC]
df_target.columns = [AA_muts[t].split('(')[1] if '(' in AA_muts[t] else AA_muts[t]  for t in df_target.columns]
df_target_sub.columns = [AA_muts[t].split('(')[1] if '(' in AA_muts[t] else AA_muts[t] for t in df_target_sub.columns]
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib

for metro in ['Ekurhuleni MM','Tshwane MM']:
    mut_freqs = {}
    for dfi in df_target.index:
        if metro == times_df.loc[dfi,'DistrictName']:
            if ('Hospital' in times_df.loc[dfi,'SiteName']) or ('Other' in times_df.loc[dfi,'SiteName']):
                continue
            mut_freqs[times_df.loc[dfi,'SampleCollectionDate']] = df_target_sub.loc[dfi].to_dict()
            print(dfi,times_df.loc[dfi,'SampleCollectionDate'],times_df.loc[dfi,'DistrictName'],times_df.loc[dfi,'SiteName'])
    df_freq = pd.DataFrame(mut_freqs).T
    df_freq = df_freq.sort_index()
    fig, ax = plt.subplots()
    for c in df_freq.columns:
        ax.plot(df_freq.index,df_freq[c],'-o',label=c)
    mean = df_freq.mean(axis=1)
    ax.plot(df_freq.index,mean,'--',color='black')
    ax.set_ylabel('SNV frequency')
    locator = mdates.MonthLocator(bymonthday=1)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    ax.legend(bbox_to_anchor=(0.5, 1.16), loc='upper center',ncols=4)
    fig.savefig(f'mut_trajectories_{metro}.png')
    plt.close()
df_target.to_csv('ba3_new_muts.csv')
df_target_sub.index = [dfi+'/'+times_df.loc[dfi,'SiteProvince']+'/'+times_df.loc[dfi,'DistrictName']  for dfi in df_target_sub.index]
df_target_sub.to_csv('ba3_new_muts_subset.csv')

# tMuts = ['S:T76P','S:R78M','S:M153T','S:A163V']
# tMutsNUC = [AA_muts_r[t] for t in tMuts]
# df_target = df[(df[tMutsNUC]>0.001).mean(axis=1)>=0.60]

# df_target = df_target[df_target.columns[df_target.sum(axis=0)>0]]
# df_target = df_target[df_target.columns[df_target.mean(axis=0)<0.95]]

# df_target.index = [dfi+'/'+times_df.loc[dfi,'SiteProvince']+'/'+times_df.loc[dfi,'DistrictName']  for dfi in df_target.index]
# df_target_sub = df_target[tMutsNUC]
# df_target.columns = [AA_muts[t] for t in df_target.columns]
# df_target_sub.columns = [AA_muts[t] for t in df_target_sub.columns]
# df_target.to_csv('ba215_new_muts.csv')
# df_target_sub.to_csv('ba215_new_muts_subset.csv')