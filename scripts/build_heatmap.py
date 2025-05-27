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

daysIncluded=60
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
lineages = ['BA.2.86','JN.1','BA.2.87.1','XBB.1.5']
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

# import plotly.graph_objects as go


# fig = px.imshow(df,x=df.columns,y=df.index,zmin=0,zmax=1.,labels=dict(x="Mutation",y="Sample",color="Frequency"),color_continuous_scale='YlOrRd')
# fig.update_layout(showlegend=False)
# fig.update_xaxes(visible=False)
# fig.update_yaxes(visible=False)
# fig.update_coloraxes(colorbar={'orientation':'h', 'thickness':20, 'y': -0.2,'len':0.4})
# fig.write_html("test_heatmap.html")

# convert to AA mutations. 

import freyja.read_analysis_utils as rutils
gff = rutils.parse_gff('NC_045512_Hu-1.gff')
AA_muts = rutils.translate_snvs(df.columns,'NC_045512_Hu-1.fasta',gff)


from plotly.subplots import make_subplots
import plotly.graph_objects as go
fullLen = df.shape[0]
df2 = df.loc[lineages]
df = df.iloc[0:df.shape[0]-len(lineages)]

tt = [times_df.loc[dfi,'SampleCollectionDate'].strftime('%Y-%m-%d') for dfi in df.index]

df = pd.concat((df,pd.Series(tt,name='SampleCollectionDate',index=df.index)),axis=1)
df = df.sort_values(by='SampleCollectionDate')
df = df.drop(columns='SampleCollectionDate')


#resorted version. 
tt = [times_df.loc[dfi,'SampleCollectionDate'].strftime('%Y-%m-%d') for dfi in df.index]
df_dates = pd.DataFrame(index=df.index,columns=df.columns)
for j0,dfi in enumerate(df.index):
    df_dates.loc[dfi,:] = tt[j0]

df_AA = pd.DataFrame(index=df.index,columns=df.columns)
for j0,dfc in enumerate(df.columns):
    df_AA.loc[:,dfc] = AA_muts[dfc]
df_AA = df_AA.fillna('')


tt = [times_df.loc[dfi,'SiteProvince'] for dfi in df.index]
df_Province = pd.DataFrame(index=df.index,columns=df.columns)
for j0,dfi in enumerate(df.index):
    df_Province.loc[dfi,:] = tt[j0]

df_custom = np.stack( (~df.isna(),df_dates,df_Province,df_AA), axis=-1)
df2_custom = ~df2.isna()

df.to_feather('mutations_wastewater.feather')

with open('mutations_wastewater_meta.npy', 'wb') as f:
    np.save(f, df_custom)

with open('mutations_wastewater_meta.npy', 'rb') as f:
    df00 = np.load(f,allow_pickle=True)

fig = make_subplots(1, 1, vertical_spacing=0.05,subplot_titles=["Mutations in Wastewater"],row_heights=[df.shape[0]/fullLen])#, (df2.shape[0]-1)/fullLen])
fig.add_trace(go.Heatmap(z=df,x=df.columns,y=df.index,zmin=0,zmax=1.,colorscale='YlOrRd', customdata = df_custom,
                         colorbar={'title':'SNV Frequency', 'orientation':'h', 'thickness':20, 'y': -0.1,'len':0.4},
                         showscale=True, hovertemplate='<b>Nuc. Mutation:%{x}</b><br><b>AA Mutation:%{customdata[3]}<br><b>Sample:%{y}<br><b>Seq. Covered:%{customdata[0]}<br><b>Frequency:%{z}<br><b>Date:%{customdata[1]}<br><b>Province:%{customdata[2]}',
                         name=""), 1, 1)
# # fig.update_layout(yaxis = dict(scaleanchor = 'x'))
# fig.add_trace(go.Heatmap(z=df2,x=df2.columns,y=df2.index,zmin=0,zmax=1.,colorscale='YlOrRd',showscale=False, customdata = df2_custom,
#                          hovertemplate='<b>Mutation:%{x}</b><br> <b>Lineage:%{y}<br> <b>Present:%{z}',
#                          name=""), 2, 1)
# fig.update_traces(yaxis = dict(scaleanchor = 'x'))
fig.update_layout(showlegend=False)
fig.update_xaxes(visible=False)
fig.update_yaxes(visible=False)

# fig['layout']['yaxis2'].update(scaleanchor = 'x')
# fig.update_coloraxes(colorbar={'orientation':'h', 'thickness':20, 'y': -0.2,'len':0.4})
fig.write_html("../SNV_heatmap.html")

