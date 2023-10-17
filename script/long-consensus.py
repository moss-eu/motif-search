import numpy as np
import pandas as pd
import sys
import csv

# nb_motif = 20
# nb_oligo = 2751
# nb_oligo_bi = 101010111111 #binary of 2751
# # search_path = "/Users/yan/Desktop/eurecom_code/motif-search/sample_data/EMP/read-medium-20x.res"
# # out_path = '/Users/yan/Desktop/eurecom_code/motif-search/sample_data/EMP/read-medium-20x.csv'
# # ref_path = "/Users/yan/Desktop/eurecom_code/motif-search/sample_data/EMP/ref-motif-only"
# search_path = "/media/ssd/ngs-data-analysis/2022_motif/mc_testbedenv_05/res/read-medium-20x.res"
# out_path = '/media/ssd/ngs-data-analysis/2022_motif/mc_testbedenv_05/res/read-medium-20x.csv'
# ref_path = "/media/ssd/ngs-data-analysis/2022_motif/mc_testbedenv_05/ref/ref-motif-only"

nb_motif = int(sys.argv[1])
nb_oligo = int(sys.argv[2])
nb_oligo_bi = int(sys.argv[3])
search_path = str(sys.argv[4])
out_path = str(sys.argv[5])
if len(sys.argv) == 6:
    ref_path = ''
else:
    ref_path = str(sys.argv[6])

cnt_threshold = 10

import collections

def get_most_common(m_list):
    length = len(m_list)
    haha = []
    i = 0
    while i < length:
        j = i
        while (j < length) and (m_list[j][0] == m_list[i][0]) and (m_list[j][1] == m_list[i][1]):
            j+=1
        k = 2
        res = m_list[i][0] + '-' + m_list[i][1]
        while k < nb_motif:
            m = i
            motifs = []
            while m < j:
                motifs.append(m_list[m][k])
                m+=1
            res += '-' + collections.Counter(motifs).most_common()[0][0] 
            k+=1
        haha.append(res)
        i = j
    return haha


# In[3]:


colnames = ['read_id', 'FLAG', 'POS', 'DIST', 'motifs']

search = pd.read_csv(search_path, sep = "\t", header = None,  names = colnames, quoting=csv.QUOTE_NONE)
search['ref_id'] = search['read_id'].str.split(',').str[0].str.split(' ').str[1]
search['read_id'] = search['read_id'].str.split(' ').str[0].str[2: ]
search['motifs'] = search['motifs'].str[1: ]
search['m0'] = search['motifs'].str.split('-').str[0]
search['m1'] = search['motifs'].str.split('-').str[1]
print("nb of reads: " + str(len(search)))

good_start = search[search['m0'] == '00000000000000000']
print("nb of reads start with 00000000000000000: " + str(len(good_start)))

good_start = good_start[good_start['m1'].str[1:].astype(int) < nb_oligo_bi]
print("nb of reads valid index: " + str(len(good_start)))

good_start = good_start[good_start['motifs'].str.split('->').str.len()==nb_motif]
print("nb of reads has " + str(nb_motif) + " motifs: " + str(len(good_start)))

bad_start = search[ ~search.isin(good_start)].dropna()
print("nb of reads not valid index: " + str(len(bad_start)))


# In[4]:


good_start = good_start.sort_values(by=["motifs"])

m_list = good_start['motifs'].str.split('->',expand=True).values.tolist()

good_start_consensus = get_most_common(m_list)

good_start_cnt = good_start.groupby(['m0', 'm1']).size().reset_index()


# In[5]:


bad_m = bad_start['motifs'].str.split('->',expand=True).values.tolist()
good_m = [s.split('-') for s in good_start_consensus]


bad_m_updated = []
for bm in bad_m:
    best_nb_match = 0
    for gm in good_m:
        nb_match = 0
        for i in range(nb_motif):
            if bm[i] == gm[i]:
                nb_match+= 1
        if nb_match > best_nb_match:
            best_nb_match = nb_match
            bm[0] = gm[0]
            bm[1] = gm[1]
    bad_m_updated.append('->'.join(bm))


# In[6]:


## filter the valid m0-m1 of rescued oligos
bad_m_df = pd.DataFrame(bad_m_updated)
bad_m_df['m0'] = bad_m_df[0].str.split('-').str[0]
bad_m_df['m1'] = bad_m_df[0].str.split('-').str[1]


# In[17]:

invalid_cnt = good_start_cnt[good_start_cnt[0] < cnt_threshold][['m0', 'm1']]
bad_m_df_filerted = pd.merge(bad_m_df, invalid_cnt, how='inner', left_on=['m0', 'm1'], right_on=['m0', 'm1'])


## merge the good reads and rescued bad reads
all_oligo = good_start['motifs'].values.tolist()
all_oligo.extend(bad_m_df_filerted[0].values.tolist())
all_oligo.sort()

all_oligo_list = [s.split('->') for s in all_oligo]


# In[23]:
consensus = get_most_common(all_oligo_list)

res = pd.DataFrame(consensus, columns = ['motifs'])
res['key_0'] = res['motifs'].str.split('-').str[0]
res['key_1'] = res['motifs'].str.split('-').str[1]
res['motifs'].to_csv(out_path, header=False, index=False)

if ref_path != '':
    ref = pd.read_csv(ref_path, sep=" ", usecols=[0, 1], names=['id', 'motifs'], header=None, quoting=csv.QUOTE_NONE)
    ref['motifs'] = ref['motifs'].str[1:]
    ref['key_0'] = ref['motifs'].str.split('>').str[0]
    ref['key_1'] = ref['motifs'].str.split('>').str[1]
    ref['motifs'] = ref['motifs'].str.split('>')
    ref['motifs'] = ref['motifs'].str.join(">")
    ref['id'] = ref['id'].str[1:]

    merged = pd.merge(ref, res, how='inner', left_on=['key_0', 'key_1'], right_on=['key_0', 'key_1'])
    motifs_x_split = merged['motifs_x'].str.split('>', expand=True)
    motifs_y_split = merged['motifs_y'].str.split('-', expand=True)

    merged['matched'] = (motifs_x_split.T == motifs_y_split.T).astype(int).sum()
    merged['num_motif_y'] = motifs_y_split.count(axis=1)

    total = merged['num_motif_y'].sum()
    match = merged['matched'].sum()
    full_recover_read = merged[merged['matched'] == 20]

    print("nb of reference oligos: " + str(len(ref)))
    print("nb of inferred oligos: " + str(len(res)))
    print("nb of fully recovered oligos: " + str(len(full_recover_read)))
    print("recovered oligo (%): " + str(len(full_recover_read) / len(ref)))
    print("recovered motifs (%): " + str(match / total))
