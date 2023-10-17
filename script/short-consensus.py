import numpy as np
import pandas as pd
import sys
import csv

min_as = 0
res_path = str(sys.argv[1])
out_path = str(sys.argv[2])
if len(sys.argv) == 3:
    ref_path = ''
else:
    ref_path = str(sys.argv[3])

def valid_dist(d):
    for dd in d:
        if int(dd) < min_as:
            return False
    return True

## load motif and ref
colnames = ['RNAME', 'FLAG', 'POS', 'DIST', 'SEQ']
raw_motif = pd.read_csv(res_path, sep='\t', header=None, names=colnames, quoting=csv.QUOTE_NONE, engine='python')

# dist filter
raw_motif['dist'] = raw_motif['DIST'].str.split('>').str[:-1]
raw_motif['valid_dist'] = raw_motif['dist'].apply(valid_dist)
infered_motif = raw_motif[raw_motif['valid_dist']]

# cnt filter
cnt = infered_motif.groupby(['SEQ']).size().reset_index(name='counts')
cnt = cnt[cnt['counts'] > 1]

valid_motif = pd.merge(cnt, infered_motif, how='inner', left_on=['SEQ'], right_on=['SEQ'])
valid_motif['SEQ'].to_csv(out_path, index=False, header=None)

if ref_path != '':
    ref = pd.read_csv(ref_path, header=None, sep=' ', names=['ref-id', 'ref-motif'], quoting=csv.QUOTE_NONE)
    merged = pd.merge(ref, valid_motif, how='outer', left_on=['ref-motif'], right_on=['SEQ'])

    TP = merged[merged['ref-motif'].notna() & merged['SEQ'].notna()]
    FP = merged[merged['ref-motif'].isna() & merged['SEQ'].notna()]
    FN = merged[merged['ref-motif'].notna() & merged['SEQ'].isna()]

    TP_unq = TP['SEQ'].nunique()
    FP_unq = FP['SEQ'].nunique()
    FN_unq = FN['ref-motif'].nunique()

    print("yyan-log:----- for dataset: " + res_path + "-----")
    print("yyan-key-log nb of motifs detect in total: " + str(len(raw_motif)))
    # print("yyan: nb of motifs detect, AS is bigger than " + str(min_as) + " (TP+FP): " + str(len(infered_motif)))
    print('yyan-key-log TP unique oligos: ' + str(TP_unq))
    print('yyan-key-log FP unique oligos: ' + str(FP_unq))
    print('yyan-key-log FN unique oligos: ' + str(FN_unq) + '\n')
