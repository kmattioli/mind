
# coding: utf-8

# In[22]:


import warnings
warnings.filterwarnings('ignore')

import argparse
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.patheffects
import matplotlib.pyplot as plt
import multiprocessing
import math
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys

from matplotlib.font_manager import FontProperties
from matplotlib import gridspec
from matplotlib import transforms
from multiprocessing import Pool
from statsmodels.sandbox.stats import multicomp

sys.path.append("utils")
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()},
                                    reload_support=True)
import cython_fnx
from utils import *
from plotting_utils import *


# ## 1. define & collect arguments

# In[23]:


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", type=int, required=False, default=12345, 
                    help="numpy.random seed to use (for reproducibility)")
parser.add_argument("-p", "--pwm_file", type=str, required=True, 
                    help="path to pwm file in MEME format")
parser.add_argument("-d", "--deletion_info_file", type=str, required=True, 
                    help="path to file containing list of deletion files to analyze (full paths)")
parser.add_argument("-l", "--seq_len", type=int, required=True, 
                    help="length of sequences used in deletion MPRA")
parser.add_argument("-f", "--offset", type=int, required=False, default=1,
                    help="# bp that the start number in deletion files is offset from 0")
parser.add_argument("-w", "--bandwidth", type=int, required=False, default=5,
                    help="bandwidth to use in moving average smoother (NOTE: should be odd)")
parser.add_argument("-p", "--peak_cutoff", type=float, required=False, default=0.5,
                    help="cutoff to use when calling peaks")
parser.add_argument("-t", "--score_type", type=str, required=True, 
                    help="either 'loss' or 'gain'")
parser.add_argument("-n", "--n_shuffles", type=int, required=False, default=1000,
                    help="# times to shuffle peak data to get null distribution")
parser.add_argument("-e", "--tfs_expressed_file", type=str, required=False, default=None, 
                    help="path to file containing list of TFs expressed in cell line of interest")
parser.add_argument("-c", "--cores", type=int, required=True,
                    help="# cores to use when computing")
parser.add_argument("-o", "--out_dir", type=str, required=True, 
                    help="directory where results will be stored")


# In[24]:


args = parser.parse_args()
seed = args.seed
pwm_file = args.pwm_file
deletion_info_file = args.deletion_info_file
seq_len = args.seq_len
offset = args.offset
bandwidth = args.bandwidth
peak_cutoff = args.peak_cutoff
score_type = args.score_type
n_shuffles = args.n_shuffles
tfs_expressed_file = args.tfs_expressed_file
cores = args.cores
out_dir = args.out_dir


# In[25]:


# defaults for debugging
# seed = 12345
# pwm_file = "inputs/0__pwm/pfm_vertebrates_meme_motifNameChanged.txt"
# deletion_info_file = "inputs/1__dels/deletion_files.txt"
# seq_len = 94
# offset = 1
# bandwidth = 5
# peak_cutoff = 0.5
# score_type = "loss"
# n_shuffles = 1000
# tfs_expressed_file = None
# cores = 4
# out_dir = "results/test"


# In[26]:


### argument assertions ###

# pwm file exists
assert os.path.exists(pwm_file), "--pwm_file path does not exist"

# deletion file exists
assert os.path.exists(deletion_info_file), "--deletion_info_file path does not exist"

# bandwidth is odd
assert bandwidth % 2 == 1, "--bandwidth should be an odd number"

# score is either loss or gain
assert score_type in ["loss", "gain"], "--score_type should be either 'loss' or 'gain'"

# tf file exists if given
if tfs_expressed_file != None:
    assert os.path.exists(tfs_expressed_file), "--tfs_expressed_file path does not exist"


# In[27]:


### set plotting defaults ###
sns.set(**PRESET)
fontsize = FONTSIZE


# ## 2. import data & make out dir if needed

# In[28]:


# set seed for reproducibility!
np.random.seed(seed)


# In[29]:


# read in pwm file
motifs, motif_lens, motif_len_map = parse_pfm(pwm_file, False)


# In[30]:


# find max motif length
max_motif_len = np.max(list(motif_lens))


# In[31]:


# read in file with paths to all of the deletion data
deletion_info = pd.read_table(deletion_info_file, sep="\t", header=None)
deletion_info.columns = ["path", "name"]
deletion_info["path"] = deletion_info["path"].map(str.strip)
deletion_info["name"] = deletion_info["name"].map(str.strip)
deletion_info = zip(list(deletion_info["path"]), list(deletion_info["name"]))


# In[32]:


# read in all of the deletion data and make sure it has the columns we need
data = {}
for path, name in deletion_info:
    
    assert os.path.exists(path), "path to deletion data %s does not exist" % path
    df = pd.read_table(path, sep="\t")
    
    assert "delpos" in df.columns, "deletion file %s does not have 'delpos' as a column name" % path
    assert "seq" in df.columns, "deletion file %s does not have 'seq' as a column name" % path
    assert "mean.log2FC" in df.columns, "deletion file %s does not have 'mean.log2FC' as a column name" % path
    assert "sd" in df.columns, "deletion file %s does not have 'sd' as a column name" % path
    assert "se" in df.columns, "deletion file %s does not have 'se' as a column name" % path
    
    data[name] = df


# In[33]:


# read in the tfs expressed file, if it exists
if tfs_expressed_file != None:
    tfs_expressed = pd.read_table(tfs_expressed_file, sep="\t", header=None)
    tfs_expressed.columns = ["tf"]
    tfs_expressed = list(tfs_expressed["tf"])
else:
    tfs_expressed = None


# In[34]:


# make subdir for results files
res_dir = "%s/files" % out_dir
if not os.path.exists(res_dir):
    os.makedirs(res_dir)
    
# make subdir for figures
figs_dir = "%s/figs" % out_dir
if not os.path.exists(figs_dir):
    os.makedirs(figs_dir)


# ## 2. create loss or gain score
# loss = looking for transcriptional activators; gain = looking for transcriptional repressors

# In[35]:


for seq in data.keys():
    df = data[seq]

    if score_type == "loss":
        df["loss_score_raw"] = df.apply(loss_score, col="mean.log2FC", axis=1)
        
        # scale scores (only if there are scores > 2, otherwise don't)
        scores = list(df["loss_score_raw"])
        scores.extend([0, 2])
        scaled_scores = scale_range(scores, 0, 2)
        del scaled_scores[-2:]
        df["loss_score_raw_scaled"] = scaled_scores
        
    elif score_type == "gain":
        df["gain_score_raw"] = df.apply(gain_score, col="mean.log2FC", axis=1)
        
        # scale scores (only if there are scores > 2, otherwise don't)
        scores = list(df["gain_score_raw"])
        scores.extend([0, 2])
        scaled_scores = scale_range(scores, 0, 2)
        del scaled_scores[-2:]
        df["gain_score_raw_scaled"] = scaled_scores


# ## 3. find peaks in the sequences and write files w/ peak info

# In[36]:


# make subdir for peak figures
peak_figs_dir = "%s/0__peaks" % figs_dir
if not os.path.exists(peak_figs_dir):
    os.makedirs(peak_figs_dir)


# In[37]:


# make subdir for peak results
peak_res_dir = "%s/0__ntd_scores" % res_dir
if not os.path.exists(peak_res_dir):
    os.makedirs(peak_res_dir)


# In[38]:


score_col = "%s_score_raw_scaled" % (score_type)
data_peaks = {}
peak_dfs = {}

for seq in data.keys():
    print("plotting: %s" % seq)
    seq_name = "%s__%s" % (seq, score_type)

    # extract bases & scores from df
    df = data[seq]
    bases = list(df["seq"])
    scaled_scores = list(df[score_col])
    yerrs = list(df["se"])
    raw_scores = list(df["mean.log2FC"])
    
    # apply moving average to scores
    scores_filt = moving_average(bandwidth, scaled_scores)
    df["filtered_score"] = scores_filt
    
    # find peaks
    widths, peak_info, df = find_peaks(peak_cutoff, seq_len, df, 
                                       scores_filt, scaled_scores, bases, offset, 
                                       max_motif_len)
    data_peaks[seq] = peak_info
    peak_dfs[seq] = df
    
    # write file
    df.to_csv("%s/%s.%s_score.txt" % (peak_res_dir, seq, score_type), sep="\t", index=False)
    
    # plot peaks
    plot_peaks((4.9, 1.4), 6, score_type, seq_len, seq_name, bandwidth, widths, raw_scores, yerrs, scores_filt, 
               scaled_scores, bases, peak_figs_dir)


# ## 4. find motifs in the peaks

# In[39]:


print("")
print("mapping motifs...")
results_dict = get_all_results(cores, data_peaks, motifs, n_shuffles, seed, parallel=True)
print("")


# ## 5. correct results for mult. hyp. & plot

# In[40]:


# make subdir for results figures
res_figs_dir = "%s/1__results" % figs_dir
if not os.path.exists(res_figs_dir):
    os.makedirs(res_figs_dir)


# In[41]:


# make subdir for results figures
motif_res_dir = "%s/1__motif_scores" % res_dir
if not os.path.exists(motif_res_dir):
    os.makedirs(motif_res_dir)


# In[42]:


# check motifs at 3 FDRs (since every peak is different)
alphas = [0.05, 0.1, 0.15]

for seq in results_dict:
    seq_results = results_dict[seq]
    n_peaks = len(seq_results)
    for p in range(n_peaks):
        name = "%s__peak%s" % (seq, p+1)
        print(name)
        peak_results = seq_results[p]
        deduped_dict = {}

        # for every motif, choose the max of either sense or antisense pwm
        for motif in peak_results:
            sense_score = peak_results[motif][0][2]
            antisense_score = peak_results[motif][1][2]
            if sense_score >= antisense_score:
                max_row = peak_results[motif][0]
                max_row.extend(["sense"])
            else:
                max_row = peak_results[motif][1]
                max_row.extend(["antisense"])
            max_row_fixed = list(max_row[0:8])
            deduped_dict[motif] = max_row_fixed
            
        df = pd.DataFrame.from_dict(deduped_dict, orient="index").reset_index()
        df.columns = ["motif", "start", "end", "score", "pval", "tile_chr", "tile_start", 
                      "tile_end", "strand"]
        
        # correct p values for multiple testing
        padj = multicomp.multipletests(list(df["pval"]), method="fdr_bh")[1]
        df["padj"] = padj
        df["neg_log_pval"] = -np.log10(df["padj"])
        
        for alpha in alphas:
            sig_df = df[df["padj"] < alpha]
            if len(sig_df) > 0:
                break
        print("ALL MOTIFS: found %s motifs at %s FDR" % (len(sig_df), alpha))
        df["fdr_cutoff"] = alpha
        df = df.sort_values(by="score", ascending=False)
        df.to_csv("%s/%s.%s__results.txt" % (res_dir, name, score_type), sep="\t", index=False)

        
        # if a list of TFs expressed was provided, filter to those only & adjust those only
        if tfs_expressed != None:
            df_sub = df[df["motif"].isin(tfs_on)]
        
            padj = multicomp.multipletests(list(df_sub["pval"]), method="fdr_bh")[1]
            df_sub["padj"] = padj
            df_sub["neg_log_pval"] = -np.log10(df_sub["padj"])
        
        
            for alpha in alphas:
                sig_df_sub = df_sub[df_sub["padj"] < alpha]
                if len(sig_df_sub) > 0:
                    break
            print("ONLY TFS EXPRESSED: found %s motifs at %s FDR" % (len(sig_df), alpha))
            df_sub["fdr_cutoff"] = alpha
            df_sub = df_sub.sort_values(by="score", ascending=False)
            df_sub.to_csv("%s/%s.%s__results.expr_filt.txt" % (res_dir, name, score_type), sep="\t", 
                          index=False)
        
        # plot the all results only
        plot_motif_results(df, 2.5, name, alpha, res_figs_dir)
        


# In[ ]:




