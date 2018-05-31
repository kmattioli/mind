import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

from matplotlib.font_manager import FontProperties
from matplotlib import gridspec
from matplotlib import transforms

sys.path.append("utils")
from utils import *

##### global variables #####
COLOR_DICT = {"A": "crimson", "C": "mediumblue", "G": "orange", "T": "forestgreen"}


PRESET = {"style": "white", "font": "Helvetica", "context": "paper", 
          "rc": {"font.size":6.5,"axes.titlesize":6.5,
                  "axes.labelsize":6.5, 'axes.linewidth':0.5,
                  "legend.fontsize":6.5, "xtick.labelsize":6.5,
                  "ytick.labelsize": 6.5}}
FONTSIZE = 6.5


class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def plot_peaks(figsize, fontsize, score_type, seq_len, seq_name, window_len, widths, scores, yerrs, scores_filt, scaled_scores, bases, figs_dir):
    
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1], hspace=0)
    peak_ax = plt.subplot(gs[0])
    seq_ax = plt.subplot(gs[1])
    
    ### peaks figure ###
    # plot peak locations
    for w in widths:
        peak_ax.axvline(x=w[0], color="gray", linestyle="solid", linewidth=0.5, zorder=1)
        peak_ax.axvline(x=w[1], color="gray", linestyle="solid", linewidth=0.5, zorder=1)
        peak_ax.axvspan(w[0], w[1], alpha=0.5, color="gainsboro", zorder=1)
    
    # plot deletion values
    xs = list(range(0, seq_len))
    peak_ax.bar(xs, scores, yerr=yerrs, color="lightgray", edgecolor="gray", linewidth=0.5, ecolor="gray", 
                error_kw={"elinewidth": 0.75})
    if score_type == "loss": 
        scores_filt = [-1*x for x in scores_filt]
    else:
        scores_filt = scores_filt
    peak_ax.plot(scores_filt, linestyle="dashed", color="black", zorder=2, label="smoothed data")
    
    # labels
    peak_ax.set_xlim((-0.75, seq_len))
    peak_ax.set_ylabel("log2(del/WT)")
    peak_ax.xaxis.set_visible(False)
    peak_ax.set_title("filtered scores and peaks (window length: %s): %s" % (window_len, seq_name))
    
    ### seq logo ###
    mpl.rcParams["font.family"] = "Arial"
    scaled_scores = scale_range(scaled_scores, 0.5, 2.0)
    
    font = FontProperties()
    font.set_size(fontsize)
    font.set_weight("bold")
    
    seq_ax.set_xticks(range(1,len(scaled_scores)+1))
    seq_ax.set_ylim((0, 2))
    seq_ax.axis("off")
    trans_offset = transforms.offset_copy(seq_ax.transData, 
                                          fig=fig, 
                                          x=1, 
                                          y=0, 
                                          units="dots")
    
    for i in range(0, len(scaled_scores)):
        score = scaled_scores[i]
        base = bases[i]
        color = COLOR_DICT[base]
        txt = seq_ax.text(i+0.5, 0, base, transform=trans_offset,fontsize=fontsize, color=color, 
                          ha="center", fontproperties=font)
        txt.set_path_effects([Scale(1.0, score)])
        fig.canvas.draw()
        trans_offset = transforms.offset_copy(seq_ax.transData, fig=fig, x=1, y=0, units='points')
    
    plt.tight_layout()
    plotname = "%s__scores_and_peaks_w%s.pdf" % (seq_name, window_len)
    fig.savefig("%s/%s" % (figs_dir, plotname), dpi="figure", bbox_inches="tight")
    plt.close()


def plot_motif_results(results, size, seq_name, fdr, figs_dir):
    xlims = (np.min(results["score"]-0.1), np.max(results["score"]+2))
    ylims = (np.min(results["neg_log_pval"]-0.1), np.max(results["neg_log_pval"]+0.2))
    sig_results = results[results["padj"] < fdr]
    non_sig_results = results[results["padj"] >= fdr]
    g = sns.JointGrid("score", "neg_log_pval", results, size=size, space=0.02)

    # score distplots
    sns.distplot(results["score"], kde=False, bins=25, color="gray", ax=g.ax_marg_x)
    g.ax_marg_x.set_xlim(xlims)
    g.ax_marg_x.set_xlabel("")

    # pval distplots
    sns.distplot(results["neg_log_pval"], kde=False, bins=25, color="gray", ax=g.ax_marg_y, vertical=True)
    g.ax_marg_y.set_ylim(ylims)
    g.ax_marg_y.set_ylabel("")

    # scatter plot
    g.ax_joint.plot(non_sig_results["score"], non_sig_results["neg_log_pval"], "o", 
                    ms=5, label="", alpha=0.5, color="gray")
    g.ax_joint.plot(sig_results["score"], sig_results["neg_log_pval"], "o", 
                    ms=5, label="FDR < %s\n[n = %s]" % (fdr, len(sig_results)), alpha=0.75, 
                    color=sns.color_palette()[1])
    g.ax_joint.set_xlim(xlims)
    g.ax_joint.set_ylim(ylims)
    g.ax_joint.set_xlabel("motif score")
    g.ax_joint.set_ylabel("-log10(adjusted pval)")
    
    # make legend
    g.ax_joint.legend(loc=4, frameon=True)
    
    plt.tight_layout()
    g.savefig("%s/%s__motif_results.pdf" % (figs_dir, seq_name), bbox_inches="tight", dpi="figure")
    plt.close()
