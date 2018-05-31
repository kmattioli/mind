import multiprocessing
import math
import numpy as np
import time
import sys

from multiprocessing import Pool

sys.path.append("utils")
import pyximport; pyximport.install()
import cython_fnx

def get_peak_results(peak_start, tile_locs, mat, motifs, n_shuf, seed):
    peak_results = {}
    for motif in motifs:
        pwm_sense = motifs[motif][0]
        pwm_anti = motifs[motif][1]
        pwm_results = []
        for pwm in [pwm_sense, pwm_anti]:
            info_content = np.sum(np.sum(pwm_to_bits(pwm)))
            results = cython_fnx.get_motif_results(pwm, info_content, mat, peak_start, n_shuf, seed)
            results.extend(tile_locs)
            pwm_results.append(results)
        peak_results[motif] = pwm_results
    return peak_results

def get_all_results(processes, data_peaks, motifs, n_shuf, seed, parallel=True):
    if parallel:
        print("using multiprocessing")
        t1 = time.time()
        mp_dict = {}
        with Pool(processes=processes) as pool:
            for seq_name in list(data_peaks.keys()):
                data = data_peaks[seq_name]
                n_peaks = len(data)
                print("seq name: %s | n peaks: %s" % (seq_name, n_peaks))
                seq_results = []
                for p in range(n_peaks):
                    peak_data = data[p]
                    peak_start = peak_data[4]
                    tile_locs = [peak_data[6], peak_data[7], peak_data[8]]
                    bases = peak_data[3]
                    scores = peak_data[2]
                    mat = cython_fnx.seq_to_mat(bases, scores)
                    args = (peak_start, tile_locs, mat, motifs, n_shuf, seed)
                    peak_results = pool.apply_async(get_peak_results, args)
                    seq_results.append(peak_results)
                mp_dict[seq_name] = seq_results
            pool.close()

            # block to get results
            print("")
            print("getting results...")
            results_dict = {}
            for key in mp_dict:
                print("%s (%s)" % (key, time.ctime()))
                peak_res = mp_dict[key]
                new_peak_res = []
                for p in peak_res:
                    res = p.get()
                    new_peak_res.append(res)
                results_dict[key] = new_peak_res
        t2 = time.time()
        print("=================")
        print("time elapsed: %s" % (t2-t1))
    
    else:
        print("serial")
        t1 = time.time()
        results_dict = {}
        for seq_name in list(data_peaks.keys()):
            data = data_peaks[seq_name]
            n_peaks = len(data)
            print("seq name: %s | n peaks: %s" % (seq_name, n_peaks))
            seq_results = []
            for p in range(n_peaks):
                peak_data = data[p]
                peak_start = peak_data[0]
                tile_locs = [peak_data[6], peak_data[7], peak_data[8]]
                bases = peak_data[3]
                scores = peak_data[2]
                mat = cython_fnx.seq_to_mat(bases, scores)
                peak_results = get_peak_results(peak_start, tile_locs, mat, motifs, n_shuf, seed)
                peak_results.extend(tile_locs)
                seq_results.append(peak_results)
            results_dict[seq_name] = seq_results
        t2 = time.time()
        print("=================")
        print("time elapsed: %s" % (t2-t1))
    return results_dict

def parse_pfm(pfm_file, pseudo):
    motif_lens = []
    motifs = {}
    switch = 0
    counter = 0
    motif_len = -1
    motif = ""
    motif_len_map = {}
    with open(pfm_file) as f:
        for line in f:
            tmp = line.strip().split()
            if line.startswith("MOTIF"):
                motif = tmp[1]
                motif_len = -1
            if line.startswith("letter"):
                motif_len = int(tmp[5])
                motif_lens.append(motif_len)
                motif_mat = np.zeros((4, motif_len))
                n_sites = line.split(" ")[7]
                counter = 0
            if line.startswith("  0.") or line.startswith("  1."):
                tmp_ar = np.asarray(tmp)
                motif_mat[:,counter] = tmp_ar
                counter += 1
            if counter == motif_len and counter != 0:
                if pseudo:
                    counts = (motif_mat * int(n_sites)) + 1
                    motif_mat = counts / np.sum(counts, axis=0)   
                reverse_mat = reverse_complement_pwm(motif_mat) 
                motifs[motif] = [motif_mat, reverse_mat]
                motif_len_map[motif] = motif_len
    motif_lens = set(motif_lens)
    return motifs, motif_lens, motif_len_map

def reverse_complement_pwm(matrix):
    reverse_matrix = np.flipud(np.fliplr(matrix))
    return reverse_matrix

def scale_range(data, minTo, maxTo):
    """
    function to scale data linearly to a new min/max value set
    
    parameters
    ----------
    data: array like, list of numbers
    minTo: float, minimum of new range desired
    maxTo: float, maximum of new range desired
    
    returns
    -------
    scaled_data: array like, new list of numbers (appropriately scaled)
    """
    minFrom = np.min(data)
    maxFrom = np.max(data)
    
    scaled_data = []
    
    for point in data:
        new_point = minTo + (maxTo - minTo) * ((point - minFrom)/(maxFrom - minFrom))
        scaled_data.append(new_point)
    
    return scaled_data

def loss_score(row, col):
    raw_score = row[col]
    if raw_score < 0:
        return abs(raw_score)
    else:
        return 0

def gain_score(row, col):
    raw_score = row[col]
    if raw_score > 0:
        return abs(raw_score)
    else:
        return 0

def pwm_to_bits(pwm):
    total_i = 2 + np.sum(np.nan_to_num(pwm * np.log2(pwm)), axis=0)
    bits = pwm * total_i
    return bits

def moving_average(bandwidth, scores):
    scores_filt = [np.nan] * int(np.floor(bandwidth/2))
    for i in range(0, len(scores)-bandwidth+1):
        scores_ = scores[i:i+bandwidth]
        avg = np.nanmean(scores_)
        scores_filt.append(avg)
    scores_filt.extend([np.nan] * int(np.floor(bandwidth/2)))
    return scores_filt

def in_peak(row, offset, widths):
    delpos = row["delpos"]
    rawpos = int(delpos)-offset
    for start, stop in widths:
        if rawpos >= start and rawpos < stop:
            return "peak"
    return "no peak"

def find_peaks(peak_cutoff, seq_len, df, scores_filt, scaled_scores, bases, offset, max_motif_len):
    # find tile start & end if given in the dataframe
    if "tile_chr" in df.columns and "tile_start" in df.columns and "tile_end" in df.columns:
        tile_chr = df["tile_chr"].iloc[0]
        tile_start = df["tile_start"].iloc[0]
        tile_end = df["tile_end"].iloc[0]
    else:
        tile_chr = "NA"
        tile_start = 0
        tile_end = seq_len
    
    # find the average standard deviation across nucleotides and calculate peak threshold
    avg_sd = df["sd"].mean()
    if avg_sd < 0:
        scaled_peak_thresh = peak_cutoff * (avg_sd)
    else:
        scaled_peak_thresh = peak_cutoff * (avg_sd * 2.5)
    
    # find the gradient of the scores
    grad = np.gradient(scores_filt)
    asign = np.sign(grad)
    
    # find all nucleotides t hat are above peak threshold
    peak_locs = [i for i, x in enumerate(scores_filt) if x > scaled_peak_thresh]
    
    # iterate through these to find peak starts and ends
    starts = []
    stops = []
    for idx in peak_locs:
        for x, y in zip(starts, stops):
            if idx >= x and idx <= y:
                continue
        peak_val = scores_filt[idx]
        val = peak_val
        l_idx = idx-1
        r_idx = idx+1
        try:
            l_val = scores_filt[l_idx]
        except IndexError:
            l_idx = 0
            l_val = scores_filt[idx]
        try:
            r_val = scores_filt[r_idx]
        except IndexError:
            r_idx = len(scores_filt)-1
            r_val = scores_filt[-1]
        r_val = scores_filt[r_idx]
        
        c_sign = asign[idx]
        l_sign = asign[l_idx]
        r_sign = asign[r_idx]
        
        the_start = False
        while l_sign == c_sign:
            c_sign = l_sign
            prev_l_sign = l_sign
            try:
                l_idx -= 1
                l_sign = asign[l_idx]
            except:
                the_start = True
                l_idx = 0
                l_sign = np.nan
        
        if l_idx not in starts:
            starts.append(l_idx)
        elif the_start:
            starts.append(l_idx)
        
        c_sign = asign[idx]
        the_end = False
        while c_sign == r_sign:
            c_sign = r_sign
            if r_idx < seq_len:
                r_idx += 1
                r_sign = asign[r_idx]
            else:
                the_end = True
                r_idx = seq_len-1
                r_sign = np.nan
        
        if r_idx not in stops:
            stops.append(r_idx)
        elif the_end:
            stops.append(the_end)
    widths = list(zip(starts, stops))
    
    # fix overlapping peaks
    # decide to merge or break them up by checking score at overlap point
    if len(widths) < 2:
        non_ov_widths = widths
    else:
        non_ov_widths = []
        for i in range(len(widths)):
            start = widths[i][0]
            end = widths[i][1]
            if i == 0:
                prev_start = -3
                prev_end = -2
            else:
                prev_start = widths[i-1][0]
                prev_end = widths[i-1][1]
            if prev_end >= start:
                score_there = scores_filt[prev_end]
                prev_start = non_ov_widths[-1][0]
                prev_end = non_ov_widths[-1][1]
                if score_there > (0.75 * scaled_peak_thresh):
                    # merge the peaks
                    non_ov_widths = non_ov_widths[:-1]
                    non_ov_widths.append((prev_start, end))
                else:
                    # just add start and end
                    non_ov_widths = non_ov_widths[:-1]
                    non_ov_widths.append((prev_start, prev_end - 1))
                    non_ov_widths.append((start + 1, end))
            else:
                non_ov_widths.append((start, end))
    
    # make sure every peak widths contain a peak idx
    # limit to peaks >= 4 
    fixed_widths = []
    for w in non_ov_widths:
        if w[1] - w[0] >= 4:
            for idx in peak_locs:
                if idx >= w[0] and idx <= w[1]:
                    if w not in fixed_widths:
                        fixed_widths.append((w[0], w[1]))
    
    # extract scores & bases for best peaks
    peak_info = []
    for w in fixed_widths:
        
        # expand peak to be max motif length if it isn't
        w_diff = max_motif_len - (w[1]-w[0])
        if w_diff > 0:
            w_add = math.ceil(w_diff/2.)
        else:
            w_add = 0
        start = w[0]-int(w_add)
        end = w[1]+int(w_add)

        if start < 0:
            start = 0
            end += int(w_add)
        if end > seq_len:
            start -= int(w_add)+(end-seq_len)
            end = int(seq_len)

        peak_scores = scaled_scores[start:end]
        peak_bases = bases[start:end]
        assert len(peak_scores) >= max_motif_len
        
        peak_info.append((w[0], w[1], peak_scores, peak_bases, start, end, tile_chr, tile_start, tile_end))

    # add peak locations to df
    df["peak"] = df.apply(in_peak, offset=offset, widths=fixed_widths, axis=1)
    
    return fixed_widths, peak_info, df
