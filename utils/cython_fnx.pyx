import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.double
ctypedef np.double_t DTYPE_t

def seq_to_mat(list bases, list scores, double min_score=1e-05):
    cdef np.ndarray[DTYPE_t, ndim=2] mat = np.zeros((4, len(bases)), dtype=DTYPE)
    cdef str base
    cdef double score
    
    for i in range(len(bases)):
        base = bases[i]
        score = scores[i]
        if score != 0:
            if base == "A":
                mat[0,i] = score
            elif base == "C":
                mat[1,i] = score
            elif base == "G":
                mat[2,i] = score
            elif base == "T":
                mat[3,i] = score
        else:
            if base == "A":
                mat[0,i] = min_score
            elif base == "C":
                mat[1,i] = min_score
            elif base == "G":
                mat[2,i] = min_score
            elif base == "T":
                mat[3,i] = min_score
    return mat


cpdef double calculate_motif_score(np.ndarray[DTYPE_t, ndim=2] mat, np.ndarray[DTYPE_t, ndim=2] pwm, double info_content):
    cdef int m = mat.shape[0]
    cdef int n = mat.shape[1]
    cdef int i = 0
    cdef int j = 0
    cdef double zero_scores = 0
    cdef double motif_score = 0
    cdef double n_base_matches = 0
    cdef double mat_max
    cdef int mat_idx
    cdef double pwm_max
    cdef int pwm_idx
    cdef double norm_motif_score
    
    # non-np matrix mult
    for i in range(m):
        for j in range(n):
            if pwm[i,j] < 0.03:
                zero_scores += mat[i,j]**2
            else:
                motif_score += mat[i,j]*pwm[i,j]
    motif_score = motif_score - zero_scores
    
    # find num of matching bases
    for j in range(n):
        mat_max = -1
        mat_idx = -1
        pwm_max = -1
        pwm_idx = -1
        for i in range(m):
            if mat[i, j] > mat_max:
                mat_max = mat[i, j]
                mat_idx = i
            if pwm[i, j] > pwm_max:
                pwm_max = pwm[i, j]
                pwm_idx = i
        if mat_idx == pwm_idx:
            n_base_matches += 1
    
    norm_motif_score = (motif_score * n_base_matches) / info_content
    return norm_motif_score


def slide_motifs_across_seq(np.ndarray[DTYPE_t, ndim=2] mat, int motif_len, np.ndarray[DTYPE_t, ndim=2] pwm, double info_content):
    cdef int n_scores = mat.shape[1]-motif_len+1
    cdef int i = 0 
    cdef np.ndarray[DTYPE_t, ndim=2] scores = np.zeros((1, n_scores), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] mat_ = np.zeros((4, motif_len), dtype=DTYPE)

    for i in range(0, n_scores):
        mat_ = mat[:, i:i+motif_len]
        norm_motif_score = calculate_motif_score(mat_, pwm, info_content)
        scores[0,i] = norm_motif_score
    return scores

def get_motif_results(np.ndarray pwm, double info_content, np.ndarray seq, int peak_start, int n_shuf, int seed):
    cdef int motif_len = pwm.shape[1]
    cdef int n_scores = seq.shape[1]-motif_len+1
    cdef np.ndarray[DTYPE_t, ndim=2] motif_scores = np.zeros((1, n_scores), dtype=DTYPE)
    cdef int i
    cdef double val
    cdef double max_score = -100
    cdef int max_idx
    cdef np.ndarray[DTYPE_t, ndim=2] shuf_motif_scores = np.zeros((1, n_scores), dtype=DTYPE)
    cdef int n = 0
    cdef double shuf_max_score
    cdef int n_shuf_motifs_with_higher_score = 0
    cdef double pval = 0.99999
    cdef int motif_start
    cdef int motif_end
    
    # set seed
    np.random.seed(seed)
    
    # first calculate score using original mat
    motif_scores = slide_motifs_across_seq(seq, motif_len, pwm, info_content)

    # non-np max and argmax
    for i in range(n_scores):
        val = motif_scores[0,i]
        if val > max_score:
            max_score = val
            max_idx = i
    
    # then find shuffled scores for that motif in that peak
    c_seq = seq.copy()
    for n in range(n_shuf):
        shuf_max_score = -100
        np.random.shuffle(c_seq)
        np.random.shuffle(c_seq.T)
        shuf_motif_scores = slide_motifs_across_seq(c_seq, motif_len, pwm, info_content)
        
        for i in range(n_scores):
            val = shuf_motif_scores[0,i]
            if val > shuf_max_score:
                shuf_max_score = val
        if shuf_max_score >= max_score:
            n_shuf_motifs_with_higher_score += 1                     

    # calculate p value based on these shuffled scores
    pval = float(n_shuf_motifs_with_higher_score+1)/(n_shuf+1)

    # find the start and end of the best hit motif in overall coords
    motif_start = peak_start + max_idx + 1
    motif_end = motif_start + motif_len
    results = [motif_start, motif_end, max_score, pval]  
    return results
