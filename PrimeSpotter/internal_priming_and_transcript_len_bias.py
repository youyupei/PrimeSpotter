import pysam
import numpy as np
import pandas as pd


# find plyA site
def find_polyA_site(txname, fasta_handle, window_size=10, 
                    minAorT=8, merge_dist=3, flank_butter=25):
    """
    return:
        number of polyA site
        length of the transcript
    """
    # get refernce names
    seq = fasta_handle.fetch(txname)
    seq = np.array([int(x=='A') for x in seq])
    # get sliding window (size = 10) sum 
    window_sum = np.convolve(seq, np.ones(window_size), 'valid')
    # find location of consecutive window sum >= minAT
    polyA_site = np.where(window_sum >= minAorT)[0]
    if len(polyA_site):
        # merge consecutive polyA when the difference is <=3
        polyA_site = np.split(polyA_site, np.where(np.diff(polyA_site) > merge_dist)[0]+1)
        return len(polyA_site), len(seq)
    else:
    
        return 0, len(seq)




if __name__ == "__main__":
    # transcript
    fasta_file_path = "/home/users/allstaff/you.yu/project/LR_multiome/results/flames_out/scm_old/transcript_assembly.fa"

    fasta = pysam.FastaFile(fasta_file_path)

    rst = {
        "transcript_id": [],
        "transcript_len": [],
        "n_polyA_site": []
    }
    for txname in fasta.references:
        n, l = find_polyA_site(txname, fasta)
        rst['transcript_id'].append(txname)
        rst["n_polyA_site"].append(n)
        rst["transcript_len"].append(l)

    df = pd.DataFrame(rst)
    df.to_csv("transcript_len_and_polyA_site.csv", index=False)



