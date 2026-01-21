from pathlib import Path
import os, glob
import numpy as np, matplotlib.pyplot as plt

os.system('vcftools --gzvcf GHIST_2025_multisweep.15.final.vcf.gz --TajimaD 10000')
os.system('vcftools --gzvcf GHIST_2025_multisweep.15.final.vcf.gz --window-pi 10000')

truth = [int(line.split()[1]) for line in open("GHIST_submissions_and_truth/groundtruth/GHIST_2025_multi_sweep_final_goldstandard.bed")]

data = np.genfromtxt(fname="out.windowed.pi", delimiter="\t", skip_header=1)

for fname in glob.glob("GHIST_submissions_and_truth/final_submissions/multi_sweep/*"):
    stem = Path(fname).stem
    uid = stem.split("_")[-1]
    if uid.isdigit():
        continue

    fid = open(fname)
    intervals = [(int(line.split()[1]), int(line.split()[2])) for line in fid if line.strip()]

    fig = plt.figure(231412, figsize=(6, 2))
    fig.clear()
    ax = fig.add_subplot(111)

    ax.plot(data[:,1]/1e6, data[:,3], '-', ms=2, color='C02', lw=0.5)
    ax.set_ylabel("Nucleotide diversity", color='C02', fontsize=12)

    for loc in truth:
        ax.axvline(x=loc/1e6, color='k', lw=0.8, ls='--', alpha=0.7)

    for left, right in intervals:
        ax.axvspan(left/1e6, right/1e6, color='red', alpha=0.3)

    ax.set_xlabel("Genomic position (Mb)")
    #fig.text(0.5, 1.0, "Multiple Sweep Detection Challenge", va='top', ha='center', fontsize=15, color='blue')
    fig.tight_layout(pad=0.1)
    fig.subplots_adjust(top=0.89)
    fig.savefig(f"submission_{uid}_multisweep.pdf")