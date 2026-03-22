import allel, os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

os.makedirs("results/paper/", exist_ok=True)

chal_type2chrom = {
    'testing':'21',
    'final':'15'
}

for chal_deets in [('multi_sweep_bgs', 'chal9', 'final'), ('multi_sweep', 'chal8', 'final'), ('multi_sweep', 'chal8', 'testing')]:
    chal, chal_num, chal_type = chal_deets
    chrom = chal_type2chrom[chal_type]
    # In a BED, the start is the equivalent of 0-based position for SNPs
    sweeps = [int(ele.split('\t')[1])+1 for ele in open(f"data/groundtruth/GHIST_2025_{chal}_{chal_type}_goldstandard.bed").readlines()]

    os.makedirs(f"results/paper/", exist_ok=True)

    for i in list(range(len(sweeps)-1,0,-1)):
        print(f"sweeps distance:",sweeps[i] - sweeps[i-1])

    # Path to VCF
    if "_bgs" in chal:
        vcf_fi = f"data/vcf/GHIST_2025_multisweep.growth_bg.{chrom}.{chal_type}.vcf.gz"
    else:
        vcf_fi = f"data/vcf/GHIST_2025_multisweep.{chrom}.{chal_type}.vcf.gz"

    callset = allel.read_vcf(
        vcf_fi,
        fields=["samples", "variants/CHROM", "variants/POS", "calldata/GT"],
        alt_number=1
    )

    samples = np.array(callset["samples"])
    chrom = np.array(callset["variants/CHROM"])
    pos   = np.array(callset["variants/POS"])
    gt    = allel.GenotypeArray(callset["calldata/GT"])

    ac = gt.count_alleles()
    gn = gt.to_n_alt()
    haps = gt.to_haplotypes()

    # Define windows (e.g., 10 kb)
    win_size = 10_000
    chrom_len = pos.max()
    win_starts = np.arange(1, chrom_len, win_size)
    win_ends = win_starts + win_size

    loc = allel.SortedIndex(pos)

    # Nucleotide diversity (π)
    pi, pi_windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=win_size)

    # Tajima's D
    tajima_d, td_windows, counts_td = allel.windowed_tajima_d(pos, ac, size=win_size)

    # Plot
    fig, axs = plt.subplots(2, 1, figsize=(6, 3), sharex=True)

    axs[0].plot(pi_windows.mean(axis=1), pi, lw=1)
    axs[0].set_ylabel("π (diversity)")

    axs[1].plot(td_windows.mean(axis=1), tajima_d, lw=1, color="darkred")
    axs[1].axhline(0, color="gray", lw=0.8, ls="--")
    axs[1].set_ylabel("Tajima's D")
    axs[1].set_xlabel(f"Position on simulated human chromosome {chrom[0]} (bp)")

    # Add figure panel labeling
    for ax, label in zip(axs, ["A", "B"]):
        ax.text(
            -0.145, 1.15, label,          # adjust -0.12 to taste
            transform=ax.transAxes,
            fontsize=12, #fontweight="bold",
            va="top", ha="left"
        )

    # Plot true sweeps
    for ax in axs:
        for sweep in sweeps:
            ax.axvline(sweep, color="black", linestyle="--")
    plt.tight_layout()

    # Get submission intervals
    fid = open(f"data/final_submissions/{chal}/{chal}_srong.bed")
    intervals = [(int(line.split()[1]), int(line.split()[2])) for line in fid if line.strip()]
    fid.close()

    if chal_type != 'testing':
        for ax in axs:
            for left, right in intervals:
                ax.axvspan(left, right, color='green', alpha=0.3)
    if chal_type == 'testing':
        plt.savefig(f"results/paper/2025_{chal_num}_{chal_type}.pdf", dpi=150)
    else:
        plt.savefig(f"results/paper/2025_{chal_num}.pdf", dpi=150)
    plt.show()



