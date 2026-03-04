import allel, glob, pickle, os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# In a BED, the start is the equivalent of 0-based position for SNPs
sweeps = [int(ele.split('\t')[1])+1 for ele in open("data/groundtruth/GHIST_2025_multi_sweep_final_goldstandard.bed").readlines()]

os.makedirs("results/multi_sweep/", exist_ok=True)

for i in list(range(len(sweeps)-1,0,-1)):
    print("sweeps distance:",sweeps[i] - sweeps[i-1])

# Path to your VCF (bgzipped + tabix indexed is ideal)
vcf_fi = "data/vcf/GHIST_2025_multisweep.15.final.vcf.gz"

callset = allel.read_vcf(
    vcf_fi,
    fields=["samples", "variants/CHROM", "variants/POS", "calldata/GT"],
    alt_number=1  # simplify to biallelic; drop multiallelic for many stats
)

samples = np.array(callset["samples"])
chrom = np.array(callset["variants/CHROM"])
pos   = np.array(callset["variants/POS"])
gt    = allel.GenotypeArray(callset["calldata/GT"])

# # If necessary, subset to a chromosome
# mask = chrom == chrom[0]  # single contig example
# pos = pos[mask]
# gt  = gt.compress(mask, axis=0)

print(gt.shape, len(pos), len(samples))

# Diploid example; if haploid adjust accordingly
ac = gt.count_alleles()                      # overall allele counts
gn = gt.to_n_alt()                           # genotype as {0,1,2}
haps = gt.to_haplotypes()                    # haplotype array (n_variants, n_haplotypes)

# Define windows (e.g., 10 kb)
win_size = 50_000
chrom_len = pos.max()
win_starts = np.arange(1, chrom_len, win_size)
win_ends = win_starts + win_size

loc = allel.SortedIndex(pos)

# Nucleotide diversity (π)
pi, pi_windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=win_size)

# Tajima's D
tajima_d, td_windows, counts_td = allel.windowed_tajima_d(pos, ac, size=win_size)

# Plot
fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

axs[0].plot(pi_windows.mean(axis=1), pi, lw=1)
axs[0].set_ylabel("π (diversity)")
axs[0].set_title(f"Windowed Sweeps, size = {win_size}")

axs[1].plot(td_windows.mean(axis=1), tajima_d, lw=1, color="darkred")
axs[1].axhline(0, color="gray", lw=0.8, ls="--")
axs[1].set_ylabel("Tajima's D")
axs[1].set_xlabel(f"Position on {chrom[0]} (bp)")

# Plot true sweeps
for ax in axs:
    for sweep in sweeps:
        ax.axvline(sweep, color="black", linestyle="--")

# plt.suptitle("Windowed Sweeps")
plt.tight_layout()
plt.savefig(f"sweep_visualization_base_{win_size}_window_size.png", dpi=150)

pickle.dump(fig, open(f"sweep_visualization_base_{win_size}_window_size.bpkl",'wb'))

plt.show()

# import pickle

# # Save to memory
# fig_bytes = pickle.dumps(fig)

# # Restore elsewhere
# fig_copy = pickle.loads(fig_bytes)
for fname in glob.glob("data/final_submissions/multi_sweep/*"):
    stem = Path(fname).stem
    uid = stem.split("_")[-1]
    if uid.isdigit():
        continue
    fig_copy = pickle.load(open(f"sweep_visualization_base_{win_size}_window_size.bpkl",'rb'))
    axs_copy = fig_copy.axes

    fid = open(fname)
    intervals = [(int(line.split()[1]), int(line.split()[2])) for line in fid if line.strip()]

    for ax in axs_copy:
        for left, right in intervals:
            ax.axvspan(left, right, color='green', alpha=0.3)
    plt.savefig(f"results/multi_sweep/sweep_visualization_{uid}_{win_size}_window_size.png", dpi=150)
    plt.show()



