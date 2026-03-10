'''
Run after MakeTabels.py
'''
import pickle, demes, demesdraw, os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

res_dd = pickle.load(open("results/results.bpkl",'rb'))

# --- layout ---
fig = plt.figure(figsize=(4.2, 4.2), constrained_layout=True)

# 2 rows, 2 cols; top row is shorter, bottom row taller
gs = GridSpec(
    nrows=2, ncols=2, figure=fig,
    height_ratios=[1, 2],   # top row small, bottom row big
    width_ratios=[1, 1]
)

# Top row: two small subplots
axA = fig.add_subplot(gs[0, 0])  # left small
axB = fig.add_subplot(gs[0, 1])  # right small
axB.set_yscale('log')
axB.set_ylabel('RRMSE')

# Bottom row: one subplot spanning both columns
axC = fig.add_subplot(gs[1, :])  # span all columns
axC.set_xlabel('Generations ago')
axC.set_ylabel('Decline fraction')

true_demes = demes.load(f"data/simulated_demes_for_paper/demes_GHIST_2025_{chal}.final.yaml")
true_demes.demes[0].name = ''
demesdraw.tubes(true_demes, seed=1234, colours='black', ax=axA)
axA.set_ylabel('Generations ago')

for uid in res_dd['bottleneck']:
    axB.axhline(res_dd['bottleneck'][uid]['RRMSE'], color='black', alpha=0.5)
    axC.scatter(res_dd['bottleneck'][uid]['generations'], res_dd['bottleneck'][uid]['post_decline_fraction'], color='black', alpha=0.5)

axC.axhline(8119.6675/27293, color='black')
axC.axvline(2025.6318740000002, color='black')


# --- Panel labels (A, B, C) ---
axA.text(-0.66, 1.15, "A", transform=axA.transAxes)
axB.text(-0.60, 1.15, "B", transform=axB.transAxes)
axC.text(-0.25, 1.05, "C", transform=axC.transAxes)

# plt.tight_layout()
plt.savefig(f"results/2025_Chal1.pdf")
plt.savefig(f"results/2025_Chal1.png")
plt.clf()
