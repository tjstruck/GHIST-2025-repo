'''
Run after MakeTabels.py
'''
import demes, demesdraw, os
import matplotlib.pyplot as plt

os.makedirs("results/paper/", exist_ok=True)

fig, axs = plt.subplots(2, 1, figsize=(4, 6))#, sharex=True)

true_demes = demes.load(f"data/demes_for_paper/demes_GHIST_2025_secondary_contact.final.yaml")

w = demesdraw.utils.separation_heuristic(true_demes)
positions_true = dict(mainland=1.3*w, island=2*w)
demesdraw.tubes(true_demes, seed=1234, ax=axs[0])#, positions=positions_true)

best_demes = demes.load(f"data/demes_for_paper/secondary_contact_srong.yaml")
w = demesdraw.utils.separation_heuristic(best_demes)
positions_best = dict(Population_1=1.5*w, Population_2=2.5*w)
positions_best['']=2*w
best_demes.demes[0].name = ''
best_demes.demes[1].ancestors = ['']
best_demes.demes[2].ancestors = ['']
demesdraw.tubes(best_demes, seed=1234, ax=axs[1])#, positions=positions_best)

axs[0].set_ylabel = 'Generations ago'
axs[1].set_ylabel = 'Generations ago'

# axs[0].tick_params(labelbottom=True)
# axs[0].set_xticks([positions_true["mainland"], positions_true["island"]], ["Mainland", "Island"])
# axs[0].set_xticks([positions_true["mainland"], positions_true["island"]], ["Population_1", "Population_2"])

for ax, label in zip(axs, ["A", "B"]):
    ax.text(
        -0.27, 1.1, label,          # adjust -0.12 to taste
        transform=ax.transAxes,
        fontsize=12, #fontweight="bold",
        va="top", ha="left"
    )

plt.tight_layout()
plt.savefig(f"results/2025_Chal2.png")
plt.savefig(f"results/paper/2025_Chal2.pdf")
plt.clf()