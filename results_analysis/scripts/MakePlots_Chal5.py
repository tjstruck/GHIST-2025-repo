'''
Run after MakeTabels.py
'''
import demes, demesdraw, os
import matplotlib.pyplot as plt

os.makedirs("results/paper/", exist_ok=True)

true_demes = demes.load(f"data/demes_for_paper/demes_GHIST_2025_admixture.final.yaml")

w = demesdraw.utils.separation_heuristic(true_demes)
anc_pop = 0.9
anc_pop_off = 0.33
# positions = dict(ancestral=0, ancient_population=0.8 * w, Modern_pop1=0.4 * w, Modern_pop2=1.2 * w, AncPop1=2 * w, AncOut2=3 * w, AncOut3=4 * w, AncOut1=5 * w, AncPop2=6 * w)
positions = dict(ancestral=0, 
                  Modern_pop1=(anc_pop-anc_pop_off) * w, ancient_population=anc_pop * w, Modern_pop2=(anc_pop+anc_pop_off) * w,
                  AncPop1=2.0 * w, AncOut2=2.7 * w, 
                  AncOut3=3.2 * w, AncOut1=3.7 * w, AncPop2=4.2 * w)
# fig, ax = plt.subplots(figsize=(6, 3))
demesdraw.tubes(true_demes, positions=positions, seed=134, labels="xticks")
ax = plt.gca()
fig = plt.gcf()
fig.set_size_inches(6, 2.5)
# ax = fig.axes[0]
# for txt in ax.texts:
#     txt.set_visible(False)
ax.set_ylabel("Generations ago")
ax.set_xticks([positions["Modern_pop1"], positions["Modern_pop2"], positions["AncPop1"], positions["AncOut2"]], ["modern1", "modern2", "ancient1", "ancient2"])
plt.tight_layout()
plt.savefig(f"results/2025_Chal5.png")
plt.savefig(f"results/paper/2025_Chal5.pdf")
plt.show()