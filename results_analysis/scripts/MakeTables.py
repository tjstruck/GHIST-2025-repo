from pathlib import Path
import os, glob, yaml, sys
import numpy as np, matplotlib.pyplot as plt

sys.path.insert(1, 'scripts')
from scoring_sweeps import ScoreSweeps
from scoring_yaml import relative_root_mean_squared_error
from table_dicts import *

table_dir = "results/latex_tables"
os.makedirs(table_dir, exist_ok=True)

res_dd = {}

for groundtruthfi in glob.glob("data/groundtruth/*"):
    chal = Path(groundtruthfi).stem.split('GHIST_2025_')[-1].split('_final_goldstandard')[0]
    res_dd[chal] = {}
    chal_res = open(f"{table_dir}/{chal}_res_table.latex", 'w')
    if 'sweep' in chal:
        keys = ['Recall', 'Size']
        chal_res.write(f"F1 & Recall & Size (x100 kb) & Competitor & Approach\\\\\n")
        chal_res.write(f"truth & {len(open(groundtruthfi).readlines())} & 1e-05 & & \\\\\n")
    else:
        keys, col_dict = demo_dict(chal)
        chal_res.write(f"RRMSE & {' & '.join([col_dict[key] for key in keys])} & Competitor & Approach\\\\\n")

        with open(groundtruthfi) as stream:
            try:
                groundtruth = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        chal_res.write(f"truth & {' & '.join([str(groundtruth['parameters'][key]) for key in keys])} & & \\\\\n")
    # print(chal, len(glob.glob(f"data/final_submissions/{chal}/*"))/2)
    for fname in glob.glob(f"data/final_submissions/{chal}/*"):
        stem = Path(fname).stem
        uid = stem.split("_")[-1]
        if uid.isdigit():
            continue
        res_dd[chal][uid] = {}
        
        if 'sweep' in chal:
            sweepscore = ScoreSweeps(fname, groundtruthfi)
            sweepscore.calculate_stats()
            res_dd[chal][uid]['Recall'] = sweepscore.true_positive_count / (sweepscore.true_positive_count + sweepscore.false_negative_count)
            res_dd[chal][uid]['Size'] = np.sum(sweepscore.interval_lengths) / 1e5
            res_dd[chal][uid]['F1'] = sweepscore.f1
            sort_by = "F1"

        else:
            with open(fname) as stream:
                try:
                    submission = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print(exc)

            # keys = list(groundtruth['parameters'].keys())
            # keys.sort()
            RRMSE = float(relative_root_mean_squared_error(np.array([groundtruth['parameters'][key] for key in keys]), np.array([submission['parameters'][key] for key in keys])))
            res_dd[chal][uid]['RRMSE'] = RRMSE
            for key in keys:
                res_dd[chal][uid][key] = float(submission['parameters'][key])
            sort_by = "RRMSE"

        # print(chal, uid, RRMSE, recall, size, f1)
    users_by_truthness = sorted(res_dd[chal], key=lambda u: res_dd[chal][u][sort_by], reverse=True)
    for uid in users_by_truthness:
        chal_res.write(f"{round(res_dd[chal][uid][sort_by], 4)} & {' & '.join([str(round(res_dd[chal][uid][key], 4)) for key in keys])} & {uid} & insert \\\\\n")
    chal_res.close()


# finals = glob.glob("data/final_submissions/*/*")

# participant = '.'.join(finals[0].split('/')[-1].split(f"{finals[0].split('/')[-2]}_")[-1].split('.')[:-1])

# # os.system('vcftools --gzvcf data/vcf/GHIST_2025_multisweep.15.final.vcf.gz --TajimaD 10000 --out results/out')
# # os.system('vcftools --gzvcf data/vcf/GHIST_2025_multisweep.15.final.vcf.gz --window-pi 10000 --out results/out')

# truth = [int(line.split()[1]) for line in open("data/groundtruth/GHIST_2025_multi_sweep_final_goldstandard.bed")]

# data = np.genfromtxt(fname="results/out.windowed.pi", delimiter="\t", skip_header=1)

# for fname in glob.glob("data/final_submissions/multi_sweep/*"):
#     stem = Path(fname).stem
#     uid = stem.split("_")[-1]
#     if uid.isdigit():
#         continue

#     fid = open(fname)
#     intervals = [(int(line.split()[1]), int(line.split()[2])) for line in fid if line.strip()]

#     fig = plt.figure(231412, figsize=(6, 2))
#     fig.clear()
#     ax = fig.add_subplot(111)

#     ax.plot(data[:,1]/1e6, data[:,3], '-', ms=2, color='C02', lw=0.5)
#     ax.set_ylabel("Nucleotide diversity", color='C02', fontsize=12)

#     for loc in truth:
#         ax.axvline(x=loc/1e6, color='k', lw=0.8, ls='--', alpha=0.7)

#     for left, right in intervals:
#         ax.axvspan(left/1e6, right/1e6, color='red', alpha=0.3)

#     ax.set_xlabel("Genomic position (Mb)")
#     #fig.text(0.5, 1.0, "Multiple Sweep Detection Challenge", va='top', ha='center', fontsize=15, color='blue')
#     fig.tight_layout(pad=0.1)
#     fig.subplots_adjust(top=0.89)
#     fig.savefig(f"results/submission_{uid}_multisweep.pdf")