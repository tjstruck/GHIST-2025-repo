#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=sim-recomb
#SBATCH --output=hpc_output/%x-%A-%a.out
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjstruck@arizona.edu
###SBATCH --partition=windfall
#SBATCH --partition=high_priority
###SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
###SBATCH --mem-per-cpu=5gb
###SBATCH --constraint=hi_mem
###SBATCH --mem-per-cpu=32gb
#SBATCH --time=120:00:00
#SBATCH --array=1-272
import sys,os
# Python-specific stuff
print('Script running\n')
if 'SLURM_SUBMIT_DIR' in os.environ:
#       # Set module search path so it will search in qsub directory
      sys.path.insert(0, os.environ['SLURM_SUBMIT_DIR'])
#       # Set current working directory to qsub directory
      os.chdir(os.environ['SLURM_SUBMIT_DIR'])
# Which process am I?
process_ii = int(os.environ.get('SLURM_ARRAY_TASK_ID',1))-1
print(process_ii)

os.system("module load bcftools htslib samtools")

import stdpopsim, demesdraw, demes, tskit, time, cyvcf2
import matplotlib.pyplot as plt

shape = 0.21
scale = 3.28e-3

chroms = ["2L", "2R", "3L", "3R"]
dems = ["SplitMig_testing_20kPop", "SplitMig_final_20kPop", "AfricanGrowth_initial_10kPop", "AfricanGrowth_final_10kPop", ]

run = []
for dem in dems:
    for chrom in chroms:
        if "SplitMig" in dems:
            recomb = 0.5
        else:
            recomb = 1
        run.append([chrom, dem, recomb])

chrom, dem, recomb = run[process_ii]

species = stdpopsim.get_species("DroMel")
model = species.get_demographic_model(dem)
contig_temp = species.get_contig(
    chrom, mutation_rate=model.mutation_rate, genetic_map="ComeronCrossoverV2_dm6"
)
print(f"{contig_temp.recombination_map.mean_rate:.4e}")
import msprime
positions = contig_temp.recombination_map.position
rates = contig_temp.recombination_map.rate*float(recomb)
new_recomb = msprime.RateMap(position=positions, rate=rates)
contig = stdpopsim.Contig(
    mutation_rate=model.mutation_rate, 
    recombination_map=new_recomb,
    genetic_map=contig_temp.genetic_map
)
contig.coordinates = (chrom, contig.coordinates[1], contig.coordinates[2])

print(f"{contig.recombination_map.mean_rate:.4e}")
recomb_mean = f"{contig.recombination_map.mean_rate:.4e}"

_, chal, extra = dem.split("_")

run_dir = f"GHIST_2025-{dem}-shape_{shape}-scale_{scale}-recombScale_{recomb}"
os.makedirs(run_dir, exist_ok=True)

if "AfricanGrowth" in dem:
    samples = {"AFR": 32}
else:
    samples = {"AFR": 12, "EUR": 23}

graph = model.model.to_demes()
demes.dump(graph, f"{run_dir}/demes_{dem}.{chal}.yaml")
demesdraw.tubes(model.model.to_demes(), title=f"{dem}, {chal} Challenge")
plt.savefig(f"{run_dir}/demes_{dem}.{chal}.png", dpi=300)

exons = species.get_annotations("FlyBase_BDGP6.32.51_CDS")
exon_intervals = exons.get_chromosome_annotations(chrom)
species.get_genetic_map("ComeronCrossoverV2_dm6")

engine = stdpopsim.get_engine("msprime")
print("running msprime")
ts = engine.simulate(
    model,
    contig,
    samples,
    seed=248,
)

ts.dump(f"{run_dir}/{dem}.{chrom}.{chal}.trees")

# import numpy as np
# position_transform = lambda x: 1 + np.array(x)
vcf_fi = open(f"{run_dir}/{dem}.{chrom}.{chal}.raw.vcf","w")
# vcf_fi.write(mts.as_vcf(position_transform=position_transform))
vcf_fi.write(ts.as_vcf())
vcf_fi.close()


# ts = tskit.load(f"{run_dir}/{dem}.{chrom}.{chal}.trees")

if "AfricanGrowth" in dem:
    fi = open(f"{run_dir}/{_}.popfile.txt","w")
    for i in range(samples['AFR']):
        fi.write(f"modern_{i+1}\tmodern\n")
    fi.close()

    fi = open(f"{run_dir}/{_}.new_samples.txt","w")
    for i in range(samples['AFR']):
        fi.write(f"tsk_{i}\tmodern_{i+1}\n")
    fi.close()
else:
    fi = open(f"{run_dir}/{_}.new_samples.txt","w")
    for i in range(samples['AFR']):
        fi.write(f"tsk_{i}\tAFR_{i+1}\n")
    for i in range(samples['AFR'],samples['AFR']+samples['EUR']):
        fi.write(f"tsk_{i}\tEUR_{i-samples['AFR']+1}\n")
    fi.close()

    fi = open(f"{run_dir}/{_}.popfile.txt","w")
    for i in range(samples['AFR']):
        fi.write(f"AFR_{i+1}\tAFR\n")
    for i in range(samples['AFR'],samples['AFR']+samples['EUR']):
        fi.write(f"EUR_{i-samples['AFR']+1}\tEUR\n")
    fi.close()

fi = open(f"{run_dir}/{_}.{chrom}.chr_map.txt","w")
fi.write(f"1\t{chrom}")
fi.close()
