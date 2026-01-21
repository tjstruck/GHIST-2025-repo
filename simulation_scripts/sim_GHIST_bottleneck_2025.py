#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=sim-bottleneck
#SBATCH --output=hpc_output/%x-%A.out
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjstruck@arizona.edu
###SBATCH --partition=windfall
#SBATCH --partition=high_priority
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
###SBATCH --mem-per-cpu=5gb
###SBATCH --constraint=hi_mem
###SBATCH --mem-per-cpu=32gb
#SBATCH --time=24:00:00
#SBATCH --array=1-2

import dadi
import msprime
import demes
import demesdraw
import matplotlib.pyplot as plt
import random
import ast
import sys,os
# Python-specific stuff
print('Script running\n')
if 'SLURM_SUBMIT_DIR' in os.environ:
#       # Set module search path so it will search in qsub directory
      sys.path.insert(0, os.environ['SLURM_SUBMIT_DIR'])
#       # Set current working directory to qsub directory
      # os.chdir(os.environ['SLURM_SUBMIT_DIR'])
# Which process am I?
process_ii = int(os.environ.get('SLURM_ARRAY_TASK_ID',1))-1
print(process_ii)

os.system("module load bcftools htslib samtools")

challenge_name = "GHIST_2025_bottleneck_NoDiscrete"
os.makedirs(challenge_name, exist_ok=True)

seq_l = "1e8"

def msprime_bottleneck(Na, p):
    (nu, T) = p
    dem = msprime.Demography()
    dem.add_population(name="modern", initial_size=Na*nu)
    dem.add_population(name="ancestral", initial_size=Na)
    dem.add_population_split(time=2*Na*T, derived=["modern"], ancestral="ancestral")
    return dem

Na_I = 21378
nu_I = 0.523
print(f"Pop-size, testing challenge: {nu_I*Na_I}")
# 10763.194

T_I = 0.1231
print(f"Generations, testing challenge: {2*Na_I*T_I}")
# 321.98036
print(f"T/nu: {T_I/nu_I}")
# 0.014957472660996356
p_I = (nu_I, T_I)


dem_I = msprime_bottleneck(Na_I, p_I)
dem_I.sort_events()

graph = dem_I.to_demes()
demes.dump(graph, f"{challenge_name}/demes_{challenge_name}.testing.yaml")
demesdraw.tubes(graph, title=f"{challenge_name.replace('_', ' ')}, Testing Challenge")
plt.savefig(f"{challenge_name}/demes_{challenge_name}.testing.png", dpi=300)
# os.system(f"demesdraw tubes --title '{challenge_name}/demes_{challenge_name}.testing.png' {challenge_name}/demes_{challenge_name}.testing.yaml")

Na_F = 27293
nu_F = 0.2975
print(f"\n\nPop-size, one-shot challenge: {nu_F*Na_F}")
# 8119.6675
T_F = 0.037109
print(f"Generations, one-shot challenge: {2*Na_F*T_F}")
# 3117.351874
print(f"T/nu: {T_F/nu_F}")
# 0.19196302521008404
p_F = (nu_F, T_F)

ns = {"modern":22} # individuals
ploidy = 2 # diploid
mut = 2.4245e-8 # mutation rate
recomb = 1.5271e-08

dem_F = msprime_bottleneck(Na_F, p_F)
dem_F.sort_events()

# os.makedirs(f"{challenge_name}", exist_ok=True)

graph = dem_F.to_demes()
demes.dump(graph, f"{challenge_name}/demes_{challenge_name}.final.yaml")
demesdraw.tubes(graph, title=f"{challenge_name.replace('_', ' ')}, One-Shot Challenge")
plt.savefig(f"{challenge_name}/demes_{challenge_name}.final.png", dpi=300)
# os.system(f"demesdraw tubes --title '{challenge_name}/demes_{challenge_name}.final.png' {challenge_name}/demes_{challenge_name}.final.yaml")

print("Demes file made")


fi = open(f"{challenge_name}/{challenge_name}.popfile.txt","w")
for i in range(1,ns['modern']+1):
    fi.write(f"wisent_{i+1}\twisent\n")
fi.close()

fi = open(f"{challenge_name}/{challenge_name}.popfile.txt","w")
for i in range(ns['modern']):
    fi.write(f"sample_{i+1}\tmodern\n")
fi.close()

fi = open(f"{challenge_name}/{challenge_name}.new_samples.txt","w")
for i in range(ns['modern']):
    fi.write(f"tsk_{i}\tsample_{i+1}\n")
fi.close()

ts = msprime.sim_ancestry(
    samples=ns, 
    demography=[dem_I, dem_F][process_ii], 
    sequence_length=ast.literal_eval(seq_l), 
    recombination_rate=recomb, 
    ploidy=ploidy,
    random_seed=12
    )

dataset = ['testing', 'final']

mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=False, random_seed=5566)
mts.num_mutations

vcf_fi = open(f"{challenge_name}/{challenge_name}.{dataset[process_ii]}.raw.vcf","w")
vcf_fi.write(mts.as_vcf())
vcf_fi.close()

mts.dump(f"{challenge_name}/{challenge_name}.{dataset[process_ii]}.trees")
