#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=sim-sec-cont
#SBATCH --output=hpc_output/%x-%A.out
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjstruck@arizona.edu
###SBATCH --partition=windfall
#SBATCH --partition=high_priority
#SBATCH --nodes=1
#SBATCH --ntasks=90
#SBATCH --cpus-per-task=1
###SBATCH --mem-per-cpu=5gb
###SBATCH --constraint=hi_mem
###SBATCH --mem-per-cpu=32gb
#SBATCH --time=72:00:00
#SBATCH --array=1-2

# import dadi
import msprime
import demes
import demesdraw
import matplotlib.pyplot as plt
import numpy as np
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

extra = ''
# extra = '_40kPop'
# extra = '_30kPop'
# extra = '_20kPop'
# extra = '_10kPop'

seq_l = "1e8"

# extra += f"_{seq_l}"
challenge_name = "GHIST_2025_secondary_contact"

os.makedirs(f"{challenge_name}{extra}", exist_ok=True)

# seq_l = ["1e8"][process_ii]



def msprime_secondary_contact_testing(Na, p):
    (nu1, nu2, nu2Cont, TAncSplit, TCont, m) = p
    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=Na)
    dem.add_population(name="mainland", initial_size=Na*nu1)
    dem.add_population(name="island", initial_size=Na*nu2)

    dem.add_population_parameters_change(time=TAncSplit*2, population="ancestral", initial_size=Na*1.5)
    dem.add_population_parameters_change(time=TAncSplit*1.9, population="ancestral", initial_size=Na*0.5)
    dem.add_population_parameters_change(time=TAncSplit*1.75, population="ancestral", initial_size=Na*2)
    # Set ancestral pop to Na some time before TAncSplit
    dem.add_population_parameters_change(time=TAncSplit*1.35, population="ancestral", initial_size=Na)

    dem.add_population_split(time=TAncSplit, derived=["mainland", "island"], ancestral="ancestral")

    # This first bit is how long the mainland population is size Na*nu1
    dem.add_population_parameters_change(time=TAncSplit*0.85, population="mainland", initial_size=Na*nu1)
    # Add in "wiggles"
    dem.add_population_parameters_change(time=TAncSplit*0.8, population="mainland", initial_size=Na*nu1*1.8)
    dem.add_population_parameters_change(time=TAncSplit*0.5, population="mainland", initial_size=Na*nu1*0.25)
    # Set back to known end size
    dem.add_population_parameters_change(time=0, population="mainland", initial_size=Na*nu1)

    # dem.set_symmetric_migration_rate(["mainland", "island"], 0)
    dem.add_migration_rate_change(time=TCont, rate=0, source="mainland", dest="island")
    dem.add_migration_rate_change(time=TCont, rate=0, source="island", dest="mainland")
    dem.add_migration_rate_change(time=0, rate=m, source="mainland", dest="island")
    dem.add_migration_rate_change(time=0, rate=m, source="island", dest="mainland")

    # This first bit is how long the island population is size Na*nu2
    dem.add_population_parameters_change(time=TAncSplit*0.5, population="island", initial_size=Na*nu2)
    # These next bits are how long the island population is size Na*nu2Cont
    # Time based on TCont
    dem.add_population_parameters_change(time=TCont*2, population="island", initial_size=Na*nu2*2.2)
    # Add a wiggle step that ends at TCont
    dem.add_population_parameters_change(time=TCont, population="island", initial_size=Na*nu2*1.2)
    # Add how long the island population is size Na*nu2Cont from TCont to TCont*scaler
    dem.add_population_parameters_change(time=TCont*0.9, population="island", initial_size=Na*nu2Cont)
    dem.add_population_parameters_change(time=TCont*0.45, population="island", initial_size=Na*nu2Cont*0.5)
    # Set back to known end size
    dem.add_population_parameters_change(time=0, population="island", initial_size=Na*nu2Cont)
    return dem

def msprime_secondary_contact_final(Na, p):
    (nu1, nu2, nu2Cont, TAncSplit, TCont, m) = p
    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=Na)
    dem.add_population(name="mainland", initial_size=Na*nu1)
    dem.add_population(name="island", initial_size=Na*nu2)

    dem.add_population_parameters_change(time=TAncSplit*3.12, population="ancestral", initial_size=Na*3.5)
    dem.add_population_parameters_change(time=TAncSplit*2.01, population="ancestral", initial_size=Na*1.25)
    dem.add_population_parameters_change(time=TAncSplit*1.8, population="ancestral", initial_size=Na*0.33)
    # Set ancestral pop to Na some time before TAncSplit
    dem.add_population_parameters_change(time=TAncSplit*1.5, population="ancestral", initial_size=Na)

    dem.add_population_split(time=TAncSplit, derived=["mainland", "island"], ancestral="ancestral")

    # This first bit is how long the mainland population is size Na*nu1
    dem.add_population_parameters_change(time=TAncSplit*0.75, population="mainland", initial_size=Na*nu1)
    # Add in "wiggles"
    dem.add_population_parameters_change(time=TAncSplit*0.55, population="mainland", initial_size=Na*nu1*1.8)
    dem.add_population_parameters_change(time=TAncSplit*0.25, population="mainland", initial_size=Na*nu1*3.5)
    # Set back to known end size
    dem.add_population_parameters_change(time=0, population="mainland", initial_size=Na*nu1)

    # dem.set_symmetric_migration_rate(["mainland", "island"], 0)
    dem.add_migration_rate_change(time=TCont, rate=0, source="mainland", dest="island")
    dem.add_migration_rate_change(time=TCont, rate=0, source="island", dest="mainland")
    dem.add_migration_rate_change(time=0, rate=m, source="mainland", dest="island")
    dem.add_migration_rate_change(time=0, rate=m, source="island", dest="mainland")

    # This first bit is how long the island population is size Na*nu2
    dem.add_population_parameters_change(time=TAncSplit*0.5, population="island", initial_size=Na*nu2)
    # These next bits are how long the island population is size Na*nu2Cont
    # Time based on TCont
    dem.add_population_parameters_change(time=TCont*2, population="island", initial_size=Na*nu2*2.2)
    # Add a wiggle step that ends at TCont
    dem.add_population_parameters_change(time=TCont, population="island", initial_size=Na*nu2*1.2)
    # Add how long the island population is size Na*nu2Cont from TCont to TCont*scaler
    dem.add_population_parameters_change(time=TCont*0.75, population="island", initial_size=Na*nu2Cont)
    dem.add_population_parameters_change(time=TCont*0.45, population="island", initial_size=Na*nu2Cont*0.5)
    # Set back to known end size
    dem.add_population_parameters_change(time=0, population="island", initial_size=Na*nu2Cont)
    return dem

ns = {"mainland":24, "island":11} # individuals
ploidy = 2 # diploid
# based on https://github.com/popsim-consortium/stdpopsim/blob/main/stdpopsim/catalog/MusMus/species.py
# mut = 5.4e-7 # mutation rate
mut = 8.4e-8
# mean recomb from stdpopsim
recomb = 5.386e-08

# Na_I = 100421
# Na_I = 40421 # ancestral pop size
Na_I = 11421
nu1_I = 1.347
nu2_I = 0.373
nu2Cont_I = 4.309

TAncSplit_I = 13345.54
TCont_I = 4297.54

m_I = 15.24e-5

TAncSplitPre1_I = 17000
TAncSplitPre2_I = 16000
TAncSplitPre3_I = 14000
NaPre1_I = 80421
NaPre2_I = 50421
NaPre3_I = 120421

popsize1 = nu1_I*Na_I
popsize2 = nu2_I*Na_I
popsize2Cont = nu2Cont_I*Na_I
print(f"Pop-size 1, testing challenge: {popsize1}")
# 135267.087
print(f"Pop-size 2 split, testing challenge: {popsize2}")
# 37457.033
print(f"Pop-size 2 contact, testing challenge: {popsize2Cont}\n")
# 432714.08900000004

TAncSplit_dadi_I = TAncSplit_I/(2*Na_I)
TCont_dadi_I = TCont_I/(2*Na_I)
print(f"dadi ancient split-time, testing challenge: {TAncSplit_dadi_I}")
# 0.06644795411318351
print(f"dadi secondary contact time, testing challenge: {TCont_dadi_I}\n")
# 0.02139761603648639

print(f"T_nu, pop1 ancient split, testing challenge: {TAncSplit_dadi_I/nu1_I}")
# 0.04933032970540721
print(f"T_nu, pop2 ancient split, testing challenge: {TAncSplit_dadi_I/nu2_I}")
# 0.17814464909700672
print(f"T_nu, pop2 secondary contact, testing challenge: {TCont_dadi_I/nu2Cont_I}\n#####\n\n\n#####")
# 0.004965796248894498

p_I = (nu1_I, nu2_I, nu2Cont_I, TAncSplit_I, TCont_I, m_I)
# p_I = (nu1_I, nu2_I, nu2Cont_I, TAncSplit_I, TCont_I, m_I, TAncSplitPre1_I, TAncSplitPre2_I, TAncSplitPre3_I, NaPre1_I, NaPre2_I, NaPre3_I)

# msprime demography model with dadi parameters
dem_I = msprime_secondary_contact_testing(Na_I, p_I)
dem_I.sort_events()

graph = dem_I.to_demes()
demes.dump(graph, f"{challenge_name}{extra}/demes_{challenge_name}{extra}.testing.yaml")
demesdraw.tubes(graph, title=f"{challenge_name.replace('_', ' ')}, Testing Challenge")
plt.savefig(f"{challenge_name}{extra}/demes_{challenge_name}{extra}.testing.png", dpi=300)

print("Demes file made")



# Na_F = 89421
# Na_F = 41532 # ancestral pop size
Na_F = 13532
nu1_F = 3.233
nu2_F = 0.09957
nu2Cont_F = 1.986

TAncSplit_F = 17880.54
TCont_F = 2573

m_F = 7.33e-5

TAncSplitPre1_F = 19000
TAncSplitPre2_F = 20000
TAncSplitPre3_F = 21000
NaPre1_F = 80421
NaPre2_F = 50421
NaPre3_F = 120421

popsize1 = nu1_F * Na_F
popsize2 = nu2_F * Na_F
popsize2Cont = nu2Cont_F * Na_F
print(f"Pop-size 1, testing challenge: {popsize1}")
# 289098.093
print(f"Pop-size 2 split, testing challenge: {popsize2}")
# 8903.64897
print(f"Pop-size 2 contact, testing challenge: {popsize2Cont}\n")
# 177590.106

TAncSplit_dadi_F = TAncSplit_F / (2 * Na_F)
TCont_dadi_F = TCont_F / (2 * Na_F)
print(f"dadi ancient split-time, testing challenge: {TAncSplit_dadi_F}")
# 0.09997953500855504
print(f"dadi secondary contact time, testing challenge: {TCont_dadi_F}\n")
# 0.014387000816363047

print(f"T_nu, pop1 ancient split, testing challenge: {TAncSplit_dadi_F / nu1_F}")
# 0.03092469378551037
print(f"T_nu, pop2 ancient split, testing challenge: {TAncSplit_dadi_F / nu2_F}")
# 1.0041130361409565
print(f"T_nu, pop2 secondary contact, testing challenge: {TCont_dadi_F / nu2Cont_F}\n#####\n\n\n#####")
# 0.007244209877322783

p_F = (nu1_F, nu2_F, nu2Cont_F, TAncSplit_F, TCont_F, m_F)
# p_F = (nu1_F, nu2_F, nu2Cont_F, TAncSplit_F, TCont_F, m_F, TAncSplitPre1_F, TAncSplitPre2_F, TAncSplitPre3_F, NaPre1_F, NaPre2_F, NaPre3_F)

# msprime demography model with dadi parameters
dem_F = msprime_secondary_contact_final(Na_F, p_F)
dem_F.sort_events()

graph = dem_F.to_demes()
demes.dump(graph, f"{challenge_name}{extra}/demes_{challenge_name}{extra}.final.yaml")
demesdraw.tubes(graph, title=f"{challenge_name.replace('_', ' ')}, One-Shot Challenge")
plt.savefig(f"{challenge_name}{extra}/demes_{challenge_name}{extra}.final.png", dpi=300)

print("Demes file made")

position_transform = lambda x: 1 + np.array(x)
# position_transform = lambda x: np.fmax(1, x)

dataset = ['testing', 'final']

ts = msprime.sim_ancestry(
    samples=ns, 
    demography=[dem_I, dem_F][process_ii], 
    sequence_length=ast.literal_eval(seq_l), 
    recombination_rate=recomb, 
    ploidy=ploidy,
    random_seed=12
    )

mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=True, random_seed=5566)
mts.num_mutations

vcf_fi = open(f"{challenge_name}{extra}/{challenge_name}{extra}.{dataset[process_ii]}.vcf","w")
vcf_fi.write(mts.as_vcf(position_transform=position_transform))#position_transform=position_transform))
vcf_fi.close()

mts.dump(f"{challenge_name}/{challenge_name}.{dataset[process_ii]}.trees")


# if process_ii == 0:

#     ts
#     #position_transform = lambda x: 1 + x





# for i in range(0,ns['mainland']):
#     fi.write(f"tsk_{i}\tmainland\n")
# for ii in range(i+1, i+ns['island']+1):
#     fi.write(f"tsk_{ii}\tisland\n")
# fi.close()

# for i in range(1,ns['mainland']+1):
#     fi.write(f"mainland_{i}\tmainland\n")
# for ii in range(1, ns['island']+1):
#     fi.write(f"island_{ii}\tisland\n")
# fi.close()


# if process_ii == 1:
#     ts = msprime.sim_ancestry(
#         samples=ns, 
#         demography=dem_F, 
#         sequence_length=ast.literal_eval(seq_l), 
#         recombination_rate=recomb, 
#         ploidy=ploidy,
#         random_seed=12
#         )
#     ts

#     mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=True, random_seed=5566)
#     mts.num_mutations

#     vcf_fi = open(f"{challenge_name}{extra}/{challenge_name}{extra}.final.vcf","w")
#     vcf_fi.write(mts.as_vcf(position_transform=position_transform))
#     vcf_fi.close()

# fi = open(f"{challenge_name}{extra}/{challenge_name}{extra}.popfile.txt","w")

# for i in range(0,ns['mainland']):
#     fi.write(f"tsk_{i}\tmainland\n")
# for ii in range(i+1, i+ns['island']+1):
#     fi.write(f"tsk_{ii}\tisland\n")
# fi.close()

# for i in range(1,ns['mainland']+1):
#     fi.write(f"mainland_{i}\tmainland\n")
# for ii in range(1, ns['island']+1):
#     fi.write(f"island_{ii}\tisland\n")
# fi.close()

# fi = open(f"{challenge_name}{extra}/{challenge_name}{extra}.popfile.txt","w")
# for i in range(ns['mainland']):
#     fi.write(f"tsk_{i}\tmainland\n")
# for i in range(ns['mainland'],ns['mainland']+ns['island']):
#     fi.write(f"tsk_{i}\tisland\n")
# fi.close()

fi = open(f"{challenge_name}{extra}/{challenge_name}{extra}.new_samples.txt","w")
for i in range(ns['mainland']):
    fi.write(f"tsk_{i}\tmainland_{i+1}\n")
for i in range(ns['mainland'],ns['mainland']+ns['island']):
    fi.write(f"tsk_{i}\tisland_{i-ns['mainland']+1}\n")
fi.close()

# if process_ii == 0:
#     os.system(f"bcftools reheader -s {challenge_name}/{challenge_name}.new_samples.txt -o {challenge_name}/{challenge_name}.testing.reheader.vcf {challenge_name}/{challenge_name}.testing.vcf")
#     os.system(f"bgzip -c {challenge_name}/{challenge_name}.testing.reheader.vcf > {challenge_name}/{challenge_name}.testing.reheader.vcf.gz")
#     os.system(f"bcftools index {challenge_name}/{challenge_name}.testing.reheader.vcf.gz")
# if process_ii == 1:
#     os.system(f"bcftools reheader -s {challenge_name}/{challenge_name}.new_samples.txt -o {challenge_name}/{challenge_name}.final.reheader.vcf {challenge_name}/{challenge_name}.final.vcf")
#     os.system(f"bgzip -c {challenge_name}/{challenge_name}.final.reheader.vcf > {challenge_name}/{challenge_name}.final.reheader.vcf.gz")
#     os.system(f"bcftools index {challenge_name}/{challenge_name}.final.reheader.vcf.gz")


# # Example from slack

# def msprime_split_mig(s1, p):
#     (nu1, nu2, T, m) = p
#     dem = msprime.Demography()
#     dem.add_population(name="A", initial_size=s1*10**nu1) # pop1 at present time
#     dem.add_population(name="B", initial_size=s1*10**nu2) # pop2 at present time
#     dem.add_population(name="C", initial_size=s1) # ancestral pop
#     dem.add_population_split(time=2*s1*T, derived=["A", "B"], ancestral="C")
#     dem.set_symmetric_migration_rate(["A", "B"], m/(2*s1))
#     return dem


# nu1 = random.random() * 4 - 2 # nu1 in log scale
# nu2 = random.random() * 4 - 2 # nu2 in log scale
# T = random.random() * 1.9 + 0.1
# m = random.random() * 9 + 1
# p = (nu1, nu2, T, m)
# p

# # msprime demographic and ancestry simulation parameters (unchanged params)
# s1 = 1e4 # ancestral pop size
# ns = {"A":10, "B":10} # sample size
# ploidy = 2 # diploid
# mut = 1e-8 # mutation rate
# # msprime demography params to alter for variance
# seq_l = 1e7
# recomb = 1e-8
# # msprime demography model with dadi parameters
# dem = msprime_split_mig(s1, p)
# dem.debug()


# # generate msprime ts
# ts = msprime.sim_ancestry(samples=ns, ploidy=ploidy, demography=dem, 
#                         sequence_length=seq_l, recombination_rate=recomb)
# # add mutations
# mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=False)
# mts.num_mutations





# # convert mts to afs 
# # specifying sample nodes in the tree sequence for ts-to-afs conversion
# s0 = [node_id for node_id in range(0,20)]
# s1 = [node_id for node_id in range(20,40)]
# sample_nodes = [s0, s1]

# afs = mts.allele_frequency_spectrum(sample_sets=sample_nodes,
#                             polarised=True, span_normalise=False)
# # convert to dadi fs
# fs = dadi.Spectrum(afs)
# # normalize
# fs_norm = fs/fs.sum()
