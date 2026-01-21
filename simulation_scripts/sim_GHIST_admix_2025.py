#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=sim-admix
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
###SBATCH --array=2
import dadi
import msprime
import demes, demesdraw
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


challenge_name = "GHIST_2025_admixture"
os.makedirs(challenge_name, exist_ok=True)


extra = ''
# seq_l = ["2.5e8"][process_ii]
seq_l = "2.5e8"

def msprime_split_isolation_then_migtation_I(Na, p):
    (nuExt, nuAnc, nuMod1_1, nuMod1_2, nuMod2_1, nuMod2_2, TSplitAnc, TSplitMod, TAdmix1, TAdmix2, TFinal, m, factor, pulse_size1, pulse_size2) = p
    dem = msprime.Demography()
    dem.add_population(name="common_ancestor", initial_size=Na/factor, initially_active=True)
    dem.add_population(name="Modern_pop1", initial_size=nuMod1_2/factor)
    dem.add_population(name="Modern_pop2", initial_size=nuMod2_2/factor)
    dem.add_population(name="AncOut1", initial_size=nuExt/factor, default_sampling_time=4523)
    dem.add_population(name="AncOut2", initial_size=nuExt/factor, default_sampling_time=3897)
    dem.add_population(name="AncPop1", initial_size=nuExt/factor, default_sampling_time=2984)
    dem.add_population(name="AncOut3", initial_size=nuExt/factor, default_sampling_time=5123)
    dem.add_population(name="AncPop2", initial_size=nuExt/factor, default_sampling_time=2765)
    dem.add_population(name="ancestral", initial_size=Na/factor)
    dem.add_population(name="ancient_population", initial_size=nuAnc/factor)

    # dem.add_population_split(time=TSplitAnc*1.8/factor, derived=["ancestral", "bad_sample"], ancestral="common")
    # dem.add_population_split(time=TSplitAnc*1.8/factor, derived=["ancestral"], ancestral="common")

    dem.add_population_parameters_change(time=TSplitAnc*1.5/factor, population="ancestral", initial_size=Na/factor*.5)
    dem.add_population_parameters_change(time=TSplitAnc*1.2/factor, population="ancestral", initial_size=Na/factor*1.5)
    dem.add_population_parameters_change(time=TSplitAnc/factor, population="ancestral", initial_size=Na/factor*2)

    # dem.add_population_split(time=TSplitAnc*1.1/factor, derived=["bad_sample"], ancestral="ancestral")

    # dem.add_population_split(time=TSplitAnc/factor, derived=["extinct_pop", "ancient_population", "bad_sample"], ancestral="ancestral")
    dem.add_population_split(time=TSplitAnc/factor, derived=["common_ancestor", "ancient_population"], ancestral="ancestral")

    dem.add_population_split(time=TSplitMod/factor, derived=["Modern_pop1", "Modern_pop2"], ancestral="ancient_population")

    dem.add_population_split(time=9500.12, derived=["AncOut1"], ancestral="common_ancestor")
    dem.add_population_split(time=8732.89, derived=["AncOut2"], ancestral="common_ancestor")
    dem.add_population_split(time=12500.78, derived=["AncPop1"], ancestral="common_ancestor")
    dem.add_population_split(time=11200.56, derived=["AncOut3"], ancestral="common_ancestor")
    dem.add_population_split(time=14200.34, derived=["AncPop2"], ancestral="common_ancestor")

    dem.add_mass_migration(time=TAdmix1/factor, source="Modern_pop1", dest="AncPop1", proportion=pulse_size1/factor)
    dem.add_mass_migration(time=TAdmix2/factor, source="Modern_pop2", dest="AncPop2", proportion=pulse_size2/factor)

    # dem.add_population_parameters_change(TFinal, initial_size=0, population="extinct_population")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="Modern_pop1", dest="Modern_pop2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="Modern_pop2", dest="Modern_pop1")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop1", dest="AncPop2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop2", dest="AncPop1")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut3", dest="AncPop2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop2", dest="AncOut3")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut1", dest="AncPop1")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop1", dest="AncOut1")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop1", dest="AncOut2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut2", dest="AncPop1")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut1", dest="AncOut2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut2", dest="AncOut1")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut1", dest="AncOut3")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut3", dest="AncOut1")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut2", dest="AncOut3")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut3", dest="AncOut2")

    dem.add_population_parameters_change(time=TAdmix1*.8, population="Modern_pop1", initial_size=nuMod1_1)
    dem.add_population_parameters_change(time=TAdmix2*.8, population="Modern_pop2", initial_size=nuMod2_1)

    # dem.add_migration_rate_change(time=Na*T, rate=0, source="modern_1", dest="modern_2")
    # dem.add_migration_rate_change(time=Na*T, rate=0, source="modern_2", dest="modern_1")
    # dem.add_migration_rate_change(time=0, rate=m/(2*Na), source="modern_1", dest="modern_2")
    # dem.add_migration_rate_change(time=0, rate=m/(2*Na), source="modern_2", dest="modern_1")
    return dem

def msprime_split_isolation_then_migtation_F(Na, p):
    (nuExt, nuAnc, nuMod1_1, nuMod1_2, nuMod2_1, nuMod2_2, TSplitAnc, TSplitMod, TAdmix1, TAdmix2, TFinal, m, factor, pulse_size1, pulse_size2) = p
    dem = msprime.Demography()
    dem.add_population(name="common_ancestor", initial_size=Na/factor, initially_active=True)
    dem.add_population(name="Modern_pop1", initial_size=nuMod1_2/factor)
    dem.add_population(name="Modern_pop2", initial_size=nuMod2_2/factor)
    dem.add_population(name="AncOut1", initial_size=nuExt/factor, default_sampling_time=4789)
    dem.add_population(name="AncOut2", initial_size=nuExt/factor, default_sampling_time=3654)
    dem.add_population(name="AncPop1", initial_size=nuExt/factor, default_sampling_time=3123)
    dem.add_population(name="AncOut3", initial_size=nuExt/factor, default_sampling_time=5342)
    dem.add_population(name="AncPop2", initial_size=nuExt/factor, default_sampling_time=2890)
    dem.add_population(name="ancestral", initial_size=Na/factor)
    dem.add_population(name="ancient_population", initial_size=nuAnc/factor)

    # dem.add_population_split(time=TSplitAnc*1.8/factor, derived=["ancestral", "bad_sample"], ancestral="common")
    # dem.add_population_split(time=TSplitAnc*1.8/factor, derived=["ancestral"], ancestral="common")

    dem.add_population_parameters_change(time=TSplitAnc*1.5/factor, population="ancestral", initial_size=Na/factor*.5)
    dem.add_population_parameters_change(time=TSplitAnc*1.2/factor, population="ancestral", initial_size=Na/factor*1.5)
    dem.add_population_parameters_change(time=TSplitAnc/factor, population="ancestral", initial_size=Na/factor*2)

    # dem.add_population_split(time=TSplitAnc*1.1/factor, derived=["bad_sample"], ancestral="ancestral")

    # dem.add_population_split(time=TSplitAnc/factor, derived=["extinct_pop", "ancient_population", "bad_sample"], ancestral="ancestral")
    dem.add_population_split(time=TSplitAnc/factor, derived=["common_ancestor", "ancient_population"], ancestral="ancestral")

    dem.add_population_split(time=TSplitMod/factor, derived=["Modern_pop1", "Modern_pop2"], ancestral="ancient_population")

    dem.add_population_split(time=7950.89, derived=["AncOut1"], ancestral="common_ancestor")
    dem.add_population_split(time=9950.34, derived=["AncOut2"], ancestral="common_ancestor")
    dem.add_population_split(time=12950.23, derived=["AncPop1"], ancestral="common_ancestor")
    dem.add_population_split(time=10875.67, derived=["AncOut3"], ancestral="common_ancestor")
    dem.add_population_split(time=6850.12, derived=["AncPop2"], ancestral="common_ancestor")


    dem.add_mass_migration(time=TAdmix1/factor, source="Modern_pop1", dest="AncPop1", proportion=pulse_size1/factor)
    # Replacing AncPop2 with AncOut2 to mix up the admixed population
    dem.add_mass_migration(time=TAdmix2/factor, source="Modern_pop2", dest="AncOut2", proportion=pulse_size2/factor)

    # dem.add_population_parameters_change(TFinal, initial_size=0, population="extinct_population")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="Modern_pop1", dest="Modern_pop2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="Modern_pop2", dest="Modern_pop1")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop1", dest="AncPop2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop2", dest="AncPop1")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut3", dest="AncPop2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop2", dest="AncOut3")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut1", dest="AncPop1")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop1", dest="AncOut1")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncPop1", dest="AncOut2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut2", dest="AncPop1")

    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut1", dest="AncOut2")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut2", dest="AncOut1")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut1", dest="AncOut3")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut3", dest="AncOut1")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut2", dest="AncOut3")
    dem.add_migration_rate_change(time=0, rate=m/factor, source="AncOut3", dest="AncOut2")

    dem.add_population_parameters_change(time=TAdmix1*.8, population="Modern_pop1", initial_size=nuMod1_1)
    dem.add_population_parameters_change(time=TAdmix2*.8, population="Modern_pop2", initial_size=nuMod2_1)

    # dem.add_migration_rate_change(time=Na*T, rate=0, source="modern_1", dest="modern_2")
    # dem.add_migration_rate_change(time=Na*T, rate=0, source="modern_2", dest="modern_1")
    # dem.add_migration_rate_change(time=0, rate=m/(2*Na), source="modern_1", dest="modern_2")
    # dem.add_migration_rate_change(time=0, rate=m/(2*Na), source="modern_2", dest="modern_1")
    return dem



# factor = 1
# Na = 32671 # ancestral pop size
# nuExt = 13249
# nuAnc = 5083

# nuMod1_1 = 2231.0
# nuMod1_2 = 9025.0

# nuMod2_1 = 1293.0
# nuMod2_2 = 6962.0

# TSplitAnc = 15090.0 # gen
# TSplitMod = 1758.0 # gen

# TAdmix1 = 1566
# TAdmix2 = 883

# pulse_size1 = 0.011
# pulse_size2 = 0.002

# TFinal = 0.1
# m = 3.14e-05

# p = (nuExt, nuAnc, nuMod1_1, nuMod1_2, nuMod2_1, nuMod2_2, TSplitAnc, TSplitMod, TAdmix1, TAdmix2, TFinal, m, factor)

factor_I = 1
Na_I = 33012.7
nuExt_I = 13501.3
nuAnc_I = 5102.8

nuMod1_1_I = 2250.45
nuMod1_2_I = 9100.67

nuMod2_1_I = 1323.89
nuMod2_2_I = 7023.34

TSplitAnc_I = 15100.56
TSplitMod_I = 1761.78

TAdmix1_I = 1270.12
TAdmix2_I = 913.45

pulse_size1_I = 0.0123
pulse_size2_I = 0.0034

TFinal_I = 0.152
m_I = 3.215e-05

p_I = (nuExt_I, nuAnc_I, nuMod1_1_I, nuMod1_2_I, nuMod2_1_I, nuMod2_2_I, TSplitAnc_I, TSplitMod_I, TAdmix1_I, TAdmix2_I, TFinal_I, m_I, factor_I, pulse_size1_I, pulse_size2_I)

factor_F = 1
Na_F = 31987.4
nuExt_F = 12999.6
nuAnc_F = 5048.9

nuMod1_1_F = 2199.78
nuMod1_2_F = 8949.56

nuMod2_1_F = 1279.34
nuMod2_2_F = 6899.12

TSplitAnc_F = 15079.89
TSplitMod_F = 1754.67

TAdmix1_F = 1159.45
TAdmix2_F = 823.78

pulse_size1_F = 0.0202
pulse_size2_F = 0.0018

TFinal_F = 0.048
m_F = 3.105e-05

p_F = (nuExt_F, nuAnc_F, nuMod1_1_F, nuMod1_2_F, nuMod2_1_F, nuMod2_2_F, TSplitAnc_F, TSplitMod_F, TAdmix1_F, TAdmix2_F, TFinal_F, m_F, factor_F, pulse_size1_F, pulse_size2_F)

ploidy = 2 # diploid
mut = 1.29e-8 # mutation rate

# 1e-8 recombination per base seems roughtly right 0.0099±0.0052 and 0.0088±0.0053 per 1MB in two species of cattle
recomb = 1.38e-08
# msprime demography model with dadi parameters

dem_I = msprime_split_isolation_then_migtation_I(Na_I, p_I)
dem_I.sort_events()
graph_I = dem_I.to_demes()
demes.dump(graph_I, f"{challenge_name}/demes_{challenge_name}.testing.yaml")
demesdraw.tubes(graph_I, title=f"{challenge_name.replace('_', ' ')}, Testing Challenge")
plt.savefig(f"{challenge_name}{extra}/demes_{challenge_name}{extra}.testing.png", dpi=300)
print("Demes file made")

dem_F = msprime_split_isolation_then_migtation_F(Na_F, p_F)
dem_F.sort_events()
graph_F = dem_F.to_demes()
demes.dump(graph_F, f"{challenge_name}/demes_{challenge_name}.final.yaml")
demesdraw.tubes(graph_F)
plt.savefig(f"{challenge_name}{extra}/demes_{challenge_name}{extra}.final.png", dpi=300)
print("Demes file made")

ns = {"Modern_pop1":21, "Modern_pop2":18, "AncOut1":4, "AncOut2":3, "AncPop1": 5, "AncOut3":3, "AncPop2": 2} # individuals

fi = open(f"{challenge_name}/{challenge_name}.popfile.txt","w")
for i in range(1,ns['Modern_pop1']+1):
    fi.write(f"Modern1.{i}\tModern1\n")
for i in range(1,ns['Modern_pop2']+1):
    fi.write(f"Modern2.{i}\tModern2\n")
for i in range(1,ns['AncOut1']+1):
    fi.write(f"AncSite1.{i}\tAnc1\n")
for i in range(1,ns['AncOut2']+1):
    fi.write(f"AncSite2.{i}\tAnc2\n")
for i in range(1,ns['AncPop1']+1):
    fi.write(f"AncSite3.{i}\tAnc3\n")
for i in range(1,ns['AncOut3']+1):
    fi.write(f"AncSite4.{i}\tAnc4\n")
for i in range(1,ns['AncPop2']+1):
    fi.write(f"AncSite5.{i}\tAnc5\n")
fi.close()

fi = open(f"{challenge_name}/{challenge_name}.new_samples.txt","w")
# for i in range(ns['mainland']):
#     fi.write(f"tsk_{i}\tmainland_{i+1}\n")
# for i in range(ns['mainland'],ns['mainland']+ns['island']):
#     fi.write(f"tsk_{i}\tisland_{i-ns['mainland']+1}\n")

for i in range(ns['Modern_pop1']):
    fi.write(f"tsk_{i}\tModern1.{i+1}\n")

for i in range(ns['Modern_pop1'],ns['Modern_pop2']+ns['Modern_pop1']):
    fi.write(f"tsk_{i}\tModern2.{i+1-ns['Modern_pop1']}\n")

for i in range(ns['Modern_pop2']+ns['Modern_pop1'],ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1']):
    fi.write(f"tsk_{i}\tAncSite1.{i+1-(ns['Modern_pop2']+ns['Modern_pop1'])}\n")

for i in range(ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1'],ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1']):
    fi.write(f"tsk_{i}\tAncSite2.{i+1-(ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1'])}\n")

for i in range(ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1'],ns['AncPop1']+ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1']):
    fi.write(f"tsk_{i}\tAncSite3.{i+1-(ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1'])}\n")

for i in range(ns['AncPop1']+ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1'],ns['AncOut3']+ns['AncPop1']+ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1']):
    fi.write(f"tsk_{i}\tAncSite4.{i+1-(ns['AncPop1']+ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1'])}\n")

for i in range(ns['AncOut3']+ns['AncPop1']+ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1'],ns['AncPop2']+ns['AncOut3']+ns['AncPop1']+ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1']):
    fi.write(f"tsk_{i}\tAncSite5.{i+1-(ns['AncOut3']+ns['AncPop1']+ns['AncOut2']+ns['AncOut1']+ns['Modern_pop2']+ns['Modern_pop1'])}\n")
fi.close()

dataset = ['testing', 'final']

print('treeseq')
ts = msprime.sim_ancestry(
    samples=ns, 
    demography=[dem_I, dem_F][process_ii], 
    sequence_length=ast.literal_eval(seq_l), 
    recombination_rate=recomb, 
    ploidy=ploidy,
    random_seed=12
    )
print("TS sim done")
mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=True, random_seed=5566)
mts.num_mutations

vcf_fi = open(f"{challenge_name}/{challenge_name}.{dataset[process_ii]}.raw.vcf","w")
vcf_fi.write(mts.as_vcf())
vcf_fi.close()

mts.dump(f"{challenge_name}/{challenge_name}.{dataset[process_ii]}.trees")

# if process_ii == 0:
    


#     os.system(f"bcftools reheader -s {challenge_name}/{challenge_name}.new_samples.txt -o {challenge_name}/{challenge_name}.testing.reheader.vcf {challenge_name}/{challenge_name}.testing.vcf")
#     os.system(f"bgzip -c {challenge_name}/{challenge_name}.testing.reheader.vcf > {challenge_name}/{challenge_name}.testing.reheader.vcf.gz")
#     os.system(f"bcftools index {challenge_name}/{challenge_name}.testing.reheader.vcf.gz")

# if process_ii == 1:
    
#     ts = msprime.sim_ancestry(
#         samples=ns, 
#         demography=dem_F, 
#         sequence_length=ast.literal_eval(seq_l), 
#         recombination_rate=recomb, 
#         ploidy=ploidy,
#         random_seed=12
#         )
#     print("TS sim done")
#     mts = msprime.sim_mutations(ts, rate=mut, discrete_genome=False, random_seed=5566)
#     mts.num_mutations

#     vcf_fi = open(f"{challenge_name}/{challenge_name}.final.raw.vcf","w")
#     vcf_fi.write(mts.as_vcf())
#     vcf_fi.close()

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
