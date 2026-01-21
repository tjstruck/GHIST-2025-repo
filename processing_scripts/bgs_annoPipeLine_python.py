#!/usr/bin/env python

###SBATCH --job-name=hpc_msprime
#SBATCH --job-name=anno-sel-tests
#SBATCH --output=hpc_output/%x-%A.out
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
#SBATCH --time=1:30:00
#SBATCH --array=1-80
###SBATCH --begin=2025-07-06T06:00:00
###SBATCH --array=1-12
###SBATCH --array=1-62

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

import sys,os,glob
import tskit, stdpopsim
import dadi
import matplotlib.pyplot as plt

#tag = sys.argv[1]
#print(tag)
species = stdpopsim.get_species("DroMel")
cds = species.get_annotations("FlyBase_BDGP6.32.51_CDS")

os.makedirs("vcf_to_sfs_files", exist_ok=True)
os.makedirs("SFS", exist_ok=True)

cds_positions_dict = {}

samples = {"AFR": 32}
fi = open(f"vcf_to_sfs_files/popfile.txt","w")
for i in range(samples['AFR']):
    fi.write(f"AFR_{i+1}\tAFR\n")
fi.close()

fi = open(f"vcf_to_sfs_files/new_samples.txt","w")
for i in range(samples['AFR']):
    fi.write(f"tsk_{i}\tAFR_{i+1}\n")
fi.close()

samples = {"AFR": 12, "EUR": 23}
fi = open(f"vcf_to_sfs_files/SplitMig.new_samples.txt","w")
for i in range(samples['AFR']):
    fi.write(f"tsk_{i}\tAFR_{i+1}\n")
for i in range(samples['AFR'],samples['AFR']+samples['EUR']):
    fi.write(f"tsk_{i}\tEUR_{i-samples['AFR']+1}\n")
fi.close()

samples = {"AFR": 12, "EUR": 23}
fi = open(f"vcf_to_sfs_files/SplitMig.popfile.txt","w")
for i in range(samples['AFR']):
    fi.write(f"AFR_{i+1}\tAFR\n")
for i in range(samples['AFR'],samples['AFR']+samples['EUR']):
    fi.write(f"EUR_{i-samples['AFR']+1}\tEUR\n")
fi.close()

chroms = list(range(1, 22)) # Human chromosomes

# ## Uncomment to clean up annotated file
# os.system('rm GHIST_2025-AfricanGrowth_testing_10kPop-shape_0.*/*.rehead*')
# os.system('rm GHIST_2025-AfricanGrowth_testing_10kPop-shape_0.*/*.anno*')
# os.system('rm GHIST_2025-AfricanGrowth_testing_10kPop-shape_0.*/*.syn.vcf')
# os.system('rm GHIST_2025-AfricanGrowth_testing_10kPop-shape_0.*/*.syntrue.vcf')
# os.system('rm GHIST_2025-AfricanGrowth_testing_10kPop-shape_0.*/*.nonsyn.vcf')
# os.system('rm GHIST_2025-AfricanGrowth_testing_10kPop-shape_0.*/*.noncoding.vcf')
# os.system('rm GHIST_2025-AfricanGrowth_testing_10kPop-shape_0.*/*.allneg.vcf')
# os.system('rm GHIST_2025-AfricanGrowth_testing_10kPop-shape_0.*/*.testing.vcf')

# Get intron positions
intron_bed= open("../anno_introns/transcripts.bed", "r").readlines()
import collections
intron_positions = collections.defaultdict(set)
# intron_positions_dict[chrom] = 0
for interval in intron_bed:
    chrom, start, stop = interval.split('\t')[:3]
    intron_positions[chrom].update(range(int(start), int(stop) + 1))
# cds_positions_dict[chrom] = len(cds_positions)

#if tag == 'vcf':
for chrom in ["2L", "2R", "3L", "3R"]:
    
    # fi_bash = open(f"anno.{chrom}.slurm", "w")
    # fi_bash.write("#!/bin/bash\n")
    # fi_bash.write(
    # f"#SBATCH --job-name=anno-{chrom}\n"\
    # "#SBATCH --output=hpc_output/%x-%A.out\n"\
    # "#SBATCH --account=rgutenk\n"\
    # "#SBATCH --qos=user_qos_rgutenk\n"\
    # "#SBATCH --mail-type=ALL\n"\
    # "#SBATCH --mail-user=tjstruck@arizona.edu\n"\
    # "###SBATCH --partition=windfall\n"\
    # "#SBATCH --partition=high_priority\n"\
    # "#SBATCH --nodes=1\n"\
    # "#SBATCH --ntasks=90\n"\
    # "#SBATCH --cpus-per-task=1\n"\
    # "###SBATCH --mem-per-cpu=5gb\n"\
    # "###SBATCH --constraint=hi_mem\n"\
    # "###SBATCH --mem-per-cpu=32gb\n"\
    # "#SBATCH --time=120:00:00\n"\
    # "###SBATCH --array=1-4\n"\
    # "module load bcftools htslib samtools\n"\
    # )

    print(chrom)
    fi = open(f"vcf_to_sfs_files/{chrom}.chr_map.txt","w")
    fi.write(f"1\t{chrom}")
    fi.close()
    # Get CDS positions
    bed_fi = open(f"vcf_to_sfs_files/GHIST_2025_Backgound_Selection_Demography_Challenges_{chrom}.bed","w")
    cds_bed = cds.get_chromosome_annotations(chrom)
    cds_positions = set()
    cds_positions_dict[chrom] = 0
    for interval in cds_bed:
        start, stop = interval
        bed_fi.write(f"{chrom}\t{start}\t{stop}\n")
        cds_positions.update(range(start, stop + 1))
    bed_fi.close()
    cds_positions_dict[chrom] = len(cds_positions)

    #runs = glob.glob(f"GHIST_2025-AfricanGrowth_final_10kPop-shape_*/*.{chrom}.testing.raw.vcf")
    runs = glob.glob(f"GHIST_2025-SplitMig*kPop-shape_0.21-scale_*-recomb_*e-09-1xScale/*.{chrom}.*.raw.vcf")
    runs = glob.glob(f"GHIST_2025-*kPop-shape_*-scale_*-recomb_*-1xScale/*.{chrom}.*.raw.vcf")
    runs = glob.glob(f"GHIST_2025-S*kPop-shape_*-scale_*-recomb_*-1xScale/*.{chrom}.*.raw.vcf")
    runs = glob.glob(f"GHIST_2025-*kPop-shape_*-scale_0-recomb_1e-10-1xScale/*.{chrom}.*.raw.vcf") + glob.glob(f"GHIST_2025-*kPop-shape_*-scale_0.00*-recomb_1e-10-1xScale/*.{chrom}.*.raw.vcf")
    runs = glob.glob(f"GHIST_2025-AfricanGrowth_final_10kPop-shape_0-scale_0-recombScale_0.1-1xScale/*.{chrom}.*.raw.vcf")
    runs = glob.glob(f"GHIST_2025-S*/*.{chrom}.*.raw.vcf")
    print(len(runs))
    # btrt
    runs=[runs[process_ii]]
    for run in runs:
        run_reheader_fi = run.replace('.raw.vcf', '.reheader.vcf')
        run_chrom_fi = run.replace('.raw.vcf', '.anno_chroms.vcf')
        run_treeseq = run.replace('.raw.vcf', '.trees')
        run_anno = run.replace('.raw.vcf', '.anno.vcf')
        #try:
        #    open(run_anno).readlines()
        #    print(f"Already annotated: {run}")
        #    continue
        #except FileNotFoundError:
        #    pass
        if not os.path.exists(run.replace('.raw.vcf', '.trees')):
            print(f"No treeseq: {run}")
            continue
        try:
            _, model, shape, scale, recomb, scaling = run.split('/')[0].split('-')
        except ValueError:
            _, model, shape, scale, recomb1, recomb2, scaling = run.split('/')[0].split('-')
            recomb = f"{recomb1}-{recomb2}"
        dd_key = f"{model}_{shape}_{scale}_{recomb}_{scaling}"

        # if "AfricanGrowth" in run or "African3Epoch_1S16" in run or "stdpopsim_dromel_growth" in run:
        if "AfricanGrowth" in run:
            sample_fi = "vcf_to_sfs_files/new_samples.txt"
        else:
            sample_fi = "vcf_to_sfs_files/SplitMig.new_samples.txt"
        chrom_map_fi = f"vcf_to_sfs_files/{chrom}.chr_map.txt"
        # else:
        #     continue

        try:
            os.system(f"bcftools reheader -s {sample_fi} -o {run_reheader_fi} {run}")
            os.system(f"bcftools annotate --rename-chrs {chrom_map_fi} -o {run_chrom_fi} {run_reheader_fi}")
            # fi_bash.write(f"bcftools reheader -s {sample_fi} -o {run_reheader_fi} {run}\n")
            # fi_bash.write(f"bcftools annotate --rename-chrs {chrom_map_fi} -o {run_chrom_fi} {run_reheader_fi}\n")
        
            # This line needs to be moved left to avoid errors and properly close the file
            # fi_bash.close()
        
        except:
            
            print(f"Error processing {run}")

        # Load the trees
        ts = tskit.load(f"{run_treeseq}")
        try:
            fi = open(f"{run_chrom_fi}").readlines()
        except FileNotFoundError:
            print(f"File not found: {run_chrom_fi}")
            continue
        fi_anno = open(f"{run_anno}",'w')

        for line in fi:
            # Fill in the header
            if line.startswith("##"):
                fi_anno.write(line)
            # Add the AA field to the header
            elif line.startswith("#CHROM"):
                fi_anno.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele. Format: AA. AA: Ancestral allele\n')
                fi_anno.write(line)
            else:
                if "-shape_0-scale_0-" in run:
                    fi_anno.write(f"{'\t'.join(line.split('\t')[:7])}\tAA={line.split('\t')[3]}\t{'\t'.join(line.split('\t')[8:])}")
                else:
                    # ts.site() takes the ID (3rd element) in the VCF which corrosponds to the site number
                    site = int(line.split("\t")[2])
                    position = int(line.split("\t")[1])
                    mut = ts.site(site).mutations[0].metadata['mutation_list'][0]['mutation_type']
                    sel = ts.site(site).mutations[0].metadata['mutation_list'][0]['selection_coeff']
                    # Check if the site is in the CDS positions
                    if position in cds_positions:
                        if sel == 0:
                            mutType = "Synonymous"
                        elif sel != 0:
                            mutType = "Nonsynonymous"
                    elif position in intron_positions[chrom]:
                        mutType = "Intronic"
                        #print(mutType)
                    else:
                        mutType = "Intergenic"
                    # Add the selection coefficient to the INFO field
                    # fi_anno.write(f"{'\t'.join(line.split('\t')[:7])}\tAA={line.split('\t')[3]};Variant_Type={mutType};SLiM_mut_ID={mut};s={sel}\t{'\t'.join(line.split('\t')[8:])}")
                    fi_anno.write(f"{'\t'.join(line.split('\t')[:7])}\tAA={line.split('\t')[3]};Variant_type={mutType}\t{'\t'.join(line.split('\t')[8:])}")
                    # fi_anno.write(f"{'\t'.join(line.split('\t')[:7])}\tAA={line.split('\t')[3]};Variant_Type={mutType};s={sel}\t{'\t'.join(line.split('\t')[8:])}")
        fi_anno.close()
