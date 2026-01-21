import os

challenges = ["GHIST_2025_bottleneck", "GHIST_2025_secondary_contact", "GHIST_2025_admixture"]
for chal in challenges:
    for chaltype in ["final", "testing"]:
        fi = open(f"{chal}/{chal}.{chaltype}.reheader.vcf").readlines()
        fi_aa = open(f"{chal}/{chal}.{chaltype}.vcf",'w')
        for line in fi:
            if line.startswith("##"):
                fi_aa.write(line)
            elif line.startswith("#CHROM"):
                fi_aa.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele. Format: AA. AA: Ancestral allele\n')
                fi_aa.write(line)
            else:
                fi_aa.write(f"{'\t'.join(line.split('\t')[:7])}\tAA={line.split("\t")[3]}\t{'\t'.join(line.split('\t')[8:])}")
        fi_aa.close()
        os.system(f"bgzip -c {chal}/{chal}.{chaltype}.vcf > {chal}/{chal}.{chaltype}.vcf.gz")
        os.system(f"bcftools index {chal}/{chal}.{chaltype}.vcf.gz")
