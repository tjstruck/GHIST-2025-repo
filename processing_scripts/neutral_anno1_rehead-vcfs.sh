
###############################
####################### MSPRIME
#################################


# mv GHIST_2025_bottleneck/GHIST_2025_bottleneck.testing.vcf GHIST_2025_bottleneck/GHIST_2025_bottleneck.testing.raw.vcf
bcftools reheader -s GHIST_2025_bottleneck/GHIST_2025_bottleneck.new_samples.txt \
-o GHIST_2025_bottleneck/GHIST_2025_bottleneck.testing.reheader.vcf \
GHIST_2025_bottleneck/GHIST_2025_bottleneck.testing.raw.vcf

# mv GHIST_2025_bottleneck/GHIST_2025_bottleneck.final.vcf GHIST_2025_bottleneck/GHIST_2025_bottleneck.final.raw.vcf
bcftools reheader -s GHIST_2025_bottleneck/GHIST_2025_bottleneck.new_samples.txt \
-o GHIST_2025_bottleneck/GHIST_2025_bottleneck.final.reheader.vcf \
GHIST_2025_bottleneck/GHIST_2025_bottleneck.final.raw.vcf


# mv GHIST_2025_admixture/GHIST_2025_admixture.testing.vcf GHIST_2025_admixture/GHIST_2025_admixture.testing.raw.vcf
bcftools reheader -s GHIST_2025_admixture/GHIST_2025_admixture.new_samples.txt \
-o GHIST_2025_admixture/GHIST_2025_admixture.testing.reheader.vcf \
GHIST_2025_admixture/GHIST_2025_admixture.testing.raw.vcf

# mv GHIST_2025_admixture/GHIST_2025_admixture.final.vcf GHIST_2025_admixture/GHIST_2025_admixture.final.raw.vcf
bcftools reheader -s GHIST_2025_admixture/GHIST_2025_admixture.new_samples.txt \
-o GHIST_2025_admixture/GHIST_2025_admixture.final.reheader.vcf \
GHIST_2025_admixture/GHIST_2025_admixture.final.raw.vcf


# mv GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.testing.vcf GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.testing.raw.vcf
bcftools reheader -s GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.new_samples.txt \
-o GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.testing.reheader.vcf \
GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.testing.raw.vcf

# mv GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.final.vcf GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.final.raw.vcf
bcftools reheader -s GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.new_samples.txt \
-o GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.final.reheader.vcf \
GHIST_2025_secondary_contact/GHIST_2025_secondary_contact.final.raw.vcf

