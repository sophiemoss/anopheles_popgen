# rename 's/^SRR/SAMN0176022_SRR/' SRR*
# rename 's/^SAMN0176022/SAMN01760622/' SAMN0176022*

# to bgzip multiple files at once using one bash command:
find /mnt/storage11/sophie/bijagos_mosq_wgs/globalmelasfq2vcf/SAMN01760622_multipleSRRs -type f -name "*.fastq" -print0 | xargs -0 -P 4 -n 1 bgzip