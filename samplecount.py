import gzip

vcf_file_path = "your_file.vcf.gz"

with gzip.open(vcf_file_path, "rt") as vcf_file:
    for line in vcf_file:
        if line.startswith("#CHROM"):
            samples = line.strip().split("\t")[9:]
            sample_count = len(samples)
            print(f"Number of samples: {sample_count}")
            break
