# running nextflow variant caller/ SNP mapping and genotyping pipeline from malariagen, as they have file for vqsr for anopheles gambiae
# only need to do this for my Bijagos samples as the malariagen samples were genotyped using this pipeline already
# create fastq entry-point and manifest

nextflow /mnt/storage11/sophie/ag1000g_genotyping_nf-v1.0.3/pipelines/snp_genotyping_vector.nf \
-c /mnt/storage11/sophie/ag1000g_genotyping_nf-v1.0.3/pipelines/snp_genotyping_vector.gambiae.config \
-profile standard \
-with-trace \
--s3_bucket mocks3bucket-ag3 \
--input_manifest_fastq /mnt/storage11/sophie/bijagos_mosq_wgs/nextflow_malariagen/pipeline_manifest.tsv \
-resume

# use -resume so that it can continue from where it failed previously

# see how it's going with cat .nextflow.log (in /mnt/storage11/sophie/bijagos_mosq_wgs/nextflow_malariagen)

# java version that appeared when running pipeline: openjdk 17.0.8-internal 2023-07-18

# conda install gatk and nextflow at the same time so that the java files are compatible
