# prepare locations.txt file
# 3R,2194206,2195205

while IFS=',' read -r chr start end; do
    echo "Processing: $chr from $start to $end" # Debugging line
    awk -v chr="$chr" -v start="$start" -v end="$end" \
        '$1 == chr && $4 <= end && $5 >= start' \
        /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3 >> X_fst_genes.txt
done < X_fst_locations.txt
