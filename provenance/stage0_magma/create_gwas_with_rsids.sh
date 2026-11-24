#!/bin/bash
# Replace custom SNP IDs with rsIDs by matching chromosome:position

GWAS_IN="$1"
REF_BIM="$2"
GWAS_OUT="$3"

echo "Creating position-to-rsID lookup..."
# Create lookup: chr:pos -> rsID
awk '{print $1":"$4, $2}' "${REF_BIM}" > /tmp/pos_to_rsid.txt

echo "Replacing SNP IDs in GWAS file..."
# Match GWAS SNPs by position and replace ID with rsID
# Output: rsID chr pos pval (needed for both annotation and analysis)
awk 'BEGIN{OFS="\t"}
NR==FNR {
    # Read lookup table
    lookup[$1] = $2
    next
}
{
    # Process GWAS file
    snp_id = $1
    chr = $2
    pos = $3
    pval = $4

    # Create position key
    pos_key = chr ":" pos

    # If we find a matching rsID, use it; otherwise skip
    if (pos_key in lookup) {
        print lookup[pos_key], chr, pos, pval
    }
}' /tmp/pos_to_rsid.txt "${GWAS_IN}" > "${GWAS_OUT}"

n_out=$(wc -l < "${GWAS_OUT}")
echo "Done. Output: ${GWAS_OUT} (${n_out} SNPs with rsIDs)"
