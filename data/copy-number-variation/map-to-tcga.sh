#!/usr/bin/env bash
for ID in $(head -n 1 KIRC.focal_score_by_genes.tsv | cut -d $'\t' -f4- | tr $'\t' $'\n'); do;
    echo $ID,$(curl "https://api.gdc.cancer.gov/v0/all?query=${ID}&size=5" | fx '.data.query.hits[0].submitter_id') >> tcga-mapping.csv;
done
