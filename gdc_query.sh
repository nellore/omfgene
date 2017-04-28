#!/usr/bin/env bash
curl --request POST --header "Content-Type: application/json" --data @gdc_filename_to_id_query.json 'https://gdc-api.nci.nih.gov/files' | gzip -9 >gdc_metadata.tsv.gz
