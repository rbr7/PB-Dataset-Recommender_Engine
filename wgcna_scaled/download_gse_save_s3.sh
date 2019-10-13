#!/bin/bash
## ./download_gse.sh <GEO search results file path>"
## download files and save to S3 bucket

# path to the file containing new gse_ids
FILE=$1

mkdir GSE_soft/

while read line; do
    echo "GSE id: $line"
    g=$(sed "s/...$/nnn/"<<<$line)
    FTP_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/$g/$line/matrix/"
    wget -r -np -nH --cut-dirs=4 $FTP_URL
    aws s3 cp matrix/ "s3://wgcna-geo-datasets/geo_series_matrix/$line/" --recursive
    rm -r matrix/
    echo "saved to s3 bucket"

done <$FILE

echo "All downloads are now complete!" 

