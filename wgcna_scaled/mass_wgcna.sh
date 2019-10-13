#!/bin/bash
echo "------------Executing MASS WGCNA---------------"

current_time=$(date "+%d.%m.%Y-%H.%M")
log_fileName="mass_run_"$current_time.log
touch $log_fileName



while getopts "t:g:f:" opt; do
  case $opt in 
    t) tod="$OPTARG" ;;
    g) gse_file="$OPTARG" ;;
    f) master_folder_path="$OPTARG" ;;
  esac
done

# type of dataset
echo "$tod"

if [ -f "module-info.txt" ]
then
  rm module-info.txt
fi

touch module-info.txt

if [ -d "modules/" ]
then 
  rm -r modules/
fi
mkdir modules/

if [ -d "dataset_logs/" ]
then
  rm -r dataset_logs/
fi

if [ -d "matrix" ]
then 
  rm -r matrix/
fi


if [ -d "soft/" ]
then 
  rm -r soft/
fi


mkdir dataset_logs/

# save file with ids into s3 folder
aws s3 cp $gse_file "s3://wgcna-geo-datasets/$master_folder_path/run_$current_time/$gse_file" 


while read line; do
    echo "GSE id: $line"
    #download the dataset matrix
    acc=$(grep -o "[0-9]*" <<< $line)
    len=$(echo $acc | wc -c)
    if [ $len -gt 3 ]
    then
        g=$(sed "s/...$/nnn/"<<<$line)
    else
        g="GSEnnn"
    fi
    FTP_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/$g/$line/matrix/"
    wget -r -np -nH --cut-dirs=4 $FTP_URL
    # run the pipeline
    Rscript run.R $line dataset_logs/ $log_fileName $tod
    #save found module to s3
    if [ -f "modules/$line.gmt" ]
    then
      aws s3 cp "modules/$line.gmt" "s3://wgcna-geo-datasets/$master_folder_path/run_$current_time/WGCNA_modules/$line.gmt"
      echo "module sent to s3 bucket"
    fi
    if [ -s "module-info.txt" ]
    then 
      aws s3 cp "module-info.txt" "s3://wgcna-geo-datasets/$master_folder_path/run_$current_time/module-info.txt"
      echo "module information sent to s3 bucket"
    fi
    #send dataset log to s3
    aws s3 cp "dataset_logs/$line.log" "s3://wgcna-geo-datasets/$master_folder_path/run_$current_time/dataset_logs/$line.log"
    # send complete run logs to s3
    aws s3 cp $log_fileName "s3://wgcna-geo-datasets/$master_folder_path/run_$current_time/$log_fileName"
    # delete the downloaded dataset and platform files
    rm -r matrix/
    rm -r soft/
done <$gse_file


echo "Bye bye"
