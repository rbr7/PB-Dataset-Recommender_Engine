while read line
do
   echo "GSE id: $line"
   python3 recommendation_metadata_extraction.py $line
done <$1
echo 'i am done'