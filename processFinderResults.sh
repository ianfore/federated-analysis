#!/bin/bash

if [ $# -ne 1 ]
then
	echo "usage $0 <chrom>"
	exit 1
fi
CHR=$1

# find total num vus in a clean file
cat $CHR-out.json | python -m json.tool > $CHR-pretty.json

grep '(' $CHR-pretty.json > ${CHR}.txt

sed -i 's/\ \ \ \ \ \ \ \ \"//' ${CHR}.txt

sed -i 's/\".*//' ${CHR}.txt

cat ${CHR}.txt | sort | uniq > $CHR-uniq.txt

numVUS=$(wc -l $CHR-uniq.txt|awk '{print $1}')

# find total num coocs
sed -n '/cooccurring/,/homozygous/p' $CHR-pretty.json > $CHR-coocs.json

grep '(' $CHR-coocs.json > $CHR-coocs.txt

sed -i 's/\ \ \ \ \ \ \ \ \"//' $CHR-coocs.txt

sed -i 's/\".*//' $CHR-coocs.txt

cat $CHR-coocs.txt | sort > $CHR-coocs-sorted.txt

numCoocs=$(wc -l $CHR-coocs-sorted.txt|awk '{print $1}')

# find total num homo
sed -n '/homozygous/,$ p' $CHR-pretty.json > $CHR-homo.json

grep '(' $CHR-homo.json > $CHR-homo.txt

sed -i 's/\ \ \ \ \ \ \ \ \"//' $CHR-homo.txt

sed -i 's/\".*//' $CHR-homo.txt

cat $CHR-homo.txt | sort > $CHR-homo-sorted.txt

numHomo=$(wc -l $CHR-homo-sorted.txt|awk '{print $1}')

# find unique and intersecting counts
# comm -123 (1 = suppress unique to f1, 2 => suppress unique to f2, 3=>suppress in both)

numOnlyCooc=$(comm -23 $CHR-coocs-sorted.txt $CHR-homo-sorted.txt | wc -l)

numOnlyHomo=$(comm -13 $CHR-coocs-sorted.txt $CHR-homo-sorted.txt | wc -l)

((numIntersection=$numVUS - $numOnlyCooc - $numOnlyHomo))

echo "total VUS = $numVUS"
echo "num only cooccuring = $numOnlyCooc"
echo "num only homozygous = $numOnlyHomo"
echo "num of intersection = $numIntersection"

# cleanup files
rm $CHR-coocs-sorted.txt
rm $CHR-homo-sorted.txt
rm $CHR-homo.txt
rm $CHR-coocs.txt
rm $CHR-homo.json
rm $CHR-coocs.json
rm $CHR-pretty.json
rm $CHR-uniq.txt
rm ${CHR}.txt
