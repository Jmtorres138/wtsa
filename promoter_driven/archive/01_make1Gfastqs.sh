echo
echo "Making 1G sized fastqs .."
echo
echo
hostname
echo
pwd
echo
echo $0
echo
date
echo
module purge
module list


infolder="/t1-data/user/hugheslab/jtorres/promoter-driven/Fastq_files"

for infile in ${infolder}/*.fastq.gz
do

basename=$(basename $infile | sed 's/.fastq.gz//')
echo $infile
echo $basename
echo $infile >> "/dev/stderr"
echo $basename >> "/dev/stderr"
echo

mkdir ${basename}

echo zcat
echo zcat >> "/dev/stderr"
date

zcat $infile | paste - - - - | awk 'BEGIN {FS="\t"}{print $1"\n"$2"\n"$3"\n"$4 >> "'${basename}/${basename}'_"NR%20".fastq" }'

echo gzip
echo gzip >> "/dev/stderr"
date
c=0
for file in ${basename}/${basename}*.fastq
do
echo ${c}
echo ${c} >> "/dev/stderr"
gzip $file
c=$((${c}+1))
done

done



echo
date
echo 'All done !'
echo


