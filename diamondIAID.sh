#!/bin/bash
# HEADER - Do Not Modify!
set -e
shopt -s expand_aliases
export LC_ALL=C
########
export PATH="/public/software/apps/bin:/public/software/apps/bin/singularity:/public/software/apps/go/bin:$PATH"           
alias fastqc="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif fastqc"
alias fastp="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif fastp"
alias bowtie2-build="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif bowtie2-build"
alias bowtie2="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif bowtie2"
alias diamond="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif diamond"

level=S
 # HEADER END
#############you can input your fasta
############ensure the diamond db
############ if you modify the database, please change the IAID all in the script



 #######makediamond db#######

mkdir -p /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID
cd /home/stu_2/bozhang/databaseIAID
for i in `ls ./ | grep .fasta`
do
a=${i%.fasta};
diamond makedb --in ./${a}.fasta -d /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID/${a};
done


##########query#########

cp  /home/stu_2/bozhang/P3.Assembly_annotation/04.cdhit/unigene_pep.fasta /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID/
cd  /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID
for i in  `ls ./ | grep .dmnd`
do
 diamond blastp  \
    --db ${i%.dmnd}.dmnd \
    --query  unigene_pep.fasta \
    --out ${i%.dmnd}_unigene.txt \
    --threads 30 \
    --outfmt 6 \
    --max-target-seqs 1 \
    --evalue 1e-5
done

###########filter identity 40, query length 70###########

cd  /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID
for i in `ls ./ | grep _unigene.txt`
do
echo $i
awk '$3 > 40' ${i} > ./${i% _unigene.txt}_new_identity40.txt
awk '$4 > 70' ./${i% _unigene.txt}_new_identity40.txt > ./${i% _unigene.txt}_new_identity_querylength.txt
##############awk '$11 > 70 ' ${i} > ./${i% _unigene.txt}_new_identity_querylength_bitscore.txt
done

#######################################calculate#######################################

cd  /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID
for i in `ls ./ | grep _new_identity_querylength.txt`
do
mkdir -p /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID/"${i%.txt_new_identity_querylength.txt}"
done

##############double loop through ##########

cd  /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID
for i in ls ./*
do 
if [ -d ${i} ]
then
cp ${i}.txt_new_identity_querylength.txt ./${i}/DATA.txt
cd ${i} || continue
for j in `ls /home/stu_2/bozhang/P3.Assembly_annotation/10.unigene/quants/ | grep quant`
do
awk 'NR==FNR{a[$1]; next} $1 in a'  ./DATA.txt /home/stu_2/bozhang/P3.Assembly_annotation/10.unigene/quants/${j} > ./temp_file.txt
awk '$5 != 0' ./temp_file.txt > ./temp_file2.txt && mv temp_file2.txt ./${j%.quant.sf}quant_result.txt
rm ./temp_file.txt
done
cd ..
fi
done


cd  /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID
for i in ls ./*
do 
if [ -d ${i} ]
then
cd ${i} || continue
for j in `ls ./ | grep quant`
do
awk '{ sum += $5 } END { print sum }' ./${j}  > ./${j%quant_result.txt}_file.txt

 sed "s/^/$(basename "${j%quant_result.txt}" _file.txt)\t/" "${j%quant_result.txt}_file.txt" > "${j%quant_result.txt}_updated.txt"

cat  ./*updated.txt> ${i}_allsample_results.txt
done
cd ..
fi
done


