#!/bin/bash
# HEADER - Do Not Modify!
set -e
shopt -s expand_aliases
export LC_ALL=C
########
alias diamond="singularity exec /home/stu_2/linguopeng/sif/run_dbcan_latest.sif  diamond"
level=S
# HEADER END

############db-pwy /home/stu_2/bozhang/P3.Assembly_annotation/11.cazy/db
#############ensure you had made 10.unigene
############# please choose the specific CAZymes and print ctrl + F in keyborad to mofify all the pwy and results
############ personal code 
############ use please reference Bo Zhang1996 in Github
############ if the # more than 5, please do not use the sentence

cd /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID/iad_unigene
mkdir tax
cp /home/stu_2/bozhang/P3.Assembly_annotation/10.unigene/pwy/unigene_contig.txt ./tax/
cut -f1 /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID/iad_unigene.txt_new_identity_querylength.txt > ./tax/iad_contig_hit.txt

cd tax
for i in `ls /home/stu_2/bozhang/P3.Assembly_annotation/10.unigene/quants/ | grep quant`
do
awk 'NR==FNR{a[$1]; next} $1 in a' /home/stu_2/bozhang/P3.Assembly_annotation/12.IAID/iad_unigene/tax/iad_contig_hit.txt /home/stu_2/bozhang/P3.Assembly_annotation/10.unigene/quants/${i} > ./temp_file.txt
######unigene select the geneid sample dont include
awk '$4 != 0' ./temp_file.txt > ./temp_file2.txt && mv temp_file2.txt ./${i%.quant.sf}quant_result.txt
rm ./temp_file.txt
done

mkdir tax_id

for i in `ls ./ | grep quant`
do
awk -F "\t" 'NR==FNR{a[$1]; next} $2 in a {print $2, $3}' ./${i} ./unigene_contig.txt > ./tax_id/${i%quant_result.txt}_tax.tsv
awk 'NR==FNR{a[$1]=$0;next}{print $0"\t"a[$1]}' ./${i} ./tax_id/${i%quant_result.txt}_tax.tsv> ./tax_id/${i%quant_result.txt}_tax_merged.tsv
awk -F'\t' '{print $1 "\t" $6}' ./tax_id/${i%quant_result.txt}_tax_merged.tsv > ./tax_id/${i%quant_result.txt}_tax_merged_CPM.tsv

#########delete the taxidnumber

sed -i -E 's/ \(taxid [0-9]+\)\t/\t/g' ./tax_id/${i%quant_result.txt}_tax_merged_CPM.tsv

#########reverse the tab

sed -i 's/ /\'$'\t''/'  ./tax_id/${i%quant_result.txt}_tax_merged_CPM.tsv

#######select the data

awk -F'\t' '{
    seen[$2] += $3
}
END {
    for (key in seen) {
        print key "\t" seen[key]
    }
}' ./tax_id/${i%quant_result.txt}_tax_merged_CPM.tsv >./tax_id/${i%quant_result.txt}_tax_merged_CPM_singlebac.tsv
done

cd  tax_id
for i in `ls ./ | grep merged_CPM_singlebac.tsv`
do
 sed "s/^/$(basename "${i}")\t/" "${i%}" > "${i%_tax_merged_CPM_singlebac.tsv}_updated.tsv"
done

cat ./*_updated.tsv > all.tsv