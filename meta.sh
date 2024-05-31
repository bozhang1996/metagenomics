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
alias samtools="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif samtools"
alias kraken2="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif kraken2"
alias bracken="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif bracken"
alias kraken-biom="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif kraken-biom"
alias biom="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif biom"
alias megahit="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif megahit"
alias quast.py="singularity  exec /home/stu_2/linguopeng/sif/MetaGenome.sif quast.py "
alias prokka="singularity  exec /home/stu_2/linguopeng/sif/MetaGenome.sif prokka "
alias seqtk="singularity  exec /home/stu_2/linguopeng/sif/MetaGenome.sif seqtk "
alias cd-hit-est="singularity  exec /home/stu_2/linguopeng/sif/MetaGenome.sif cd-hit-est "
alias perl="singularity  exec /home/stu_2/linguopeng/sif/MetaGenome.sif perl "
alias salmon="singularity  exec /home/stu_2/linguopeng/sif/MetaGenome.sif salmon"
alias emapper.py="singularity  exec /home/stu_2/linguopeng/sif/emapper.sif emapper.py"
alias emapperx.R="singularity  exec /home/stu_2/linguopeng/sif/emapper.sif emapperx.R"
alias diamond="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif diamond"
alias salmon="singularity exec /home/stu_2/linguopeng/sif/MetaGenome.sif salmon"
alias Rscript="singularity  exec /home/stu_2/linguopeng/sif/MetaGenome.sif Rscript"
level=S
# # HEADER END
# cd /home/stu_2/BOZHANG/
# mkdir -p P1.data_filter/01.quality
# cd P1.data_filter/01.quality
# mkdir ./qc
# #########fastqc#######
# for i in `ls /home/stu_2/BOZHANG/cleanreads | grep _clean_r1.gz`
# do
# fastqc  --outdir ./qc  --threads 16  /home/stu_2/BOZHANG/cleanreads/${i%_clean_r1.gz*}_clean_r1.gz  /home/stu_2/BOZHANG/cleanreads/${i%_clean_r1.gz*}_clean_r2.gz 
# done
# #########fastp#######
# mkdir ./clean_data
# for i in `ls /home/stu_2/BOZHANG/cleanreads | grep _clean_r1.gz`
# do
# fastp  --thread  16 -i /home/stu_2/BOZHANG/cleanreads/${i%_clean_r1.gz*}_clean_r1.gz -I  /home/stu_2/BOZHANG/cleanreads/${i%_clean_r1.gz*}_clean_r2.gz  -o clean_data/${i%_clean_r1.gz*}_1.fq.gz -O clean_data/${i%_clean_r1.gz*}_2.fq.gz -j  clean_data/${i%_clean_r1.gz*}.fastp.json -h  clean_data/${i%_clean_r1.gz*}.fastp.html
# done
# #########bowtie2#######
# # bowtie2-build  genome.fa genome.db
# # bowtie2-build  c57.fa c57.db 
# # bowtie2-build  hg38.fa hg38.db
# cd ..
# mkdir -p 02.contaminant && cd 02.contaminant
# #########sam#######
# mkdir -p sam
# for i in `ls /home/stu_2/BOZHANG/cleanreads | grep _clean_r1.gz`
# do
# bowtie2 --threads 24 -x /home/stu_2/linguopeng/database/genome/c57.db  -1 ../01.quality/clean_data/${i%_clean_r1.gz*}_1.fq.gz -2 ../01.quality/clean_data/${i%_clean_r1.gz*}_2.fq.gz -S ./sam/${i%_clean_r1.gz*}.sam 2>./sam/${i%_clean_r1.gz*}.map.log 
# done
# #########bam#######
# mkdir -p bam
# for i in `ls /home/stu_2/BOZHANG/cleanreads | grep _clean_r1.gz`
# do
# samtools  view -@ 24 -f 12 ./sam/${i%_clean_r1.gz*}.sam >./bam/${i%_clean_r1.gz*}.unmap.bam 
# done
# rm -rf sam
# #########fastq#######
# mkdir clean_unmap
# for i in `ls /home/stu_2/BOZHANG/cleanreads | grep _clean_r1.gz`
# do
# samtools fastq -1 ./clean_unmap/${i%_clean_r1.gz*}_1.clean.fq.gz -2 ./clean_unmap/${i%_clean_r1.gz*}_2.clean.fq.gz -s ./clean_unmap/${i%_clean_r1.gz*}_singleton.clean.fq.gz ./bam/${i%_clean_r1.gz*}.unmap.bam 
# done
# #######kraken2#######
# cd /home/stu_2/BOZHANG/P1.data_filter/02.contaminant
# cd ../../
# mkdir -p P2.Taxonomic_profiling/1.taxon && cd P2.Taxonomic_profiling/1.taxon
# mkdir kraken
# for i in `ls /home/stu_2/BOZHANG/cleanreads | grep _clean_r1.gz`
# do
# kraken2 --threads 24 --paired --db /home/stu_2/linguopeng/database/k2 --report kraken/${i%_clean_r1.gz*}.kreport --output kraken/${i%_clean_r1.gz*}.kraken /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap/${i%_clean_r1.gz*}_1.clean.fq.gz /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap/${i%_clean_r1.gz*}_2.clean.fq.gz 
# done
# ##########bracken#######
# cd /home/stu_2/BOZHANG/P2.Taxonomic_profiling/1.taxon
# # mkdir out_$level
# # for i in `ls /home/stu_2/BOZHANG/cleanreads | grep _clean_r1.gz`
# # do
# # bracken -d /home/stu_2/linguopeng/database/k2 -i  kraken/${i%_clean_r1.gz*}.kreport -o out_$level/${i%_clean_r1.gz*}.bracken.$level -w out_$level/${i%_clean_r1.gz*}.kreport  -l  $level -t 24
# # done
# kraken-biom  kraken/*.kreport --max D  -o  ./out_$level/$level.biom  
# biom  convert -i  ./out_$level/$level.biom -o  ./out_$level/$level.count.tsv.tmp  --to-tsv --header-key taxonomy
# sed 's/; g__\([^;]\+\); s__/; g__\1; s__\1 /'   ./out_$level/$level.count.tsv.tmp >  ./out_$level/$level.taxID.count.tsv
# sed  '/^#/! s/^[0-9]\+\t\(.*[A-Za-z]\+__\([^;]\+\)\)$/\2\t\1/'  ./out_$level/$level.taxID.count.tsv  >  ./out_$level/$level.taxName.count.tsv
# sed -e '1d' -e '2s/^#//' ./out_$level/$level.taxName.count.tsv | awk -F "\t" -v OFS="\t" '{NF--; print}' | sed 's/\t$//' > ./out_$level/$level.count.tsv
# Rscript  /home/stu_2/linguopeng/script/draw_taxonBarplot.R    out_$level/$level.count.tsv  10   out_$level/$level.count.out
# mkdir -p /home/stu_2/BOZHANG/P3.Assembly_annotation/01.megahit && cd /home/stu_2/BOZHANG/P3.Assembly_annotation/01.megahit
# for i in `ls /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap | grep _1.clean.fq.gz`
# do
# megahit \
   # -1 /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap/${i%_1.clean.fq.gz*}_1.clean.fq.gz \
   # -2 /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap/${i%_1.clean.fq.gz*}_2.clean.fq.gz \
   # --min-contig-len 500 \
   # --tmp-dir  ./ \
   # --memory  0.8 \
   # --num-cpu-threads 24 \
   # --out-dir ${i%_clean_r1.gz*}_megahit \
   # --out-prefix  ${i%_clean_r1.gz*}	
# done
# mkdir /home/stu_2/BOZHANG/P3.Assembly_annotation/02.quast  && cd /home/stu_2/BOZHANG/P3.Assembly_annotation/02.quast
# for i in `ls /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap | grep _1.clean.fq.gz`
# do
# ln -s /home/stu_2/BOZHANG/P3.Assembly_annotation/01.megahit/${i%_1.clean.fq.gz*}_1.clean.fq.gz_megahit/${i%_1.clean.fq.gz*}_1.clean.fq.gz.contigs.fa  ./${i%_1.clean.fq.gz*}.megahit.fa
# quast.py ./${i%_1.clean.fq.gz*}.megahit.fa
# done
# mkdir -p /home/stu_2/BOZHANG/P3.Assembly_annotation/03.prokka  
# cd /home/stu_2/BOZHANG/P3.Assembly_annotation/03.prokka 
# for i in `ls /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap | grep _1.clean.fq.gz`
# do
# ln -s /home/stu_2/BOZHANG/P3.Assembly_annotation/01.megahit/${i%_1.clean.fq.gz*}_1.clean.fq.gz_megahit/${i%_1.clean.fq.gz*}_1.clean.fq.gz.contigs.fa ./${i%_1.clean.fq.gz*}.contigs.fa 
# prokka --outdir ${i%_1.clean.fq.gz*}_prokka --prefix ${i%_1.clean.fq.gz*} --addgenes --addmrna --locustag ${i%_1.clean.fq.gz*} --kingdom Bacteria --metagenome --cpus 24 ./${i%_1.clean.fq.gz*}.contigs.fa 
# done
# mkdir -p /home/stu_2/BOZHANG/P3.Assembly_annotation/04.cdhit 
# cd /home/stu_2/BOZHANG/P3.Assembly_annotation/04.cdhit 
# cat  ../03.prokka/*_prokka/*.ffn > all.trans.fa
# cat  ../03.prokka/*_prokka/*.faa > all.pep.fa
# sed -n "s/^>\(\S\+\).*$/\1/p"  all.pep.fa > all.cds.id
# seqtk subseq  all.trans.fa all.cds.id > all.cds.fa
# cd-hit-est -i all.cds.fa  -o  all.cds.cdhit -c 0.95 -aS 0.9  -M 204800 -T 30 
# cp  all.cds.cdhit  unigene_cds.fasta
# awk '$1 ~/^>/{print $1}'  all.cds.cdhit | sed 's/^>//' > unigene.id
# seqtk  subseq  all.pep.fa unigene.id > unigene_pep.fasta
# awk '{if($1~/^>/){printf $1 $2} else if($NF ~ /*$/){print "\t"$3}}' all.cds.cdhit.clstr |sed 's/>//g; s/...$//'|awk '{print $2"\t"$1}'  > map_id.txt
# perl /home/stu_2/linguopeng/script/map_data_ids  map_id.txt  unigene_cds.fasta 
# perl /home/stu_2/linguopeng/script/map_data_ids  map_id.txt  unigene_pep.fasta
# cd  /home/stu_2/BOZHANG/P3.Assembly_annotation && mkdir 05.abundance && cd  05.abundance
# ln  -s   ../04.cdhit/unigene_cds.fasta
# cd  /home/stu_2/BOZHANG/P3.Assembly_annotation/05.abundance
# salmon index -t  unigene_cds.fasta  -i  unigene_index -p 30
#for i in `ls /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap | grep _1.clean.fq.gz`
#do
#salmon  quant --validateMappings  --meta -p 30   -i  unigene_index  -l IU   -1 /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap/${i%_1.clean.fq.gz*}_1.clean.fq.gz -2  /home/stu_2/BOZHANG/P1.data_filter/02.contaminant/clean_unmap/${i%_1.clean.fq.gz*}_2.clean.fq.gz  -o quants/${i%_1.clean.fq.gz*}.quant
#done
# mkdir /home/stu_2/BOZHANG/P3.Assembly_annotation/06.eggnog
# cd /home/stu_2/BOZHANG/P3.Assembly_annotation/06.eggnog
# emapper.py -i unigene_pep.fasta -o my --itype proteins -m diamond --cpu 30 --override 1>emapper.log 2>&1 
# emapperx.R my.emapper.annotations unigene_pep.fasta
# mkdir /home/stu_2/BOZHANG/P3.Assembly_annotation/07.uniprot
cd /home/stu_2/BOZHANG/P3.Assembly_annotation/07.uniprot
# ln -s ../04.cdhit/unigene_pep.fasta
# diamond blastp  \
 #   --db /home/stu_2/linguopeng/database/uniprot/uniref90 \
  #  --query  unigene_pep.fasta \
   # --out unigene_pep.uniref90.m6 \
   # --threads 30 \
   # --outfmt 6 \
   # --max-target-seqs 1 \
   # --evalue 1e-5
## 基于idmapping 提取GO注释信息
 cd /home/stu_2/BOZHANG/P3.Assembly_annotation/07.uniprot
# perl /home/stu_2/linguopeng/script/uniref90_idmapping.pl \
  #  unigene_pep.uniref90.m6 \
 #   /home/stu_2/linguopeng/database/uniprot/idmapping_selected.tab \
 #   > unigene_pep.uniref90.GOanno
Rscript  /home/stu_2/linguopeng/script/GOmapperx.R  unigene_pep.uniref90.GOanno 
# cd /home/stu_2/BOZHANG/P3.Assembly_annotation/05.abundance
# quants=` ls quants/ | awk '{print "quants/"$1 }' | tr '\n'  ' ' `
# names=` ls quants/ | sed 's/.quant$//'| tr '\n'   ' '  `
# salmon quantmerge --quants $quants --names  $names  --column tpm -o unigenens.tpm
# salmon quantmerge --quants $quants --names  $names  --column numreads -o unigenens.count
# sed '1s/^Name\t//' unigenens.count > unigenens.count.matrix

#ln -s  /home/stu_2/BOZHANG/P3.Assembly_annotation/06.eggnog/R_Library/ .
#ln -s ../06.eggnog/my.emapper.annotations
Rscript \
        ../../../linguopeng/script/enrich_analysis.R \
        ./DE_out/unigenens.count.matrix.FH_vs_SC.DESeq2.DE_results \
        ./unigene.annotations \
        ./DE_out/unigenens.count.matrix.FH_vs_SC.DESeq2.DE_results.enrich