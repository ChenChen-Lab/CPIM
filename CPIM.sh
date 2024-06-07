fq1=$1
fq2=$2
sa=$3

##run metaphlan
/datapool/software/anaconda3/envs/phylophlan/bin/metaphlan --input_type fastq --index mpa_v30_CHOCOPhlAn_201901 --bowtie2out metaphlan/${sa}.res ${fq1},${fq2} -o metaphlan/${sa}.mp --nproc 10 -t rel_ab_w_read_stats
##run kraken
/datapool/software/anaconda3/bin/kraken2  -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  ${fq1} ${fq2}  --threads 10 --output kraken/${sa}.kraken --use-names --report ${sa}.label --use-mpa-style

##select candi species
grep 's__' metaphlan/${sa}.mp |awk -F"|" '{print $NF}' >>${sa}.candi
grep 's__' kraken/${sa}.label |awk -F"|" '{print $NF}' >>${sa}.candi
sort -u ${sa}.candi|while read candi_sp;do grep ${candi_sp}  /datapool/db/complete_genome/refseq_complete.list >>${sa}.candi.ref;
cut -f 3 ${sa}.candi.ref|while read ref;do cat ${ref} >>${sa}.genome.fa

##mapping and calculate coverage and uniformity
bwa index ${sa}.genome.fa
samtools sort ${sa}.bam >${sa}.sorted.bam
samtools mpileup ${sa}.sorted.bam -f ${sa}.genome.fa >${sa}.mpileup
genomeCoverageBed -ibam ${sa}.sorted.bam >${sa}.gcb
awk '$2==0' ${sa}.gcb |awk '{print $0"\t"1-$5"\t"$4-$3}' >${sa}.gcb0
samtools view ${sa}.bam|cut -f 3|sort |uniq -c|awk '{print $2"\t"$1}' >${sa}.rc
my_join.pl -F 1 -f 1 -a ${sa}.rc -b ${sa}.gcb0|cut -f 1-2,4-|awk '{print $0"\t"$2*150/$8*0.2}' >${sa}.rcc
cut -f 1,9 ${sa}.rcc |while read ch de;do a=`awk -F"\t" -v ch=$ch -v de=$de '$1==ch && $4 >de' ${sa}.mpileup|wc -l`;echo $ch $a;done |awk '{print $1"\t"$2}'|my_join.pl -F 1 -f 1 -a ${sa}.rcc -b - >${sa}.rc_unif
my_join.pl -F 3 -f 1 -a ${sa}.sp_geno -b  ${sa}.rc_unif|awk -v sa=${sa} '{ge[$2]+=$8;cov[$2]+=$11;uni[$2]+=$14}END{for (i in ge){print sa"\t"i"\t"ge[i]"\t"cov[i]"\t"uni[i]"\t"cov[i]/ge[i]"\t"uni[i]/ge[i]"\t"uni[i]/cov[i]}}' >${sa}.cov_uni





