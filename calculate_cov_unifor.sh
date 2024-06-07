sa=$1

samtools sort ${sa}.bam >${sa}.sorted.bam
samtools mpileup ${sa}.sorted.bam -f ${sa}.genome.fa >${sa}.mpileup
genomeCoverageBed -ibam ${sa}.sorted.bam >${sa}.gcb
awk '$2==0' ${sa}.gcb |awk '{print $0"\t"1-$5"\t"$4-$3}' >${sa}.gcb0
samtools view ${sa}.bam|cut -f 3|sort |uniq -c|awk '{print $2"\t"$1}' >${sa}.rc
my_join.pl -F 1 -f 1 -a ${sa}.rc -b ${sa}.gcb0|cut -f 1-2,4-|awk '{print $0"\t"$2*150/$8*0.2}' >${sa}.rcc
cut -f 1,9 ${sa}.rcc |while read ch de;do a=`awk -F"\t" -v ch=$ch -v de=$de '$1==ch && $4 >de' ${sa}.mpileup|wc -l`;echo $ch $a;done |awk '{print $1"\t"$2}'|my_join.pl -F 1 -f 1 -a ${sa}.rcc -b - >${sa}.rc_unif 
my_join.pl -F 3 -f 1 -a ${sa}.sp_geno -b  ${sa}.rc_unif|awk -v sa=${sa} '{ge[$2]+=$8;cov[$2]+=$11;uni[$2]+=$14}END{for (i in ge){print sa"\t"i"\t"ge[i]"\t"cov[i]"\t"uni[i]"\t"cov[i]/ge[i]"\t"uni[i]/ge[i]"\t"uni[i]/cov[i]}}' >${sa}.cov_uni
