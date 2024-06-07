
kraken_njy=read.delim("njy_kraken",header=F)
kraken_njy_filter=subset(kraken_njy,V3>1)
mat_njy=acast(V2~V1,value="V3",data =kraken_njy_filter ,fill = 0)
mat_njy_rarefy=t(rrarefy(t(mat_njy),min(apply(mat_njy,2 , sum))))


shannon_div=data.frame(shannon=diversity(t(mat_njy_rarefy)))
simpson_div=data.frame(simpson=diversity(t(mat_njy_rarefy),index="simpson"))

ggplot(shannon_div,aes(x=1,y=shannon))+geom_violin()+geom_jitter()+theme_classic()
ggplot(simpson_div,aes(x=1,y=simpson))+geom_violin()+geom_jitter()+theme_classic()

genus_reads_count=data.frame(count=apply(mat_njy_rarefy,1, sum))
genus_reads_count$genus=rownames(genus_reads_count)
genus_reads_count=genus_reads_count[order(genus_reads_count$count,decreasing = T),]
top10_genus=head(genus_reads_count,10)$genus
mat_njy_rarefy_top10=data.frame(t(mat_njy_rarefy[top10_genus,]))
mat_njy_rarefy_top10$other=as.numeric(apply(mat_njy_rarefy,2,sum)-apply(mat_njy_rarefy_top10,1,sum))
mt_top10_mat=melt(as.matrix(mat_njy_rarefy_top10))

shannon_div$ExpID=rownames(shannon_div)
shannon_div=merge(shannon_div,sampleID)


ggplot(shannon_div,aes(x=Library     ,y=shannon,col=Library     ))+geom_violin()+geom_jitter()+theme_classic()

simpson_div$ExpID=rownames(simpson_div)

new_sample=subset(sampleID,Library =="BGI")
top1_count=mat_njy_rarefy[top10_genus[1],]
top1_count=top1_count[order(top1_count)]

p1=ggplot(mt_top10_mat,aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+scale_fill_brewer(palette = "Set3")+xlim(names(top1_count))
p2=ggplot(sampleID,aes(x=ExpID,y=1,fill=Library))+geom_bar(stat="identity")+scale_fill_brewer(palette = "Set3")+xlim(names(top1_count))

top1_count=top1_count[names(top1_count) %in% new_sample$ExpID ]







sample_list=read.delim("all_sample_list")


all_kraken_sp_count=read.delim("all_kraken_sp_count",header=F)
colnames(all_kraken_sp_count)=c("ExpID","all_kraken_sp_count")
all_meta_sp_count=read.delim("all_metaphlan_sp_count")
colnames(all_meta_sp_count)=c("ExpID","all_meta_sp_count")

all_kraken_bac_reads=read.delim("all_kraken_bacterial_count",header = F)
colnames(all_kraken_bac_reads)=c("ExpID","kindom","all_kraken_bac_reads")
all_metaphlan_bac_reads=read.delim("all_metaphlan_bacterial_count")
colnames(all_metaphlan_bac_reads)=c("ExpID","clade_name", "clade_taxid","relative_abundance","coverage", "estimated_number_of_reads_from_the_clade")

all_sample_count=read.delim("all.sample.count",header = F,sep = " ")
colnames(all_sample_count)=c("ExpID","RawReads","FilterReads")

kraken_sp_top1=read.delim("all_kraken_sp_top1",header = F)
 metaphlan_sp_top1=read.delim("all_metaphlan_sp_top1",header = F)
 colnames(kraken_sp_top1)=c("ExpID","Kraken_top1","Kraken_top1_reads")
 colnames(metaphlan_sp_top1)=c("ExpID","metaphlan_top1","clade_taxid","Metaphlan_top1_abundance","Metaphlan_top1_coverage", "Metaphlan_top1_reads")

sample_anno=merge(sample_list,all_sample_count,all.x=T)
sample_anno=merge(sample_anno,all_kraken_bac_reads[,c(1,3)],all.x=T)
sample_anno=merge(sample_anno,all_kraken_sp_count,all.x=T)
sample_anno=merge(sample_anno,kraken_sp_top1,all.x=T)
sample_anno=merge(sample_anno,all_metaphlan_bac_reads[,c(1:6)],all.x=T)
sample_anno=merge(sample_anno,all_meta_sp_count,all.x=T)
sample_anno=merge(sample_anno,metaphlan_sp_top1[,c(1,2,6)],all.x=T)

write.csv(sample_anno,"sample_anno.csv")

kraken_sp_all=read.delim("filtered_kraken_bacterial_species",header = F)
colnames(kraken_sp_all)=c("ExpID","sp","readsCount","mean")
kraken_sp_all$mean=as.numeric(tapply(kraken_sp_all$readsCount,kraken_sp_all$ExpID,mean)[as.character(kraken_sp_all$ExpID)])
kraken_sp_all$sd=as.numeric(tapply(kraken_sp_all$readsCount,kraken_sp_all$ExpID,sd)[as.character(kraken_sp_all$ExpID)])
kraken_sp_all$z_score=(kraken_sp_all$readsCount-kraken_sp_all$mean)/kraken_sp_all$sd
kraken_sp_all$sp_count=as.numeric(table(kraken_sp_all$ExpID)[as.character(kraken_sp_all$ExpID)])
kraken_sp_all_candi=subset(kraken_sp_all,readsCount>50&(z_score>1|sp_count<3))
kraken_sp_all_candi$sp2=kraken_sp_all_candi$sp
kraken_sp_all_candi$sp2=gsub(" ","_",kraken_sp_all_candi$sp2)
kraken_sp_all_candi$sp2=gsub("d__","k__",kraken_sp_all_candi$sp2)
kraken_sp_all_candi$tools="Kraken"
sample_anno$kraken_candi=as.numeric(table(kraken_sp_all_candi$ExpID)[as.character(sample_anno$ExpID)])
sample_anno[is.na(sample_anno$kraken_candi),"kraken_candi"]=0
kraken_sp_all=merge(kraken_sp_all,sample_anno[,c(1,3,5,6)],by="ExpID")

metaphlan_sp_all=read.delim("all_metaphlan_bacterial_species",header = F)
colnames(metaphlan_sp_all)=c("ExpID","sp","tax_id","perc","coverage","readsCount")
metaphlan_sp_all$mean=as.numeric(tapply(metaphlan_sp_all$readsCount,metaphlan_sp_all$ExpID,mean)[as.character(metaphlan_sp_all$ExpID]))
metaphlan_sp_all$sd=as.numeric(tapply(metaphlan_sp_all$readsCount,metaphlan_sp_all$ExpID,sd)[as.character(metaphlan_sp_all$ExpID)])
metaphlan_sp_all$z_score=(metaphlan_sp_all$readsCount-metaphlan_sp_all$mean)/metaphlan_sp_all$sd
metaphlan_sp_all$sp_count=as.numeric(table(metaphlan_sp_all$ExpID)[as.character(metaphlan_sp_all$ExpID)])
metaphlan_sp_all_candi=subset(metaphlan_sp_all,readsCount>50&(z_score>1|sp_count<3))
metaphlan_sp_all_candi$sp2=metaphlan_sp_all_candi$sp
metaphlan_sp_all_candi$sp2=gsub("s__Pseudomonas_aeruginosa_group","s__Pseudomonas_aeruginosa",metaphlan_sp_all_candi$sp2)
metaphlan_sp_all_candi$tools="MetaPhlan"

sample_anno$metaphlan_candi=as.numeric(table(metaphlan_sp_all_candi$ExpID)[as.character(sample_anno$ExpID)])
sample_anno[is.na(sample_anno$metaphlan_candi),"metaphlan_candi"]=0


sp_all_candi=rbind(kraken_sp_all_candi[,c("ExpID","sp2","tools","readsCount")],metaphlan_sp_all_candi[,c("ExpID","sp2","tools","readsCount")])
sp_all_candi$tmp=paste(sp_all_candi$ExpID,sp_all_candi$sp2)
write.csv(acast(sp_all_candi,tmp~tools ,value.var="readsCount",fill=0),"sp_candi.csv")



kraken_sp_all_ascites=subset(kraken_sp_all,sampleType=="ascites"&readsCount >10)
kraken_sp_all_ascites=merge(kraken_sp_all_ascites,ases_group,by.x="ExpID",by.y="SampleID",all.x = "T")
kraken_sp_all_ascites[is.na(kraken_sp_all_ascites$Group),"Group"]="BA"
ascites_sample_sp=table(kraken_sp_all_ascites$sp)
ascites_sample_sp=ascites_sample_sp[order(ascites_sample_sp,decreasing=T)]
write.csv(table(kraken_sp_all_ascites[,c("sp","Group")]),"ascites_group_count.csv")
write.csv(table(subset(kraken_sp_all_ascites,z_score >1)[,c("sp","Group")]),"ascites_group_count.csv")



all_sample_cov_uni=read.delim("all_sample_cov_uni_dep",header = F)
colnames(all_sample_cov_uni)=c("ExpID","candiSp","readCount","genomeSize","coveredGenome","highCover","Coverage","Uniformity","averageDepth")
gplot(all_sample_cov_uni,aes(x=Coverage,y=Uniformity,col=log10(averageDepth)))+geom_point()+theme_classic()+xlab("Coverage")+ylab("Uniformity")+scale_color_gradient2(low = "blue",high="red",mid="yellow",midpoint = 2)+facet_wrap(~sampleType )

colnames(sample_anno)[16]="all_metaphlan_bac_reads"

ggplot(sample_anno,aes(x=all_metaphlan_bac_reads,y=all_metaphlan_bac_reads/RawReads ))+geom_point()+scale_x_log10()+scale_y_log10()+facet_wrap(~sampleType)

ggplot(metaphlan_sp_all,aes(x=readsCount,y=perc/100 ))+geom_point()+scale_x_log10()+theme_classic()+facet_wrap(~sampleType)+scale_y_log10()

ggplot(metaphlan_sp_all,aes(x=readsCount,y=z_score  ))+geom_point()+scale_x_log10()+theme_classic()+facet_wrap(~sampleType)+geom_hline(yintercept = 0.5)

ggplot(sample_anno,aes(x=all_kraken_bac_reads,y=all_kraken_bac_reads/RawReads ))+geom_point()+scale_x_log10()+scale_y_log10()+facet_wrap(~sampleType)+theme_classic()

ggplot( subset(all_sample_cov_uni,IfInclude=="Y"),aes(x=Coverage,y=Uniformity,col=log10(averageDepth)))+geom_point()+theme_classic()+xlab("Coverage")+ylab("Uniformity")+scale_color_gradient2(low = "blue",high="red",mid="yellow",midpoint = 2)+facet_wrap(~sampleType )

all_sample_sp=merge(all_sample_cov_uni,candi_sp,by.x=c("ExpID","candiSp"),by.y=c("ExpID","Species2"))

kraken_ases=read.csv("kraken_sp_ascites.csv")
select_genus=c("g__Klebsiella","g__Pseudomonas","g__Escherichia","g__Enterococcus","g__Staphylococcus")
kraken_ases_genus=subset(kraken_ases,genus %in% select_genus)

kraken_genus_count_mt=melt(tapply(kraken_ases_genus$readsCount,kraken_ases_genus[,c("ExpID","genus")],sum))
kraken_genus_count_mt[is.na(kraken_genus_count_mt$value),"value"]=0
#kraken_genus_count_mt=subset(kraken_genus_count_mt,!is.na(value))
kraken_genus_count_mt$sum=tapply(kraken_ases_genus$readsCount,kraken_ases_genus[,c("ExpID")],sum)[kraken_genus_count_mt$ExpID]
kraken_genus_count_mt$perc=kraken_genus_count_mt$value/kraken_genus_count_mt$sum

ases_culture=read.csv("ases_culture.csv")


kraken_genus_count_mt_culture_pos=subset(kraken_genus_count_mt,ExpID %in% ases_culture$ExpID)
kraken_genus_count_mt_culture_pos=merge(kraken_genus_count_mt_culture_pos,ases_culture[,c(1,3,4)],all.x = T)
kraken_genus_count_mt_culture_pos[is.na(kraken_genus_count_mt_culture_pos$culture),"culture"]="N"

roc_all=roc(kraken_genus_count_mt_culture_pos$culture,kraken_genus_count_mt_culture_pos$perc)

coords_sg=data.frame()
data_sg=data.frame()
for (sg in select_genus){
  roc_sg=roc(subset(kraken_genus_count_mt_culture_pos,genus==sg)$culture,subset(kraken_genus_count_mt_culture_pos,genus==sg)$perc)
  add_data_sg=data.frame(coords(roc_sg))
  add_data_sg$genus=sg
  data_sg=rbind(data_sg,add_data_sg)
  coords_sg=rbind(coords_sg,c(as.numeric(auc(roc_sg)), as.numeric(coords(roc_sg, "best", ret=c("threshold", "specificity", "sensitivity"))),sg))
}
colnames(coords_sg)=c("auc","threshold","specificity", "sensitivity","genus")
coords_sg$sensitivity=as.numeric(coords_sg$sensitivity)
coords_sg$specificity=as.numeric(coords_sg$specificity)
ggplot()+geom_line(data=data_sg,aes(y=sensitivity ,x=specificity ,col=genus))+theme_classic()+geom_point(data=coords_sg,aes(y=sensitivity ,x=specificity ,col=genus))+geom_text(data=coords_sg,aes(y=sensitivity ,x=specificity ,label=paste(genus,threshold ,sep="\n")))


kraken_genus_count_mt_culture_pos=merge(kraken_genus_count_mt_culture_pos,coords_sg[,c("genus","threshold")])
kraken_genus_count_mt_culture_pos$ifPos="N"
kraken_genus_count_mt_culture_pos[kraken_genus_count_mt_culture_pos$perc>kraken_genus_count_mt_culture_pos$threshold,"ifPos"]="Y"

kraken_genus_count_mt=merge(kraken_genus_count_mt,coords_sg[,c("genus","threshold")])
kraken_genus_count_mt$ifPos="N"
kraken_genus_count_mt[kraken_genus_count_mt$perc>kraken_genus_count_mt$threshold,"ifPos"]="Y"

write.csv(subset(kraken_genus_count_mt,ifPos=="Y"),file="aes_genus_filtred.csv")


culture_blood=read.csv("blood_culture_result.csv")
blood_kraken=read.csv("kraken_blood.csv",row.names = "X")
colnames(blood_kraken)=c("ExpID","Patient")

kraken_genus_count_mt=melt(tapply(blood_kraken$readsCount,blood_kraken[,c("Patient","genus")],sum))
kraken_genus_count_mt[is.na(kraken_genus_count_mt$value),"value"]=0
kraken_genus_count_mt$sum=tapply(blood_kraken$readsCount,blood_kraken[,c("Patient")],sum)[kraken_genus_count_mt$Patient]
kraken_genus_count_mt$perc=kraken_genus_count_mt$value/kraken_genus_count_mt$sum

colnames(blood_culture)[4]="culture"
blood_kraken_selected=subset(kraken_genus_count_mt,Patient%in% blood_culture$Patient & genus %in% blood_culture$genus)
blood_kraken_selected=merge(blood_kraken_selected,blood_culture[,c(1,2,4)],all.x="T")
blood_kraken_selected[is.na(blood_kraken_selected$culture),"culture"]="N"

tbl_blood=table(blood_kraken_selected[,c("genus","culture")])
select_genus_bl=rownames(tbl_blood[tbl_blood[,2]>0,])
blood_kraken_selected=subset(blood_kraken_selected,genus %in% select_genus_bl)
tbl_blood=table(blood_kraken_selected[,c("genus","culture")])

blood_kraken_selected$genus=as.character(blood_kraken_selected$genus)
tbl_blood=table(blood_kraken_selected[,c("genus","culture")])

coords_sg=data.frame()
 data_sg=data.frame()
 for (sg in select_genus_bl){
     roc_sg=roc(subset(blood_kraken_selected,genus==sg)$culture,subset(blood_kraken_selected,genus==sg)$perc)
     add_data_sg=data.frame(coords(roc_sg))
     add_data_sg$genus=sg
     data_sg=rbind(data_sg,add_data_sg)
     coords_sg=rbind(coords_sg,c(as.numeric(auc(roc_sg)), as.numeric(coords(roc_sg, "best", ret=c("threshold", "specificity", "sensitivity"))[1,]),sg))
 }

 coords_sg_select=subset(coords_sg,auc>0.5)
 coords_sg_select=coords_sg_select[-2,]
 
blood_kraken_selected=merge(blood_kraken_selected,coords_sg_select[,c("genus","threshold")])
blood_kraken_selected$ngs="N"
blood_kraken_selected[blood_kraken_selected$perc>blood_kraken_selected$threshold,"ngs"]="Y"
table(blood_kraken_selected[,c("culture","ngs")])

ggplot()+geom_line(data=subset(data_sg,genus %in% coords_sg[coords_sg$auc>0.5,"genus"]),aes(y=sensitivity ,x=specificity ,col=genus))+theme_classic()+geom_point(data=coords_sg_select,aes(y=sensitivity ,x=specificity ,col=genus))

kraken_genus_count_mt_sg=merge(kraken_genus_count_mt,coords_sg_select[,c("genus","threshold")])
kraken_genus_count_mt_sg$ngs="N"
kraken_genus_count_mt_sg[kraken_genus_count_mt_sg$perc>kraken_genus_count_mt_sg$threshold,"ngs"]="Y"
kraken_genus_count_mt_sg$genus=as.character(kraken_genus_count_mt_sg$genus)
table(kraken_genus_count_mt_sg[,c("genus","ngs")])

table(subset(blood_culture,Patient%in% kraken_genus_count_mt$Patient & genus %in% kraken_genus_count_mt$genus )$genus)

kraken_sp_new=read.delim("new_blodd.species.count",header = F)
colnames(kraken_sp_new)=c("ExpID","sp","readsCount")

kraken_sp_new$mean=as.numeric(tapply(kraken_sp_new$readsCount,kraken_sp_new$ExpID,mean)[as.character(kraken_sp_new$ExpID)])
kraken_sp_new$sd=as.numeric(tapply(kraken_sp_new$readsCount,kraken_sp_new$ExpID,sd)[as.character(kraken_sp_new$ExpID)])
kraken_sp_new$z_score=(kraken_sp_new$readsCount-kraken_sp_new$mean)/kraken_sp_new$sd

kraken_sp_new_filtered=read.delim("filtered_ecmo_kraken",header = F)
colnames(kraken_sp_new_filtered)=c("ExpID","sp","readsCount","mean")
kraken_sp_new_filtered$mean=as.numeric(tapply(kraken_sp_new_filtered$readsCount,kraken_sp_new_filtered$ExpID,mean)[as.character(kraken_sp_new_filtered$ExpID)])
kraken_sp_new_filtered$sd=as.numeric(tapply(kraken_sp_new_filtered$readsCount,kraken_sp_new_filtered$ExpID,sd)[as.character(kraken_sp_new_filtered$ExpID)])
kraken_sp_new_filtered$z_score=(kraken_sp_new_filtered$readsCount-kraken_sp_new_filtered$mean)/kraken_sp_new_filtered$sd


kraken_sp_all=merge(kraken_sp_all,sample_anno[,c(1,3,5,6)],by="ExpID")

kraken_sp_all_blood=subset(kraken_sp_all,sampleType  =="blood")
kraken_sp_all_blood=rbind(kraken_sp_all_blood[,1:6],kraken_sp_new_filtered[,1:6])
kraken_sp_all_blood=merge(kraken_sp_all_blood,sample_blood,by.x="ExpID",by.y="sampleID")
tbl_blood_sp=table(subset(kraken_sp_all_blood,sample_type!="healthy"&z_score >1)$sp)
tbl_blood_sp=tbl_blood_sp[order(tbl_blood_sp,decreasing = T)]
write.csv(tbl_blood_sp,"tbl_blood_sp.csv")

tbl_blood_sp=table(subset(kraken_sp_all_blood,sample_type!="healthy"&z_score>1&readsCount >100)$sp)

select_sp=c("d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Pseudomonadales|f__Pseudomonadaceae|g__Pseudomonas|s__Pseudomonas aeruginosa","d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Klebsiella|s__Klebsiella pneumoniae","d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Xanthomonadales|f__Xanthomonadaceae|g__Stenotrophomonas|s__Stenotrophomonas maltophilia","d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Pseudomonadales|f__Moraxellaceae|g__Acinetobacter|s__Acinetobacter baumannii","d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Enterobacter|s__Enterobacter cloacae","d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Salmonella|s__Salmonella enterica")

kraken_sp_all_blood$sum=as.numeric(tapply(kraken_sp_all_blood$readsCount,kraken_sp_all_blood$ExpID,sum)[as.character(kraken_sp_all_blood$ExpID)])
kraken_sp_all_blood$sp_perc=kraken_sp_all_blood$readsCount/kraken_sp_all_blood$sum
kraken_sp_all_blood$sp_perc_bac=kraken_sp_all_blood$readsCount/kraken_sp_all_blood$bact_reads  

blood_cutoff=tapply(subset(kraken_sp_all_blood,sample_type=="healthy")$sp_perc_bac,subset(kraken_sp_all_blood,sample_type=="healthy")$sp,mean)
kraken_sp_all_blood$mean_healthy=as.numeric(tapply(subset(kraken_sp_all_blood,sample_type=="healthy")$sp_perc_bac,subset(kraken_sp_all_blood,sample_type=="healthy")$sp,mean)[kraken_sp_all_blood$sp])
kraken_sp_all_blood$sd_healthy=as.numeric(tapply(subset(kraken_sp_all_blood,sample_type=="healthy")$sp_perc_bac,subset(kraken_sp_all_blood,sample_type=="healthy")$sp,sd)[kraken_sp_all_blood$sp])
kraken_sp_all_blood[is.na(kraken_sp_all_blood$mean_healthy),"mean_healthy"]=0
kraken_sp_all_blood[is.na(kraken_sp_all_blood$sd_healthy),"sd_healthy"]=0

kraken_sp_all_blood_1_1_amr=merge(kraken_sp_all_blood_1_1[,1:2],amr_uniq,by.x="ExpID",by.y="V1",all.x=T)
kraken_sp_all_blood_1_1_amr=merge(kraken_sp_all_blood_1_1_amr,amr_anno_uniq,by.x="ExpID",by.y="V1",all.x=T)


kraken_sp_all_blood_identified_amr=merge(kraken_sp_all_blood_identified[,c(1:3,13)],amr_uniq,by.x="ExpID",by.y="V1",all.x=T)
kraken_sp_all_blood_identified_amr=merge(kraken_sp_all_blood_identified_amr,amr_anno_uniq,by.x="ExpID",by.y="V1",all.x=T)

length(unique(subset(kraken_sp_all_blood,sample_type!="healthy"&z_score  >1&sp_perc_bac >mean_healthy&readsCount     >100)$ExpID))
sample_blood$sample_sourc="blood"
ggplot(kraken_sp_all_blood,aes(x=readsCount      )  )+geom_density()+scale_x_log10()+geom_vline(xintercept = 50)+theme_classic()

tbl_blood_sp_100=table(subset(kraken_sp_all_blood,sample_type!="healthy"&z_score  >1&sp_perc_bac >mean_healthy&readsCount     >100)$sp)
tbl_blood_sp_100=tbl_blood_sp_100[order(tbl_blood_sp_100,decreasing = T)]
write.csv(tbl_blood_sp_100,"tbl_blood_sp.csv")

tbl_blood_sp_10=table(subset(kraken_sp_all_blood,sample_type!="healthy"&z_score  >0&sp_perc_bac >mean_healthy)$sp)
tbl_blood_sp_10=tbl_blood_sp_10[order(tbl_blood_sp_10,decreasing = T)]
write.csv(tbl_blood_sp_10,"tbl_blood_sp_10.csv")

tbl_blood_sp_100=table(subset(kraken_sp_all_blood,sample_type!="healthy"&z_score  >1&sp_perc_bac >mean_healthy&readsCount     >50)$sp)
tbl_blood_sp_100=tbl_blood_sp_100[order(tbl_blood_sp_100,decreasing = T)]
write.csv(tbl_blood_sp_100,"tbl_blood_sp_50.csv")

tbl_blood_sp_100=table(subset(kraken_sp_all_blood,sample_type!="healthy"&z_score  >1&sp_perc_bac >mean_healthy+sd_healthy)$sp)
tbl_blood_sp_100=tbl_blood_sp_100[order(tbl_blood_sp_100,decreasing = T)]
write.csv(tbl_blood_sp_100,"tbl_blood_sp_sd.csv")




kraken_sp_all_blood$sp_perc_raw=kraken_sp_all_blood$readsCount/kraken_sp_all_blood$rawReads
kraken_sp_all_blood$mean_healthy_raw=as.numeric(tapply(subset(kraken_sp_all_blood,sample_type=="healthy")$sp_perc_raw,subset(kraken_sp_all_blood,sample_type=="healthy")$sp,mean)[kraken_sp_all_blood$sp])
kraken_sp_all_blood$sdhealthy_raw=as.numeric(tapply(subset(kraken_sp_all_blood,sample_type=="healthy")$sp_perc_raw,subset(kraken_sp_all_blood,sample_type=="healthy")$sp,sd)[kraken_sp_all_blood$sp])
kraken_sp_all_blood_1_1=subset(kraken_sp_all_blood,sample_type!="healthy"&z_score  >1&sp_perc_bac >mean_healthy+sd_healthy)
write.csv(subset(kraken_sp_all_blood,sample_type!="healthy"&z_score  >1&sp_perc_bac >mean_healthy+sd_healthy),"blood_sp_final.csv")

sample_aese=subset(sample_anno,sampleType=="ascites")[,1:7]
sample_aese=merge(sample_aese,ases_group[,1:2],by.x="ExpID",by.y="SampleID",all.x = "T")
kraken_sp_all_ascites$sp_perc=kraken_sp_all_ascites$readsCount/kraken_sp_all_ascites$all_kraken_bac_reads
kraken_sp_all_ascites$mean_healthy=as.numeric(tapply(subset(kraken_sp_all_ascites,Group=="NC")$sp_perc,subset(kraken_sp_all_ascites,Group=="NC")$sp,mean)[kraken_sp_all_ascites$sp])
kraken_sp_all_ascites$sd_healthy=as.numeric(tapply(subset(kraken_sp_all_ascites,Group=="NC")$sp_perc,subset(kraken_sp_all_ascites,Group=="NC")$sp,sd)[kraken_sp_all_ascites$sp])
sample_aese$sp_identified=as.numeric(table(kraken_sp_ascites_filterd$ExpID)[sample_aese$ExpID])
aese_ptc=head(subset(ases_group,Group=="DP"),11)
aese_ptc$Group="PTC"
aese_ntc=subset(ases_group,Group=="NC")
aese_ntc$Group="NTC"
aese_training=rbind(aese_ptc,aese_ntc)

kraken_sp_all_ascites$z_score_healthy=(kraken_sp_all_ascites$sp_perc-kraken_sp_all_ascites$mean_healthy)/kraken_sp_all_ascites$sd_healthy

calculate_backgroud=function(kraken_sp_all_ascites){
  kraken_sp_all_ascites$sp_perc=kraken_sp_all_ascites$readsCount/kraken_sp_all_ascites$all_kraken_bac_reads
  kraken_sp_all_ascites$mean_ntc=as.numeric(tapply(subset(kraken_sp_all_ascites,Group=="NTC")$sp_perc,subset(kraken_sp_all_ascites,Group=="NTC")$sp,mean)[kraken_sp_all_ascites$sp])
  kraken_sp_all_ascites$sd_ntc=as.numeric(tapply(subset(kraken_sp_all_ascites,Group=="NTC")$sp_perc,subset(kraken_sp_all_ascites,Group=="NTC")$sp,sd)[kraken_sp_all_ascites$sp])
  return(kraken_sp_all_ascites)
}



kraken_sp_all_csf=subset(kraken_sp_all,sampleType=="cerebrospinal")
colnames(sample_csf)[2]="Group"
kraken_sp_all_csf=merge(kraken_sp_all_csf,sample_csf[,c("ExpID","all_kraken_bac_reads","Group")],all.x = T)
kraken_sp_all_csf=calculate_backgroud(kraken_sp_all_csf)

kraken_sp_all_csf$z_score_ntc=(kraken_sp_all_csf$sp_perc-kraken_sp_all_csf$mean_ntc)/kraken_sp_all_csf$sd_ntc
kraken_sp_all_csf[is.na(kraken_sp_all_csf$Group),"Group"]=""
kraken_sp_all_csf[kraken_sp_all_csf$Group=="","Group"]="to_test"
kraken_sp_all_csf[is.na(kraken_sp_all_csf$z_score_ntc),"z_score_ntc"]=0

kraken_sp_csf_filtred=subset(kraken_sp_all_csf,z_score>1&sp_perc>mean_ntc+sd_ntc)

write_tbl=function(dataset_select,filename){
  tbl_blood_sp_100=table(dataset_select$sp)
tbl_blood_sp_100=tbl_blood_sp_100[order(tbl_blood_sp_100,decreasing = T)]
write.csv(tbl_blood_sp_100,filename)
}
write_tbl(subset(kraken_sp_all_ascites,Group !="NC"&z_score  >1&sp_perc  >mean_healthy),"tbl_ascite_sp_mean.csv")




length(unique(subset(kraken_sp_all_ascites,Group !="NC"&z_score  >0&sp_perc  >mean_healthy)$ExpID))
write.csv(table(subset(kraken_sp_all_ascites,Group !="NC"&z_score  >0&sp_perc  >mean_healthy)[,c("sp","Group")]),"ascite_sp_group_mean.csv")
write.csv(subset(kraken_sp_all_ascites,Group !="NC"&z_score  >0&sp_perc  >mean_healthy),"ascite_sp_final.csv")
kraken_sp_ascites_filterd=subset(kraken_sp_all_ascites,z_score  >0&sp_perc  >mean_healthy)
kraken_sp_ascites_train_filterd=subset(kraken_sp_ascites_filterd,ExpID %in% aese_training$SampleID)
table(unique(kraken_sp_ascites_train_filterd[,c("ExpID","Group")])$Group)

kraken_sp_all_ascites[is.na(kraken_sp_all_ascites$z_score_healthy),"z_score_healthy"]=0

kraken_sp_ascites_filterd1=subset(kraken_sp_all_ascites,z_score  >0&sp_perc  >mean_healthy+sd_healthy)
kraken_sp_ascites_train_filterd1=subset(kraken_sp_ascites_filterd1,ExpID %in% aese_training$SampleID)
table(unique(kraken_sp_ascites_train_filterd1[,c("ExpID","Group")])$Group)

kraken_sp_all_ascites[kraken_sp_all_ascites$Group=="NC","Group"]="NTC"
kraken_sp_all_ascites[kraken_sp_all_ascites$Group=="DP","Group"]="PTC"
colnames(kraken_sp_all_ascites)[17]="z_score_ntc"

select_best_zh=function(kraken_sp_all_ascites,z){
ascites_train=subset(kraken_sp_all_ascites,z_score>=z&Group%in%c("PTC","NTC"))
to_test_ascite_z=unique(ascites_train$z_score_ntc)
to_test_ascite_z=data.frame(z_h=to_test_ascite_z[order(to_test_ascite_z)])
to_test_ascite_z$PTC=as.numeric(apply(to_test_ascite_z,1,function(a){ length(unique(subset(ascites_train,z_score_ntc>a[1]&Group=="PTC")$ExpID))}))
 to_test_ascite_z$NTC=as.numeric(apply(to_test_ascite_z,1,function(a){ length(unique(subset(ascites_train,z_score_ntc>a[1]&Group=="NTC")$ExpID))}))
 to_test_ascite_z$sense=to_test_ascite_z$PTC/max(to_test_ascite_z$PTC)
to_test_ascite_z$spec=1-to_test_ascite_z$NTC/max(to_test_ascite_z$NTC)
 to_test_ascite_z$auc=to_test_ascite_z$sense+to_test_ascite_z$spec-1
return(to_test_ascite_z)
}

to_test_ascite_z[to_test_ascite_z$auc==max(to_test_ascite_z$auc),]
 
to_test_csf_z=select_best_zh(kraken_sp_all_csf)
to_test_csf_z[to_test_csf_z$auc==max(to_test_csf_z$auc),]
dim(subset(csf_group,ExpID %in% unique(subset(kraken_sp_all_csf,z_score>0&z_score_ntc>3.189433)$ExpID)))
table(subset(csf_group,ExpID %in% unique(subset(kraken_sp_all_csf,z_score>0&z_score_ntc>3)$ExpID))$group)

final_csf_sp=subset(kraken_sp_all_csf,z_score>0&z_score_ntc>3)
final_csf_sp=merge(final_csf_sp,csf_group[,c("ExpID","group","isloate","control")])

kraken_sp_all_blood$z_score_ntc=(kraken_sp_all_blood$sp_perc_bac-kraken_sp_all_blood$mean_healthy)/kraken_sp_all_blood$sd_healthy
kraken_sp_all_blood=merge(kraken_sp_all_blood,blood_group,all.x = "T")
kraken_sp_all_blood[is.na(kraken_sp_all_blood$z_score),"z_score"]=1
to_test_blood_z=select_best_zh(kraken_sp_all_blood)
to_test_blood_z[to_test_blood_z$auc==max(to_test_blood_z$auc),]
ggplot(to_test_blood_z,aes(x=1-sense,y=spec))+geom_point()+geom_line()+theme_classic()

blood_16S_add=read.csv("blood_sample_16S_add.csv")

kraken_sp_all_blood_identified=subset(kraken_sp_all_blood,z_score>=1&z_score_ntc>2.86)
kraken_sp_all_blood_identified=merge(kraken_sp_all_blood_identified,blood_16S_add[,c("ExpID","PatientID","patient_no")],all.x = T)
amr_anno_uniq=merge(amr_anno_uniq,amr_uniq,by=c("V1"))
kraken_sp_all_blood_identified=merge(kraken_sp_all_blood_identified,amr_anno_uniq,all.x = T,by.x="ExpID",by.y="V1")
write.csv(kraken_sp_all_blood_identified,"blood_identified_final.csv")

kraken_mat=acast(kraken_sp_all, sp~ExpID,value.var = "readsCount",fill = 0)
perc_mat=t(t(kraken_mat)/apply(kraken_mat,2,sum))
generate_umap_data=function(mat,samp,clu){
  umap_genus=umap(t(log10(mat+0.000001)))
  umap_genus_data=data.frame(umap_genus$layout)
  umap_genus_data$ExpID=rownames(umap_genus_data)
  umap_genus_data=merge(umap_genus_data,samp)
  if(clu>0){
    dist_umap_genus=hclust(dist(umap_genus$layout),method="ward.D2")
    umap_genus_data$cluster=cutree(dist_umap_genus,clu)[umap_genus_data$ExpID]
  }
  return(umap_genus_data)
}
umap_kraken_all=generate_umap_data(perc_mat,sample_anno[,1:7],0)
ggplot(umap_kraken_all,aes(x=X1,y=X2,col=sampleType ))+geom_point()

generate_pca_data=function(mat,samp,clu){
pca_genus=prcomp(t(mat))
pca_genus_data=data.frame(pca_genus$x[,1:2])
pca_genus_data$ExpID=rownames(pca_genus_data)
pca_genus_data=merge(pca_genus_data,samp)
if(clu>0){
  dist_pca_genus=hclust(dist(pca_genus),method="ward.D2")
  pca_genus_data$cluster=cutree(dist_pca_genus,clu)[pca_genus_data$ExpID]
}
return(pca_genus_data)
}


sample_csf=merge(sample_csf,csf_group[,c("ExpID","control")],all.x=T)

csf_amr_mat=read.csv("csf/csf_amr_count.csv",row.names = "X")
csf_candi_sp_mat=cbind(csf_candi_sp,csf_amr_mat[as.character(csf_candi_sp$ExpID),])

csf_candi_sp=read.csv("csf/csf_candi_sp.csv")
amr_uniq=read.delim("all_sample.resfined_uniq",header = F)
amr_anno_uniq=read.delim("all_sample.res_anno_uniq",header = F)
csf_candi_sp_amr=merge(csf_candi_sp,amr_uniq,by.x="ExpID",by.y="V1")
csf_candi_sp_amr=merge(csf_candi_sp_amr,amr_anno_uniq,by.x="ExpID",by.y="V1")
write.csv(csf_candi_sp_amr,"csf_candi_sp_amr.csv")

blood_culture_all=read.csv("blood_culture_all.csv",row.names = "X")
blood_culture_all=merge(blood_culture_all,patient_ecmo_year,all.x="T")
print_table(subset(blood_culture_all,!is.na(blood_culture_all$year)),"Genus","year")
blood_culture_all$c_bact=paste(blood_culture_all[,2],blood_culture_all[,1],sep="_")
blood_culture_all_ecmo=subset(blood_culture_all,!is.na(blood_culture_all$year))

drug_group=read.csv("drug_group.csv")
rownames(drug_group)=drug_group$drug
drug_group=drug_group[colnames(blood_culture_all)[5:46],]
mat_blood_culture=as.matrix(blood_culture_all[,5:46])
rownames(mat_blood_culture)=paste(blood_culture_all[,2],blood_culture_all[,1],sep="_")
mat_blood_culture[mat_blood_culture!="R"]=0
mat_blood_culture[mat_blood_culture=="R"]=1
mat_blood_culture_drug=t(apply(mat_blood_culture, 1, function(a){tapply(as.numeric(a),drug_group$drug_group,sum)}))
combined_blood_culture_drug=data.frame(mat_blood_culture_drug)
combined_blood_culture_drug$c_bact=rownames(combined_blood_culture_drug)
combined_blood_culture_drug=merge(combined_blood_culture_drug,blood_culture_all[,c(1:3,47,48)])

mt_blood_drug_group=melt(mat_blood_culture_drug)
colnames(mt_blood_drug_group)[1]="c_bact"
mt_blood_drug_group=merge(mt_blood_drug_group,blood_culture_all[,c(1:3,47,48)])
mt_blood_drug_group=subset(mt_blood_drug_group,value>0)
mt_blood_drug_group$if_r=1

mt_blood_drug_group=subset(mt_blood_drug_group,!is.na(mt_blood_drug_group$year))


print_table=function(select,group,cluster){
  to_print_table=data.frame(matrix(0,nrow=nlevels(select[,group])+1,ncol=nlevels(select[,cluster])+1))
  rownames(to_print_table)=c("Total",levels(select[,group]))
  colnames(to_print_table)=c("Total",levels(select[,cluster]))
  tbl_count=table(select[,c(cluster)])
  tbl_perc=table(select[,c(cluster)])/nrow(select)*100
  to_print_table[1,1]=paste0(nrow(select),"(100)")
  for (i in 1:nlevels(select[,cluster])){
    to_print_table[1,(i+1)]=paste0(tbl_count[i],"(",round(tbl_perc[i],0),")")
  }
  
  for(j in 1:nlevels(select[,group])){
    select_s=select[select[,group]==levels(select[,group])[j],]
    tbl_count_s=table(select_s[,c(cluster)])
    to_print_table[(j+1),1]=paste0(nrow(select_s),"(",round(nrow(select_s)/nrow(select)*100,0),")")
    for (i in 1:nlevels(select[,cluster])){
      to_print_table[(j+1),(i+1)]=paste0(tbl_count_s[i],"(",round(tbl_count_s[i]/tbl_count[i]*100,0),")")
    }
  }
  return(to_print_table)
  
}

print_amr_prevelance=function(Genus_s){
  tbl_select=subset(mt_blood_drug_group,Genus==Genus_s)
  count_mat=t(table(tbl_select[,c("Var2","year")]))
  perc_mat= count_mat/as.numeric(table(subset(blood_culture_all_ecmo,Genus==Genus_s)$year)[rownames(count_mat)])
  to_print_table=count_mat
  
  for(j in 1:nrow(to_print_table)){
    rownames(to_print_table)[j]=paste0(rownames(to_print_table)[j],"(",table(subset(blood_culture_all,Genus==Genus_s)$year)[rownames(to_print_table)[j]],")")
    for (i in 1:ncol(to_print_table)){
      to_print_table[j,i]=paste0(count_mat[j,i],"(",round(perc_mat[j,i]*100,1),"%)")
    }
  }
  return(to_print_table)
}
tmp= table(blood_culture_all_ecmo$Genus)
tmp=tmp[order(tmp,decreasing = T)]
write.csv(rbind(print_amr_prevelance(names(tmp)[2]),print_amr_prevelance(names(tmp)[3]),print_amr_prevelance(names(tmp)[4]),print_amr_prevelance(names(tmp)[5]),print_amr_prevelance(names(tmp)[6])),"prevelence.csv")

umap_ecmo_cluster=read.csv("../2022_ECMO/ECMO_2024/umap_genus_mat_all_va.csv")
blood_amr_identified=read.csv("ECMO_identified_amr_final.csv")

blood_16S_add$amr_identified=as.numeric(table(blood_amr_identified$ExpID)[blood_16S_add$ExpID])
blood_16S_add[is.na(blood_16S_add$amr_identified),"amr_identified"]=0
blood_16S_add=merge(blood_16S_add,umap_ecmo_cluster[,c("ExpID3","cluster.x")],by.x="ExpID_16S",by.y="ExpID3",all.x = T)

table(blood_16S_add$cluster.x)
table(subset(blood_16S_add,amr_identified>0)$cluster.x)

ases_final_result=read.csv("腹水筛选得到的结果.csv")
ases_final_result=merge(ases_final_result,amr_anno_uniq,by.x="ExpID",by.y="V1")

