library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(UpSetR)
library(reshape2)
library(ggplot2)
library(ggpubr)
cutoff=0.25
stratagy='log2FC_delt_fixed' ## targets' definition: |log2FC| > 0.25 between sh vs. NTC 

## read KD defined targets
df.shs_merged = data.frame()
for(gene in c('PABPN1','PCBP2','U2AF2')){
  df.x = read_rds(paste0('06.PDUI_analys/03.',gene,'.shs_merged.PDUI_melt.rds'))## two sh-replicates were agragated
  df.shs_merged = rbind(df.shs_merged, df.x)
}
saveRDS(df.shs_merged, file='01.all_shs.merged.rds')
###############################################################################
##################### Figure1g ################################################
###############################################################################
df.shs_merged <- readRDS('01.all_shs.merged.rds')
df.shs_merged.signif = subset(df.shs_merged,Is_target=='Yes')%>%
  mutate(value=1,directionality = ifelse(log2FC < 0,'shortened','lengthened'))%>% unique()
saveRDS(df.shs_merged.signif, file='01.all_shs.merged.signif.rds')

## summarize the targets sharing across regulators
df.shs_merged.d = dcast(Gene ~shs_merged , data= df.shs_merged.signif, value.var = 'value', fun.aggregate = mean)
df.shs_merged.d[is.na(df.shs_merged.d)] = 0
pdf(file='03.genes.targets.upsetR.pdf', height=4.5, width=6.5)
upset(df.shs_merged.d, nintersects = 1000)
dev.off()
## summarize the targets number
df.sum = df.shs_merged.signif%>%
  group_by(shs_merged, directionality)%>%
  summarise(N=n())
df.sum$shs_merged = factor(df.sum$shs_merged, levels=c( "U2AF2-1&U2AF2-2","PABPN1-1&PABPN1-2","PCBP2-1&PCBP2-2"  ))
ggplot(data = df.sum, aes(y=shs_merged, x=N,fill=directionality))+
  geom_bar(stat='identity')+
  theme_classic()+
  scale_fill_manual(values=c('#ea8460','#4770a2'))
write.csv(df.sum, file='03.genes.targets.direc.sum.csv')
ggsave(file='03.genes.targets.bar.pdf', height = 2.5, width=5)

###############################################################################
##################### Figure S1 k-l ###########################################
###############################################################################
## include SF3B4
df.rep_merged = data.frame()
for(gene in c('PABPN1','U2AF2','PCBP2')){
  df.x = read_rds(paste0('06.PDUI_analys/03.',gene,'.reps_merged.PDUI_melt.rds')) ## two technical replicates were agragated
  df.rep_merged = rbind(df.rep_merged, df.x)
}

saveRDS(df.rep_merged, file=paste0('02.genes.reps_merged.','_',cutoff,'.PDUI_melt.rds'))

all_combos = combn(unique(df.filter$shID_merged),2)
# sorted_combos <- apply(all_combos, 2, sort)
# unique_sorted_combos <- unique(sorted_combos)

df.jaccard = data.frame()
for(i in seq(ncol(all_combos))){
  shs = unlist(all_combos[,i])
  df.y = subset(df.filter, shID_merged %in% shs)%>%
    filter(Is_target == 'Yes')
  df.sum.x = df.y %>%
    group_by(Gene)%>%
    summarise(N.sh=n())%>%
    group_by(N.sh)%>%
    summarise(N=n())
  df.sum.y = data.frame(sh1=unlist(all_combos[1,i]),
                        sh2= unlist(all_combos[2,i]),
                        intersection=df.sum.x[[which(df.sum.x$N.sh==2),'N']],
                        Jaccard_index = df.sum.x[[which(df.sum.x$N.sh==2),'N']]/sum(df.sum.x$N))
  
  df.jaccard = rbind(df.jaccard, df.sum.y)
  df.sum.y = data.frame(sh2=unlist(all_combos[1,i]),
                        sh1= unlist(all_combos[2,i]),
                        intersection=df.sum.x[[which(df.sum.x$N.sh==2),'N']],
                        Jaccard_index = df.sum.x[[which(df.sum.x$N.sh==2),'N']]/sum(df.sum.x$N))
  df.jaccard = rbind(df.jaccard, df.sum.y)
}
df.jaccard= df.jaccard %>%
  mutate(Is_same_gene = ifelse(str_split_i(sh1,'\\-',1)==str_split_i(sh2,'\\-',1),'Yes','No'))
ggplot(df.jaccard, aes(x=sh2, y=sh1, fill=Jaccard_index))+
  geom_tile()+
  geom_text(aes(x=sh2, y=sh1, label=ifelse(Is_same_gene=='Yes','*','')))+
  scale_fill_distiller(palette = "Spectral")+
  theme_classic()+
  ggtitle('Jaccard index')+
  theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file='04.genes.Jaccard_index.pdf', height = 3.5, width=4.8)
ggplot(df.jaccard, aes(x=sh2, y=sh1, fill=intersection ))+
  geom_tile()+
  geom_text(aes(x=sh2, y=sh1, label=ifelse(Is_same_gene=='Yes','*','')))+
  scale_fill_distiller(palette = "Spectral")+
  theme_classic()+
  ggtitle('Intersection size')+
  theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file='04.genes.Intersection.pdf', height = 3.5, width=4.8)

write.csv(df.jaccard, file='04.genes.jaccard.csv')

###############################################################################
##################### Figure S1 m-n ###########################################
###############################################################################

df.rep_merged = read_rds(paste0('02.genes.reps_merged.','_',cutoff,'.PDUI_melt.rds'))
df.sc = read_rds('0.single_cell_res/250909.single_cell.targets.rds') ## regression analysis based on single cell data
df.sc_eclip_motif = df.sc %>%
  filter(group.single %in% c("cor_eCLIP", "cor_motif"))%>%
  select(Regulator, phenotype_id)%>%
  unique()%>%
  mutate(group.single='cor_motif_eCLIP')
df.sc_only = df.sc %>%
  filter(group.single =='cor_only')
df.sc_eclip_motif =merge(df.sc_only%>% select(-group.single), df.sc_eclip_motif)
df.sc = rbind(df.sc, df.sc_eclip_motif)
df.sc = subset(df.sc, group.single %in% c('cor_only','cor_motif_eCLIP')&method.sc =='spearman')
saveRDS(df.sc, file='03.single_cell.Allcts.spearman.rds')

df.rep_merged = df.rep_merged %>%
  mutate(Regulator=str_split_i(shID_merged,'\\-',1))
df.merge = merge(df.rep_merged, df.sc, by.x=c('Regulator','Gene'), by.y=c('Regulator','phenotype_id'))
saveRDS(df.merge, file='04.sc_ALLcts_KD.spearman.merged_genes.rds')
df.kd_signif = subset(df.merge, Is_target=='Yes')

## integrate AD-DAPA information with PCBP2-KD
df.dapa = read_rds('../re-DAPA_V6/01.spearman_wilcoxon.fdr_0.1.pairwise_regression.rds')%>%
  filter(Is_DAPA=='Yes',celltype %in% c('Exc','Inh')) ## read AD-DAPA results
df.kd_signif = df.kd_signif %>%
  mutate(Is_AD_DAPA = ifelse(Gene %in% unique(df.dapa$phenotype_id),'Yes','No'),
         gene_name = str_split_i(Gene, '\\|', 2))
ggplot(df.kd_signif%>% filter(group.single=='cor_only', Regulator=='PCBP2',celltype=='Exc'), aes(x=log2FC, y=as.numeric(Estimate)))+
  geom_point(size=0.25,aes(col=Is_AD_DAPA))+
  stat_cor(method='pearson')+
  geom_smooth(method='lm')+
  facet_wrap(~shID_merged,scales='free', ncol=2)+
  theme_classic()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_vline(xintercept = 0, linetype='dashed')+
  ggtitle('cor_only')+
  scale_color_manual(values=c('grey','red4'))
geom_text_repel(data=df.kd_signif%>% filter(group.single=='cor_only', Regulator=='PCBP2', abs(delt_PDUI)>0.1)%>%
                  mutate(gene_name= str_split_i(Gene,'\\|',2))%>%
                  group_by(gene_name)%>%
                  arrange(abs(delt_PDUI))%>%slice_head(n=1),
                aes(x=log2FC, y=as.numeric(Estimate), label=gene_name))
ggsave(file='06.PCBP2.DAPA.KD_sc.compare.log2FC.point.cor_only.pdf', height=2.5, width=2*2.6)
ggplot(df.kd_signif%>% filter(group.single=='cor_only', Regulator=='PCBP2',celltype=='Exc'), aes(x=log2FC, y=as.numeric(Estimate),col=Is_AD_DAPA))+
  geom_point(size=0.25,aes(col=Is_AD_DAPA))+
  stat_cor(method='pearson')+
  geom_smooth(method='lm')+
  facet_wrap(~shID_merged,scales='free', ncol=2)+
  theme_classic()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_vline(xintercept = 0, linetype='dashed')+
  ggtitle('cor_only')+
  scale_color_manual(values=c('grey','red4'))
ggsave(file='06.PCBP2.DAPA.KD_sc.compare.log2FC.point.colByDAPA.cor_only.pdf', height=2.5, width=2*2.8)

ggplot(df.kd_signif%>% filter(group.single=='cor_only',celltype=='Exc'), aes(x=log2FC, y=as.numeric(Estimate)))+
  geom_point(size=0.25,col='#54307f')+
  stat_cor(method='pearson')+
  geom_smooth(method='lm')+
  facet_wrap(~shID_merged,scales='free', ncol=2)+
  theme_classic()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_vline(xintercept = 0, linetype='dashed')+
  ggtitle('cor_only')
ggsave(file='05.genes.KD_sc.compare.point.cor_only.pdf', height=2*3, width=2*2)

## prepare the supplementary tables
df.kd_signif= read_rds('01.all_shs.merged.signif.rds')%>%
  mutate(Regulator=str_split_i(shs_merged,'\\-',1))
df.kd_signif1 = df.kd_signif %>% filter(Is_target=='Yes')%>%select(-value,-shs_merged,-Is_target)
write.csv(df.kd_signif1, file='01.all_shs.FinalForStable.csv') 
df.merge = merge(df.kd_signif, df.sc, by.x=c('Regulator','Gene'), by.y=c('Regulator','phenotype_id'))
df.kd_signif = df.merge
## summarize the consistency of PCBP2-targets that is AD-DAPA 
df.pcbp2 = subset(df.kd_signif, Regulator=='PCBP2' & group.single=='cor_only' & celltype=='Exc')%>%
  mutate(Is_AD_DAPA = ifelse(Gene %in% df.dapa$phenotype_id, 'Yes', 'No'),
         Is_cons = ifelse(log2FC*as.numeric(Estimate)<0,'Yes', 'No'))
df.pcbp2.sum = df.pcbp2 %>%
  group_by(Is_cons, Is_AD_DAPA)%>%
  summarise(N=n())

ggplot(df.merge%>% filter(group.single=='cor_only',celltype=='Exc'), aes(x=log2FC, y=as.numeric(Estimate)))+
  geom_point(size=0.15,col='#54307f')+
  stat_cor(method='pearson')+
  geom_smooth(method='lm')+
  facet_wrap(~Regulator,scales='free')+
  theme_classic()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_vline(xintercept = 0, linetype='dashed')+
  ggtitle('cor_only')
ggsave(file='05.genes.KD_sc.compare.shs_merged.point.cor_only.pdf', height=2.2, width=2*3)

df.kd_signif = df.merge
df.kd_signif = df.kd_signif %>%
  mutate(gene_name = str_split_i(Gene, '\\|', 2),
         Is_AD_DAPA = ifelse( Gene %in% unique(df.dapa$phenotype_id),'Yes','No'))
ggplot(df.kd_signif%>% filter(group.single=='cor_only', Regulator=='PCBP2',celltype=='Exc'), aes(x=log2FC, y=as.numeric(Estimate)))+
  geom_point(size=0.25,aes(col=Is_AD_DAPA))+
  stat_cor(method='pearson')+
  geom_smooth(method='lm')+
  facet_wrap(~shs_merged,scales='free', ncol=2)+
  theme_classic()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_vline(xintercept = 0, linetype='dashed')+
  ggtitle('cor_only')+
  scale_color_manual(values=c('grey','red4'))
ggsave(file='06.PCBP2.DAPA.KD_sc.compare.shs_merged.log2FC.point.cor_only.pdf', height=2.5, width=3.6)

ggplot(df.kd_signif%>% filter(group.single=='cor_only', Regulator=='PCBP2',celltype=='Exc'), aes(x=log2FC, y=as.numeric(Estimate),col=Is_AD_DAPA))+
  geom_point(size=0.25,aes(col=Is_AD_DAPA))+
  stat_cor(method='pearson')+
  geom_smooth(method='lm')+
  facet_wrap(~shs_merged,scales='free', ncol=2)+
  theme_classic()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_vline(xintercept = 0, linetype='dashed')+
  ggtitle('cor_only')+
  scale_color_manual(values=c('grey','red4'))
ggsave(file='06.PCBP2.DAPA.KD_sc.compare.shs_merged.log2FC.point.corByDAPA.cor_only.pdf', height=2.5, width=3.6)
## AD_DAPA overalapped targets at gene level 
df.kd_signif = df.kd_signif %>%
  mutate( gene_name = str_split_i(Gene, '\\|', 2),
         Is_AD_DAPA = ifelse(gene_name %in% unique(df.dapa$Gene),'Yes','No'))
targets.list=c('TOLLIP',
               'SSR3',
               'GLYR1', 
               'PNKD',
               'ATL2'
)
df.targets = subset(df.kd_signif, str_split_i(Gene, '\\|',2) %in% targets.list & group.single =='cor_only'& Regulator=='PCBP2')%>%
  mutate(gene_name =str_split_i(Gene, '\\|',2))

ggplot(df.kd_signif%>% filter(group.single=='cor_only', Regulator=='PCBP2',celltype=='Exc'), aes(x=log2FC, y=as.numeric(Estimate),col=Is_AD_DAPA))+
  geom_point(size=0.25,aes(col=Is_AD_DAPA))+
  stat_cor(method='pearson')+
  geom_smooth(method='lm')+
  facet_wrap(~shs_merged,scales='free', ncol=2)+
  theme_classic()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_vline(xintercept = 0, linetype='dashed')+
  ggtitle('cor_only')+
  scale_color_manual(values=c('grey','red4'))+
  geom_point(data=df.kd_signif%>% filter(group.single=='cor_only', Regulator=='PCBP2',Is_AD_DAPA=='Yes',
                                         gene_name %in% targets.list,
                                         log2FC*as.numeric(Estimate)<0)%>%
               mutate(gene_name= str_split_i(Gene,'\\|',2))%>%
               group_by(gene_name)%>%
               arrange(abs(log2FC))%>%slice_head(n=1),
             aes(x=log2FC, y=as.numeric(Estimate)),size=1)+
  geom_text_repel(data=df.kd_signif%>% filter(group.single=='cor_only', Regulator=='PCBP2',Is_AD_DAPA=='Yes',
                                              gene_name %in% targets.list,
                                              log2FC*as.numeric(Estimate)<0)%>%
                    mutate(gene_name= str_split_i(Gene,'\\|',2))%>%
                    group_by(gene_name)%>%
                    arrange(abs(log2FC))%>%slice_head(n=1),
                  aes(x=log2FC, y=as.numeric(Estimate), label=gene_name),
                  box.padding = 0.5,           # 标签周围的填充
                  point.padding = 0.3,         # 点周围的填充
                  arrow = arrow(length = unit(0.01, "npc")),
                  min.segment.length = 0,      # 显示所有连线
                  segment.color = "gray5",
                  segment.size = 0.5)
ggsave(file='06.PCBP2.DAPA.KD_sc.compare.shs_merged.log2FC.point.corByDAPA.GENE.add_anti.cor_only.pdf', height=3.5, width=4.6)


###############################################################################
##################### Figure S1 o #############################################
###############################################################################
## reps_merged
df.reps_merged = read_rds('04.sc_ALLcts_KD.spearman.merged_genes.rds')%>%
  mutate(Is_consistent = ifelse(log2FC*as.numeric(Estimate)<0,'Yes','No'))

df.sum.reps = df.reps_merged %>%
  group_by(shID_merged,Is_consistent,group.single,Is_target)%>%
  summarise(N=n())
df.cons = df.sum.reps %>%
  filter(Is_target=='Yes')%>%
  group_by( shID_merged,group.single  )%>%
  mutate(N.total= sum(N), perct_cons = round(N/N.total*100,2))%>%
  filter(Is_consistent=='Yes')

df.fisher = data.frame()
for(gp in c('cor_motif_eCLIP', 'cor_only')){
  for(sh in unique(df.sum.reps$shID_merged)){
    df.x = subset(df.sum.reps , group.single==gp & shID_merged==sh)
    if(nrow(df.x)>1){
      a=ifelse(is.na(df.x[which(df.x$Is_consistent =='Yes'&df.x$Is_target=='Yes'),'N']%>% as.numeric() ),0,df.x[which(df.x$Is_consistent =='Yes'&df.x$Is_target=='Yes'),'N']%>% as.numeric())## KD & cons
      b=df.x[which(df.x$Is_consistent =='No'&df.x$Is_target=='Yes'),'N']%>% as.numeric() ## KD & not_cons
      c=df.x[which(df.x$Is_consistent =='Yes'&df.x$Is_target=='No'),'N']%>% as.numeric() ## not_KD & cons
      d=df.x[which(df.x$Is_consistent =='No'&df.x$Is_target=='No'),'N']%>% as.numeric() ## not_KD & not_cons
      a = ifelse(is.na(a),0, a)
      b = ifelse(is.na(b),0, b)
      c = ifelse(is.na(c),0, c)
      d = ifelse(is.na(d),0, d)
      df.mat = matrix(c(a,c,b,d), nrow=2)
      df.y = data.frame(df.x[1,c('group.single','shID_merged')],
                        OR=fisher.test(df.mat%>% apply(2,as.numeric))$estimate, pval=fisher.test(df.mat%>% apply(2,as.numeric))$p.value,
                        pva.hyper = dhyper(a,a+c, b+d, a+b))
      if(nrow(df.fisher) > 0){
        df.fisher = rbind(df.fisher, df.y)
      }else{
        df.fisher = df.y
      }
    }
  }
}

df.reps.sum = merge(df.fisher, df.cons, by=c('group.single','shID_merged'))
df.reps.sum = df.reps.sum %>%
  mutate(Regulator = str_split_i(shID_merged, '\\-',1))%>%
  filter(group.single=='cor_only')
ct_group='Exc'
group='shs_merged'
df.shs_sum = data.frame()
df.shs_cons= data.frame()
for(gene in c('PABPN1','PCBP2','U2AF2')){
  df.x = read.csv(paste0('07.compare_sc/',stratagy,'/',gene,'.',group,'.',ct_group,'.compare_sc_',cutoff,'.fisher.csv'))[,-1]
  df.shs_sum = rbind(df.shs_sum, df.x)
  df.y=read.csv(paste0('07.compare_sc/',stratagy,'/',gene,'.',group,'.compare_sc_',cutoff,'.Is_target_consistency.csv'))[,-1]
  df.shs_cons= rbind(df.shs_cons, df.y)
}


df.shs_sum = df.shs_sum%>%
  mutate(Regulator = str_split_i(shs_merged, '\\-',1))%>%
  filter(index.sc =='spearman cor_only')
df.shs_cons = df.shs_cons %>%
  mutate(Regulator = str_split_i(shs_merged, '\\-',1))%>%
  filter(index.sc =='spearman cor_only')
df.shs_cons = read_rds(paste0('05.KD_sc.cons.cutoffs.merged.spearman_only.filter.',stratagy,'.',group_ct,'.rds'))%>%
  filter(cutoff.kd==cutoff)
df.shs_sum = merge(df.shs_sum, df.shs_cons)
names(df.shs_sum)[3] ='shID_merged'

names.list = intersect(names(df.shs_sum), names(df.reps.sum))
df.sum.merge = rbind(df.reps.sum %>% select(names.list), df.shs_sum %>% select(names.list))
df.sum.merge$reps = c('rep1','rep2','rep1','rep2','rep1','rep2','intersect','intersect','intersect')
library(ggrepel)
p1=ggplot(data=df.sum.merge%>% filter(group.single=='cor_only', reps=='intersect'), aes(y=Regulator, x=perct_cons, fill=reps))+
  geom_bar(stat='identity',position=position_dodge(preserve='single'))+
  scale_fill_manual(values=c('#be302d'))+
  theme_classic()+
  geom_text_repel(aes(y=Regulator, x=perct_cons+2,label= paste(round(perct_cons,2),'%', shID_merged)),
                  position =  position_dodge2(width=0.9))+
  geom_vline(xintercept = 50, linetype='dashed')+
  ggtitle('cor_only')+
  theme(legend.position = 'none')

p2=ggplot(data=df.sum.merge%>% filter(group.single=='cor_only', reps=='intersect'), aes(y=Regulator, x=OR, fill=reps))+
  geom_bar(stat='identity',position=position_dodge(preserve='single'))+
  scale_fill_manual(values=c('#be302d'))+
  theme_classic()+
  geom_text_repel(aes(y=Regulator, x=OR, label= case_when(pva.hyper<0.001 ~'***',
                                                          pva.hyper< 0.01 ~ '**',
                                                          pva.hyper < 0.05 ~ '*',
                                                          .default = '')),
                  position =  position_dodge2(width=0.9))+
  geom_vline(xintercept = 1, linetype='dashed')+
  ggtitle('cor_only')+
  theme(legend.position = 'none')

ggsave(plot=p1,file='07.genes.KD_sc.compare.cor_only.cons.pdf', height=3.8, width=5.5)
ggsave(plot=p2,file='07.genes.KD_sc.compare.cor_only.fisher.pdf', height=3.8, width=5)
