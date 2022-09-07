#example for HPV6
#we present one example for E1 gene

library(stringr)
library(dplyr)
setwd("F:/HPV/L324_trim/pileup_apobec")

#we get all pileup files within the directory
filenames = dir(pattern="*.pileup")

#we import the annotation file
gff_file = read.table(file = "coordinates3.txt", header=T)

#we define each category with the different genbank identifiers
#HIGH RISk
high_risk=c("K02718.1", "X05015.1", "J04353.1", "M12732.1",
            "X74477.1","M62849.1", "X74479.1", "M62877.1",
            "X74481.1","X74483.1","D90400.1", "X77858.1",
            "U31794.1", "DQ080079.1")

#LOW RISk
low_risk=c("X00203.1", "M14119.1", "M73236.1", "AJ620205.1",
           "U31788.1", "U37488.1", "U31793.1", "X94164.1",
           "AJ620209.1", "AF436128.1", "AY057438.1")

#INTERMEDIATE RISK
intermediate_risk=c("X74482.1", "D21208.1", "U21941.1", "X94165.1",
                    "AB027021.1")

#UNKNOWN RISK
unknown_risk=c("X74475.1", "X74476.1", "AY395706.1", "AB027020.1",
               "AF436130.1", "AF151983.1", "AF293960.1", "AF349909.1",
               "AJ400628.2", "DQ080081.1", "GQ244463.1")

#import a file that indicates the HPV types for each genbank accession number
ref_types = read.table("HPV_ref_type.txt", header = T, sep = "\t")

#summary_table initialization
summary_apobec=data.frame(sample=character(),genome=character(),mut_apo=integer(), nb_tcw=integer(), mut_non_apo=integer(), nb_ct=integer(),
                          ratio_apo=double(), ratio_non_apo=double(), fold_change=double(),pval=double())

#for each sample
for(my_file in filenames){
  print(my_file)
  #we open the file as a tsv format
  file_name = sub("\\..*", "", my_file)
  my_data = readr::read_tsv(my_file, col_names = F)
  my_data = my_data[,c(1:5)]
  
  #we import a file that indicates the different HPV observed for each sample 
  my_references = read.table("reference_observed.txt", header = T, sep = "\t")
  my_references = my_references[my_references$sample==file_name,]
  my_references = c(my_references$references)
  
  final_df = data.frame(1)
  
  #for each type in the sample
  for(one_ref in my_references){

    #we focus on a specific reference, and the E1 gene
    gff_sub = gff_file[gff_file$V1==one_ref,]
    gff_sub = gff_sub[gff_sub$V10=="E1",]
    
    print(one_ref)
    my_data_focus = my_data[my_data$X1==one_ref,]
    ref_genome = my_data_focus$X1[1]
    
    #we set all variables to 0
    lacking_data=0
    my_data_focus$prop_lack=0
    my_data_focus$match=0
    my_data_focus$mut_A=0
    my_data_focus$mut_T=0
    my_data_focus$mut_G=0
    my_data_focus$mut_C=0
    my_data_focus$insert=0
    my_data_focus$del=0
    
    #we count the number and type of observation at each position of the genome
    for(i in 1:nrow(my_data_focus)){
      #we consider a lacking position when the deapth of coverage is less than 10X
      if(my_data_focus$X4[i]<10){
        lacking_data=lacking_data+1
      } else{
        my_data_focus$prop_lack[i] = 0
        my_data_focus$match[i] = str_count(my_data_focus$X5[i], "\\.")
        my_data_focus$match[i] = my_data_focus$match[i] + str_count(my_data_focus$X5[i], "\\,")
        my_data_focus$mut_A[i] = str_count(my_data_focus$X5[i], "A|a")
        my_data_focus$mut_T[i] = str_count(my_data_focus$X5[i], "T|t")
        my_data_focus$mut_G[i] = str_count(my_data_focus$X5[i], "G|g")
        my_data_focus$mut_C[i] = str_count(my_data_focus$X5[i], "C|c")
        my_data_focus$insert[i] = str_count(my_data_focus$X5[i], "\\+")
        my_data_focus$del[i] = str_count(my_data_focus$X5[i], "\\*")
      }
      
    }
    
    #we calculate the percentage of each observation
    my_data_focus$prop_match= my_data_focus$match/my_data_focus$X4*100
    my_data_focus$prop_mut_A = my_data_focus$mut_A/my_data_focus$X4*100
    my_data_focus$prop_mut_T = my_data_focus$mut_T/my_data_focus$X4*100
    my_data_focus$prop_mut_G = my_data_focus$mut_G/my_data_focus$X4*100
    my_data_focus$prop_mut_C = my_data_focus$mut_C/my_data_focus$X4*100
    my_data_focus$prop_insert = my_data_focus$insert/my_data_focus$X4*100
    my_data_focus$prop_del = my_data_focus$del/my_data_focus$X4*100
    
    #for each position, we define the consensus base
    for(i in 1:nrow(my_data_focus)){
      #if the depth is less than 10X, it is a lacking position
      if(my_data_focus$X4[i]<10){
        my_data_focus$prop_lack[i] = 100
        my_data_focus$prop_match[i] = 0
        my_data_focus$prop_mut_A[i] = 0
        my_data_focus$prop_mut_T[i] = 0
        my_data_focus$prop_mut_G[i] = 0
        my_data_focus$prop_mut_C[i] = 0
        my_data_focus$prop_insert[i] = 0
        my_data_focus$prop_del[i] = 0
      }
      
      my_data_focus$obs[i] = sum(my_data_focus[i,c(14:20)]>50)
      
      #if the proportion of insertion is > 50%, we add the additional bases
      if(my_data_focus$prop_insert[i]>50){
        
        test = my_data_focus$X5[i]
        test2 = str_split(test, "\\+")
        test3=as.data.frame(test2, col.names = "observ")
        for(j in 1:nrow(test3)){
          test3$observ[j]=gsub("\\,","",test3$observ[j])
          test3$observ[j]=gsub("\\.","",test3$observ[j])
          test3$observ[j]=gsub("\\$","",test3$observ[j])
          test3$observ[j]=toupper(test3$observ[j])
        }
        test3$number=1
        #we check the different insertion detected
        test4 = test3 %>%
          group_by(observ) %>%
          summarise(total=sum(number)) %>%
          top_n(1, total)
        #we consider the insertion the most present
        insert_seq = test4$observ[1]
        insert_seq = gsub('[[:digit:]]+', '', insert_seq)
        my_data_focus$prop_match[i] = 0
        my_data_focus$prop_mut_A[i] = 0
        my_data_focus$prop_mut_T[i] = 0
        my_data_focus$prop_mut_G[i] = 0
        my_data_focus$prop_mut_C[i] = 0
        my_data_focus$prop_insert[i] = test4$total[1]/my_data_focus$X4[i]*100
        my_data_focus$prop_del[i] = 0
        #my_data_focus$prop_remo_del[i] = 0
      }
      
      #if the proportion of bases is >50% of match, we consider a match
      #else, it is a mutation
      if(my_data_focus$prop_match[i]>50){
        my_data_focus$sequence[i]=my_data_focus$X3[i]
      }else if(my_data_focus$prop_mut_A[i]>50){
        my_data_focus$sequence[i]="A"
      }else if(my_data_focus$prop_mut_T[i]>50){
        my_data_focus$sequence[i]="T"
      }else if(my_data_focus$prop_mut_G[i]>50){
        my_data_focus$sequence[i]="G"
      }else if(my_data_focus$prop_mut_C[i]>50){
        my_data_focus$sequence[i]="C"
      }else if(my_data_focus$prop_insert[i]>50){
        my_data_focus$sequence[i]=paste(my_data_focus$X3[i], insert_seq, sep="")
      }else if(my_data_focus$prop_del[i]>50){
        my_data_focus$sequence[i]="-"
      }else{
        my_data_focus$sequence[i]="N"
      }
    }
    
    #we generate a sub_data table
    my_sub_data = my_data_focus[,c(1,2,3,22)]
    #we initialize if mutations is possible apobec or not
    my_sub_data$apobec="No"
    my_sub_data$apobec_possible="No"
    
    possible_apobec=0
    possible_non_apobec=0
    real_apobec=0
    non_apobec_mut=0
    
    #for each position
    for(my_position in 3:nrow(my_sub_data)){
      if(my_position>=gff_sub$V4 && my_position<=gff_sub$V5){
        #we define the codon from the reference and the sample
        pos1_ref = my_sub_data$X3[my_position-2]
        pos2_ref = my_sub_data$X3[my_position-1]
        pos3_ref = my_sub_data$X3[my_position]
        pos1_seq = my_sub_data$sequence[my_position-2]
        pos2_seq = my_sub_data$sequence[my_position-1]
        pos3_seq = my_sub_data$sequence[my_position]
        motif_ref = paste(pos1_ref, pos2_ref, pos3_ref, "")
        motif_ref = gsub(" ", "", motif_ref, fixed = T)
        motif_seq = paste(pos1_seq, pos2_seq, pos3_seq, "")
        motif_seq = gsub(" ", "", motif_seq, fixed = T)
        #if the codon from the reference is TCW
        if(motif_ref=="TCA" || motif_ref=="TCT"){
          possible_apobec=possible_apobec+1
          #if the sequence is not well covered, we cannot exclude that an apobec mutation exists
          if(pos1_seq=="N" || pos2_seq=="N" || pos3_seq=="N"){
            possible_apobec=possible_apobec-1
          }else{
            #if the sequence has a c>t mutation in TCW, it is apobec mutation
            if(motif_seq=="TTA" || motif_seq=="TTT"){
              real_apobec = real_apobec+1
            }
          }
        #if the codon from the reference is not TCW
        }else{
          #if we have a c nucleotide not in TCW
          if(pos2_ref=="C" && motif_ref!="TCT" && motif_ref!="TCA"){
            possible_non_apobec=possible_non_apobec+1
            #if the sequence is not well covered, we cannot exclude that the nucleotide is a c
            if(pos1_seq=="N" || pos2_seq=="N" || pos3_seq=="N"){
              possible_non_apobec=possible_non_apobec-1
            }
            #if there if a mutation, add c>t mutation not in TCW motif
            if(pos2_seq=="T"){
              non_apobec_mut = non_apobec_mut + 1
            }
          }
        }
      }
    }
    
    #we perform a fisher test to check the enrichment of apobec mutations
    my_fisher = fisher.test(matrix(c((real_apobec),((possible_apobec)-(real_apobec)), non_apobec_mut, ((possible_non_apobec)-(non_apobec_mut))),2,2, byrow=T))
    my_p_val=my_fisher$p.value
    #we calcutate the ratios of c>t mutations compared to the number of TCW and no TCW motifs, respectively
    ratio_apo = (real_apobec/possible_apobec*100)
    ratio_non_apo=(non_apobec_mut/possible_non_apobec*100)
    fold_change = ratio_apo/ratio_non_apo
    #we add in the summarized table the different statistics
    summary_apobec[nrow(summary_apobec)+1,] = list(file_name, one_ref,real_apobec, possible_apobec,
                                                   non_apobec_mut,possible_non_apobec,
                                                   ratio_apo, ratio_non_apo, fold_change, my_p_val)
  }
}

#import a file that indicates the HPV types for each genbank accession number
HPV_ref = read.table(file = "HPV_ref_type.txt", header=T, sep="\t")

#we add the risk category and the corresponding HPV type in the dataframe
for(i in 1:nrow(summary_apobec)){
  for(j in 1:nrow(HPV_ref)){
    if(summary_apobec$genome[i]==HPV_ref$genome[j]){
      summary_apobec$type_HPV[i]=HPV_ref$HPV_type[j]
    }
  }
  if(summary_apobec$genome[i] %in% high_risk){
    summary_apobec$risk[i]="high"
  }else if(summary_apobec$genome[i] %in% low_risk){
    summary_apobec$risk[i]="low"
  }else if(summary_apobec$genome[i] %in% intermediate_risk){
    summary_apobec$risk[i]="intermediate"
  }else{
    summary_apobec$risk[i]="unknown"
  }
}

#write the table into a file, here for E1
write.table(summary_apobec, "new_details_apobec_E1.txt", sep="\t", row.names=F, col.names = T)

#we consider a ratioC>T >=2 as a gene sequence enriched in APOBEC mutations
summary_apobec_significant = summary_apobec[summary_apobec$fold_change>=2,]
summary_apobec_significant2 = summary_apobec_significant[complete.cases(summary_apobec_significant),]

####
####
#When the analysis was performed for all the genes
####
####

#we import the results for all the genes
E1 = read.table("new_details_apobec_E1.txt", sep="\t", header=T)
E1$gene="E1"
E2 = read.table("new_details_apobec_E2.txt", sep="\t", header=T)
E2$gene="E2"
E4 = read.table("new_details_apobec_E4.txt", sep="\t", header=T)
E4$gene="E4"
E6 = read.table("new_details_apobec_E6.txt", sep="\t", header=T)
E6$gene="E6"
E7 = read.table("new_details_apobec_E7.txt", sep="\t", header=T)
E7$gene="E7"
L1 = read.table("new_details_apobec_L1.txt", sep="\t", header=T)
L1$gene="L1"
L2 = read.table("new_details_apobec_L2.txt", sep="\t", header=T)
L2$gene="L2"
#we combine the results of the different genes
all_genes = rbind(E1,E2,E4,E6,E7,L1,L2)

#we focus on a specific risk, here low risk
all_genes_low = all_genes[all_genes$risk=="low",]
#we consider only significant APOBEC mutations, i.e. with ratioC>T >=2
all_genes_low = all_genes_low[all_genes_low$fold_change>=2,]
all_genes_low = all_genes_low[complete.cases(all_genes_low),]
#we calculate the number of observations
all_genes_low2 = all_genes_low %>%
  group_by(gene)%>%
  summarize(total=n())
#for E4 and E7, we had 0 sample with APOBEC mutations we thus update this information
all_genes_low2[nrow(all_genes_low2)+1,] = list("E4",0)
all_genes_low2[nrow(all_genes_low2)+1,] = list("E7",0)
all_genes_low2$group="low"
#we calculate the proportion of samples having APOBEC mutations for each gene, among all samples with APOBEC mutations
all_genes_low2$repartition=all_genes_low2$total/sum(all_genes_low2$total)*100

#we perform the same analyses for other risk categories
all_genes_high = all_genes[all_genes$risk=="high",]
all_genes_high = all_genes_high[all_genes_high$fold_change>=2,]
all_genes_high = all_genes_high[complete.cases(all_genes_high),]
all_genes_high2 = all_genes_high %>%
  group_by(gene)%>%
  summarize(total=n())
all_genes_high2[nrow(all_genes_high2)+1,] = list("E2",1)
all_genes_high2$group="high"
all_genes_high2$repartition=all_genes_high2$total/sum(all_genes_high2$total)*100

all_genes_inter = all_genes[all_genes$risk=="intermediate",]
all_genes_inter = all_genes_inter[all_genes_inter$fold_change>=2,]
all_genes_inter = all_genes_inter[complete.cases(all_genes_inter),]
all_genes_inter2 = all_genes_inter %>%
  group_by(gene)%>%
  summarize(total=n())
all_genes_inter2$group="intermediate"
all_genes_inter2$repartition=all_genes_inter2$total/sum(all_genes_inter2$total)*100

all_genes_unknown = all_genes[all_genes$risk=="unknown",]
all_genes_unknown = all_genes_unknown[all_genes_unknown$fold_change>=2,]
all_genes_unknown = all_genes_unknown[complete.cases(all_genes_unknown),]
all_genes_unknown2 = all_genes_unknown %>%
  group_by(gene)%>%
  summarize(total=n())
all_genes_unknown2[nrow(all_genes_unknown2)+1,] = list("E4",1)
all_genes_unknown2[nrow(all_genes_unknown2)+1,] = list("E7",1)
all_genes_unknown2$group="unknown"
all_genes_unknown2$repartition=all_genes_unknown2$total/sum(all_genes_unknown2$total)*100

#we combine the proportions for all categories
combine=rbind(all_genes_low2, all_genes_high2, all_genes_inter2, all_genes_unknown2)

library(ggplot2)
#plot the proportion
p = ggplot(combine)+
  geom_bar(aes(x = gene, y = total, fill=group),
           stat="identity", position=position_dodge())

#show the plot
p
