#example for HPV11
setwd("D:/HPV11")
#we import the alignment of HPV11 sequences
#note that each sequence is on one line
my_sequences <- readLines("hpv11_aligned_clean.fasta")
type="HPV11"
risk="low"

#initialization of dataframes
apo_mut = data.frame(id_seq=character(), type=character(),risk=character(),position=integer())
apo_mut_noCT = data.frame(id_seq=character(), type=character(),risk=character(),position=integer(),mut=character())
tcw_motif = data.frame(id_seq=character(), type=character(),risk=character(),position=integer())
total_apo = data.frame(id_seq=character(), type=character(),risk=character(),c_to_t=integer())

####First line = reference sequence
for(my_seq in 1:length(my_sequences)){
  nb_apo = 0
  #even numbers correspond to sequences
  if(my_seq%%2 == 0){
    #we read the reference genome
    if(my_seq == 2){
      ref_split <- strsplit(my_sequences, "")[[2]]
      for(i in 3:length(ref_split)){
        #we check if we have a TCW motif
        if(ref_split[i-2] == "t"){
          if(ref_split[i-1] == "c"){
            if(ref_split[i] == "a" || ref_split[i] == "t"){
              tcw_motif[nrow(tcw_motif) + 1,] = list(my_sequences[my_seq-1], risk, type, i)
            }
          }
        }
      }
    #we look for other sequences than the reference
    }else{
      seq_split <- strsplit(my_sequences, "")[[my_seq]]
      for(i in 3:length(ref_split)){
        #we check if we have a TCW motif
        if(ref_split[i-2] == "t"){
          if(ref_split[i-1] == "c"){
            if(ref_split[i] == "a" || ref_split[i] == "t"){
              #we check if we have a c>t mutation in TCW motif
              if(seq_split[i-1] == "t"){
                apo_mut[nrow(apo_mut) + 1,] = list(my_sequences[my_seq-1], risk, type, i)
                nb_apo = nb_apo+1
              }else{
                #this is not a c>t mutation in TCW
                #apo_mut_noCT[nrow(apo_mut_noCT) + 1,] = list(my_sequences[my_seq-1], risk, type, i, seq_split[i-1])
              }
            }
          }
        }
      }
    }
    #we count the total number of c>t in TCW motif
    total_apo[nrow(total_apo) + 1,] = list(my_sequences[my_seq-1], risk, type, nb_apo)
  }
}

my_apo = total_apo[c(2:nrow(total_apo)),]

#table_freq = my_apo %>%
#  group_by(risk, type, c_to_t) %>%
#  summarize(total=n())

#table_freq$prop = table_freq$total/sum(table_freq$total)*100

#we save the table as a file
write.table(my_apo, paste(type,".txt"), sep="\t", row.names = F, col.names = T)




####
####
#When the analysis was performed for some types (here HPV6, 11, 16 and 18), we can regroup the files
#and perform some comparisons
####
####

setwd("D:/all_HPV")

library(ggplot2)
library(data.table)
library(dplyr)

#we get all .txt files (here, one file for each type)
temp = list.files(pattern="*.txt")
myfiles = lapply(temp, read.delim)

#we combine the different files as a dataframe
df_HPV = rbindlist(myfiles, fill=F, idcol=NULL)

#we focus on the main types
HPV_type=c("HPV6", "HPV11","HPV16", "HPV18")

#we initialize a table that will store the statistic results for each replicate
result_stat=data.frame(replicat=integer(),pval=double(), median_1=double(), median_2=double())

#we perform here 100 replicates
for(replicat in 1:100){
  print(replicat)
  #we initialize a sub_table
  fuse_HPV = data.frame()
  #for each type, we randomly retrived 50 sequences
  for(a_type in HPV_type){
    df_sub_HPV = df_HPV[df_HPV$risk==a_type,]
    df_draw = sample_n(df_sub_HPV,50)
    fuse_HPV = rbind(fuse_HPV, df_draw)
  }

  #we prepare sub tables
  df_6 = fuse_HPV[fuse_HPV$risk=="HPV6",]
  df_11 = fuse_HPV[fuse_HPV$risk=="HPV11",]
  df_16 = fuse_HPV[fuse_HPV$risk=="HPV16",]
  df_18 = fuse_HPV[fuse_HPV$risk=="HPV18",]
  #we calculate the medians
  median_6 = median(df_6$c_to_t)
  median_11 = median(df_11$c_to_t)
  median_16 = median(df_16$c_to_t)
  median_18 = median(df_18$c_to_t)
  #we perform the test
  #here, we focus in a comparison of hpv 6 and 16
  test = wilcox.test(df_6$c_to_t,df_16$c_to_t)
  #we get the pvalue of the test
  my_pval = test$p.value
  #we store the result in the dataframe
  result_stat[nrow(result_stat) + 1,] = list(replicat, my_pval, median_6, median_16)
}

#we count the number of significant tests
result_stat = result_stat[result_stat$pval<=0.05,]

#we combine the sequences based on the last replicates
combine_sub = rbind(df_6, df_11, df_16, df_18)

#we plot for each type the number of C>T mutations for each type, based on the last replicate
plot1 = ggplot(combine_sub, aes(x=risk, y=c_to_t, fill=type))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  scale_x_discrete(limits=c("HPV6","HPV11","HPV16","HPV18"))+
  #geom_hline(yintercept = 2, color="red", size=.5, linetype="dashed")+
  scale_fill_manual(values=c("#c00000","#18b410"))+
  ylim(0,8)

#show the plot
plot1
