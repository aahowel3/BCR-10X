#provides VDJ gene family assignments NS CDR3 sequences but not coordinates
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/")
data=read.csv("filtered_contig_annotations.csv")

#remove cell barcodes with incomplete chains
library(tidyverse)
full_length=data %>%
  group_by(barcode) %>%
  filter(n() == 2) %>%
  ungroup()

####################################################################################################################
#updated graphs by V gene not entire region - load JSON file with VDJ region details
library(jsonlite)
myData <- jsonlite::fromJSON(txt="all_contig_annotations.json")
#reduce down to only viable cell barcodes
myData=myData[myData$contig_name %in% full_length$contig_id,]
myData=myData[order(match(myData$contig_name,full_length$contig_id)),]

#convert subcolumns to regular columns
#you will need this for the for loops
myData$fwr1_start = myData$fwr1$start
myData$fwr1_stop = myData$fwr1$stop
myData$fwr2_start = myData$fwr2$start
myData$fwr2_stop = myData$fwr2$stop
myData$fwr3_start = myData$fwr3$start
myData$fwr3_stop = myData$fwr3$stop
myData$fwr4_start = myData$fwr4$start
myData$fwr4_stop = myData$fwr4$stop
myData$cdr1_start = myData$cdr1$start
myData$cdr1_stop = myData$cdr1$stop
myData$cdr2_start = myData$cdr2$start
myData$cdr2_stop = myData$cdr2$stop

###########################################################################################################
#for each bam qname entry  - find corresponding ref seq 
#cigar orient them, json feature extract L+V region, count mutations
library(pkgmaker)
library(tidyverse)

pair_bam=read.csv("concat_ref_pairwise_MDtag.out",header=FALSE)
bamtable <- data.frame(matrix(pair_bam$V1, ncol = 4, byrow = TRUE))
diff_table <- separate(bamtable, X1, into = c("qname", "tags"), sep = "^\\S*\\K\\s+")
diff_table=diff_table[diff_table$qname %in% full_length$contig_id,]
diff_table=diff_table[order(match(diff_table$qname,full_length$contig_id)),]

#apply strsplit to every row in the two designated columns
fun <- function(x, character = FALSE) {
  strsplit(x, "")[[1]]
}
diff_table$query = sapply(diff_table$X2, fun)
diff_table$ref = sapply(diff_table$X4, fun)

#######################################################################################
#convert wide format sequences into variable-length long format mutations against ref genome per position
minime2=diff_table[c("qname", "X2", "X4", "query","ref")]
#data format is the same copy over chain column
minime2$chain=full_length$chain
minime2$chain_type=ifelse(minime2$chain=="IGH","heavy","light")

#numbering schema for query column
#cannot be the same as cellranger listed position as insertions will change numbering
for(r in 1:length(minime2$query)) {
  minime2$schema[[r]]=ifelse(minime2$query[[r]] != "-" , cumsum(minime2$query[[r]] != "-" ), NA)
}

minime2$schema2 <- sapply(minime2$schema, paste, collapse = ",")
minime2$query2 <- sapply(minime2$query, paste, collapse = ",")
minime2$ref2 <- sapply(minime2$ref, paste, collapse = ",")

#convert from wide to long
tall_df_1 = minime2 %>%
  separate_rows(c(query2,schema2,ref2), sep=',')
#suequential numbering by condiiton
tall_df_1 = tall_df_1 %>% group_by(qname) %>% mutate(POS = seq_len(n()))


#add indexing column to reference which labeled sections to extract 
tall_df_1 = tall_df_1 %>%   ungroup(qname)  %>%   
  mutate(index = 1:length(tall_df_1$POS))

#if you get an error here its likely bc you're trying to do this operation on a truncated dataset with only mutations showing
#which messes up the indexing asks and breaks the loop
tall_df_1$label = 0
#for CDR3 label ONLY will need to repeat function
for(contig in unique(tall_df_1$qname)) { 
  start_num = myData[which(myData$contig_name==contig),"cdr3_start"]
  end_num = myData[which(myData$contig_name==contig),"cdr3_stop"]
  start_num = as.numeric(start_num + 1)
  end_num = as.numeric(end_num)
  start <- tall_df_1[which(tall_df_1$schema2 == start_num & tall_df_1$qname == contig), "index"]
  end   <- tall_df_1[which(tall_df_1$schema2 == end_num & tall_df_1$qname == contig), "index"]
  for(i in tall_df_1$index) {
    if(tall_df_1$index[i] <= as.numeric(end) & tall_df_1$index[i] >= as.numeric(start)) {
      tall_df_1$label[i] = "CDR3"
    }
    if(tall_df_1$index[i] <= as.numeric(start) & tall_df_1$index[i] >= as.numeric(end)) {
      next   
    }
  } 
}

for(contig in unique(tall_df_1$qname)) { 
  start_num = myData[which(myData$contig_name==contig),"fwr1_start"]
  end_num = myData[which(myData$contig_name==contig),"fwr1_stop"]
  start_num = start_num + 1
  start <- tall_df_1[which(tall_df_1$schema2 == start_num & tall_df_1$qname == contig), "index"]
  end   <- tall_df_1[which(tall_df_1$schema2 == end_num & tall_df_1$qname == contig), "index"]
  for(i in tall_df_1$index) {
    if(tall_df_1$index[i] <= as.numeric(end) & tall_df_1$index[i] >= as.numeric(start)) {
      tall_df_1$label[i] = "FWR1"
    }
    if(tall_df_1$index[i] <= as.numeric(start) & tall_df_1$index[i] >= as.numeric(end)) {
      next   
    }
  } 
}

for(contig in unique(tall_df_1$qname)) { 
  start_num = myData[which(myData$contig_name==contig),"fwr2_start"]
  end_num = myData[which(myData$contig_name==contig),"fwr2_stop"]
  start_num = start_num + 1
  start <- tall_df_1[which(tall_df_1$schema2 == start_num & tall_df_1$qname == contig), "index"]
  end   <- tall_df_1[which(tall_df_1$schema2 == end_num & tall_df_1$qname == contig), "index"]
  for(i in tall_df_1$index) {
    if(tall_df_1$index[i] <= as.numeric(end) & tall_df_1$index[i] >= as.numeric(start)) {
      tall_df_1$label[i] = "FWR2"
    }
    if(tall_df_1$index[i] <= as.numeric(start) & tall_df_1$index[i] >= as.numeric(end)) {
      next   
    }
  } 
}

##
for(contig in unique(tall_df_1$qname)) { 
  start_num = myData[which(myData$contig_name==contig),"fwr3_start"]
  end_num = myData[which(myData$contig_name==contig),"fwr3_stop"]
  start_num = start_num + 1
  start <- tall_df_1[which(tall_df_1$schema2 == start_num & tall_df_1$qname == contig), "index"]
  end   <- tall_df_1[which(tall_df_1$schema2 == end_num & tall_df_1$qname == contig), "index"]
  for(i in tall_df_1$index) {
    if(tall_df_1$index[i] <= as.numeric(end) & tall_df_1$index[i] >= as.numeric(start)) {
      tall_df_1$label[i] = "FWR3"
    }
    if(tall_df_1$index[i] <= as.numeric(start) & tall_df_1$index[i] >= as.numeric(end)) {
      next   
    }
  } 
}

for(contig in unique(tall_df_1$qname)) { 
  start_num = myData[which(myData$contig_name==contig),"fwr4_start"]
  end_num = myData[which(myData$contig_name==contig),"fwr4_stop"]
  start_num = start_num + 1
  start <- tall_df_1[which(tall_df_1$schema2 == start_num & tall_df_1$qname == contig), "index"]
  end   <- tall_df_1[which(tall_df_1$schema2 == end_num & tall_df_1$qname == contig), "index"]
  for(i in tall_df_1$index) {
    if(tall_df_1$index[i] <= as.numeric(end) & tall_df_1$index[i] >= as.numeric(start)) {
      tall_df_1$label[i] = "FWR4"
    }
    if(tall_df_1$index[i] <= as.numeric(start) & tall_df_1$index[i] >= as.numeric(end)) {
      next   
    }
  } 
}

for(contig in unique(tall_df_1$qname)) { 
  start_num = myData[which(myData$contig_name==contig),"cdr1_start"]
  end_num = myData[which(myData$contig_name==contig),"cdr1_stop"]
  start_num = start_num + 1
  start <- tall_df_1[which(tall_df_1$schema2 == start_num & tall_df_1$qname == contig), "index"]
  end   <- tall_df_1[which(tall_df_1$schema2 == end_num & tall_df_1$qname == contig), "index"]
  for(i in tall_df_1$index) {
    if(tall_df_1$index[i] <= as.numeric(end) & tall_df_1$index[i] >= as.numeric(start)) {
      tall_df_1$label[i] = "CDR1"
    }
    if(tall_df_1$index[i] <= as.numeric(start) & tall_df_1$index[i] >= as.numeric(end)) {
      next   
    }
  } 
}

for(contig in unique(tall_df_1$qname)) { 
  start_num = myData[which(myData$contig_name==contig),"cdr2_start"]
  end_num = myData[which(myData$contig_name==contig),"cdr2_stop"]
  start_num = start_num + 1
  start <- tall_df_1[which(tall_df_1$schema2 == start_num & tall_df_1$qname == contig), "index"]
  end   <- tall_df_1[which(tall_df_1$schema2 == end_num & tall_df_1$qname == contig), "index"]
  for(i in tall_df_1$index) {
    if(tall_df_1$index[i] <= as.numeric(end) & tall_df_1$index[i] >= as.numeric(start)) {
      tall_df_1$label[i] = "CDR2"
    }
    if(tall_df_1$index[i] <= as.numeric(start) & tall_df_1$index[i] >= as.numeric(end)) {
      next   
    }
  } 
}


#truncate to plot only mutated positions 
#DO NOT SUBSET TALLDF1 UNTIL AFTER APPLYING LABLES
tall_df_1$mut=ifelse(tall_df_1$query2 != tall_df_1$ref2,"mu","invar")
tall_df_2=tall_df_1[tall_df_1$mut == "mu", ]
tall_df_2=tall_df_2[tall_df_2$label != "0", ]

#interchange if you wanna look at heavy or light
tall_df_2=tall_df_2[tall_df_2$chain_type == "heavy", ]
tall_df_2=tall_df_2[tall_df_2$chain_type == "light", ]

#plot frequency of mutations based on functional region
g <- ggplot(tall_df_2, aes(qname))
g + geom_bar(aes(fill = label)) + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_x_discrete(limits = unique(tall_df_2$qname))

#subset even further just CDR3 - but need to do prop of mut/length otherwise its not normalized
tall_df_2=tall_df_2[tall_df_2$label == "CDR3", ]

#use frequency of those entries bybarcode name to get a number of mutations table
#and then glue full sequence CDR3 seq onto it 
#and do CDR3 mutations PROPRTIONAL to the CDR3 length 
#more accurate idea of the number of mutations across clonotypes
proportion=as.data.frame(table(tall_df_2$qname))
proportion$Var1 = as.character(proportion$Var1)
proportion = proportion[order(match(proportion$Var1,myData$contig_name)),]
#onyl get entries from proportion that appear in full dataset
df.2.sub <- myData[myData$contig_name %in% proportion$Var1,]
proportion$length=nchar(df.2.sub$cdr3_seq)
proportion$prop = proportion$Freq / proportion$length

ggplot(proportion, aes(x=Var1,y=prop)) +
  geom_bar(stat = "identity") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_x_discrete(limits = proportion$Var1)


###################################################################################################
##plot CDR3 legnth dist as AAs
##quick quick plots just to share in report 
full_length$CDR3_length=nchar(full_length$cdr3)
#tweak alt column chain2 just to get IGH to plot first so its length is easier to see
full_length$chain2=ifelse(full_length$chain=="IGH","zIGH", full_length$chain)
full_length$unique=duplicated(full_length$raw_clonotype_id)
truncated=full_length[full_length$unique == "FALSE", ]

#CDR3 AA lenght plot 
ggplot(full_length, aes(x=barcode,y=CDR3_length, fill=chain)) + 
  geom_bar(position="stack", stat = "identity") + 
  scale_x_discrete(limits=full_length$barcode,breaks=truncated$barcode, labels=levels(as.factor(full_length$raw_clonotype_id))) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 

#v gene usage
ggplot(full_length, aes(v_gene)) + 
  geom_bar() +
  #  geom_bar(position="stack", stat = "identity") + 
  #  scale_x_discrete(limits=full_length$barcode,breaks=truncated$barcode, labels=levels(as.factor(full_length$raw_clonotype_id))) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 