#This creates an output folder that will hold the figures and GDC data.
#You should change the path to your own.
dir.create("C:/Users/ME/Documents/QBIO/sp24_cw/qbio490_sp24_final_PAAD/outputs")

#Below sets the output directory as our worjing directory, allowing us to use relative
#paths to save things.
setwd("C:/Users/ME/Documents/QBIO/sp24_cw/qbio490_sp24_final_PAAD/outputs")

#This section loads the packages used, and installs them if needed.
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(checkBuilt = TRUE, 
                       lib = "C:/Users/ME/AppData/Local/R/win-library/4.3", 
                       version = "3.18")
  library(TCGAbiolinks)
}
if (!require("SummarizedExperiment", quietly = TRUE)){
  BiocManager::install("SummarizedExperiment")
  library(SummarizedExperiment)
}
if (!require("DESeq2", quietly = TRUE)){
  BiocManager::install("DESeq2")
  library(DESeq2)
}
if (!require("EnhancedVolcano", quietly = TRUE)){
  BiocManager::install("EnhancedVolcano")
  library(EnhancedVolcano)
}
if (!require("maftools", quietly = TRUE)){
  BiocManager::install("maftools")
  library(maftools)
}
if (!require(survival)) {
  install.packages("survival")
  library(survival)
}

if (!require(survminer)) {
  install.packages("survminer")
  library(survminer)
}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

#This downloads and loads the TCGA PAAD clinical information to the output folder.
clin_query <- GDCquery(project = "TCGA-PAAD",
                       data.category = "Clinical",
                       data.type = "Clinical Supplement",
                       data.format = 'BCR Biotab')
GDCdownload(clin_query) #You can comment this code after running it once.
#You only need to load the data afterwards.
clinical.BCRtab.all <- GDCprepare(clin_query)
#Below you grab the sections of the clinical data we explored from our research question.
clinic <- clinical.BCRtab.all$clinical_patient_paad[-c(1,2),]
clinic = clinic[clinic$tumor_grade!="G4"&clinic$tumor_grade!="GX",]
#Removed the stage 4 and stage unknown samples. No conclusions can be drawn from the
#two G4 samples, and the one GX isn't helpful either.
rad <- clinical.BCRtab.all$clinical_radiation_paad[-c(1,2),] #Radiation treatment data.
drug <- clinical.BCRtab.all$clinical_drug_paad[-c(1,2),] #Drug treatment data.
clinic$Tumor_Sample_Barcode <- clinic$bcr_patient_barcode #Just adds a column with the same information but
# a more informative name.

#Same process with the mutation (MAF) data.
maf_query <- GDCquery(project = "TCGA-PAAD", data.category = "Simple Nucleotide Variation", access = "open", data.type = "Masked Somatic Mutation", workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query)#Again, comment after the first time.
maf <- GDCprepare(maf_query) 
maf_object <- read.maf(maf = maf, clinicalData = clinic, isTCGA = TRUE)

#And the gene expression data.
rna_query <- GDCquery(project ="TCGA-PAAD", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")

GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)
rna_clinical <- data.frame(rna_se@colData) #Clinical information of the RNA samples.
rna_clinical$Tumor_Sample_Barcode <- rna_clinical$patient

#The following is a function that borrows survminer functions to 
#plot Kaplan-Meier survival plots and calculate the significance of
#survival relationships.
KM_plot_aid <- function(data,giveTitle="", split = "tumor_grade", clinical_df = F, just_pval = F, verbose = F){
  #Data is the clinical data you are plotting from, split is the name of the column
  #with the categories you are comparing across. 
  surv_data = as.data.frame(data) #The added columns don't stay with data.
  if(grepl("-",split)|grepl("/",split)){#"-" messes with the formula, so I will add a column correcting for that and rename split.
    clean_split =gsub("-", "_", split)
    clean_split =gsub("/", "_slash_", clean_split)
    surv_data[clean_split]=surv_data[split]
    split=clean_split
  }#Formulas are sensitive to some types of names, this cleans for that.
  if(clinical_df){
    surv_data$survival_time <- ifelse(is.na(surv_data$last_contact_days_to)|
                                        (!is.na(surv_data$death_days_to)&
                                           (surv_data$death_days_to>surv_data$last_contact_days_to)),                        
                                      surv_data$death_days_to,surv_data$last_contact_days_to )
  }#survival time is how long the patient survived, in days.
  else{
    surv_data$survival_time <- ifelse(is.na(surv_data$days_to_last_follow_up)|
                                        (!is.na(surv_data$days_to_death)&
                                           (surv_data$days_to_death>surv_data$days_to_last_follow_up)),                        
                                      surv_data$days_to_death,surv_data$days_to_last_follow_up )
  }#The clinical and RNA clinical data have different column names.
  if(verbose)
    print(dim(surv_data))
  
  surv_data = surv_data[!is.na(surv_data$survival_time),]
  surv_data$survival_time <- as.numeric(surv_data$survival_time)
  surv_data$death_event <- ifelse(surv_data$vital_status=="Alive", F, T)
  #death_event is whether a patient is listed as alive or dead.
  if(verbose){#If you want to verify there are no nulls here. Didn't see issues with the PAAD data.
    print(paste("Time ",sum(is.na(surv_data$survival_time))))
    print(paste(" Death ",sum(is.na(surv_data$death_event))))
    print(paste(" Grade ",sum(is.na(surv_data[split]))))
  }
  surv_object = Surv(surv_data$survival_time, surv_data$death_event)
  formula = as.formula(surv_object ~ .)
  fit <- surv_fit(update(formula, paste("~",split)), data = surv_data)
  #survminer's function to calculate find the relationship between the split
  #and survival.
  if(just_pval){
    return(as.numeric(surv_pvalue(fit, data=surv_data)["pval"]))
  }#There is the option of just returning the p value.
  survival_plot <- ggsurvplot(fit, title = giveTitle,
                              data = surv_data,
                              pval=TRUE,
                              ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                              legend = 'right')
  KM <- survival_plot$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                axis.text = element_text(size=16),
                                                legend.title = element_text(size=14),
                                                legend.text = element_text(size=12))
  return(KM) #Or you can return the plot.
}

#Below makes the Co-oncoplots.
G1_maf = subsetMaf(maf=maf_object, tsb =clinic$Tumor_Sample_Barcode[clinic$tumor_grade=="G1"] )

G2_maf = subsetMaf(maf=maf_object, tsb =clinic$Tumor_Sample_Barcode[clinic$tumor_grade=="G2"] )

G3_maf = subsetMaf(maf=maf_object, tsb =clinic$Tumor_Sample_Barcode[clinic$tumor_grade=="G3"] )
jpeg("./Co-oncoplot_PAAD_G1vG2.jpeg")
co_onco_plot1 = coOncoplot(m1 = G1_maf, 
           m2 = G2_maf, titleFontSize = .7,geneNamefont = .4,
           m1Name = "Grade 1 Tumours", 
           m2Name = "Grade 2 Tumours",  
           borderCol = NA)
dev.off()
jpeg("./Co-oncoplot_PAAD_G2vG3.jpeg")
coOncoplot(m1 = G2_maf,
           m2 = G3_maf, titleFontSize = .7,geneNamefont = .4,
           m1Name = "Grade 2 Tumours", 
           m2Name = "Grade 3 Tumours",
           borderCol = NA)
dev.off()
#Below makes the KM plot of tumor grade and the potential mutation biomarkers.
jpeg("./KM_plot_PAAD_Survival by Tumor Grade.jpeg")
KM_plot_aid(clinic, "Pancreatic Cancer Survival by Tumor Grade", clinical_df =T)
dev.off()
jpeg("./KM_plot_PAAD_Survival by KRAS Gene.jpeg")
clinic["KRAS_mut"] = ifelse(clinic$Tumor_Sample_Barcode %in% maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol=="KRAS"],"Mutated", "Wild-type")
KM_plot_aid(clinic, "Survival by KRAS Gene", split = "KRAS_mut", clinical_df =T)
dev.off()
jpeg("./KM_plot_PAAD_Survival by TP53 Gene.jpeg")
clinic["TP53_mut"] = ifelse(clinic$Tumor_Sample_Barcode %in% maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol=="TP53"],"Mutated", "Wild-type")
KM_plot_aid(clinic, "Survival by TP53 Gene", split = "TP53_mut", clinical_df =T)
dev.off()
jpeg("./KM_plot_PAAD_Survival by SMAD4 Gene.jpeg")
clinic["SMAD4_mut"] = ifelse(clinic$Tumor_Sample_Barcode %in% maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol=="SMAD4"],"Mutated", "Wild-type")
KM_plot_aid(clinic, "Survival by SMAD4 Gene", split = "SMAD4_mut", clinical_df =T)
dev.off()
jpeg("./KM_plot_PAAD_Survival by TTN Gene.jpeg")
clinic["TTN_mut"] = ifelse(clinic$Tumor_Sample_Barcode %in% maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol=="TTN"],"Mutated", "Wild-type")
KM_plot_aid(clinic, "Survival by TTN Gene", split = "TTN_mut", clinical_df =T)
dev.off()
jpeg("./KM_plot_PAAD_Survival by CDKN2A Gene.jpeg")
clinic["CDKN2A_mut"] = ifelse(clinic$Tumor_Sample_Barcode %in% maf_object@data$Tumor_Sample_Barcode[maf_object@data$Hugo_Symbol=="CDKN2A"],"Mutated", "Wild-type")
KM_plot_aid(clinic, "Survival by CDKN2A Gene", split = "CDKN2A_mut", clinical_df =T)
dev.off()
#Below makes the co-lollipop plots of the mutations associated with survival.
jpeg("./KM_plot_PAAD_Survival by CDKN2A Gene.jpeg")
lollipopPlot2(m1 =G1_maf , 
              m2 = G3_maf, gene = "KRAS",     
              m1_name = "G1", 
              m2_name = "G3")
dev.off()
lollipopPlot2(m1 =G1_maf , 
              m2 = G3_maf, gene = "TP53",     
              m1_name = "G1", 
              m2_name = "G3")
dev.off()
lollipopPlot2(m1 =G1_maf , 
              m2 = G3_maf, gene = "CDKN2A",     
              m1_name = "G1", 
              m2_name = "G3")
dev.off()

#The following is for the analysis of the transcriptome.
rna_clinical$Radiation = ifelse(rna_clinical$Tumor_Sample_Barcode %in% rad$bcr_patient_barcode, 1, 0)
rna_clinical$Chemo = ifelse((rna_clinical$Tumor_Sample_Barcode %in% drug$bcr_patient_barcode[drug$pharmaceutical_therapy_type =="Chemotherapy"]), 1, 0)
#These one-hot encoded columns allow you to control for those treatments
#I only controlled for those two b/c these ones had enough diversity to 
#work with. 
list_columns <- sapply(rna_clinical, is.list)
rna_clinical <- rna_clinical[, !list_columns]
#Above removes nested columns that would scramble analysis.

#rna_genes holds the metadata of the genes, we will use it for labelling.
#Each row is a bit of RNA.
rna_genes<- rna_se@rowRanges@elementMetadata
rna_genes= as.data.frame(rna_genes)
rownames(rna_genes) = rna_genes$gene_id

#rna_counts is the RNA piece counts by the tumor samples.
rna_counts <- as.data.frame(rna_se@assays@data$unstranded)
rownames(rna_counts) = rownames(rna_genes)
colnames(rna_counts) = rownames(rna_clinical)

#Below is data preprocessing to remove genes for which not enough was
#picked up to make comparisons.
higher_count_mask = rowSums(rna_counts)>=20
quantile_valid_mask = vector(mode="logical", length=length(higher_count_mask))
for (i in 1:nrow(rna_counts)){
  gene_expression=rna_counts[i,]
  quantile_valid_mask[i] = quantile(gene_expression, 0.75, na.rm = TRUE)>quantile(gene_expression, 0.25, na.rm = TRUE)
}

rna_counts=rna_counts[higher_count_mask&quantile_valid_mask,]#filter out low expression and no quantiles.
rna_genes=rna_genes[higher_count_mask&quantile_valid_mask,]#filter out low expression

#This grabs the tumor grades from the clinical data.
for (bar in rna_clinical$Tumor_Sample_Barcode[rna_clinical$Tumor_Sample_Barcode %in% clinic$Tumor_Sample_Barcode]){
  rna_clinical$tumor_grade[rna_clinical$Tumor_Sample_Barcode==bar] = clinic$tumor_grade[clinic$Tumor_Sample_Barcode==bar]
}
DeSeq_grades_mask = rna_clinical$tumor_grade!="Not Reported"&rna_clinical$tumor_grade!="G2"
# Excluded G2 tumors for simplicity's sake. Making and interpreting two group DeSeqs is much easier.
rna_clinical_DeSeq = rna_clinical[DeSeq_grades_mask,]
rna_counts_DeSeq <- rna_counts[,DeSeq_grades_mask]

#Below is the setup of all the controls and the contrast.
rna_clinical_DeSeq$tumor_grade = factor(rna_clinical_DeSeq$tumor_grade)
rna_clinical_DeSeq$vital_status=factor(rna_clinical_DeSeq$vital_status)
rna_clinical_DeSeq$gender=factor(rna_clinical_DeSeq$gender)
rna_clinical_DeSeq$age_diagnosed_categ = factor(ifelse(rna_clinical_DeSeq$age_at_diagnosis/365.25>60,"Older", ifelse(rna_clinical_DeSeq$age_at_diagnosis/365.25<30,"Younger", "Middle")))
rna_clinical_DeSeq$Radiation=factor(rna_clinical_DeSeq$Radiation)
rna_clinical_DeSeq$Chemo=factor(rna_clinical_DeSeq$Chemo)

#The design controls each of the below except the last, the contrast.
design = ~ vital_status+age_diagnosed_categ+Radiation+Chemo+ gender+ tumor_grade

#Below the DESeq analysis occurs.
dds <- DESeqDataSetFromMatrix(countData = rna_counts_DeSeq,colData = rna_clinical_DeSeq, design= design)
dds_obj <- DESeq(dds)

#Below the volcano plot showing the results is made.
#plot the results.
resultsNames(dds_obj)
results <- results(dds_obj, format = "DataFrame", contrast = c("tumor_grade", "G3" , "G1"))
results <- data.frame(results)
results$gene_name <- rna_genes$gene_name
results$"-log10(padj)" <- -log10(results$padj) 
row.names(results) <- rna_genes$gene_id
results= results[!is.na(results$log2FoldChange)&!is.na(results$padj),]
#Catches nulls that might try to slip through. 
jpeg("./RNA_Volcano_PAAD_by_tumour_grade.jpeg")
EnhancedVolcano(results,
                lab = results$gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Tumour Grade: G3 vs G1',
                axisLabSize = 10,
                titleLabSize = 14, subtitleLabSize = 0,
                pointSize = 1.0,
                labSize = 4.0)
dev.off()

#Below tells you that the mutated genes were not DE'ed
paste("KRAS expression levels were significantly different between G1 and G3 samples: ", results$padj[results$gene_name=="KRAS"]<.05, sep="")
paste("TP53 expression levels were significantly different between G1 and G3 samples: ", results$padj[results$gene_name=="TP53"]<.05, sep="")
paste("SMAD4 expression levels were significantly different between G1 and G3 samples: ", results$padj[results$gene_name=="SMAD4"]<.05, sep="")
paste("TTN expression levels were significantly different between G1 and G3 samples: ", results$padj[results$gene_name=="TTN"]<.05, sep="")
paste("CDKN2A expression levels were significantly different between G1 and G3 samples: ", results$padj[results$gene_name=="CDKN2A"]<.05, sep="")

#Below grabs the results that significant and G3 is 2x and up from G1 or .5x and down.
Sig_to_Us_results = results[(results$log2FoldChange<=-1|results$log2FoldChange>=1)&results$padj<= .05,]
Sig_to_Us_results$direction_of_change = ifelse(Sig_to_Us_results$log2FoldChange>0, "Up Regulated", "Down Regulated")

#Below is a function that divides the patients into expression levels for a given gene.
#High and low are the top and bottom quartiles.
assign_expression_category <- function(gene_expression) {
  ifelse(gene_expression >= quantile(gene_expression, 0.75, na.rm = TRUE), "high",
         ifelse(gene_expression <= quantile(gene_expression, 0.25, na.rm = TRUE), "low", "medium"))
}

#Below plots KM survival plots of some of the genes most down and up regulated
#in G3 to G1 tumours.
#Biggest fold change down.
jpeg("./KM_plot_PAAD_DEFA5 effect on survival.jpeg")
rna_clinical$DEFA5_level =t(assign_expression_category(rna_counts[rna_se@rowRanges$gene_id[rna_se@rowRanges$gene_name == "DEFA5"],]))
KM_plot_aid(rna_clinical[rna_clinical$DEFA5_level!="medium",], "DEFA5 effect on survival", split = "DEFA5_level")
dev.off()
jpeg("./KM_plot_PAAD_AC025756.1 effect on survival.jpeg")
rna_clinical$AC025756.1_level =t(assign_expression_category(rna_counts[rna_se@rowRanges$gene_id[rna_se@rowRanges$gene_name == "AC025756.1"],]))
KM_plot_aid(rna_clinical[rna_clinical$AC025756.1_level!="medium",], "AC025756.1 effect on survival", split = "AC025756.1_level")
dev.off()
jpeg("./KM_plot_PAAD_APOA2 effect on survival.jpeg")
rna_clinical$APOA2_level =t(assign_expression_category(rna_counts[rna_se@rowRanges$gene_id[rna_se@rowRanges$gene_name == "APOA2"],]))
KM_plot_aid(rna_clinical[rna_clinical$APOA2_level!="medium",], "APOA2 effect on survival", split = "APOA2_level")
dev.off()

#Up regulated
jpeg("./KM_plot_PAAD_CASP14 effect on survival.jpeg")
rna_clinical$CASP14_level =t(assign_expression_category(rna_counts[rna_se@rowRanges$gene_id[rna_se@rowRanges$gene_name == "CASP14"],]))
KM_plot_aid(rna_clinical[rna_clinical$CASP14_level!="medium",], "CASP14 effect on survival", split = "CASP14_level")
dev.off()
jpeg("./KM_plot_PAAD_MIR205HG effect on survival.jpeg")
rna_clinical$MIR205HG_level =t(assign_expression_category(rna_counts[rna_se@rowRanges$gene_id[rna_se@rowRanges$gene_name == "MIR205HG"],]))
KM_plot_aid(rna_clinical[rna_clinical$MIR205HG_level!="medium",], "MIR205HG effect on survival", split = "MIR205HG_level")
dev.off()
jpeg("./KM_plot_PAAD_KRT16P6 effect on survival.jpeg")
rna_clinical$KRT16P6_level =t(assign_expression_category(rna_counts[rna_se@rowRanges$gene_id[rna_se@rowRanges$gene_name == "KRT16P6"],]))
KM_plot_aid(rna_clinical[rna_clinical$KRT16P6_level!="medium",], "KRT16P6 effect on survival", split = "KRT16P6_level")
dev.off()

#Below finds the significance of every gene in Sig_to_Us_results to survival.
clutter_rna_clin = as.data.frame(rna_clinical)
survival_effect = vector(length=nrow(Sig_to_Us_results))
count = 1
for(i in 1:(length(Sig_to_Us_results$gene_name))){
  name = Sig_to_Us_results$gene_name[i]
  gene_expression = rna_counts[rna_se@rowRanges$gene_id[rna_se@rowRanges$gene_name == name],]
  if(nrow(gene_expression)>1){#Some genes have multiple bits of RNA in the data.
    ID = rownames(Sig_to_Us_results)[i] # I say bits b/c it's not all mRNA.
    gene_expression = rna_counts[ID,]
    name = capture.output(cat(c(name,".",ID),sep=""))
  }#Formulas are particular about the char's in inputs. They don't like "()".
  clutter_rna_clin[paste(name,sep="","_level")]= t(assign_expression_category(gene_expression))
  survival_effect[count]= KM_plot_aid(clutter_rna_clin[clutter_rna_clin[paste(name,sep="","_level")]!="medium",], "", split = paste(name,sep="","_level"),just_pval = T)
  count=count+1
  clutter_rna_clin <- clutter_rna_clin[, 1:ncol(rna_clinical)]
}
Sig_to_Us_results$survival_relationship_sig = survival_effect

#Below captures all the genes with expression levels with a significant association with survival.
results_Surv_effect <- Sig_to_Us_results[Sig_to_Us_results$survival_relationship_sig<0.05,]

#None of the potential mutation biomarkers have expression levels that affect survival.
paste("KRAS expression levels affected survival: ", "KRAS" %in% results_Surv_effect$gene_name, sep="")
paste("TP53 expression levels affected survival: ", "TP53" %in% results_Surv_effect$gene_name, sep="")
paste("SMAD4 expression levels affected survival: ", "SMAD4" %in% results_Surv_effect$gene_name, sep="")
paste("TTN expression levels affected survival: ", "TTN" %in% results_Surv_effect$gene_name, sep="")
paste("CDKN2A expression levels affected survival: ", "CDKN2A" %in% results_Surv_effect$gene_name, sep="")

#Below writes out those 470 gene entries for your viewing pleasure.
ordered_log <- order(results_Surv_effect$direction_of_change)
results_Surv_effect <- results_Surv_effect[ordered_log, ]
write.csv(results_Surv_effect, file="/PAAD_genes_affected_Surv.csv")