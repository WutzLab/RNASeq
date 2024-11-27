directory <- "../../"
file_list <- grep("_counts", list.files(directory), value=TRUE)
ordered_files <- c(
"HATX3_counts",          # WT Ctrl
"HATX8_counts",
"HATX12_counts",
"HATX3_DOX_counts",      # WT Dox
"HATX8_DOX_counts",
"HATX12_DOX_counts",
"U1_10_counts",          # Ubn1 Ctrl
"U1_13_counts",
"U1_36_counts",
"U1_10_DOX_counts",      # Ubn1 Dox
"U1_13_DOX_counts",
"U1_36_DOX_counts",
"U9_counts",             # Ubn2 Ctrl
"U14_counts",
"Ue38_counts",
"U9_DOX_counts",         # Ubn2 Dox
"U14_DOX_counts",
"Ue38_DOX_counts",
"U12_10_counts",         # Ubn1/2 Ctrl
"U12_27_counts",
"U12_82_counts",
"U12_10_DOX_counts",     # Ubn1/2 Dox
"U12_27_DOX_counts",
"U12_82_DOX_counts",
"H23_counts",            # Hira Ctrl
"H18_counts",
"HE46_counts",
"H23_DOX_counts",        # Hira Dox
"H18_DOX_counts",
"HE46_DOX_counts")

fls <- ordered_files
condition <- ifelse(grepl("DOX", fls), "Dox", "Ctrl")
Ctrl_files <- ordered_files[condition == "Ctrl"]
Dox_files <- ordered_files[condition == "Dox"]

condition <- ifelse(grepl("DOX", fls), "Dox", "Ctrl")
genotype <- ifelse(grepl("HATX", fls), "WT",
            ifelse(grepl("U1_", fls), "Ubn1",
            ifelse(grepl("U12_", fls), "Ubn12",
            ifelse(grepl("U", fls), "Ubn2",
            ifelse(grepl("H", fls), "Hira", "NA")))))

WT_files <- ordered_files[genotype == "WT"]
Ubn1_files <- ordered_files[genotype == "Ubn1"]
Ubn2_files <- ordered_files[genotype == "Ubn2"]
Ubn12_files <- ordered_files[genotype == "Ubn12"]
Hira_files <- ordered_files[genotype == "Hira"]

concatstr <- function(..., sep='') {  # joins two or more strings
   paste(..., sep=sep, collapse=sep)
}

library("DESeq2")

for(gt in c("WT","Ubn1","Ubn2","Ubn12","Hira")) {
   print(concatstr("Results for ",gt," +Dox/-Dox"))
   fls <- ordered_files[genotype == gt]
   condition <- ifelse(grepl("DOX", fls), "Dox", "Ctrl")
   celline <- unlist(ifelse(grepl("DOX", fls), strsplit(fls, split="_DOX_counts"), strsplit(fls, split="_counts")))

   clonenumber <- c()
   for (i in 1:length(fls)) {
      if(gt == "WT") {
            clnum <- grep(celline[i], c("HATX3", "HATX8", "HATX12"))
         } else if(gt == "Ubn1") {
            clnum <- grep(celline[i], c("U1_10", "U1_13", "U1_36"))
         } else if(gt == "Ubn2") {
            clnum <- grep(celline[i], c("U9", "U14", "Ue38"))
         } else if(gt == "Ubn12") {
            clnum <- grep(celline[i], c("U12_10", "U12_27", "U12_82"))
         } else if(gt == "Hira") {
            clnum <- grep(celline[i], c("H23", "HE46", "H18"))
         } else {
            clnum <- 0
         }
      clonenumber <- c(clonenumber, clnum)
      }
   samplename <- c()
   for (i in 1:length(fls)) {
      samplename <- c(samplename, concatstr(gt, "_", as.character(clonenumber[i]), condition[i]))
   }

   cat("genotype\tclone number\tcondition\tcell line\tsample name\tfilename\n")
   for (i in 1:length(fls)) {
      cat(concatstr(gt, concatstr("#",as.character(clonenumber[i])), condition[i], celline[i], samplename[i], fls[i], sep='\t\t'))
      cat("\n")
   }

   print(concatstr(as.character(length(fls)),"filenames of type", typeof(fls), sep=" "))

   sampleTable <- data.frame(sampleName = samplename,
                             fileName = fls,
                             condition = condition)
   sampleTable$condition <- factor(sampleTable$condition)
   ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                          directory = directory,
                                          design= ~ condition)
   ddsHTSeq$genotype <- factor(ddsHTSeq$condition, levels = c("Ctrl", "Dox"))
   smallestGroupSize <- 3                                          # prefiltering
   keep <- rowSums(counts(ddsHTSeq) >= 10) >= smallestGroupSize
   ddsHTSeq <- ddsHTSeq[keep,]
   ddsHTSeq <- DESeq(ddsHTSeq)
   res05 <- results(ddsHTSeq, contrast=c("condition","Dox","Ctrl"), alpha=0.05)
   resOrdered <- res05[order(res05$pvalue),]
   print(concatstr("./",gt,"_Dox_vs_Ctrl_results.csv"))
   write.csv(as.data.frame(resOrdered), 
             file=concatstr("./",gt,"_Dox_vs_Ctrl_results.csv"))
   resSig <- subset(resOrdered, padj <= 0.05)
   print(concatstr("./",gt,"_Dox_vs_Ctrl_results05.csv"))
   write.csv(as.data.frame(resSig), 
             file=concatstr("./",gt,"_Dox_vs_Ctrl_results05.csv"))
   resultsNames(res05)
   summary(res05)
   print(concatstr(as.character(sum(res05$padj <= 0.05, na.rm=TRUE)),"significantly differentially expressed genes (adj.p-value <= 0.05).", sep=" "))
   resSig_1_5foldUP <- subset(resSig, log2FoldChange >= (log(1.5)/log(2)))
   print(concatstr(as.character(sum(resSig_1_5foldUP$padj < 0.05, na.rm=TRUE)),"1.5-fold or more upregulated genes (adj.p-value < 0.05).", sep=" "))
   resSig_2foldUP <- subset(resSig, log2FoldChange >= 1)
   print(concatstr(as.character(sum(resSig_2foldUP$padj < 0.05, na.rm=TRUE)),"2-fold or more upregulated genes (adj.p-value < 0.05).", sep=" "))
   resSig_1_5foldDOWN <- subset(resSig, log2FoldChange <= (-log(1.5)/log(2)))
   print(concatstr(as.character(sum(resSig_1_5foldDOWN$padj < 0.05, na.rm=TRUE)),"1.5-fold or more downregulated genes (adj.p-value < 0.05).", sep=" "))
   resSig_2foldDOWN <- subset(resSig, log2FoldChange <= -1)
   print(concatstr(as.character(sum(resSig_2foldDOWN$padj < 0.05, na.rm=TRUE)),"2-fold or more downregulated genes (adj.p-value < 0.05).", sep=" "))

   plotMA(res05, ylim=c(-3,3))
   # plotCounts(ddsHTSeq, gene=which.min(res05$padj), intgroup="genotype")
}



# resLFC <- lfcShrink(ddsHTSeq, coef="genotype_Ubn12_vs_WT") # effect size adjust log fold change

# res05_Ubn12 <- results(ddsHTSeq, contrast=c("genotype","Ubn12","WT"), alpha=0.05)
# resultsNames(res05_Ubn12)
# summary(res05_Ubn12)
# sum(res05_Ubn12$padj < 0.05, na.rm=TRUE)
# resOrdered <- res05_Ubn12[order(res05_Ubn12$pvalue),]
# resSig <- subset(resOrdered, padj < 0.1)
# write.csv(as.data.frame(resOrdered), file="Ubn12_vs_WT_results.csv")
# plotMA(res05_Ubn12, ylim=c(-3,3))
# plotCounts(ddsHTSeq, gene=which.min(res05_Ubn12$padj), intgroup="genotype")


# SAMPLE ORDER IN PYTHON ANALYSIS
# genotypes = ["WT", "Ubn1", "Ubn2", "Ubn1/2", "Hira", "Ubn1_2"]
# conditions = ["Ctrl", "Dox"]
# cellLines = { "WT": ["HATX3", "HATX8", "HATX12"],
#              "Ubn1": ["U1_10", "U1_13", "U1_36"],
#              "Ubn2": ["U9", "U14", "Ue38"],
#              "Ubn1/2": ["U12_10", "U12_27", "U12_82"],
#              "Ubn1_2": ["U12_8"],
#              "Hira": ["H23", "HE46", "H18"] }
