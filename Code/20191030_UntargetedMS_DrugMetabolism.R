# Untargeted Mass Spectrometry - Data Analysis - Drug Metabolism
# AUTHOR: Alan K Jarmusch (ajarmusch@ucsd.edu)
# DATE UPDATED: 2019-10-30

# R Libraries
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(viridis))
suppressMessages(library(cowplot))
suppressMessages(library(pcaMethods))

# Data Formatting ---------------------------------------------------------

# INPUT: data to reformat from MZmine2 output
df <- read.csv("GNPS_quant.csv", header=TRUE, sep=",")
  df <- df[,-length(df)]
  # Optional: Cut input data based on RT
    #df <- subset(df, df[,3] >= 0.5 & df[,3] <= 6.0)
  # Optional: Cut input data based on m/z
    #df <- subset(df, df[,2] >= 100 & df[,2] <= 600)

# bind m/z and RT into a feature ID
MS1_feature <- paste(df[,1],df[,2],df[,3],sep="_")
  df <- df[,-c(1:10)]
  df <- cbind(MS1_feature, df)

# reformat file names output from MZmine2 to metadata format
colnames <- read.csv("GNPS_quant.csv", header=FALSE, nrow=1)
  colnames <- colnames[,-length(df)]
  colnames <- colnames[,-c(1:10)]
  colnames <- gsub(" Peak area", "", as.matrix(colnames))
  colnames <- cbind("MS1_feature",colnames)

# final reformat step
colnames(df) <- colnames

# transpose
t_df <- transpose(df)

# calculate sample sum, sample mean, sample standard deviation, and MS features with non-zero peak areas
# format for calculations
df_test <- as.matrix(as.data.frame(lapply(t_df[-1,], as.numeric)))
  # calculate
  sum <- c("sum",rowSums(df_test,na.rm=TRUE))
  mean <- c("mean",rowMeans(df_test, na.rm=TRUE))
  sd <- c("sd",rowSds(df_test, na.rm=TRUE))
  nonzero_features <- c("nonzero_features",length(df_test[1,]) - rowCounts(df_test, value = 0 , na.rm=TRUE))
  # add calculated values to data
  out <- cbind(sum,mean,sd,nonzero_features,t_df)
  out <- cbind(rownames(out),out)
  colnames(out) <- out[1,]
  out <- out[-1,]
  out[,1] <- colnames[-1]
  colnames(out)[1:5] <- c("filename","sum","mean","sd","nonzero_features")

# OUTPUT: formatted table
write.table(out,"20190710_DrugMetabolism_Formatted.txt", sep="\t", row.names=FALSE)

rm(list=ls())

# Data Quality Assessment -------------------------------------------------

# INPUT: Formatted Data
df <- read.delim("20190710_DrugMetabolism_Formatted.txt", sep="\t", header=TRUE)

# Compute QC/QA information on standards
sulfamethizole <- df[,grep("X9_271.031631288304_2.45372310691824", colnames(df))]
sulfamethizole_split <- df[,grep("X5896_271.031858170132_2.45018746189025", colnames(df))]
sulfamethazine <- df[,grep("X158_279.088420881166_2.47601433480608", colnames(df))]
sulfadimethoxine <- df[,grep("X159_285.020546247722_2.63319565217391", colnames(df))]
sulfachloropyridazine <- df[,grep("X160_311.081214930407_3.06229951907131", colnames(df))]
number_sulfamethizole <- grep("X9_271.031631288304_2.45372310691824", colnames(df))
number_sulfamethizole_split <- grep("X5896_271.031858170132_2.45018746189025", colnames(df))
number_sulfamethazine <- grep("X158_279.088420881166_2.47601433480608", colnames(df))
number_sulfadimethoxine <- grep("X159_285.020546247722_2.63319565217391", colnames(df))
number_sulfachloropyridazine <- grep("X160_311.081214930407_3.06229951907131", colnames(df))

# Create table of just standards
IS_plot <- data.frame(df[,1],sulfamethizole,sulfamethizole_split,sulfamethazine,sulfadimethoxine,sulfachloropyridazine)
  IS_plot <- data.frame(lapply(IS_plot, as.character), stringsAsFactors=FALSE)
  colnames(IS_plot)[1] <- "filename" 

# Load Metadata
metadata <- read.delim("Metadata.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

# Merge standards data with metadata
IS_plot_merged <- merge(metadata, IS_plot, by="filename")

  # Troubleshoot merge
  # troubleshoot <- IS_plot[!IS_plot$filename %in% IS_plot_merged$MS_filename_mzXML,]

# Merge standards data with specific values from MS1 data matrix
IS_plot_merged <- merge(df[,1:5], IS_plot_merged, by="filename")

# OUTPUT: write matrix to file (.csv)
write.csv(IS_plot_merged,"20190710_DrugMetabolism_ISinfo.csv", row.names = FALSE)

# PLOT
IS_plot_merged <- fread("20190710_DrugMetabolism_ISinfo.csv", sep=",", header = TRUE)

# Wide to Long format
IS_plot_merged_long <- gather(IS_plot_merged, standard, peak_area, 220:ncol(IS_plot_merged))

ggplot_IS <- ggplot(data=subset(IS_plot_merged_long,IS_plot_merged_long$ATTRIBUTE_sample_type != "not applicable"), aes(x=as.factor(ATTRIBUTE_MS_SampleType_SampleSub), y=as.numeric(peak_area)))+
  facet_wrap(~standard, scales = "free_y")+
  geom_jitter(cex=2, pch=20)+
  geom_boxplot(alpha=0, width = 0.5, outlier.alpha = 0)+
  theme_minimal()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
        axis.text=element_text(colour="black",size=8),
        axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
        axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
        axis.title=element_text(colour="black",size=8),
        strip.text.x = element_text(colour="black",size=8),
        aspect.ratio=1,
        legend.position="none") +
  labs(x="Sample_Type", y="Peak Area of Standard")
print(ggplot_IS)

ggplot_IS <- ggplot(data=subset(IS_plot_merged_long,IS_plot_merged_long$ATTRIBUTE_MS_SampleType_SampleSub == "blood_plasma"),
                    aes(x=as.Date(ATTRIBUTE_MS_dateacquired, format = "%m/%d/%Y"), y=as.numeric(peak_area)))+
  facet_wrap(~standard, scales = "free")+
  geom_jitter(width=0.3, cex=1, pch=16)+
  theme_minimal()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
        axis.text=element_text(colour="black",size=8),
        axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
        axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
        axis.title=element_text(colour="black",size=8),
        strip.text.x = element_text(colour="black",size=8),
        aspect.ratio=0.5,
        legend.position="none") +
  labs(x="Sample_Type", y="Peak Area of Standard")
print(ggplot_IS)

# Remove Standards Prior to Further Analysis (remove # of columns = # internal standards)
dim(df)
df <- df[,-c(number_sulfamethizole,number_sulfamethizole_split,number_sulfamethazine,number_sulfadimethoxine,number_sulfachloropyridazine)]
dim(df)

# OUTPUT: merge normalized data matrix with metdata and output
out_unnormalize <- merge(metadata, df, by="filename")
write.csv(out_unnormalize,"20190710_DrugMetabolism_MS1DataMatrix_NoNorm_NoStds_wMetadata.csv", row.names=FALSE)

rm(list=ls())


# Unsupervised Multivariate Statistical Analysis - PCA --------------------

# Urine
  # INPUT: Data
  df <- fread("20190710_DrugMetabolism_MS1DataMatrix_NoNorm_NoStds_wMetadata.csv", header=TRUE, sep=",")
    df <- subset(df,df$ATTRIBUTE_MS_SampleType_SampleSub == "urine")
  
  # Order study day
  df$ATTRIBUTE_Study_Day <- factor(df$ATTRIBUTE_Study_Day, levels = c('0','1','2','3','4','5','6','7','8','9','15','38') )
  df$ATTRIBUTE_Time_Point_Mins <- factor(df$ATTRIBUTE_Time_Point_Mins, levels = c('0','5','30','60','120','240','300','360','480','not applicable','FALSE'))
  df$ATTRIBUTE_Study_Day_Text <- factor(df$ATTRIBUTE_Study_Day_Text, levels = c('enrollment','pre-antibiotic PK','antibiotic Rx','post-antibiotic PK',
                                                                                '1 week post-antibiotic Rx','1 month post-antibiotic Rx'))
  # Normalize by creatinine
  df <- cbind(df[,1:219],sweep(df[,-c(1:219)], df$X283_114.062927408424_0.302783602150537, MARGIN=1, "/"))
    df <- df[complete.cases(df), ]
  
  # PCA
  mypca <- pca(df[,-c(1:219)], method="nipals", center=TRUE, nPcs = 5, scale = "pareto")
    df <- cbind(scores(mypca),df)
    loadings <- as.data.frame(loadings(mypca))
    rownames <- rownames(loadings)
    loadings <- cbind(rownames,loadings)
    rownames(loadings) <- NULL
    head(loadings)
    test_loading_names <- loadings
    test_loading_names <- cbind(substr(test_loading_names$rownames, 0, 9),loadings)
    colnames(test_loading_names)[1] <- "mz"
    test_loading_names <- 
      separate(test_loading_names, col="mz", remove= TRUE, "_", into=c("feature_number","mz"))
    head(test_loading_names)
  
  # Define Variable to Plot
  plot_value <- df$ATTRIBUTE_Time_Point_Mins
  
  # Plot
  score_plot <- ggplot(df, aes(as.numeric(PC1), as.numeric(PC2)))+
    facet_wrap(~ATTRIBUTE_Study_Day_Text, ncol=1)+
    geom_point(aes(colour=as.factor(plot_value)), pch=16, cex=1, alpha=1.0) +
    scale_color_brewer(palette = "Dark2")+
    theme_minimal()+
    theme(panel.border = element_rect(colour ="black", size=0.5, linetype = "solid", fill = NA),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(colour="black", size=6),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x=element_text(colour="black",size=6),
          aspect.ratio=1.0,
          legend.position="right",
          legend.title=element_text(colour="black", size=6),
          legend.text=element_text(colour="black", size=6) ) +
    labs(x="PC1", y="PC2")
  print(score_plot)
  ggsave("Plot_Urine_PCA.pdf", width = 4, height = 4.0, units = "in", dpi=300, useDingbats=FALSE)
  
  # Plot
  density_score_plot <- ggplot(df, aes(as.numeric(PC1)))+
    facet_wrap(~ATTRIBUTE_Study_Day_Text, ncol=1)+
    geom_density(aes(colour=as.factor(plot_value),fill=as.factor(plot_value)), alpha=0.5, size=0.5) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal()+
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(colour="black", size=6),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x=element_text(colour="black",size=6),
          aspect.ratio=1.5,
          legend.position="none",
          legend.title=element_text(colour="black", size=6),
          legend.text=element_text(colour="black", size=6) ) +
    labs(x="PC1 score values", y="Density")
  print(density_score_plot)
  ggsave("Plot_Urine_Density.pdf", width = 2.0, height = 4.0, units = "in", dpi=300, useDingbats=FALSE)

  rm(list=ls())
  
# Plasma (Blood)

  # INPUT: Data
  df <- fread("20190710_DrugMetabolism_MS1DataMatrix_NoNorm_NoStds_wMetadata.csv", header=TRUE, sep=",")
    df <- subset(df,df$ATTRIBUTE_MS_SampleType_SampleSub == "blood_plasma")
  
  # Order study day
  df$ATTRIBUTE_Study_Day <- factor(df$ATTRIBUTE_Study_Day, levels = c('0','1','2','3','4','5','6','7','8','9','15','38') )
  df$ATTRIBUTE_Time_Point_Mins <- factor(df$ATTRIBUTE_Time_Point_Mins, levels = c('0','5','30','60','120','240','300','360','480','not applicable','FALSE'))
  df$ATTRIBUTE_Study_Day_Text <- factor(df$ATTRIBUTE_Study_Day_Text, levels = c('enrollment','pre-antibiotic PK','antibiotic Rx','post-antibiotic PK',
                                                                                '1 week post-antibiotic Rx','1 month post-antibiotic Rx'))
  # Normalization - none
  
  # PCA
  mypca <- pca(df[,-c(1:219)], method="nipals", center=TRUE, nPcs = 5, scale = "pareto")
    df <- cbind(scores(mypca),df)
    loadings <- as.data.frame(loadings(mypca))
    rownames <- rownames(loadings)
    loadings <- cbind(rownames,loadings)
    rownames(loadings) <- NULL
    head(loadings)
    test_loading_names <- loadings
    test_loading_names <- cbind(substr(test_loading_names$rownames, 0, 9),loadings)
    colnames(test_loading_names)[1] <- "mz"
    test_loading_names <- 
      separate(test_loading_names, col="mz", remove= TRUE, "_", into=c("feature_number","mz"))
    head(test_loading_names)
  
  # Define Variable to Plot
  plot_value <- df$ATTRIBUTE_Time_Point_Mins
  
  # Plot
  score_plot <- ggplot(df, aes(as.numeric(PC1), as.numeric(PC2)))+
    facet_wrap(~ATTRIBUTE_Study_Day_Text)+
    geom_point(aes(colour=as.factor(plot_value)), pch=16, cex=1.0, alpha=1.0) +
    theme_minimal()+
    theme(panel.border = element_rect(colour ="black", size=0.5, linetype = "solid", fill = NA),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(colour="black", size=6),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x=element_text(colour="black",size=6),
          aspect.ratio=1.0,
          legend.position="right",
          legend.title=element_text(colour="black", size=6),
          legend.text=element_text(colour="black", size=6) ) +
    labs(x="PC1", y="PC2")
  print(score_plot)
  
  # Plot
  density_score_plot <- ggplot(df, aes(as.numeric(PC1)))+
    facet_wrap(~ATTRIBUTE_Study_Day_Text)+
    geom_density(aes(colour=as.factor(plot_value),fill=as.factor(plot_value)), alpha=0.3, size=0.5) +
    theme_minimal()+
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(colour="black", size=6),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x=element_text(colour="black",size=6),
          aspect.ratio=1.0,
          legend.position="right",
          legend.title=element_text(colour="black", size=6),
          legend.text=element_text(colour="black", size=6) ) +
    labs(x="PC1 score values", y="Density")
  print(density_score_plot)

  rm(list=ls())
  
# Drug Levels Over Study Day or Time -------------------------------------

# INPUT: Data
df <- fread("20190710_DrugMetabolism_MS1DataMatrix_NoNorm_NoStds_wMetadata.csv", header=TRUE, sep=",")

# order study day
df$ATTRIBUTE_Study_Day <- factor(df$ATTRIBUTE_Study_Day, levels = c('0','1','2','3','4','5','6','7','8','9','15','38') )
df$ATTRIBUTE_Time_Point_Mins <- factor(df$ATTRIBUTE_Time_Point_Mins, levels = c('0','5','30','60','120','240','300','360','480','not applicable','FALSE'))
df$ATTRIBUTE_Study_Day_Text <- factor(df$ATTRIBUTE_Study_Day_Text, levels = c('enrollment','pre-antibiotic PK','antibiotic Rx','post-antibiotic PK',
                                                                              '1 week post-antibiotic Rx','1 month post-antibiotic Rx'))
# Identify Chemcials Associated with MS1 Features
  # Chemicals for Standards
      stercobilin <- df$X517_595.346390005393_3.12385043310876
      creatinine <- df$X283_114.062927408424_0.302783602150537
  # Internal standards
      IS_caffeine_d3 <- df$X328_198.107181109308_2.17063853080569
      IS_dextrorphan_d3 <- df$X354_261.203493607764_2.58340415695416
      IS_dextromethorphan_d3 <- df$X360_275.219843209508_3.02965592369478
      IS_omeprazole_d3 <- df$X807_349.140854149502_2.69188895131086
      IS_midazolam_d4 <- df$X372_330.111152477801_2.97450968325791 
  # Caffeine metabolites
      CAF_caffeine <- df$X972_195.089410390392_2.23677936743266
      CAF_methylxanthine <- df$X306_167.054350979357_0.513631927710843
      CAF_paraxanthine <- df$X316_181.072688941673_1.71319848484849
      CAF_methyluric <- df$X319_183.051179307203_0.570414741784038
      CAF_dimethyluric <- df$X327_197.068050226177_1.41345342247692
      CAF_theobromine <- df$X317_181.073012468502_1.10669328985507
      #CAF_AFMU <- not detected
  # Dextromethorphan metabolites
      DEX_dextromethorphan <- df$X1885_272.200753827816_3.11473963223787
      DEX_dextrorphan <- df$X836_258.184394650555_2.53843409465021
      DEX_dextrorphan_GLUC <- df$X379_434.217522717808_2.13214033101045
  # Midazolam metabolites
      MDZ_midazolam <- df$X5546_326.087324176636_3.01630319598521      
      MDZ_hydroxyMDZ <- df$X1877_342.082565409694_3.04860479148181
      MDZ_hydroxyMDZ_GLUC <- df$X858_518.113651968449_2.25715646135266
      MD_dihyroxyMDZ_GLUC <- df$X1880_534.108018180799_2.69321526294978
  # Omeprazole
      OME1_Omeprazole <- df$X5458_346.122664218162_2.7420736347359
      OME3_Omeprazole_sulfide <- df$X5112_330.126192609595_2.95065243759177
      OME4_Omeprazole_sulfone <- df$X5459_362.117961220971_3.03408651445966
      OME5_Hydroxyomeprazole <- df$X1046_362.118023279698_2.62215116245278
      OME7_Hydroxyomeprazole_sulfide <- df$X1052_346.122640331923_2.64336876520681
      OME11_Carboxyomeprazole <- df$X1047_376.097099204851_2.59171429324895
      OME13_Carboxyomeprazole_sulfide <- df$X854_360.100744014923_2.60551224007562
      OME1920_O_desmethylomeprazole_sulfide <- df$X850_316.112425327301_2.58288125
      OME29_Omeprazole_sulfide_4_O_glucuronide <- df$X1064_492.144319901645_2.3099519016864
      OME30_Omeprazole_5_O_glucuronide <- df$X1051_538.149987539353_2.57202220475192
      #OME9_Hydroxyomeprazole_sulfone <- not detected in FBMN
      #OME15_Carboxyomeprazole_sulfone <- not detected in FBMN
      #OME17_5_O_desmethylomeprazole <- not detected in FBMN
      #OME22_Didesmethylomeprazole_sulfide <- not detected in FBMN
      #OME27_Omeprazole_sulfide_O_glucuronide <- not detected in FBMN

# Create a matrix    
drugs <- cbind(df[,1:219],
               stercobilin,creatinine,IS_caffeine_d3,IS_dextrorphan_d3,IS_dextromethorphan_d3,
               IS_omeprazole_d3,IS_midazolam_d4,CAF_caffeine,CAF_methylxanthine,CAF_paraxanthine,
               CAF_methyluric,CAF_dimethyluric,CAF_theobromine,DEX_dextromethorphan,DEX_dextrorphan,
               DEX_dextrorphan_GLUC,MDZ_midazolam,MDZ_hydroxyMDZ,MDZ_hydroxyMDZ_GLUC,MD_dihyroxyMDZ_GLUC, OME1_Omeprazole,
               OME3_Omeprazole_sulfide,OME4_Omeprazole_sulfone,OME5_Hydroxyomeprazole,
               OME7_Hydroxyomeprazole_sulfide,OME11_Carboxyomeprazole,OME13_Carboxyomeprazole_sulfide,
               OME1920_O_desmethylomeprazole_sulfide, OME29_Omeprazole_sulfide_4_O_glucuronide,OME30_Omeprazole_5_O_glucuronide)

# remove include/exclude samples.
dim(drugs)
drugs <- subset(drugs,drugs$Include_Exclude == "include")
drugs <- subset(drugs,drugs$filename != "DM000088109_RD2_01_29084.mzXML")
dim(drugs)

# Assemble normalized matrix for all sample types
  # Fecal
    # subset drugs
    drugs_fecal <- subset(drugs, drugs$ATTRIBUTE_MS_SampleType_SampleSub == "fecal")
    # normalization by MS1 feature
    drugs_fecal_normal <- cbind(drugs_fecal[,1:219], sweep(drugs_fecal[,-c(1:219)], drugs_fecal$stercobilin, MARGIN=1, "/"))
    # wide to long format
    drugs_fecal_normal_long <- gather(drugs_fecal_normal, drug, peak_area, 220:ncol(drugs_fecal_normal))
    drugs_fecal_normal_long <- subset(drugs_fecal_normal_long, drugs_fecal_normal_long$drug != "stercobilin" & drugs_fecal_normal_long$drug != "creatinine")

  # Urine
    # subset drugs
    drugs_urine <- subset(drugs, drugs$ATTRIBUTE_MS_SampleType_SampleSub == "urine")
    # normalization by MS1 feature
    drugs_urine_normal <- cbind(drugs_urine[,1:219],sweep(drugs_urine[,-c(1:219)], drugs_urine$creatinine, MARGIN=1, "/"))
    # wide to long format
    drugs_urine_normal_long <- gather(drugs_urine_normal, drug, peak_area, 220:ncol(drugs))
    drugs_urine_long_plot <- subset(drugs_urine_normal_long, drugs_urine_normal_long$drug != "stercobilin" & drugs_urine_normal_long$drug != "creatinine")  

  # Plasma (Blood)
    # subset drugs
    drugs_blood <- subset(drugs, drugs$ATTRIBUTE_MS_SampleType_SampleSub == "blood_plasma")
    # normalization by MS1 feature - none
    # wide to long format
    drugs_blood_long <- gather(drugs_blood, drug, peak_area, 220:ncol(drugs))
    drugs_blood_long_plot <- subset(drugs_blood_long, drugs_blood_long$drug != "stercobilin" & drugs_blood_long$drug != "creatinine")  

# OUTPUT: Normalized Data for Comparison with Microbiome Data
out <- rbind(drugs_blood,drugs_urine_normal,drugs_fecal_normal)
write.table(out,"20190710_MSdata_drugs_normalized.tsv", sep="\t", row.names = FALSE)

# Plot - Drug level of time by subject
  # Plasma (blood)
    blood <- ggplot(data=drugs_blood_long_plot, aes(x=as.factor(ATTRIBUTE_Time_Point_Mins), y=as.numeric(peak_area)))+
      #facet_wrap(~drug, scales = "free")+
      facet_wrap(~drug, scales = "free")+
      geom_jitter(aes(group=drug), colour="black", width=0.15, cex=0.5, pch=20)+
      geom_boxplot(colour="black", fill="red", position = position_dodge(), 
                   alpha=0.25, width = 0.5, outlier.alpha = 0)+
      # optional scale y-axis
        #scale_y_log10()+
      theme_minimal()+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
            axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
            axis.text=element_text(colour="black",size=6),
            axis.text.x=element_text(angle=0, hjust=0.5, vjust=1.25),
            axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
            axis.title=element_text(colour="black",size=6),
            strip.text.x = element_text(colour="black",size=6),
            aspect.ratio=1,
            legend.position="bottom") +
      labs(x="Time (min)", y="Log10(Peak Area)")
    print(blood)
    ggsave("Curves_Blood.pdf", width = 10, height = 10, units = "in", dpi=300, useDingbats=FALSE)
  
  # Urine
  urine <- ggplot(data=drugs_urine_normal_long, aes(x=as.factor(ATTRIBUTE_Time_Point_Mins), y=as.numeric(peak_area)))+
    facet_wrap(~drug, scales = "free")+
    geom_jitter(aes(group=drug), colour="black", width=0.15, cex=0.5, pch=20)+
    geom_boxplot(colour="black", fill="orange", position = position_dodge(), 
                 alpha=0.25, width = 0.5, outlier.alpha = 0)+
    # optional scale y-axis
      #scale_y_log10()+
    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=1.25),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1,
          legend.position="bottom") +
    labs(x="Time (min)", y="Log10(Peak Area)")
  print(urine)
  ggsave("Curves_Urine.pdf", width = 10, height = 10, units = "in", dpi=300, useDingbats=FALSE)
  
  # Fecal
  fecal <- ggplot(data=drugs_fecal_normal_long, aes(x=as.factor(ATTRIBUTE_Study_Day), y=as.numeric(peak_area)))+
    facet_wrap(~drug, scales = "free")+
    geom_jitter(aes(group=drug), colour="black", width=0.15, cex=0.5, pch=20)+
    geom_boxplot(colour="black", fill="brown", position = position_dodge(), 
                 alpha=0.25, width = 0.5, outlier.alpha = 0)+
    # optional scale y-axis
      #scale_y_log10()+
    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=90, hjust=0.5, vjust=1.25),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1,
          legend.position="bottom") +
    labs(x="Time (min)", y="Peak Area")
  print(fecal)
  ggsave("Curves_Fecal.pdf", width = 10, height = 10, units = "in", dpi=300, useDingbats=FALSE)


# Drugs Levels Over Time or Day - Omeprazole Focus ------------------------
# Selection of Metabolite to Highlight in Manuscript Figure
drugs_blood_long_plot <- subset(drugs_blood_long_plot,
                                  drugs_blood_long_plot$drug == "OME1_Omeprazole" |
                                  drugs_blood_long_plot$drug == "OME5_Hydroxyomeprazole" |
                                  drugs_blood_long_plot$drug == "OME3_Omeprazole_sulfide" |
                                  drugs_blood_long_plot$drug == "OME7_Hydroxyomeprazole_sulfide")
drugs_blood_long_plot$drug <- factor(drugs_blood_long_plot$drug, 
                                     levels = c('OME1_Omeprazole','OME5_Hydroxyomeprazole',
                                                'OME3_Omeprazole_sulfide','OME7_Hydroxyomeprazole_sulfide'))
drugs_urine_normal_long <- subset(drugs_urine_normal_long,
                                  drugs_urine_normal_long$drug == "OME1_Omeprazole" |
                                  drugs_urine_normal_long$drug == "OME5_Hydroxyomeprazole" |
                                  drugs_urine_normal_long$drug == "OME3_Omeprazole_sulfide" |
                                  drugs_urine_normal_long$drug == "OME7_Hydroxyomeprazole_sulfide")
drugs_urine_normal_long$drug <- factor(drugs_urine_normal_long$drug, 
                                     levels = c('OME1_Omeprazole','OME5_Hydroxyomeprazole',
                                                'OME3_Omeprazole_sulfide','OME7_Hydroxyomeprazole_sulfide'))
drugs_fecal_normal_long <- subset(drugs_fecal_normal_long,
                                  drugs_fecal_normal_long$drug == "OME1_Omeprazole" |
                                    drugs_fecal_normal_long$drug == "OME5_Hydroxyomeprazole" |
                                    drugs_fecal_normal_long$drug == "OME3_Omeprazole_sulfide" |
                                    drugs_fecal_normal_long$drug == "OME7_Hydroxyomeprazole_sulfide")
drugs_fecal_normal_long$drug <- factor(drugs_fecal_normal_long$drug, 
                                       levels = c('OME1_Omeprazole','OME5_Hydroxyomeprazole',
                                                  'OME3_Omeprazole_sulfide','OME7_Hydroxyomeprazole_sulfide'))

# Plots
  # Plasma (blood)
  blood <- ggplot(data=drugs_blood_long_plot, aes(x=as.factor(ATTRIBUTE_Time_Point_Mins), y=as.numeric(peak_area)))+
    facet_wrap(~drug, scales = "free", ncol=4)+
    #geom_jitter(aes(group=drug), colour="black", width=0.15, cex=0.25, pch=20)+
    geom_boxplot(aes(colour=ATTRIBUTE_Study_Day_Text, fill=ATTRIBUTE_Study_Day_Text), alpha=0.3, width = 0.5, outlier.alpha = 0)+
    # optional scale y-axis
      scale_y_log10()+
    scale_color_manual(values = c("#2c7bb6", "#d7191c"))+
    scale_fill_manual(values = c("#2c7bb6", "#d7191c"))+
    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=90, hjust=0, vjust=0.25),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1.0,
          legend.position="none") +
    labs(x="Time (min)", y="Peak Area")
  print(blood)
  ggsave("OME_Curves_Blood.pdf", useDingbats=FALSE)
  
  # Urine
  urine <- ggplot(data=drugs_urine_normal_long, aes(x=as.factor(ATTRIBUTE_Time_Point_Mins), y=as.numeric(peak_area)))+
    facet_wrap(~drug, scales = "free", ncol=4)+
    #geom_jitter(aes(group=drug), colour="black", width=0.15, cex=0.25, pch=20)+
    geom_boxplot(aes(colour=ATTRIBUTE_Study_Day_Text, fill=ATTRIBUTE_Study_Day_Text),
                 alpha=0.3, width = 0.5, outlier.alpha = 0)+
    # optional scale y-axis
      scale_y_log10()+
    scale_color_manual(values = c("#2c7bb6", "#d7191c"))+
    scale_fill_manual(values = c("#2c7bb6", "#d7191c"))+
    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=90, hjust=0, vjust=0.25),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1,
          legend.position="bottom") +
    labs(x="Time (min)", y="Normalized Peak Area")
  print(urine)
  ggsave("OME_Curves_Urine.pdf", useDingbats=FALSE)
  
  # Feces 
  fecal <- ggplot(data=drugs_fecal_normal_long,
                  aes(x=as.factor(ATTRIBUTE_Study_Day), y=as.numeric(peak_area)))+
    facet_wrap(~drug, scales = "free", ncol=4)+
    #geom_jitter(aes(group=drug), colour="black", width=0.15, cex=0.25, pch=20)+
    geom_boxplot(colour="black", fill="grey50", position = position_dodge(), 
                 alpha=0.3, width = 0.5, outlier.alpha = 0)+
    # optional scale y-axis
      scale_y_log10()+
    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
          axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
          axis.text=element_text(colour="black",size=6),
          axis.text.x=element_text(angle=90, hjust=0, vjust=0.25),
          axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
          axis.title=element_text(colour="black",size=6),
          strip.text.x = element_text(colour="black",size=6),
          aspect.ratio=1,
          legend.position="bottom") +
    labs(x="Study Day", y="Normalized Peak Area")
  print(fecal)
  ggsave("OME_Curves_Fecal.pdf", useDingbats=FALSE)
  