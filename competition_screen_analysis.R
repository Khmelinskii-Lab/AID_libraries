# Load remotes package for installing specific package versions
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library(remotes)

# Install packages with specific versions used for the Manuscript
install_version("ggVennDiagram", version = "1.5.0")
install_version("ggpubr", version = "0.6.0")
install_version("skimr", version = "2.1.5")
install_version("readr", version = "2.1.4")
install_version("dplyr", version = "1.1.3")
install_version("ggpmisc", version = "0.5.4-1")
install_version("scales", version = "1.2.1")
install_version("lubridate", version = "1.9.3")
install_version("tidyr", version = "1.3.0")
install_version("ggplot2", version = "3.4.3")
install_version("ggpp", version = "0.5.4")
install_version("limma", version = "3.56.2")
install_version("forcats", version = "1.0.0")
install_version("tibble", version = "3.2.1")
install_version("ggrepel", version = "0.9.3")
install_version("lessR", version = "4.2.9")
install_version("stringr", version = "1.5.0")
install_version("tidyverse", version = "2.0.0")
install_version("usethis", version = "2.2.2")
install_version("moderndive", version = "0.5.5")
install_version("purrr", version = "1.0.2")
install_version("knitr", version = "1.44")

# Bioconductor package installation
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.17")  
BiocManager::install("limma", version = "3.56.2")

# Install from GitHub
remotes::install_github("kassambara/ggpubr")

# Load packages
library(ggplot2)
library(dplyr)
library(knitr)
library(tidyverse)
library(skimr)
library(moderndive)
library(stringr)
library(lessR)
library(limma)
library(scales)
library(ggpubr)
library(devtools)
library(ggrepel)
library(ggpmisc)
library(ggVennDiagram)

#Please change the directory path to where the code and files needed to run it are stored
path_input_dataset<-"~/Downloads/AID_libraries_1_competition_screen/"
setwd(path_input_dataset)

#--------------------------------------------------------------------------------------------------------------
#Fig1b bottom panel. Distribution of relative fitness effect in the genome-wide competition assay
#In this file the columns Winner_5, Winner_20 & Winner_50 correspond to the ORF classification based on the thresholds for 5%, 20% and 50% fitness difference
#In the columns Winner_5, Winner_20 & Winner_50; "N" stands for N-winners/C-losers, "C" stands for C-winners/N-losers & "neutral" stands for neutral ORFs
combined_plates_final_NoP59_t0_Winner_Corrected_doublings<-read.csv("combined_plates_final_NoP59_t0_Winner_Corrected_doublings.csv") 
Mean_Slope<-mean(combined_plates_final_NoP59_t0_Winner_Corrected_doublings$Slope_lm, na.rm=T)
SD_Slope<-sd(combined_plates_final_NoP59_t0_Winner_Corrected_doublings$Slope_lm, na.rm=T)

#Plot option 1
plot1b<-ggplot(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, aes(x=Slope_lm,fill =Winner_5)) +
  geom_histogram(binwidth=0.01, bins=1000) +
  scale_fill_manual(values=c("green","red","gray"), labels = c("C-winners/N-losers", "N-winners/C-losers", "neutrals"))+
  theme(axis.text.x = element_text(angle = 90))+
  geom_vline(xintercept=Mean_Slope, size=1.5, color="black", linetype="dotted")+
  scale_y_continuous(trans='log10')+
  xlab("Slope_lm")+ylab("log10 count")+
  ggtitle("lm slope distribution t0")+
  scale_x_continuous(breaks = seq(-2, 1.5, 0.5))+
  theme_classic()

cairo_pdf(filename = "Fig1b_bottom_panel.pdf", 
          width = 11, height = 8, onefile = TRUE)
plot1b
dev.off()

#Plot option 2 to visualize fitness effect with frequency of 1 ORF
hist_data = hist(combined_plates_final_NoP59_t0_Winner_Corrected_doublings$Slope_lm, breaks = 1000, plot = F)
hist_data$counts = log10(hist_data$counts) + 1
pdf("Fig1b_bottom_panel_v2.pdf", width = 11, height = 8)
plot(hist_data, ylim = c(-1,4))
dev.off()

#Count of ORFs by loser type for a given threshold column
#Function to summarize loser types 
count_loser_types <- function(data, threshold_column) {
  data %>%
    group_by(Loser_type = case_when(
      .data[[threshold_column]] == 'C' ~ 'N-losers',
      .data[[threshold_column]] == 'N' ~ 'C-losers',
      .data[[threshold_column]] == 'neutral' ~ 'Neutrals'
    ), .add = TRUE) %>%
    summarize(Loser_type_Count = n()) %>%
    mutate(Threshold = threshold_column)
}

#Applying the function for each threshold column and combining the results
Hits_table <- bind_rows(
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_5'),
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_20'),
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_50')
)

# Reordering columns
Hits_table <- Hits_table[, c(3, 1:2)]
print(Hits_table)
write.csv(Hits_table, paste0(path_input_dataset,"combined_plates_final_NoP59_t0_Winner_Corrected_doublings_Hits.csv"))

#Fig1b top panel. For visualizing the linear regression of all ORFs
#Examples shown in Fig1b can be searched by their ORF name: YNL031C, YER179W and YBR162W-A
all_days_competition<-read.csv("combined_plates_final_NoP59.csv")

all_days_competition<-subset(all_days_competition,!(all_days_competition$Slope_lm=='NA'))
write.csv(all_days_competition, paste0(path_input_dataset,"combined_plates_final_NoP59_noNA.csv"))

plate_well_vector<-unique((all_days_competition$Plate_Well))
for (i in 1:length(plate_well_vector)) {
  #i<-3  
  a<-subset(all_days_competition, all_days_competition$Plate_Well==plate_well_vector[i])
  ORF_name<-unique(a[,c("ORF")])
  Slope<-unique(round(a[,c("Slope_lm")],3))
  pdf(paste0(path_input_dataset,"Linear_Regression_Competition/Plate_Well_",plate_well_vector[i],"_ORF_",ORF_name,".pdf"))
  print(
    
    ggplot(a,aes(x=TimePoint,y = Log_Normalized_Ratio)) +
      geom_point(color='#2980B9', size = 4)+
      geom_abline(intercept=0, slope=lm(Log_Normalized_Ratio ~ 0+ TimePoint, data = a)$coefficients[1], color='#2C3E50', size=1.1) + 
      xlab("Time point")+ylab("Log_Normalized-Ratio C-pop/ N-pop ")+
      theme_bw()+
      ggtitle(paste0("lm regression intercept 0 of Plate-well : ",plate_well_vector[i]," ORF : ",ORF_name, "\nSlope : ",Slope))
    
  )
  dev.off()
  print(i)
}

#--------------------------------------------------------------------------------------------------------------
#Fig1c. ORFs affected by N' or C'terminal tagging as a function of a relative fitness threshold
combined_plates_final_NoP59_t0_Winner_Corrected_doublings<-read.csv("combined_plates_final_NoP59_t0_Winner_Corrected_doublings.csv") 

#Count of ORFs by loser type and essentiality for a given threshold column

#Function to summarize loser types
count_loser_types <- function(data, threshold_column, essentiality) {
  data %>%
    filter(ORF_essential == essentiality) %>%
    group_by(Loser_type = case_when(
      .data[[threshold_column]] == 'C' ~ 'N-losers',
      .data[[threshold_column]] == 'N' ~ 'C-losers',
      .data[[threshold_column]] == 'neutral' ~ 'Neutrals'
    ), .add = TRUE) %>%
    summarize(Loser_type_Count = n()) %>%
    mutate(Threshold = threshold_column, Essentiality = essentiality)
}

# Applying the function for each threshold column and essentiality and combining the results
Hits_essentiality_table <- bind_rows(
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_5', 'essential'),
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_20', 'essential'),
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_50', 'essential'),
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_5', 'non-essential'),
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_20', 'non-essential'),
  count_loser_types(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, 'Winner_50', 'non-essential')
)

# Reordering columns
Hits_essentiality_table <- Hits_essentiality_table[, c(1:4)]
print(Hits_essentiality_table)
write.csv(Hits_essentiality_table, paste0(path_input_dataset,"combined_plates_final_NoP59_t0_Winner_Corrected_doublings_Hits_Essentiality.csv"))

#Fig1c. Plot 
combined_plates_final_NoP59_t0_Winner_Corrected_doublings<-read.csv("combined_plates_final_NoP59_t0_Winner_Corrected_doublings.csv") 

#Essential ORFs
essential<-subset(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, combined_plates_final_NoP59_t0_Winner_Corrected_doublings$ORF_essential=="essential")
empty_dataframe<-as.data.frame(matrix(ncol=4))
min_slope<-min(combined_plates_final_NoP59_t0_Winner_Corrected_doublings$abs_Slope_lm)
max_slope<-max(combined_plates_final_NoP59_t0_Winner_Corrected_doublings$abs_Slope_lm)
interval<-seq(min_slope,max_slope,0.05)

for(j in 1:39){
  #j<-2 
  a<-c(as.numeric(interval[j]),
       length(unique(subset(essential,essential$abs_Slope_lm == interval[j])[,c("Plate_Well")])),
       length(unique(subset(essential,essential$abs_Slope_lm < interval[j])[,c("Plate_Well")])),
       length(unique(subset(essential,essential$abs_Slope_lm > interval[j])[,c("Plate_Well")])))
  empty_dataframe<-as.data.frame(rbind(empty_dataframe,a))
}

empty_dataframe<-empty_dataframe[2:nrow(empty_dataframe),]
names(empty_dataframe)<-c("Interval","Equal","less","more")
threshold_abs_slope_essential<-empty_dataframe
essential_ORFs<-as.numeric(length(essential$abs_Slope_lm))
threshold_abs_slope_essential$less_percentage<-(threshold_abs_slope_essential$less/essential_ORFs)*100
threshold_abs_slope_essential$more_percentage<-(threshold_abs_slope_essential$more/essential_ORFs)*100

#Non-essential ORFs
non_essential<-subset(combined_plates_final_NoP59_t0_Winner_Corrected_doublings, combined_plates_final_NoP59_t0_Winner_Corrected_doublings$ORF_essential=="non-essential")
empty_dataframe<-as.data.frame(matrix(ncol=4))
interval<-seq(min_slope,max_slope,0.05)

for(j in 1:39){
  #j<-2 
  a<-c(as.numeric(interval[j]),
       length(unique(subset(non_essential,non_essential$abs_Slope_lm == interval[j])[,c("Plate_Well")])),
       length(unique(subset(non_essential,non_essential$abs_Slope_lm < interval[j])[,c("Plate_Well")])),
       length(unique(subset(non_essential,non_essential$abs_Slope_lm > interval[j])[,c("Plate_Well")])))
  empty_dataframe<-as.data.frame(rbind(empty_dataframe,a))
}

empty_dataframe<-empty_dataframe[2:nrow(empty_dataframe),]
names(empty_dataframe)<-c("Interval","Equal","less","more")
threshold_abs_slope_non_essential<-empty_dataframe
non_essential_ORFs<-as.numeric(length(non_essential$abs_Slope_lm))
threshold_abs_slope_non_essential$less_percentage<-(threshold_abs_slope_non_essential$less/non_essential_ORFs)*100
threshold_abs_slope_non_essential$more_percentage<-(threshold_abs_slope_non_essential$more/non_essential_ORFs)*100

total_ORFs<-as.numeric(length(combined_plates_final_NoP59_t0_Winner_Corrected_doublings$abs_Slope_lm))
empty_dataframe<-as.data.frame(matrix(ncol=4))
min_slope<-min(combined_plates_final_NoP59_t0_Winner_Corrected_doublings$abs_Slope_lm)
max_slope<-max(combined_plates_final_NoP59_t0_Winner_Corrected_doublings$abs_Slope_lm)
interval<-seq(min_slope,max_slope,0.05)

for(j in 1:39){
  #j<-2 
  a<-c(as.numeric(interval[j]),
       length(unique(subset(combined_plates_final_NoP59_t0_Winner_Corrected_doublings,combined_plates_final_NoP59_t0_Winner_Corrected_doublings$abs_Slope_lm == interval[j])[,c("Plate_Well")])),
       length(unique(subset(combined_plates_final_NoP59_t0_Winner_Corrected_doublings,combined_plates_final_NoP59_t0_Winner_Corrected_doublings$abs_Slope_lm < interval[j])[,c("Plate_Well")])),
       length(unique(subset(combined_plates_final_NoP59_t0_Winner_Corrected_doublings,combined_plates_final_NoP59_t0_Winner_Corrected_doublings$abs_Slope_lm > interval[j])[,c("Plate_Well")])))
  empty_dataframe<-as.data.frame(rbind(empty_dataframe,a))
}

empty_dataframe<-empty_dataframe[2:nrow(empty_dataframe),]
names(empty_dataframe)<-c("Interval","Equal","less","more")
threshold_abs_slope<-empty_dataframe
threshold_abs_slope$less_percentage<-(threshold_abs_slope$less/total_ORFs)*100
threshold_abs_slope$more_percentage<-(threshold_abs_slope$more/total_ORFs)*100

threshold_abs_slope$Equal_AllORFs<-threshold_abs_slope$Equal
threshold_abs_slope$less_AllORFs<-threshold_abs_slope$less
threshold_abs_slope$more_AllORFs<-threshold_abs_slope$more
threshold_abs_slope$less_percentage_AllORFs<-threshold_abs_slope$less_percentage
threshold_abs_slope$more_percentage_AllORFs<-threshold_abs_slope$more_percentage
threshold_abs_slope$Equal<-NULL
threshold_abs_slope$less<-NULL
threshold_abs_slope$more<-NULL
threshold_abs_slope$less_percentage<-NULL
threshold_abs_slope$more_percentage<-NULL

threshold_abs_slope_essential$Equal_essentialORFs<-threshold_abs_slope_essential$Equal
threshold_abs_slope_essential$less_essentialORFs<-threshold_abs_slope_essential$less
threshold_abs_slope_essential$more_essentialORFs<-threshold_abs_slope_essential$more
threshold_abs_slope_essential$less_percentage_essentialORFs<-threshold_abs_slope_essential$less_percentage
threshold_abs_slope_essential$more_percentage_essentialORFs<-threshold_abs_slope_essential$more_percentage
threshold_abs_slope_essential$Equal<-NULL
threshold_abs_slope_essential$less<-NULL
threshold_abs_slope_essential$more<-NULL
threshold_abs_slope_essential$less_percentage<-NULL
threshold_abs_slope_essential$more_percentage<-NULL

threshold_abs_slope_non_essential$Equal_nonessentialORFs<-threshold_abs_slope_non_essential$Equal
threshold_abs_slope_non_essential$less_nonessentialORFs<-threshold_abs_slope_non_essential$less
threshold_abs_slope_non_essential$more_nonessentialORFs<-threshold_abs_slope_non_essential$more
threshold_abs_slope_non_essential$less_percentage_nonessentialORFs<-threshold_abs_slope_non_essential$less_percentage
threshold_abs_slope_non_essential$more_percentage_nonessentialORFs<-threshold_abs_slope_non_essential$more_percentage
threshold_abs_slope_non_essential$Equal<-NULL
threshold_abs_slope_non_essential$less<-NULL
threshold_abs_slope_non_essential$more<-NULL
threshold_abs_slope_non_essential$less_percentage<-NULL
threshold_abs_slope_non_essential$more_percentage<-NULL

threshold_final<-left_join(threshold_abs_slope, threshold_abs_slope_essential, by = c("Interval" = "Interval")) 
threshold_final<-left_join(threshold_final, threshold_abs_slope_non_essential, by = c("Interval" = "Interval")) 

Fig1c<-ggplot(threshold_final, aes(x=Interval)) + 
  geom_line(aes(y = less_percentage_AllORFs), color = "black", size=1.5) + 
  geom_line(aes(y = more_percentage_AllORFs), color="black", size=1.5)+
  geom_line(aes(y = less_percentage_essentialORFs), color = "blue", size=1.5) + 
  geom_line(aes(y = more_percentage_essentialORFs), color="blue", size=1.5)+
  geom_line(aes(y = less_percentage_nonessentialORFs), color = "purple",size=1.5) + 
  geom_line(aes(y = more_percentage_nonessentialORFs), color="purple",  size=1.5)+
  theme(axis.text.x = element_text(angle = 90))+
  theme_classic()+
  #scale_y_continuous(trans='log10')+
  scale_y_log10()+
  xlab("Interval of Abs slope")+ylab("log10 ORF-Hits percentage")+
  ggtitle("Threshold of Abs slopes for hits")

cairo_pdf(filename = "Fig1c_Plot.pdf", 
          width = 11, height = 8, onefile = TRUE)
Fig1c
dev.off()

#_______________________________________________________________________________________________
#Fig1d. Subcellular localization analysis using Hu et al, 2003 (Global analysis of protein localization in budding yeast)
combined_plates_final_signals_SWAT_SP_MTS<-read.csv("combined_plates_final_signals_SWAT_SP_MTS.csv")
Localization_Huh <- read.csv("20220912_GFPLocalizationLibrary 1.csv")

combined_plates_final_signals_SWAT_SP_MTS_subcellular<-left_join(combined_plates_final_signals_SWAT_SP_MTS, Localization_Huh, by = c("ORF" = "ORF")) 
#write.csv(combined_plates_final_signals_SWAT_SP_MTS_subcellular, paste0(path_input_dataset,"combined_plates_final_signals_SWAT_SP_MTS_subcellular.csv"))
#combined_plates_final_signals_SWAT_SP_MTS_subcellular<-read.csv("combined_plates_final_signals_SWAT_SP_MTS_subcellular.csv")

#Filter out rows where "GFP.visualized" is "not visualized"
subcellular_visualized<-subset(combined_plates_final_signals_SWAT_SP_MTS_subcellular,!(GFP.visualized.=="not visualized")) #only GFP visualized proteins to be analized
#Subset the dataset into "non-essential" and "essential" based on the "ORF_essential" column
subcellular_visualized_Essential <- subset(subcellular_visualized, ORF_essential == "essential")
subcellular_visualized_NonEssential <- subset(subcellular_visualized, ORF_essential == "non-essential")

#5% Essential & Non-Essentials N, C & neutrals 
# Function to count occurrences of a certain value in a column based on TRUE/FALSE entries in other columns
count_entries <- function(data, target_value, count_column, logical_columns) {
  true_count <- numeric(length(logical_columns))
  false_count <- numeric(length(logical_columns))
  
  for (i in seq_along(logical_columns)) {
    true_count[i] <- sum(data[[count_column]] == target_value & data[[logical_columns[i]]] == "TRUE")
    false_count[i] <- sum(data[[count_column]] == target_value & data[[logical_columns[i]]] == "FALSE")
  }
  
  return(list("TRUE" = true_count, "FALSE" = false_count))
}

# List of columns to iterate over for counting
columns_to_count <- c("Ambiguous", "Mitochondrion", "Vacuole", "Spindle.pole", "Cell.periphery", 
                      "Punctate.composite", "Vacuolar.membrane", "ER", "Nuclear.periphery", 
                      "Endosome", "Bud.neck", "Microtubule", "Golgi", "Late.Golgi", "Peroxisome", 
                      "Actin", "Nucleolus", "Cytoplasm", "ER.to.Golgi", "Early.Golgi", 
                      "Lipid.particle", "Nucleus", "Bud")

# Function to perform counting for both essential and non-essential subsets
# To asses fitness thresholds 20% or 50%, change the "Winner_5" in the function for "Winner_20" or "Winner_50"
perform_counting <- function(data, target_value) {
  count_results <- list()
  
  for (col in columns_to_count) {
    counts <- count_entries(data, target_value, "Winner_5", col)
    count_results[[col]] <- counts
  }
  
  return(count_results)
}

# Count occurrences of "N", "neutral" and "C" for essential subset
essential_N_counts <- perform_counting(subcellular_visualized_Essential, "N")
essential_neutral_counts <- perform_counting(subcellular_visualized_Essential, "neutral")
essential_C_counts <- perform_counting(subcellular_visualized_Essential, "C")

# Count occurrences of "N", "neutral" and "C" for non-essential subset
non_essential_N_counts <- perform_counting(subcellular_visualized_NonEssential, "N")
non_essential_neutral_counts <- perform_counting(subcellular_visualized_NonEssential, "neutral")
non_essential_C_counts <- perform_counting(subcellular_visualized_NonEssential, "C")

# Create empty dataframe to store output
output_table <- data.frame(
  Protein_Essentiality = character(),
  Winner_5 = character(),
  Protein_Essentiality_Count = numeric(),
  Subcellular_Localization = character(),
  Subcellular_TRUE_Count = numeric(),
  Subcellular_FALSE_Count = numeric(),
  stringsAsFactors = FALSE
)

# Function to add counts to output table
add_counts_to_output <- function(output_df, counts, essentiality, target_value) {
  for (category in names(counts)) {
    true_counts <- counts[[category]]$"TRUE"
    false_counts <- counts[[category]]$"FALSE"
    
    output_df <- rbind(output_df, data.frame(
      Protein_Essentiality = rep(essentiality, length(true_counts)),
      Winner_5 = rep(target_value, length(true_counts)),
      Protein_Essentiality_Count = true_counts + false_counts,
      Subcellular_Localization = rep(category, length(true_counts)),
      Subcellular_TRUE_Count = true_counts,
      Subcellular_FALSE_Count = false_counts,
      stringsAsFactors = FALSE
    ))
  }
  
  return(output_df)
}

# Add counts for essential subset
output_table <- add_counts_to_output(output_table, essential_N_counts, "essential", "N")
output_table <- add_counts_to_output(output_table, essential_neutral_counts, "essential", "neutral")
output_table <- add_counts_to_output(output_table, essential_C_counts, "essential", "C")

# Add counts for non-essential subset
output_table <- add_counts_to_output(output_table, non_essential_N_counts, "non-essential", "N")
output_table <- add_counts_to_output(output_table, non_essential_neutral_counts, "non-essential", "neutral")
output_table <- add_counts_to_output(output_table, non_essential_C_counts, "non-essential", "C")
write.csv(output_table, paste0(path_input_dataset,"Fig1d_5.csv"))

# Read the CSV file
Fig1d_5 <- read.csv("Fig1d_5.csv")

# Calculate Subcellular_TRUE_Percentage
unique_combinations <- unique(paste(Fig1d_5$Subcellular_Localization, Fig1d_5$Protein_Essentiality))

for (combination in unique_combinations) {
  split_comb <- strsplit(combination, " ")[[1]]
  loc <- split_comb[1]
  ess <- split_comb[2]
  subset_table <- subset(Fig1d_5, Subcellular_Localization == loc & Protein_Essentiality == ess)
  total_true_counts <- sum(subset_table$Subcellular_TRUE_Count)
  subset_table$Subcellular_TRUE_Percentage <- with(subset_table, (Subcellular_TRUE_Count / total_true_counts) * 100)
  Fig1d_5[Fig1d_5$Subcellular_Localization == loc & Fig1d_5$Protein_Essentiality == ess, "Subcellular_TRUE_Percentage"] <- subset_table$Subcellular_TRUE_Percentage
}

Fig1d_5
write.csv(Fig1d_5, paste0(path_input_dataset,"Fig1d_5_PlotTable.csv"))
Fig1d_5<- read.csv("Fig1d_5_PlotTable.csv")

subcellular_5<-ggplot(Fig1d_5, aes(x = Subcellular_TRUE_Percentage, y =Subcellular_Localization, fill = Winner_5)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Protein_Essentiality, scales = "free") +
  labs(title = "% ORFs with specific subcellular localization 5% threshold",
       x = "% ORFs", y = "Protein localization") +
  scale_fill_manual(values = c("magenta", "green","gray")) +
  theme_minimal()

cairo_pdf(filename = "Fig1d_5_Plot.pdf", 
          width = 11, height = 8, onefile = TRUE)
subcellular_5
dev.off()

#_______________________________________________________________________________________________
#Fig1e. Subcellular localization analysis using Hu et al, 2003 (Global analysis of protein localization in budding yeast)
#Proteins with differential localization of N- and C-terminally GFP tagged variants (Huh et al., 2003; Weill et al., 2018), stratified by gene essentiality and differential fitness according to b. N-loser, neutral and C-loser ORFs defined at two differential fitness thresholds. Number of ORFs in each group is indicated.

#5% 
# Function to subset data based on specified conditions
combined_plates_final_signals_SWAT_SP_MTS_localization<-read.csv("combined_plates_final_signals_SWAT_SP_MTS_localization.csv")#this is the newest complete file

generate_summary_table <- function(data, essentiality) {
  # Subset the data based on essentiality
  if (essentiality == "Non-Essential") {
    subset_data <- subset(data, !(ORF_essential == "essential"))
  } else if (essentiality == "Essential") {
    subset_data <- subset(data, !(ORF_essential == "non-essential"))
  }
  
  # Create empty vectors to store counts
  counts <- c()
  same_counts <- c()
  diff_counts <- c()
  nnop1pr_counts <- c()
  nnop1pr_same_counts <- c()
  nnop1pr_diff_counts <- c()
  
  # Group by the 'Winner_5' column
  groups <- unique(subset_data$Winner_5)
  
  # Iterate over groups
  for (group in groups) {
    # Filter data for the current group
    group_data <- subset(subset_data, Winner_5 == group)
    
    # Count N vs C endogenous localization
    n_vs_c_endogenous <- table(group_data$N_C_localization_endogenous_v2)
    n_vs_c_endogenous_same <- n_vs_c_endogenous["Same"]
    n_vs_c_endogenous_diff <- n_vs_c_endogenous["Different"]
    
    # Count NNop1pr vs C endogenous localization
    nnop1pr_vs_c_endogenous <- table(group_data$N_C_localization_NNop1pr_C_v2)
    nnop1pr_vs_c_endogenous_same <- nnop1pr_vs_c_endogenous["Same"]
    nnop1pr_vs_c_endogenous_diff <- nnop1pr_vs_c_endogenous["Different"]
    
    # Count total entries for NNop1pr vs C endogenous localization
    total_nnop1pr <- sum(nnop1pr_vs_c_endogenous)
    
    # Count total entries for the current group
    total_entries <- nrow(group_data)
    
    # Append counts to vectors
    counts <- c(counts, total_entries)
    same_counts <- c(same_counts, n_vs_c_endogenous_same)
    diff_counts <- c(diff_counts, n_vs_c_endogenous_diff)
    nnop1pr_counts <- c(nnop1pr_counts, total_nnop1pr)
    nnop1pr_same_counts <- c(nnop1pr_same_counts, nnop1pr_vs_c_endogenous_same)
    nnop1pr_diff_counts <- c(nnop1pr_diff_counts, nnop1pr_vs_c_endogenous_diff)
  }
  
  # Combine counts into a data frame
  summary_table <- data.frame(Data = groups,
                              Protein_Essentiality = essentiality,
                              NvsC_endogenous_Localization_Count = counts,
                              NvsC_endogenous_Localization_Same = same_counts,
                              NvsC_endogenous_Localization_Different = diff_counts,
                              NNop1prvsC_Localization_Count = nnop1pr_counts,
                              NNop1prvsC_Localization_Same = nnop1pr_same_counts,
                              NNop1prvsC_Localization_Different = nnop1pr_diff_counts)
  
  return(summary_table)
}

# Define a list to store the summary tables for each subset
summary_tables <- list()

# Iterate over essentiality levels
essentialities <- c("Non-Essential", "Essential")
for (ess in essentialities) {
  summary_tables[[length(summary_tables) + 1]] <- generate_summary_table(combined_plates_final_signals_SWAT_SP_MTS_localization, ess)
}

# Combine summary tables into one
All_5_table <- do.call(rbind, summary_tables)
All_5_table$Data_type<-'All'
All_5_table$ID<-paste(All_5_table$Data_type,All_5_table$Data,All_5_table$Protein_Essentiality,sep="_")

#Endogenous
generate_summary_table <- function(data) {
  # Define function to subset data based on missing, below threshold, and ambiguous values
  subset_data <- function(data, column, value) {
    data <- subset(data, !(column == value))
  }
  
  # Subsetting the data 
  data <- subset_data(data, data$N_NATIVEpr.GFP.localization, "missing")
  data <- subset_data(data, data$C_Localization.combined, "missing")
  data <- subset_data(data, data$N_NATIVEpr.GFP.localization, "below threshold")
  data <- subset_data(data, data$C_Localization.combined, "below threshold")
  data <- subset_data(data, data$N_NATIVEpr.GFP.localization, "ambiguous")
  data <- subset_data(data, data$C_Localization.combined, "ambiguous")
  
  # Group by the 'Winner_5' column
  groups <- unique(data$Winner_5)
  
  # Create empty vectors to store counts
  counts <- c()
  same_counts <- c()
  diff_counts <- c()
  group_info <- c()
  essentiality_info <- c()
  
  # Iterate over groups
  for (group in groups) {
    # Filter data for the current group
    group_data <- subset(data, Winner_5 == group)
    
    # Create subsets for Non-Essential and Essential within the current group
    subset_non_essential <- subset(group_data, !(ORF_essential == "essential"))
    subset_essential <- subset(group_data, !(ORF_essential == "non-essential"))
    
    # Count occurrences of "Same" and "Different" for Non-Essential subset within the current group
    n_same_non_essential <- nrow(subset_non_essential[subset_non_essential$N_C_localization_endogenous_v2 == 'Same', ])
    n_different_non_essential <- nrow(subset_non_essential[subset_non_essential$N_C_localization_endogenous_v2 == 'Different', ])
    
    # Count occurrences of "Same" and "Different" for Essential subset within the current group
    n_same_essential <- nrow(subset_essential[subset_essential$N_C_localization_endogenous_v2 == 'Same', ])
    n_different_essential <- nrow(subset_essential[subset_essential$N_C_localization_endogenous_v2 == 'Different', ])
    
    # Append counts to vectors
    counts <- c(counts, n_same_non_essential + n_different_non_essential, n_same_essential + n_different_essential)
    same_counts <- c(same_counts, n_same_non_essential, n_same_essential)
    diff_counts <- c(diff_counts, n_different_non_essential, n_different_essential)
    group_info <- c(group_info, rep(group, 2))
    essentiality_info <- c(essentiality_info, "Non-Essential", "Essential")
  }
  
  # Combine counts into a data frame
  summary_table <- data.frame(
    Data = group_info,
    Protein_Essentiality = essentiality_info,
    NvsC_endogenous_Localization_Count = counts,
    NvsC_endogenous_Localization_Same = same_counts,
    NvsC_endogenous_Localization_Different = diff_counts
  )
  
  # Return summary table
  return(summary_table)
}

endogenous_table <- generate_summary_table(combined_plates_final_signals_SWAT_SP_MTS_localization)
endogenous_table$Data_type<-'Subset'
endogenous_table$ID<-paste(endogenous_table$Data_type,endogenous_table$Data,endogenous_table$Protein_Essentiality,sep="_")

#non endogenous
generate_summary_table <- function(data) {
  # Define function to subset data based on missing, below threshold, and ambiguous values
  subset_data <- function(data, column, value) {
    data <- subset(data, !(column == value))
  }
  
  # Subsetting the data 
  data <- subset_data(data, data$NOP1pr.GFP.localization, "missing")
  data <- subset_data(data, data$C_Localization.combined, "missing")
  data <- subset_data(data, data$NOP1pr.GFP.localization, "below threshold")
  data <- subset_data(data, data$C_Localization.combined, "below threshold")
  data <- subset_data(data, data$NOP1pr.GFP.localization, "ambiguous")
  data <- subset_data(data, data$C_Localization.combined, "ambiguous")
  
  # Group by the 'Winner_5' column
  groups <- unique(data$Winner_5)
  
  # Create empty vectors to store counts
  counts <- c()
  same_counts <- c()
  diff_counts <- c()
  group_info <- c()
  essentiality_info <- c()
  
  # Iterate over groups
  for (group in groups) {
    # Filter data for the current group
    group_data <- subset(data, Winner_5 == group)
    
    # Create subsets for Non-Essential and Essential within the current group
    subset_non_essential <- subset(group_data, !(ORF_essential == "essential"))
    subset_essential <- subset(group_data, !(ORF_essential == "non-essential"))
    
    # Count occurrences of "Same" and "Different" for Non-Essential subset within the current group
    n_same_non_essential <- nrow(subset_non_essential[subset_non_essential$N_C_localization_NNop1pr_C_v2 == 'Same', ])
    n_different_non_essential <- nrow(subset_non_essential[subset_non_essential$N_C_localization_NNop1pr_C_v2 == 'Different', ])
    
    # Count occurrences of "Same" and "Different" for Essential subset within the current group
    n_same_essential <- nrow(subset_essential[subset_essential$N_C_localization_NNop1pr_C_v2 == 'Same', ])
    n_different_essential <- nrow(subset_essential[subset_essential$N_C_localization_NNop1pr_C_v2 == 'Different', ])
    
    # Append counts to vectors
    counts <- c(counts, n_same_non_essential + n_different_non_essential, n_same_essential + n_different_essential)
    same_counts <- c(same_counts, n_same_non_essential, n_same_essential)
    diff_counts <- c(diff_counts, n_different_non_essential, n_different_essential)
    group_info <- c(group_info, rep(group, 2))
    essentiality_info <- c(essentiality_info, "Non-Essential", "Essential")
  }
  
  # Combine counts into a data frame
  summary_table <- data.frame(
    Data = group_info,
    Protein_Essentiality = essentiality_info,
    NNop1prvsC_Localization_Count = counts,
    NNop1prvsC_Localization_Same = same_counts,
    NNop1prvsC_Localization_Different = diff_counts
  )
  
  # Return summary table
  return(summary_table)
}

nonendogenous_table <- generate_summary_table(combined_plates_final_signals_SWAT_SP_MTS_localization)
colnames(nonendogenous_table) <- c("Data","Protein_Essentiality","NNop1prvsC_Localization_Count","NNop1prvsC_Localization_Same","NNop1prvsC_Localization_Different")
nonendogenous_table$Data_type<-'Subset'
nonendogenous_table$ID<-paste(nonendogenous_table$Data_type,nonendogenous_table$Data,nonendogenous_table$Protein_Essentiality,sep="_")

# Combine all dataframes
subset_5_table <- left_join(endogenous_table, nonendogenous_table, by="ID")
subset_5_table$Data.y<-NULL
subset_5_table$Data_type.y<-NULL
subset_5_table$Protein_Essentiality.y<-NULL
subset_5_table$Data_type<-subset_5_table$Data_type.x
subset_5_table$Data_type.x<-NULL
subset_5_table$ID2<-subset_5_table$ID
subset_5_table$ID<-NULL
subset_5_table$ID<-subset_5_table$ID2
subset_5_table$ID2<-NULL
names(subset_5_table)[1:2] <- c("Data","Protein_Essentiality")
threshold5_table<- rbind(All_5_table,subset_5_table)
threshold5_table<- threshold5_table[, c("Data_type", "Data","Protein_Essentiality","ID","NvsC_endogenous_Localization_Count","NvsC_endogenous_Localization_Different","NvsC_endogenous_Localization_Same","NNop1prvsC_Localization_Count","NNop1prvsC_Localization_Different","NNop1prvsC_Localization_Same")]
threshold5_table$Data[threshold5_table$Data=="N"]<-'C-loser'
threshold5_table$Data[threshold5_table$Data=="C"]<-'N-loser'
#write.csv(threshold5_table, paste0(path_input_dataset,"Fig1e_5.csv"))
#threshold5_table<-read.csv("Fig1e_5.csv")

# Calculate percentages and add new columns
threshold5_table$NvsC_endogenous_Localization_Different_Percentage <- (threshold5_table$NvsC_endogenous_Localization_Different / threshold5_table$NvsC_endogenous_Localization_Count) * 100
threshold5_table$NvsC_endogenous_Localization_Same_Percentage <- (threshold5_table$NvsC_endogenous_Localization_Same / threshold5_table$NvsC_endogenous_Localization_Count) * 100
threshold5_table$NNop1prvsC_Localization_Different_Percentage <- (threshold5_table$NNop1prvsC_Localization_Different / threshold5_table$NNop1prvsC_Localization_Count) * 100
threshold5_table$NNop1prvsC_Localization_Same_Percentage <- (threshold5_table$NNop1prvsC_Localization_Same / threshold5_table$NNop1prvsC_Localization_Count) * 100
threshold5_table$X<-NULL
threshold5_table$Category<-NULL
threshold5_table$ID<-NULL
colnames(threshold5_table)[2]<-"Winner_5"
threshold5_table$NvsC_endogenous_Localization_Count<-NULL
threshold5_table$NvsC_endogenous_Localization_Different<-NULL
threshold5_table$NvsC_endogenous_Localization_Same<-NULL
threshold5_table$NNop1prvsC_Localization_Count<-NULL
threshold5_table$NNop1prvsC_Localization_Different<-NULL
threshold5_table$NNop1prvsC_Localization_Same<-NULL
write.csv(threshold5_table, paste0(path_input_dataset,"Fig1e_5_PlotTable.csv"))
Fig1e_5_PlotTable<-read.csv("Fig1e_5_PlotTable.csv")
Fig1e_5_PlotTable$X<-NULL

# Reshape the data
Fig1e_5_PlotTable <- pivot_longer(Fig1e_5_PlotTable, 
                                  cols = starts_with("NvsC") | starts_with("NNop1prvsC"), 
                                  names_to = "Category", 
                                  values_to = "Percentage")

# Reorder the columns
Fig1e_5_PlotTable <- Fig1e_5_PlotTable[, c("Data_type", "Protein_Essentiality", "Winner_5", "Category", "Percentage")]

# Write the formatted data to a new CSV file
write.csv(Fig1e_5_PlotTable, "Fig1e_5_PlotTable_LongFormat.csv", row.names = FALSE)

#PLOT for Endogenous Subset data
Fig1e_5_PlotTable_Subset<-subset(Fig1e_5_PlotTable,!(Fig1e_5_PlotTable$Data_type=="All"))
Fig1e_5_PlotTable_Subset_Endogenous<-subset(Fig1e_5_PlotTable_Subset,!(Fig1e_5_PlotTable_Subset$Category=="NNop1prvsC_Localization_Different_Percentage"))
Fig1e_5_PlotTable_Subset_Endogenous<-subset(Fig1e_5_PlotTable_Subset_Endogenous,!(Fig1e_5_PlotTable_Subset_Endogenous$Category=="NNop1prvsC_Localization_Same_Percentage"))

e_5<-ggplot(Fig1e_5_PlotTable_Subset_Endogenous, aes(x = Percentage, y = Winner_5, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Protein_Essentiality, scales = "free") +
  labs(title = "% ORFs tagged Endogenously with different Localization 5% threshold",
       x = "% ORFs", y = "Loser type") +
  scale_fill_manual(values = c("red", "gray")) +
  theme_minimal()

#PLOT for Non-Endogenous Subset data
Fig1e_5_PlotTable_Subset<-subset(Fig1e_5_PlotTable,!(Fig1e_5_PlotTable$Data_type=="All"))
Fig1e_5_PlotTable_Subset_NonEndogenous<-subset(Fig1e_5_PlotTable_Subset,!(Fig1e_5_PlotTable_Subset$Category=="NvsC_endogenous_Localization_Different_Percentage"))
Fig1e_5_PlotTable_Subset_NonEndogenous<-subset(Fig1e_5_PlotTable_Subset_NonEndogenous,!(Fig1e_5_PlotTable_Subset_NonEndogenous$Category=="NvsC_endogenous_Localization_Same_Percentage"))

ne_5<-ggplot(Fig1e_5_PlotTable_Subset_NonEndogenous, aes(x = Percentage, y = Winner_5, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Protein_Essentiality, scales = "free") +
  labs(title = "% ORFs tagged Non-Endogenously with different Localization 5% threshold",
       x = "% ORFs", y = "Loser type") +
  scale_fill_manual(values = c("red", "gray")) +
  theme_minimal()

cairo_pdf(filename = "Fig1e_5_Plots.pdf", 
          width = 11, height = 8, onefile = TRUE)
e_5
ne_5
dev.off()

#20% 
# Function to subset data based on specified conditions
combined_plates_final_signals_SWAT_SP_MTS_localization<-read.csv("combined_plates_final_signals_SWAT_SP_MTS_localization.csv")#this is the newest complete file

generate_summary_table <- function(data, essentiality) {
  # Subset the data based on essentiality
  if (essentiality == "Non-Essential") {
    subset_data <- subset(data, !(ORF_essential == "essential"))
  } else if (essentiality == "Essential") {
    subset_data <- subset(data, !(ORF_essential == "non-essential"))
  }
  
  # Create empty vectors to store counts
  counts <- c()
  same_counts <- c()
  diff_counts <- c()
  nnop1pr_counts <- c()
  nnop1pr_same_counts <- c()
  nnop1pr_diff_counts <- c()
  
  # Group by the 'Winner_20' column
  groups <- unique(subset_data$Winner_20)
  
  # Iterate over groups
  for (group in groups) {
    # Filter data for the current group
    group_data <- subset(subset_data, Winner_20 == group)
    
    # Count N vs C endogenous localization
    n_vs_c_endogenous <- table(group_data$N_C_localization_endogenous_v2)
    n_vs_c_endogenous_same <- n_vs_c_endogenous["Same"]
    n_vs_c_endogenous_diff <- n_vs_c_endogenous["Different"]
    
    # Count NNop1pr vs C endogenous localization
    nnop1pr_vs_c_endogenous <- table(group_data$N_C_localization_NNop1pr_C_v2)
    nnop1pr_vs_c_endogenous_same <- nnop1pr_vs_c_endogenous["Same"]
    nnop1pr_vs_c_endogenous_diff <- nnop1pr_vs_c_endogenous["Different"]
    
    # Count total entries for NNop1pr vs C endogenous localization
    total_nnop1pr <- sum(nnop1pr_vs_c_endogenous)
    
    # Count total entries for the current group
    total_entries <- nrow(group_data)
    
    # Append counts to vectors
    counts <- c(counts, total_entries)
    same_counts <- c(same_counts, n_vs_c_endogenous_same)
    diff_counts <- c(diff_counts, n_vs_c_endogenous_diff)
    nnop1pr_counts <- c(nnop1pr_counts, total_nnop1pr)
    nnop1pr_same_counts <- c(nnop1pr_same_counts, nnop1pr_vs_c_endogenous_same)
    nnop1pr_diff_counts <- c(nnop1pr_diff_counts, nnop1pr_vs_c_endogenous_diff)
  }
  
  # Combine counts into a data frame
  summary_table <- data.frame(Data = groups,
                              Protein_Essentiality = essentiality,
                              NvsC_endogenous_Localization_Count = counts,
                              NvsC_endogenous_Localization_Same = same_counts,
                              NvsC_endogenous_Localization_Different = diff_counts,
                              NNop1prvsC_Localization_Count = nnop1pr_counts,
                              NNop1prvsC_Localization_Same = nnop1pr_same_counts,
                              NNop1prvsC_Localization_Different = nnop1pr_diff_counts)
  
  return(summary_table)
}

# Define a list to store the summary tables for each subset
summary_tables <- list()

# Iterate over essentiality levels
essentialities <- c("Non-Essential", "Essential")
for (ess in essentialities) {
  summary_tables[[length(summary_tables) + 1]] <- generate_summary_table(combined_plates_final_signals_SWAT_SP_MTS_localization, ess)
}

# Combine summary tables into one
All_20_table <- do.call(rbind, summary_tables)
All_20_table$Data_type<-'All'
All_20_table$ID<-paste(All_20_table$Data_type,All_20_table$Data,All_20_table$Protein_Essentiality,sep="_")

#Endogenous
generate_summary_table <- function(data) {
  # Define function to subset data based on missing, below threshold, and ambiguous values
  subset_data <- function(data, column, value) {
    data <- subset(data, !(column == value))
  }
  
  # Subsetting the data 
  data <- subset_data(data, data$N_NATIVEpr.GFP.localization, "missing")
  data <- subset_data(data, data$C_Localization.combined, "missing")
  data <- subset_data(data, data$N_NATIVEpr.GFP.localization, "below threshold")
  data <- subset_data(data, data$C_Localization.combined, "below threshold")
  data <- subset_data(data, data$N_NATIVEpr.GFP.localization, "ambiguous")
  data <- subset_data(data, data$C_Localization.combined, "ambiguous")
  
  # Group by the 'Winner_20' column
  groups <- unique(data$Winner_20)
  
  # Create empty vectors to store counts
  counts <- c()
  same_counts <- c()
  diff_counts <- c()
  group_info <- c()
  essentiality_info <- c()
  
  # Iterate over groups
  for (group in groups) {
    # Filter data for the current group
    group_data <- subset(data, Winner_20 == group)
    
    # Create subsets for Non-Essential and Essential within the current group
    subset_non_essential <- subset(group_data, !(ORF_essential == "essential"))
    subset_essential <- subset(group_data, !(ORF_essential == "non-essential"))
    
    # Count occurrences of "Same" and "Different" for Non-Essential subset within the current group
    n_same_non_essential <- nrow(subset_non_essential[subset_non_essential$N_C_localization_endogenous_v2 == 'Same', ])
    n_different_non_essential <- nrow(subset_non_essential[subset_non_essential$N_C_localization_endogenous_v2 == 'Different', ])
    
    # Count occurrences of "Same" and "Different" for Essential subset within the current group
    n_same_essential <- nrow(subset_essential[subset_essential$N_C_localization_endogenous_v2 == 'Same', ])
    n_different_essential <- nrow(subset_essential[subset_essential$N_C_localization_endogenous_v2 == 'Different', ])
    
    # Append counts to vectors
    counts <- c(counts, n_same_non_essential + n_different_non_essential, n_same_essential + n_different_essential)
    same_counts <- c(same_counts, n_same_non_essential, n_same_essential)
    diff_counts <- c(diff_counts, n_different_non_essential, n_different_essential)
    group_info <- c(group_info, rep(group, 2))
    essentiality_info <- c(essentiality_info, "Non-Essential", "Essential")
  }
  
  # Combine counts into a data frame
  summary_table <- data.frame(
    Data = group_info,
    Protein_Essentiality = essentiality_info,
    NvsC_endogenous_Localization_Count = counts,
    NvsC_endogenous_Localization_Same = same_counts,
    NvsC_endogenous_Localization_Different = diff_counts
  )
  
  # Return summary table
  return(summary_table)
}

endogenous_table <- generate_summary_table(combined_plates_final_signals_SWAT_SP_MTS_localization)
endogenous_table$Data_type<-'Subset'
endogenous_table$ID<-paste(endogenous_table$Data_type,endogenous_table$Data,endogenous_table$Protein_Essentiality,sep="_")

#non endogenous
generate_summary_table <- function(data) {
  # Define function to subset data based on missing, below threshold, and ambiguous values
  subset_data <- function(data, column, value) {
    data <- subset(data, !(column == value))
  }
  
  # Subsetting the data 
  data <- subset_data(data, data$NOP1pr.GFP.localization, "missing")
  data <- subset_data(data, data$C_Localization.combined, "missing")
  data <- subset_data(data, data$NOP1pr.GFP.localization, "below threshold")
  data <- subset_data(data, data$C_Localization.combined, "below threshold")
  data <- subset_data(data, data$NOP1pr.GFP.localization, "ambiguous")
  data <- subset_data(data, data$C_Localization.combined, "ambiguous")
  
  # Group by the 'Winner_20' column
  groups <- unique(data$Winner_20)
  
  # Create empty vectors to store counts
  counts <- c()
  same_counts <- c()
  diff_counts <- c()
  group_info <- c()
  essentiality_info <- c()
  
  # Iterate over groups
  for (group in groups) {
    # Filter data for the current group
    group_data <- subset(data, Winner_20 == group)
    
    # Create subsets for Non-Essential and Essential within the current group
    subset_non_essential <- subset(group_data, !(ORF_essential == "essential"))
    subset_essential <- subset(group_data, !(ORF_essential == "non-essential"))
    
    # Count occurrences of "Same" and "Different" for Non-Essential subset within the current group
    n_same_non_essential <- nrow(subset_non_essential[subset_non_essential$N_C_localization_NNop1pr_C_v2 == 'Same', ])
    n_different_non_essential <- nrow(subset_non_essential[subset_non_essential$N_C_localization_NNop1pr_C_v2 == 'Different', ])
    
    # Count occurrences of "Same" and "Different" for Essential subset within the current group
    n_same_essential <- nrow(subset_essential[subset_essential$N_C_localization_NNop1pr_C_v2 == 'Same', ])
    n_different_essential <- nrow(subset_essential[subset_essential$N_C_localization_NNop1pr_C_v2 == 'Different', ])
    
    # Append counts to vectors
    counts <- c(counts, n_same_non_essential + n_different_non_essential, n_same_essential + n_different_essential)
    same_counts <- c(same_counts, n_same_non_essential, n_same_essential)
    diff_counts <- c(diff_counts, n_different_non_essential, n_different_essential)
    group_info <- c(group_info, rep(group, 2))
    essentiality_info <- c(essentiality_info, "Non-Essential", "Essential")
  }
  
  # Combine counts into a data frame
  summary_table <- data.frame(
    Data = group_info,
    Protein_Essentiality = essentiality_info,
    NNop1prvsC_Localization_Count = counts,
    NNop1prvsC_Localization_Same = same_counts,
    NNop1prvsC_Localization_Different = diff_counts
  )
  
  # Return summary table
  return(summary_table)
}

nonendogenous_table <- generate_summary_table(combined_plates_final_signals_SWAT_SP_MTS_localization)
colnames(nonendogenous_table) <- c("Data","Protein_Essentiality","NNop1prvsC_Localization_Count","NNop1prvsC_Localization_Same","NNop1prvsC_Localization_Different")
nonendogenous_table$Data_type<-'Subset'
nonendogenous_table$ID<-paste(nonendogenous_table$Data_type,nonendogenous_table$Data,nonendogenous_table$Protein_Essentiality,sep="_")

# Combine all dataframes
subset_20_table <- left_join(endogenous_table, nonendogenous_table, by="ID")
subset_20_table$Data.y<-NULL
subset_20_table$Data_type.y<-NULL
subset_20_table$Protein_Essentiality.y<-NULL
subset_20_table$Data_type<-subset_20_table$Data_type.x
subset_20_table$Data_type.x<-NULL
subset_20_table$ID2<-subset_20_table$ID
subset_20_table$ID<-NULL
subset_20_table$ID<-subset_20_table$ID2
subset_20_table$ID2<-NULL
names(subset_20_table)[1:2] <- c("Data","Protein_Essentiality")
threshold20_table<- rbind(All_20_table,subset_20_table)
threshold20_table<- threshold20_table[, c("Data_type", "Data","Protein_Essentiality","ID","NvsC_endogenous_Localization_Count","NvsC_endogenous_Localization_Different","NvsC_endogenous_Localization_Same","NNop1prvsC_Localization_Count","NNop1prvsC_Localization_Different","NNop1prvsC_Localization_Same")]
threshold20_table$Data[threshold20_table$Data=="N"]<-'C-loser'
threshold20_table$Data[threshold20_table$Data=="C"]<-'N-loser'
#write.csv(threshold20_table, paste0(path_input_dataset,"Fig1e_20.csv"))
#threshold20_table<-read.csv("Fig1e_20.csv")

# Calculate percentages and add new columns
threshold20_table$NvsC_endogenous_Localization_Different_Percentage <- (threshold20_table$NvsC_endogenous_Localization_Different / threshold20_table$NvsC_endogenous_Localization_Count) * 100
threshold20_table$NvsC_endogenous_Localization_Same_Percentage <- (threshold20_table$NvsC_endogenous_Localization_Same / threshold20_table$NvsC_endogenous_Localization_Count) * 100
threshold20_table$NNop1prvsC_Localization_Different_Percentage <- (threshold20_table$NNop1prvsC_Localization_Different / threshold20_table$NNop1prvsC_Localization_Count) * 100
threshold20_table$NNop1prvsC_Localization_Same_Percentage <- (threshold20_table$NNop1prvsC_Localization_Same / threshold20_table$NNop1prvsC_Localization_Count) * 100
threshold20_table$X<-NULL
threshold20_table$Category<-NULL
threshold20_table$ID<-NULL
colnames(threshold20_table)[2]<-"Winner_20"
threshold20_table$NvsC_endogenous_Localization_Count<-NULL
threshold20_table$NvsC_endogenous_Localization_Different<-NULL
threshold20_table$NvsC_endogenous_Localization_Same<-NULL
threshold20_table$NNop1prvsC_Localization_Count<-NULL
threshold20_table$NNop1prvsC_Localization_Different<-NULL
threshold20_table$NNop1prvsC_Localization_Same<-NULL
write.csv(threshold20_table, paste0(path_input_dataset,"Fig1e_20_PlotTable.csv"))
Fig1e_20_PlotTable<-read.csv("Fig1e_20_PlotTable.csv")
Fig1e_20_PlotTable$X<-NULL

# Reshape the data
Fig1e_20_PlotTable <- pivot_longer(Fig1e_20_PlotTable, 
                                   cols = starts_with("NvsC") | starts_with("NNop1prvsC"), 
                                   names_to = "Category", 
                                   values_to = "Percentage")

# Reorder the columns
Fig1e_20_PlotTable <- Fig1e_20_PlotTable[, c("Data_type", "Protein_Essentiality", "Winner_20", "Category", "Percentage")]

# Write the formatted data to a new CSV file
write.csv(Fig1e_20_PlotTable, "Fig1e_20_PlotTable_LongFormat.csv", row.names = FALSE)

#PLOT for Endogenous Subset data
Fig1e_20_PlotTable_Subset<-subset(Fig1e_20_PlotTable,!(Fig1e_20_PlotTable$Data_type=="All"))
Fig1e_20_PlotTable_Subset_Endogenous<-subset(Fig1e_20_PlotTable_Subset,!(Fig1e_20_PlotTable_Subset$Category=="NNop1prvsC_Localization_Different_Percentage"))
Fig1e_20_PlotTable_Subset_Endogenous<-subset(Fig1e_20_PlotTable_Subset_Endogenous,!(Fig1e_20_PlotTable_Subset_Endogenous$Category=="NNop1prvsC_Localization_Same_Percentage"))

e_20<-ggplot(Fig1e_20_PlotTable_Subset_Endogenous, aes(x = Percentage, y = Winner_20, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Protein_Essentiality, scales = "free") +
  labs(title = "% ORFs tagged Endogenously with different Localization 20% threshold",
       x = "% ORFs", y = "Loser type") +
  scale_fill_manual(values = c("red", "gray")) +
  theme_minimal()

#PLOT for Non-Endogenous Subset data
Fig1e_20_PlotTable_Subset<-subset(Fig1e_20_PlotTable,!(Fig1e_20_PlotTable$Data_type=="All"))
Fig1e_20_PlotTable_Subset_NonEndogenous<-subset(Fig1e_20_PlotTable_Subset,!(Fig1e_20_PlotTable_Subset$Category=="NvsC_endogenous_Localization_Different_Percentage"))
Fig1e_20_PlotTable_Subset_NonEndogenous<-subset(Fig1e_20_PlotTable_Subset_NonEndogenous,!(Fig1e_20_PlotTable_Subset_NonEndogenous$Category=="NvsC_endogenous_Localization_Same_Percentage"))

ne_20<-ggplot(Fig1e_20_PlotTable_Subset_NonEndogenous, aes(x = Percentage, y = Winner_20, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Protein_Essentiality, scales = "free") +
  labs(title = "% ORFs tagged Non-Endogenously with different Localization 20% threshold",
       x = "% ORFs", y = "Loser type") +
  scale_fill_manual(values = c("red", "gray")) +
  theme_minimal()

cairo_pdf(filename = "Fig1e_20_Plots.pdf", 
          width = 11, height = 8, onefile = TRUE)
e_20
ne_20
dev.off()

#_______________________________________________________________________________________________
#FigS1a. Tagging efficiency of Reference ORFs tagged with mNG endogenously, parallel to the Halo libraries SGA

#Install packages if not present (needs internet connection)
library(BiocManager)
BiocManager::install_version("flowCore", version = "2.12.2")
BiocManager::install_version("flowStats", version = "4.12.0")
BiocManager::install_version("flowViz", version = "1.64.0")
BiocManager::install_version("ggcyto", version = "1.28.1")
BiocManager::install_version("openCyto", version = "2.12.0")

# Load packages
library(tidyverse)
library(flowCore)
library(flowStats)
library(flowViz)
library(ggcyto)
library(openCyto) 
library(gridExtra) 
library(vioplot)

# Load custom functions needed
cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix(NA, mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
} #Function to bind columns together

###Tagging efficiency of C-tagged Reference ORFs
# Data import
#Change source to subfolder of choice; where your ".fcs" files are located.
sub <- "./23.05.2022_N-Ceffect-ScreenI_mCh-mNG_C-ORFs-mNG_Ctrls"

# Find file names of .fcs files in the current working directory
# Need the "./" before the folder if entering a folder, change the number at "^number(.*)..."
FileNames <- list.files(path=sub, pattern = "^(.*).fcs$", recursive = FALSE) 
FileNames
#Remove a well that was done as an extra test.
FileNames2 <- FileNames[c(1:14)]
FileNames2

#Prepare this before hand, labels for strains, has to be unique
StrainNames <- read.csv("Strainlabels_mNG_23.05.2022_Screen.csv", header = TRUE, sep = ",") 
SN <- as.character(StrainNames$ORF)
SN

# Compiling all fcs files in directory into one dataset
#flowData <- read.flowSet(FileNames, path = sub)
flowData <- read.flowSet(FileNames2, path = sub)
sampleNames(flowData)

#Changing the sample names from the file name
sampleNames(flowData) <- SN 
sampleNames(flowData)


# Data prep

#Showing the parameters measured for the first file in the folder
pData(parameters(flowData[[1]])) 

for (i in SN) {
  #flowData@frames[[i]]@parameters@data["$P10",1] <- "Ratio"
  #flowData@frames[[i]]@description[["$P10N"]] <- "Ratio"
  #colnames(flowData@frames[[i]]@exprs)[10] <- "Ratio"
  
  flowData@frames[[i]]@parameters@data[, "desc"] <- c("Size-A", "Size-H", "Size-W",
                                                      "Complexity-A", "Complexity-H",
                                                      "Singularity-W", "mNG", "mCH",
                                                      "DAPI", "NA")
  
}
pData(parameters(flowData[[1]])) #Check that the parameters are changed correctly

#Check that the table has correct number of experiments/samples
flowData

# Transformation of data

tData <- transform(flowData, transformList(colnames(flowData)[7:8], logicleTransform())) #logicle transformation for fluorescence data only
#print(splom(flowData[[4]])) #to compare, will show the fluorescence graphs squished to side (original data)
#print(splom(tData[[4]])) #to compare, should show better distribution of the fluorescence graphs (transformed data)
#print(xyplot(`mNG` ~ `Size-A`, data = tData))


# Gating the data
## Setup

# Setting up the channels/axes
chnl1 <- c("FSC-A", "SSC-A") #This has to refer to column names for a flow set, not the description
chnl2 <- c("SSC-A", "SSC-W") #This has to refer to column names for a flow set, not the description
chnl3 <- c("FSC-A", "BL488nm 530/30-A") #This has to refer to column names for a flow set, not the description
chnl4 <- c("FSC-A", "YG561nm 610/20-A") #This has to refer to column names for a flow set, not the description
chnl5 <- "SSC-A" #This has to refer to column names for a flow set, not the description
chnl6 <- "SSC-W" #This has to refer to column names for a flow set, not the description

# Setting up the gates (use next paragraph to view)
c_sample <- 14
Flive <- openCyto:::.boundary(tData[[c_sample]], channels = chnl1, min = c(0.35e5, 0.15e5), max=c(2e5, 1.75e5)) #adjust to remove cell debris
Fsingle <- openCyto:::.boundary(tData[[c_sample]], channels = chnl2, min = c(0.2e5, 0.62e5), max=c(1.85e5, 0.85e5)) #adjust to select single cells
Fgreen <- openCyto:::.boundary(tData[[c_sample]], channels = chnl3, min = c(0.35e5, 2), max=c(2e5, Inf)) #adjust to select above background
Fred <- openCyto:::.boundary(tData[[c_sample]], channels = chnl4, min = c(0.35e5, 1), max=c(2e5, Inf)) #adjust to select above background
FBGG <- openCyto:::.boundary(tData[[c_sample]], channels = chnl3, min = c(0.35e5, 0), max=c(2e5, 2)) #adjust to select above background
FBGR <- openCyto:::.boundary(tData[[c_sample]], channels = chnl4, min = c(0.35e5, 0), max=c(2e5, 1)) #adjust to select above background

# For viewing the gates
p1 <- ggcyto(data = tData[[c_sample]], aes(x = UQ(chnl1[1]), y = UQ(chnl1[2]))) +
  geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 3e5))) # change chnlnumber according to the gate
p2 <- ggcyto(data = tData[[c_sample]], aes(x = UQ(chnl2[1]), y = UQ(chnl2[2]))) +
  geom_hex(bins = 100)  + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0.5e5, 2.5e5))) # change chnlnumber according to the gate
p3 <- ggcyto(data = tData[[c_sample]], aes(x = UQ(chnl3[1]), y = UQ(chnl3[2]))) +
  geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 5))) # change chnlnumber according to the gate
p4 <- ggcyto(data = tData[[c_sample]], aes(x = UQ(chnl4[1]), y = UQ(chnl4[2]))) +
  geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 5))) # change chnlnumber according to the gate

#view plot with gate; change the geom_gate as needed
grid.arrange(as.ggplot(p1 + geom_gate(Flive)),
             as.ggplot(p2 + geom_gate(Fsingle)),
             as.ggplot(p3 + geom_gate(FBGG)),
             as.ggplot(p4 + geom_gate(FBGR)),
             nrow = 2, ncol = 2)

#p1 #view the plot without gate
## Once complete, transfer details to .csv file manually

#Only run to save
pdf(file = paste0("Output/",(format(Sys.time(), "%Y%m%d_%H%M")), "_Gate_plots.pdf"), width = 12, height = 10, onefile = T)

for (id in 1:length(SN)) {
  # For viewing the gates
  p1 <- ggcyto(data = tData[[id]], aes(x = UQ(chnl1[1]), y = UQ(chnl1[2]))) +
    geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 3e5))) # change chnlnumber according to the gate
  p2 <- ggcyto(data = tData[[id]], aes(x = UQ(chnl2[1]), y = UQ(chnl2[2]))) +
    geom_hex(bins = 100)  + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0.5e5, 2.5e5))) # change chnlnumber according to the gate
  p3 <- ggcyto(data = tData[[id]], aes(x = UQ(chnl3[1]), y = UQ(chnl3[2]))) +
    geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 5))) # change chnlnumber according to the gate
  p4 <- ggcyto(data = tData[[id]], aes(x = UQ(chnl4[1]), y = UQ(chnl4[2]))) +
    geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 5))) # change chnlnumber according to the gate
  
  #view plot with gate; change the geom_gate as needed
  grid.arrange(as.ggplot(p1 + geom_gate(Flive)),
               as.ggplot(p2 + geom_gate(Fsingle)),
               as.ggplot(p3 + geom_gate(FBGG)),
               as.ggplot(p4 + geom_gate(FBGR)),
               nrow = 2, ncol = 2)
}

dev.off()

## Applying the Gates

# Automated gating using a template for samples
GT = gatingTemplate(x = "./Gates_mNG.csv", verbose=TRUE) #gating template for automation, template file needs to be in directory or sub-directory
#plot(GT) #show the gating sequence and all possible/gated nodes and types of gating method
GS = GatingSet(tData) #Construct gating set with tData file number [a:z] or dont put [] if all data

gt_gating(GT, GS) #apply the gates in template (GT) to the flowdata (GS)
warn = warnings() #save any warning messages. MUST CHECK!

plot(GS) #Shows gates applied
ggcyto::autoplot(GS[[1]]) #need double square bracket; plot all gates for the first file; check label name
ggcyto::autoplot(GS[[length(tData)]]) #plot the last to make sure no files missing

nodelist = gs_get_pop_paths(GS, path = "full") #getting the list of full path of gates applied
nodelist #Showing the list. Don't run this if there is a lot of gates

#getGate(GS[[1]], nodelist[5]) #Retrieve a certain gate for GS[[filenumber]] at the node[number]. Don't run without filenumber if there is a lot of files
#plotGate(GS[[1]], nodelist[3]) #Brings up a specific plot of the file in GS[[file.number]] at node[number]

#Density plots for the transformed and gated data
#DPred = densityplot(~ `mCH`, getData(GS, nodelist[3]), overlap=-0.5, xlab='Intensity', ylab='ORF') #value after ~ is the fluorescence data, node[number] corresponds to the gate in list
#DPred #show the densityplot
#DPgreen = densityplot(~ `mNG`, getData(GS, nodelist[3]), overlap=-0.5, xlab='Intensity', ylab='ORF') #value after ~ is the fluorescence data, node[number] corresponds to the gate in list
#DPgreen #show the densityplot
#grid.arrange(DPred,DPgreen, ncol=2) #arrange the two density plots side by side 
#The lattice plots do not use the par settings in general. They have their own set of settings from Grid graphics. 
#Use gridExtra package with grid.arrange().
#densityplot(~ `Singulari`, getData(GS[[1]], nodelist[3]))

# Data analysis

## Data counts

Countdat = gs_pop_get_count_fast(GS) #Obtain the counts for the gated and parent populations
Countdat$Percent = round((Countdat$Count/Countdat$ParentCount)*100, digits=1)
CDatF = Countdat %>% dplyr::filter(grepl("redFP", Population) | grepl("greenFP", Population)) #Only fluorescence percentage

setwd("Output")
write.csv(Countdat, file=paste((format(Sys.time(), "%Y_%m_%d_%H%M")),"_CountData.csv", sep = ""), row.names=FALSE) # Give a name for your output .csv file (numbered by date and time)
write.csv(CDatF, file=paste((format(Sys.time(), "%Y_%m_%d_%H%M")),"_FluoData.csv", sep = ""), row.names=FALSE) # Give a name for your output .csv file (numbered by date and time)

setwd(path_input_dataset)


## Data visualisation

oripar = par() # get original par
par(mar = c(5.1, 13, 4.1, 2.1)) #Adjust the left margin bigger

## For mCherry ##
#nodelist need to check whether 3("/livecells/singlecells")
#bracket outside of "as.data.frame" should correspond to the column with the fluo values
Vred = as.data.frame(exprs(gh_pop_get_data(GS[[1]], nodelist[3])))[ c(8) ] 
for (i in 2:(length(flowData))) {
  fdata = as.data.frame(exprs(gh_pop_get_data(GS[[i]], nodelist[3])))[ c(8) ] #nodelist need to check whether 3("/livecells/singlecells")
  Vred=cbindPad(Vred, fdata)
}
colnames(Vred)=StrainNames[,1] #in square brackets [row,col]; Changing the column names to the ORFs

mcherry = CDatF %>% dplyr::filter(grepl("redFP", Population))

par(mar = c(5.1, 13, 4.1, 2.1)) #Adjust the left margin bigger
(plot(0,yaxt='n',pch='',ylab='',xlab=expression(Log[10]*"(mCherry Intensity)(a.u.)"), xlim = c(0,5), ylim =c(0,length(flowData))) 
  + axis(2, seq_len(length(flowData)), SN, las=1) 
  + (grid(nx = NULL, ny = NA, col = "lightgray") 
     + axis(4, seq_len(length(flowData)), mcherry$Percent, las=1, line=-5, tick=0)) 
  + title(ylab='ORF', mgp=c(11,1,1)) 
  + title(main=paste("Fluorescence Intensity of mCherry", sep = ""))) #Change title as per needed

for (u in 1:ncol(Vred)) {
  Vred1 = as.numeric(Vred[ which(Vred[,u] > 0), u])
  vioplot(Vred1, h = 0.05, ylim=c(0,length(flowData)), col = "red", add = TRUE, at=u, horizontal=TRUE)}
abline(v=2, col="red", lty=2, lwd=3) #Change this to threshold value

setwd("Output")
#Export image as A3 size
dev.copy2pdf(device=pdf, file=paste((format(Sys.time(), "%Y_%m_%d_%H%M")),"Ctagged_mCherryVioplot.pdf", sep = ""), width=11.7, height=16.5)

setwd(path_input_dataset)

## For GFP ##
Vgreen = as.data.frame(exprs(gh_pop_get_data(GS[[1]], nodelist[3])))[ c(7) ] #nodelist need to check whether 3 or 5 (singlecell or fluo gate)
for (i in 2:(length(flowData))) {
  fdata = as.data.frame(exprs(gh_pop_get_data(GS[[i]], nodelist[3])))[ c(7) ] #nodelist need to check whether 3 or 7 (singlecell or fluo gate)
  Vgreen=cbindPad(Vgreen, fdata)
}
colnames(Vgreen)=StrainNames[,1] #in square brackets [row,col]

GFP=CDatF %>% dplyr::filter(grepl("greenFP", Population))

par(mar = c(5.1, 13, 4.1, 2.1)) #Adjust the left margin bigger
(plot(0,yaxt='n',pch='',ylab='',xlab=expression(Log[10]*"(mNG Intensity)(a.u.)"), xlim = c(0,5), ylim =c(0,length(flowData))) 
  + axis(2, seq_len(length(flowData)), SN, las=1) 
  + (grid(nx = NULL, ny = NA, col = "lightgray") 
     + axis(4, seq_len(length(flowData)), GFP$Percent, las=1, line=-5, tick=0)) 
  + title(ylab='ORF', mgp=c(11,1,1)) 
  + title(main=paste("Fluorescence Intensity of mNG", sep = ""))) #Change title as per needed

for (u in 1:ncol(Vgreen)) {
  Vgreen1 = as.numeric(Vgreen[ which(Vgreen[,u] > 0), u])
  vioplot(Vgreen1, h = 0.05, ylim=c(0,length(flowData)), col = "green", add = TRUE, at=u, horizontal=TRUE)}
abline(v=2, col="red", lty=2, lwd=3) #Change this to threshold value

setwd("Output")
#Export image as A3 size
dev.copy2pdf(device=pdf, file=paste((format(Sys.time(), "%Y_%m_%d_%H%M")),"Ctagged_mNGVioplot.pdf", sep = ""), width=11.7, height=16.5)

setwd(path_input_dataset)
par(oripar) #reset par


###Tagging efficiency of N-tagged Reference ORFs
# Data import
#Change source to subfolder of choice; where your ".fcs" files are located.
sub <- "./23.05.2022_N-Ceffect-ScreenI_mCh-mNG_N-ORFs-mNG_Ctrls"

# Find file names of .fcs files in the current working directory
# Need the "./" before the folder if entering a folder, change the number at "^number(.*)..."
FileNames <- list.files(path=sub, pattern = "^(.*).fcs$", recursive = FALSE) 
FileNames
#Remove a well that was done as an extra test.
FileNames2 <- FileNames[c(1:14)]
FileNames2

#Prepare this before hand, labels for strains, has to be unique
StrainNames <- read.csv("Strainlabels_mNG_23.05.2022_Screen.csv", header = TRUE, sep = ",") 
SN <- as.character(StrainNames$ORF)
SN

# Compiling all fcs files in directory into one dataset
#flowData <- read.flowSet(FileNames, path = sub)
flowData <- read.flowSet(FileNames2, path = sub)
sampleNames(flowData)

#Changing the sample names from the file name
sampleNames(flowData) <- SN 
sampleNames(flowData)

# Data prep

#Showing the parameters measured for the first file in the folder
pData(parameters(flowData[[1]])) 

for (i in SN) {
  #flowData@frames[[i]]@parameters@data["$P10",1] <- "Ratio"
  #flowData@frames[[i]]@description[["$P10N"]] <- "Ratio"
  #colnames(flowData@frames[[i]]@exprs)[10] <- "Ratio"
  
  flowData@frames[[i]]@parameters@data[, "desc"] <- c("Size-A", "Size-H", "Size-W",
                                                      "Complexity-A", "Complexity-H",
                                                      "Singularity-W", "mNG", "mCH",
                                                      "DAPI", "NA")
  
}
pData(parameters(flowData[[1]])) #Check that the parameters are changed correctly

#Check that the table has correct number of experiments/samples
flowData



# Transformation of data

tData <- transform(flowData, transformList(colnames(flowData)[7:8], logicleTransform())) #logicle transformation for fluorescence data only
#print(splom(flowData[[4]])) #to compare, will show the fluorescence graphs squished to side (original data)
#print(splom(tData[[4]])) #to compare, should show better distribution of the fluorescence graphs (transformed data)
#print(xyplot(`mNG` ~ `Size-A`, data = tData))

# Gating the data

## Setup

# Setting up the channels/axes
chnl1 <- c("FSC-A", "SSC-A") #This has to refer to column names for a flow set, not the description
chnl2 <- c("SSC-A", "SSC-W") #This has to refer to column names for a flow set, not the description
chnl3 <- c("FSC-A", "BL488nm 530/30-A") #This has to refer to column names for a flow set, not the description
chnl4 <- c("FSC-A", "YG561nm 610/20-A") #This has to refer to column names for a flow set, not the description
chnl5 <- "SSC-A" #This has to refer to column names for a flow set, not the description
chnl6 <- "SSC-W" #This has to refer to column names for a flow set, not the description

# Setting up the gates (use next paragraph to view)
c_sample <- 14
Flive <- openCyto:::.boundary(tData[[c_sample]], channels = chnl1, min = c(0.35e5, 0.15e5), max=c(2e5, 1.75e5)) #adjust to remove cell debris
Fsingle <- openCyto:::.boundary(tData[[c_sample]], channels = chnl2, min = c(0.2e5, 0.62e5), max=c(1.85e5, 0.85e5)) #adjust to select single cells
Fgreen <- openCyto:::.boundary(tData[[c_sample]], channels = chnl3, min = c(0.35e5, 2), max=c(2e5, Inf)) #adjust to select above background
Fred <- openCyto:::.boundary(tData[[c_sample]], channels = chnl4, min = c(0.35e5, 1), max=c(2e5, Inf)) #adjust to select above background
FBGG <- openCyto:::.boundary(tData[[c_sample]], channels = chnl3, min = c(0.35e5, 0), max=c(2e5, 2)) #adjust to select above background
FBGR <- openCyto:::.boundary(tData[[c_sample]], channels = chnl4, min = c(0.35e5, 0), max=c(2e5, 1)) #adjust to select above background

# For viewing the gates
p1 <- ggcyto(data = tData[[c_sample]], aes(x = UQ(chnl1[1]), y = UQ(chnl1[2]))) +
  geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 3e5))) # change chnlnumber according to the gate
p2 <- ggcyto(data = tData[[c_sample]], aes(x = UQ(chnl2[1]), y = UQ(chnl2[2]))) +
  geom_hex(bins = 100)  + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0.5e5, 2.5e5))) # change chnlnumber according to the gate
p3 <- ggcyto(data = tData[[c_sample]], aes(x = UQ(chnl3[1]), y = UQ(chnl3[2]))) +
  geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 5))) # change chnlnumber according to the gate
p4 <- ggcyto(data = tData[[c_sample]], aes(x = UQ(chnl4[1]), y = UQ(chnl4[2]))) +
  geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 5))) # change chnlnumber according to the gate

#view plot with gate; change the geom_gate as needed
grid.arrange(as.ggplot(p1 + geom_gate(Flive)),
             as.ggplot(p2 + geom_gate(Fsingle)),
             as.ggplot(p3 + geom_gate(FBGG)),
             as.ggplot(p4 + geom_gate(FBGR)),
             nrow = 2, ncol = 2)

#p1 #view the plot without gate
## Once complete, transfer details to .csv file manually



#Only run to save
pdf(file = paste0("Output/",(format(Sys.time(), "%Y%m%d_%H%M")), "_Gate_plots.pdf"), width = 12, height = 10, onefile = T)

for (id in 1:length(SN)) {
  # For viewing the gates
  p1 <- ggcyto(data = tData[[id]], aes(x = UQ(chnl1[1]), y = UQ(chnl1[2]))) +
    geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 3e5))) # change chnlnumber according to the gate
  p2 <- ggcyto(data = tData[[id]], aes(x = UQ(chnl2[1]), y = UQ(chnl2[2]))) +
    geom_hex(bins = 100)  + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0.5e5, 2.5e5))) # change chnlnumber according to the gate
  p3 <- ggcyto(data = tData[[id]], aes(x = UQ(chnl3[1]), y = UQ(chnl3[2]))) +
    geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 5))) # change chnlnumber according to the gate
  p4 <- ggcyto(data = tData[[id]], aes(x = UQ(chnl4[1]), y = UQ(chnl4[2]))) +
    geom_hex(bins = 100) + ggcyto_par_set(limits = list(x = c(0, 3e5), y = c(0, 5))) # change chnlnumber according to the gate
  
  #view plot with gate; change the geom_gate as needed
  grid.arrange(as.ggplot(p1 + geom_gate(Flive)),
               as.ggplot(p2 + geom_gate(Fsingle)),
               as.ggplot(p3 + geom_gate(FBGG)),
               as.ggplot(p4 + geom_gate(FBGR)),
               nrow = 2, ncol = 2)
}

dev.off()


## Applying the Gates


# Automated gating using a template for samples
GT = gatingTemplate(x = "./Gates_mNG.csv", verbose=TRUE) #gating template for automation, template file needs to be in directory or sub-directory
#plot(GT) #show the gating sequence and all possible/gated nodes and types of gating method
GS = GatingSet(tData) #Construct gating set with tData file number [a:z] or dont put [] if all data

gt_gating(GT, GS) #apply the gates in template (GT) to the flowdata (GS)
warn = warnings() #save any warning messages. MUST CHECK!

plot(GS) #Shows gates applied
ggcyto::autoplot(GS[[1]]) #need double square bracket; plot all gates for the first file; check label name
ggcyto::autoplot(GS[[length(tData)]]) #plot the last to make sure no files missing

nodelist = gs_get_pop_paths(GS, path = "full") #getting the list of full path of gates applied
nodelist #Showing the list. Don't run this if there is a lot of gates

#getGate(GS[[1]], nodelist[5]) #Retrieve a certain gate for GS[[filenumber]] at the node[number]. Don't run without filenumber if there is a lot of files
#plotGate(GS[[1]], nodelist[3]) #Brings up a specific plot of the file in GS[[file.number]] at node[number]

#Density plots for the transformed and gated data
#DPred = densityplot(~ `mCH`, getData(GS, nodelist[3]), overlap=-0.5, xlab='Intensity', ylab='ORF') #value after ~ is the fluorescence data, node[number] corresponds to the gate in list
#DPred #show the densityplot
#DPgreen = densityplot(~ `mNG`, getData(GS, nodelist[3]), overlap=-0.5, xlab='Intensity', ylab='ORF') #value after ~ is the fluorescence data, node[number] corresponds to the gate in list
#DPgreen #show the densityplot
#grid.arrange(DPred,DPgreen, ncol=2) #arrange the two density plots side by side 
#The lattice plots do not use the par settings in general. They have their own set of settings from Grid graphics. 
#Use gridExtra package with grid.arrange().
#densityplot(~ `Singulari`, getData(GS[[1]], nodelist[3]))


# Data analysis

## Data counts

Countdat = gs_pop_get_count_fast(GS) #Obtain the counts for the gated and parent populations
Countdat$Percent = round((Countdat$Count/Countdat$ParentCount)*100, digits=1)
CDatF = Countdat %>% dplyr::filter(grepl("redFP", Population) | grepl("greenFP", Population)) #Only fluorescence percentage

setwd("Output")
write.csv(Countdat, file=paste((format(Sys.time(), "%Y_%m_%d_%H%M")),"_CountData.csv", sep = ""), row.names=FALSE) # Give a name for your output .csv file (numbered by date and time)
write.csv(CDatF, file=paste((format(Sys.time(), "%Y_%m_%d_%H%M")),"_FluoData.csv", sep = ""), row.names=FALSE) # Give a name for your output .csv file (numbered by date and time)

setwd(path_input_dataset)


## Data visualisation

oripar = par() # get original par
par(mar = c(5.1, 13, 4.1, 2.1)) #Adjust the left margin bigger

## For mCherry ##
#nodelist need to check whether 3("/livecells/singlecells")
#bracket outside of "as.data.frame" should correspond to the column with the fluo values
Vred = as.data.frame(exprs(gh_pop_get_data(GS[[1]], nodelist[3])))[ c(8) ] 
for (i in 2:(length(flowData))) {
  fdata = as.data.frame(exprs(gh_pop_get_data(GS[[i]], nodelist[3])))[ c(8) ] #nodelist need to check whether 3("/livecells/singlecells")
  Vred=cbindPad(Vred, fdata)
}
colnames(Vred)=StrainNames[,1] #in square brackets [row,col]; Changing the column names to the ORFs

mcherry = CDatF %>% dplyr::filter(grepl("redFP", Population))

par(mar = c(5.1, 13, 4.1, 2.1)) #Adjust the left margin bigger
(plot(0,yaxt='n',pch='',ylab='',xlab=expression(Log[10]*"(mCherry Intensity)(a.u.)"), xlim = c(0,5), ylim =c(0,length(flowData))) 
  + axis(2, seq_len(length(flowData)), SN, las=1) 
  + (grid(nx = NULL, ny = NA, col = "lightgray") 
     + axis(4, seq_len(length(flowData)), mcherry$Percent, las=1, line=-5, tick=0)) 
  + title(ylab='ORF', mgp=c(11,1,1)) 
  + title(main=paste("Fluorescence Intensity of mCherry", sep = ""))) #Change title as per needed

for (u in 1:ncol(Vred)) {
  Vred1 = as.numeric(Vred[ which(Vred[,u] > 0), u])
  vioplot(Vred1, h = 0.05, ylim=c(0,length(flowData)), col = "red", add = TRUE, at=u, horizontal=TRUE)}
abline(v=2, col="red", lty=2, lwd=3) #Change this to threshold value

setwd("Output")
#Export image as A3 size
dev.copy2pdf(device=pdf, file=paste((format(Sys.time(), "%Y_%m_%d_%H%M")),"Ntagged_mCherryVioplot.pdf", sep = ""), width=11.7, height=16.5)

setwd(path_input_dataset)

## For GFP ##
Vgreen = as.data.frame(exprs(gh_pop_get_data(GS[[1]], nodelist[3])))[ c(7) ] #nodelist need to check whether 3 or 5 (singlecell or fluo gate)
for (i in 2:(length(flowData))) {
  fdata = as.data.frame(exprs(gh_pop_get_data(GS[[i]], nodelist[3])))[ c(7) ] #nodelist need to check whether 3 or 7 (singlecell or fluo gate)
  Vgreen=cbindPad(Vgreen, fdata)
}
colnames(Vgreen)=StrainNames[,1] #in square brackets [row,col]

GFP=CDatF %>% dplyr::filter(grepl("greenFP", Population))

par(mar = c(5.1, 13, 4.1, 2.1)) #Adjust the left margin bigger
(plot(0,yaxt='n',pch='',ylab='',xlab=expression(Log[10]*"(mNG Intensity)(a.u.)"), xlim = c(0,5), ylim =c(0,length(flowData))) 
  + axis(2, seq_len(length(flowData)), SN, las=1) 
  + (grid(nx = NULL, ny = NA, col = "lightgray") 
     + axis(4, seq_len(length(flowData)), GFP$Percent, las=1, line=-5, tick=0)) 
  + title(ylab='ORF', mgp=c(11,1,1)) 
  + title(main=paste("Fluorescence Intensity of mNG", sep = ""))) #Change title as per needed

for (u in 1:ncol(Vgreen)) {
  Vgreen1 = as.numeric(Vgreen[ which(Vgreen[,u] > 0), u])
  vioplot(Vgreen1, h = 0.05, ylim=c(0,length(flowData)), col = "green", add = TRUE, at=u, horizontal=TRUE)}
abline(v=2, col="red", lty=2, lwd=3) #Change this to threshold value

setwd("Output")
#Export image as A3 size
dev.copy2pdf(device=pdf, file=paste((format(Sys.time(), "%Y_%m_%d_%H%M")),"Ntagged_mNGVioplot.pdf", sep = ""), width=11.7, height=16.5)

setwd(path_input_dataset)

par(oripar) #reset par

#_______________________________________________________________________________________________
#FigS1b. 
#Reproducibility of relative fitness estimation with the competition assay.
#Relative fitnesses of two replicates for 92 ORFs, determined as in Fig.1a. r, pearson correlation coefficient. 

##Correlation between P1 and P59

combined_plates_final_P1P59_t0<-read.csv("combined_plates_final_P1P59_t0.csv")

reproducibility<-ggplot(combined_plates_final_P1P59_t0, aes(P1_Slope_lm , P59_Slope_lm))+
  geom_point()+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  xlab("Fitness effect Replicate 1")+ylab("Fitness effect Replicate 2")+
  #stat_poly_eq(use_label(c("eq", "R2"))) +
  stat_poly_line() +
  stat_poly_eq() +
  ggtitle("Correlation of fitness effect of replicates")+
  theme_classic()

cairo_pdf(filename = "FigS1b_Screen_Reproducibility.pdf", 
          width = 11, height = 8, onefile = TRUE)
reproducibility
dev.off()

#_______________________________________________________________________________________________
#FigS1b. 
#Noise determination for flow cytometry
combined_plates_final_noise_t0<-read.csv("combined_plates_final_noise_t0.csv")

#General plot
plotS1b<-ggplot(combined_plates_final_noise_t0, aes(x=Slope_lm, fill=ORF)) +
  geom_histogram(binwidth=.01, bins=100, colour="black") +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Slope_lm")+ylab("count")+
  ggtitle("Fitness effect distribution for different strains at t0")+
  xlim(-1,1)+
  theme_classic()

cairo_pdf(filename = "FigS1b_NoiseCytometer.pdf", 
          width = 11, height = 8, onefile = TRUE)
plotS1b
dev.off()

#Individual SD and Mean relative fitness per strain-pair
Donors_t0<-subset(combined_plates_final_noise_t0, combined_plates_final_noise_t0$ORF == "Donors")
Mean_Slope_Donors_t0<-mean(Donors_t0$Slope_lm, na.rm=T)
sd_Slope_Donors_t0<-sd(Donors_t0$Slope_lm, na.rm=T)

YPL093W_t0<-subset(combined_plates_final_noise_t0, combined_plates_final_noise_t0$ORF == "YPL093W")
Mean_Slope_YPL093W_t0<-mean(YPL093W_t0$Slope_lm, na.rm=T)
sd_Slope_YPL093W_t0<-sd(YPL093W_t0$Slope_lm, na.rm=T)

YLR024C_t0<-subset(combined_plates_final_noise_t0, combined_plates_final_noise_t0$ORF == "YLR024C")
Mean_Slope_YLR024C_t0<-mean(YLR024C_t0$Slope_lm, na.rm=T)
sd_Slope_YLR024C_t0<-sd(YLR024C_t0$Slope_lm, na.rm=T)

YER177W_t0<-subset(combined_plates_final_noise_t0, combined_plates_final_noise_t0$ORF == "YER177W")
Mean_Slope_YER177W_t0<-mean(YER177W_t0$Slope_lm, na.rm=T)
sd_Slope_YER177W_t0<-sd(YER177W_t0$Slope_lm, na.rm=T)

#_______________________________________________________________________________________________
#FigS1d & FigS1e. 
#Relative fitness of Halo-tagged strains determined according to Fig.1a and in the triple competition assay
combined_plates_final_Triple<- read.csv("combined_plates_final_Triple_ValidationSet.csv")
combined_plates_final_Triple$X<-NULL

Screen_vs_Triple_sub<-read.csv("Screen_vs_Triple_sub.csv")
Screen_vs_Triple_sub$X<-NULL

plotS1d<-ggplot(Screen_vs_Triple_sub, aes(Screen_vs_Triple_sub$Slope_lm , Screen_vs_Triple_sub$Slope_CN_Triple))+
  geom_point()+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  xlab("Slopes Screen")+ylab("Slopes Triple competition C vs N")+
  #stat_poly_eq(use_label(c("eq", "R2"))) +
  stat_poly_line() +
  stat_poly_eq() +
  ggtitle("Correlation of Screen Validation set and Triple competition C vs N comparison")+
  xlim(-1.5,1)+
  ylim(-1.5,1)+
  theme_classic()

cairo_pdf(filename = "FigS1d.pdf", 
          width = 11, height = 8, onefile = TRUE)
plotS1d
dev.off()

#Pairwise relative fitness of the three strains for each ORF in the triple competition assay
combined_plates_final_Triple$ORF<- factor(combined_plates_final_Triple$ORF, levels=unique(combined_plates_final_Triple$ORF))

plotS1e<-ggplot(combined_plates_final_Triple, aes(fill=Slope_type, y=Slope, x=ORF)) +
  geom_bar(stat = "identity",position = "dodge")+
  xlab("ORF") +
  ylab("Slope of linear regression of pairwise comparisons") +
  #geom_text(aes(label=C_Localization_Percentage), position=position_dodge(width=0.9), vjust=-0.25,size=1.9)+
  #theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic()

cairo_pdf(filename = "FigS1e.pdf", 
          width = 11, height = 8, onefile = TRUE)
plotS1e
dev.off()

#___________________________________________________
#FigS1f.
#Protein localization signals analysis
combined_plates_final_signals_SWAT_SP_MTS<-read.csv("combined_plates_final_signals_SWAT_SP_MTS.csv")
Mean_Slope<-mean(combined_plates_final_signals_SWAT_SP_MTS$Slope_lm)

# plots not included in the manuscript
ggplot(combined_plates_final_signals_SWAT_SP_MTS, aes(x=Slope_lm)) +
  facet_wrap(~Signal_Protein_Terminal, ncol=1)+
  geom_histogram(binwidth=.01, bins=1000, colour="black", fill="black") +
  geom_vline(xintercept=Mean_Slope, size=1.5, color="red", linetype="dotted")+
  scale_y_log10()+
  xlab("Slope")+ylab("count log10")+
  ggtitle("ORFs with different signals at N, C or both N-C terminal")+
  theme_classic()

ggplot(combined_plates_final_signals_SWAT_SP_MTS, aes(x=Slope_lm)) +
  facet_wrap(~Signal_Protein_Terminal_Summary, ncol=1)+
  geom_histogram(binwidth=.01, bins=1000, colour="black", fill="black") +
  geom_vline(xintercept=Mean_Slope, size=1.5, color="red", linetype="dotted")+
  scale_y_log10()+
  xlab("Slope")+ylab("count log10")+
  ggtitle("ORFs with different signals at N, C or both N-C terminal")+
  theme_classic()

###Localization signals analysis count
combined_plates_final_signals_SWAT_SP_MTS<-read.csv("combined_plates_final_signals_SWAT_SP_MTS.csv")
combined_plates_final_signals_SWAT_SP_MTS_essentials<-subset(combined_plates_final_signals_SWAT_SP_MTS,!(ORF_essential=="non-essential"))
combined_plates_final_signals_SWAT_SP_MTS_nonessentials<-subset(combined_plates_final_signals_SWAT_SP_MTS,!(ORF_essential=="essential"))

#General count of protein localization signals for all 5247 ORFs
count_occurrences <- function(combined_plates_final_signals_SWAT_SP_MTS, ORF_essential, essential, Signal_Protein_Terminal_Summary) {
  subset_data <- combined_plates_final_signals_SWAT_SP_MTS[combined_plates_final_signals_SWAT_SP_MTS[[ORF_essential]] == essential, ]
  count_N <- nrow(subset_data[subset_data[[Signal_Protein_Terminal_Summary]] == 'N', ])
  count_C <- nrow(subset_data[subset_data[[Signal_Protein_Terminal_Summary]] == 'C', ])
  count_N_C <- nrow(subset_data[subset_data[[Signal_Protein_Terminal_Summary]] == 'N_C', ])
  count_no_signal <- nrow(subset_data[subset_data[[Signal_Protein_Terminal_Summary]] == 'no_signal', ])
  
  count_df <- data.frame(
    Protein_Essentiality = essential,
    Protein_Essentiality_Count = nrow(subset_data),
    N_Termini_signals = count_N,
    C_Termini_signals = count_C,
    N_and_C_Termini_signals = count_N_C,
    no_signals = count_no_signal
  )
  
  return(count_df)
}

# For non-essential
count_nonessential <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "non-essential", "Signal_Protein_Terminal_Summary")

# For essential
count_essential <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "essential", "Signal_Protein_Terminal_Summary")

# Combine the results
General_FigS1f <- rbind(count_nonessential, count_essential)

# Calculate percentages
General_FigS1f$N_Termini_signals_Percentage <- (General_FigS1f$N_Termini_signals / General_FigS1f$Protein_Essentiality_Count) * 100
General_FigS1f$C_Termini_signals_Percentage <- (General_FigS1f$C_Termini_signals / General_FigS1f$Protein_Essentiality_Count) * 100
General_FigS1f$N_and_C_Termini_signals_Percentage <- (General_FigS1f$N_and_C_Termini_signals / General_FigS1f$Protein_Essentiality_Count) * 100
General_FigS1f$no_signals_Percentage <- (General_FigS1f$no_signals / General_FigS1f$Protein_Essentiality_Count) * 100

#5% Threshold. Count of protein localization signals for all 5247 ORFs
count_occurrences <- function(combined_plates_final_signals_SWAT_SP_MTS, ORF_essential, essential, winner, N, Signal_Protein_Terminal_Summary) {
  subset_data <- combined_plates_final_signals_SWAT_SP_MTS[combined_plates_final_signals_SWAT_SP_MTS[[ORF_essential]] == essential & combined_plates_final_signals_SWAT_SP_MTS[[winner]] == N, ]
  print(subset_data)
  count_N <- sum(subset_data[[Signal_Protein_Terminal_Summary]] == 'N')
  count_C <- sum(subset_data[[Signal_Protein_Terminal_Summary]] == 'C')
  count_N_C <- sum(subset_data[[Signal_Protein_Terminal_Summary]] == 'N_C')
  count_no_signal <- sum(subset_data[[Signal_Protein_Terminal_Summary]] == 'no_signal')
  
  count_df <- data.frame(
    Protein_Essentiality = essential,  # Renamed to keep the desired column name
    Winner_5 = N,
    Protein_Essentiality_Count = nrow(subset_data),
    N_Termini_signals = count_N,
    C_Termini_signals = count_C,
    N_and_C_Termini_signals = count_N_C,
    no_signals = count_no_signal,
    stringsAsFactors = FALSE
  )
  
  count_df$N_Termini_signals_Percentage <- (count_df$N_Termini_signals / count_df$Protein_Essentiality_Count) * 100
  count_df$C_Termini_signals_Percentage <- (count_df$C_Termini_signals / count_df$Protein_Essentiality_Count) * 100
  count_df$N_and_C_Termini_signals_Percentage <- (count_df$N_and_C_Termini_signals / count_df$Protein_Essentiality_Count) * 100
  count_df$no_signals_Percentage <- (count_df$no_signals / count_df$Protein_Essentiality_Count) * 100
  
  return(count_df)
}

#For essential and N-winner
count_essential_N <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "essential", "Winner_5", "N", "Signal_Protein_Terminal_Summary")

#For non-essential and N-winner
count_nonessential_N <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "non-essential", "Winner_5", "N", "Signal_Protein_Terminal_Summary")

#For essential and C-winner
count_essential_C <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "essential", "Winner_5", "C", "Signal_Protein_Terminal_Summary")

#For non-essential and C-winner
count_nonessential_C <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "non-essential", "Winner_5", "C", "Signal_Protein_Terminal_Summary")

#For essential and neutral
count_essential_neutral <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "essential", "Winner_5", "neutral", "Signal_Protein_Terminal_Summary")

#For non-essential and neutral
count_nonessential_neutral <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "non-essential", "Winner_5", "neutral", "Signal_Protein_Terminal_Summary")

# Combine the results
General_FigS1f_5 <- rbind(count_essential_N, count_nonessential_N,count_essential_C,count_nonessential_C,count_essential_neutral,count_nonessential_neutral )
General_FigS1f_5$Winner_5[General_FigS1f_5$Winner_5=="C"]<- 'N-loser'
General_FigS1f_5$Winner_5[General_FigS1f_5$Winner_5=="N"]<- 'C-loser'
write.csv(General_FigS1f_5, paste0(path_input_dataset,"General_FigS1f_5.csv"))
General_FigS1f_5<-read.csv("General_FigS1f_5.csv")

#PLOT
# Prepare data for plotting
prepare_data <- function(data, protein_type) {
  subset_data <- subset(data, Protein_Essentiality == protein_type)
  categories <- c("N_Termini_signals_Percentage", "C_Termini_signals_Percentage", "N_and_C_Termini_signals_Percentage", "no_signals_Percentage")
  percentages <- subset_data[, categories]
  percentages$Category <- rownames(percentages)
  percentages$Protein_Type <- protein_type
  percentages$Winner_5 <- subset_data$Winner_5  # Add Winner_5 column
  return(percentages)
}

# Prepare data for essential and non-essential rows
essential_data <- prepare_data(General_FigS1f_5, "essential")
nonessential_data <- prepare_data(General_FigS1f_5, "non-essential")

# Combine the data
combined_data <- rbind(essential_data, nonessential_data)

# Reshape data into long format
combined_data_long <- gather(combined_data, key = "Category", value = "Percentage", -Category, -Protein_Type, -Winner_5)
combined_data_long_sub<-subset(combined_data_long,!combined_data_long$Category=="no_signals_Percentage")
write.csv(combined_data_long, paste0(path_input_dataset,"General_FigS1f_5_PlotTable.csv"))

# Plot using ggplot2
plotS1f_5<-ggplot(combined_data_long_sub, aes(x = Percentage, y = Winner_5, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Protein_Type, scales = "free") +
  labs(title = "Percentage of ORFs by Winner Type 5% threshold, abs(a) 0.175",
       x = "% ORFs", y = "Loser type") +
  scale_fill_manual(values = c("magenta", "cyan", "green", "red")) +
  theme_minimal()

cairo_pdf(filename = "FigS1f_alpha0.175.pdf", 
          width = 11, height = 8, onefile = TRUE)
plotS1f_5
dev.off()

#20% Threshold. Count of protein localization signals for all 5247 ORFs
count_occurrences <- function(combined_plates_final_signals_SWAT_SP_MTS, ORF_essential, essential, winner, N, Signal_Protein_Terminal_Summary) {
  subset_data <- combined_plates_final_signals_SWAT_SP_MTS[combined_plates_final_signals_SWAT_SP_MTS[[ORF_essential]] == essential & combined_plates_final_signals_SWAT_SP_MTS[[winner]] == N, ]
  print(subset_data)
  count_N <- sum(subset_data[[Signal_Protein_Terminal_Summary]] == 'N')
  count_C <- sum(subset_data[[Signal_Protein_Terminal_Summary]] == 'C')
  count_N_C <- sum(subset_data[[Signal_Protein_Terminal_Summary]] == 'N_C')
  count_no_signal <- sum(subset_data[[Signal_Protein_Terminal_Summary]] == 'no_signal')
  
  count_df <- data.frame(
    Protein_Essentiality = essential,  # Renamed to keep the desired column name
    Winner_20 = N,
    Protein_Essentiality_Count = nrow(subset_data),
    N_Termini_signals = count_N,
    C_Termini_signals = count_C,
    N_and_C_Termini_signals = count_N_C,
    no_signals = count_no_signal,
    stringsAsFactors = FALSE
  )
  
  count_df$N_Termini_signals_Percentage <- (count_df$N_Termini_signals / count_df$Protein_Essentiality_Count) * 100
  count_df$C_Termini_signals_Percentage <- (count_df$C_Termini_signals / count_df$Protein_Essentiality_Count) * 100
  count_df$N_and_C_Termini_signals_Percentage <- (count_df$N_and_C_Termini_signals / count_df$Protein_Essentiality_Count) * 100
  count_df$no_signals_Percentage <- (count_df$no_signals / count_df$Protein_Essentiality_Count) * 100
  
  return(count_df)
}

#For essential and N-winner
count_essential_N <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "essential", "Winner_20", "N", "Signal_Protein_Terminal_Summary")

#For non-essential and N-winner
count_nonessential_N <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "non-essential", "Winner_20", "N", "Signal_Protein_Terminal_Summary")

#For essential and C-winner
count_essential_C <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "essential", "Winner_20", "C", "Signal_Protein_Terminal_Summary")

#For non-essential and C-winner
count_nonessential_C <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "non-essential", "Winner_20", "C", "Signal_Protein_Terminal_Summary")

#For essential and neutral
count_essential_neutral <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "essential", "Winner_20", "neutral", "Signal_Protein_Terminal_Summary")

#For non-essential and neutral
count_nonessential_neutral <- count_occurrences(combined_plates_final_signals_SWAT_SP_MTS, "ORF_essential", "non-essential", "Winner_20", "neutral", "Signal_Protein_Terminal_Summary")

# Combine the results
General_FigS1f_20 <- rbind(count_essential_N, count_nonessential_N,count_essential_C,count_nonessential_C,count_essential_neutral,count_nonessential_neutral )
General_FigS1f_20$Winner_20[General_FigS1f_20$Winner_20=="C"]<- 'N-loser'
General_FigS1f_20$Winner_20[General_FigS1f_20$Winner_20=="N"]<- 'C-loser'
write.csv(General_FigS1f_20, paste0(path_input_dataset,"General_FigS1f_20.csv"))

#PLOT
# Prepare data for plotting
prepare_data <- function(data, protein_type) {
  subset_data <- subset(data, Protein_Essentiality == protein_type)
  categories <- c("N_Termini_signals_Percentage", "C_Termini_signals_Percentage", "N_and_C_Termini_signals_Percentage", "no_signals_Percentage")
  percentages <- subset_data[, categories]
  percentages$Category <- rownames(percentages)
  percentages$Protein_Type <- protein_type
  percentages$Winner_20 <- subset_data$Winner_20  # Add Winner_20 column
  return(percentages)
}

# Prepare data for essential and non-essential rows
essential_data <- prepare_data(General_FigS1f_20, "essential")
nonessential_data <- prepare_data(General_FigS1f_20, "non-essential")

# Combine the data
combined_data <- rbind(essential_data, nonessential_data)

# Reshape data into long format
combined_data_long <- gather(combined_data, key = "Category", value = "Percentage", -Category, -Protein_Type, -Winner_20)
combined_data_long_sub<-subset(combined_data_long,!combined_data_long$Category=="no_signals_Percentage")

# Plot using ggplot2
plotS1f_20<-ggplot(combined_data_long_sub, aes(x = Percentage, y = Winner_20, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Protein_Type, scales = "free") +
  labs(title = "Percentage of ORFs by Winner Type 20% threshold",
       x = "% ORFs", y = "Loser type") +
  scale_fill_manual(values = c("magenta", "cyan", "green", "red")) +
  theme_minimal()

cairo_pdf(filename = "FigS1f_alpha0.611.pdf", 
          width = 11, height = 8, onefile = TRUE)
plotS1f_20
dev.off()
