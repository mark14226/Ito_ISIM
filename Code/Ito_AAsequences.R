# From data list, compile data frame of ALL information we have
TotalRep.df <- do.call(rbind.data.frame, Sample.aa.list)
# Determine the total number of detected amino acid sequences, thus the complete TCR repertoire across samples
length(unique(TotalRep.df$aminoAcid))

# Set coverage filter
# Set threshold
coverage = 5
# Consolidate 
FilterRep.df <- subset(TotalRep.df, TotalRep.df$count >= 5)







## Generate data frame of total counts across all samples (FOR CLONES WITH AT LEAST 1 SAMPLE WITH ENOUGH COVERAGE)
TotCountDF <- Sample.aa.list[[1]] %>% 
  dplyr::select(aminoAcid, count) %>%
  dplyr::mutate(Sample = names(Sample.aa.list[1]))
# Loop across samples and add each to data frame
for(i in 2:length(Sample.aa.list)){
    counts <- Sample.aa.list[[i]] %>% 
    dplyr::select(aminoAcid, count) %>%
    dplyr::mutate(Sample = names(Sample.aa.list[i]))
    TotCountDF <- merge(TotCountDF, counts, all = TRUE)
  }

# Generate wide version of data frame 
TotCountDF.wide <- dcast(TotCountDF, aminoAcid ~ Sample, value.var = "count")
TotCountDF.wide[is.na(TotCountDF.wide)] <- 0
rownames(TotCountDF.wide) <- TotCountDF.wide$aminoAcid
TotCountDF.wide$aminoAcid <- NULL
# Save to file
write.csv(TotCountDF.wide, paste0(mount, wd, "/Results/TotalRep_CountTable.csv"))


# Filter for aa with certain coverage in at least 1 sample
# Set threshold
coverage = 5
# Consolidate 
FilterCountDF.wide <- TotCountDF.wide[apply(TotCountDF.wide, 1, function(x) { max(x) >= 5 } ) ,]

















## Generate data frame of total frequencies across all samples (FOR CLONES WITH AT LEAST 1 SAMPLE WITH ENOUGH COVERAGE)
TotFreqDF <- Sample.aa.list[[1]] %>% 
  dplyr::select(aminoAcid, frequencyCount) %>%
  dplyr::mutate(Sample = names(Sample.aa.list[1]))
# Loop across samples and add each to data frame
for(i in 2:length(Sample.aa.list)){
  frequencies <- Sample.aa.list[[i]] %>% 
    dplyr::select(aminoAcid, frequencyCount) %>%
    dplyr::mutate(Sample = names(Sample.aa.list[i]))
  TotFreqDF <- merge(TotFreqDF, frequencies, all = TRUE)
}

# Generate wide version of data frame 
TotFreqDF.wide <- dcast(TotFreqDF, aminoAcid ~ Sample, value.var = "frequencyCount")
TotFreqDF.wide[is.na(TotFreqDF.wide)] <- 0
rownames(TotFreqDF.wide) <- TotFreqDF.wide$aminoAcid
TotFreqDF.wide$aminoAcid <- NULL
# Save to file
write.csv(TotFreqDF.wide, paste0(mount, wd, "/Results/TotalRep_FrequencyTable.csv"))


# Filter for aa with certain coverage in at least 1 sample
# Set threshold
MinFreq = 0.01
# Consolidate 
FilterFreqDF.wide <- TotFreqDF.wide[apply(TotFreqDF.wide, 1, function(x) { max(x) >= MinFreq } ) ,]
