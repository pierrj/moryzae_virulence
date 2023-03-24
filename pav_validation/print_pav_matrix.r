library(data.table)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

## read in orthogroups

df <- data.frame(fread(args[1], na.strings = ''))

df <- subset(df, select = -c(Orthogroup))

## change names to match lineage info spreadsheet

str_split_vector <- function(x, string){
  output_vector <- c()
  for ( item in x ) {
    output <- strsplit(item, string)[[1]][1]
    output_vector <- c(output_vector, output)
  }
  return(output_vector)
}


# fix col names (assuming the first string is the genome name, breaks otherwise)
colnames(df) <- str_split_vector(colnames(df), '_')

## change colnames to reflect genome names

colnames(df) <- paste(colnames(df), args[2], sep='')

# convert to binary
df[is.na(df)] <- 0

df[df != 0] <- 1

df[] <- lapply(df, function(x) as.numeric(x))


## add back OG names

df_ogs <- data.frame(fread('Orthogroups.tsv', na.strings = ''))

rownames(df) <- df_ogs$Orthogroup


# read in validated PAVs

validated_pavs <- data.frame(fread(args[3], header = FALSE))

colnames(validated_pavs) <- c('genome', 'og', 'actuallyfound')

# keep only those that were thought to be pav and were verified to not be in genome
validated_pavs <- validated_pavs[validated_pavs$actuallyfound == 'no',]

# drop col
validated_pavs <- subset(validated_pavs, select = -c(actuallyfound))

# change to dataframe with only "actuallyfound: no" filled in
validated_pavs <- validated_pavs %>% 
            group_by(og) %>% 
            mutate(id=row_number()) %>% 
            pivot_wider(names_from=og, values_from=genome) %>% 
            ungroup()

validated_pavs <- subset(validated_pavs, select = -c(id))

# turn to df of 1s and 0s
validated_pavs <- validated_pavs %>%
                      pivot_longer(cols = everything(), names_to = 'Groups', 
                           values_drop_na = TRUE) %>%
                      distinct %>%
                      mutate(new =1) %>% 
                      pivot_wider(names_from =value, values_from = new,  
                             values_fill = list(new = 0))


## flip sign
validated_pavs[validated_pavs == 1] <- 2 # so you dont double switch

validated_pavs[validated_pavs == 0] <- 1

validated_pavs[validated_pavs == 2] <- 0


validated_pavs <- data.frame(validated_pavs)

rownames(validated_pavs) <- validated_pavs$Groups

validated_pavs <- subset(validated_pavs, select = -c(Groups))

# add missing orthogroups/rows

missing_ogs <- setdiff(rownames(df),rownames(validated_pavs))

rows_to_add <- data.frame(matrix(ncol = length(colnames(df)), nrow = length(missing_ogs)))

colnames(rows_to_add) <- colnames(df)
  
rownames(rows_to_add) <- missing_ogs

rows_to_add[is.na(rows_to_add)] <- 1

# bind together

validated_pavs <- rbind(validated_pavs, rows_to_add)

# sort

validated_pavs <- validated_pavs[order(rownames(validated_pavs)),]

validated_pavs <- validated_pavs[,order(colnames(validated_pavs))]

## add the two dataframes together

df <- as.data.frame(as.matrix(validated_pavs) + as.matrix(df))


# finally change 2s to 1s, 0 means a gene is really missing

df[df == 2] <- 1


## write pav matrix

output.file <- file(args[4], "wb") ## to get unix line endings

write.table(df, file = output.file, quote = FALSE, col.names = TRUE, row.names = TRUE)

close(output.file)