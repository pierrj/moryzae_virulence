library(data.table)

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


# convert to binary
df[is.na(df)] <- 0

df[df != 0] <- 1

df[] <- lapply(df, function(x) as.numeric(x))


## add back OG names

df_ogs <- data.frame(fread(args[1], na.strings = ''))

rownames(df) <- df_ogs$Orthogroup

## change colnames to reflect genome names

colnames(df) <- paste(colnames(df), args[2], sep='')

## stack data frame and add back orthogroup information

df_stacked <- stack(df)

colnames(df_stacked) <- c('pav', 'genome')

df_stacked$og <- rep(rownames(df), ncol(df))

## subset to missing

df_stacked <- df_stacked[df_stacked$pav == 0,]

df_stacked <- subset(df_stacked, select = -c(pav))

## write to tsv

output.file <- file(args[3], "wb") ## to get unix line endings

write.table(df_stacked, file = output.file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

close(output.file)