library(dplyr)
library(ggplot2)
library(tidyr)

# Read combined data
df <- readRDS("~/BDA Project/data/clinical_demog_data_combined.rds")

# Create mapping of old to new column names
col_name_mapping <- list(
  "#" = "ID",
  "Date of Evaluation" = "DATE_OF_EVALUATION",
  "GCS (Neurotrax)" = "GCS_NEUROTRAX",
  "EFI (Exe. Func. Index)" = "EFI_EXEC_FUNC_INDEX",
  "Gender(1-female, 0-male)" = "GENDER",
  "Age" = "AGE",
  "Year Fall" = "YEAR_FALL",
  "6 Months Fall" = "SIX_MONTHS_FALL",
  "yr almost" = "YEAR_ALMOST",
  "GDS" = "GDS",
  "ABC Tot %" = "ABC_TOTAL_PERCENT",
  "SF-36" = "SF36",
  "PASE" = "PASE",
  "MMSE" = "MMSE",
  "MoCa" = "MOCA",
  "FAB" = "FAB",
  "TMTa" = "TMT_A",
  "TMTb" = "TMT_B",
  "TUG" = "TUG",
  "FSST" = "FSST",
  "BERG" = "BERG",
  "DGI" = "DGI",
  "DGI stairs" = "DGI_STAIRS",
  "base(velocity)" = "BASE_VELOCITY",
  "s3(velocity)" = "S3_VELOCITY",
  "feet close eyes open" = "FEET_CLOSE_EYES_OPEN",
  "feet close eyes closed" = "FEET_CLOSE_EYES_CLOSED",
  "tandem_eyes_open" = "TANDEM_EYES_OPEN",
  "tandem_eyes_closed" = "TANDEM_EYES_CLOSED",
  "FALLER" = "FALLER"
)

# Function to clean column names
clean_column_names <- function(df, mapping) {
  # Rename columns based on mapping
  names(df) <- sapply(names(df), function(x) {
    if (x %in% names(mapping)) {
      mapping[[x]]
    } else {
      x
    }
  })
  return(df)
}

# Example usage:
df <- clean_column_names(df, col_name_mapping)

# Change all versions of NA values to NA value
df <- df %>%
  mutate(across(-1, ~ replace(., . %in% c(
    "N/A", "n/a", "NA", "", "na"
  ), NA)))


# An example how to approach column type conversion

# An example how to convert column types
convert_special_numeric <- function(df, col_name, value_mappings = list("at least 2" = 2)) {
  # Convert column name to symbol for tidy evaluation
  col_sym <- sym(col_name)
  
  # First convert text values to numeric where possible
  df_converted <- df
  df_converted[[col_name]] <- as.character(df_converted[[col_name]])  # Convert to character first
  
  # Apply mappings if they exist
  if (length(value_mappings) > 0) {
    for (key in names(value_mappings)) {
      df_converted[[col_name]][df_converted[[col_name]] == key] <- value_mappings[[key]]
    }
  }
  
  df_converted[[col_name]] <- gsub(",", ".", df_converted[[col_name]])
  
  # Convert to numeric and handle NA values safely
  numeric_values <- suppressWarnings(as.numeric(df_converted[[col_name]]))
  na_rows <- is.na(numeric_values)
  
  # Calculate mean of non-NA values
  non_na_values <- numeric_values[!na_rows]
  non_na_col_mean <- if(length(non_na_values) > 0) mean(non_na_values, na.rm = TRUE) else NA
  
  print(paste0("COLUMN: ", col_name, " | NA REPLACED: ", sum(na_rows), " | NON NA MEAN: ", non_na_col_mean))
  
  na_stats <- data.frame(
    column = col_name,
    na_count = sum(na_rows),
    non_na_mean = non_na_col_mean,
    stringsAsFactors = FALSE
  )
  
  # Replace NA values with mean
  if(!is.na(non_na_col_mean)) {
    numeric_values[na_rows] <- non_na_col_mean
  }
  df_converted[[col_name]] <- numeric_values
  
  return(list(
    data = df_converted,
    na_stats = na_stats
  ))
}

process_columns <- function(df, column_mappings) {
  result <- df
  all_na_stats <- data.frame(
    column = character(),
    na_count = numeric(),
    non_na_mean = numeric(),
    stringsAsFactors = FALSE
  )
  
  # First process columns with mappings
  for (col_name in names(column_mappings)) {
    if (col_name %in% names(result)) {
      processed <- convert_special_numeric(result, col_name, column_mappings[[col_name]])
      result <- processed$data
      all_na_stats <- rbind(all_na_stats, processed$na_stats)
      message(sprintf("Processed mapped column: %s", col_name))
    }
  }
  
  # Find columns with NA values
  cols_with_na <- names(result)[sapply(result, function(x) any(is.na(x)))]
  
  # Process remaining columns that have NAs and can be converted to numeric
  remaining_cols <- setdiff(cols_with_na, names(column_mappings))
  
  for (col_name in remaining_cols) {
    processed <- convert_special_numeric(result, col_name, list())
    result <- processed$data
    all_na_stats <- rbind(all_na_stats, processed$na_stats)
    message(sprintf("Processed non mapped column: %s", col_name))
  }
  
  # Process non NA columns
  non_na_cols <- setdiff(names(df), union(remaining_cols, names(column_mappings)))
  
  for (col_name in non_na_cols) {
    if (!col_name %in% c("ID", "DATE_OF_EVALUATION")) {  # Changed from list() to c()
      result[[col_name]] <- as.numeric(gsub(",", ".", result[[col_name]]))
      message(sprintf("Processed non nan column: %s", col_name))
    }
  }
  
  return(list(
    data = result,
    na_stats = all_na_stats
  ))
}

# Example usage:
# df_processed <- process_columns(df, column_mappings)

# Example usage:
# df_processed <- process_columns(df, column_mappings)

# List column types
# sapply(df, class)

col_names <- colnames(df)

# col_names
# 1                      ID
# 2      DATE_OF_EVALUATION
# 3           GCS_NEUROTRAX
# 4     EFI_EXEC_FUNC_INDEX
# 5                  GENDER
# 6                     AGE
# 7               YEAR_FALL
# 8         SIX_MONTHS_FALL
# 9             YEAR_ALMOST
# 10                    GDS
# 11      ABC_TOTAL_PERCENT
# 12                   SF36
# 13                   PASE
# 14                   MMSE
# 15                   MOCA
# 16                    FAB
# 17                  TMT_A
# 18                  TMT_B
# 19                    TUG
# 20                   FSST
# 21                   BERG
# 22                    DGI
# 23             DGI_STAIRS
# 24          BASE_VELOCITY
# 25            S3_VELOCITY
# 26   FEET_CLOSE_EYES_OPEN
# 27 FEET_CLOSE_EYES_CLOSED
# 28       TANDEM_EYES_OPEN
# 29     TANDEM_EYES_CLOSED
# 30                 FALLER

# List column unique values
print(unique(df[sym("FALLER")]), n = dim(df)[1])

col_name <- "Gender(1-female, 0-male)"

# Example usage:
column_mappings <- list(
  "YEAR_FALL" = list(
    "at least 2" = 2
  ),
  "SIX_MONTHS_FALL" = list(
    "yes, but can't tell" = 1
  ),
  "YEAR_ALMOST" = list(
    "at least 1 per week" = 1,
    "more than once" = 1,
    "don't remember" = 0,
    "can't tell" = 0,
    "yes, no number" = 1,
    "more than 12" = 12,
    "once a week" = 52,
    "once a month" = 12,
    "2-3 a week" = 52 * 2.5,
    "many times" = 52
  )
)
df_processed <- process_columns(df, column_mappings)

saveRDS(df_processed$data, "data/clinical_demog_clean.rds")
saveRDS(df_processed$na_stats, "data/clinical_demog_clean_na_stats.rds")

# Calculate and visualise NA value counts

na_counts_col <- data.frame(
  column = names(colSums(is.na(df))),
  count = colSums(is.na(df))
)

na_counts_row <- data.frame(
  ID = df$ID,
  count = rowSums(is.na(df))
)

ggplot(na_counts_col, aes(x = count, y = reorder(column, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Missing Values by Column",
       x = "Count of NA Values",
       y = "Column Name") +
  theme(axis.text.y = element_text(size = 8))

ggplot(na_counts_row[na_counts_row$count > 0,], aes(x = count, y = reorder(ID, count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Missing Values by Subject",
       x = "Count of NA Values",
       y = "SUBJECT") +
  theme(axis.text.y = element_text(size = 8))