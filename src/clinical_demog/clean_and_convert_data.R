library(dplyr)
library(ggplot2)
library(tidyr)

# Read combined data
df <- readRDS("~/BDA Project/data/clinical_demog_data_combined.rds")

# Rename # column to ID
df <- rename(df, "ID" = "#")

# Change all versions of NA values to NA value
df <- df %>%
  mutate(across(-1, ~ replace(., . %in% c(
    "N/A", "n/a", "NA", ""
  ), NA)))


# An example how to approach column type conversion

# An example how to convert column types
convert_special_numeric <- function(df, col_name, value_mappings = list("at least 2" = "2")) {
  # Convert column name to symbol for tidy evaluation
  col_sym <- sym(col_name)
  
  # Create case_when statements dynamically from the mappings
  cases <- lapply(names(value_mappings), function(key) {
    quo(!!sym(col_name) == !!key ~ !!value_mappings[[key]])
  })
  
  # Add the default case
  cases[[length(cases) + 1]] <- quo(TRUE ~ !!sym(col_name))
  
  df[!is.na(df[[col_name]]), ] %>%
    mutate(
      !!col_sym := case_when(!!!cases)
    ) %>%
    mutate(
      !!col_sym := as.numeric(!!col_sym)
    )
}

# List column types
sapply(df, class)

col_name <- "Year Fall"

# List column unique values
unique(df[sym(col_name)])

# Example usage:
mappings <- list(
  "at least 2" = "2",
  "more than 3" = "3",
  "about 5" = "5"
)
result <- convert_special_numeric(df, col_name, mappings)


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