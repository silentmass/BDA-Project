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