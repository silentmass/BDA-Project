# Load required library
library(readr)
library(dplyr)

get_gender_col <- function(df) {
  grep("^Gender", names(df), value = TRUE)
}



# Read and combine the sheets
# Get sheet names
file_path_prefix <- "data/ClinicalDemogData_COFL/"

# Read and combine the sheets

controls <- read_delim(
  paste0(file_path_prefix, "Controls-Table 1.csv"),
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE,
)
controls <- controls %>%
  mutate_at(
    vars(
      "Year Fall",
      "6 Months Fall",
      "ABC Tot %",
      "FAB",
      "TMTa",
      "TMTb",
      "FSST",
      "feet close eyes open",
      "feet close eyes closed"
    ),
    as.character
  )
controls$FALLER <- 0

fallers <- read_delim(
  paste0(file_path_prefix, "Fallers-Table 1.csv"),
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

fallers$FALLER <- 1
fallers <- fallers %>%
  rename("Gender(1-female, 0-male)" = "Gender(0-male,1-female)")

variable_legend <- read_delim(
  paste0(file_path_prefix, "Legend-Table 1.csv"),
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)



combined_df <- bind_rows(controls, fallers)

col_names <- colnames(combined_df)

# saveRDS(combined_df, "data/clinical_demog_data_combined.rds")
# saveRDS(variable_legend, "data/clinical_demog_variable_legend.rds")