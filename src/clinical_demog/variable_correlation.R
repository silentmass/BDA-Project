install.packages("corrplot")

##checking the variable correlations

data <- readRDS("data/clinical_demog_clean.rds")

selected_variables <- c(
  "FALLER",
  "GDS",
  "AGE",
  "GCS_NEUROTRAX",
  "EFI_EXEC_FUNC_INDEX",
  "GENDER",
  "ABC_TOTAL_PERCENT",
  "SF36",
  "MMSE",
  "MOCA",
  "FAB",
  "TUG",
  "FSST",
  "BERG",
  "DGI",
  "TMT_A",
  "TMT_B",
  "BASE_VELOCITY",
  "S3_VELOCITY",
  "PASE",
  "FEET_CLOSE_EYES_OPEN",
  "FEET_CLOSE_EYES_CLOSED",
  "TANDEM_EYES_OPEN",
  "TANDEM_EYES_CLOSED"
)
new_data <- data
new_data[, selected_variables] <- lapply(data[, selected_variables], scale)
correlation_matrix <- cor(new_data[, selected_variables])
library(corrplot)
corrplot(correlation_matrix, method = "circle", type = "lower")
