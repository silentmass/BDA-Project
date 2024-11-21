# Analysis helper functions

# Function to check convergence metrics
check_convergence <- function(model) {
  # Get summary
  summary <- summary(model)
  
  # Extract Rhat values
  rhat <- summary$fixed[, "Rhat"]
  
  # Extract bulk and tail ESS
  bulk_ess <- summary$fixed[, "Bulk_ESS"]
  tail_ess <- summary$fixed[, "Tail_ESS"]
  
  # Print diagnostics
  cat("\nConvergence Diagnostics:\n")
  cat("----------------------\n")
  cat("Rhat values (should be close to 1):\n")
  print(rhat)
  cat("\nBulk Effective Sample Sizes:\n")
  print(bulk_ess)
  cat("\nTail Effective Sample Sizes:\n")
  print(tail_ess)
  
  # Check for convergence issues
  rhat_issues <- any(rhat > 1.05, na.rm = TRUE)
  ess_issues <- any(bulk_ess < 1000 | tail_ess < 1000, na.rm = TRUE)
  
  if(rhat_issues) {
    cat("\nWarning: Some Rhat values are > 1.05, indicating potential convergence issues.\n")
  }
  if(ess_issues) {
    cat("\nWarning: Some parameters have low effective sample sizes (< 1000).\n")
  }
  if(!rhat_issues && !ess_issues) {
    cat("\nNo major convergence issues detected.\n")
  }
}

evaluate_fall_model <- function(model, data, threshold = 0.5) {
  # Get predictions
  predictions <- posterior_predict(model)
  pred_probs <- colMeans(predictions)
  
  # ROC curve
  roc_obj <- roc(data$FALLER, pred_probs)
  roc_plot <- ggroc(roc_obj) +
    labs(title = "ROC Curve for Fall Classification",
         subtitle = paste("AUC =", round(auc(roc_obj), 3))) +
    theme_minimal()
  
  # Confusion matrix
  pred_class <- ifelse(pred_probs > threshold, 1, 0)
  conf_mat <- table(Predicted = pred_class, Actual = data$FALLER)
  
  # Calculate metrics
  sensitivity <- conf_mat[2,2] / sum(conf_mat[,2])
  specificity <- conf_mat[1,1] / sum(conf_mat[,1])
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  
  # Return results as list
  return(list(
    roc_plot = roc_plot,
    confusion_matrix = conf_mat,
    metrics = list(
      sensitivity = round(sensitivity, 3),
      specificity = round(specificity, 3),
      accuracy = round(accuracy, 3),
      auc = round(auc(roc_obj), 3)
    )
  ))
}

print_fall_results <- function(results) {
  cat("\n=== Fall Classification Model Results ===\n\n")
  
  # Print metrics
  cat("Performance Metrics:\n")
  cat(sprintf("Sensitivity: %.3f\n", results$metrics$sensitivity))
  cat(sprintf("Specificity: %.3f\n", results$metrics$specificity))
  cat(sprintf("Accuracy: %.3f\n", results$metrics$accuracy))
  cat(sprintf("AUC: %.3f\n\n", results$metrics$auc))
  
  # Print confusion matrix
  cat("Confusion Matrix:\n")
  print(results$confusion_matrix)
  cat("\n")
  
  # Print ROC plot
  print(results$roc_plot)
}

format_variables <- function(vars) {
  # Split into groups of 4
  chunks <- split(vars, ceiling(seq_along(vars)/4))
  # Join variables in each chunk with commas, then join chunks with newlines
  paste(sapply(chunks, paste, collapse=", "), collapse="\n")
}