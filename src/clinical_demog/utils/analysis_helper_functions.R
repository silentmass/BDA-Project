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