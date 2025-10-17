# Function to compute contingency table and chi-square test for two variables
contingency_chisq <- function(x, y, id) {
    # Create contingency table
    contingency_table <- table(c(x, TRUE, FALSE, FALSE, TRUE), c(y, TRUE, TRUE, FALSE, FALSE)) - 1

    # Perform chi-square test (with error handling for small tables)
    chisq_result <- tryCatch(
        {
            test <- chisq.test(contingency_table, correct = FALSE)
            list(
                id = id,
                contingency_table = contingency_table,
                chisq_statistic = test$statistic,
                p_value = test$p.value
            )
        },
        error = function(e) {
            list(
                id = id,
                contingency_table = contingency_table,
                chisq_statistic = NA_real_,
                p_value = NA_real_
            )
        }
    )
    return(chisq_result)
}
