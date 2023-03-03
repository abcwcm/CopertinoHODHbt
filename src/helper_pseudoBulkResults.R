summarizeTestsPerLabel <- function (results, ...) 
{
    if (!is.matrix(results)) {
        results <- decideTestsPerLabel(results, ...)
    }
    output <- list()
    available <- sort(unique(as.vector(results)), na.last = TRUE)
    for (i in available) {
        output[[as.character(i)]] <- if (is.na(i)) {
            colSums(is.na(results))
        }
        else {
            colSums(results == i, na.rm = TRUE)
        }
    }
    do.call(cbind, output)
}
