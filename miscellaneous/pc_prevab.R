# prevab_pc: a function for extracting scores from comp1 of a PC analysis carried
# out on a prevab table, as a measure of inmportance

prevab_pc <- function(prevab_df, verbose = verbose_output){
  prevab_mat <- prevab_df |>
    select(label, relAbundance, relprev) |>
    column_to_rownames("label") |>
    as.matrix() |>
    scale()
  prevabpc <- princomp(prevab_mat)
  # if(verbose){
  #  summary(prevabpc)
  #  plot(prevabpc)
  #  biplot(prevabpc)
  #  cat("The variance explained by components of the relprev/relab matrix is: ",
  #      prevabpc$sdev^2/sum(prevabpc$sdev^2))
  #  }
  scores_df <- as.data.frame(prevabpc$scores)
  return(scores_df)
}
