# scratch.R

lData.2 = lapply(lData, function(df) {
  df = df[is.finite(df$log2.fold_change) & df$status == 'OK',]
  fc = df$log2.fold_change
  fc[fc > 5] = 5
  fc[fc < -5] = -5
  df$log2.fold_change. = fc
  df2 = data.frame(gene_id=as.character(df$test_id), log2.fold_change=df$log2.fold_change., p_value=df$p_value, q_value = df$q_value)
  rownames(df2) = NULL
  return(df2)
})

n = names(lData.2)

sapply(n, function(x){
  n2 = paste('Results/', x, '.csv', sep='')
  write.csv(lData.2[[x]], file=n2, row.names = F)
})
