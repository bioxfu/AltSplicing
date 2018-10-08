#VS <- c('mdg10nramp1_Mn_vs_nramp1_Mn', 'mdg10nramp1_vs_mdg10nramp1_Mn', 'mdg10nramp1_vs_nramp1', 'nramp1_vs_nramp1_Mn')
VS <- commandArgs(T)

dfms <- list()

for (i in 1:length(VS)) {
  dfms[[i]] <- read.table(paste0('tables/', VS[i], '_event_compare.out'), header = T, sep = '\t', quote = '')
}

if (length(VS) == 4) {
  m1 <- merge(dfms[[1]], dfms[[2]], by.x = 1, by.y = 1, all = T)
  m2 <- merge(m1, dfms[[3]], by.x = 1, by.y = 1, all = T)
  m3 <- merge(m2, dfms[[4]], by.x = 1, by.y = 1, all = T)
  m <- m3[c(1:7,11:13, 17:19, 23:25)]
  colnames(m) <- sub('\\..+', '', colnames(m))
  colnames(m)[-c(1:4)] <- paste0(rep(VS, each=3), '.', colnames(m)[-c(1:4)])
  m <- m[m[,6] < 0.05 | m[,9] < 0.05 | m[,12] < 0.05 | m[,15] < 0.05, ]
}

write.table(m, 'tables/merge_AS_out', row.names = F, sep = '\t', quote = F)
