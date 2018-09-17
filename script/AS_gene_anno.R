# set parameters
argv <- commandArgs(T)
anno_file <- argv[1]
AS_file <- argv[2]
output <- paste0(AS_file, '_gene_anno.tsv')

AS <- read.table(AS_file, header = T)
anno <- read.table(anno_file, sep='\t', header = T, quote = '', row.names = 1)
AS_table <- merge(AS, anno, by.x = 2, by.y = 0, all.x = T)

write.table(AS_table, output, row.names=F, sep='\t', quote=F)
