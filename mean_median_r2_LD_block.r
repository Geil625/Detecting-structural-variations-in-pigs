rm(list=ls())
pig=... # the pig individual under consideration

setwd(paste0('.../',pig))
system('mkdir temp')
blocks <- readLines('LD.blocks')
res <- c()
for (i in 1:length(blocks))
{
  snps <- blocks[i]
  snps <- unlist(strsplit(snps,split = ' '))
  snps <- snps[-1]
  write.table(snps,file = paste0('temp/snps_',pig,'_',i,'.txt'),quote = F,row.names = F,col.names = F)
  system(paste0('plink --bfile ',pig,' --extract temp/snps_',pig,'_',i,'.txt --make-bed --out temp/',pig,'_',i))
  system(paste0('plink --bfile temp/',pig,'_',i,' --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out temp/',pig,'_',i))
  r2 <- read.table(paste0('temp/',pig,'_',i,'.ld'),header = T,as.is = T)[,7]
  ave <- mean(r2,na.rm=T)
  med <- median(r2,na.rm=T)
  snps <- paste(snps,collapse = '|')
  system(paste0('rm temp/*',pig,'_',i,'*'))
  res <- rbind(res,c(ave,med,snps))
}
colnames(res) <- c('average','median','snps')
write.table(res,file = 'blocks_r2.txt',quote = F,row.names = F,sep = '\t')
