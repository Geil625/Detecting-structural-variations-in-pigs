### this script is to remove variants that could not parsed by SURVIVOR ("SVTYPE=NA") and rearrange/recalculate the SV length
rm(list=ls())
pigs=c('AW','BKS','BMA','LW','PTR')

for (pig in pigs)
{
  setwd('...')
  input <- file('.../SURVIVOR_output.vcf','r')
  output <- file('.../SURVIVOR_output_rearranged.vcf','w')
  
  while (TRUE)
  {
    oneline <- readLines(con = input,n = 1)
    if (length(oneline)==0) break
    if (grepl('#',oneline))
    {
      if (grepl('^#CHROM',oneline)) 
      {
        oneline <- unlist(strsplit(oneline,split = '\t'))
        oneline <- oneline[-11]
        oneline[10] <- 'SAMPLE'
        oneline <- paste(oneline,collapse = '\t')
      }
      writeLines(oneline,con = output)
      next
    }
    oneline <- unlist(strsplit(oneline,split = '\t'))
    oneline <- oneline[-11]
    info <- oneline[8]
    info <- unlist(strsplit(info,split = ';'))
    svtype <- info[grep('^SVTYPE=',info)]
    svtype <- substring(svtype,8)
    if (svtype != 'NA')
    {
      oneline[4] <- 'N'
      if (svtype == 'INS') {
        oneline[5] <- paste0('<',svtype,'>')
        endpos <- oneline[2]
        info[grep('^END=',info)] <- paste0('END=',endpos)
      } else if (svtype == 'TRA') {
        info <- info[-(grep('^SVLEN=',info))]
        endpos <- info[grep('^END=',info)]
        endpos <- substring(endpos,5)
        if (!grepl('[',oneline[5],fixed=T) & !grepl(']',oneline[5],fixed=T)) print (c(oneline[1],oneline[2]))
      } else {
        oneline[5] <- paste0('<',svtype,'>')
        endpos <- info[grep('^END=',info)]
        endpos <- substring(endpos,5)
        if (svtype == 'DEL') {
          svlen <- -abs(as.integer(oneline[2])-as.integer(endpos))
        } else {
          svlen <- abs(as.integer(oneline[2])-as.integer(endpos))
        }
        info[grep('^SVLEN=',info)] <- paste0('SVLEN=',svlen)
      }
      oneline[3] <- paste(oneline[1],oneline[2],endpos,svtype,sep = '_')
      info <- paste(info,collapse = ';')
      oneline[8] <- info
      oneline <- paste(oneline,collapse = '\t')
      writeLines(oneline,con = output)
    }
  }
  close(input)
  close(output)
}
