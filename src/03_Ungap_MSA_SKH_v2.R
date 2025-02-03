pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }

pkgTest("rphast")


library(rphast)
args = commandArgs(trailingOnly=TRUE)
mymsa <- read.msa(args[1])
#ref = read.table(args[2],header=F)[,1]
ref = args[2]
Tripsicum_index <- grep(ref,mymsa[[2]])
out=args[3]
ungap_index <- which(unlist(strsplit(mymsa[[1]][Tripsicum_index], split=""))!="-")
missing_index <- c()
for(x in 1:length(mymsa[[1]]))
{
mymsa[[1]][x] <- paste0(unlist(strsplit(mymsa[[1]][x],split=""))[ungap_index],collapse="")
if(!(length(unlist(strsplit(mymsa[[1]][x],split="")))==sum(unlist(strsplit(mymsa[[1]][x],split=""))== "-")))
	{
	missing_index <- c(missing_index, x)
	}
}
print(missing_index)

##Remove the UTR copy of the reference
#missing_index <- missing_index[-2]
mymsa[[1]] <- mymsa[[1]][missing_index]
mymsa[[2]] <- mymsa[[2]][missing_index]
write.msa(mymsa,out,format="FASTA")
