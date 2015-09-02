# load initial data
file<-readData("FastaSeqsPau",include.unknown=T)

# useful stuff
get.sum.data(file)->summary
summary

# first we need to apply the neutrality stats
file<-neutrality.stats(file)

# then we apply F_ST stats
file<-F_ST.stats(file)

# we call Hudson program, thetaID needs to be defined for any module
ms.file1 <- MS(file, niter=100,thetaID="Tajima", neutrality=TRUE)
ms.file1 # we can see available slots
ms.file1@locus

ms.file2 <- MS(file, niter=100,thetaID="Tajima", F_ST=TRUE)
ms.file2
ms.file2@locus

# It appears that from both, the neutrality or F_ST simulations
# we should be able to recover number of haplotypes with
ms.file1@locus[[1]]@haplotypes 
# but the result is <0 x 0 matrix>

# if we try
ms.file1@locus[[1]]@haplotypes[1] 
# the result is [1] NA

# if we try 
ms.file1@locus[[1]]@haplotypes[[1]] 
# the result is Error in ms.file1@locus[[1]]@haplotypes[[1]] : subscript out of bounds
