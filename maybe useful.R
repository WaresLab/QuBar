ms.output <- function(filename = "", outfile = paste(filename, ".out", sep = ""), out.opt = "medium", do.amova = FALSE, graph = "none", splitpops = FALSE, lumppops = FALSE, excludepops = FALSE, ...){
  if(do.amova) library(ade4)	#make sure ade4 is loaded so we can use the amova function
  
  ## ms.header parses the header information from the ms output.  
  ms.header <- function(ms.outfile){
    header <- ms.outfile[1:(grep("//", ms.outfile)[1] - 1)]
    header.list <- c()
    if(any(header == "-I")) {		
      npops <<- as.numeric(header[grep("-I", header)[1] + 1])
      nseqs <<- as.numeric(header[grep("-I", header)[1] + (2:(npops + 1))])
      header.list <- c(header.list, npops = npops, nseqs = nseqs)
      if(any(header == "-g")) {  # removed as.numeric() Dec2011
        gps <- header[grep("-g", header) + 1]
        g <- header[grep("-g", header) + 2]
        header.list <- c(header.list, g = g)
      }
      
      if(any(header == "-m")){  
        m <- c()
        for(i in 1:length(grep("-m", header)))  {
          mlab <- paste("m", header[grep("-m", header) + 1][i], header[grep("-m", header) + 2][i], sep="")
          m[i] <- header[grep("-m", header) + 3][i]
          names(m)[i] <- mlab
        }
        header.list <- c(header.list, m)	
      }
      
      if(length(grep("-", header[which(header == "-I") + npops + 2])) == 0) {  # removed as.numeric() Dec2011
        fourNm <- header[which(header == "-I") + npops + 2]
        header.list <- c(header.list, fourNm = fourNm)
      }	
      
    } else {   #this is else to "if any header=I" so should mean there was no I in the header
      npops <<- 1
      nseqs <<- header[2]
      header.list <<- c(header.list, npops = npops, nseqs = nseqs)
    } 
    
    # The rest of the header should be parsed the same way whether I specified or not
    
    nsims <<- as.numeric(header[3])
    
    if(any(header == "-T")) ms.trees <<- TRUE
    
    if(any(header == "-s")) {  # removed as.numeric() Dec2011
      s <- header[grep("-s", header)[1] + 1]
      header.list <- c(header.list, s = s)
    }
    
    if(any(header == "-s") && any(header == "t")) st <<- TRUE   # because the documentation implies the output is different in this case, although I don't see it so far
    
    if(any(header == "-G")) {  # removed as.numeric() Dec2011
      G <<- header[grep("-G", header)[1] + 1] 
      header.list <- c(header.list, G = G)
    }
    
    if(any(header == "-n")){
      nps <- header[grep("-n", header) + 1]
      n <- header[grep("-n", header) + 2]
      names(n) <- paste("n", nps, sep="")
      header.list <- c(header.list, n)
    }
    
    if(any(header == "-r")) {   # removed as.numeric() Dec2011
      r <- header[which(header == "-r") + 1]
      nsites <<- header[which(header == "-r") + 2]
      header.list <- c(header.list, r = r)
    }
    if(any(header == "-c")) {  # removed as.numeric() Dec2011
      f <- header[which(header == "-c") + 1]
      L <- header[which(header == "-c") + 2]
      header.list <- c(header.list, f = f, L = L)
    }
    if(any(header == "-ma")){  # removed as.numeric() Dec2011
      ms <- header[which(header == "-ma") + (1:(npops^2))]  #generates warning messages for nonnumerics... not ideal, but not fatal
      ms <- matrix(ms, nrow = npops, ncol = npops, byrow = TRUE)
      header.list <- c(header.list, ms = ms)
    }
    if(length(grep("-e", header)) != 0){	#should parse -e headers correctly but there are a lot of options...
      if(any(header == "-eG")){  # removed as.numeric() Dec2011
        eGt <- header[which(header == "-eG") + 1]
        eGa <- header[which(header == "-eG") + 2]
        header.list <- c(header.list, eGt = eGt, eGa = eGa)
      }
      if(any(header == "-eg")){  # removed as.numeric() Dec2011
        egt <- header[which(header == "-eg") + 1]
        egi <- header[which(header == "-eg") + 2]
        ega <- header[which(header == "-eg") + 3]
        header.list <- c(header.list, egt = egt, egi = egi, ega = egi)
      }
      if(any(header == "-eN")){  # removed as.numeric() Dec2011
        eNt <- header[which(header == "-eN") + 1]
        eNx <- header[which(header == "-eN") + 2]
        header.list <- c(header.list, eNt = eNt, eNx = eNx)
      }
      if(any(header == "-en")){  # removed as.numeric() Dec2011  #tbs problem
        ent <- header[which(header == "-en") + 1]
        eni <- header[which(header == "-en") + 2]
        enx <- header[which(header == "-en") + 3]
        header.list <- c(header.list, ent = ent, eni = eni, enx = enx)
      }
      if(any(header == "-eM")){  # removed as.numeric() Dec2011
        eMt <- header[which(header == "-eM") + 1]
        eMx <- header[which(header == "-eM") + 2]
        header.list <- c(header.list, eMt = eMt, eMx = eMx)
      }
      if(any(header == "-em")){  # removed as.numeric() Dec2011
        emt <- header[which(header == "-em") + 1]
        emi <- header[which(header == "-em") + 2]
        emj <- header[which(header == "-em") + 3]
        emx <- header[which(header == "-em") + 4]
        header.list <- c(header.list, emt = emt, emi = emi, emj = emj, emx = emx)
      }
      if(any(header == "-ema")){
        emat <- as.numeric(header[which(header == "-ema") + 1])
        emanpop <- as.numeric(header[which(header == "-ema") + 2])	
        emaM <- matrix(as.numeric(header[which(header == "-ema") + seq(3, length.out = emanpop^2)]), nrow = emanpop, ncol = emanpop, byrow = TRUE)
        header.list <- c(header.list, emat = emat, emanpop = emanpop, emaM = emaM)
      }
      if(any(header == "-es")){  # removed as.numeric() Dec2011
        est <- header[which(header == "-es") + 1]
        esi <- header[which(header == "-es") + 2]
        esp <- header[which(header == "-es") + 3]
        header.list <- c(header.list, est = est, esi = esi, esp = esp)
      }
      if(any(header == "-ej")){  # removed as.numeric() Dec2011
        ejt <- header[which(header == "-ej") + 1]
        eji <- header[which(header == "-ej") + 2]
        ejj <- header[which(header == "-ej") + 3]
        header.list <- c(header.list, ejt = ejt, eji = eji, ejj = ejj)
      }
    }
    ## adjust the population numbers if you want to
    ## combine, split, or exclude some populations
    ## WHY ARE THESE IN THE HEADER FUNCTION?
    if(lumppops) {
      nseqs <<- lump.pops(npops, nseqs, ...)
      npops <<- length(nseqs)
    }
    
    if(splitpops) {
      nseqs <<- split.pops(npops, nseqs, ...)
      npops <<- length(nseqs)
    }
    ## Implementation of excludepops is slightly different than 
    ##  lumppops and splitpops.  Exclusion all happens here in
    ##  the header function, while lumping and splitting 
    ##  call functions that are defined separately.
    if(excludepops){
      epop <- as.numeric(strsplit(readline("Exclude Which population?"), split = ",")[[1]])[1]
      slashes <- grep("//", ms.outfile)
      ssites <- as.numeric(ms.outfile[grep("segsites:", ms.outfile) + 1])
      if(epop == npops){ #excluding last pop
        # remove population data
        ms.outfile <- ms.outfile[1:(length(ms.outfile) - 1)]
        for(i in (length(slashes) - 1):1){
          ms.outfile <- c(ms.outfile[1:(slashes[i] + 3 + ssites[i] + sum(nseqs[1:(epop - 1)]))], ms.outfile[(slashes[i] + 3 + ssites[i] + sum(nseqs) + 1):length(ms.outfile)])
        }
        nseqs <<- nseqs[1:(epop - 1)]
      }
      if(epop == 1){ #excluding first pop
        # remove population data
        for(i in length(slashes):1){
          ms.outfile <- c(ms.outfile[1:(slashes[i] + 3 + ssites[i])], ms.outfile[(slashes[i] + 3 + ssites[i] + nseqs[1] + 1):length(ms.outfile)])
        }
        nseqs <<- nseqs[2:length(nseqs)]
      }
      if(epop>1 && epop<npops){ #excluding interior pop
        # remove population data
        for(i in length(slashes):1){
          ms.outfile <- c(ms.outfile[1:(slashes[i] + 3 + ssites[i] + sum(nseqs[1:(epop - 1)]))], ms.outfile[(slashes[i] + 3 + ssites[i] + sum(nseqs[1:(epop)]) + 1):length(ms.outfile)])
        }
        nseqs <<- c(nseqs[1:(epop - 1)], nseqs[(epop + 1):length(nseqs)])
      }
      
      ms.in <<- ms.outfile
      npops <<- length(nseqs)
    }
    
    pops <<- c()
    pops[1] <<- 1
    if (npops>1){
      for (i in 2:npops){
        pops[i] <<- pops[i-1] + nseqs[i-1]
      }
    }
    
    return(header.list)
  }	
  ## end of the header parsing function
  
  
  
  ## lump.pops() and split.pops() are functions to combine or split populations.
  ## If lumping or splitting is specified on the command line they are called
  ##  from within the ms.header() function.
  lump.pops <- function(npops, nseqs, lpop = c(0), lumpno = c(0)){
    if(lpop[1] == 0) lpop <- as.numeric(strsplit(readline("First pops in each lump?"), split = ",")[[1]])
    lumps <- length(lpop)
    if(lumpno[1] == 0) lumpno <- as.numeric(strsplit(readline("Number of pops in each lump?"), split = ",")[[1]])
    lumpees <- c()
    for(i in 1:lumps) lumpees[[i]] <- c(lpop[i]:(lpop[i] + lumpno[i] - 1))
    for(i in 1:npops) if(i %in% lpop) nseqs[i] <- sum(nseqs[i:(i + lumpno[which(lpop == i)] - 1)])
    nseqs <- nseqs[which(1:npops %in% c(lpop, setdiff(1:npops, unlist(lumpees))))]
  }
  
  split.pops <- function(npops, nseqs, spop = c(0), newpops = c(0), newns = c(0)){
    if(spop[1] == 0) spop <-  as.numeric(strsplit(readline("Split which pops?"), split = ",")[[1]])
    splits <- length(spop)
    if(newpops[1] == 0) newpops <- as.numeric(strsplit(readline("Split into how many?"), split = ",")[[1]])
    if(newns[1] == 0) newns <- as.numeric(strsplit(readline("New pop ns?"), split = ",")[[1]])
    
    newns.list <- c()
    for(i in 1:splits) {
      if(i == 1) newns.list[[i]] <- newns[1:(newpops[i])]
      if(i>1) newns.list[[i]] <- newns[(sum(newpops[1:(i - 1)]) + 1):(sum(newpops[1:(i - 1)]) + newpops[i])]
    }	#newns.list should now be a list of all the new subpop ns
    
    
    newnpops <- npops + sum(newpops) - splits
    newnseqs <- c()
    for(i in 1:npops){
      if(i %in% spop) {
        newnseqs[i] <- newns.list[which(spop == i)]
      } else {
        newnseqs[i] <- nseqs[i]
      }
    }
    newnseqs <- unlist(newnseqs)	 
    return(newnseqs)
  }
  
  ##
  ## Functions to manipulate simulated data and calculate
  ##  various population genetic summary statistics.
  ##  Function definitions only here; functions will get 
  ##  called in the main script farther down.
  
  ## Distance matrices
  ## Comparisons with pop i in sim k are
  ##  absdist[[k]][pops[i]:(pops[i] + nseqs[i] - 1), pops[i]:(pops[i] + nseqs[i] - 1)]  
  ## Comparisons between pop i and pop j in sim k are
  ##  absdist[[k]][pops[i]:(pops[i] + nseqs[j] - 1), pops[j]:(pops[j] + nseqs[j] - 1)]  
  ##  or generalize; for i = j this is within a pop, and for i != j it's between pops
  
  
  dist.mat <- function(ms.sim, nseqs){
    if(is.na(ms.sim[1])) {
      a <- matrix(0, nrow = sum(nseqs), ncol = 1)
    } else {
      if (nchar(ms.sim[1])>1) {
        a <- t(sapply(strsplit(ms.sim, split = ""), as.matrix))
      } else {
        if (nchar(ms.sim[1]) == 1) {
          a <- t(t(sapply(strsplit(ms.sim, split = ""), as.matrix)))
        }
      }
    }
    absdist <- as.matrix(dist(a, method = "manhattan"))
    return(absdist)
  }
  
  sq.dist.mat <- function(ms.sim, nseqs){
    if(is.na(ms.sim[1])) {
      a <- matrix(0, nrow = sum(nseqs), ncol = 1)
    } else {
      if (nchar(ms.sim[1])>1) {
        a <- t(sapply(strsplit(ms.sim, split = ""), as.matrix))
      } else {
        if (nchar(ms.sim[1]) == 1) {
          a <- t(t(sapply(strsplit(ms.sim, split = ""), as.matrix)))
        }
      }
    }
    sqdist <- as.matrix(dist(a, method = "manhattan")^2)
    return(sqdist)
  }
  
  
  ht.dist.mat <- function(ms.sim, nseqs){
    if(is.na(ms.sim[1])) {
      a <- matrix(0, nrow = sum(nseqs), ncol = 1)
    } else {
      if (nchar(ms.sim[1])>1) {
        a <- t(sapply(strsplit(ms.sim, split = ""), as.matrix))
      } else {
        if (nchar(ms.sim[1]) == 1) {
          a <- t(t(sapply(strsplit(ms.sim, split = ""), as.matrix)))
        }
      }
    }
    u <- unique(a)
    htdist <- as.matrix(dist(u, method = "manhattan"))
    return(htdist)
  }
  
  ##
  ## Functions for single simulations
  ## These return arrays with nsims cols and npops rows (for pi) or npops x npops rows (for d)
  ##
  pi.sim <- function(distmat, nseqs, pops, npops){
    pi.pops <- c()
    d <- array(distmat, dim = c(sum(nseqs), sum(nseqs)))
    for (i in 1:npops){
      pi.pops[i] <- sum(d[pops[i]:(pops[i] + nseqs[i] - 1), pops[i]:(pops[i] + nseqs[i] - 1)])/(nseqs[i] * (nseqs[i] - 1))
    }
    return(pi.pops)
  }
  
  d.sim <- function(distmat, pops, nseqs, npops){
    d.pops <- array(dim = c(npops, npops))
    d <- array(distmat, dim = c(sum(nseqs), sum(nseqs)))
    for (i in 1:(npops - 1)){
      for (j in (i + 1):npops){
        d.pops[i,j] <- mean(d[pops[i]:(pops[i] + nseqs[i] - 1), pops[j]:(pops[j] + nseqs[j] - 1)])
      }
    }
    return(d.pops)
  }
  
  d.to.dist <- function(d.all.col, npops){
    d.dist <- as.dist(t(matrix(d.all.col, nrow = npops, ncol = npops)))
    return(d.dist)
  }
  
  phist.pair <- function(dxy, pix, piy) (dxy - (mean(c(pix, piy)))) / dxy
  
  phist.inds <- function(ind.v, d.mat, pi.vec) phist.pair(d.mat[ind.v[1],ind.v[2]], pi.vec[ind.v[1]], pi.vec[ind.v[2]])
  
  phist.pair.sim <- function(d.long.v, pi.v, npops){  ## I think this works but it's ugly
    phist.v <- c()
    d.mat <- matrix(d.long.v, nrow = npops, ncol = npops) ## upper triangular
    phist.v <- apply(combn(1:npops,2), 2, phist.inds, d.mat=d.mat, pi.vec=pi.v)
    names(phist.v) <- apply(combn(1:npops,2), 2, paste, sep="", collapse="")
    return(phist.v)
  }
  
  phist.pair.simno <- function(no, d.long.v, pi.v, npops) phist.pair.sim(d.long.v[,no], pi.v[,no], npops)
  
  ## Segregating sites
  pop.segs <- function(ms.sim, npops, nseqs, pops){
    segs <- array(dim = npops)
    if(is.na(ms.sim[1])) {
      segs[1:npops] <- 0
    } else {
      for (j in 1:npops){
        seg <- c()
        for(i in 1:nchar(ms.sim[1])){
          site <- sapply(ms.sim[pops[j]:(pops[j] + nseqs[j] - 1)], substr, i, i)
          seg[i] <- !all(site == site[1])   #will be FALSE if constant or TRUE if segregating; interpret as "does site i segregate in pop j?"
        }
        segs[j] <- length(which(seg))
      }
    }
    return(segs)
  }
  
  ## Find and count haplotypes
  ## for a single simulation ms.sim
  ## apply to rows(1) of sims; returns a list whose elements are arrays, all different sizes
  haplotypes <- function(ms.sim, npops, pops, nseqs){
    nhap <- length(unique(ms.sim))
    haps <- unique(ms.sim)
    ht <- array(dim = c(nhap, npops))
    for (j in 1:nhap){
      for (k in 1:npops){
        ht[j,k] <- length(which(ms.sim[(pops[k]:(pops[k] + nseqs[k] - 1))] == haps[j]))
      }
    }
    return(ht)
  }
  
  ht.all <- function(ht.mat){
    rowsums <- function(num.array) apply(num.array, 1, sum)
    sapply(ht.mat, rowsums)
  }
  
  nhaps <- function(ht.matrix){
    nhap.pop <- function(ht.matrix.column){
      nhap <- length(which(ht.matrix.column != 0))
    }
    hts.pop <- apply(ht.matrix, 2, nhap.pop)
  }
  
  hap.freqs <- function(ht.matrix, nseqs){
    ht.freq <- ht.matrix/nseqs
  }
  
  
  ## Tajima's D
  #using http://ocw.mit.edu/NR/rdonlyres/Health-Sciences-and-Technology/HST-508Fall-2005/15E2F713-CE88-4D7D-ACE1-D18004F26E40/0/tajimad1.pdf
  #that is, D = ( pi - S/a1 ) / sqrt( e1 * S + (e2 * S * (S - 1)))
  tajima.d <- function(pi.all, S.all, nseqs, pops, npops){
    
    a1 <- c()
    a2 <- c()
    for (i in 1:npops){
      a1[i] <- sum(1/(1:(nseqs[i] - 1)))
      a2[i] <- sum(1/((1:(nseqs[i] - 1))^2))
    }
    b1 <- (nseqs + 1)/(3 * (nseqs - 1))
    b2 <- (2 * (nseqs^2 + nseqs + 3))/(9 * nseqs * (nseqs - 1))
    c1 <- b1 - (1/a1)
    c2 <- b2 - ((nseqs + 2)/(a1*nseqs)) + (a2/(a1^2))
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    
    D <- (pi.all - (S.all/a1))/sqrt(e1 * S.all + e2 * S.all * (S.all - 1))
  }
  
  ## Haplotype diversity for a single haplotype matrix
  hap.div <- function(ht.freq.matrix, nseqs){
    if(sum(ht.freq.matrix) == 0) Hd.pops <- c(rep(0, times = ncol(ht.freq.matrix)))
    if(sum(ht.freq.matrix) > 0) Hd.pops <- (nseqs/(nseqs - 1)) * (1 - apply(ht.freq.matrix^2, 2, sum))
    return(Hd.pops)
  }
  
  ## Haplotype counts
  sumhaps <- function(ht.matrix){
    apply(ht.matrix, 1, sum)
  }
  
  ## Haplotype diversity for a list of haplotype matrices
  hap.div.all <- function(htm.list, nseqs){
    pophaps <- sapply(htm.list, sumhaps)
    div.by.nseqs <- function(ht.counts, nseqs) ht.counts/sum(nseqs)
    all.freqs <- sapply(pophaps, div.by.nseqs, sum(nseqs))
    square <- function(numbers) numbers^2
    squares <- sapply(all.freqs, square)
    Hd.list <- function(sq.freq.list, nseqs) {
      if(sum(sq.freq.list) == 0) h <- 0
      if(sum(sq.freq.list) != 0) h <- (sum(nseqs)/(sum(nseqs) - 1)) * (1 - sum(sq.freq.list))
      return(h)
    }
    hapdiv.all <- sapply(squares, Hd.list, nseqs)
    return(hapdiv.all)
  }
  
  ## AMOVA and phiST
  #should return NA for sims with no seg sites and an amova for everything else.  Apply to absdist dim 2 with samples = fake samples matrix
  ms.amova <- function(distances, samples, nseqs, struct = NULL){
    if(any(distances != 0)) amova(samples, sqrt(as.dist(matrix(distances, nrow = sum(nseqs), ncol = sum(nseqs)))), struct) else NA
  }
  
  get.phist <- function(amova){
    if(all(is.na(amova))) NA else phist.amova <- amova$statphi[[1]][1]
  }
  
  ##
  ## End of function definitions
  ##
  
  ##
  ## Now the main script!
  ##
  
  if(outfile != "") {
    #sink(file = outfile, split = TRUE)
    message(paste("Writing text output to", outfile))
  }
  message(paste("Output mode is", out.opt, "."))
  # if(graph != "screen" & graph != "none") message(paste("Writing graphical output to", graph))
  if(filename == "") filename <- readline("filename with output from ms? (just name, and path if necessary; no quotation marks)  ")
  if(filename == "") stop("exiting--no filename entered", call. = FALSE)
  ms.in <- scan(filename, what = "c", quiet = TRUE)
  message(paste("Reading", filename, "and calculating.... please be patient; I have a lot of pairwise comparisons to do!\n"))
  ## strip the tbs arguments after //
  slashes <- grep("//", ms.in)
  starting.spots <- grep("segsites:", ms.in)
  for(i in length(slashes):1){
    ms.in <- c(ms.in[1:slashes[i]], ms.in[starting.spots[i]:length(ms.in)])
  }
  headerlist <<- ms.header(ms.in)		## parses the header
  segsites <- as.numeric(ms.in[grep("segsites:", ms.in) + 1])  ## identifies segregating sites (from the output file)
  delims <- which(ms.in == "//")   ## finds the slashes in the output file that delimit simulations
  ## pull out the simulated data and store it in an array
  sims <- array(dim = c(nsims, sum(nseqs)))	
  #if(st)    ##documentation implies there's an extra output line in this case, but actual output file seems to have a blank line so the total is the same
  for (i in 1:nsims) {
    if(segsites[i] != 0) {
      for (j in 1:(sum(nseqs))) {
        sims[i,j] <- ms.in[(delims[i] + 3 + segsites[i] + 1):(delims[i] + 3 + segsites[i] + 1 + sum(nseqs) - 1)][j]
      }
    }
  }
  ## Calculate pop gen statistics on the simulated data
  absdist <<- apply(sims, 1, dist.mat, nseqs)
  pi.all <- apply(absdist, 2, pi.sim, nseqs, pops, npops)
  if(npops == 1) pi.all <- t(pi.all)
  pi.tot <- apply(absdist, 2, pi.sim, sum(nseqs), 1, 1)
  if(npops>1) d.all <<- apply(absdist, 2, d.sim, pops, nseqs, npops)
  if(npops>1) d.all.short <- d.all[which(!is.na(d.all[,1])),] 
  S.all <- apply(sims, 1, pop.segs, npops, nseqs, pops)
  if(npops == 1) S.all <- t(S.all)
  ht <<- apply(sims, 1, haplotypes, npops, pops, nseqs)
  ht.freqs <<- sapply(ht, hap.freqs, nseqs)
  ht.dists <<- apply(sims, 1, ht.dist.mat, nseqs) 
  ht.pops <<- sapply(ht, nhaps)
  if(npops == 1) ht.pops <- t(ht.pops)
  tajD.all <- tajima.d(pi.all, S.all, nseqs, pops, npops)
  if(npops == 1) tajD.all <- t(t(tajD.all))  ## What am I trying to accomplish here?
  Hd.pops <- sapply(ht.freqs, hap.div, nseqs)
  if(npops == 1) Hd.pops <- t(Hd.pops)
  Hd.all <- hap.div.all(ht, nseqs)
  #pi.among.all <<- apply(absdist, 2, pi.among, nseqs)
  if(npops > 1) phi.st <- (apply(d.all, 2, mean, na.rm = TRUE) - apply(pi.all, 2, mean))/(apply(d.all, 2, mean, na.rm = TRUE))
  ## Pairwise phist here  ## comes out as ncomps x nsims array, rownames are comparisons, no colnames
  if(npops > 2) phi.st.pairs <- sapply(1:nsims, phist.pair.simno, d.long.v=d.all, pi.v=pi.all, npops=npops)
  if(do.amova){
    amova.fake.ht <- array(dim = c(sum(nseqs), npops), data = 0)
    for (i in 1:npops) amova.fake.ht[pops[i]:(pops[i] + nseqs[i] - 1), i] <- 1
    amova.fake.ht <- as.data.frame(amova.fake.ht)
    amova.all <- apply(absdist, 2, ms.amova, samples = data.frame(amova.fake.ht), nseqs = nseqs)
    phi.st.amova <- unlist(sapply(amova.all, get.phist), use.names = F)
    phi.st <- cbind(phi.st, phi.st.amova)
  }
  ## Store the calculated statistics that are vectors in a data frame
  df.1 <- data.frame(S.total = segsites, S = t(S.all), pi.total = pi.tot, pi = t(pi.all), nhap.total = sapply(ht, nrow), nhap = t(ht.pops), Hd.total = Hd.all, Hd = t(Hd.pops), tajD = t(tajD.all), pi.s = apply(pi.all, 2, mean))
  ## Some of these only exist if npops > 1 or npops > 2, and for some the structure of the statistic depends on npops
  if(npops > 1) df.1 <- cbind(df.1, phi.st)
  if(npops > 1) {
    d.cols <- c()
    for (i in 2:(npops)) {  ## new labeling loop Dec2011 to fix problem JP pointed out of mislabeled columns
      for (j in 1:(i-1)){
        d.cols <- c(d.cols, paste("d.", j, i, sep=""))
      }
    }   
    if (npops > 2) d.df <- data.frame(t(d.all.short))  ## creates a data frame
    if (npops == 2) d.df <- data.frame(d.all.short)  ## necessary because with only 2 pops there's a vector of d values, but with more than 2 there's a matrix
    names(d.df) <- d.cols
    numeric.output <<- data.frame(df.1, d.df)
    if (npops > 2) {  ## if > 2 do the pairwise phi.st calculations and put them in
      phist.df <- data.frame(phist=t(phi.st.pairs))
      numeric.output <<- data.frame(df.1, d.df, phist.df)
    }    
  } else {  ## if npops == 1
    numeric.output <<- df.1
  }
  #if(npops == 2) numeric.output$d.all <<- d.all
  if(outfile == "") {
    writeLines("Done!  You now have access to the following variables:
               Variables that are not per simulation:
               headerlist has the information from the ms command line (the \"header\" in the input file)           
               nsims is the number of simulations
               npops is the number of populations
               nseqs is a vector of the number of seqs per population
               pops is a vector of the first seq in each population (i.e. where in the total list of seqs each pop starts)
               
               Variables that list nicely are in a data frame named numeric.output:
               S.total is the number of segregating sites across all populations (straight from ms)
               pi.n is pi for population n (e.g. pi.1, pi.2, etc)
               pi.total is pi if you consider all seqs to come from the same population
               S.n is the number of segregating sites for population n (e.g. S.1, S.2...)
               nhap.total is the number of unique haplotypes in each simulation across all populations
               nhap.n is the number of unique haplotypes in population n
               tajD.n is Tajima's D for population n
               pi.t is pi across populations
               pi.s is the mean of the population pis
               Hd.n is haplotype diversity in population n
               Hd.total is haplotype diversity considering all seqs to come from the same population
               phi.st is Phi.st calculated the \"simple\" way from pi.t and pi.s
               phist.xy is the pairwise Phi.st for populations x and y
               phi.st.amova is Phi.st calculated by the amova
               d.xy is divergence between pops x and y (Nei's Dxy, not Da)
               
               Variables that don't list nicely (generally because they have unpredictable numbers of values, or sizes, per simulation):
               d.all has the divergences, in nsims columns and n x n rows IF n != 2; fold into a matrix if you want
               ht has the haplotype counts per pop, in a list of nsims little nhaps by npops arrays
               ht.freqs is like ht only frequency (per pop) rather than count
               ht.dists is a list of nsims distance matrices, each nhaps x nhaps
               amova.all is a list of nsims amovas
               absdist has the distance matrix; each column is one simulation, unfolded into a single line (so dim(absdist) == c(sum(nseqs)^2, nsims))
               
               Useful things to remember:
               \t*call functions with their arguments in parentheses, e.g. mean(pi.all)
               \t*pull out individual entries in lists, arrays, or vectors with square brackets, e.g. pi.all[1,1]
               \t*pull out rows or columns of arrays with square brackets and one number missing, e.g. pi.all[,1] or pi.all[1,]
               \t*use the colon to get ranges, e.g. 1:10 will give you integers 1,2,3,...,10.  Useful with subsets, e.g. pi.all[1,1:10]
               \t*you can always use an expression rather than a number, e.g. pi.all[1, 1:nseqs[1]] or even pi.all[1,pops[1]:(pops[1] + nseqs[1] - 1)]
               \t*just type the name of a variable to see the current value of that variable
               \t*use help(command) to see a description, options, and syntax for a command")
  }	
  ##
  ## UN-COMMENT THESE LINES IF YOU WANT THE MEANS/MEDIANS/STDEVS, SOME PLOTS, ETC:
  ##
  #sample stats and output
  #pi.mean <- apply(pi.all, 2, mean)
  #pi.median <- apply(pi.all, 2, median)
  #pi.stdev <- apply(pi.all, 2, sd)
  #S.mean <- mean(segsites)
  #S.median <- median(segsites)
  #S.stdev <- sd(segsites)
  #d.mean <- mean(d.all)
  #d.median <- median(d.all)
  #d.stdev <- sd(d.all)
  #ms.stats <<- data.frame(pi.mean, pi.median, pi.stdev, S.mean, S.median, S.stdev, d.mean, d.median, d.stdev)
  #writeLines("OUTPUT:\n\n")
  #print(ms.stats)
  #if (graph != "none"){
  #	if (graph == "screen") x11(50,20) else (postscript(file = graph, width = 10, height = 7.5))
  #		par(mfcol = c(npops, 4))
  #		for (i in 1:npops) hist(pi.all[,i], nclass = 50, col = c("red", "blue", "orange", "purple", "yellow", "black", "gray")[i])
  #		plot(segsites, pi.all[,1], col = "red")
  #		plot(segsites, pi.all[,2], col = "blue")
  #		hist(segsites, nclass = (max(segsites) - min(segsites) + 1))
  #		hist(d.all, nclass = 50)
  #		plot(pi.all[,1], pi.all[,2], col = "green")
  #		plot(segsites, d.all, col = "orange")
  #		if(graph != "screen") dev.off()
  #	}
  
  ## If an outfile was specified, write the output to it:
  if(outfile != ""){
    if(out.opt %in% c("medium", "m", "long", "l")) write.table(format(numeric.output, digits=6), outfile, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE, na = "NA", quote=FALSE)
    ## Note changes in v4.1: NA entries are now "NA" in the output, rather than blank;
    ##   precision set to 6 sig figs; quoting of strings (including headers) removed
    if(out.opt %in% c("long", "l")) {
      sink(paste(outfile, ".ht", sep = ""))
      print(ht)
      print(ht.freqs)
      print(ht.dists)
      sink()
      sink(paste(outfile, ".amova", sep = ""))
      print(amova.all)
      sink()
    }
  }
  ## The main object returned is a data frame containing most of the calculated statistics.
  ## However, this (and quite a few other objects) are created by the script via the '<<-' operator,  
  ##  so the script can be called from the R prompt without assigning the output to anything.
  ##  This is inelegant but designed to be easier for non-R-experts.  
  invisible(numeric.output)
  }


output<-ms.output("ms.out",haplotypes=T)

ms.header <- function(ms.outfile){
  header <- ms.outfile[1:(grep("//", ms.outfile)[1] - 1)]
  header.list <- c()
  if(any(header == "-I")) {		
    npops <<- as.numeric(header[grep("-I", header)[1] + 1])
    nseqs <<- as.numeric(header[grep("-I", header)[1] + (2:(npops + 1))])
    header.list <- c(header.list, npops = npops, nseqs = nseqs)
    if(any(header == "-g")) {  # removed as.numeric() Dec2011
      gps <- header[grep("-g", header) + 1]
      g <- header[grep("-g", header) + 2]
      header.list <- c(header.list, g = g)
    }
  }
}
ms.header(ms.in)

ms.in <- scan("ms.out", what = "c", quiet = TRUE)
message(paste("Reading", filename, "and calculating.... please be patient; I have a lot of pairwise comparisons to do!\n"))
## strip the tbs arguments after //
slashes <- grep("//", ms.in)

# source for this
# https://datadryad.org/bitstream/handle/10255/dryad.36550/ms.out.v4.2.R?sequence=1