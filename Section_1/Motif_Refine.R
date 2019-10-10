suppressWarnings(suppressMessages({
  
  elapsed <- system.time({
    
    args <- commandArgs(trailingOnly=TRUE)
    
    libpath <- '/Library/Frameworks/R.framework/Versions/3.6/Resources/library/'
    strt <- Sys.time()
    #system("printf '%*s\\n' \"${COLUMNS:-$(tput cols)}\" '' | tr ' ' =")
    cat("\n\nImporting necessary libraries\n")
    
    library(data.table, lib.loc = libpath);cat("data.table ")
    library(stringi, lib.loc = libpath);cat("stringi ")
    library(stringr, lib.loc = libpath);cat("stringr ")
    library(plyr, lib.loc = libpath);cat("plyr ")  # count conflict with
    library(nnet, lib.loc = libpath);cat("nnet ")
    #library(qpcR, lib.loc = libpath);cat("qpcR ")  # this
    library(seqinr, lib.loc = libpath);cat("seqinr ")
    library("Biostrings");cat("Biostrings ")
    library(parallel, lib.loc = libpath);cat("parallel ")
    library(R.utils);cat("R.utils ")
    library(ggplot2);cat("ggplot2 ")
    library(gridExtra);cat("gridExtra ")
    library(reshape2);cat("reshape2 ")
    library(foreach);cat("foreach ")
    #library(doMC);cat("doMC ")
    library(RSQLite);cat("RSQLite ")
    library(gplots);cat("gplots ")
    
    #Feeding it with Parameter_*.txt file
	## READ SETTINGS FROM TEXT FILE ##########################################
	parameters <- args[1]
  	lines <- readLines(parameters)
	cl1 <- c(lines[grep("PARAMETERS", lines)+1],lines[grep("STRAIN_FOLDERS", lines)-1])
    cl_max <- as.numeric(sub(".*=","",cl1[cl1 != ""]))  
    strain1 <- c(lines[(grep("STRAIN_FOLDERS", lines)+1):(grep("OUTPUT_DIRECTORY", lines)-1)])
    directories <- sub(".*=","",strain1[strain1 != ""])
    motiftable <- read.table(parameters, skip=grep("MOTIFS", lines), header=T,sep=",")
    outdir1 <- c(lines[(grep("OUTPUT_DIRECTORY", lines)+1):(grep("MOTIFS", lines)-1)])
    outdir <- sub(".*=","",outdir1[outdir1 != ""])
    
    #Feeding it with Modifications file
    # core max 8
    ### BIG CSV FILE VARS ####################################################
    cov_cut <- c(10,20,30)  # coverage cutoffs                 
    cov_cut_spec <- 10                  # cov_cut for ranksum & percentage
    size <- c(100,500,1000)     # random size samples              
    size_spec <- 100                     # sample size for ranksum table    
    percentage <- 50                    # percentage threshold             
    ## CHECKED OR UNCHECKED VARIABLES ########################################
    rthres <- 3                        # ratio threshold
    ##########################################################################
    
    cat("\n")
    
    #Checking if outdir exists otherwise create it
    unlink(outdir, recursive = T)
    ifelse(!dir.exists(outdir), dir.create(outdir), FALSE)
    
    #Checks the genomes present
    #directories <- #system("find . -maxdepth 2 | grep -E '[0-9]{6}$'", intern=T)[c(4,16,49,53,54,55)] # changed for troubleshooting
    namesg2 <- as.numeric(gsub(".*/","",directories))
    genomes <- directories
    # separate by genome type
    gentypes <- gsub(".*/","",sub("/[^/]*$", "", genomes))
    genunique <- unique(gentypes)
    cat("Genomes detected:", "\n")
    for(unic in 1:length(genunique)){
      unic2 <- namesg2[which(gentypes == genunique[unic])]
      unic3 <- unname(tapply(unic2,cumsum(c(TRUE,diff(unic2)!=1)), FUN= function(x) 
        if(length(x)>1) paste(range(x), collapse='-') else x))
      cat(paste(" ", genunique[unic]),"\n")
      cat(paste("  +", unic3, collapse = "\n"), "\n")
    }
    
    cat("\n")
    
    # functions
    #drawline <- function(){system("printf '%*s\\n' \"${COLUMNS:-$(tput cols)}\" '' | tr ' ' =")}
    # imported mgsub function from qdap
	
	cbind.na <- function (..., deparse.level = 1) 
{
    na <- nargs() - (!missing(deparse.level))
    deparse.level <- as.integer(deparse.level)
    stopifnot(0 <= deparse.level, deparse.level <= 2)
    argl <- list(...)
    while (na > 0 && is.null(argl[[na]])) {
        argl <- argl[-na]
        na <- na - 1
    }
    if (na == 0) 
        return(NULL)
    if (na == 1) {
        if (isS4(..1)) 
            return(cbind2(..1))
        else return(matrix(...))
    }
    if (deparse.level) {
        symarg <- as.list(sys.call()[-1L])[1L:na]
        Nms <- function(i) {
            if (is.null(r <- names(symarg[i])) || r == "") {
                if (is.symbol(r <- symarg[[i]]) || deparse.level == 
                  2) 
                  deparse(r)
            }
            else r
        }
    }
    if (na == 0) {
        r <- argl[[2]]
        fix.na <- FALSE
    }
    else {
        nrs <- unname(lapply(argl, nrow))
        iV <- sapply(nrs, is.null)
        fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
        if (deparse.level) {
            if (fix.na) 
                fix.na <- !is.null(Nna <- Nms(na))
            if (!is.null(nmi <- names(argl))) 
                iV <- iV & (nmi == "")
            ii <- if (fix.na) 
                2:(na - 1)
            else 2:na
            if (any(iV[ii])) {
                for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i))) 
                  names(argl)[i] <- nmi
            }
        }
        nRow <- as.numeric(sapply(argl, function(x) NROW(x)))
        maxRow <- max(nRow, na.rm = TRUE)
        argl <- lapply(argl, function(x) if (is.null(nrow(x))) 
            c(x, rep(NA, maxRow - length(x)))
        else rbind.na(x, matrix(, maxRow - nrow(x), ncol(x))))
        r <- do.call(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
    }
    d2 <- dim(r)
    r <- cbind2(argl[[1]], r)
    if (deparse.level == 0) 
        return(r)
    ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
    ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
    if (ism1 && ism2) 
        return(r)
    Ncol <- function(x) {
        d <- dim(x)
        if (length(d) == 2L) 
            d[2L]
        else as.integer(length(x) > 0L)
    }
    nn1 <- !is.null(N1 <- if ((l1 <- Ncol(..1)) && !ism1) Nms(1))
    nn2 <- !is.null(N2 <- if (na == 2 && Ncol(..2) && !ism2) Nms(2))
    if (nn1 || nn2 || fix.na) {
        if (is.null(colnames(r))) 
            colnames(r) <- rep.int("", ncol(r))
        setN <- function(i, nams) colnames(r)[i] <<- if (is.null(nams)) 
            ""
        else nams
        if (nn1) 
            setN(1, N1)
        if (nn2) 
            setN(1 + l1, N2)
        if (fix.na) 
            setN(ncol(r), Nna)
    }
    r
}
	
    mgsub <- function (pattern, replacement, text.var, leadspace = FALSE, 
                       trailspace = FALSE, fixed = TRUE, trim = TRUE, order.pattern = fixed, 
                       ...) 
    {
      if (leadspace | trailspace) 
        replacement <- spaste(replacement, trailing = trailspace, 
                              leading = leadspace)
      if (fixed && order.pattern) {
        ord <- rev(order(nchar(pattern)))
        pattern <- pattern[ord]
        if (length(replacement) != 1) 
          replacement <- replacement[ord]
      }
      if (length(replacement) == 1) 
        replacement <- rep(replacement, length(pattern))
      for (i in seq_along(pattern)) {
        text.var <- gsub(pattern[i], replacement[i], text.var, 
                         fixed = fixed, ...)
      }
      if (trim) 
        text.var <- gsub("\\s+", " ", gsub("^\\s+|\\s+$", "", 
                                           text.var, perl = TRUE), perl = TRUE)
      text.var
    }
    # reverse function
    myreverse <- function(x,...){unlist(lapply(strsplit(as.vector(x),""),function(z)paste(rev(z),collapse="")))}
    
    # shorten motif
    shortenme <- function(x){
      r <- rle(unlist(strsplit(x, "")))
      m <- c(rbind(r$values, r$lengths))
      paste(ifelse(m == 1, "", m), collapse="") 
    }
    
    
    # prettify
    prettify <- theme(panel.background = element_rect(fill = NA,color="gray"), 
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(size=.1, color="black",linetype="dotted"), 
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_line(size=.1, color="black"),
                      legend.position="bottom")
    
    # TMUX SPLITTER
#     nom_tmux <- as.numeric(args[2]) #12
#     nomnom <- as.numeric(args[1]) #i
    
    namesg3 <- namesg2#split(namesg2, ceiling(seq_along(namesg2)/floor((length(namesg2)/nom_tmux))))[[nomnom]]  
    #split(namesg2, countme5)[[nomnom]]  
    genomes3 <- genomes#split(genomes, ceiling(seq_along(genomes)/floor((length(genomes)/nom_tmux))))[[nomnom]]  
    #split(genomes, countme5)[[nomnom]]
    
    ## SEPARATE BY GENOME ###########################################
    for(set in 1:length(namesg3)){                                  #
      csv_file <- paste0(genomes3[set],"/modifications.csv.gz")      #
      fasta_file <- paste0(genomes3[set],"/genome.fasta")            #
      motif_file <- paste0(genomes3[set],"/motif_summary.csv")       #
      
      name <- namesg3[set]
       
      cat(paste("Working on", name), "\n")
      dog2s <- 0
      
      # FASTA FILE
      cat(paste(" + Reading genome.fasta"), "\n")
      gene <- readDNAStringSet(fasta_file)  # DNA String class
      
      # MOTIF SUMMARY CSV
      cat(paste(" + Reading motif_summary.csv"), "\n")
      #motifsfiles <- list.files(pattern = "motif_summary.csv$", recursive = TRUE)
      #allmotifs <- do.call(rbind,lapply(motifsfiles,function(z)read.csv(z)))[,1:4]
      
      ##final motifs dataframe
      # mMS <- c("CTTTANNNNNNNTG", "CANNNNNNNTAAAG", "GCMGAAG", "CAAAAA", "TAGNNNNNRTGAA", 
               # "TTCAYNNNNNCTA", "GAAGNNNNNRAAT", "ATTYNNNNNCTTC", "CATCG", "CGATG", 
               # "GCAYNNNNNNTTC", "GAANNNNNNRTGC", "CAGNNNNNNTTC", "GAANNNNNNCTG", 
               # "GATC", "GACNNNNTTTA", "RGACNNNNNRCT", "AGYNNNNNGTCY", "CAGNNNNNNCRT", 
               # "AYGNNNNNNCTG", "CTGCAG", "GCNGCAGCNNB", "GCWGCV", "ATCNNNNNNCTC", 
               # "GAGNNNNNNGAT", "TGANNNNNNCTC", "GAGNNNNNNTCA", "CAAYNNNNNCTT", 
               # "AAGNNNNNRTTG", "CCANNNNNNNTAYC", "GRTANNNNNNNTGG", "CAYCNNNNTTTA", 
               # "TAAANNNNGRTG", "TAACTG", "TTRAAYG", "GACNNNNNACT", "AGTNNNNNGTC", 
               # "GAANNNNNNNCTC", "GAGNNNNNNNTTC", "TCANNNNNGTC", "GACNNNNNTGA", 
               # "RGAAAGR", "AYCNNNNRTGT", "ACAYNNNNGRT", "RGACNNNNNNTAG", "CTANNNNNNGTCY", 
               # "AYCNNNNNCTRC", "GYAGNNNNNGRT", "CAAAA", "CATCC", "GGATG", "GCATC", 
               # "GATGC","GGSCCM")
      # mCP <- c(4L, 1L, 5L, 5L, 1L, 3L, 2L, 0L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 
               # 1L, 2L, 0L, 1L, 0L, 4L, 1L, 1L, 0L, 1L, 2L, 1L, 2L, 1L, 2L, 3L, 
               # 1L, 2L, 3L, 4L, 1L, 0L, 2L, 1L, 2L, 1L, 4L, 0L, 2L, 2L, 2L, 0L, 
               # 2L, 4L, 1L, 2L, 2L, 1L, 1L)
      # mMT <- c("m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", 
               # "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", 
               # "m6A", "m6A", "m6A", "m4C", "m4C", "m6A", "m6A", "m6A", "m6A", 
               # "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m4C", "m6A", "m6A", 
               # "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", 
               # "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m6A", "m4C")
      # mF <- c(0.9807136, 0.97782063, 0.97613066, 0.9661038, 0.98872787, 0.9742351, 
              # 0.9871668, 0.9832182, 0.92561984, 0.90082645, 1, 0.9929578, 0.9959547, 
              # 0.99433655, 0.9913669, 0.6087613, 1, 0.99862826, 0.9922179, 0.98346305, 
              # 0.9705882, 0.5314286, 0.38787377, 0.70463735, 0.6946492, 0.60175145, 
              # 0.5942452, 0.4642397, 0.46134022, 0.45861602, 0.45771143, 1, 
              # 0.9989407, 0.9842437, 0.9928469, 0.992236, 0.98913044, 0.99148697, 
              # 0.9801362, 0.994363, 0.9921082, 0.9880056, 0.6639828, 0.66317767, 
              # 0.99667776, 0.9817276, 0.9961046, 0.99499166, 0.06592127, 0.06959002, 
              # 0.047552202, 0.13474387, 0.11183582, 0.22867298)
      
      # motiftable <- data.frame(motifString = mMS, centerPos = mCP, modificationType = mMT, fraction = mF)
      
      #motiftable <- read.table(motif_file, header = T, sep=",")[,1:4]
      
      motiftable <- rbind(motiftable[!motiftable$modificationType == "modified_base",],motiftable[motiftable$modificationType == "modified_base",])
      motofle <- as.vector(motiftable$motifString)
      modtype <- as.vector(motiftable$modificationType)
      modnom <- as.vector(motiftable$centerPos)
      percmeth <- as.vector(motiftable$fraction)
      urder2 <- which(!modtype %in% c("modified_base"))
      urder3 <- which(modtype %in% c("modified_base"))
      methvec <- gsub("\\d|[[:lower:]]","",modtype)
      urder <- modtype %in% c("modified_base")
      if(!length(urder3) == 0){
        for(mayo in 1:length(urder3)){
          if(grepl("A",substr(motofle[urder3][mayo],modnom[urder3][mayo]+1,modnom[urder3][mayo]+1))){
            urder[urder3][mayo] <- FALSE
            methvec[urder3][mayo] <- "A"
          }else{
            if(grepl("C",substr(motofle[urder3][mayo],modnom[urder3][mayo]+1,modnom[urder3][mayo]+1))){
              urder[urder3][mayo] <- FALSE
              methvec[urder3][mayo] <- "C"
            }else{
              urder[urder3][mayo] <- TRUE
            }
          }
        }  
      }
      
      nodataatpos <- vector()
      
      if(nrow(motiftable) == 0){
        cat("Error motif_summary.csv is empty. Moving onto next genome.", "\n")
        next
      }
      
      
      # CSV FILE
      cat(paste(" + Reading modifications.csv"), "\n")
      csv <- fread(gunzip(csv_file, remove = F, temporary = T, overwrite = T), sep = ",", header = T, verbose = F, drop = c(6, 7, 8, 11, 12, 13))  # Set csv parameters
      
      moditype <- modtype[!urder]
      moditype_c <- modtype[urder]
      motif_f <- motofle[!urder]
      motif_mod <- motofle[urder]
      perc_meth <- percmeth[!urder]
      perc_meth_c <- percmeth[urder]
      
      methvec <- methvec[!urder]
      
      motif_r <- myreverse(motif_f)
      motif_ref <- data.frame(sym = c("W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"), 
                              bases = c("(A|T)", "(C|G)", "(A|C)", "(G|T)", "(A|G)", "(C|T)", "(C|G|T)", "(A|G|T)", "(A|C|T)", "(A|C|G)", "(A|C|G|T)"))
      motif_ref_cut <- transform(motif_ref, bases = motif_ref$bases <- 
                                   gsub("[A-Za-z ']*(?:(?:~~|&&)[A-Za-z ']*)*\\K(?:[^A-Za-z ']|\\z)", "", motif_ref$bases, perl = TRUE))
      dir <- paste0(genomes3[set],"/")
      nom_matrix <- matrix(1:(length(cov_cut) * length(size)), nrow = length(cov_cut), ncol = length(size), byrow = T)
      csv_list <- list()
      final_motifs <- vector()
      final_motifs6 <- vector()
      # create list of csv_specific
      for (i in 1:length(cov_cut)) {
        csv1 <- csv[base == "A"]
        csv2 <- csv1[coverage > cov_cut[i], ]
        # cat(paste("  + coverage:", cov_cut[i], "| size: "))
        for (j in 1:length(size)) {
          # percentage cutoff
          if (nrow(csv2) * 100/nrow(csv1) < percentage) {
            jui <- data.table(data.frame(matrix(, nrow = 0, ncol = ncol(csv)+8)))
            setnames(jui, names(jui), c(names(csv), "name", "cov_cutoff", "samp_size", "percent", "scoremean", "scoremedian", "ipdratiomean", "ipdratiomedian"))
            csv_list[[nom_matrix[i, j]]] <- jui
          } else {
            csv_list[[nom_matrix[i, j]]] <- data.table(csv2[sample(1:nrow(csv2), size[j]), ], 
                                                       name = name, cov_cutoff = cov_cut[i], 
                                                       samp_size = size[j], percent = nrow(csv2) * 100/nrow(csv1), 
                                                       scoremean = mean(csv2[["score"]]), 
                                                       scoremedian = median(csv2[["score"]]), 
                                                       ipdratiomean = mean(csv2[["ipdRatio"]]), 
                                                       ipdratiomedian = median(csv2[["ipdRatio"]]))
          }
          # cat(paste(size[j],""))
        }
        # cat("\n")
      }
      csv_listo <- do.call(rbind, lapply(csv_list, function(z) z[1, (ncol(z)-7):ncol(z), with = F]))
      csv_listo <- csv_listo[complete.cases(csv_listo),]
      
      # create data folder
      dir2 <- paste0(dir, "data/")
      unlink(dir2, recursive = T)
      ifelse(!dir.exists(dir2), dir.create(dir2), FALSE)
      
      # create results folder
      dir4 <- paste0(dir, "results/")
      unlink(dir4, recursive = T)
      ifelse(!dir.exists(dir4), dir.create(dir4), FALSE)
      
      # WRITE DATA/COV_SIZE.CSV
      wfilename <- paste0(name, "_covsize.csv")
      filename <- paste0(dir2, name, "_covsize.csv")
      cat(paste(" + Writing to", wfilename), "\n")
      write.csv(csv_listo, filename, row.names = F)
      
      # Create lists and vectors
      contig_vec <- vector()  # contig names vector
      control <- vector()
      r_cutoff <- vector()
      dog2_list <- list()
      positive_list <- list()
      negative_list <- list()
      timegat <- list()
      
      # take longest contig
      gag <- which.is.max(width(gene))
      cat(" + Analyzing Motifs", "\n")
      
      # change genome_f and csv to match contig
      cat(paste("  + Reading largest contig:", names(gene)[gag]), "\n")
      genome_f <- toString(gene[[gag]])
      csv2 <- csv[csv$refName == names(gene)[gag]]
      
      # R and F prep
      csv2 <- csv2[complete.cases(csv2), ]
      csv2 <- csv2[order(csv2$tpl),]
      param <- c(csv2$tpl[1], csv2$tpl[nrow(csv2)])  # create reverse genome
      genome_r <- chartr("GATC", "CTAG", genome_f)  # chart
      csv_r <- csv2[(csv2$strand == 1), ]  # separate csv r by strandness
      csv_f <- csv2[(csv2$strand == 0), ]  # separate csv f by strandness
      # cluster prep
      csv_sp <- csv[csv$coverage >= cov_cut_spec, ]
      csv_sp_f <- csv_f[csv_f$coverage >= cov_cut_spec, ]
      csv_sp_r <- csv_r[csv_r$coverage >= cov_cut_spec, ]
      
      exportme2 <- list(data.table(csv_sp_f), data.table(csv_sp_r), genome_f)
      
      # create cluster
      starterp <- Sys.time()
      cat(paste0("  + Making cluster of ", cl_max, " cores"), "\n")
      cat("   + Exporting global variables/libraries", "\r")
      cl <- makeCluster(cl_max)
      clusterExport(cl, c("exportme2"))
      clusterEvalQ(cl, {
        library(stringi)
        library(IRanges)
        library(data.table)
      })
      enderp <- Sys.time()
      elapsederp <- signif(enderp - starterp, 3)
      flush.console()
      cat(paste("  + Exporting global variables/libraries:", elapsederp, "sec"), "\n")
      j <- 1
      jmax <- length(motif_f)
      refyon <- vector()
      allvarsbwhile <- c(ls(),"allvarsbwhile")
      while (j <= jmax) {
        refyon[j] <- 0
        r_cutoff[j] <- 0
        # create data/motif folder
        dir3 <- paste0(dir2, motif_f[j], "/")
        ifelse(!dir.exists(dir3), dir.create(dir3), FALSE)
        meth <- methvec[j]
        
        cat(paste("  + Motif:", motif_f[j]), "\n")
        motif_r <- myreverse(motif_f[j])  # Reverse motif_f
        mot_f <- stri_replace_all_fixed(motif_f[j], motif_ref$sym, motif_ref$bases, vectorize_all = F)  
        # create and translate forward and reverse motifsf
        mot_r <- paste(rev(stri_extract_all(mot_f, regex = "\\([^)]+\\)|.")[[1]]), collapse = "")
        
        # automatic methlyation detector
        methloc <- function(csv_t, genome_t, mot_t) {
          # weird error skip
          csv_t <- csv_t[order(csv_t$tpl)]
          param2 <- c(csv_t$tpl[1], csv_t$tpl[nrow(csv_t)])
          loc_pre <- data.frame(rbind.fill(stri_locate_all_regex(genome_t, mot_t)))  # forward[1] and reverse[2] same result
          loc_pre <- loc_pre[!(loc_pre$start < param2[1] | loc_pre$start > param2[2]), ]  # shave loc_pre
          loc_locate <- stri_locate_all_regex(substring(genome_t, loc_pre[1:nrow(loc_pre), 1], loc_pre[1:nrow(loc_pre), 2]), meth)
          
          # Expand the sequences
          loc_pre_list <- data.frame(index = loc_pre[, 1])
          loc_pre_list <- split(loc_pre_list, rownames(loc_pre_list))
          binded <- binded <- Map(function(z,...)cbind(z,...,row.names = NULL), loc_locate, index = loc_pre_list)
          binded2 <- Map(cbind, binded, lapply(lapply(lapply(binded, 
                                                             function(z) z[, "start"] - 1 + z[, "index"]), as.data.frame), setNames, nm = c("tpl")))
          binded2_list <- data.frame(id = 1:sum(sapply(binded2, nrow)))
          binded2_list <- split(binded2_list, rep(1:length(sapply(binded2, nrow)), sapply(binded2, nrow)))
          binded2 <- Map(cbind, binded2, binded2_list)
          binded2.a <- do.call("rbind", binded2)
          binded3 <- data.frame(subset(merge(csv_t, binded2.a, by = "tpl"), 
                                       select = c("base", "start", "score", "tpl", "id")))
          binded3 <- binded3[order(binded3$id), ]
          remove_me <- binded2.a[which(!binded2.a$tpl %in% binded2.a$tpl[binded2.a$tpl %in% binded3$tpl]), ]
          binded2 <- lapply(binded2, function(z) z[!(z$id %in% remove_me$id), ])
          vectornom <- sapply(binded2, nrow)
          score_loc <- lapply(split(binded3, rep(1:length(vectornom), vectornom)), data.frame)
          names(score_loc) <- NULL
          score_loc2 <- plyr:::count(do.call(rbind, lapply(score_loc, function(z) z[which.is.max(z$score), ]$start))[, 1])
          # count plyr and qcdR conflict
          transl <- as.numeric()
          transl <- score_loc2[which.is.max(score_loc2$freq), 1] - 1
          transl
        }
        
        # cat("   + centerPos Checking","\n")
        # transl2 <- methloc(csv_f, genome_f, mot_f)
        transl_r <- transl_f <- motiftable$centerPos[j]
        
        if (is.na(transl_f)) 
          next  # skip contig if transl_f return NA
        
        # create locations and motif table
        locations_f <- data.frame(data.frame(matrix(unlist(stri_locate_all_regex(genome_f, mot_f)), ncol = 2))[1] + transl_f, 
                                  motif = unlist(stri_extract_all_regex(genome_f, mot_f)))
        locations_r <- data.frame(data.frame(matrix(unlist(stri_locate_all_regex(genome_r, mot_r)), ncol = 2))[2] - transl_r, 
                                  motif = myreverse(unlist(stri_extract_all_regex(genome_r, mot_r))))
        # shave motif locations with param
        loc_f_s <- locations_f[!(locations_f$X1 < param[1] | locations_f$X1 > param[2]), ]
        loc_r_s <- locations_r[!(locations_r$X2 < param[1] | locations_r$X2 > param[2]), ]
        # merge locations and scores
        merged_f <- merge(data.frame(tpl = loc_f_s$X1, motif = loc_f_s$motif), csv_f, all = F)
        merged_r <- merge(data.frame(tpl = loc_r_s$X2, motif = loc_r_s$motif), csv_r, all = F)
        # merge forward and reverse motif locations 
        dat <- data.table(rbind(merged_f, merged_r))
        
        ## FIND DATA FOR NEGATIVE CONTROLS AND PLOT NEXT TO DISTRO FOR SCORE AND IPDRATIO
        mot_wrong <- unlist(strsplit(motif_f[j],"")) 
        mot_wrong_pos <- which(mot_wrong %in% c("A","T","C","G"))
        mot_wrong2 <- chartr("AGTC","TCAG",mot_wrong[mot_wrong_pos])
        mot_wrong[mot_wrong_pos] <- mot_wrong2
        motif_w <- paste(mot_wrong, collapse = "")
        mot_f_w <- stri_replace_all_fixed(motif_w, motif_ref$sym, motif_ref$bases, vectorize_all = F)  
        # create and translate forward and reverse motifs
        mot_r_w <- paste(rev(stri_extract_all(mot_f_w, regex = "\\([^)]+\\)|.")[[1]]), collapse = "")
        locations_f_w <- data.frame(data.frame(matrix(unlist(stri_locate_all_regex(genome_f, mot_f_w)), ncol = 2))[1] + transl_f, 
                                    motif = unlist(stri_extract_all_regex(genome_f, mot_f_w)))
        locations_r_w <- data.frame(data.frame(matrix(unlist(stri_locate_all_regex(genome_r, mot_r_w)), ncol = 2))[2] - transl_r, 
                                    motif = myreverse(unlist(stri_extract_all_regex(genome_r, mot_r_w))))
        # shave motif locations with param
        loc_f_s_w <- locations_f_w[!(locations_f_w$X1 < param[1] | locations_f_w$X1 > param[2]), ]
        loc_r_s_w <- locations_r_w[!(locations_r_w$X2 < param[1] | locations_r_w$X2 > param[2]), ]
        # merge locations and scores
        merged_f_w <- merge(data.frame(tpl = loc_f_s_w$X1, motif = loc_f_s_w$motif), csv_f, all = F)
        merged_r_w <- merge(data.frame(tpl = loc_r_s_w$X2, motif = loc_r_s_w$motif), csv_r, all = F)
        # merge forward and reverse motif locations
        datw <- data.table(rbind(merged_f_w, merged_r_w))
        
        motif_sh <- shortenme(motif_f[j])
        motif_sh_w <- shortenme(motif_w)
        
        #         if(nrow(dat)==0 | nrow(datw)==0){
        #           cat(paste("   - Cannot create distrib.png, low count"), "\n")
        #         }else{
        #           hist1 <- ggplot(dat, aes(score)) + 
        #             geom_histogram(binwidth = (max(dat$score)-min(dat$score))/30, aes(fill=..count..)) + 
        #             scale_fill_gradient("Count", low = "green", high = "red") +
        #             scale_x_continuous(limits = c(0, max(dat$score))) +
        #             theme_bw() +
        #             theme(legend.position="none") + 
        #             ggtitle(paste("Score Methyl:", motif_sh))
        #           hist2 <- ggplot(dat, aes(ipdRatio)) + 
        #             geom_histogram(aes(fill=..count..), binwidth = (max(dat$ipdRatio)-min(dat$ipdRatio))/30) + 
        #             scale_fill_gradient("Count", low = "green", high = "red") + 
        #             scale_x_continuous(limits = c(0, max(dat$ipdRatio))) +
        #             theme_bw() +
        #             theme(legend.position="none", axis.title.y=element_blank()) + 
        #             ggtitle(paste("IPDRatio Methyl:", motif_sh))
        #           hist1w <- ggplot(datw, aes(score)) + 
        #             geom_histogram(binwidth = (max(dat$score)-min(dat$score))/30, aes(fill=..count..)) + 
        #             scale_fill_gradient("Count", low = "green", high = "red") + 
        #             #scale_x_continuous(limits = c(0, max(dat$score))) +
        #             theme_bw() +
        #             theme(legend.position="none") + 
        #             ggtitle(paste("Score Control:", motif_sh_w))
        #           hist2w <- ggplot(datw, aes(ipdRatio)) + 
        #             geom_histogram(aes(fill=..count..), binwidth = (max(dat$ipdRatio)-min(dat$ipdRatio))/30) + 
        #             scale_fill_gradient("Count", low = "green", high = "red") + 
        #             #scale_x_continuous(limits = c(0, max(dat$ipdRatio))) +
        #             theme_bw() +
        #             theme(legend.position="none", axis.title.y=element_blank()) + 
        #             ggtitle(paste("IPDRatio Control:", motif_sh_w))
        #           
        #           png(paste0(dir2,motif_f[j],"_distrib.png"), width=2000, height=1350, res=120)
        #           grid.arrange(hist1, hist2, hist1w, hist2w, ncol=2)
        #           dev.off()
        #           cat(paste("   - Writing to distrib.png"), "\n")
        #         }
        
        csv_mot <- list()
        levels <- as.vector(unique(dat$motif))
        s_r <- list()
        
        tryCatch({
          cat(paste("   + Analyzing child motifs for coverages:", paste(cov_cut, collapse = ", ")), "\n")
          for (u in 1:length(levels)) {
            mot_spec <- dat[motif == levels[u]]
            for (v in 1:length(cov_cut)) {
              csv_lista <- lapply(lapply(csv_list, function(z) z[, c(5, 6, 9, 10), with = F]), as.data.frame)
              csv_listb <- lapply(csv_lista, function(z) z[(z$samp_size == size_spec), ])
              csv_listc <- do.call(rbind, csv_listb)
              levels2 <- cov_cut_spec  
              for (h in 1:length(levels2)) {
                csv_listd <- csv_listc[csv_listc$cov_cutoff == levels2[h], ][, 1:2]
                p_score_val <- -log10(wilcox.test(mot_spec$score, csv_listd$score, alternative = c("greater"))$p.value)
                p_ipdratio_val <- -log10(wilcox.test(mot_spec$ipdRatio, csv_listd$ipdRatio, alternative = c("greater"))$p.value)
                s_r_sub <- data.frame(p_score_val, p_ipdratio_val)
                names(s_r_sub) <- c(paste0("pval_s.", levels2[h]), paste0("pval_r.", levels2[h]))
                s_r[[h]] <- s_r_sub
              }
            }
            csv_mot[[u]] <- data.frame(motif = levels[u], count = nrow(mot_spec), 
                                       mean_s = mean(mot_spec$score), mean_r = mean(mot_spec$ipdRatio), 
                                       do.call(cbind, s_r))
          }
          csv_mot_listo <- do.call(rbind, csv_mot)
          
          # WRITE RANKSUM.CSV
          wfilename2 <- paste0(name, "_ranksum.csv")
          filename2 <- paste0(dir3, motif_f[j], "_ranksum.csv")
          cat(paste("   + Writing to", wfilename2), "\n")
          write.csv(csv_mot_listo, filename2, row.names = F)
        },error=function(e){cat("   + Error: Could not write to ranksum.csv", "\n")})
        
        tryCatch({
          
          ###### Phase V ####### 
          num_deg <- regmatches(mot_f, gregexpr("(?<=\\().*?(?=\\))", mot_f, perl = T))[[1]]  
          # degenerate translations 'A|C|T' 'A|C|T' 'A|C|G'
          new_mots <- unlist(strsplit(motif_f[j], ""))
          nom3 <- rep(4, length(new_mots))
          nom_matrix2 <- do.call(cbind.na, split(1:sum(nom3), rep(1:length(nom3), nom3)))
          colnames(nom_matrix2) <- NULL
          
          new_list_o <- list()
          new_list_o2 <- list()
          new_mot_ref <- rbind(motif_ref, data.frame(sym = c("A", "C", "G", "T"), bases = c("A", "C", "G", "T")))
          cat(paste0("   + Methy Checking child motifs from pos. 1-", nchar(motif_f[j])), "\n")
          for (q in 1:nchar(motif_f[j])) {
            cat <- 1:nchar(motif_f[j])
            x <- motif_f[j]
            pat2 <- paste0("(", paste(as.vector(new_mot_ref$sym), collapse = "|"), ")")
            rep <- as.vector(new_mot_ref$bases)
            g <- (1:nchar(x))[cat[!cat %in% q]]  #gregexpr(pat2, x)[[1]][cat[!cat %in% q]]
            cat2 <- as.vector(new_mot_ref$sym)
            k <- unlist(strsplit(x, ""))
            for (a in 1:length(g)) {
              k[g[a]] <- rep[match(substr(x, g[a], g[a]), cat2)]
            }
            new_mots2 <- paste0(k, collapse = "")
            #      cat(paste(" + Checking child motif:", new_mots2));threedots()
            num_deg2 <- c("A", "C", "G", "T")  #unlist(strsplit(num_deg[q],'[|]'))
            for (s in 1:length(num_deg2)) {
              new_mots2.a <- strsplit(new_mots2, "(?<=\\))|(?<=[[:alpha:]])(?=[[:alpha:]\\(])", perl = TRUE)[[1]]
              # unlist(stri_extract_all_regex(new_mots2, '\\([ACGT\\|]*\\)|[ACGT]'))
              new_mots2.a[q] <- num_deg2[s]
              new_mots3 <- paste0(new_mots2.a, collapse = "")
              mottleme <- unlist(strsplit(motif_f[j],""))
              mottleme[q] <- num_deg2[s]
              mottleme2 <- paste(mottleme,collapse="")
              # cat(paste(num_deg2[s], ""), "\n")
              # new_mot_ref2 <- new_mot_ref[new_mot_ref$sym == new_mots[q], ] 
              # new_mots3 <- stri_replace_all_fixed(new_mots2, toString(new_mot_ref2$sym), num_deg2[s], vectorize_all = F)
              nom_brak <- str_count(new_mots3, "[\\(]")
              nom_brak_ref <- str_count(motif_f[j], paste0("(", paste(as.vector(motif_ref$sym), collapse = "|"), ")"))
              # NEW DATA DATA FOR LISTED DATAFRAME new_listo <- subset(csv_mot_listo, grepl(new_mots3, motif))
              o_motif <- new_mots3
              o_motif_n <- motif_f[j]
              o_pos <- gregexpr(pat2, x)[[1]][cat[cat %in% q]]
              o_deg_y.or.n <- ifelse(nom_brak < nom_brak_ref, "Y", "N")
              new_list_o[[nom_matrix2[s, q]]] <- data.frame(motif_ = mottleme2, motif = o_motif, motif_n = o_motif_n, 
                                                            pos = o_pos, deg_y.or.n = o_deg_y.or.n)
            }
          }
          new_list_o <- lapply(new_list_o, function(z) replace(z, is.na(z), 0))  # replace NaN's with 0
          #########################################################################################################
          ## PRUNING PROCESS determine checked or unchecked, read in dog <- new_list_o
          dog1.me <- do.call(rbind, new_list_o)  # all checked at this point # count>=5 & mean_r<1.5 & -log10 p <3
          dog1.me2 <- dog1.me[1]
          dog1 <- dog1.me[2:ncol(dog1.me)]          ## NEW DATA prep CSV
          csv_sp_2 <- csv_sp[order(csv_sp$tpl)][complete.cases(csv_sp[order(csv_sp$tpl)]), ]
          param_sp <- c(csv_sp_2$tpl[1], csv_sp_2$tpl[nrow(csv_sp_2)])
          dat2.A.1 <- csv_sp[csv_sp$base == meth, ]
          dat2.A <- dat2.A.1[sample(nrow(dat2.A.1), size_spec), ]
          
          exportme <- list(transl_f, size_spec, param_sp)
          ## PARALLELIZE
          clstart1 <- Sys.time()
          clusterExport(cl, c("exportme"))
          clend.a <- Sys.time()
          cat("   + Merging csv and dat", "\r")
          dat2 <- clusterApply(cl, dog1$motif, function(z) {
            # VARS
            csv_sp_f <- exportme2[[1]]
            csv_sp_r <- exportme2[[2]]
            genome_f <- exportme2[[3]]
            transl_f <- exportme[[1]]
            size_spec <- exportme[[2]]
            param_sp <- exportme[[3]]
            transl_r <- transl_f
            genome_r <- chartr("GATC", "CTAG", genome_f)  # chart
            myreverse <- function(x,...){unlist(lapply(strsplit(as.vector(x),""),function(z)paste(rev(z),collapse="")))}
            ar <- paste(rev(stri_extract_all(z, regex = "\\([^)]+\\)|.")[[1]]), collapse = "")
            
            locations_f2 <- data.frame(data.frame(matrix(unlist(stri_locate_all_regex(genome_f, z)), ncol = 2))[1] + transl_f, 
                                       motif = unlist(stri_extract_all_regex(genome_f, z)))
            locations_r2 <- data.frame(data.frame(matrix(unlist(stri_locate_all_regex(genome_r, ar)), ncol = 2))[2] - transl_r, 
                                       motif = myreverse(unlist(stri_extract_all_regex(genome_r, ar))))
            loc_f_s2 <- locations_f2[!(locations_f2$X1 < param_sp[1] | locations_f2$X1 > param_sp[2]), ]
            loc_r_s2 <- locations_r2[!(locations_r2$X2 < param_sp[1] | locations_r2$X2 > param_sp[2]), ]
            merged_f2 <- merge(data.frame(tpl = loc_f_s2$X1, motif = loc_f_s2$motif), csv_sp_f, all = F)
            merged_r2 <- merge(data.frame(tpl = loc_r_s2$X2, motif = loc_r_s2$motif), csv_sp_r, all = F)
            # merge forward and reverse motif locations
            data.table(rbind(merged_f2, merged_r2))
          })
          
          
          clend2 <- Sys.time() - clend.a
          cat(paste("   + Merging csv and dat", signif(clend2, 3), "sec"), "\n")
          
          # time gather
          tot_size <- as.numeric(object.size(exportme2)) + as.numeric(object.size(exportme))
          timegat[[j]] <- data.frame(size = tot_size, nom_core = cl_max, 
                                     export = elapsederp, merge = clend2)
          j2 <- FALSE
          dog1.a_list <- list()
          
          quantthresh <- 1-perc_meth[j]/2 # or .5
          
          for(jal in 1:NROW(dat2)){
            if(nrow(dat2[[jal]]) == 0){
              o_p_val <- 0
              o_quant_r <- 0
              o_quant_s <- 0
              o_coverage <- 0
            }else{
              o_p_val <- -log10(wilcox.test(dat2[[jal]]$ipdRatio, dat2.A$ipdRatio, alternative = c("greater"))$p.value)#/nrow(dat2[[jal]])
              o_quant_r <- as.numeric(quantile(dat2[[jal]]$ipdRatio, probs = (quantthresh), type = 4))
              o_quant_s <- as.numeric(quantile(dat2[[jal]]$score, probs = (quantthresh), type = 4))
              o_coverage <- mean(dat2[[jal]]$coverage)
            }
            o_totcount <- nrow(dat2[[jal]])
            dog1.a_list[[jal]] <- data.frame(totcount = o_totcount, quant_s = o_quant_s, 
                                             quant_r = o_quant_r, p_val = o_p_val,
                                             mean_cov = o_coverage)
          }
          
          # dat and datw
          # dat and datw
          if(!nrow((dat))==0){
            final_exp <- data.frame(motif_n = motif_f[j], totcount = nrow(dat), 
                                    quant_s = as.numeric(quantile(dat$score, probs = (quantthresh), type = 4)),
                                    quant_r = as.numeric(quantile(dat$ipdRatio, probs = (quantthresh), type = 4)),
                                    p_val_s = -log10(wilcox.test(dat$score, dat2.A$score, alternative = c("greater"))$p.value),
                                    p_val_r = -log10(wilcox.test(dat$ipdRatio, dat2.A$ipdRatio, alternative = c("greater"))$p.value),
                                    mean_cov = mean(dat$coverage))
          }else{
            final_exp <- data.frame(motif_n = motif_f[j], totcount = nrow(dat), 
                                    quant_s = 0,
                                    quant_r = 0,
                                    p_val_s = 0,
                                    p_val_r = 0,
                                    mean_cov = 0)
          }
          
          if(!nrow(datw)==0){
            final_con <- data.frame(motif_n = motif_w, totcount = nrow(datw), 
                                    quant_s = as.numeric(quantile(datw$score, probs = (quantthresh), type = 4)),
                                    quant_r = as.numeric(quantile(datw$ipdRatio, probs = (quantthresh), type = 4)),
                                    p_val_s = -log10(wilcox.test(datw$score, dat2.A$score, alternative = c("greater"))$p.value),
                                    p_val_r = -log10(wilcox.test(datw$ipdRatio, dat2.A$ipdRatio, alternative = c("greater"))$p.value),
                                    mean_cov = mean(datw$coverage))
          }else{
            final_con <- data.frame(motif_n = motif_w, totcount = nrow(datw), 
                                    quant_s = 0,
                                    quant_r = 0,
                                    p_val_s = 0,
                                    p_val_r = 0,
                                    mean_cov = 0)
          }
          write.csv(final_exp, paste0(dir2, motif_f[j], "_exp.csv"), row.names = F)
          write.csv(final_con, paste0(dir2, motif_f[j], "_con.csv"), row.names = F)
          #cat("WITNESS----")
          write.csv(do.call(rbind, dat2), paste0(dir2, motif_f[j], "_expdat2.csv"), row.names = F)
          #cat("ME---------")
          #write.csv(final_con, paste0(dir2, motif_f[j], "_condat.csv"), row.names = F)
          
          cat(paste("   + Writing to exp/con.csv\n"))
          
          dog1.a <- do.call(rbind, dog1.a_list)
          dog2 <- cbind(dog1, dog1.a, id = 1:nrow(dog1))
          dog3 <- cbind(dog1.me2,dog2[c(2,3,5,6,7,8,9)])
          ifelse(!dir.exists(dir3), dir.create(dir3), FALSE)
          # Specific Level 1
          for(dutu in 1:nrow(dog2)){
            filenamer <- paste0("pos", sprintf("%02d",as.numeric(dog2$pos[dutu])), ".csv")
            write.csv(dat2[[dutu]], paste0(dir3, filenamer), row.names = F)
          }
          cat(paste("   + Writing to", paste0(motif_f[j], "_pos1-",nrow(dog2)/4, ".csv")), "\n")
          
          slurm <- data.table(do.call(rbind, dat2))
          
          patdur <- "(?<=\\))|(?<=[[:alpha:]])(?=[[:alpha:]\\(])"
          slurm2 <- list()
          for(slurpy in 1:nrow(dog1)){
            burn <- strsplit(as.character(dog1$motif[slurpy]), patdur, perl=TRUE)[[1]]
            #busedat <- slurm[motif %like% dog1$motif[slurpy]]
            busedat <- slurm[grep(dog1$motif[slurpy], slurm$motif, perl=TRUE),]
            if(!nrow(busedat)==0){
              slurm2[[slurpy]] <- cbind(busedat, 
                                        base2=burn[as.numeric(dog1$pos[slurpy])], 
                                        pos=as.numeric(dog1$pos[slurpy]))
            }else{
              #               dum <- data.table(matrix(ncol = 10, nrow = 0))
              #               setnames(dum, c(colnames(slurm),"base2","pos"))
              dum <- data.table(tpl = NA, motif = NA, refName = slurm$refName[1],strand=0, 
                                base = burn[as.numeric(dog1$pos[slurpy])], score = 0, ipdRatio = 0,
                                base2 = burn[as.numeric(dog1$pos[slurpy])], pos = as.numeric(dog1$pos[slurpy]))
              slurm2[[slurpy]] <- dum
            }
          }
          
          slurm2 <- lapply(slurm2, as.data.frame)
          
          # Specific Level 2
          for(dutu2 in 1:NROW(slurm2)){
            filenamer2 <- paste0("pos", formatC(as.numeric(slurm2[[dutu2]]$pos[1]), width = 2, format = "d", flag = "0"), 
                                 slurm2[[dutu2]]$base2[1],".csv")
            write.csv(slurm2[[dutu2]], paste0(dir3, filenamer2), row.names = F)
          }
          cat(paste("   + Writing to", paste0(motif_f[j], "_pos1-",nrow(dog2)/4, "AGTC.csv")), "\n")
          dog2_list[[j]] <- dog3
          
          checkif <- read.csv(paste0(dir,"motif_summary.csv"), header=T, sep=",")
          perc_meth[j] <- ifelse(motif_f[j] %in% checkif$motifString,
                                 checkif$fraction[which(checkif$motifString %in% motif_f[j])],.50)
          
          #determine cutoff of quant_r
          if(motiftable$modificationType[j] != "m4C"){
            r_cut <- 4 * perc_meth[j]  + (1-perc_meth[j])# Changed to 2 from 3
          }else{
            r_cut <- 3 * perc_meth[j]  + (1-perc_meth[j])
          }
          
          r_cutoff[j] <- r_cut
          
          
          
          j2 <- TRUE
          
          #           checkif <- read.csv(paste0(dir,"motif_summary.csv"), header=T, sep=",")
          #           perc_meth[j] <- ifelse(motif_f[j] %in% checkif$motifString,
          #                                  checkif$fraction[which(checkif$motifString %in% motif_f[j])],.50)
          #           
          #           #determine cutoff of quant_r
          #           if(motiftable$modificationType[j] != "m4C"){
          #             r_cut <- 4 * perc_meth[j]  + (1-perc_meth[j])# Changed to 2 from 3
          #           }else{
          #             r_cut <- 3 * perc_meth[j]  + (1-perc_meth[j])
          #           }
          #           
          #           r_cutoff[j] <- r_cut
          
          # HARD CODE PRUNING Version 2
          removedog <- dog2[!(dog2$quant_r >= r_cut),] # dog2[!(dog2$quant_r >= 3 | dog2$quant_s >= 10 | dog2$p_val >= 5),] # changed to "or" from "and"
          
          # NEXT STEP: Unchecking bases
          cat("   + Methy Unchecking child motifs", "\n") 
          if(nrow(removedog) == nrow(dog2)){
            cat("   + Error: removing all child motifs", "\n")
            final_motifs[j] <- NA
            checkme <- NA
            checkme2 <- NA
            final_motif3 <- NA
            final_motif2 <- NA            
            final_motifs6[j] <- NA
            control[j] <- NA
            contig_vec[j] <- names(gene)[gag]
          }else{
            if (nrow(removedog) == 0) {
              bigdog <- dog2
            } else {
              bigdog <- dog2[-removedog$id, ]
            }
            #  | quant_r < quanthres2)
            
            bigdog2.a <- bigdog
            positive_list[[j]] <- bigdog
            negative_list[[j]] <- removedog
            # check if motif is experimental or control
            if (is.null(bigdog2.a)) {
              final_motifs[j] <- motif_f[j]
              final_motifs6[j] <- motif_f[j]
              control[j] <- "control"
              contig_vec[j] <- names(gene)[gag]
            } else {
              # put motif back together into parent
              cat("   + Reconstructing Motif", "\n")
              findpos <- function(a, b) {
                red <- gsub(" *\\(.*?\\) *", ".", toString(a))
                matrix(substr(red, b, b))
              }
              bigdog2.a$deg <- apply(bigdog2.a, 1, function(z) findpos(z["motif"], z["pos"]))
              bigdog2 <- bigdog2.a[order(bigdog2.a$pos), ]
              vectornom2 <- plyr:::count(bigdog2$pos)$freq
              bigdog3 <- split(bigdog2, rep(1:length(vectornom2), vectornom2))
              bigcat <- data.frame(matrix(unlist(lapply(bigdog3, function(z) list(paste(sort(z$deg), collapse = ""), 
                                                                                  unique(z$pos)))), ncol = 2, byrow = T))
              bigcat$X1 <- mgsub(motif_ref_cut$bases, as.vector(motif_ref_cut$sym), bigcat$X1)
              final_motif <- unlist(strsplit(motif_f[j], ""))
              cat("   + Replacing with references", "\n")
              for (car in 1:nrow(bigcat)) {
                final_motif[as.numeric(as.character(bigcat$X2[car]))] <- bigcat$X1[car]
              }
              
              cat("   + Double Checking", "\n")
              # DOUBLE CHECK
              checkme <- unlist(strsplit(motif_f[j],""))
              checkme2 <- final_motif
              
              final_motif3 <- checkme2 #unshaved
              final_motif2 <- gsub("^N+|N+$", "", paste0(final_motif3, collapse = ""))  # remove N's from head and tail ends
              final_motifs[j] <- final_motif2
              final_motifs6[j] <- paste0(final_motif3, collapse = "")
              control[j] <- "experim"
              contig_vec[j] <- names(gene)[gag]
              refyon[j] <- 1
            } 
          }
          # WRITING TO DEBUG_GRAPH 2 NEW
          cat(paste("   + Writing to debug_graph2.png"), "\n")
          bases <- c("A","C","G","T")
          deglist<- list()
          prettify <- theme(panel.background = element_rect(fill = NA,color="gray"), 
                            panel.grid.major.y = element_blank(),
                            panel.grid.major.x = element_line(size=.1, color="black",linetype="dotted"), 
                            panel.grid.minor.y = element_blank(),
                            panel.grid.minor.x = element_line(size=.1, color="black"),
                            legend.position="bottom")
          
          
          for(newy in 1:nchar(motif_f[j])){
            for(newy2 in 1:4){
              locpos <- bases[newy2]
              namepos <- paste0(dir3, "/pos", sprintf("%02d",newy), locpos, ".csv")
              dataworky10 <- read.csv(namepos, sep = ",", header = T)
              deglist[[((newy-1)*4+newy2)]] <- data.frame(base = locpos, pos = newy, values = c(dataworky10$score, dataworky10$ipdRatio), 
                                                          type = c(rep("score", nrow(dataworky10)), rep("ipdRatio", nrow(dataworky10))),
                                                          labels = paste0(substr(motif_f[j],newy,newy), "\U2192", 
                                                                          substr(final_motifs6[j],newy,newy)),
                                                          count = nrow(dataworky10),
                                                          county = c(rep(max(dataworky10$score), nrow(dataworky10)), rep(max(dataworky10$ipdRatio), nrow(dataworky10))))
            }
          }
          dataw <- do.call(rbind, deglist)
          valueme2 <- data.frame(valueme2 = 0.1 + max(as.numeric(as.vector(sapply(dataw[dataw$type=="ipdRatio",]$values, function (x) rep(x,3))))))
          dataworky <- cbind(dataw, valueme2)
          dataworky2 <- data.frame()
          dataworky2 <- dataworky[dataworky$base=="A",]
          
          write.csv(dataw, paste0(dir3,"/debug2_data.csv"), row.names = F)
          dummyline <- data.frame(type = c("ipdRatio", "score"), value = c(r_cut, 0))
          
          graphme666 <- ggplot(data=transform(dataw, plt_labels = paste("Position ", sprintf("%02d", pos), ", ", labels, sep="")),aes(x = base, y = values, color = type, group = base)) + 
            geom_violin() + 
            facet_grid(type ~ plt_labels, scales="free_y") + 
            theme_gray() %+replace% prettify + 
            geom_hline(data = dummyline, aes(yintercept = value),colour="#BB0000", linetype="dashed") + 
            geom_text(aes(y=county, label=count),size = 5, colour = "grey25") +
            ggtitle(paste0(motif_f[j]," to ",final_motifs[j],", ",max(dataw$pos)," positions"))
          filename666 <- paste0(dir3,"/",name,"_",as.character(motif_f[j]),"_debug_graph2.png")
          ggsave(graphme666, filename = filename666, width=length(unique(dataw$pos))*2.8, height=8, limitsize = F)
          
        },error=function(e){
          cat("   + Error: Not enough observations for this motif\n")
        })
        if(j2 == FALSE){ # there is error
          final_motifs[j] <- NA
          final_motifs6[j] <- NA
          control[j] <- NA
          contig_vec[j] <- NA
          if(file.exists(paste0(dir3, "/pos01A.csv"))){
            # WRITING TO DEBUG_GRAPH 2 NEW
            cat(paste("   + Writing to debug_graph2.png"), "\n")
            bases <- c("A","C","G","T")
            deglist<- list()
            prettify <- theme(panel.background = element_rect(fill = NA,color="gray"), 
                              panel.grid.major.y = element_blank(),
                              panel.grid.major.x = element_line(size=.1, color="black",linetype="dotted"), 
                              panel.grid.minor.y = element_blank(),
                              panel.grid.minor.x = element_line(size=.1, color="black"),
                              legend.position="bottom")
            
            
            for(newy in 1:nchar(motif_f[j])){
              for(newy2 in 1:4){
                locpos <- bases[newy2]
                namepos <- paste0(dir3, "/pos", sprintf("%02d",newy), locpos, ".csv")
                dataworky10 <- read.csv(namepos, sep = ",", header = T)
                deglist[[((newy-1)*4+newy2)]] <- data.frame(base = locpos, pos = newy, values = c(dataworky10$score, dataworky10$ipdRatio), 
                                                            type = c(rep("score", nrow(dataworky10)), rep("ipdRatio", nrow(dataworky10))),
                                                            labels = paste0(substr(motif_f[j],newy,newy), "\U2192", 
                                                                            substr(final_motifs6[j],newy,newy)),
                                                            count = nrow(dataworky10),
                                                            county = c(rep(max(dataworky10$score), nrow(dataworky10)), rep(max(dataworky10$ipdRatio), nrow(dataworky10))))
              }
            }
            dataw <- do.call(rbind, deglist)
            valueme2 <- data.frame(valueme2 = 0.1 + max(as.numeric(as.vector(sapply(dataw[dataw$type=="ipdRatio",]$values, function (x) rep(x,3))))))
            dataworky <- cbind(dataw, valueme2)
            dataworky2 <- data.frame()
            dataworky2 <- dataworky[dataworky$base=="A",]
            
            cat(dir3)
            write.csv(dataw, paste0(dir3,"/debug2_data.csv"), row.names = F)
            dummyline <- data.frame(type = c("ipdRatio", "score"), value = c(r_cut, 0))
            
            graphme666 <- ggplot(data=transform(dataw, plt_labels = paste("Position ", sprintf("%02d", pos), ", ", labels, sep="")),aes(x = base, y = values, color = type, group = base)) + 
              geom_boxplot() + 
              facet_grid(type ~ plt_labels, scales="free_y") + 
              theme_gray() %+replace% prettify + 
              geom_hline(data = dummyline, aes(yintercept = value),colour="#BB0000", linetype="dashed") +
              geom_text(aes(y=county, label=count), size = 5, colour = "grey25") +
              ggtitle(paste0(motif_f[j]," to ",final_motifs[j],", ",max(dataw$pos)," positions"))
            filename666 <- paste0(dir3,"/",name,"_",as.character(motif_f[j]),"_debug_graph2.png")
            ggsave(graphme666, filename = filename666, width=length(unique(dataw$pos))*2.8, height=8, limitsize = F)
            
          }}
        j <- j + 1
        cat(paste0("   + j incremented to ", j), "\n")
        rm(list=ls()[!ls()%in%allvarsbwhile])
      }
      #######################################################################################################
      pos1 <- motiftable$centerPos[!urder]
      pos2 <- motiftable$centerPos[urder]
      pos <- c(pos1, pos2)
      controlf <- c(control, rep("experim",length(moditype_c)))
      contig_vecf <- rep(contig_vec[1],length(pos))
      modis <- c(moditype, moditype_c)
      fractionsp <- c(perc_meth, perc_meth_c)
      modifs_b <- c(motif_f, motif_mod)
      modifs_a <- c(final_motifs, motif_mod)
      modifs_u <- c(final_motifs6, motif_mod)
      r_cutoffs <- c(r_cutoff, rep(0,length(moditype_c)))
      refyonp <- c(refyon,rep(0,length(perc_meth_c)))
      if(!length(modifs_b)==length(modifs_a)){
        while(!length(modifs_a) == length(modifs_b)){
          modifs_a[length(modifs_b)] <- NA
        }
      }
      if(!length(modifs_b)==length(modifs_u)){
        while(!length(modifs_u) == length(modifs_b)){
          modifs_u[length(modifs_b)] <- NA
        }
      }
      final_table <- data.frame(before = modifs_b, after = modifs_a, unshaved = modifs_u, modificationType = modis, 
                                centerPos = pos, fraction = fractionsp, refyon = refyonp, r_cutoffs = r_cutoffs)#, "nodataforpos" = nodataatpos
      cat(paste(" + Deduping for strain ", name),"\n")
      # DEDUPING
      ulgy <- vector()
      final_table1.a <- final_table[complete.cases(final_table[,"after"]),]
      final_table1.b <- final_table[!complete.cases(final_table[,"after"]),]
      final_motifs9 <- as.vector(final_table1.a$after)
      if(!length(final_motifs9)==0){
        #           for(ulg in 1:length(final_motifs9)){
        #             if(!final_motifs9[ulg] == ""){
        #               if(sum(grepl(final_motifs9[ulg],final_motifs9)) != 1){
        #                 ulgy[ulg] <- F
        #               }else{ulgy[ulg] <- T}
        #             }else{
        #               ulgy[ulg] <- T
        #             }
        #           }
        #           final_table1.c <- final_table1.a[ulgy,]
        #           final_table2 <- rbind(final_table1.c,final_table1.b)
        #           final_table.deduped <- final_table2#[!duplicated(final_table2[,c(2)]),]  # prune duplicate motifs
        final_table.deduped <- final_table
        dog2s <- do.call(rbind, dog2_list)
        positive <- do.call(rbind, positive_list)
        negative <- do.call(rbind, negative_list)
        
        # Writing to positive and negative
        filename97 <- paste0(dir4, name, "_neg.csv")
        filename98 <- paste0(dir4, name, "_pos.csv")
        cat("  - Writing to pos/neg.csv \n")
        write.csv(negative, filename97, row.names = F)
        write.csv(positive, filename98, row.names = F)
        
        # WRITE MOTIF_POSTMOTIFS.CSV
        filename3 <- paste0(dir4, name, "_postmotifs.csv")
        wfilename3 <- paste0(name, "_postmotifs.csv")
        cat(paste("  - Writing to", wfilename3), "\n")
        write.csv(final_table.deduped, filename3, row.names = F)
        
        # WRITE MOTIF_PREMOTIFS.CSV
        filename4 <- paste0(dir4, name, "_premotifs.csv")
        wfilename4 <- paste0(name, "_premotifs.csv")
        cat(paste("  - Writing to", wfilename4), "\n")
        write.csv(dog2s, filename4, row.names = F)
      }
      
      if(!dog2s == 0){
        # DATA FOR GRAPHS
        #         uniquerly <- unique(dog2s$motif_n)
        #         
        #         # DEBUG PLOT V.1
        #         if(!all(is.na(final_table$after))){
        #           cat(paste("  - Writing to debug_graph.png"), "\n")
        #           for(jagery in 1:length(uniquerly)){
        #             workwmelist <- list()
        #             mymotif <- uniquerly[jagery]
        #             gatter <- dog2s[dog2s$motif_n == mymotif,]
        #             gattersplit <- split(gatter,gatter$pos)
        #             numberme <- matrix(1:(NROW(gattersplit)*4), ncol=4, byrow = T)
        #             for(jagery2 in 1:NROW(gattersplit)){
        #               workwme <- gattersplit[[jagery2]]
        #               for(jagery3 in 1:4){
        #                 worky <- workwme[jagery3,]
        #                 baseref <- c("A","C","G","T")
        #                 workwmelist[[numberme[jagery2,jagery3]]] <- data.frame(base = baseref[jagery3],
        #                                                                        pos = jagery2,
        #                                                                        type = c("quant_s","quant_r","p_val"), 
        #                                                                        values = c(workwme[jagery3,"quant_s"],
        #                                                                                   workwme[jagery3,"quant_r"],
        #                                                                                   workwme[jagery3,"p_val"]),
        #                                                                        labels = workwme[jagery3,"totcount"])
        #               }
        #             }
        #             dataworky <- do.call(rbind, workwmelist)
        #             
        #             write.csv(dataworky, paste0(dir, "data/", uniquerly[jagery],"/debug_data.csv"), row.names = F)
        #             
        #             valueme2 <- data.frame(valueme2 = 0.1 + max(as.numeric(as.vector(sapply(dataworky[dataworky$type=="quant_r",]$values, function (x) rep(x,3))))))
        #             dataworky <- cbind(dataworky, valueme2)
        #             dataworky2 <- data.frame()
        #             dataworky2 <- dataworky[dataworky$type=="quant_r",]
        #             prettify <- theme(panel.background = element_rect(fill = NA,color="gray"), 
        #                               panel.grid.major.y = element_blank(),
        #                               panel.grid.major.x = element_line(size=.1, color="black",linetype="dotted"), 
        #                               panel.grid.minor.y = element_blank(),
        #                               panel.grid.minor.x = element_line(size=.1, color="black"),
        #                               legend.position="bottom")
        #             
        #             graphme1 <- ggplot(transform(dataworky, type=factor(type,levels=c("quant_r","quant_s","p_val"))), 
        #                                aes(x = base, y = values, color = type, group = type)) + 
        #               geom_line(size=1) + 
        #               facet_wrap(type ~ pos, scales="free", nrow = 3) + 
        #               theme_gray() %+replace% prettify + 
        #               geom_text(data=dataworky2, aes(y=valueme2, label=labels), aes=1, colour = "dark gray",size=3) + 
        #               ggtitle(paste0(uniquerly[jagery]," to ",final_motifs[jagery],", ",dataworky[nrow(dataworky),2]," positions"))
        #             
        #             filename91 <- paste0(dir, "data/", uniquerly[jagery],"/",name,"_",as.character(uniquerly[jagery]),"_debug_graph.png")
        #             ggsave(graphme1, filename = filename91, width=length(unique(dataworky$pos))*2.8, height=6, limitsize = F)
        #           }
        #         }
        
        #         # DEBUG V. 2 MEMORY INTENSIVE
        #         if(!all(is.na(final_table$after))){
        #           cat(paste("  - Writing to debug_graph2.png"), "\n")
        #           finalelist <- list()
        #           for(jagery in 1:length(uniquerly)){
        #             workwmelist <- list()
        #             mymotif <- uniquerly[jagery]
        #             gatter <- dog2s[dog2s$motif_n == mymotif,]
        #             gattersplit <- split(gatter,gatter$pos)
        #             #numberme <- matrix(1:(NROW(gattersplit)*4), ncol=4, byrow = T)
        #             for(jagery2 in 1:NROW(gattersplit)){
        #               workwme <- gattersplit[[jagery2]]
        #               deglist<- list()
        #               for(jagery3 in 1:length(workwme$motif_)){
        #                 locpos <- substr(workwme$motif_[jagery3], workwme$pos[jagery3],workwme$pos[jagery3])
        #                 namepos <- paste0(dir2, mymotif, "/pos", sprintf("%02d", workwme$pos[jagery3]), locpos, ".csv")
        #                 dataworky10 <- read.csv(namepos, sep = ",", header = T)
        #                 dataworky11 <- data.frame(base = locpos, pos = jagery2, values = c(dataworky10$score, dataworky10$ipdRatio), 
        #                                           type = c(rep("score", nrow(dataworky10)), rep("ipdRatio", nrow(dataworky10))),
        #                                           labels = paste0(substr(uniquerly[jagery],jagery2,jagery2), "\U2192", 
        #                                                           substr(final_motifs6[jagery],jagery2,jagery2)))
        #                 deglist[[jagery3]] <- dataworky11
        #               }
        #               deglist2 <- cbind(do.call(rbind, deglist))
        #               workwmelist[[jagery2]] <- deglist2
        #             }
        #             dataw <- do.call(rbind, workwmelist)
        #             valueme2 <- data.frame(valueme2 = 0.1 + max(as.numeric(as.vector(sapply(dataw[dataw$type=="ipdRatio",]$values, function (x) rep(x,3))))))
        #             dataworky <- cbind(dataw, valueme2)
        #             dataworky2 <- data.frame()
        #             dataworky2 <- dataworky[dataworky$base=="A",]
        #             
        #             write.csv(dataw, paste0(dir, "data/", uniquerly[jagery],"/debug2_data.csv"), row.names = F)
        #             
        #             prettify <- theme(panel.background = element_rect(fill = NA,color="gray"), 
        #                               panel.grid.major.y = element_blank(),
        #                               panel.grid.major.x = element_line(size=.1, color="black",linetype="dotted"), 
        #                               panel.grid.minor.y = element_blank(),
        #                               panel.grid.minor.x = element_line(size=.1, color="black"),
        #                               legend.position="bottom")
        #             
        #             r_cutg <- final_table[final_table$before==as.character(uniquerly[jagery]),]$r_cutoffs
        #             
        #             graphme666 <- ggplot(data=transform(dataw, plt_labels = paste("Position ", sprintf("%02d", pos), ", ", labels, sep="")),aes(x = base, y = values, color = type, group = base)) + 
        #               geom_boxplot() + 
        #               facet_grid(type ~ plt_labels, scales="free_y") + 
        #               theme_gray() %+replace% prettify + geom_hline(aes(yintercept = r_cutg),colour="#BB0000", linetype="dashed") +
        #               #geom_text(data=dataworky2, aes(y=valueme2, label=labels), aes=1, colour = "dark gray",size=3) + 
        #               ggtitle(paste0(uniquerly[jagery]," to ",final_motifs[jagery],", ",max(dataw$pos)," positions"))
        #             filename666 <- paste0(dir, "data/", uniquerly[jagery],"/",name,"_",as.character(uniquerly[jagery]),"_debug_graph2.png")
        #             ggsave(graphme666, filename = filename666, width=length(unique(dataw$pos))*2.8, height=8, limitsize = F)
        #           }
        #         }
        
        premotif <- dog2s
        allmotifs <- unique(premotif$motif_n)
        for(layo in 1:length(allmotifs)){
          suballmotifs <- premotif[(premotif$motif_n == allmotifs[layo]),]
          allpositions <- unique(suballmotifs$pos)
          for(layo2 in 1:length(allpositions)){
            suballpositions <- suballmotifs[(suballmotifs$pos == allpositions[layo2]),]
            baseone <- rep(c("A","C","G","T"),4)
            basetwo <- unlist(strsplit(stri_dup(c("A","C","G","T")[1:4],4),""))
            
            fileA <- paste0(dir, "data/", allmotifs[layo],"/pos", sprintf("%02d",layo2),"A.csv")
            fileC <- paste0(dir, "data/", allmotifs[layo],"/pos", sprintf("%02d",layo2),"C.csv")
            fileG <- paste0(dir, "data/", allmotifs[layo],"/pos", sprintf("%02d",layo2),"G.csv")
            fileT <- paste0(dir, "data/", allmotifs[layo],"/pos", sprintf("%02d",layo2),"T.csv")
            
            ipdA <- read.csv(fileA, header = TRUE, sep=",")$ipdRatio
            ipdC <- read.csv(fileC, header = TRUE, sep=",")$ipdRatio
            ipdG <- read.csv(fileG, header = TRUE, sep=",")$ipdRatio
            ipdT <- read.csv(fileT, header = TRUE, sep=",")$ipdRatio
            
            data=c(ipdA,ipdC,ipdG,ipdT)
            
            # grouping
            group = rep( c('A','C','G','T') , c(length(ipdA),length(ipdC),length(ipdG),length(ipdT)) )
            
            # get the pvalue for each pair
            pval = TukeyHSD(aov(data~group))$group[,4]
            
            valuedat <- as.vector(c(1,pval[1],pval[2],pval[3],pval[1],1,pval[4],pval[5],pval[2],pval[4],1,pval[6],pval[3],pval[5],pval[6],1))
            
            heatdata <- data.frame(Bases1 = baseone, 
                                   Bases2 = basetwo,
                                   value = valuedat)
            # heatplot <- qplot(x=Bases1, y=Bases2, data=heatdata, fill=value, geom="tile")
            heatplot <- ggplot(heatdata, aes(Bases1, Bases2)) + geom_tile(aes(fill=value), color="white") + 
              scale_fill_gradient(low = "aliceblue", high ="steelblue") + 
              geom_text(aes(fill = heatdata$value, label = signif(heatdata$value,4))) + 
              theme_gray()%+replace%prettify + theme(legend.position="none", axis.title.x=element_blank(),
                                                     axis.title.y=element_blank())
            
            filename66 <- paste0(dir, "data/", allmotifs[layo],"/pos", layo2, "_heatplot.png")
            ggsave(heatplot, filename = filename66, width= 4, height=4, limitsize = F)
            
          }
        }
      }
      
      # Time gather
      timegath <- do.call(rbind, timegat)
      write.csv(timegath, paste0(dir4, "times.csv"), row.names = F)
      
      # Stop cluster
      stopCluster(cl)
      
      # Stop connections
      showConnections()
      closeAllConnections()
    }
    
    # BREAK
     
    #cat(paste("Total time elapsed:",signif(as.numeric(elapsed[1])/60,3),"min."), "\n")
    
    cat("Combining tables into one csv")
    
    files2 <- list.files(pattern = "\\_postmotifs.csv$", recursive = T)
    names2 <- gsub("/", "", stri_extract_first_regex(files2, "[A-z]+/+[0-9]{6}"))
    myfiles <- lapply(files2, function(z)read.csv(z))
    logical <- unlist(lapply(myfiles, function(z) nrow(z) != 0))
    names2 <- names2[logical]
    if(any(logical == FALSE)){
      falsey <- which(logical==FALSE)
      for(ful in 1:length(falsey)){
        myfiles[[falsey[ful]]] <- NULL
      }
    }
    finale96 <- mapply(cbind, myfiles, "ID" = names2, SIMPLIFY=F)
    #finale96 <- finale95[sapply(finale95, function(x) dim(x)[1]) > 0]
    finale <- data.table(do.call(rbind, finale96))
    finale2 <- setcolorder(finale, c(ncol(finale), 1:(ncol(finale)-1)))
    filename5 <- "combined.csv"
    sameyon <- vector()
    b_unique <- vector()
    a_unique <- vector()
    for(hurly in 1:nrow(finale)){
      sameyon[hurly] <- ifelse(as.character(finale$after[hurly]) == as.character(finale$before[hurly]), 1, 0)
      bef <- as.vector(finale$before)
      bef[hurly] <- NA
      bef <- bef[!is.na(bef)]
      af <- as.vector(finale$after)
      af[hurly] <- NA
      af <- af[af != ""]
      af <- af[!is.na(af)]
      b_unique[hurly] <- ifelse(as.character(finale$before[hurly]) %in% bef,0,1)
      if((!is.na(as.character(finale$after[hurly]))) && (as.character(finale$after[hurly])!="")){
        if((as.character(finale$after[hurly]) %in% af)){
          a_unique[hurly] <- 0
        }else{
          a_unique[hurly] <- 1
        }
      }else{
        a_unique[hurly] <- NA
      }
    }
    finale <- data.frame(finale, same = sameyon,b_unique=b_unique,a_unique=a_unique)
    
    cat(paste("   - Writing", filename5, "to current directory"), "\n")
    write.csv(finale, paste0(outdir,"/",filename5), row.names = F)
    
    ###################GENERATE TABLES AND PLOTS##############################
    
    #     finale <- read.csv("combined.csv",header=T,sep=",") 
    #     
    #     
    #     # Fix finale table;fill in missing
    #     for(i in 1:length(unique(finale$ID))){
    #       myid <- as.character(unique(finale$ID)[i])
    #       boolme <- !(motiftable$motifString %in% finale[finale$ID == myid,]$before)
    #       if(nrow(finale[finale$ID == myid,]) != nrow(motiftable)){
    #         finale <- rbind(finale,data.frame(ID = myid, 
    #                                           before = motiftable$motifString[boolme],
    #                                           after = NA,
    #                                           unshaved = NA,
    #                                           modificationType = motiftable$modificationType[boolme],
    #                                           centerPos = motiftable$centerPos[boolme],
    #                                           fraction = motiftable$fraction[boolme],
    #                                           refyon = 0,
    #                                           r_cutoffs = NA,
    #                                           same = NA,
    #                                           b_unique = 0,
    #                                           a_unique =NA))
    #       }
    #     }
    #     
    #     #prep finale for table 1
    #     finale1 <- finale
    #     finale1$after <- as.character(finale1$after)
    #     finale1$unshaved <- as.character(finale1$unshaved)
    #     for(j in 1:nrow(finale1)){
    #       newnameme <- as.character(finale1$ID[j])
    #       newmotiflist <- as.vector(read.csv(paste0("Cdiff/",substr(newnameme,6,nchar(newnameme)),"/motif_summary.csv"),header=T,sep=",")$motifString)
    #       newmotifme <- as.character(finale1$after[j])
    #       newmotifme2 <- as.character(finale1$unshaved[j])
    #       if(!(is.na(newmotifme) || newmotifme == "") && (newmotifme %in% newmotiflist)){ # already in and is motif
    #         finale1$after[j] <- paste0(".",newmotifme)
    #         finale1$unshaved[j] <- paste0(".",newmotifme2)
    #       }else{
    #         if((is.na(newmotifme) || newmotifme == "") && (newmotifme %in% newmotiflist)){ # already in and isn't motif
    #           finale1$after[j] <- ".-"
    #           finale1$unshaved[j] <- ".-"
    #         }else{
    #           if(!(is.na(newmotifme) || newmotifme == "") && !(newmotifme %in% newmotiflist)){ # isn't in and is motif
    #             finale1$after[j] <- paste0(".+",newmotifme)
    #             finale1$unshaved[j] <- paste0(".+",newmotifme2)
    #           }else{
    #             finale1$after[j] <- ".0" # isn't in and isn't motif
    #             finale1$unshaved[j] <- ".0"
    #           }
    #         }
    #       }
    #     }
    #     
    #     #order
    #     motiftable <- motiftable[order(motiftable$motifString),]
    #     finale <- finale[order(finale$ID,finale$before),]
    #     finale1 <- finale1[order(finale1$ID,finale1$before),] #for table 1
    #     
    #     # Table 1
    #     table1 <- matrix(finale1$after, nrow(motiftable))
    #     table1u <- matrix(finale1$unshaved, nrow(motiftable))
    #     rownames(table1) <- motiftable$motifString
    #     colnames(table1) <- as.vector(unique(finale1$ID))
    #     rownames(table1u) <- motiftable$motifString
    #     colnames(table1u) <- as.vector(unique(finale1$ID))
    #     write.csv(as.data.frame(table1),"table1.csv")
    #     write.csv(as.data.frame(table1u),"table1_N.csv")
    #     
    #     # Table 2
    #     table2 <- matrix(finale$after, nrow(motiftable))
    #     rownames(table2) <- motiftable$motifString
    #     colnames(table2) <- as.vector(unique(finale1$ID))
    #     for(sag in 1:ncol(table2)){
    #       for(sag2 in 1:nrow(table2)){
    #         table2[sag2,sag] <- ifelse(rownames(table2)[sag2]==as.vector(table2[sag2,sag]),1,0)
    #       }
    #     }
    #     write.csv(as.data.frame(table2),"table2.csv")
    #     
    #     
    #     # Table 3
    #     table3 <- matrix(finale$unshaved, nrow(motiftable))
    #     rownames(table3) <- motiftable$motifString
    #     colnames(table3) <- as.vector(unique(finale$ID))
    #     referrer <- data.frame(sym = c("A","G","T","C","W","S","M","K","R","Y","B","D","H","V","N"),
    #                            bases = c("","","","","AT","CG","AC","GT","AG","CT","SKYCGT",
    #                                      "RWKAGT","MWYACT","MRSACG","VHDBYRKMSWACGT"))
    #     for(sag in 1:ncol(table3)){
    #       for(sag2 in 1:nrow(table3)){
    #         reference <- rownames(table3)[sag2]
    #         tester <- as.vector(table3[sag2,sag])
    #         bool <- vector()
    #         for(sag3 in 1:nchar(reference)){
    #           basern <- unlist(strsplit(reference,""))[sag3]
    #           testern <- unlist(strsplit(tester,""))[sag3]
    #           basern2 <- as.vector(referrer[referrer$sym == basern,]$bases)
    #           bool[sag3] <- grepl(testern, basern2)
    #         }
    #         
    #         table3[sag2,sag] <- ifelse(any(bool),1,0)
    #       }
    #     }
    #     write.csv(as.data.frame(table3),"table3.csv")
    #     
    #     # Table 4
    #     table4 <- matrix(finale$unshaved, nrow(motiftable))#table3
    #     rownames(table4) <- motiftable$motifString
    #     colnames(table4) <- as.vector(unique(finale$ID))
    #     referrer2 <- data.frame(sym = c("A","G","T","C","W","S","M","K","R","Y","B","D","H","V","N"),
    #                             bases = c("WMRDHVN","SKRBDVN","WKYBDHN","SMYBHVN","DHN","BVN",
    #                                       "HVN","BDN","DVN","BHN","","","","",""))
    #     for(sag in 1:ncol(table4)){
    #       for(sag2 in 1:nrow(table4)){
    #         reference <- rownames(table4)[sag2]
    #         tester <- as.vector(table4[sag2,sag])
    #         bool <- vector()
    #         for(sag3 in 1:length(reference)){
    #           basern <- unlist(strsplit(reference,""))[sag3]
    #           testern <- unlist(strsplit(tester,""))[sag3]
    #           basern2 <- as.vector(referrer2[referrer2$sym == basern,]$bases)
    #           bool[sag3] <- grepl(testern, basern2)
    #         }
    #         table4[sag2,sag] <- ifelse(any(bool),1,0)
    #       }
    #     }
    #     write.csv(as.data.frame(table4),"table4.csv")
    #     
    #     #merge table 2,3,4 into table 5
    #     mtable <- table2
    #     for(r in 1:nrow(table2)){
    #       for(c in 1:ncol(table2)){
    #         if(!is.na(table2[r,c]) && table2[r,c]=="1"){
    #           mtable[r,c] <- "same"
    #         }
    #         if(!is.na(table2[r,c]) && table3[r,c]=="1"){
    #           mtable[r,c] <- "under" #SWITCHED WITH OVER
    #         }
    #         if(!is.na(table2[r,c]) && table4[r,c]=="1"){
    #           mtable[r,c] <- "over"  #SWITCHED WITH UNDER
    #         }
    #         if(!is.na(table2[r,c]) && mtable[r,c]=="0"){
    #           mtable[r,c] <- "diff"
    #         }
    #         if(is.na(mtable[r,c])){
    #           mtable[r,c] <- "NA"
    #         }
    #       }  
    #     }
    #     write.csv(as.data.frame(mtable),"table5.csv")
    #     
    #     #merge into ggplot heatmap
    #     mtable2 <- melt(mtable)
    #     names(mtable2) <- c("x","y","color")
    #     color_palette <- c("#ffffff","#8dd3c7","#ffffb3","#bebada","#fb8072")
    #     gm2 <- ggplot(mtable2, aes(y, x, fill=color)) + geom_tile(color="black") +
    #       theme(axis.text.x = element_text(angle = 90),
    #             axis.text=element_text(size=10)) + 
    #       xlab("") + ylab("") +
    #       scale_fill_manual(values = color_palette, name = "")
    #     ggsave(gm2,filename="ggplot2.png",width=13,height=9,units="in",dpi=600)
    #     
    #     #ggplot1
    #     a <- c("A","C","G","T","W","S","M","K","R","Y","B","D","H","V","N")
    #     b <- c("A","C","G","T","(A|T)","(C|G)","(A|C)","(G|T)","(A|G)","(C|T)","(C|G|T)","(A|G|T)","(A|C|T)","(A|C|G)","(A|C|G|T)")
    #     
    #     registerDoMC(16)  # number of CPU cores  
    #     
    #     expdir <-  system("find . -maxdepth 2 | grep -E '[0-9]{6}$'", intern=T)
    #     expdatlist <- list()
    #     expdatlist <- foreach(expy=1:length(expdir)) %dopar% {
    #       specexpdir <- list.files(paste0(expdir[expy],"/data"))
    #       specexpdir <- as.vector(paste0(expdir[expy],"/data/",specexpdir[grepl("expdat2",specexpdir)]))
    #       dater2 <- lapply(lapply(specexpdir,fread), function(z) setkey(z,"ipdRatio"))
    #       
    #       if(all(lapply(dater2,nrow)==0)){
    #         datery <- data.table(matrix(nrow=0,ncol=3))
    #         setnames(datery,names(datery),c("Motif","Strain","ipdRatio"))
    #       }else{
    #         dater2_10 <- dater2[[10]]
    #         dater2_11 <- dater2[[11]]
    #         dater2[[10]] <-  dater2_11
    #         dater2[[11]] <-  dater2_10
    #         dater3 <- Map(cbind, dater2, "index" = c(1:NROW(dater2))) #IDK WHY CAAAA and CAAAAA switch
    #         dater4 <- lapply(dater3, function(z) 
    #           setkey(z,motif)[str_detect(motif,mgsub(a,b,motiftable$motifString[index]))])
    #         dater5 <- lapply(lapply(dater4,function(z) z[ ,names(z)[!names(z)%in%c("ipdRatio","index")] := NULL]), function(z) z[!duplicated(z), ]) 
    #         dater <- dater5[sapply(dater5, function(z) dim(z)[1]) > 0]
    #         modmotiftab <- motiftable[c(1:nrow(motiftable))%in%unique(do.call(rbind, dater)$index),]
    #         datery <- mapply(cbind, "Motif"=modmotiftab$motifString, "Strain" = namesg2[expy], dater, SIMPLIFY=F)
    #       }
    #       cat(paste0(expy,"."))
    #       datery
    #     } 
    #     
    #     for(expyr in 1:NROW(expdatlist)){
    #       if(unlist(lapply(expdatlist, function(z)all(is.na(z))))[expyr])
    #       {expdatlist[[expyr]]=NULL}
    #     }
    #     dftab5 <- rbindlist(lapply(expdatlist,rbindlist))
    #     write.csv(dftab5,"ggplot1dat.csv",row.names = F)
    #     
    #     gm <- ggplot(dftab5, aes(x=0, y=ipdRatio)) + 
    #       geom_boxplot(outlier.size = 0.5) +
    #       facet_grid(Motif ~ Strain) +
    #       theme_bw() + theme(strip.text.y = element_text(size = 4),
    #                          axis.text.x = element_blank()) +
    #       ylim(0, 8)
    #     ggsave(gm,filename="ggplot1.png", width = 50, height = 50, units = "in",limitsize=F)
    #     
    #     #ggplot3
    #     
    #     matg3v2 <- matrix(, nrow = nrow(table4), ncol = ncol(table4)) 
    #     dimnames(matg3v2) <- list(rownames(table4), colnames(table4))
    #     matg3 <- dftab5[, mean(ipdRatio), list(Motif, Strain)]
    #     for(rowi in 1:nrow(matg3)){
    #       matg3v2[which(rownames(matg3v2)==matg3$Motif[rowi]),
    #               which(colnames(matg3v2)==paste0("Cdiff0",matg3$Strain[rowi]))] <- matg3$V1[rowi]
    #     }
    #     clusterdat <- matg3v2
    #     dat983 <- melt(matg3v2)
    #     dat983[is.na(dat983)]<-0
    #     dat983v2 <- dat983
    #     setnames(dat983v2,names(dat983v2),c("Motif","Strain","ipdRatio")) 
    #     dat983v2$ipdRatio <- signif(dat983v2$ipdRatio,2)
    #     dat983v2$ipdRatio[dat983v2$ipdRatio>6] <- ">6"
    #     dat983$ipdRatio[dat983$ipdRatio>6] <- 6
    #     setnames(dat983,names(dat983),c("Motif","Strain","ipdRatio")) 
    #     gm983 <- ggplot(dat983, aes(x=Strain, y=Motif)) + 
    #       geom_tile(aes(fill=ipdRatio),colour = "black") +
    #       geom_text(aes(fill=dat983$ipdRatio, label=dat983v2$ipdRatio),size=1) + 
    #       theme(axis.text.x = element_text(angle = 90),
    #             axis.text=element_text(size=10)) + 
    #       xlab("") + ylab("") + scale_fill_gradient(low = "white", high = "red") 
    #     ggsave(gm983,filename="ggplot3.png",width=13,height=9,units="in",dpi=600)
    #     
    #     # Motif Strain IPDRatio
    #     # A         A1        1
    #     # A         A1        2
    #     # A         A1        3
    #     # A         B1        1
    #     # A         B1        2
    #     # B .....
    #     
    #     #ggplot4 cluster heatmap
    #     mtscaled <- as.matrix(scale(clusterdat))
    #     png(file = "ggplot4.png",width=13,height=9,units="in",res=1200)
    #     heatmap.2(mtscaled,margins=c(5,8),col=colorRampPalette(c('white', 'red'))(12))
    #     dev.off()
    
    
    
    ################################################################
  })
  
  # clean
  debugfolder <- paste0(outdir,"/debug_plots/")
  tarfile <- 'debug_plots.tgz'
  if(dir.exists(debugfolder) & file.exists(tarfile)){
    unlink(debugfolder, recursive = TRUE)
    unlink(tarfile)
  }
  
  # Create zip file of all debug plots
  cat(paste("   - Gzipping all debug_plots to debug_plots.tgz"), "\n")
  debuglist <- list.files(pattern = "\\_debug_graph.*.png$", recursive = TRUE)
  dir.create(debugfolder, showWarnings = FALSE)
  file.copy(debuglist, debugfolder)
  tar(tarfile,debugfolder,compression='gzip')
  
  # BREAK AND END
   
  #cat(paste("Total time elapsed:",signif(as.numeric(elapsed[1])/60,3),"min."), "\n")
}))