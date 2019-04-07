#' Plot CNV segmentation with log2 ratio of each bin in the background;
#' can zoom in to specified genomic region per chromosome
#' For help: plot.l2r.ms_zoom.R -h
#'
#' Author: Maxwell Sherman & Semin Lee (?)
#' Revised: 08-15-16
#' Revision notes: generic fns from R/fnUtils.R
Args <- commandArgs()

bin_dir = Args[6]
seg_file = Args[7]
bin_size = Args[8]
mouse = Args[9]
title = Args[10]

bin.dir.tmr = bin_dir
seg.file    = seg_file
bin.size    = bin_size
lamda       = 3
mouse       = mouse
title       = title
loc       = NULL
chrm        = 'TRUE'
xy          = TRUE
k           = 5 # integer width of median window, must be odd


#library("Cairo")

chrFromSpecies <- function(species) {
    # Returns chromosome numbers based on species
    #
    # Args:
    #   species     str     species for which to retrieve chrm numbers
    #
    # Returns:
    #   chrl        vector  chromosome numbers
    if (species == 'human') {
        chrl <- c(1:22, "X")
    } else if (species == 'mouse') {
        chrl <- c(1:19, "X")
    } else {
        stop("species not recognized. Don't know what to use for chromosomes")
    }
    names(chrl) <- chrl
    return(chrl)
}

loadBinData <- function(bin.dir.tmr, loc=NULL, species='human', chrm=TRUE) {
    # Plots CNV segments with log2 ratio of each bin in the backgroun
    #
    #  Args:
    #      bin.dir.tmr     str     path to directory of bin data
    #      loc             df      location: chrom   start   end
    #      species         str     species
    #      chrm            bool    does file name include "chrm"
    #
    #  Returns:
    #      bins            df      chr start   end log2ratio   color   subject

    write(sprintf("Processing: %s", bin.dir.tmr), stderr())

    if (! is.null(loc)) {
        chrl <- as.character(loc$chrom)
    } else {
        chrl <- chrFromSpecies(species)
    }
    names(chrl) <- chrl
    # Read bin files
    # start   end     obs     expected        var
    # 61762   78721   12      16.35           72.4539
    cntl.tmr <- lapply(chrl, function(c) {
           if (suppressWarnings(all(is.na(as.numeric(c))))) {
               if (chrm) {
                    pattern = sprintf("[.]chrm_%s[.]", c)
               } else {
                    # pattern = sprintf("-%s[.]", c)
                    pattern = sprintf("[.]%s[.]", c)
               }
           } else {
               if (chrm) {
                    pattern = sprintf("[.]chrm_%02d[.]", as.numeric(c))
               } else {
                    # pattern = sprintf("-%02d[.]", as.numeric(c))
                    pattern = sprintf("[.]%02d[.]", as.numeric(c))
               }
           }
           print(pattern)
           file <- file.path(bin.dir.tmr,
                             list.files(path=bin.dir.tmr, pattern=pattern)[1])
           print(file)

           # cnt <- scan(file, list(1,1,1,1.1,1.1), sep="\t", skip=1)
           # names(cnt) <- c("start", "end", "obs", "expected", "var")
           cnt <- read.table(file, sep="\t", header=T, stringsAsFactors=F)
           return(cnt)
        }  # End function
    )      # End lapply
    # Total tumor and normal cnt
    # In BICseq2, observed counts treated as tumor, expected as normal
    tc <- 0
    nc <- 0
    # Get the ratio of total tumor and normal read counts
    tc.chr <- sapply(chrl, function(c) {sum(cntl.tmr[[c]]$obs)})
    nc.chr <- sapply(chrl, function(c) {sum(cntl.tmr[[c]]$expected)})
    tc <- sum(as.numeric(tc.chr))
    nc <- sum(as.numeric(nc.chr))
    # Calculate [1/(tc/nc)] for entire genome for normalization purposes
    R <- nc / tc
    # Create data frame of bin edges and copy number ratio
    # chrom    start    end    ratio
    #   1          0    10000   0
    #   1      10001    20000   0.5
    bins <- data.frame(do.call(rbind,
                              lapply(chrl,
                                     function(c) {
                                       r <- (cntl.tmr[[c]]$obs + 0.00000001) /
                                                (cntl.tmr[[c]]$expected + 0.0000001) * R
                                       ec <- sub("chr", "", c)
                                       if (ec[1] != 'X') {
                                            label <- sprintf("%02d", as.numeric(ec))
                                       } else {
                                            label <- as.character(ec)
                                       }
                                       df <- data.frame(ec, cntl.tmr[[c]]$start,
                                                    cntl.tmr[[c]]$end, log2(r), label)

                                       if (! is.null(cntl.tmr[[c]]$gc)) {df$gc <- cntl.tmr[[c]]$gc}
                                       return(df)
                                     }  # End function
                                    )   # End lapply
                             )          # End do.call
                      )                 # End data.frame

    colnames(bins)[1:5] <- c("chrom", "start", "end", "ratio", "chrm.label")
    bins$chrom <- as.character(bins$chrom)
    bins$chrm.label <- as.character(bins$chrm.label)
    bins$copy <- rep(0, nrow(bins))
    bins$copy[bins$ratio > 0.3] <- 1
    bins$copy[bins$ratio < -0.3] <- -1
    file <- list.files(bin.dir.tmr)[1]
    subject <- tail(unlist(strsplit(bin.dir.tmr, '/', fixed=TRUE)), n=2)[1]
    print(subject)
    bins$subject <- subject
    if (! is.null(loc)) {
        bins <- bins[(bins$chrom==loc$chrom) & (bins$start>loc$start) & (bins$end<loc$end), ]
    }
    write(sprintf("Done processing: %s", bin.dir.tmr), stderr())
    return(bins)
}

plot.l2r <- function(bin.dir.tmr, seg.file, bin.size, lamda,
                     mouse, title, xy=1, endrule="median", k=3, loc=NULL, chrm=TRUE) {
    ## Plots CNV segments with log2 ratio of each bin in the backgroun
    ##
    ##  Args:
    ##      bin.dir.tmr     str     path to directory of bin data
    ##      seg.file        str     path to segmentation file
    ##      bin.size        int     size of bins
    ##      lamda           int     lambda smoothing factor
    ##      mouse           bool    is this mouse data?
    ##      title           str     title of plot
    ##      xy              bool    plot X chrm?
    ##      endrule         str     what to do with endpoints
    ##      k               int     number of neighboring points to include in
    ##                              running median
    ##
    ##  Returns:
    ##      --              --      Saves plot

    if (mouse == "0") {
        species <- "mouse"
    } else {
        species <- "human"
    }
    bins <- loadBinData(bin.dir.tmr, loc=loc, species=species, chrm=chrm)
    segments <- read.table(seg.file, sep="\t", header=T, stringsAsFactors=F)
    #chrom   start   end     binNum  observed        expected        log2.copyRatio
    #chr1    10027   990262  216     7163    4923.08 0.153223364203753
    segments <- segments[,c(1,2,3,7)]
    colnames(segments) <- c("chrom", "start", "end", "ratio")
    segments$chrom <- sub("chr", "", segments$chrom)
    # Plot
    pngf <- plotRatioNSeg(seg.file, bins, segments, title, bin.size, lamda, endrule=endrule, k=k, loc=loc)
    print(paste("done plotting", pngf))
}

plotRatioNSeg <- function(seg.file, bins, segments, sampleName=NULL, bin.size, lamda,
                         save = TRUE, endrule = c("none", "median", "keep", "constant"),
                         k = 3, indexOnly = FALSE, loc=NULL) {

    ## Actual plotting routine to plot CNV segs and bin ratios
    endrule     <- match.arg(endrule)
    graphList   <- list()
    locations   <- alignGenes(bins[, c("chrom", "start")], indexOnly = indexOnly)
    adjustments <- getAdjustments(bins[, c("chrom", "start")], indexOnly = indexOnly)
    num.segs    <- dim(segments)[1]
    print(adjustments)

    if (is.null(loc)) {
        xlim <- c(0, max(locations) + 10)
    } else {
        xlim <- c(adjustments[loc$chrom] + loc$start, adjustments[loc$chrom] + loc$end)
    }
    print(loc)
    print(xlim)

    if (save) {
        seg.dir <- dirname(seg.file)
        graphList[[sampleName]] <- file.path(seg.dir, paste(c(head(strsplit(basename(seg.file), ".", fixed=T)[[1]], n=-1), "l2r.pdf"), collapse="."))
        pdf(file = graphList[[sampleName]], width = 18, height = 5)
    }
    # Title generation
    main <- paste(sampleName, " ( b", bin.size, " : k", k, " : l", lamda,
                  " : s", num.segs, " )", sep="")
    # Create plot and axis options
    graphics::plot(0, 0, type="n", main=main, xlab="Chromosome", ylab=" Log2 ratio",
                   ylim=c(-2.2, 2.2), axes=FALSE, xlim=xlim, cex.lab=1.2)

    # Create axis (on the left of plot) and draw a box arond it
    axis(2)
    box()

    # Add colored boxes around odd numbered chromosomes
    HighLightChrom(adjustments, -2, 2)
    # Plot bin ratios
    if (endrule != "none") {
        nas <- which(is.na(as.numeric(bins[, "ratio"])))
        if (length(nas) != 0) {
            points(locations[-nas], runmed(as.vector(bins[-nas, "ratio"]),
                                           k=k, endrule=endrule),
                   cex = 0.3, pch = 16, col='grey')
        } else {
            points(locations, runmed(as.vector(bins[, "ratio"]), k = k, endrule = endrule),
                   cex = 0.3, pch = 16, col='black')
        }
    } else {
        points(locations, as.numeric(as.vector(bins[, "ratio"])),
               cex=0.3, pch=16, col='grey80')
    }

    # Add lines at -1, 0, +1
    lines(c(min(locations), max(locations)), rep(0, 2), lwd = 2, col = "blue")
    lines(c(min(locations), max(locations)), rep(1, 2), col = "blue")
    lines(c(min(locations), max(locations)), rep(-1, 2), col = "blue")

    # Plot segments

    if (!indexOnly) {
        drawSegs(segments, adjustments)  # This is default
    } else {
        for (index in 1:nrow(temp)) {
            positions = range(locations[which(as.character(bins[, "chrom"]) ==
                                              as.character(temp[index, "chrom"]) &
                                              as.numeric(bins[, "start"]) <
                                              as.numeric(segments[index, "end"]) &
                                              as.numeric(bins[, "end"]) >
                                              as.numeric(segments[index, "start"]))])
            lines(c(positions[1], positions[2]), rep(segments[index, "ratio"], 2), col = "red", lwd = 3)
        }
    }

   # # Label odd chromosomes
    #markChrom(adjustments, -5)
    markChrom(adjustments, -2)

    # Save plot
    if(save){
        dev.off()
    }
    if(save){
        return(graphList)
    }else{
        return(invisible())
    }
}

alignGenes = function(positions, indexOnly = FALSE){
    if(indexOnly){
        return(1:nrow(positions))
    }else{
        adjustments = getAdjustments(positions, indexOnly)
        for(chrom in names(adjustments)){
            positions[positions[, 1] == chrom, 2] =
                as.numeric(positions[positions[, 1] == chrom, 2]) +
                as.numeric(adjustments[chrom])
        }
        return(as.numeric(positions[, 2]))
    }
}

# This function gets the values that can be used to adjust chromosomal
# locations to align genes on different chromosomes so that they  appear
# in sequence from chromosome one to Y along a single strand
getAdjustments = function(positions, indexOnly = FALSE){
    if(indexOnly){
        temp = split.data.frame(positions, factor(positions[, 1]))
        temp = unlist(lapply(temp, FUN = function(x) nrow(x) + 1))
    }else{
        temp = split.data.frame(positions, factor(positions[, 1]))
        temp = unlist(lapply(temp, FUN = function(x) max(as.numeric(x[, 2]))))
    }
    # Chromosomes to 30 for other organisms (e. g. fish - 25)
    chroms = sort(as.numeric(names(temp)[names(temp) %in% 1:30]))
    if(any(names(temp) %in% "X")){
        chroms = c(chroms, "X")
    }
    if(any(names(temp) %in% "Y")){
        chroms = c(chroms, "Y")
    }
    adjustments = 0
    for(index in 1:length(chroms) - 1) {
        adjustments = c(adjustments, (adjustments[length(adjustments)]+
                                       temp[as.character(chroms[index])]))
    }
    names(adjustments) = chroms
    return(adjustments)
}

HighLightChrom <- function(adjustments, min, max) {
    for (index in 1:length(adjustments)) {
        if (index %% 2 == 1) {
            polygon(c(adjustments[index], adjustments[index + 1],
                      adjustments[index + 1], adjustments[index]),
                    c(min, min, max, max), col = "wheat1", border = "white")
        }
    }
    return(invisible())
}

drawSegs = function(segdata, adjustments, seqOnly = TRUE){
    drawSegLine = function(segLocs){
        toAdd = adjustments[as.character(as.vector(segLocs["chrom"]))]
        lines(c(as.numeric(segLocs["start"]) + toAdd,
                as.numeric(segLocs["end"]) + toAdd),
              rep(segLocs["ratio"], 2), col = "red", lwd = 3)
    }
    junck = apply(segdata, 1, FUN = drawSegLine)
    return(invisible())
}

markChrom = function(adjustments, min){
    chromLocs = NULL
    chromNames = NULL
    for(i in 1:length(adjustments) - 1){
        if(i %% 2 == 1){
            chromLocs = c(chromLocs, mean(c(adjustments[i], adjustments[i + 1])))
            chromNames = c(chromNames, names(adjustments)[i])
        }
    }
    text(chromLocs, rep(min - 0.125, length(chromLocs)), chromNames, cex = 1.2)
    return(invisible())
}

ParseLocation <- function(loc.str) {
    if (is.character(loc.str)) {
        spl.chr <- strsplit(loc.str, ':')
        chrom <- as.numeric(spl.chr[[1]][1])

        spl.loc <- strsplit(spl.chr[[1]][2], '-')
        start <- as.numeric(spl.loc[[1]][1])
        end <- as.numeric(spl.loc[[1]][2])

        print('returning df')
        return(data.frame(chrom, start, end))
    } else {
        print('returning null')
        return(NULL)
    }
}


plot.l2r(bin.dir.tmr, seg.file, bin.size, lamda, mouse, title, xy=xy, k=k, loc=loc, chrm=chrm)
