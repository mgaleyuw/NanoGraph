#' Read in metadata.txt
#'
#'Reads in metadata.txt produced by NanoR and returns a Data.Frame object
#'Optionally, reads below a quality threshold can be dropped. By default all
#'reads are returned
#'
#'@param infile Path to the input file
#'@param qualityCutoff reads with quality scores >= to this cutoff are kept
#'@return a DataFrame of the infile
#'@export
loadTable <- function(infile, qualityCutoff=0){
  df <- utils::read.table(infile, sep="\t", header=TRUE, row.names=1)
  return(df[which(df$Quality >= 10),c("Relative.Time", "Length.of.Read", "Quality")])
}




#' Calculate read length distribution and coverage then save to .csv
#' 
#' Uses metadata DataFrame to calculate read distribution statistics and saves
#' to file. Optionally append results to existing file where row 
#' name = sample name.
#' 
#' @param indf DataFrame of read metadata
#' @param outputname filename (.csv) for export, including path is not saved to
#' working directory
#' @param append TRUE or FALSE to append summary statistics to an existing file
#' @return a vector containing results
#' @export
createTextSummary <- function(indf, outputname=NA, append=FALSE){
  #basic read length statistics
  quantiles <- summary(indf$Length.of.Read)
  
  #N50,N90
  N50 <- calcN50(indf)
  N90 <- calcN50(indf, nfactor=0.9)
  
  #coverage based on HGCh38.p13
  p13Tot=3272116950
  p13TotNonN=3110748599
  readSum <- sum(indf$Length.of.Read)
  coverageAll <- readSum/p13Tot
  coverageNonN <- readSum/p13TotNonN
  
  #The other definition of coverage can't be calculated pre alignment https://www.ecseq.com/support/ngs/how-to-calculate-the-coverage-for-a-sequencing-experiment
  #still need to look into how to use RSamtools on an index file, for this problem and also for filtering out only T-LRS reads.
  
  if(! is.na(outputname)){
    writeline=paste(coverageAll, coverageNonN, quantiles["Min."], quantiles["Max."], quantiles["Mean"], quantiles["Median"],  N50, N90,sep="\t")
    if(append){
      write(writeline,file=outputname.txt,append=TRUE)
    } else {
      outfile=file(outputname, "w")
      header=paste("coverage", "coverageNonN", "minRead", "maxRead", "meanRead", "medianReadLength",  "N50", "N90", sep="\t")
      writeLines(header, con=outfile)
      writeLines(writeline, con=outfile)
      close(outfile)
    }
  } else {
    return( c("coverage"=coverageAll, "coverageNonN"=coverageNonN, "minRead"=quantiles["Min."], "maxRead"=quantiles["Max."], "meanRead"=quantiles["Mean"],  "medianReadLength"=quantiles["Median"], "N50"=N50, "N90"=N90) )
  }
}

#'Calculate N50 for cumulative reads
#'
#'Uses metadata DataFrame and returns the N50 or any other cumulative percentile
#'
#'@param indf DataFrame of read metadata
#'@param nfactor cumulative percentile, default=0.5
#'@return a tibble of results
#'@export
calcN50 <- function(indf, nfactor=0.5){
  totalLength <- sum(indf$Length.of.Read)
  halfLength <- totalLength*nfactor
  summdf <- indf %>% dplyr::arrange(desc(Length.of.Read)) %>% dplyr::mutate("CumulativeKB" = cumsum(Length.of.Read)) %>% summarise(minL=nth(Length.of.Read, which.min(abs(CumulativeKB-halfLength))))
  return(as.vector(summdf[1]))
}

#'Plot bar chart of read counts or cumulative read KB in bins
#'
#'Uses metadata DataFrame to plot the distribution of reads in ggplot. Can plot
#'read quality in stacked bars, manipulate the width of bins.
#'
#'@param indf DataFrame of read metadata
#'@param outname Filename for plot to save. Can be .pdf or .png. Set to NA to 
#'skip plot saving.
#'@param binNum Number of Bins
#'@param counts TRUE or FALSE, default FALSE, plot counts, not cumulative KB
#'@param skinnyBin TRUE or FALSE, default FALSE, make bars skinny for overlay
#'@param outwidth width dimension to save plot
#'@param outheight height dimension to save plot
#'@param outunits unit of width and height, default 'in'
#'@param plotQual TRUE or FALSE, default FALSE, plot quality distribution as fill
#'@param qual.bins a vector of quality cutoffs for use with plotQual
#'@param barColor an r color or color hex code
#'@param indf2 DataFrame of read metadata for a second sample
#'@param colorMapVar name of a column in df to map fill to
#'@param dodgeAdjust a number to adjust the spacing of color mapped bars. Default is 10kb.
#'@param sampleNames a vector of samplenames for plotting two datasets
#'@param binWidth a number, representing width of bins. overrides skinnybin.
#'@return a ggplot2 plot
#'@export
plotReadCountDistributions <- function(indf, outname, binNum=50, counts=FALSE, skinnyBin=FALSE, outwidth=8, outheight=7, outunits="in", plotQual=FALSE, qual.bins=c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60), barColor="#8ec462", indf2=NA, barColorRamp=c("#8ec462", "firebrick1"), dodgeAdjust=10000, sampleNames=c("Sample 1", "Sample 2"), binWidth=NA) {
  ggplot2::theme_set(ggplot2::theme_classic(base_size=12))
  if(counts){
    bw=10000
    if(skinnyBin){
      bw=5000
    }
    if(!is.na(binWidth)){
      bw=binWidth
    }
    if(is.na(indf2)){
      g <- plotCountBases(indf, binNum, plotQual=plotQual, qual.bins=qual.bins, binWidth=bw, barColor=barColor)
    } else {
      indf$Dataset <- sampleNames[1]
      indf2$Dataset <- sampleNames[2]
      df <- rbind(indf, indf2)
      g <- plotCountBases(df, binNum, plotQual=FALSE, qual.bins=qual.bins, binWidth=bw, barColorRamp = barColorRamp, dodgeAdjust = dodgeAdjust, colorMapVar = "Dataset")
    }
  } else {
    bw=5000
    if(skinnyBin){
      bw=2500
    }
    if(!is.na(binWidth)){
      bw=binWidth
    }
    if(is.na(indf2)){
    g <- plotCumulativeBases(indf, binNum, plotQual=plotQual, qual.bins=qual.bins, binWidth=bw, barColor=barColor)
    } else {
     indf$Dataset <- sampleNames[1]
     indf2$Dataset <- sampleNames[2]
     df <- rbind(indf, indf2)
     g <- plotCumulativeBases(df, binNum, plotQual=FALSE, qual.bins=qual.bins, binWidth=bw, barColorRamp = barColorRamp, dodgeAdjust = dodgeAdjust, colorMapVar = "Dataset")
     
    }
    }
  if(! is.na(outname)){
    ggplot2::ggsave(outname, dpi=300, width=outwidth, height=outheight, units=outunits)
  }
  return(g)
}


#'Plot histogram of read counts
#'
#'Uses metadata DataFrame to plot the distribution of reads in ggplot. Can plot
#'read quality in stacked bars, manipulate the width of bins.
#'
#'@param df DataFrame of read metadata
#'@param binNum Number of Bins
#'@param plotQual TRUE or FALSE, default FALSE, plot quality distribution as fill
#'@param qual.bins a vector of quality cutoffs for use with plotQual
#'@param binWidth the width of each bar in b
#'@param barColor the name or hex string of a color if not plotting quality
#'@param barColorRamp a vector of rcolors to be used when mapping two datasets
#'@param colorMapVar name of a column in df to map fill to
#'@param dodgeAdjust a number to adjust the spacing of color mapped bars. Default is 10kb.
#'@return a ggplot2 plot
#'@export
plotCountBases <- function(df, bin.num, plotQual=TRUE, qual.bins=c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60), binWidth=8000, barColor="#8ec462", barColorRamp=c("#8ec462", "violet"), colorMapVar=NA, dodgeAdjust=10000){
  bins <- seq(from=min(df$Length.of.Read), to=max(df$Length.of.Read), length.out=bin.num)
  df$bincode <- .bincode(df$Length.of.Read, bins, right=FALSE)
  if(plotQual==FALSE){
    if(is.na(colorMapVar)){
      df.countKB <- df %>% dplyr::group_by(bincode) %>% dplyr::summarize(Count=n())
      df.countKB <- data.frame(df.countKB)
      
      names(bins) <- seq(1, length(bins))
      df.countKB$bin.label <- bins[as.vector(df.countKB$bincode)]
      
      
      g <- ggplot2::ggplot(df.countKB, ggplot2::aes(x=bin.label, y=Count)) + ggplot2::geom_bar(stat="identity", width=binWidth, fill=barColor) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) + ggplot2::scale_x_continuous(n.breaks=20, breaks=waiver(),  label=scales::label_bytes("kB")) + ggplot2::xlab("Read Length") + ggplot2::ylab("Read Count") + ggplot2::scale_y_continuous(label= scales::label_bytes("MB"))
    } else {
      df.countKB <- df %>% dplyr::group_by(across({{ colorMapVar }}), bincode) %>% dplyr::summarize(Count=n())
      df.countKB <- data.frame(df.countKB)
      
      names(bins) <- seq(1, length(bins))
      df.countKB$bin.label <- bins[as.vector(df.countKB$bincode)]
      
      
      g <- ggplot2::ggplot(df.countKB, ggplot2::aes_string(x="bin.label", y="Count", fill=colorMapVar)) + ggplot2::geom_bar(stat="identity", width=binWidth, position=position_dodge(width=dodgeAdjust)) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) + ggplot2::scale_x_continuous(n.breaks=20, breaks=waiver(),  label=scales::label_bytes("kB")) + ggplot2::xlab("Read Length") + ggplot2::ylab("Read Count") + ggplot2::scale_y_continuous(label= scales::label_bytes("MB")) + scale_fill_manual(values=barColorRamp)
      
    }
  } else {
    df$Quality.Bin <- factor(.bincode(df$Quality, qual.bins, right=FALSE))
    df$Quality.Bin <- factor(df$Quality.Bin, levels=sort(as.numeric(levels(df$Quality.Bin)), decreasing=TRUE))
    levels(df$Quality.Bin) <- paste(">", qual.bins[as.numeric(levels(df$Quality.Bin))])
    
    df.countKB <- df %>% dplyr::group_by(bincode, Quality.Bin) %>% dplyr::summarize(Count=n())
    df.countKB <- data.frame(df.countKB)
    
    names(bins) <- seq(1, length(bins))
    df.countKB$bin.label <- bins[as.vector(df.countKB$bincode)]
    
    g <- ggplot2::ggplot(df.countKB, ggplot2::aes(x=bin.label, y=Count, fill=Quality.Bin)) + ggplot2::geom_bar(stat="identity", width=binWidth) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) + ggplot2::scale_x_continuous(n.breaks=20, breaks=waiver(),  label=scales::label_bytes("kB")) + ggplot2::xlab("Read Length") + ggplot2::ylab("Read Count") + ggplot2::scale_y_continuous(label= scales::label_bytes("MB")) + ggplot2::scale_fill_viridis_d(name="Quality", direction=-1)
  }
}


#'Plot bar chart of cumulative read count
#'
#'Uses metadata DataFrame to plot the distribution of reads in ggplot. Can plot
#'read quality in stacked bars, manipulate the width of bins.
#'
#'@param df DataFrame of read metadata
#'@param bin.num Number of Bins
#'@param plotQual TRUE or FALSE, default FALSE, plot quality distribution as fill
#'@param qual.bins a vector of quality cutoffs for use with plotQual
#'@param binWidth the width of each bar in b
#'@param barColor the name or hex string of a color if not plotting quality
#'@param barColorRamp a vector of rcolors to be used when mapping two datasets
#'@param colorMapVar name of a column in df to map fill to
#'@param dodgeAdjust a number to adjust the spacing of color mapped bars. Default is 10kb.
#'@return a ggplot2 plot
#'@export
plotCumulativeBases <- function(df, bin.num, plotQual=TRUE, qual.bins=c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60), binWidth=6000, barColor="#8ec462", barColorRamp=c("#8ec462", "firebrick1"), colorMapVar=NA, dodgeAdjust=10000){
  bins <- seq(from=min(df$Length.of.Read), to=max(df$Length.of.Read), length.out=bin.num)
  df$bincode <- .bincode(df$Length.of.Read, bins, right=FALSE)
  
  if(plotQual==FALSE){
    if(is.na(colorMapVar)){
      df.cumuKB <- df %>% dplyr::group_by(bincode) %>% dplyr::mutate("CumulativeKB" =  cumsum(Length.of.Read)) %>% dplyr::group_by(bincode) %>% dplyr::slice_max(CumulativeKB)
      df.cumuKB <- data.frame(df.cumuKB)
      
      names(bins) <- seq(1, length(bins))
      df.cumuKB$bin.label <- bins[as.vector(df.cumuKB$bincode)]
      g <- ggplot2::ggplot(df.cumuKB, ggplot2::aes(x=bin.label, y=CumulativeKB)) + ggplot2::geom_bar(stat="identity", width=binWidth, fill=barColor) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) + ggplot2::scale_x_continuous(n.breaks=20, breaks=waiver(),  label=scales::label_bytes("kB")) + ggplot2::xlab("Read Length") + ggplot2::ylab("Number of Bases") + ggplot2::scale_y_continuous(label= scales::label_bytes("MB"))
    } else {
      
      df.cumuKB <- df %>% group_by(across({{ colorMapVar }}), bincode) %>% dplyr::mutate("CumulativeKB" =  cumsum(Length.of.Read)) %>% dplyr::group_by(across({{ colorMapVar }}), bincode) %>% dplyr::slice_max(CumulativeKB)
      df.cumuKB <- data.frame(df.cumuKB)
      
      names(bins) <- seq(1, length(bins))
      df.cumuKB$bin.label <- bins[as.vector(df.cumuKB$bincode)]
      
      g <- ggplot2::ggplot(df.cumuKB, ggplot2::aes_string(x="bin.label", y="CumulativeKB", fill=colorMapVar)) + ggplot2::geom_bar(stat="identity", width=binWidth, position = ggplot2::position_dodge(width = dodgeAdjust)) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) + ggplot2::scale_x_continuous(n.breaks=20, breaks=waiver(),  label=scales::label_bytes("kB")) + ggplot2::xlab("Read Length") + ggplot2::ylab("Number of Bases") + ggplot2::scale_y_continuous(label= scales::label_bytes("MB")) + ggplot2::scale_fill_manual(values=barColorRamp)
    }
  } else { 
    df$Quality.Bin <- factor(.bincode(df$Quality, qual.bins, right=FALSE))
    df$Quality.Bin <- factor(df$Quality.Bin, levels=sort(as.numeric(levels(df$Quality.Bin)), decreasing=TRUE))
    levels(df$Quality.Bin) <- paste(">", qual.bins[as.numeric(levels(df$Quality.Bin))])
    
    
    df.cumuKB <- df %>% dplyr::group_by(bincode, Quality.Bin) %>% dplyr::mutate("CumulativeKB" = cumsum(Length.of.Read)) %>% dplyr::group_by(bincode, Quality.Bin) %>% dplyr::slice_max(CumulativeKB)
    
    names(bins) <- seq(1, length(bins))
    df.cumuKB$bin.label <- bins[as.vector(df.cumuKB$bincode)]
    
    
    g <- ggplot2::ggplot(df.cumuKB, ggplot2::aes(x=bin.label, y=CumulativeKB, fill=Quality.Bin)) + ggplot2::geom_bar(stat="identity") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) + ggplot2::scale_x_continuous(n.breaks=20, breaks=waiver(), label=scales::label_bytes("kB")) + ggplot2::xlab("Read Length") + ggplot2::ylab("Number of Bases") + ggplot2::scale_fill_viridis_d(name="Quality", direction=-1) + ggplot2::scale_y_continuous(label= scales::label_bytes("MB"))
  }
  
  return(g)
}