#' Read in metadata.txt
#'
#'Reads in metadata.txt produced by NanoR and returns a Data.Frame object
#'Optionally, reads below a quality threshold can be dropped. By default all
#'reads are returned
#'
#'@param infile Path to the input file
#'@param qualityCutoff reads with quality scores >= to this cutoff are kept
#'@return a DataFame of the infile
#'@export
loadTable <- function(infile, qualityCutoff=0){
  df <- utils::read.table(infile, sep="\t", header=TRUE, row.names=1)
  return(df[which(df$Quality >= 10),c("Relative.Time", "Length.of.Read", "Quality")])
}


#' Calculate read length distribution and coverage then save to .csv
#' 
#' Uses metadata DataFrame to calculate read distribution statistics and saves
#' to file. Not implemented yet: appending results to existing file where row 
#' name = sample name.
#' 
#' @param indf DataFrame of read metadata
#' @param outputname filename (.csv) for export, including path is not saved to
#' working directory
#' @param append TRUE or FALSE to append summary statistics to an existing file
#' @export
createTextSummary <- function(indf, outputname, append=FALSE){
  #basic read length statistics
  quantiles <- summary(indf$Length.of.Read)
  
  #coverage based on HGCh38.p13
  p13Tot=3272116950
  p13TotNonN=3110748599
  readSum <- sum(indf$Length.of.Read)
  coverageAll <- readSum/p13Tot
  coverageNonN <- readSum/p13TotNonN
  
  #The other definition of coverage can't be calculated pre alignment https://www.ecseq.com/support/ngs/how-to-calculate-the-coverage-for-a-sequencing-experiment
  #still need to look into how to use RSamtools on an index file, for this problem and also for filtering out only T-LRS reads.
  
df <- data.frame("coverage"=coverageAll, "coverageNonN"=coverageNonN, "minRead"=quantiles["Min."], "maxRead"=quantiles["Max."], "meanRead"=quantiles["Mean"], "N25"=quantiles["1st Qu."], "N50"=quantiles["Median"], "N75"=quantiles["3rd Qu."])
utils::write.csv(df, outputname)
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
#'@return a ggplot2 plot
#'@export
plotReadCountDistributions <- function(indf, outname, binNum=50, counts=FALSE, skinnyBin=FALSE, outwidth=8, outheight=7, outunits="in", plotQual=FALSE, qual.bins=c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)) {
  ggplot2::theme_set(ggplot2::theme_classic(base_size=12, base_family="Avenir"))
  if(counts){
    bw=10000
    if(skinnyBin){
      bw=5000
    }
    g <- plotCountBases(indf, binNum, plotQual=plotQual, qual.bins=qual.bins, binWidth=bw)
  } else {
    bw=5000
    if(skinnyBin){
      bw=2500
    }
    g <- plotCumulativeBases(indf, binNum, plotQual=plotQual, qual.bins=qual.bins, binWidth=bw)
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
#'@return a ggplot2 plot
#'@export
plotCountBases <- function(df, binNum, plotQual=TRUE, qual.bins=c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60), binWidth=8000, barColor="#8ec462"){
  bins <- seq(from=min(df$Length.of.Read), to=max(df$Length.of.Read), length.out=bin.num)
  df$bincode <- .bincode(df$Length.of.Read, bins, right=FALSE)
  if(plotQual==FALSE){
    df.countKB <- df %>% dplyr::group_by(bincode) %>% dplyr::summarize(Count=n())
    df.countKB <- data.frame(df.countKB)
    
    names(bins) <- seq(1, length(bins))
    df.countKB$bin.label <- bins[as.vector(df.countKB$bincode)]
    
    
    g <- ggplot2::ggplot(df.countKB, ggplot2::aes(x=bin.label, y=Count)) + ggplot2::geom_bar(stat="identity", width=binWidth, fill=barColor) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) + ggplot2::scale_x_continuous(n.breaks=20, breaks=waiver(),  label=scales::label_bytes("kB")) + ggplot2::xlab("Read Length") + ggplot2::ylab("Read Count") + ggplot2::scale_y_continuous(label= scales::label_bytes("MB"))
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
#'@return a ggplot2 plot
#'@export
plotCumulativeBases <- function(df, bin.num, plotQual=TRUE, qual.bins=c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60), binWidth=6000, barColor="#8ec462"){
  bins <- seq(from=min(df$Length.of.Read), to=max(df$Length.of.Read), length.out=bin.num)
  df$bincode <- .bincode(df$Length.of.Read, bins, right=FALSE)
  
  if(plotQual==FALSE){
    df.cumuKB <- df %>% dplyr::group_by(bincode) %>% dplyr::mutate("CumulativeKB" =  cumsum(Length.of.Read)) %>% dplyr::group_by(bincode) %>% dplyr::slice_max(CumulativeKB)
    df.cumuKB <- data.frame(df.cumuKB)
    
    names(bins) <- seq(1, length(bins))
    df.cumuKB$bin.label <- bins[as.vector(df.cumuKB$bincode)]
    
    
    g <- ggplot2::ggplot(df.cumuKB, ggplot2::aes(x=bin.label, y=CumulativeKB)) + ggplot2::geom_bar(stat="identity", width=binWidth, fill=barColor) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) + ggplot2::scale_x_continuous(n.breaks=20, breaks=waiver(),  label=scales::label_bytes("kB")) + ggplot2::xlab("Read Length") + ggplot2::ylab("Number of Bases") + ggplot2::scale_y_continuous(label= scales::label_bytes("MB"))
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
