#' Print Genome wide copynumber plot
#'
#' This function will take a dataframe with sample, chr, start.pos, end.pos and call.
#' It will output a full plot of all the segments
#' @title plotCopynumber
#' @param ReturnClass Class of type CNV. REQUIRED
#' @param segment_scale whether to plot CNV as discrete colors or continuous (heatmap)
#' @keywords Copynumber CNV
#' @export
#' @examples
#' Plot_Copynumber(CNVvault_Object)

plotCopynumber <- function(ReturnClass, genome="hg19", segment_scale = c("scna", "loh", "discrete"), rescale = TRUE, outliers = TRUE) {

        segment_scale = match.arg(segment_scale)
        print(segment_scale)

        if(segment_scale == "scna"){
                if(rescale){
                        max_scale = max(ReturnClass@Segments$calls, na.rm = TRUE)*0.5

                        if(outliers){
                                ReturnClass@Segments$calls[ReturnClass@Segments$calls > max_scale] <- max_scale
                        }


                        segment_scale <- scale_fill_gradientn(colors = c("blue", "white", "red"), values=scales::rescale(c(-2,0, max_scale)),
                                                              limits=c(-2,max_scale), breaks = -2:4, na.value = "black")
                } else {
                        segment_scale <- scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "black")
                }
        } else if(segment_scale == "loh"){
                segment_scale <- scale_fill_gradientn(colors = c("gray", "red"), values=scales::rescale(c(0.5,1)),
                                                      limits=c(0.5, 1), na.value = "black")
        } else if(segment_scale == "discrete"){
                segment_scale <- NULL
        }

        ReturnClass@Plot$plot <- ggplot(ReturnClass@Segments, aes(xmin=start, xmax=end, ymin=Ystart, ymax=Yend, fill=calls)) +
                segment_scale +
                # gghighlight::gghighlight(sampleID == "reference") +
                geom_rect() +
                ## Draw lines between the chromosomes
                geom_vline(xintercept = ReturnClass@Chr_Starts,
                           colour="grey",
                           alpha=0.5) +
                ## draw lines between the samples
                geom_hline(yintercept = 0:ReturnClass@NumberOfSamples,
                           colour="grey",
                           alpha=0.5) +
                ## Start x at 0. break halfway through each chromosome (this means the label is centered)
                 scale_x_continuous(expand = c(0, 0),
                                    breaks = ReturnClass@Chr_Starts[-length(ReturnClass@Chr_Starts)] +
                                            (ReturnClass@Chr_Sizes/2),
                                    ## set the labels to the chr names
                                    labels = gsub("chr","",
                                                  names(ReturnClass@Chr_Starts[-length(ReturnClass@Chr_Starts)])),
                                    position = "top") +
                ## Break y between each sample
                scale_y_continuous(expand=c(0,0),
                                   breaks = 0.5:(ReturnClass@NumberOfSamples-0.5), labels=levels(ReturnClass@Segments$sampleID)) +
                ## set the colours for the calls and remove the name
                # scale_fill_manual(name="", values=c("Gain"="#F8766D", "Loss"="#619CFF","CN-LOH"="#00BA38")) +
                ## Set axis labels
                labs(x="Chromosome", y="Sample") +
                theme_classic() +
                theme(axis.text.x = element_text(size=8,
                                                 face="bold"),
                      axis.text.y = element_text(size=10,
                                                 face="bold"),
                      axis.ticks = element_blank(),
                      legend.position = "right")
        print(ReturnClass@Plot$plot)
        ReturnClass
}
