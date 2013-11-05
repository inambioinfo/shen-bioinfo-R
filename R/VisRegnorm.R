VisRegnorm <- function(gene, spikein, x.col=1, y.col=2) {
# Visualize regression-based normalization.
# Args:
#   gene: Gene data table.
#   spikein: Spikein data table.
#   x.col: X sample position or name.
#   y.col: Y sample position or name.
    
    stopifnot(ncol(gene) >= 2 && ncol(spikein) >= 2)
    stopifnot(ncol(gene) == ncol(spikein))
    
    # Subset data.
    gene.sub <- gene[, c(x.col, y.col)]
    spikein.sub <- spikein[, c(x.col, y.col)]
    gene.sub.log <- log10(gene[, c(x.col, y.col)] + 1)
    spikein.sub.log <- log10(spikein[, c(x.col, y.col)] + 1)
    sample.names <- colnames(gene[, c(x.col, y.col)])
    sample.names <- paste("log10(", sample.names, " + 1)")

    # Function to plot gene and spikein data points.
    PlotDatpoint <- function(title) {
        smoothScatter(gene.sub.log, transformation=function(x) x^.18, 
                      xlab=sample.names[1],
                      ylab=sample.names[2],
                      nrpoints=0,
                      main=title)  # Gene.
        points(spikein.sub.log, pch=17)  # Spikein.
    }
    # Function to draw legend.
    DrawLegend <- function() {
        # Legend.
        leg.col <- c(blues9[7], "black", "red", "black")
        leg.pch <- c(20, 17, NA, NA)
        leg.lty <- c(0, 0, 1, 1)
        leg.lwd <- c(0, 0, 2, 2)
        legend("topleft", c("Gene", "Spikein", "Gene-fit", "Spikein-fit"),
               col=leg.col, pch=leg.pch, lty=leg.lty, lwd=leg.lwd)
    }
    
    oldpar <- par(mfrow=c(1, 2))
    
    # Plot Loess.
    # Gene x value min and max.
    PlotDatpoint("Loess")
    gene.minmax <- gene.sub.log[c(order(gene.sub.log[, 1])[1], 
                                  order(gene.sub.log[, 1], decreasing=T)[1]), ]
    spikein.aug <- rbind(spikein.sub.log, gene.minmax)
    loess.fit.gene <- loess.smooth(gene.sub.log[, 1], gene.sub.log[, 2])
    lines(loess.fit.gene, col="red", lwd=2)  # Gene-fit.
    loess.fit.spikein <- loess.smooth(spikein.aug[, 1], spikein.aug[, 2])
    lines(loess.fit.spikein, lwd=2)  # Spikein-fit.
    DrawLegend()
    
    # Plot sizeFactor.
    PlotDatpoint("sizeFactor")
    # Functiont calculate size factors.
    getSizeFactors <- function(dat) {
        library(DESeq)
        cds <- newCountDataSet(dat, c("x", "y"))
        cds <- estimateSizeFactors(cds)
        sizeFactors(cds)
    }
    gene.sf <- getSizeFactors(gene.sub)
    spikein.sf <- getSizeFactors(spikein.sub)
    # Function to draw lines based on size factor normalization.
    DrawSFLine <- function(xr, x.sf, y.sf, col, evaluation=50) {
        xp <- seq(log10(xr[1] + 1), log10(xr[2] + 1), length.out=evaluation)
        yp <- log10(y.sf / x.sf * (10^xp - 1) + 1)
        lines(xp, yp, col=col, lwd=2)
    }
    gsr <- range(gene.sub[, 1], spikein.sub[, 1])  # gene-spikein combined range.
    DrawSFLine(gsr, gene.sf[1], gene.sf[2], col="red")
    DrawSFLine(gsr, spikein.sf[1], spikein.sf[2], col="black")
    DrawLegend()
    

    par(oldpar)
}




















