plotEstimates <- function(list.complete,
                          return.estimates = "ALL",
                          attribute = NULL,
                          col = "auto",
                          plot.type = "p") {
    
    # select output
    if(return.estimates=="ALL") { return.estimates <- ncol(list.complete[[1]]) }
    
    # define colors
    if(col=="auto") { colorsmetric <- rainbow(length(list.complete)) }
    
    # prepare for plotting
    estimates.total <- length(return.estimates)-1
    estimates.count <- nrow(list.complete[[1]])
    
    # plot observed estimators
    if(estimates.count<=500) { lwd.by.iteration <- 1}
    if(estimates.count>500 && estimates.count<=1000) { lwd.by.iteration <- .7}
    if(estimates.count>1000 && estimates.count<=2000) { lwd.by.iteration <- .5}
    if(estimates.count>2000 && estimates.count<=4000) { lwd.by.iteration <- .3}
    if(estimates.count>4000) { lwd.by.iteration <- 0.1}
    
    # set plot window
    if(estimates.total<6) { plot.panels <- c(estimates.total,1) }
    if(estimates.total==6) { plot.panels <- c(3,2) }
    if(estimates.total==7) { plot.panels <- c(3,3) }
    if(estimates.total==8) { plot.panels <- c(2,4) }
    if(estimates.total==9) { plot.panels <- c(3,3) }
    if(estimates.total==10) { plot.panels <- c(5,2) }
    if(estimates.total>10 && estimates.total<17) { plot.panels <- c(4,4) }
    if(estimates.total==12) { plot.panels <- c(3,4) }
    if(estimates.total==15) { plot.panels <- c(3,5) }
    if(estimates.total>16 && estimates.total<18) { plot.panels <- c(4,5) }
    if(estimates.total==18) { plot.panels <- c(3,6) }
    if(estimates.total>18) { plot.panels <- c(5,5) }
    png(paste0("selected_network_estimates_complete_by_",tolower(attribute),".png"), type='cairo', width=plot.panels[2]*4,height=plot.panels[1]*4, units='in', res=200)
    par(mfrow=plot.panels, oma = c(4, 1, 1, 1))
    labels.plot1 <- 1:length(list.complete[[1]]$sample)
    labels.plot2 <- paste0(round(list.complete[[1]]$sample, 2)*100,"%")
    for(i in return.estimates[-1]) { 
        plot(as.numeric(list.complete[[1]][,i]), xlab="", ylab="", col=NULL, cex=0.5, xaxt="n", main=paste(colnames(list.complete[[1]])[i]), type=plot.type, lwd=lwd.by.iteration,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2)
        for(u in 1:length(list.complete)) {
            if(plot.type=="l") { lines(as.numeric(list.complete[[u]][,i]), xlab="", ylab="", col=colorsmetric[u], cex=0.5, xaxt="n", lwd=lwd.by.iteration,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2, type="l") }
            if(plot.type=="p") { points(as.numeric(list.complete[[u]][,i]), xlab="", ylab="", col=colorsmetric[u], cex=0.5, xaxt="n", lwd=lwd.by.iteration,cex.lab=1.6, cex.axis=1.6, cex.main=2.5, cex.sub=2, type="p") }
            
        }
        axis(1, at=labels.plot1, labels=labels.plot2)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", names(list.complete), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty="n", pch=19, col=colorsmetric, cex=3)
    dev.off() 
    return(print(paste0("Iteration completed. Plot saved at ",getwd(),"/selected_network_estimates_complete_by_",tolower(attribute),".png")))
}