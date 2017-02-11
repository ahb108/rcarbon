
reScale <- function(x, type="simple", crng=NULL, na.rm=TRUE){

    types <- c("simple","normal", "custom")
    if (!type %in% types){
        stop("The rescale type you have chosen is not currently an option.")
    }
    if (na.rm){ x <- na.omit(x) }
    if (type=="normal"){
        res <- (x-mean(x))/sd(x)
    } else if (type=="custom"){
        if (is.null(crng)){
            stop("For custom type you need to specify a crng.")
        } else {
            xrange <- range(x)
            mfac <- (crng[2] - crng[1])/(xrange[2] - xrange[1])
            res <- crng[1] + (x - xrange[1]) * mfac
        }
    } else {
       res <- (x-min(x))/(max(x) - min(x))
   }
    return(res)
}

gaussW <- function(x, bw){
    exp(-(x^2)/(2*(bw^2)))
}

runMean <- function(x, n, edge="NA"){
    res <- x
    tmp <- filter(res,rep(1/n,n), sides=2)
    if (edge == "fill"){
        res[!is.na(tmp)] <- tmp[!is.na(tmp)]
    } else {
        res <- tmp
    }
    return(res)
}

quickMarks <- function(x, verbose=TRUE){

    if (!"calDates" %in% class(x)){
        stop("Input must be of class \"calDates\"")
    }
    df <- as.data.frame(matrix(ncol=8,nrow=length(x$grid)), stringsasFactors=TRUE)
    names(df) <- c("DateID","CRA","Error","qMed","q95s","q95e","q68s","q68e")
    print("Extracting approximate values...")
    if (length(x$grid)>1 & verbose){
        flush.console()
        pb <- txtProgressBar(min=1, max=length(x$grid), style=3)
    }
    for (a in 1:length(x$grid)){
        if (length(x$grid)>1 & verbose){ setTxtProgressBar(pb, a) }
        tmp <- x$grid[[a]]
        tmp <- tmp[tmp$PrDens>0,]
        tmp <- tmp[with(tmp, order(-PrDens)), ]
        tmp$Cumul <- cumsum(tmp$PrDens)
        tmp$Cumul <- tmp$Cumul/max(tmp$Cumul)
        tmp1 <- tmp[tmp$Cumul <= 0.95,]
        df[a,"q95s"] <- min(tmp1$calBP)
        df[a,"q95e"] <- max(tmp1$calBP)
        wdth <- max(tmp1$calBP)-min(tmp1$calBP)
        df[a,"q68s"] <- min(tmp1$calBP) + (wdth*0.25)
        df[a,"q68e"] <- max(tmp1$calBP) - (wdth*0.25)
        df[a,"qMed"] <- round(mean(c(df[a,"q95s"],df[a,"q95e"])),0)
        df[a,"DateID"] <- as.character(x$metadata$DateID[a])
        df[a,"CRA"] <- x$metadata$CRA[a]
        df[a,"Error"] <- x$metadata$Error[a]
    }
    if (length(x)>1 & verbose){ close(pb) }
    class(df) <- append(class(df),"quickMarks")
    return(df)
}

