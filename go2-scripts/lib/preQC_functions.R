inform=function(...){
    cat(paste0("[ 🛈 ] ", paste(...), "\n"))
    flush.console()
}

error=function(...){
    cat(paste0("[ ✗ ] ", paste(...), "\n"))
    stop(call.=F)
}


warn=function(...){
    cat(paste0("[ ⚠ ] ", paste(...), "\n"), call.=F)
    flush.console()
}

addCCinDF=function(sumStatsDF, sumStatsInfo){
    ncases=sumStatsInfo$NCASES
    ncontrols=sumStatsInfo$NCONTROLS
    inform("Added", ncases, "cases and", ncontrols, "controls in summary statistics.")
    sumStatsDF[,c("NCASES", "NCONTROLS", "N"):=list(ncases, ncontrols, ncases+ncontrols)]
    return(sumStatsDF)
}

readFnCheckerOutput=function(fn){
    d=fread(fn)
    expectedHeader=c("PATH","STUDY","PHENO",'ETHNICITY',"NCASES","NCONTROLS","SOFTWARE","PANEL","BUILD","DATE","ANALYST","GC","FF","ISERROR")
    if(any(colnames(d) != expectedHeader)){
        error("Error in file checker output: Expected header (", paste(expectedHeader, collapse=", "), "), received (", paste(colnames(d), collapse=", "), ")")
    }
    return(d)
}

discardFilenameError=function(outFilenameChecker, strict=T){
    errors=outFilenameChecker[ISERROR=="YES"]
    if(nrow(errors)){
        if(strict){
            error("Strict mode is on and there are", nrow(errors), "errors left in your filename checker output. Fix errors, rerun filename checks and rerun this program.")
        }else{
            warn(nrow(errors), "files/rows in your filename checker output have errors. Strict mode is off, they will be ignored.")
            return(outFilenameChecker[ISERROR=="NO"])
        }
    }else{
        inform("Read filename info from", nrow(outFilenameChecker), "files. No error detected.")
        return(outFilenameChecker)
    }
}


checkColnamesBeforeORConversion=function(data){
    mandatory=c("CHR", "POS", "EA", "NEA", "P", "INFO")
    curheader=colnames(data)
    if(! all(mandatory %in% curheader)) stop(paste("The following mandatory columns could not be found or matched:\n",
                                                           paste(setdiff(mandatory, curheader), sep=",")))
    if( !(all(c("BETA", "SE") %in% curheader) | all(c("OR", "OR_L95", "OR_U95") %in% curheader))) stop(paste("Neither BETA,SE or OR,OR_L95,OR_U95 was found or matched:",
                                                                                          paste(setdiff(mandatory, curheader), sep=",")))
    cat("[ ✔ ] All required columns found.\n")
    return(1)
}

fixChrCol=function(data){
    # checks if non-numeric chr codes are present
    # if yes, checks for trailing "chr" and removes it

    # if all are numeric, replaces 23 and 25 by "X"
    if(! "CHR" %in% colnames(data)) stop("CHR is not among the columns")
    if(class(data$CHR) %in% c("numeric", "integer")){
        data[,CHRNUM:=CHR] #it is easier and faster to do matching on integers
        inform("All chromosome codes are numeric.")
        if(nrow(data[CHR>22])){
            nonAutoChrCodes=unique(data[CHR>22]$CHR)
            inform("Non-autosomal chromosome codes detected:", paste(nonAutoChrCodes, collapse=","))
            if(!all(nonAutoChrCodes %in% c(23,25))) error("Only 23 and 25 are supported as chrX numerical codes.")
            data[CHRNUM==25,CHRNUM:=23]
            data[,CHR:=as.character(CHR)]
            numconv=nrow(data[CHR %in% c("23", "25")])
            if(numconv){
                data[CHR %in% c("23", "25"),CHR:="X"]
                inform(numconv, "rows on X chromosome converted from numerical code.")
            }
        }else{
        inform("All positions are autosomal (no X detected).")
        }
        data[,CHR:=as.character(CHR)]
    }else{
        # column was not cast to numeric: there are strings
        # Step 1: check for "chr"
        # if there is a chr, it should be everywhere unless the file is very broken
        # we sample 1000 random lines and look for chr
        subset=tolower(data[sample(nrow(data), 1000)]$CHR)
        if(any(grepl("chr", subset))){
            inform("Removing \"chr\" prefix in chromosome names.")
            data[,CHR:=sub("^chr", "", CHR)]
        }

        nonIntChrCodes=unique(data$CHR)
        nonIntChrCodes=nonIntChrCodes[! nonIntChrCodes %in% as.character(c(1:22, 23, 25))]
        if(length(nonIntChrCodes)){
            inform("Non-integer chromosome codes found:", paste(nonIntChrCodes, collapse=","))
            toRemove=setdiff(nonIntChrCodes, "X")
            nRemove=nrow(data[CHR %in% toRemove])
            if(nRemove){
                data=data[! CHR %in% toRemove]
                warn( nRemove, "rows removed due to unparsable non-numeric code.")
            }
        }
        numconv=nrow(data[CHR %in% c("23", "25")])
        if(numconv){
            data[CHR %in% c("23", "25"),CHR:="X"]
            inform( numconv, "rows on X chromosome converted from numerical code.")
        }

        rem=data[!(CHR %in% c(as.numeric(1:22), "X"))]
        if(nrow(rem)){
            print(head(rem))
            error("This should not be possible, but there are irregular chromosomes remaining (see above).")
        }
        data[,CHRNUM:=CHR]
        data[CHRNUM=="X", CHRNUM:=23]
        data[,CHRNUM:=as.integer(CHRNUM)]
    }
    #  - the CHR col is character
    #  - the CHRNUM col contains a numeric code for the chr
    #  - there are no exotic chromosomes left
    #  - the "chr" prefix has been removed if present
    # The error behaviour is currently asymmetrical:
    #  - if there were non-autosomal-or-X INT codes for chr (24, 32, ...) there was an error
    #  - if there were non-autosomal-or-X string codes (chrY, mt, ...) they were removed
    return(data)
}

# checkChrCol
# checks if all chromosomes are present
# checks range of chromosomes
# plots number of variants per chr
# plots histogram of chr vars

checkChrCol=function(data){
    # Check if all autosomes are present
    presentChr=unique(data$CHRNUM)
    missingChr=setdiff(1:22, presentChr)
    if(length(missingChr)){
        stop(paste("The following autosomes are missing:", paste(missingChr, collapse=","), "\n"))
    }
    else{
        inform("All autosomes are present.\n")
    }

    # check range of chromosomes
}

checkChrColXchrFiles=function(data){
    # Check if only Chr
    presentChr=unique(data$CHRNUM)
    autoChr=setdiff(23, presentChr)
    
    if(presentChr!=23){
        stop(paste("The following autosomes are present:", paste(autoChr,  collapse=","), "\n"))
    }else{
        # Only X chromosome present
        inform("Only X chromosome present.\n")
    }
    # check range of chromosomes
}


checkHeader=function(sumStatsPath, colNameMappingTable){
    d=fread(sumStatsPath)
    inform("Read file", sumStatsPath, ".")
    curheader=colnames(d)
    patterns=strsplit(colNameMappingTable$V2, ",", fixed=T)
    existing=c()
    replacement=c()
    for(i in 1:nrow(colNameMappingTable)){
        pat=patterns[[i]]
        pat=c(colNameMappingTable$V1[i], pat)
        matched=F
        for( p in pat ){
            if(grepl("*", p, fixed=T)){
                #cat(paste("Found a pattern to search: ", p, "\n"))
                matches=grepl(p, curheader, ignore.case=T)
            }else{
                #cat(paste("Found an exact match to search: ", p, "\n"))
                matches=toupper(curheader)==toupper(p)
            }
            if(sum(matches)>1) stop(paste("Expression", p, "matched several columns : ", paste(curheader[matches], collapse=" ")))
            if(sum(matches)==0) next
            matched=T
            if(curheader[matches] == colNameMappingTable$V1[i]) next
            existing=c(existing, curheader[matches])
            replacement=c(replacement, colNameMappingTable$V1[i])
        }
        if(!matched) warn(paste("Could not find any of", colNameMappingTable$V2[i], "in header:\n", paste(curheader, collapse=", ")))
    }

    if(length(existing)){
        cat(paste("Found the following entries to replace/rename:\n"))
        replacement=data.table(ex=existing, post=replacement)
        replacement=unique(replacement)
        print(cbind(replacement$ex, replacement$post))
        setnames(d, replacement$ex, replacement$post)
    }
    return(d)
}

checkAndConvertORtoBeta=function(sumStatsDF, sumStatsInfo, softList, strict=F){
    isBeta=all(c("BETA", "SE") %in% colnames(sumStatsDF))
    isOR=all(c("OR", "OR_L95", "OR_U95") %in% colnames(sumStatsDF))

    ## get column names for type conversion
    coltocheck=NULL
    toconvert=NULL
    
    if(isOR){
      ORcol=c("P", "OR", "OR_L95", "OR_U95")
      types=sumStatsDF[,lapply(.SD, class),.SDcols=ORcol]
      toconvert=names(types)[types!="numeric"]
      coltocheck=ORcol
    }
    
    if(isBeta){
      ORcol=c("P", "BETA", "SE")
      types=sumStatsDF[,lapply(.SD, class),.SDcols=ORcol]
      toconvert=names(types)[types!="numeric"]
      coltocheck=c(coltocheck, ORcol)
    }
    
    # required if both OR and Beta are there
    toconvert=unique(toconvert)
    coltocheck=unique(coltocheck)
    
    # Type conversion block
    if(length(toconvert)){
      warn("Columns", paste(toconvert, collapse=", "), "are not numeric and will be converted.")
      sumStatsDF[,c(toconvert) := lapply(.SD, as.numeric), .SDcols=toconvert]
      prevN=nrow(sumStatsDF)
      torm=sumStatsDF[,lapply(.SD, function(x) {!is.finite(x)}), .SDcols=toconvert]
      torm=rowSums(torm)
      sumStatsDF=sumStatsDF[torm<1,]
      warn(prevN-nrow(sumStatsDF), "rows were lost during the conversion.")
    }else{
      inform("Checking columns", paste(coltocheck, collapse=", "), "for non-finite values.")
      prevN=nrow(sumStatsDF)
      torm=sumStatsDF[,lapply(.SD, function(x) {!is.finite(x)}), .SDcols=coltocheck]
      torm=rowSums(torm)
      sumStatsDF=sumStatsDF[torm<1,]
      warn(prevN-nrow(sumStatsDF), "rows were lost during the conversion.")
      
    }
    
    if(isOR & !isBeta){
        inform("Only OR are present. Converting OR to Betas.")
        sumStatsDF[,c("BETA", "SE") := list(log(OR), ((log(OR_U95)-(log(OR)))/1.95996))]
        return(sumStatsDF)
    }
    if(isOR & isBeta){
        inform("Dataframe contains both betas and OR, no conversion needed: using betas.")
        return(sumStatsDF)
    }
    #implied: if !isOR (we only have betas)
    if(!is.null(softList)){
      if(toupper(sumStatsInfo$SOFTWARE) %in% toupper(softList)){
        if(sumStatsInfo$FF=="FFNo"){
          inform("Software", toupper(sumStatsInfo$SOFTWARE), "is listed as a QT model and fudging will be applied.")
          sumStatsDF[,fudgefactor := 1/((NCASES/(NCONTROLS+NCASES))*(1-(NCASES/(NCASES+NCONTROLS))))]
          sumStatsDF[,c("BETA", "SE") := list(BETA*fudgefactor, SE*fudgefactor)]
        }else{
          inform("Software", toupper(sumStatsInfo$SOFTWARE), "is listed as a QT model but fudging has already been applied.")
        }
      }else{
        # the list of software is there but ours is not among them and only betas are present
        if(sumStatsInfo$FF=="FFNo") {
          if(strict) error("Only betas are present, but the software", sumStatsInfo$SOFTWARE, "is not in the linear list and fudging has not been applied. Exiting because of strict mode.")
          inform("Only betas are present, but the software", sumStatsInfo$SOFTWARE, "is not in the linear list. Fudging will NOT be applied, double-check this is what you want.")
        }else{
          inform("Only betas are present, the software", sumStatsInfo$SOFTWARE, "is not in the linear list, but fudging has already been applied.")
        }
      }
      return(sumStatsDF)
    }

    
    if(sumStatsInfo$FF=="FFYes"){
        inform("Only betas are present and fudge factor has been applied. Nothing to be done.")
        return(sumStatsDF)
    }else{
        if(strict) error("Only betas are present, but no linear software list has been supplied and fudging has not been applied. Exiting because of strict mode.")
        inform("Only betas are present and fudge factor has not been applied. Applying FF to betas.")
        sumStatsDF[,fudgefactor := 1/((NCASES/(NCONTROLS+NCASES))*(1-(NCASES/(NCASES+NCONTROLS))))]
        sumStatsDF[,c("BETA", "SE") := list(BETA*fudgefactor, SE*fudgefactor)]
        return(sumStatsDF)
    }
}

checkOrFillEAF=function(sumStatsDF, sumStatsInfo){
    if("EAF" %in% colnames(sumStatsDF)){
        inform("EAF present in summary stats. Not recalculating.")
        return(sumStatsDF)
    }
    if(sumStatsInfo$SOFTWARE == "SNPTEST"){
        inform("EAF not present in summary stats and SNPTEST format detected. Calculating EAF.")
        colnames(sumStatsDF)=toupper(colnames(sumStatsDF))
        sumStatsDF[,EAF := (2*(ALL_BB)+ALL_AB)/(2*(ALL_AB + ALL_BB + ALL_AA))]
        return(sumStatsDF)
    }
    error("EAF column not present in summary statistics.")
}

checkAbnormalBetas=function(sumStatsDF){
  numberAbnormal=nrow(sumStatsDF[EAF>0.05 & EAF<0.95 & abs(BETA)>0.5])
  if(numberAbnormal){
    warn(numberAbnormal, "common variants out of ",nrow(sumStatsDF),"(",round(100*numberAbnormal/nrow(sumStatsDF)),"%) have effects above 0.5.")
    print(head(sumStatsDF[EAF>0.05 & EAF<0.95 & abs(BETA)>0.5]))
  }
}

createAlphabeticCPTID=function(sumStatsDF){
    # this is an expensive function
    inform("Generating alphabetic IDs... this will take time.")
    tomerge=dcast(melt(sumStatsDF[,InternalID:=.I][,c("InternalID", "EA", "NEA")], id.var="InternalID")[order(value), .SD, by="InternalID"][,N:=paste0("Col", 1:.N) , .(InternalID)], InternalID~N, value.var=c("variable", "value"))[,alph:=paste(value_Col1, value_Col2, sep="_")][,c("InternalID", "alph")]
    sumStatsDF=merge(sumStatsDF, tomerge, by="InternalID")[,CPTID:=paste0(CHR, ":", POS, "_", alph)]
    return(sumStatsDF)
}

createDiagnosticsPlots=function(sumStatsDF, reportFilename, sumStatsBuild, chromLengths){
    #     if(names(dev.cur())!="pdf"){
    #         # We are not printing to a PDF right now, open one:
        pdf(reportFilename, width=25, height=10)
    #     }
    chrRanges=sumStatsDF[,as.list(c(range(POS), .N)), by="CHRNUM"]
    setnames(chrRanges, c("chr", "start", "end", "N"))
    chrRanges[,size_mb:=round((end-start)/1000000)]
    setkey(chrRanges, "chr")

    chroms=unique(sumStatsDF$CHRNUM)
    check_chr23= 23 %in% chroms

    # Chromosome sizes and num of variants
        if(check_chr23){

        #plt1 Chromosome sizes
        maxall=max(c(chromLengths, chrRanges$size_mb*1000000))
        plot(chrRanges$chr, chrRanges$size_mb, xlab="chromosome", ylab="Size (Mbp)", xaxt="n")
        points(chrRanges$chr, chromLengths[1:23]/1000000, pch="_", col="tomato2", cex=1.5)
        axis(1, 1:23)
        abline(v=1:23, lty=2)

        #plt2 No.of variants per Chromosome
        plot(chrRanges$chr, chrRanges$N, xlab="chromosome", ylab="Number of variants", xaxt="n")
        axis(1, 1:23)
        abline(v=1:23, lty=2)

    }else{
        #plt1 Chromosome sizes
        plot(chrRanges$chr, chrRanges$size_mb, xlab="chromosome", ylab="Size (Mbp)", xaxt="n")
        points(chrRanges$chr, chromLengths[1:22]/1000000, pch="_", col="tomato2", cex=1.5)
        axis(1, 1:22)
        abline(v=1:22, lty=2)

        #plt2 No.of variants per Chromosome
        plot(chrRanges$chr, chrRanges$N, xlab="chromosome", ylab="Number of variants", xaxt="n")
        axis(1, 1:22)
        abline(v=1:22, lty=2)
    }

    # coverage plot
    #plt3
    if(check_chr23){
        allstart=min(chrRanges$start)-1
        allend=max(chrRanges$end)+1
        allrange=seq(allstart, allend, by=5000000)
        allrange=c(allrange, allend)
        sumStatsDF[,chrChunk:=as.numeric(cut(POS, breaks=allrange))]
        varnum=sumStatsDF[,.N, by=.(CHRNUM, chrChunk)]
        varnum[,prop:=N/max(N), by=.(CHRNUM)]
        options(repr.plot.width=15, repr.plot.height=10)

        plot(0, xlim=c(0, max(varnum$chrChunk)+1), ylim=c(1,24), type="n", bty="n", xaxt="n", xlab="position", yaxt="n", ylab="chromosome")
        abline(v=seq(0, max(varnum$chrChunk)+1), lty=2, col="lightgray")
        rect(ybottom=varnum$CHRNUM, ytop=varnum$CHRNUM+varnum$prop*0.8, xleft=varnum$chrChunk-1, xright=varnum$chrChunk, col="gray")
        text(x=-0.7, y=1:23+0.3, 1:23)
        axis(1, at=seq(0, 50, by=2), paste0(seq(0, 25), "M"))
    }else{
        allstart=min(chrRanges$start)-1
        allend=max(chrRanges$end)+1
        allrange=seq(allstart, allend, by=5000000)
        allrange=c(allrange, allend)
        sumStatsDF[,chrChunk:=as.numeric(cut(POS, breaks=allrange))]
        varnum=sumStatsDF[,.N, by=.(CHRNUM, chrChunk)]
        varnum[,prop:=N/max(N), by=.(CHRNUM)]
        options(repr.plot.width=15, repr.plot.height=10)

        plot(0, xlim=c(0, max(varnum$chrChunk)+1), ylim=c(1,23), type="n", bty="n", xaxt="n", xlab="position", yaxt="n", ylab="chromosome")
        abline(v=seq(0, max(varnum$chrChunk)+1), lty=2, col="lightgray")
        rect(ybottom=varnum$CHRNUM, ytop=varnum$CHRNUM+varnum$prop*0.8, xleft=varnum$chrChunk-1, xright=varnum$chrChunk, col="gray")
        text(x=-0.7, y=1:22+0.3, 1:22)
        axis(1, at=seq(0, 50, by=2), paste0(seq(0, 25), "M"))
    }


    # qq plot (p-values)
    ##plt4 -Autosomes only
    library(manqq)
    if(check_chr23){
        newStats= sumStatsDF[sumStatsDF$CHRNUM != 23]

        # Autosomes only
        lambda=fastqq(newStats$P)
        mtext("Uniform QQ-Plot, p-values (Autosomes only)")
        if(lambda > 1.2 | lambda < 0.8) {warn("lambda = ", lambda, "deviates from 1.")}
        # X-chrom only
        sumStatsXchr= sumStatsDF[sumStatsDF$CHRNUM == 23 ]
        lambda=fastqq(sumStatsXchr$P)
        mtext("Uniform QQ-Plot, p-values (X-Chromosome only)")

    }else{
        lambda=fastqq(sumStatsDF$P)
        mtext("Uniform QQ-Plot, p-values (Autosomes only)")
        if(lambda > 1.2 | lambda < 0.8) {warn("lambda = ", lambda, "deviates from 1.")}
    }

    # qq plot (beta-se/NORM)
    ##plt5. QQ-plot sample quantiles
    if(check_chr23){
        ## Autosomes only
        newStats= sumStatsDF[sumStatsDF$CHRNUM != 23]
        smol=newStats[sample(.N, 1000000)]
        qqnorm(smol[ ,BETA/SE], pch=".", cex=1.6, col="gray5", main="Normal QQ-Plot, beta/se (what METAL will see)\n 1M random variants - Autosomes only")
        abline(a=0, b=1, col="tomato2", lty=2)

        # beta vs EAF
        plot(smol$EAF, smol$BETA, pch=".", cex=1.6, col="gray5", xlab="EAF", ylab="Beta", main="beta vs EAF, 1M random variants - Autosomes only")

        ## X-Chromosome only
        sumStatsXchr= sumStatsDF[sumStatsDF$CHRNUM == 23 ]
        smol=sumStatsXchr[sample(.N, 50000)]
        qqnorm(smol[ ,BETA/SE], pch=".", cex=1.6, col="gray5", main="Normal QQ-Plot, beta/se (what METAL will see)\n 50k random variants - X-Chromosome only")
        abline(a=0, b=1, col="tomato2", lty=2)

        # beta vs EAF
        plot(smol$EAF, smol$BETA, pch=".", cex=1.6, col="gray5", xlab="EAF", ylab="Beta", main="beta vs EAF, 50k random variants - X-Chromosome only")

    }else{
        smol=sumStatsDF[sample(.N, 1000000)]
        qqnorm(smol[ ,BETA/SE], pch=".", cex=1.6, col="gray5", main="Normal QQ-Plot, beta/se (what METAL will see)\n 1M random variants Autosomes only")
        abline(a=0, b=1, col="tomato2", lty=2)

        ##plt6. beta vs EAF
        # beta vs EAF
        plot(smol$EAF, smol$BETA, pch=".", cex=1.6, col="gray5", xlab="EAF", ylab="Beta", main="beta vs EAF, 1M random variants - Autosomes only")
    }

    # manhattan
    smol=sumStatsDF[,.(CHRNUM, POS, P)]
    setnames(smol, c("chr", "pos", "p"))
    print(head(smol))
    forgetme=fastmanh(smol, build=sumStatsBuild, maxpeaks=10, no_distance=T, no_annot=T)

    # pvalues
    ##plt7. hist betas for all variants
    if (check_chr23){
        newStats= sumStatsDF[sumStatsDF$CHRNUM != 23]
        hist(newStats$BETA, breaks=101, main="Histogram of betas (all variants) - Autosomes only", xlab="")

        #X-chrom only
        sumStatsXchr= sumStatsDF[sumStatsDF$CHRNUM == 23]
        hist(sumStatsXchr$BETA, breaks=101, main="Histogram of betas (all variants) - X-Chromosome only", xlab="")

    }else{
        hist(sumStatsDF$BETA, breaks=101, main="Histogram of betas (all variants) - Autosomes only", xlab="")
    }

    # allele tables and info
    sumStatsDF[,c("lA", "lB"):=lapply(.SD, str_count), .SDcols=c("EA", "NEA")]
    maxAL=max(c(sumStatsDF$lA, sumStatsDF$lB))
    sentence=paste("\nMaximum allele length is", maxAL)
    report=sentence
    inform(sentence)
    if(maxAL==1){
      sentence="No full-length indels are included in allele codes.\n"
      inform(sentence)
      report=c(report, sentence)
    }else{
        sublarge=sumStatsDF[lA==maxAL | lB==maxAL]
        sentence=paste(nrow(sublarge), "positions have the maximum allele length.")
        inform(sentence)
        report=c(report, sentence)
        sentence=paste("Total number of indels defined as l(A)>1 is", nrow(sumStatsDF[lA>1 | lB>1]), ",", round(nrow(sumStatsDF[lA>1 | lB>1])/nrow(sumStatsDF)*100), "% of total.")
        report=c(report, sentence)
        inform(sentence)
        sentence="Example positions with largest allele length:"
        report=c(report, sentence)
        report= c(report,strwrap(captureOutput(print(head(sublarge))),width = 250))
    }
    #print(colnames(sumStatsDF))
    biallelic=sumStatsDF[lA==1 & lB==1]
    report=c(report, "TABLE OF ALLELES IN BIALLELIC SNPs\n==================================\n",
    captureOutput(table(biallelic[,c("EA", "NEA")])))
    textplot(report, fixed.width=TRUE, mar=c(0, 0, 3, 0) + 0.1, halign='left', valign="top", cex=0.85,
        cspace=1, lspace=1)
    dev.off()

}
