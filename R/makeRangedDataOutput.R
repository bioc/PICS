makeRangedDataOutput<-function(obj, type="fixed", filter=list(delta=c(0,Inf),se=c(0,Inf),sigmaSqF=c(0,Inf),sigmaSqR=c(0,Inf),score=c(0,Inf)),length=100)
{
  nSe<-3
  if(type!="wig")
  {
    mu<-mu(obj)
    delta<-delta(obj)
    se<-se(obj)
    seF<-seF(obj)
    seR<-seR(obj)
    sf<-sigmaSqF(obj)
    sr<-sigmaSqR(obj)
    chromosome<-chromosome(obj)
    score<-score(obj)

    ## Filter regions with small deltas over a region
    if(!is.null(filter))
    {
      # indK<-unlist(sapply(obj@List,function(x,dmin){k<-K(x);if(k==0){return(NULL);};if(any(delta(x)<dmin & se(x)<20)){return(rep(FALSE,k));}else{return(rep(TRUE,k));}},dmin=filter$delta[1]))
      ### Filter based on delta
      ind1<-delta>filter$delta[1] & delta<filter$delta[2]
      ind2<-sf>filter$sigmaSqF[1] & sf<filter$sigmaSqF[2]
      ind3<-sr>filter$sigmaSqR[1] & sr<filter$sigmaSqR[2]
      ind5<-is.finite(score) & score>filter$score[1] & score<filter$score[2]
	    ind<-ind1&ind2&ind3&ind5
      
      if(!is.null(filter$se))
      {
      	ind4<-se>filter$se[1] & se<filter$se[2] & seF>filter$se[1] & seF<filter$se[2] & seR>filter$se[1] & seR<filter$se[2]
      	ind<-ind&ind4
      }
    }
    else
    {
      ind<-is.finite(score)
    }
  }
  if(type=="bed")
  {
    score<-(score(obj))[ind]
    ord<-order(-score)
    # Order the score    
    score<-score[ord]
    start<-(mu-delta/2-nSe*seF)[ind]
    end<-(mu+delta/2+nSe*seR)[ind]
    start<-start[ord]
    end<-end[ord]
    # chrom<-(paste("chr", chromosome, sep=""))[ord]    
    chrom<-(paste("",chromosome, sep=""))[ind]
    chrom<-chrom[ord]
    # strand<-NULL
  }
  else if(type=="ci")
  {
    score<-(score(obj))[ind]
    ord<-order(-score)
    score<-score[ord]
    start<-(mu-nSe*se)[ind]
    end<-(mu+nSe*se)[ind]
    start<-start[ord]
    end<-end[ord]
    chrom<-(paste("",chromosome, sep=""))[ind]
    chrom<-chrom[ord]
    # strand<-NULL
  }
  else if(type=="fixed")
  {
    score<-(score(obj))[ind]
    ord<-order(-score)
    score<-score[ord]
    start<-(mu-length)[ind]
    end<-(mu+length)[ind]
    start<-start[ord]
    end<-end[ord]
    chrom<-(paste("",chromosome, sep=""))[ind]
    chrom<-chrom[ord]
    strand<-NULL
  }
  else if(type=="wig")
  {
    temp<-wigDensity(obj,strand="*",step=10,sum=TRUE,filter=filter,scale=TRUE)
    chrom<-temp$chr
    start<-temp$x
    end<-temp$x+9
    score<-temp$density
    strand="*"
    # temp<-lapply(obj@List,"density","+",step=20,filter=filter)
    # startp<-unlist(lapply(temp,function(x){head(x$x,-1)}))
    # scorep<-unlist(lapply(temp,function(x){head(x$y,-1)}))
    # lx<-unlist(lapply(temp,function(x){if(is.null(x)){return(0);}else{return(length(x$x)-1);}}))
    # chr<-unlist(lapply(obj@List,function(x){if(is(x,"picsError")){return(NA);}else{return(x@chr);}}))
    # chrom<-(paste("", chr, sep=""))
    # chromp<-rep(chrom,lx)
    # endp<-startp+9
    
    # temp<-lapply(obj@List,"density","-",step=20,filter=filter)
    # startn<-unlist(lapply(temp,function(x){head(x$x,-1)}))
    # scoren<-unlist(lapply(temp,function(x){head(x$y,-1)}))
    # lx<-unlist(lapply(temp,function(x){if(is.null(x)){return(0);}else{return(length(x$x)-1);}}))
    # chr<-unlist(lapply(obj@List,function(x){if(is(x,"picsError")){return(NA);}else{return(x@chr);}}))
    # chrom<-(paste("", chr, sep=""))
    # chromn<-rep(chrom,lx)
    # endn<-startn+9
    
    # start<-c(startp,startn)
    # end<-c(endp,endn)
    # chrom<-c(chromp,chromn)
    # strand<-c(rep("+",length(startp)),rep("-",length(startn)))
    # score<-c(scorep,scoren)
  }
  
  if(length(score)==0) 
  {
  	gr=NULL
  }else
  {
	  ranges<-IRanges(as.integer(start),as.integer(end))
	  if(type=="bed" | type=="ci" | type=="fixed")
	  {
		names(ranges)<-paste("pics",1:(length(score)),sep="")
		gr<-RangedData(ranges, score, space = chrom)
#		gr<-GRanges(seqnames=chrom,ranges=ranges,strand="*",score)
	  }
	  else
	  {
		  
		  gr<-RangedData(ranges, score, space = chrom,  strand=strand)
		  #Create the new GRanges objec to return
#		  gr<-GRanges(seqnames=chrom, ranges=ranges, score, strand=strand)
	  }
  } 	
  return(gr)
}


##
# Make a RangedData object easier to understand that can be exported to different file formats.
#   type : fixed : fixed ranges around the predicted center
#          bed : ranges determined by the predicted delta and se
#          wig : density for windows of size 10
#   length : half the fixed length of the TF or nuclosome
##
MRO<-function(obj, type="fixed",filter=list(delta=c(0,Inf),se=c(0,Inf),sigmaSqF=c(0,Inf),sigmaSqR=c(0,Inf),score=c(0,Inf)),
		dataType="TF", length=NULL)
{
	if(is.null(length))
	{
		if(tolower(dataType)=="chip-seq" | tolower(dataType)=="tf")
		{
			length<-100
		}
		else if(tolower(dataType)=="mnase" | tolower(dataType)=="h")
		{
			length<-73
		}
	}

	mu<-mu(obj)
	delta<-delta(obj)
	score<-score(obj)
	se<-se(obj)
	seF<-seF(obj)
	seR<-seR(obj)
	sigmaSqF<-sigmaSqF(obj)
	sigmaSqR<-sigmaSqR(obj)
	chromosome<-chromosome(obj)
	
	nSe<-3

	## Filter regions with small deltas over a region
	if(!is.null(filter))
	{
		# indK<-unlist(sapply(obj@List,function(x,dmin){k<-K(x);if(k==0){return(NULL);};if(any(delta(x)<dmin & se(x)<20)){return(rep(FALSE,k));}else{return(rep(TRUE,k));}},dmin=filter$delta[1]))
		### Filter based on delta
		ind1<-delta>filter$delta[1] & delta<filter$delta[2]
		ind2<-sigmaSqF>filter$sigmaSqF[1] & sigmaSqR<filter$sigmaSqF[2]
		ind3<-sigmaSqR>filter$sigmaSqR[1] & sigmaSqR<filter$sigmaSqR[2]
		ind5<-is.finite(score) & score>filter$score[1] & score<filter$score[2]
		ind<-ind1&ind2&ind3&ind5
		
		if(!is.null(filter$se))
		{
			ind4<-se>filter$se[1] & se<filter$se[2] & seF>filter$se[1] & seF<filter$se[2] & seR>filter$se[1] & seR<filter$se[2]
			ind<-ind&ind4
		}
	}
	else
	{
		ind<-is.finite(score)
	}

	if(type=="fixed")
	{
		score<-score[ind]
		ord<-order(-score)
		score<-score[ord]
		start<-(mu-length)[ind]
		end<-(mu+length)[ind]
		start<-start[ord]
		end<-end[ord]
		chrom<-chromosome[ind]
		chrom<-chrom[ord]
	}
	else if(type=="bed")
	{
		score<-score[ind]
		ord<-order(-score)
		score<-score[ord]
		start<-(mu-delta/2-nSe*seF)[ind]
		end<-(mu+delta/2+nSe*seR)[ind]
		start<-start[ord]
		end<-end[ord]
		chrom<-chromosome[ind]
		chrom<-chrom[ord]
	}
	else if(type=="wig")
	{
		#if data.frame: cannot call density.
		temp<-density(obj,strand="*",step=10,sum=TRUE,filter=filter,scale=TRUE)
		chrom<-temp$chr
		start<-temp$x
		end<-temp$x+9
		score<-temp$density
		strand="*"
	}
	else
	{
		stop("Unknown type")
	}
		
	if(length(score)==0) 
	{
		gr=NULL
	}else
	{
		ranges<-IRanges(as.integer(start),as.integer(end))
		if(type=="bed" | type=="fixed")
		{
			names(ranges)<-paste("pics",1:(length(score)),sep="")
			gr<-RangedData(ranges, score, space = chrom)
		}
		else if(type=="wig")
		{
			#get the coverage
			cov<-coverage(ranges)$chrom
			#get the max coverage
			#if>1
			#for each base where coverage>1
			which(cov@values>1)
			#get all the ranges where they are, get the score, average it
			#
		}
		else
		{		
			gr<-RangedData(ranges, score, space = chrom,  strand=strand)
		}
	} 	
	return(gr)
}
