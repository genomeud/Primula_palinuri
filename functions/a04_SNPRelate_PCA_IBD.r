# Copyright:	    Fabio Marroni and Paloma Perez 2024
# Aim:              Compute PCA and IBD using library SNPRelate
# To add:		     
# Suggestions: 
# Fixes:  


# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

# Required input files are missing fileterd vcf and population file. Last option "to_do" needs to be filled with (pca_analysis,ibd_analysis,FST) depending on the analysis to perform
option_list = list(
  make_option(c("-p", "--spp"), type="character",
  default="primula_palinuri",help="spp name [default= %default]", metavar="character"),
  make_option(c("-V", "--vcffile"), type="character",
  default="populations.snps.filtered.recode_MIS_filt.vcf",help="stacks vcf file [default= %default]", metavar="character"),
  make_option(c("-R", "--to_remove"), type="character",
  default=NULL,help="samples separated by comma to remove from the analysis [default= %default]", metavar="character"),
  make_option(c("-P", "--popfile"), type="character",
  default="raw_reads/sample_ID_and_populations.txt",help="Population file path [default= %default]", metavar="character"),
  make_option(c("-G", "--gdsfile"), type="character",
  default="SNPRelate/primula_palinuri.cov5.info50.gds",help="SNPRelate gds file [default= %default]", metavar="character"),
  make_option(c("-O", "--outpath"), type="character",
  default="SNPRelate/",help="Output directory [default= %default]", metavar="character"),
  make_option(c("-W", "--within"), type="character",
  default="SNPRelate/IBD_within.pdf",help="Within population IBD graph file [default= %default]", metavar="character"), 
  make_option(c("-M", "--maf"), type="numeric",
  default=0.05,help="Minimum MAF to include a SNP in analysis [default= %default]", metavar="character"), 
  make_option(c("-S", "--missingness"), type="numeric",
  default=0.25,help="Maximum missingness to include a SNP in analysis [default= %default]", metavar="character"), 
  make_option(c("-L", "--LD"), type="numeric",
  default=0.2,help="Threshold for LD pruning [default= %default]", metavar="character"), 
  make_option(c("-I", "--ibdmatfile"), type="character",
  default="SNPRelate/primula_palinuri.cov5.info50.ibd", help="IBD matrix output file [default= %default]", metavar="character"),
  make_option(c("-s","--step"), action="store",dest="to_do",default="pca_analysis", type='character', help="Step of the pipeline (pca_analysis,ibd_analysis) [%default]")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

  if (is.null(opt$spp)) {
  stop("WARNING: No spp name specified with '-p' flag.")
} else {  cat ("spp name ", opt$spp, "\n")
  spp <- opt$spp  
  }
  
if (is.null(opt$vcffile)) {
  stop("WARNING: No vcf file specified with '-V' flag.")
} else {  cat ("vcf file ", opt$vcffile, "\n")
  vcffile <- opt$vcffile  
  }

if (is.null(opt$to_remove)) {
  stop("WARNING: No remove samples file specified with '-R' flag.")
} else {  cat ("remove samples ", opt$to_remove, "\n")
  to_remove <- opt$to_remove  
  }
  
  if (is.null(opt$popfile)) {
  stop("WARNING: No popfile specified with '-P' flag.")
} else {  cat ("popfile is", opt$popfile, "\n")
  popfile <- opt$popfile  
  }
  
if (is.null(opt$gdsfile)) {
  stop("WARNING: No gds file specified with '-G' flag.")
} else {  cat ("gds file ", opt$gdsfile, "\n")
  gdsfile <- opt$gdsfile  
  }
  
if (is.null(opt$outpath)) {
  stop("WARNING: No outpath specified with '-O' flag.")
} else {  cat ("outpath ", opt$outpath, "\n")
  outpath <- opt$outpath  
  }
  
  if (is.null(opt$within)) {
  stop("WARNING: No within file specified with '-W' flag.")
} else {  cat ("within is", opt$within, "\n")
  within <- opt$within  
  }

  if (is.null(opt$maf)) {
  stop("WARNING: No maf specified with '-M' flag.")
} else {  cat ("maf is", opt$maf, "\n")
  maf <- opt$maf  
  }

  if (is.null(opt$missingness)) {
  stop("WARNING: No missingness specified with '-S' flag.")
} else {  cat ("missingness is", opt$missingness, "\n")
  missingness <- opt$missingness  
  }

  if (is.null(opt$LD)) {
  stop("WARNING: No LD specified with '-L' flag.")
} else {  cat ("LD is", opt$LD, "\n")
  LD <- opt$LD  
  }

  if (is.null(opt$ibdmatfile)) {
  stop("WARNING: No ibdmatfile specified with '-I' flag.")
} else {  cat ("ibdmatfile is", opt$ibdmatfile, "\n")
  ibdmatfile <- opt$ibdmatfile  
  }
  
  if (is.null(opt$to_do)) {
  stop("WARNING: No step specified with '--step' flag.")
} else {  cat ("Step to run is", opt$to_do, "\n")
  to_do <- opt$to_do  
  }
dir.create(dirname(ibdmatfile),showWarnings=F,recursive=T)

################
# PCA ANALYSIS #
################

run_pca<-function(spp,vcffile,popfile)
{
    library(gdsfmt)
    library(SNPRelate)
	threads=4	
	dir.create(outpath,recursive=T,showWarnings=F)
	outsuffix<-spp
    vcf.fn<-vcffile
    snpgdsVCF2GDS(vcf.fn, gdsfile, method="biallelic.only")
    # it is possible to read an already created gds object:
	(genofile <- snpgdsOpen(gdsfile))	
	sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    # exclude pops with few individuals 
    to_remove<-strsplit(to_remove, ",")[[1]]
    cat("Removed samples...\n")
    print(to_remove)
    sample.id<-sample.id[! sample.id %in% to_remove]
    # calculate PCA 
    pca<-snpgdsPCA(genofile,sample.id=sample.id, num.thread=threads,autosome.only=FALSE)
    # create output dataframe with the eigenvalues 
    tab <- data.frame(sample.id=pca$sample.id,EV1 = pca$eigenvect[,1], EV2=pca$eigenvect[,2],EV3 = pca$eigenvect[,3], EV4=pca$eigenvect[,4],stringsAsFactors=F)
    # read population information 
    pop_code <- read.table(popfile,stringsAsFactors=F,header=F,sep="\t")
    colnames(pop_code)<-c("sample.id","group_assignment")
    # merge the two datasets 
    total<-merge(tab,pop_code,by="sample.id",all.x=T,sort=F)
    # # store the dataframe with PCA coordinates 
    write.table(total,paste(outpath,"PCA.all_snp.",outsuffix,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)

    ########
    # PLOT #
    ########
    cex.point=0.9
    cex.axis=0.5
    cex.cross=0.7
    cex.sylvestrys=0.5
    cex.legend=0.7
    cex.text=1
    cex.main=1

    out_color=paste(outpath,"PCA.all_snp.",outsuffix,"_low_cov.cov5.info50.jpeg",sep="")
    jpeg(out_color,width=16,height=8,units="cm",res=300,type="cairo")
    par(mgp=c(1.25,0.3,0),oma=c(0.2,0.2,0.05,0.2),mar=c(2,3,1,1))

    #################
    # ASSIGN COLORS #
    #################
    col4="gray68"
    col5="orchid2"
    col6="seagreen"
    col7="tomato2"
    col8="hotpink4"
    col9="steelblue"
    col10="lightgoldenrod"
    
    total$color_group[total$group_assignment=="LAM"]<-col4
    total$color_group[total$group_assignment=="FIU"]<-col5
    total$color_group[total$group_assignment=="ID"]<-col6
    total$color_group[total$group_assignment=="PC"]<-col7
    total$color_group[total$group_assignment=="SGP"]<-col8
    total$color_group[total$group_assignment=="CIM"]<-col9
    total$color_group[total$group_assignment=="PPA"]<-col10
    
    color_list<-c(col4,col5,col6,col7,col8,col9,col10)
    
    test<-total
 	########
	# AXES #
	########
	# set axes -> round up to nearest --> pay attention to sign
	if (min(test$EV1)<0){
		x_le<-floor(min(test$EV1)*100)/100
	} else {
		x_le<-ceiling(min(test$EV1)*100)/100
	}
	if (max(test$EV1)<0){
		x_ri<-floor(max(test$EV1)*100)/100
	} else {
		x_ri<-ceiling(max(test$EV1)*100)/100
	}
	if (max(test$EV2)<0){
		y_up<-floor(max(test$EV2)*100)/100
	} else {
		y_up<-ceiling(max(test$EV2)*100)/100
	}
	if (min(test$EV2)<0){
		y_low<-floor(min(test$EV2)*100)/100
	} else {
		y_low<-ceiling(min(test$EV2)*100)/100
	}

	############
	# VARIANCE #
	############
	# calculate % variance of each eigenvector 
    # variance proportion (%)
    pc.percent <- pca$varprop*100
    head(round(pc.percent, 2))

	pc1=1
	pc2=2
	plot(test$EV1,test$EV2,col=test$color_group,cex=cex.cross,pch=c(0),xlim=c(x_le,x_ri),ylim=c(y_low,y_up),xlab="",ylab="",font.lab=2,las=1,cex.axis=cex.axis,xaxt="n",yaxt="n",tck=-0.01,main=c(expression(paste(italic("PCA")))),cex.main=cex.main)
    
	########
	# AXES #
	########
	#draw the x axis , since tick labels are to far from axis
	axis(1,at=axTicks(1),labels=round((axTicks(1)),2),cex.axis=cex.axis+.2,lwd=1,line=0,tck=-0.01,mgp=c(1.6,0.01,0),las=0,font.lab=2)
	# add the x axis of PC1 for Myles as text
	mtext(paste("PC",pc1," (",round(pc.percent[1],2), "%)",sep=""),1,font=2,cex=1,line=1.2)
	
	#draw the y axis , since tick labels are to far from axis
	axis(2,at=round(axTicks(2),2),,labels=round((axTicks(2)),2),cex.axis=cex.axis+.2,lwd=1,line=0,tck=-0.01,mgp=c(1.6,0.3,0),las=0,font.lab=2,las=1)
	# add the y axis of PC1 for Myles as text
	mtext(paste("PC",pc2," (",round(pc.percent[2],2), "%)",sep=""),2,font=2,cex=1,line=1.8)

	##########
	# LEGEND #
	##########
    par(family="Calibri")
    legend("bottomright",legend=c(expression(paste(italic("LAM"))),expression(paste(italic("FIU"))),expression(paste(italic("ID"))),expression(paste(italic("PC"))),expression(paste(italic("SGP"))),expression(paste(italic("CIM"))),expression(paste(italic("PPA")))),text.col=c(c(color_list)[1:7]),ncol=1,x.intersp=0.7,text.font=c(4),cex=cex.legend,pch=c(0),col=c(c(color_list)[1:7]),bty = "n",bg =NA)
    dev.off()

    snpgdsClose(genofile)

}

################
# IBD ANALYSIS #
################

run_ibd<-function(gdsfile,to_remove,popfile,within,ibdmatfile,maf,missingness,LD)
{
    library(data.table)
    library(RColorBrewer)
    library(gdsfmt)
    library(SNPRelate)
    library(dplyr)
    library(multcompView)

    #Load gds file (has to be previously saved)
    # genofile <- snpgdsOpen(gdsfile)
    genofile <- snpgdsOpen(gdsfile)
    #Prune for LD (otherwise IBD estimates might be inflated and running time longer)
    # snpgdsLDpruning Recursively removes SNPs within a sliding window
    # autosome.only: if ‘TRUE’, use autosomal SNPs only
    snpset <- snpgdsLDpruning(genofile, ld.threshold=LD,autosome.only=F)   
    snpset.id<-unlist(snpset)
   
    sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    
    # exclude pops with few individuals 
    to_remove<-strsplit(to_remove, ",")[[1]]
    print(to_remove)
    sample.id<-sample.id[! sample.id %in% to_remove]
    
    #Compute IBD and write IBD matrix just in case...
    # missing.rate: to use the SNPs with "<= missing.rate"
    # Maximum likelihood estimation (MLE) for the Identity-By-Descent (IBD) Analysis
    # Description:
    # Calculate the three IBD coefficients (k0, k1, k2) for non-inbred individual pairs by Maximum Likelihood Estimation.
    ibd <- snpgdsIBDMLE(genofile,sample.id=sample.id, maf=maf, missing.rate=missingness, num.thread=8,snp.id=snpset.id,autosome.only=F)
    # snpgdsIBDSelection Return a data frame with IBD coefficients
    ibd.coeff <- snpgdsIBDSelection(ibd)
    ibd.coeff<-ibd.coeff[order(ibd.coeff$kinship,decreasing=T),]

    write.table(ibd.coeff,ibdmatfile,row.names=F,quote=F,sep="\t")

    #Assign populations to individuals
    mypop<-fread(popfile,data.table=F)
    # remove samples to exclude from pop file!!!
    mypop<-mypop[! mypop$sample_ID %in% to_remove,]  

    # ibd coeff between same pops and all
    ibd.coeff$POP2<-ibd.coeff$POP1<-NULL
    ibd.coeff$POP1<-mypop$population[match(ibd.coeff$ID1,mypop$sample_ID)]
    ibd.coeff$POP2<-mypop$population[match(ibd.coeff$ID2,mypop$sample_ID)]

    #Only select within population comparisons
    intra<-ibd.coeff[ibd.coeff$POP1==ibd.coeff$POP2,]

    # create also vs all population
    intra2<-ibd.coeff
    intra2$POP1<-"all"

    intra<-rbind(intra,intra2)
    #Plot within pop relatedness

    intra$POP1<-factor(intra$POP1,levels=c("LAM","FIU","ID","PC","SGP","CIM","PPA","all"))
    color<-c("gray68","orchid2","seagreen","tomato2","hotpink4","steelblue","lightgoldenrod","white")
 
    # ---------------------- #
    # inbreeding coefficient #
    # ---------------------- #
    # # calculate the inbdreeding coefficient

    indinb.coeff.mom.weir <- snpgdsIndInb(genofile, sample.id=sample.id, maf=maf, missing.rate=missingness, snp.id=snpset.id,autosome.only=F,method="mom.weir")
    indinb.coeff.mom.visscher <- snpgdsIndInb(genofile, sample.id=sample.id, maf=maf, missing.rate=missingness,snp.id=snpset.id,autosome.only=F,method="mom.visscher")
    indinb.coeff.mle <- snpgdsIndInb(genofile, maf=maf, sample.id=sample.id, missing.rate=missingness, snp.id=snpset.id,autosome.only=F,method="mle")

    # assign now the value to each individual / population <- i don't understand the order of samples -
    # create a df with sample id + inbred + population
    out<-data.frame(sample.id=indinb.coeff.mle$sample.id,indinb.coeff.mle=indinb.coeff.mle$inbreeding,stringsAsFactors=F)
    out$indinb.coeff.mle<-as.numeric(out$indinb.coeff.mle)
    out$POP1<-mypop$population[match(out$sample.id,mypop$sample_ID)]
    #Plot within pop relatedness
    out$POP1<-factor(out$POP1,levels=c("LAM","FIU","ID","PC","SGP","CIM","PPA"))

    # create a new distribution which include all samples 
    out2=out
    out2$POP1="all"
    out<-rbind(out,out2)

    ################
    # PLOT drawing # 
    ################

    # draw a merged plot that integrate kinship + inbreeding
    res=1200
    cex.text=0.6
    cex_text=1.6
    cex.axis=0.75
    cex.significance=1
    cex.lab=1
    ylas=1
    pos=1

    outfile<-paste(dirname(within),"/kinship_within.inbreeding_coefficient.v2.jpeg",sep="")
    jpeg(outfile,width=16,height=12,units="cm",res=res, type="cairo")
    par(mar=c(0.8,3,1.2,0.1), mgp=c(1.6,0.5,0),tck=-0.03,oma=c(1.5,0.1,0.3,0.1),mfrow=c(2,1))

    ###################
    # PLOT A # whitin #
    ###################

    posi<-boxplot(intra$kinship~intra$POP1,col=color,cex.axis=0.8,ylab="",xlab="",font=2,yaxt="n",xaxt="n",ylim=c(-0.02,0.5))
    axis(2,line=0,las=ylas,cex.axis=cex.axis,mgp=c(0,0.5,0))
    mtext("Kinship coefficient",2,font=2,line=1.8,cex=cex.lab)
    # add the panel letter
    mtext(LETTERS[pos],2,line=2,at=0,padj=-8.5,adj=1,cex=cex_text,las=2,font=2)
    pos=pos+1

    # add the pairwise wilcoxon test 
    # create a dataframe with the values from the matrix     
    pp<-pairwise.wilcox.test(intra$kinship,intra$POP1,p.adjust.method="none", paired = FALSE, exact=FALSE)
    mymat<-tri_to_squ(pp$p.value)
    myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
    print(myletters)
    print(posi$stats)

    # # extract the order (median) of groups
    ordine<-posi$names[order(posi$stats[3,])] 
    lista<-as.list(strsplit(myletters$Letters,split=""))
    # create a new order based on the median average (obtained from boxplot data)
    intra$ordinato<-0
    for(tt in 1:length(ordine)){
        intra$ordinato[intra$POP1==ordine[tt]]<-tt
    }

    # now need to reassign the corresponding data to the names 
    pp1<-pairwise.wilcox.test(intra$kinship,intra$ordinato,p.adjust.method="none", paired = FALSE)
    # the colum and row have now the ordinato information 
    colnames(pp1$p.value)<-ordine[-length(ordine)]
    rownames(pp1$p.value)<-ordine[-1]
    # transfor matrix triangular en cuadrada
    mymat<-tri_to_squ(pp1$p.value)
    myletters[1]$Letters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)$Letters
    # order back as in boxplot in order to add the text to each box
    # myletters<-myletters[1]$Letters[order(names(myletters[1]$Letters))]
    myletters<-myletters[1]$Letters
    # myletters<-myletters[order(match(as.character(unique(comp1$V2)),names(myletters)))]
    myletters<-myletters[order(match(names(myletters),as.character(levels(intra$POP1))))]
    print(myletters)
    # for the y positions use 0
    # stat$stats[5,]+(letter_p_shift*(stat$stats[5,])/100)
    text(seq(1,length(posi$names),1),par('usr')[3]-((55*par('usr')[3])/100),myletters,cex=cex.significance,font=2,adj=c(0.5,0.5) )

    #############################
    # PLOT B # inbreeding (MLE) #
    #############################
    posi<-boxplot(out$indinb.coeff.mle~out$POP1,col=color,cex.axis=0.8,ylab="",xlab="",lab.font=2,yaxt="n",xaxt="n",ylim=c(0,1))
    axis(2,line=0,las=ylas,cex.axis=cex.axis,mgp=c(0,0.5,0))
    mtext("Inbreeding coefficient",2,font=2,line=1.8,cex=cex.lab)
    axis(1,at=c(1:6),tick=T,labels=FALSE,adj=0.5,padj=0.5,font=c(4),cex.axis=cex.axis)
    text(c(1:7),-0.18,c("LAM","FIU","ID","PC","SGP","CIM","PPA"),adj=0.5,padj=0.5,font=4,cex=cex.axis,xpd=NA,las=1)
    text(c(8),-0.18,c("all"),adj=0.5,padj=0.5,font=2,cex=cex.axis,xpd=NA,las=1)
    # add the panel letter
    mtext(LETTERS[pos],2,line=2,at=0,padj=-9.3,adj=1,cex=cex_text,las=2,font=2)

    # add the pairwise wilcoxon test 
    # create a dataframe with the values from the matrix 
    pp<-pairwise.wilcox.test(out$indinb.coeff.mle,out$POP1,p.adjust.method="none", paired = FALSE)
    mymat<-tri_to_squ(pp$p.value)
    myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
    print(myletters)
    print(posi$stats)

    # # extract the order (median) of groups
    ordine<-posi$names[order(posi$stats[3,])]
    lista<-as.list(strsplit(myletters$Letters,split=""))
    # create a new order based on the median average (obtained from boxplot data)
    out$ordinato<-0
    for(tt in 1:length(ordine)){
        out$ordinato[out$POP1==ordine[tt]]<-tt
    }
    # now need to reassign the corresponding data to the names 
    pp1<-pairwise.wilcox.test(out$indinb.coeff.mle,out$ordinato,p.adjust.method="none", paired = FALSE)
    # the colum and row have now the ordinato information 
    colnames(pp1$p.value)<-ordine[-length(ordine)]
    rownames(pp1$p.value)<-ordine[-1]
    mymat<-tri_to_squ(pp1$p.value)
    myletters[1]$Letters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)$Letters
    # order back as in boxplot in order to add the text to each box
    myletters<-myletters[1]$Letters
    myletters<-myletters[order(match(names(myletters),as.character(levels(out$POP1))))]
    print(par('usr')[3])
    text(seq(1,length(posi$names),1),par('usr')[3]-((125*par('usr')[3])/100),myletters,cex=cex.significance,font=2,adj=c(0.5,0.5) )
    dev.off()

}

tri_to_squ<-function(x)
{
    rn<-row.names(x)
    cn<-colnames(x)
    an<-unique(c(cn,rn))
    myval<-x[!is.na(x)]
    mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
    for(ext in 1:length(cn))
    {
     for(int in 1:length(rn))
     {
     if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
     mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
     mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
     }  
    }
    return(mymat)
}

if(to_do=="pca_analysis")
{
    run_pca(spp=spp,vcffile=vcffile,popfile=popfile)
} else if (to_do=="ibd_analysis")
{
    run_ibd(gdsfile=gdsfile,to_remove=to_remove,popfile=popfile,within=within,ibdmatfile=ibdmatfile,maf=maf,missingness=missingness,LD=LD)
}
