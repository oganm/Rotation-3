# Specificity_index functions, version 1.0
# Created by Joe Dougherty, jdougherty_at_mail.rockefeller.edu or
# jodoc_a_ucla.edu
# Last updated 8/24/09
# The permutation testing is computationally intensive. Expect it to run for a while.

# Variable descriptions:
# summary_in - is a dataframe with expresion values for genes in rows, and samples
#or cell types in columns (at this point replicate arrays have been averaged, so
#one column per cell type)
#summary_fc is a matched array (same genes and samples) but with NA's for any
#genes that should be excluded for a particular cell type.
# We often filter prior to specificity index to remove genes below a particular
# background threshold.


# Define function
specificity_index<-function(summary_in){
    bts<-50  #set amount of distributions to average for permutation testing
    p_max<- .1 # set the maximum pvalue to be calculated
    e_min<- 50 #set the minimum expression value for a gene to be included
    #make a clean version of summary_in...
    dat_clean<-summary_in
    #...for each sample, remove values with less than
    # 50 absolute expression...
    dat_clean[summary_in<e_min]<-NA
    # for each sample remove those values that are designated to be filtered.
    #dat_clean[is.na(summary_fc)]<-NA
    #remove extra variable
    #rm(summary_fc)
    
    #make variables to catch the output of this program
    lng <- length(summary_in[1,])
    datComb<-array(NA,c(length(summary_in[,1]),2*lng)) #will be 2x the length of input - one column for rank, one for pvalue
    colnames(datComb)<-c(rep(colnames(summary_in), each=2))
    
    #for each sample(cell type
    for(j in 1:lng) {
        
        #.. determine which samples are not you
        notme =c(1:(j-1),(j+1):lng)
        
        if(j==1){
            notme=c(2:lng)}
        
        if(j==lng){
            notme=c(1:(lng-1))}
        
        #make a matrix of those genes that are rankable for this cell type.
        TmpMat<-summary_in[!is.na(dat_clean[,j]),]     # in other words, those that are NOT NA in the clean data for cell type j
        smps<-length(TmpMat[,1])   #  a variable for how many genes this is
        
        #for each other cell type sample with replacement from the values that are
        #not NA for this particular cell type.
        
        avg_ip_pre<-array(NA,c(smps,bts))  #blank a variable to hold the simulated distributions
        
        for(k in 1:bts){     #for 1 to the number of sampled distributions to create
            TmpMatSmp<-array(NA, c(smps, lng))    #blank a variable, number of (unflltered) genes by number of cell types
            for(i in 1:lng)  {  #for each cell type, sample with replacement from gene expression values
                TmpMatSmp[,i]<-sample(TmpMat[,i],smps, replace=TRUE)
            }
            #calcuate the log base 2 fold change for the current cell type (j) vs all others
            z<-log2(TmpMatSmp[,j]/TmpMatSmp[,notme])
            #blank variable
            rankz<-array(0,dim(z))
            z<-z*-1   #to invert the data, so low ranks end up being 'good', and NAs will go to the bottom
            #Calculate the rank of each gene within a each cell/cell comparison
            for(i in 1:(lng-1)) {rankz[,i]=rank(z[,i], na.last="keep")}
            #find the average rank for each gene across all samples
            avg_ip<-c()    #blank a variable
            # A distribution of the averages  each "gene"  rank
            for(i in 1:smps) {avg_ip[i]<-mean(rankz[i,])}
            avg_ip_pre[,k]<-avg_ip
        }      #end for k
        
        avg_ip_dist<-sort(avg_ip_pre)
        
        #now make the actual distribution of average rank IPs
        # blank a variable Z original (zO) to catch it
        
        zO<-array(NA,c(length(dat_clean[,1]), lng-1))
        
        #... then calculate fold change for this cell type versus all those others
        #(log base two of ratio )
        # Notice that because we are using cleaned data, we will only
        # calculate this when the numerator is not NA, but the demoninator
        # is not filtered.
        zO<-log2(dat_clean[,j]/summary_in[,notme])
        #blank an array to convert into from the data frame.
        rankzO<-array(0,dim(zO))
        zO<-zO*-1
        #Calculate the rank of each gene within a sample
        for(i in 1:(lng-1)) {rankzO[,i]=rank(zO[,i], na.last="keep")}
        # Move the values from the data frame to the array, so I can take a mean
        #blank a variable to catch the means
        avg_ipO<-c()
        #average the IPs for the actual values.
        for(i in 1:(length(rankzO[,1]))) {avg_ipO[i]<-mean(rankzO[i,], na.rm=FALSE)}
        
        # Make a picture of the situation.
        fnm = paste("hist_rnks",colnames(summary_in)[j],".png",sep="_")
        png(filename=fnm, width=960, height=960)
        par(mfrow=c(2,2))
        hist(avg_ip_dist, freq=FALSE , xlim=c(0,smps), main="Sampled distribution")
        hist(avg_ipO, freq=FALSE, xlim=c(0,smps), main="Actual distribution")
        dev.off()
        
        # make a variable to catch the Pvals
        pvals<-c()
        lng2<-length(avg_ipO)
        #for rankable genes, see where they would fall in the bootstraped distribution
        #use that to calculate a P value.
        # only calculate for p<p_max, to save time (usually p<.1)
        #so find cutoff p< p_max
        ctoff<-avg_ip_dist[(smps*bts)*p_max]
        
        for (i in 1:lng2)
        {
            if ((!is.na(avg_ipO[i]))&(avg_ipO[i]<ctoff))
            {
                avg_ip_dist[1]<-avg_ipO[i]
                pvals[i]<-((rank(avg_ip_dist, na.last="keep")[1])/(smps*bts))
            }else{
                pvals[i]<-NA
            }
        }
        datComb[,((j*2)-1)]<-avg_ipO
        datComb[,(j*2)]<-pvals
        
    } # toend J loop (each sample)
    datComb #return this
}#end function


# Future directions:
# 1) Allow users to set bts, p_max, e_min and make histogram output optional
# 2) Index the search for the pvalue, to make it faster.
# 3) Clean the code a bit more, and simplify
# 4) Get histograms to always match on bar width