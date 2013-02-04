function [ m, v, p, c] = fitVariance( samples )
    %Fits the distribution of variances for a given dataset
    %to an underlying distribution
    tic 
    szSamples=size(samples);
    
    samplesN=getNormalizedSamples( samples );    
    
    mSample=mean(samplesN,2);
    vSample=transpose(var(transpose(samplesN)))-mSample;
    npStd=sqrt(vSample)./mSample;
    
    c=corr(npStd.^2, mSample);
    
    %Do quick and dirty to get approximate of v  
    
    %The measured mean should be the best approximation
    %of the true mean, only leaving the variance to be found
    m=mean(npStd(mSample>100 & npStd.^2>0));  
    
    %The maximum variance is expected be less than the observed variance
    %because the chi square stretches out the values when the variance 
    %is measured. In samples these seemed very small.
   
    vStart=0.001;
    vEnd=var(npStd(mSample>100))*1.25;
    vInc=(vEnd-vStart)/5;   
    
    vtest=vStart:vInc:vEnd;
    
    ctr=1;    
    ms=0;
    vs=0;
    ks=0;
    ps=0;
    trialRun=1;
    p=0;
    v=0;
    
    %since there is some flux in the Ks, we don't use a binary search
    while p<0.1 & trialRun<6

        ctr=1;    
        ms=0;
        vs=0;
        ks=0;
        ps=0;    
                
        parfor j=1:length(vtest)
        	[h p k]=testVarFit( npStd, mSample, size(samples,2), m, vtest(j), 3 ); %We bootstrap it 3 times to avoid random fluctuations in the simulation
            ms(j)=m;
           	vs(j)=vtest(j);
            ks(j)=k;
            ps(j)=p;                
        end
        
        m=mean(ms(ks==min(ks)));
        v=mean(vs(ks==min(ks)));
        p=mean(ps(ks==min(ks)));  
        
        startV=max(v-vInc, 0.001);
        endV=v+vInc;
        vInc=vInc/2;
        vtest=startV:vInc:endV;            
        trialRun=trialRun+1;
    end            
    

end

