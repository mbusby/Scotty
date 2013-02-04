function [ pwrBias ] = getPwrBias(powerByReadDepth,powerByReadDepthNoPoisson, thresh )
%{
    There is more power to call differential expression in genes 
    sequenced to greater depths than at lower depths because measurements at 
    low depths have relatively lower Poisson counting noise. When
    the bilogical and non-Poisson technical variance remains
    constant, there will be an increasing power to call the same
    effect size measurement differentially expressed as read
    count increases.  This increase will asymptote at a predictable 
    point: the power to call differential expression without poisson
    noise.  We refer to this asymptote value as "ultimate power".  While
    this power will never actually be reached, with reasonable rates of 
    biological variance at counts > ~1000 power is >99% of the ultimate 
    power.  At counts below 10, there is only roughly half of the total 
    to detect changes as is seen at the higher counts.

    In datasets that have low sequenceing or few reps the power
    to detect DE is low..

    In general, adding either deeper sequencing or more reps decreases
    the bias.  The relationships, however, aren't quite predictable:   
    3 reps sometimes has more bias than 2.  Generally, bias is highest 
    below ~10.  However, the "inflection point" of "10" can be lower
    if there are more replicates.  While you will always get less bias with
    deeper sequencing, sometimes it might be better to add a greater
    number of reps instead.  

    To demonstrate this effect, we defined a metric to quantify the bias in
    a dataset.  While there were many ways I could think of to do the, the
    metric we came up with that seemed to best display the effect we are 
    trying to get at was to look at the chunk of the genes that are really
    measured with low enough counts that there is a substantial bias against
    calling them differentially expressed.  So there are two parameters to 
    define here: 1) What is "substantial"? and the second is 2) how big is 
    the chunk.  We (somewhat arbitrarily) defined "substantically biased" 
    as having less than 80% of the ultimate power to detect differential 
    gene expression.  I chose 80% somewhat empirically: generally
    counts of 10 and 20 are used as points where data can be considered
    to be continuous, and the Poisson noise is therefore low. I looked at 
    charts and these had about 80% of ultimate power bias. 

    We left the 80% cut off as a parameter that could be altered.  Even if
    bias is removed to this point by deeper sequencing, an 80% difference
    in technial measurement bias is still a pretty big artifact, and it
    really will still need to be accounted for in downstream analyses.  If
    you are sequencing a species with a small transcriptome, such as yeast,
    it is worth considering whether you should just sequence deeply enough
    (say to 90% bias, min of ~100 reads) and then just ignore the 
    measurement bias all together.
 %}
    
    pwrBiasMatrix=powerByReadDepth./powerByReadDepthNoPoisson;
   
    passedThresh=pwrBiasMatrix>=thresh;
    pwrBias=sum(sum(passedThresh),2)/(size(pwrBiasMatrix,1)*size(pwrBiasMatrix,2));
    
    

end

