function [ pwr ] = getPowerTTest( m1, v1, m2, v2, n1, n2, al )
     %{
     Input the mean and variance of the data in ordinal space

     Find the power of the test using the formulas developed by 
     Chow, Shao, and Wang (J Biopharm Stat, 2002). 

     If one of the means is zero then the mean it is set to 0.001;
     Variance cannot be zero.
    
     mx1, mx2 are the means of the two distributions
     vx1, vx2 are the variances in log space
     n1, n2 are the sample sizes
     al is alpha, the false postive rate (this is the p value cut off)
    
     %}
  
     
     %Calc degrees of freedom
     %We tested but didn't use Scatterswait's because is performed poorly
     %when we tried it against simulation, particularly when there
     %was a difference in variance between the two samples.
     %df=  (v1./n1 + v2/n2).^2 ./ (((v1/n1).^2./(n1-1)) +  ((v2/n2).^2./(n2-1)) );  
     df = (n1+n2)-2;
     
     b=nctcdf(tinv((1-(al./2)), df), df, abs(m1-m2)./sqrt( v1./n1 + v2./n2))  - nctcdf(-tinv((1-(al./2)), df), df, abs(m1-m2)./sqrt( v1./n1 + v2./n2)) ;       
        
     pwr=1-b;
     
     %Test code to demonstrate it works:
    %{ 
    clear
    ctr=1;
    n1=5;
    n2=5;
    for i=1:1:50
        x(ctr)=i;
        al=0.01;
        mx1=i;%Mean in log space, sample 1
        fc=2;
        mx2=mx1*fc; %Mean in log space, sample 2

        od1=.5;
        od2=.5;

        %Total Variance
        vx1=mx1+((mx1)*od1)^2; 
        vx2=mx2+((mx2)*od2)^2;      

        %Variance without poisson noise 
        %(Just the overdispersion, poisson variance = the count)
        vx1NoPoisson=((mx1)*od1)^2;
        vx2NoPoisson=((mx2)*od2)^2;  

        pwrCalc(ctr) = getPowerTTest( mx1, vx1, mx2, vx2, n1, n2, al );

        %Simulate real test cases
        for j=1:1:100
            [mu1 sigma1]=getMuSigmaLognormal(mx1, vx1NoPoisson);
            [mu2 sigma2]=getMuSigmaLognormal(mx2, vx2NoPoisson);

            ls01=lognrnd(mu1, sigma1, n1,1);
            ls02=lognrnd(mu2, sigma2, n2,1);
            
            %ls01=normrnd(mx1, sqrt(vx1NoPoisson), n1,1);     
            %ls02=normrnd(mx2, sqrt(vx2NoPoisson), n2,1);    
     
            %Make random data drawn from a Poisson Lognormal
            ls1=poissrnd(ls01);
            ls2=poissrnd(ls02);

            %ns1=log(ls1);
            %ns2=log(ls2);
            [h p]=ttest2(ls1, ls2);
            ps(j)=p;
        end

        pwrMeas(ctr)=sum(ps<al)/length(ps);
        ctr=ctr+1;    
    end   


    figure
    plot(x,pwrCalc, 'r')
    hold on
    plot(x,pwrMeas, 'go')
    legend('Predicted', 'Measured in simulation')
    xlabel('Read Depth (Sample 1)')
    ylabel('Genes Detected as Differentially Expressed')

    %}
     

end

