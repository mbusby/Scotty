function [ clusteringChart ] = clusterSamples( data, sampleNames )
%Clusters data for quality control purposes
    
  
    %To avoid results dominated by Poisson noise,
    %only cluster samples with greater than 10 reads per gene
    %on average    
    data=data(median(data,2)>10 ,:);
    X=transpose(data);
   
    Y = pdist(X, 'spearman' );
    Z = linkage(Y,'average');
    T = cluster(Z,5);
          
    clusteringChart=figure('Visible','on');
    [H,T,perm] = dendrogram(Z, 'labels', sampleNames, 'orientation','left' );
    xlabel('Distance (1 - Spearman Rank Correlation)', 'FontSize', 14);
    topX=max(Y)/2*1.05;
    mY=max(Y)
    
    %axis([0 topX 0 size(data,2)+1]);

end

