function [Z_pos,Z_neg]=module_degree_zscore_sign(A,Ci)
%MODULE_DEGREE_ZSCORE       Within-module degree z-score
%
%   Z=module_degree_zscore(A,Ci);
%
%   The within-module degree z-score is a within-module version of degree
%   centrality.
%
%   Inputs:     A,      binary/weighted, directed/undirected connection matrix
%               Ci,     community affiliation vector
%
%   Output:     Z,      within-module degree z-score.
%
%   Note: The output for directed graphs is the "out-neighbor" z-score.
%
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%
%   Mika Rubinov, UNSW, 2008-2010


n=length(A);                        %number of vertices
Z_pos=zeros(n,1);
Z_neg=zeros(n,1);

A_pos=(A>0).*A;
A_neg=(-A>0).*-A;

for i=1:max(Ci)
    Koi_pos=sum(A_pos(Ci==i,Ci==i),2);
    Koi_neg=sum(A_neg(Ci==i,Ci==i),2);
    
    Z_pos(Ci==i)=(Koi_pos-mean(Koi_pos))./std(Koi_pos);
    Z_neg(Ci==i)=(Koi_neg-mean(Koi_neg))./std(Koi_neg);
end

Z_pos(isnan(Z_pos))=0;
Z_neg(isnan(Z_neg))=0;