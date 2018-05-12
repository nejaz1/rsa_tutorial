function [c,cv,d,dv]=crossval_correlation(Y,varargin)
% Computes the correlation between individuals and the group mean 
%     based on leave-one-out-crossvalidation and in a non-crossvalidated 
%     fashion 
% Based on
% INPUT:
%   D   : NxK distances for N subjects and K measurements
% VARARGIN:
%   type: type of correlation to be used 
%           - Pearson (default)
%           - Spearman
%           - Kendall
% OUTPUT:
%   c   : averaged cross-validated correlation
%   cv  : 1xN vector of cross-validated correlations
%   d   : averaged non-crossvalidaed correlations 
%   dv  : 1xN vector of non-cross-validated correlations

[N,K] = size(Y);
cv    = zeros(1,N);
partition = 1:N;
type = 'Pearson';
vararginoptions(varargin,{'type'});

part = unique(partition);  % unique partitions
for i=1:length(part)
    trainI = (partition~=part(i));
    testI  = (partition==part(i));   
    cv(i)  = corr(mean(Y(trainI,:))',Y(testI,:)','type',type);
    dv(i)  = corr(mean(Y)',Y(testI,:)','type',type);
end;
c=mean(cv);
d=mean(dv); 
