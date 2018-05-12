function [pd,muK,Sw]=distance_mahalanobis(y,c,rn,varargin);
% function [C,muK,Sw]=mva_distance(y,c,varargin);
% Calculates a numclass x numclass distance matrix
% Based on
% INPUT:
%   y: PxN matrix (P voxels, N trials)
%   c: category
%   rn: crossvalidation set (ignored)
% VARARGIN:
%   'regular',0.01              : Regularization for the Sw matrix
%   'norm','none'               : Normalization of the mean patterns
%   'pooled','yes'              : Pooled or not pooled Sw? 
% OUTPUT:
% C: pdist shape of distance matrix (use squareform to get numclass x
%                                       numclass matrix)
regular=0.01; 
norm='none'; 
vararginoptions(varargin,{'distance','regular','norm'});
classes=unique(c);
numclasses=length(classes); % Number of classes in the data set
[P,N]=size(y);

muK=zeros([P numclasses]);      % means
Sw=zeros(P,P);          % Within class variability
for i=1:numclasses
    j=find(classes(i)==c);
    n(i) = length(j);      % number of sampels per category
    muK(:,i) = sum(y(:,j),2)/n(i);                             % get the Cluster means
    res = bsxfun(@minus,y(:,j),muK(:,i));
    Sw = Sw+(res*res')/n(i);                         % Estimate common covariance matrix
end;
Sw=Sw+eye(P)*mean(diag(Sw))*regular;        % Regularize the Sw matrix 

% calculate the mahalabobis distance 
pd=pdist(muK','mahalanobis',Sw)';
