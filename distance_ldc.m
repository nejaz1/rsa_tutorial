function [d,conw,sigma]=distance_ldc(Y,X,C,partition,varargin);
% function [d,conw]=distance_ldc(Y,X,C,partition,varargin);
% Calculates the crossvalidated squared Eucledian distances 
% On pre-whitened data this gives you the numerator of the LDT distance
% Trial on 
% INPUT:
%  Y          : Raw data in N x P matrix
%  X          : design matrix (NxQ matrix)
%  C          : CxQ labels / contrast denoting which betas are which categories
%  partition  : 1xQ integer value that indicates the partition for crossvalidation (typically run number)
% VARARGIN:
%  'C0',C0                  : Nuisance contrast to be removed before classification (often the
%                               run effect) (QxK) of nuisance effects
% OUTPUT:
%   d         : Average distance values for each contrast across folds (1xC)
%   ts        : t-values for each individual fold (numPartxC) 
%   ps        : p-values for each individual fold (numPartxC) 
% Get the unique Partitions

[part,~,partition]=unique(partition);
numPart   = length(part); 
[N,P]=size(Y); 

for i=1:numPart
    Ia     =  (partition~=part(i));                         % Indices for training set
    Ib     =  (partition==part(i));                         % Indices for testing set
    
    Ya     =  Y(Ia,:);
    Yb     =  Y(Ib,:);
    Xa     =  X(Ia,:);
    Xb     =  X(Ib,:);
    
    betaA =  pinv(Xa)*Ya;                    %-Parameter estimates
    betaB =  pinv(Xb)*Yb;                    %-Parameter estimates
    w      =  (C*betaA);                        % Calculate the classifier w=(muK1-muK2)'*inv(Sw)
    conw(i,:) =  sum((C*betaB).*w,2)/P;                 % project beta's on contrast in questions
end;
d=mean(conw,1); 

