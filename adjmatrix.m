function thresh = adjmatrix(WavCorrMatrix,cost,varargin)
% This function takes a symmetric matrix of association (e.g. a wavelet 
% correlation matrix), "WavCorrMatrix" and binarizes it at some cost, 
% outputing the adjacency matrix, "thresh". The cost can be a value between 
% 0 and 1, where a cost of 0.1 means that the 10% strongest connections are 
% set to 1 while the 90% weakest connections are set to 0. This function 
% uses the absolute values of the associations. If you would like to change 
% this, edit line 14.

% Modified by JB 4/13/2012: takeabs argument required boolean

% By default, do not take absolute value
if nargin > 2
    takeabs = varargin{1};
else
    takeabs = false;
end

M = numel(WavCorrMatrix(1,:));
% Checks to see if you diagonals have been set to zero yet.
if WavCorrMatrix(1,1)~=0;
    WavCorrMatrix = tril(WavCorrMatrix,-1)+triu(WavCorrMatrix,+1);
end

if takeabs
    [Wavsort ix] = sort(reshape(tril(abs(WavCorrMatrix),-1),1,M*M),'descend');
else
    [Wavsort ix] = sort(reshape(tril(WavCorrMatrix,-1),1,M*M),'descend');
end

thresh = zeros(M,M);
thresh(ix(1:round(cost*M*(M-1)/2))) = ones;
thresh = thresh + thresh';