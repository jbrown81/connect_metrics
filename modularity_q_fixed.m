function [q_fixed] = modularity_q_fixed(W,Ci0)

% calculate Q for a network (unthresholded) given a fixed community
% assignment

M=Ci0;

qtype='sta';

n=length(W);                                                %number of nodes/modules

[dum dum M] = unique(M(:).');                           %align module indices

W0= W.*(W>0);                                               %positive weights matrix
W1=-W.*(W<0);                                               %negative weights matrix
s0=sum(W0(:));                                              %positive sum of weights
s1=sum(W1(:));                                              %negative sum of weights
Knm0=zeros(n,n);                                            %positive node-to-module degree
Knm1=zeros(n,n);                                            %negative node-to-module degree
for m=1:max(M)                                              %loop over modules
    Knm0(:,m)=sum(W0(:,M==m),2);
    Knm1(:,m)=sum(W1(:,M==m),2);
end
Kn0=sum(Knm0,2);                                            %positive node degree
Kn1=sum(Knm1,2);                                            %negative node degree
Km0=sum(Knm0,1);                                            %positive module degree
Km1=sum(Knm1,1);                                            %negative module degree

switch qtype
    case 'smp';  d0 = 1/s0;       d1 = 1/s1;                %dQ = dQ0/s0 - dQ1/s1;
    case 'gja';  d0 = 1/(s0+s1);  d1 = 1/(s0+s1);           %dQ = (dQ0 - dQ1)/(s0+s1);
    case 'sta';  d0 = 1/s0;       d1 = 1/(s0+s1);           %dQ = dQ0/s0 - dQ1/(s0+s1);
    case 'pos';  d0 = 1/s0;       d1 = 0;                   %dQ = dQ0/s0;
    case 'neg';  d0 = 0;          d1 = 1/s1;                %dQ = -dQ1/s1;
    otherwise; error('qtype unknown');
end
if ~s0                                                      %adjust for absent positive weights
    s0=1;
    d0=0;
end
if ~s1                                                      %adjust for absent negative weights
    s1=1;
    d1=0;
end

[dum dum M]=unique(M(:).');                                 %realign module indices
%M=M'; % added by JB, 08/2014
%compute modularity
m = M(ones(1,n),:);
Q0 = (W0-(Kn0*Kn0.')/s0).*(m==m.');
Q1 = (W1-(Kn1*Kn1.')/s1).*(m==m.');
q_fixed = d0*sum(Q0(:)) - d1*sum(Q1(:));