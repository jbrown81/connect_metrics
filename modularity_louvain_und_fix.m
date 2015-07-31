function [q_fix]=modularity_louvain_und_fix(W,ci_fix)
%MODULARITY_LOUVAIN_UND     Optimal community structure and modularity for
%fixed community assignment

W=+W;                                       %convert from logical
n=length(W);                                %number of nodes
s=sum(W(:));                                %weight of edges

n=max(ci_fix);                              %new number of modules
W1=zeros(n);                            %new weighted matrix
for i=1:n
    for j=i:n
        w=sum(sum(W(ci_fix==i,ci_fix==j)));     %pool weights of nodes in same module
        W1(i,j)=w;
        W1(j,i)=w;
    end
end
W=W1;

q_fix=sum(diag(W))/s-sum(sum((W/s)^2));