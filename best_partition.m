function [best_ci,max_q]=best_partition(W,n_runs,sign)

% best_partition.m
% Written by Jesse Brown
% determines optimal modularity partitioning based on multiple iterations
% of modularity_louvain_und function from Brain Connectivity Toolbox

%   Input:      W,      undirected (weighted or binary) connection matrix.
%               n,      number of iterations
%               sign,   true if W contains positive and negative
%               weights

in_mat=W; % matrix on which to determine module stability
[r,c]=size(W);

ci=zeros(r,n_runs);
VIn=zeros(n_runs,n_runs);
MIn=zeros(n_runs,n_runs);

if sign
    for i=1:n_runs
        [ci(:,i) q(i)]=modularity_louvain_und_sign(in_mat);
    end
else
    for i=1:n_runs
        [ci(:,i) q(i)]=modularity_louvain_und(in_mat);
    end

end

% calculate the similarity of each partition
for j=1:n_runs
    for k=1:n_runs
        [VIn(j,k) MIn(j,k)]=partition_distance(ci(:,j),ci(:,k));
    end
end

% meta modularity: cluster the different partitions based on their
% similarity
meta_ci=modularity_louvain_und(MIn);
[p ind]=sort(meta_ci');

num=length(unique(p)); % number of different configurations detected (approximate)

clear p_agg q_mean
for i=1:num
    count=sum(p==i);
    cur=ind(find(p==i)); % get set of identical partitions

    sub_MIn=MIn(cur,cur); % get the mean and standard deviation of their similarity
    sub_MIn_mean=mean(sub_MIn(:));
    sub_MIn_std=std(sub_MIn(:));

    % aggregate measure of the goodness of a parition
    % takes the number of times this paritioning was achieved
    % multiplies by the mean similiarity of the paritions
    % divides by the variance of the partition similarity
    % CURRENTLY BASED ON PARITION CONSISTENCY
    
    % SWITCHED THIS, CHANGE BACK IF NEEDED
    %p_agg(i)=(count*sub_MIn_mean)/sub_MIn_std;
    p_agg(i)=(count*mean(q(cur)))/std(q(cur));
    
    q_mean(i) = mean(q(cur));
end

% Find the set of partitions with the best score
[m best_p_num]=max(p_agg);

best_ci_nums=ind(find(p==best_p_num));
best_cis=ci(:,best_ci_nums);

% Among those, choose the one with the maximum q-value
max_q=max(q(best_ci_nums));
max_q_num=find(q==max_q);
best_ci=ci(:,max_q_num(1));