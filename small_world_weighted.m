function [cc_cpl_ratio,cc_ratio,cpl_ratio] = small_world_weighted(connect_mat,sim_num,real_mcc,real_cpl)

% small_world_weighted.m
% JB 1/2011
% Input a real (weighted) connectivity matrix
% Calculate simulated networks with equivalent degree
% Calculate ratio of real to simulated clustering coefficient,
% characteristic path length

num_regions=length(connect_mat);

rand_network_iter=5;
sim_mcc=zeros(1,sim_num);
sim_cpl=zeros(1,sim_num);

for i=1:sim_num
    sim=randmio_und_connected(connect_mat,rand_network_iter);
    sim_mcc(i)=mean(clustering_coef_wu(sim));
    
    sim_weight_distance_mapping=1./sim;
    sim_weight_distance_mapping(find(isnan(sim_weight_distance_mapping)))=Inf;
    sim_weight_distance_mapping(find(isinf(sim_weight_distance_mapping)))=0;
    
    sim_distance_mat=distance_wei(sim_weight_distance_mapping);
    sim_inv_distance_mat=(1./sim_distance_mat).*(~eye(num_regions)); % zero diagonal
    sim_inv_distance_mat(find(isnan(sim_inv_distance_mat)))=0; % set Nan's to zero
    
    sim_cpl(i)=charpath(sim_distance_mat);  
end

cc_ratio=mean(real_mcc./sim_mcc);
cpl_ratio=mean(real_cpl./sim_cpl);
cc_cpl_ratio=cc_ratio/cpl_ratio;