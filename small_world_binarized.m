function [cc_cpl_ratio,cc_ratio,cpl_ratio] = small_world_binarized(connect_mat,sim_num,real_cc,real_cpl)

% small_world.m
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
    sim=randmio_und(connect_mat,rand_network_iter);
    sim_cc(i)=mean(clustering_coef_bu(sim));
    sim_dist=distance_bin(sim);
    sim_cpl(i)=charpath(sim_dist);
end

cc_ratio=mean(real_cc./sim_cc);
cpl_ratio=mean(real_cpl./sim_cpl);
cc_cpl_ratio=cc_ratio/cpl_ratio;