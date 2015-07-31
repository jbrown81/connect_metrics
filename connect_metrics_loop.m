% connect_metrics_loop.m
% Jesse Brown 02/28/2014

% Call connect_metrics_func in a loop with different sparsities and store
% results

%sparseness_range=.05:.05:.5;
sparseness_range=.05:.05:.1;
count=0;

for i=sparseness_range
    count=count + 1;
    [strengths_all(:,:,count),cc_all(:,:,count),mcc_all(:,count),cpl_all(:,count),bc_all(:,:,count),eloc_all(:,:,count),meloc_all(:,count),ereg_all(:,:,count),eglob_all(:,count),small_worldness_all(:,count),cc_ratio_all(:,count),cpl_ratio_all(:,count),modularity_all(:,count),participation_thresh_all(:,:,count),Ci_all(:,:,count),n_comps_all(:,count),min_edge_weights_all(:,count)]=connect_metrics_func(i,false);
    disp(i)
end