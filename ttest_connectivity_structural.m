function [mcc_p,cpl_p,meloc_p,eglob_p,cc_ratio_p,cpl_ratio_p,sw_p,mod_p,strength_ps,cc_ps,bc_ps,eloc_ps,ereg_ps] = ttest_connectivity_structural(strengths,cc,mcc,cpl,bc,ebc,eloc,meloc,ereg,eglob,small_worldness,cc_ratio,cpl_ratio,modularity,participation_thresh,num_regions,num_subjs,num_group1)
% ttest_connectivity.m
% JB 1/2011
% Run two sample, one-tailed t-tests for connectivity values from connect_metrics.m

%% User Variables
% These first two values are passed in as arguments currently
%num_subjs=20; % total number of subjects
%num_group1=20; % number of subjects in group 1
start_group2=num_group1+1;
global_alpha=.01; % alpha for global measures
node_alpha=.001; % alpha for node measures
edge_alpha=.0001; % alpha for edge measures
direction='both'; % right for group x > y, left for y > x
plot_metrics=true;
fdr_correct=false; % 1 to do False Discovery Rate correction
               % only set up for regional metrics
fdr_q_level=.05;

%% Statistical Tests
% Regional measures
for i=1:num_regions
    [strength_ttest(i,1),strength_ttest(i,2)]=ttest2(strengths(i,1:num_group1),strengths(i,start_group2:num_subjs),node_alpha,direction);
    [cc_ttest(i,1),cc_ttest(i,2)]=ttest2(cc(i,1:num_group1),cc(i,start_group2:num_subjs),node_alpha,direction);
    [bc_ttest(i,1),bc_ttest(i,2)]=ttest2(bc(i,1:num_group1),bc(i,start_group2:num_subjs),node_alpha,direction);
    [eloc_ttest(i,1),eloc_ttest(i,2)]=ttest2(eloc(i,1:num_group1),eloc(i,start_group2:num_subjs),node_alpha,direction);
    [ereg_ttest(i,1),ereg_ttest(i,2)]=ttest2(ereg(i,1:num_group1),ereg(i,start_group2:num_subjs),node_alpha,direction);
    [participation_thresh_ttest(i,1),participation_thresh_ttest(i,2)]=ttest2(participation_thresh(i,1:num_group1),participation_thresh(i,start_group2:num_subjs),node_alpha,direction);
    
    %for j=1:num_regions
        %[ebc_ttest(i,j),ebc_pval(i,j)]=ttest2(ebc(i,j,1:num_group1),ebc(i,j,start_group2:num_subjs),edge_alpha,direction);
        %[estrength_ttest(i,j),estrength_pval(i,j)]=ttest2(connect_mats(i,j,1:num_group1),connect_mats(i,j,start_group2:num_subjs),edge_alpha,direction);
    %end
end

% Get regional p-values
strength_ps=strength_ttest(:,2);
cc_ps=cc_ttest(:,2);
bc_ps=bc_ttest(:,2);
eloc_ps=eloc_ttest(:,2);
ereg_ps=ereg_ttest(:,2);

if fdr_correct
    strength_ttest(:,1)=strength_ttest(:,2)<=fdr(strength_ttest(:,2),fdr_q_level);
    cc_ttest(:,1)=cc_ttest(:,2)<=fdr(cc_ttest(:,2),fdr_q_level);
    bc_ttest(:,1)=bc_ttest(:,2)<=fdr(bc_ttest(:,2),fdr_q_level);
    eloc_ttest(:,1)=eloc_ttest(:,2)<=fdr(eloc_ttest(:,2),fdr_q_level);
    ereg_ttest(:,1)=ereg_ttest(:,2)<=fdr(ereg_ttest(:,2),fdr_q_level);
end

% Mean clustering coefficient
[mcc_ttest,mcc_p]=ttest2(mcc(1:num_group1),mcc(start_group2:num_subjs),global_alpha,direction);
disp(sprintf('mean clustering coefficient p-value: %2.4f',mcc_p))

% Characteristic path length
% Note: directionality is switched here because we are interested in
% whether a group has a lower CPL, indicating greater efficiency
[cpl_ttest,cpl_p]=ttest2(cpl(start_group2:num_subjs),cpl(1:num_group1),global_alpha,direction);
disp(sprintf('characteristic path length p-value: %2.4f',mcc_p))

% Mean local efficiency
[meloc_ttest,meloc_p]=ttest2(meloc(1:num_group1),meloc(start_group2:num_subjs),global_alpha,direction);
disp(sprintf('mean local efficiency p-value: %2.4f',mcc_p))

% Global efficiency
[eglob_ttest,eglob_p]=ttest2(eglob(1:num_group1),eglob(start_group2:num_subjs),global_alpha,direction);
disp(sprintf('global efficiency p-value: %2.4f',mcc_p))

% Normalized cc and cpl
[cc_ratio_ttest,cc_ratio_p]=ttest2(cc_ratio(start_group2:num_subjs),cc_ratio(1:num_group1),global_alpha,direction);
disp(sprintf('normalized clustering coefficient p-value: %2.4f',mcc_p))
[cpl_ratio_ttest,cpl_ratio_p]=ttest2(cpl_ratio(start_group2:num_subjs),cpl_ratio(1:num_group1),global_alpha,direction);
disp(sprintf('normalized characteristic path length p-value: %2.4f',mcc_p))

% Small worldness
[sw_ttest,sw_p]=ttest2(small_worldness(1:num_group1),small_worldness(start_group2:num_subjs),global_alpha,direction);
disp(sprintf('small worldness p-value: %2.4f',mcc_p))

% Modularity
[mod_ttest,mod_p]=ttest2(modularity(1:num_group1),modularity(start_group2:num_subjs),global_alpha,direction);
disp(sprintf('modularity Q p-value: %2.4f',mcc_p))

% Unthresholded modularity
[Q_unthresh_ttest,Q_unthresh_p]=ttest2(Q_unthresh(1:num_group1),Q_unthresh(start_group2:num_subjs),global_alpha,direction);
disp(sprintf('unthresholded modularity Q p-value: %2.4f',mcc_p))

%% Plotting
if plot_metrics == 1
    subplot(3,2,1)
    imagesc(strength_ttest)
    title('Node strength t-test')
    subplot(3,2,2)
    imagesc(cc_ttest)
    title('Node clustering coefficient t-test')
    subplot(3,2,3)
    imagesc(bc_ttest)
    title('Node betweenness centrality t-test')
    subplot(3,2,4)
    imagesc(eloc_ttest)
    title('Node local efficiency t-test')
    subplot(3,2,5)
    imagesc(ereg_ttest)
    title('Node regional efficiency t-test')
    subplot(3,2,6)
    imagesc(participation_thresh_ttest)
    title('Node participation t-test')
    %subplot(3,2,5)
    %imagesc(ebc_ttest)
    %title('Edge betweenness t-test')
    %subplot(3,2,6)
    %imagesc(ebc_pval)
    %title('Edge betweenness p-values')  
end