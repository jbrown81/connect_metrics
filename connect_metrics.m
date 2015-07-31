% connect_metrics.m
% Jesse Brown 2011-2015
% jesse.brown@ucsf.edu
% Connectivity metric calculation based on Brain Connectivity Toolbox
% http://sites.google.com/a/brain-connectivity-toolbox.net/bct/Home
% Set up to run for list of subjects
% Input files are space or tab-delimited, SYMMETRIC connectivity matrices

%% User Variables
% The name of the directory where each subject's connectivity matrix is stored
file_path='/data/mridata/jbrown/psp_nc/fc_mats';
files={'/data5/patientNIC/12813/12813_20110809/rsfmri/processedfmri_TRCNnSFmDI/matrix/rMT_27_fc_mat_covars.txt' '/data5/patientNIC/12813/12813_20120305/rsfmri/processedfmri_TRCNnSFmDI/matrix/rMT_27_fc_mat_covars.txt'};
addpath(genpath('/data/mridata/jbrown/scripts/connect_metrics_share'));

structural=false; % Set to true for a structural connectivity analysis, false otherwise
functional=true;% Set to true for a functional connectivity analysis, false otherwise
keep_mats=true; % true to keep each subject's matrix (connect_mat, length_mat, stat_mat)

% Main variables
num_regions=27; % the number of regions
do_connect_metrics=true; % true to calculate network metrics, false to just load matrices
sparseness=.2; % 0 to do no thresholding, otherwise choose a threshold % between 0 and 1 (eg .2 means keep 20% strongest edges)
binarize=false; % false to use weighted matrix, true to use tresholded binarized matrix
small_world_iters=10; % number of simulated random networks to generate, 10 or more is best
mod_iters=20; % number of different modularity partitions to compute to find highest Q, 20 or more is best
ci_fixed=[]; % pre-specified modularirty community index for each node; leave empty if unknown

% Variables ONLY for structural connectivity matrix, usually from DTI tractography
if structural
    length_scaling=false; % scale weight distance mapping by average length between regions
    if length_scaling
        lengthmat_files={'p1_lengthmat.txt' 'p2_lengthmat.txt' 'c1_lengthmat.txt' 'c2_lengthmat.txt'}; % average fiber length matrix filenames
    end
    weight_stat_scaling=false; % scale weight by statistic (eg FA)
    if weight_stat_scaling
        statmat_files={'p1_statmat.txt' 'p2_statmat.txt' 'c1_statmat.txt' 'c2_statmat.txt'}; % average FA/ADC/... matrix filenames
    end
    fiber_density=false; % false to use raw fiber counts, true to scale by region volume
    if fiber_density % if fiber_density, list of files where each file is single column with volume for each ROI on rows
        roi_volume_files={'p1_vols.txt' 'p2_vols.txt' 'c1_vols.txt' 'c2_vols.txt'};
    end
    fiber_thresh=0; % only keep connections with greater than this number of fibers
end

% Variables ONLY for functional connectivity matrix, usual from fMRI
% NOTE: if abs_value and min_shift are both set to false, and sparseness
% value is high enough (eg enough edges are kept) to include negative weights, 
% program will crash
if functional
    abs_value=false; % if true, take abs value of connectivity matrix
    min_shift=false; % shift distribution to be all positive
    r_to_z=true; % perform r-to-Z transformation
    apply_distance_thresh=false; % zero out connections between ROIs with Euclidean distance less than specified distance apart
    distance_thresh=20; % the distance cutoff in mm
    centers_file='/Users/jessebrown/Desktop/centers_149.txt'; % a file with the center of mass coordinate in mm
    if apply_distance_thresh
        roi_centers=load(centers_file);
        centers_dist_mat=centers_dist(roi_centers);
    end
    network_deconv=false; % apply network deconvolution to raw correlation matrix (Feizi et al Nature Biotech 2013, http://www.nature.com/nbt/journal/v31/n8/full/nbt.2635.html)
end

% For statistics; associated script ttest_connectivity.m will run two-group
% t-test on resultant connectivity metrics if desired
do_twogroup_statistics=false;
num_group1=17; % number of subjects in group 1


%% Variable initialization
[a num_subjs]=size(files);
cc=zeros(num_regions,num_subjs);
strengths=zeros(num_regions,num_subjs);
bc=zeros(num_regions,num_subjs);
ebc=zeros(num_regions,num_regions,num_subjs);
eloc=zeros(num_regions,num_subjs);
ereg=zeros(num_regions,num_subjs);


%% Set up graphs (connectivity matrices)
for i=1:num_subjs
    % Load files   
    connectmat_file=sprintf('%s',files{i});
    connect_mat=load(connectmat_file);
    connect_mat_raw=connect_mat;

    % zero  matrices
    connect_mat(find(isnan(connect_mat)))=0; % set NaNs to zero
    connect_mat=connect_mat.*~eye(num_regions,num_regions); % set diagonal elements to zero
    
    %% Set up connectivity matrices (graphs)
    % Structure-specific matrix setup
    if structural
        % Structure specific files
        if length_scaling
            load(sprintf('%s/%s/%s.txt',file_path,files{i},lengthmat_file));
            length_mat=eval(lengthmat_file);
            length_mat=length_mat.*(connect_mat>fiber_thresh);
        end
        
        if weight_stat_scaling
            load(sprintf('%s/%s/%s.txt',file_path,files{i},statmat_file));
            stat_mat=eval(statmat_file);
            stat_mats(:,:,i)=stat_mat;
            stat_mat=stat_mat.*(connect_mat>fiber_thresh);
        end
        
        if fiber_thresh > 0
            connect_mat=connect_mat.*(connect_mat>fiber_thresh);
        end
    
        % Create fiber density matrix
        % Divide fiber count for each region by volume
        if fiber_density
            vols=load(sprintf('%s/%s',file_path,roi_volume_files{i}));
            vols=vols(1:num_regions);
            volsums_mat=zeros(num_regions,num_regions);
            for v1=1:num_regions
                for v2=1:num_regions
                    volsums_mat(v1,v2)=mean(vols(v1)+vols(v2));
                end
            end
            connect_mat=connect_mat./volsums_mat;
            connect_mat(find(isnan(connect_mat)))=0;
        end
    
        % Scale connection matrix by statistic matrix
        if weight_stat_scaling
            connect_mat = connect_mat.*stat_mat;
        end
    end
       
    % Get the percent connectedness for the unthresholded network
    percent_connected(i)=length(find(connect_mat(:)))/(num_regions*(num_regions-1));
    
    % Create thresholded weighted matrix
    if sparseness > 0
        if structural
            [connect_mat,min_val_orig(i),max_val_orig(i),edge_count_orig(i),min_val_thr(i),max_val_thr(i),edge_count_thr(i)]=threshmat(connect_mat,sparseness,false);
            if length_scaling
                length_mat=length_mat.*threshmat(connect_mat,sparseness);
                % Get total length of fiber for all connections
                total_fiber_length(i)=sum(sum(connect_mat.*length_mat));
            end
            if weight_stat_scaling
                stat_mat=stat_mat.*threshmat(connect_mat,sparseness);
                % Total cost account for fiber number, length, and FA
                total_fiber_cost(i)=sum(sum(connect_mat.*length_mat.*stat_mat));
            end
        end
        
        if functional
            % prepare edges
            if r_to_z
                connect_mat=fisherz(connect_mat);
            end
            
            if apply_distance_thresh
                connect_mat=connect_mat.*(centers_dist_mat>distance_thresh);
            end
            
            if network_deconv
                connect_mat=ND(connect_mat);
                connect_mats_ND_unthresh(:,:,i)=connect_mat;
            end

            connect_mat_unthr=connect_mat;
            
            % apply thresholding
            if abs_value
                [connect_mat,min_val_orig(i),max_val_orig(i),edge_count_orig(i),min_val_thr(i),max_val_thr(i),edge_count_thr(i)]=threshmat(connect_mat,sparseness,true);
            elseif min_shift
                % shift the weight distribution so the min value is zero
                connect_mat=connect_mat+abs(min(connect_mat(:)));
                [connect_mat,min_val_orig(i),max_val_orig(i),edge_count_orig(i),min_val_thr(i),max_val_thr(i),edge_count_thr(i)]=threshmat(connect_mat,sparseness,false);
            else
                [connect_mat,min_val_orig(i),max_val_orig(i),edge_count_orig(i),min_val_thr(i),max_val_thr(i),edge_count_thr(i)]=threshmat(connect_mat,sparseness,false);
            end
        end
    end

    % Binarize the weighted matrix
    connect_mat_bin=double(connect_mat~=0);
        
    % Keep individual raw matrices in memory
    if keep_mats
        connect_mats(:,:,i)=connect_mat;
        connect_mats_raw(:,:,i)=connect_mat_raw;
        if structural
            if length_scaling
                len_mats(:,:,i)=length_mat;
            end
            if weight_stat_scaling
                stat_mats(:,:,i)=stat_mat;
            end
        end
    end

    
    %% Calculate connectivity metrics
    if do_connect_metrics
    if binarize
        % Strength (same for binarized/weighted matrix)
        strengths(:,i)=strengths_und(connect_mat_bin)';
        
        % Clustering coefficient
        cc(:,i)=clustering_coef_bu(connect_mat_bin);
        
        % Mean clustering coefficient
        mcc(i)=mean(cc(:,i));
        
        % Nodal betweenness centrality
        bc(:,i)=betweenness_bin(connect_mat_bin);
        
        % Edge betweenness centrality
        %ebc(:,:,i)=edge_betweenness_bin(connect_mat_bin);
        
        % Create distance matrix       
        distance_mat=distance_bin(connect_mat_bin);
        inv_distance_mat=1./distance_mat.*(~eye(num_regions)); % zero diagonal
        inv_distance_mat(find(isnan(inv_distance_mat)))=0; % set Nan's to zero
        
        % Characteristic path length, along with efficiency, node
        % eccentricity, radius, and diameter
        [cpl(i),efficiency(i),ecc(:,i),radius(i),diameter(i)]=charpath(distance_mat);
                
        % Nodal local efficiency 
        eloc(:,i)=efficiency_bin(connect_mat_bin,1);
        
        % Mean local efficiency
        meloc(i)=mean(eloc(:,i));
        
        % Regional efficiency
        ereg(:,i)=sum(inv_distance_mat)/num_regions;
        
        % Global efficiency, excluding main diagonal
        eglob(i)=efficiency_bin(connect_mat_bin);
        
        % Small worldness
        [cc_cpl_ratio,cc_ratio(i),cpl_ratio(i)]=small_world_binarized(connect_mat_bin,small_world_iters,mcc(i),cpl(i));
        small_worldness(i)=cc_cpl_ratio;
        
        % Modularity (Louvain)
        %[Ci modularity(i)]=modularity_louvain_und(connect_mat_bin);
        [Ci(:,i) modularity(i)]=best_partition(connect_mat,10,false);
        
        % Participation
        participation_thresh(:,i)=participation_coef(connect_mat_bin,Ci(:,i));
        
        % Within-module Z-score
        modz_thresh(:,i)=module_degree_zscore(connect_mat_bin,Ci(:,i));
        
        % Participation and Within-module Z-score (based on pre-specified modularity partition)
        if ci_fixed
            [partic_thr_fix(:,i)]=participation_coef(connect_mat_bin, ci_fixed);
            [within_modz_thr_fix(:,i)]=module_degree_zscore(connect_mat_bin, ci_fixed);
        end
        
        % Unthresholded weighted measures for weighted functional networks
        if functional
            %[Ci0]=modularity_louvain_und_sign(connect_mat_unthr);
            [Ci0]=best_partition(connect_mat_unthr,mod_iters,true);
            [Ci_unthresh(:,i) Q_unthresh(i)] = modularity_finetune_und_sign(connect_mat_unthr,'sta',Ci0);
            [strengths_pos(:,i) strengths_neg(:,i)]=strengths_und_sign(connect_mat_unthr);
            [participation_pos(:,i) participation_neg(:,i)]=participation_coef_sign(connect_mat_unthr, Ci_unthresh(:,i));
            [diversity_pos(:,i) diversity_neg(:,i)]=diversity_coef_sign(connect_mat_unthr, Ci_unthresh(:,i));
            modZ(:,i)=module_degree_zscore(connect_mat_unthr,Ci_unthresh(:,i));
            [within_modz_pos(:,i) within_modz_neg(:,i)]=module_degree_zscore_sign(connect_mat_unthr, Ci_unthresh(:,i));
            [without_modz_pos(:,i) without_modz_neg(:,i)]=without_module_degree_zscore_sign(connect_mat_unthr, Ci_unthresh(:,i));
            
            if ci_fixed
                [partic_pos_unthr_fix(:,i) partic_neg_unthr_fix(:,i)]=participation_coef_sign(connect_mat_unthr, ci_fixed);
                [within_modz_pos_unthr_fix(:,i) within_modz_unthr_neg_fix(:,i)]=module_degree_zscore_sign(connect_mat_unthr, ci_fixed);
            end
        end
        
        % Assortativity is correlation coefficient between the degrees of
        % all nodes on two opposite ends of a link.
        assortativity(i)=assortativity_bin(connect_mat_bin,0);
        
        % Number of components
        [comps,comp_sizes]=get_components(connect_mat_bin);
        n_comps=length(comp_sizes);
    else
        % Strength (same for binarized/weighted matrix)
        strengths(:,i)=strengths_und(connect_mat)';
        
        % Clustering coefficient
        cc(:,i)=clustering_coef_wu(connect_mat);
        
        % Mean clustering coefficient
        mcc(i)=mean(cc(:,i));

        % Inverse weight matrix, to create a weight to distance mapping
        % Typical mapping is weight inversion
        inv_connect_mat=1./connect_mat;
        
        % If scale weight distance mapping by physical length
        if structural
            if length_scaling
                weight_distance_mapping=(inv_connect_mat).*length_mat;
                % Multiplying by lengths changes Inf's to Nan's, so reset them to
                % Inf, which distance_wei function can handle
                weight_distance_mapping(find(isnan(weight_distance_mapping)))=Inf;
            else
                weight_distance_mapping=(inv_connect_mat);
            end
        else
            weight_distance_mapping=(inv_connect_mat);
            weight_distance_mapping(find(isinf(weight_distance_mapping)))=0;
        end
        
        % Create distance matrix
        % all_paths structure stores the shortest path between any pair of nodes
        %distance_mat=distance_wei(weight_distance_mapping);
        [distance_mat,n_steps_mat,all_paths]=distance_wei_paths(weight_distance_mapping);
        inv_distance_mat=(1./distance_mat).*(~eye(num_regions)); % zero diagonal
        inv_distance_mat(find(isnan(inv_distance_mat)))=0; % set Nan's to zero
        
        % Characteristic path length, along with efficiency, node
        % eccentricity, radius, and diameter
        [cpl(i),efficiency(i),ecc(:,i),radius(i),diameter(i)]=charpath(distance_mat);
        
        cpl(i)=charpath(distance_mat);
        
        % Nodal Betweenness centrality
        bc(:,i)=betweenness_wei(weight_distance_mapping);
    
        % Edge betweenness centrality  
        %ebc(:,:,i)=edge_betweenness_wei(weight_distance_mapping);
        
        % Nodal local efficiency
        if structural && length_scaling
            eloc(:,i)=local_efficiency_wu(connect_mat,length_mat);
        else
            eloc(:,i)=efficiency_wei(connect_mat,1);
        end
            
        % Mean local efficiency
        meloc(i)=mean(eloc(:,i));
        
        % Regional efficiency for a weighted matrix just the row-wise mean of 
        % the inverse distance matrix
        ereg(:,i)=sum(inv_distance_mat)/(num_regions-1);
    
        % Global efficiency is grand mean of inverse distance matrix, 
        % excluding main diagonal
        eglob(i)=efficiency_wei(connect_mat);
        
        % Small worldness
        if structural
            if length_scaling
                [cc_cpl_ratio,cc_ratio(i),cpl_ratio(i)]=small_world_weighted_distance(connect_mat,length_mat,small_world_iters,mcc(i),cpl(i));
            else
                [cc_cpl_ratio,cc_ratio(i),cpl_ratio(i)]=small_world_weighted(connect_mat,small_world_iters,mcc(i),cpl(i));
            end
        else
            [cc_cpl_ratio,cc_ratio(i),cpl_ratio(i)]=small_world_weighted(connect_mat,small_world_iters,mcc(i),cpl(i));
        end
        small_worldness(i)=cc_cpl_ratio;
        
        % Modularity (Louvain)
        %[Ci modularity(i)]=modularity_louvain_und(connect_mat);
        [Ci(:,i) modularity(i)]=best_partition(connect_mat,mod_iters,false);
        
        % Participation coefficent (based on modularity)
        participation_thresh(:,i)=participation_coef(connect_mat,Ci(:,i));
        
        % Within-module Z-score (based on modularity)
        modz_thresh(:,i)=module_degree_zscore(connect_mat,Ci(:,i));
        
        % Participation and Within-module Z-score (based on pre-specified modularity partition)
        if ci_fixed
            [partic_thr_fix(:,i)]=participation_coef(connect_mat, ci_fixed);
            [within_modz_thr_fix(:,i)]=module_degree_zscore(connect_mat, ci_fixed);
        end
        
        % Unthresholded weighted measures for weighted functional networks
        if functional
            %[Ci0]=modularity_louvain_und_sign(connect_mat_unthr);
            [Ci0]=best_partition(connect_mat_unthr,mod_iters,true);
            [Ci_unthresh(:,i) Q_unthresh(i)] = modularity_finetune_und_sign(connect_mat_unthr,'sta',Ci0);
            [strengths_pos(:,i) strengths_neg(:,i)]=strengths_und_sign(connect_mat_unthr);
            [participation_pos(:,i) participation_neg(:,i)]=participation_coef_sign(connect_mat_unthr, Ci_unthresh(:,i));
            [diversity_pos(:,i) diversity_neg(:,i)]=diversity_coef_sign(connect_mat_unthr, Ci_unthresh(:,i));
            modZ(:,i)=module_degree_zscore(connect_mat_unthr,Ci_unthresh(:,i));
            [within_modz_pos(:,i) within_modz_neg(:,i)]=module_degree_zscore_sign(connect_mat_unthr, Ci_unthresh(:,i));
            [without_modz_pos(:,i) without_modz_neg(:,i)]=without_module_degree_zscore_sign(connect_mat_unthr, Ci_unthresh(:,i));
            
            if ci_fixed
                [partic_pos_unthr_fix(:,i) partic_neg_unthr_fix(:,i)]=participation_coef_sign(connect_mat_unthr, ci_fixed);
                [within_modz_pos_unthr_fix(:,i) within_modz_unthr_neg_fix(:,i)]=module_degree_zscore_sign(connect_mat_unthr, ci_fixed);
            end
        end
        
        % Assortativity is correlation coefficient between the degrees of
        % all nodes on two opposite ends of a link.
        assortativity(i)=assortativity_wei(connect_mat_bin,0);
        
        % Number of components
        [comps,comp_sizes]=get_components(connect_mat);
        n_comps(i)=length(comp_sizes);
    end
    end
    
    mat_dim=num_regions*num_regions;
    disp(sprintf('matrix %d/%d done; %d/%d edges intact; final weight range %2.2f-%2.2f',i,num_subjs,edge_count_thr(i),mat_dim,min_val_thr(i),max_val_thr(i)));
end


%% Run between-group statistics
if do_twogroup_statistics
    if structural
        [mcc_p,cpl_p,meloc_p,eglob_p,cc_ratio_p,cpl_ratio_p,sw_p,mod_p,strength_ps,cc_ps,bc_ps,eloc_ps,ereg_ps]=ttest_connectivity_structural(strengths,cc,mcc,cpl,bc,ebc,eloc,meloc,ereg,eglob,small_worldness,cc_ratio,cpl_ratio,modularity,participation_thresh,num_regions,num_subjs,num_group1);
    end
    if functional
        [mcc_p,cpl_p,meloc_p,eglob_p,cc_ratio_p,cpl_ratio_p,sw_p,mod_p,strength_ps,cc_ps,bc_ps,eloc_ps,ereg_ps]=ttest_connectivity_functional(strengths,cc,mcc,cpl,bc,ebc,eloc,meloc,ereg,eglob,small_worldness,cc_ratio,cpl_ratio,modularity,participation_thresh,Q_unthresh,strengths_pos,strengths_neg,participation_pos,participation_neg,num_regions,num_subjs,num_group1);
    end
end