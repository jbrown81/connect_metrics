function [out_mat,min_val_orig,max_val_orig,edge_count_orig,min_val_thr,max_val_thr,edge_count_thr] = threshmat(in_mat,sparseness,varargin)

% By default, do not take absolute value
if nargin > 2
    takeabs = varargin{1};
else
    takeabs = false;
end

prct_cutoff=100-(sparseness*100);

if takeabs
    in_mat=abs(in_mat);
end

cutoff=prctile(in_mat(:),prct_cutoff);
min_val_orig=min(in_mat(:));
max_val_orig=max(in_mat(:));
edge_count_orig=sum(in_mat(:)~=0);

% makes sure floor values (ie 0) don't get set to 1
out_mat=(double(in_mat>=cutoff).*double(in_mat~=min_val_orig).*in_mat);
min_val_thr=min(out_mat(find(out_mat~=0)));
max_val_thr=max(out_mat(:));
edge_count_thr=sum(out_mat(:)~=0);