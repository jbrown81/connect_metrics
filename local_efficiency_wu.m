function[eloc]=local_efficiency_wu(G,length_mat);
% Jesse Brown, 8/2010
% For each node, find subgraph
% For subgraph, find characteristic path length
% Inverse is local efficiency
% Optional: include length_mat to scale distance matrix by fiber length

[r,c]=size(G);
eloc=zeros(r,1);

if nargin==2
    for x=1:r
        ind=find(G(x,:)>0);
        if length(ind)==0;
            eloc(x)=0;
        else
            % ind=[x ind]; % Not including i, following Latora and Marchiori 2001
            li=length(ind);
            submat=zeros(li,li);
            sublenmat=zeros(li,li);
            for sub1=1:li
                for sub2=1:li
                    submat(sub1,sub2)=G(ind(sub1),ind(sub2));
                    sublenmat(sub1,sub2)=length_mat(ind(sub1),ind(sub2));
                end
            end
            weight_distance_mapping=(1./submat).*sublenmat;
            distance_mat=distance_wei(weight_distance_mapping);
            inv_distance_mat=zero_diagonal(1./distance_mat);
            eloc(x)=sum(sum(inv_distance_mat))/(li*(li-1));
        end
    end
else
    for x=1:r
        ind=find(G(x,:)>0);
        if length(ind)==0;
            eloc(x)=0;
        else
            % ind=[x ind]; % Not including i, following Latora and Marchiori 2001
            li=length(ind);
            submat=zeros(li,li);
            for sub1=1:li
                for sub2=1:li
                    submat(sub1,sub2)=G(ind(sub1),ind(sub2));
                end
            end
            weight_distance_mapping=1./submat;
            distance_mat=distance_wei(weight_distance_mapping);
            inv_distance_mat=zero_diagonal(1./distance_mat);
            eloc(x)=sum(sum(inv_distance_mat))/(li*(li-1));
        end
    end
end