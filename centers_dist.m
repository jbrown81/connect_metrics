function[distance_mat]=centers_dist(centers)

% Find the euclidean distance between all coordinates in a connectivity
% network

[r,c]=size(centers);

for i=1:r
    for j=1:r
        x=centers(i,:);
        y=centers(j,:);
        distance_mat(i,j)=sqrt(((x(1)-y(1))^2) + ((x(2)-y(2))^2) + ((x(3)-y(3))^2));
    end
end