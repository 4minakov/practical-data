% k-means algorithm
% make synthetic data
%
data0 = randn(300,2);
xc0 = [0 0; -3 2; 2 3];
%xc0 = [0 0; -2 2;];
% true cluster centers 
data = [ ...
[data0(:,1) + xc0(1,1), data0(:,2) + xc0(1,2)];...
[data0(:,1) + xc0(2,1), data0(:,2) + xc0(2,2)];...
[data0(:,1) + xc0(3,1), data0(:,2) + xc0(3,2)];...
];
%data = data0;

figure(1)
plot(data(:,1),data(:,2),'r.', 'MarkerSize',12), hold on
plot(xc0(:,1),xc0(:,2),'k+', 'MarkerSize',12,'LineWidth',2)
title('Data');
%
% number of points and dimensions of the data
[np,ndim] = size(data); 
% assume number of clusters
k = 3;
% choose trial centers of clusters
ind = randperm(np,k);
xc  = data(ind,:);
% initialize labels for data points
labels_cur = zeros(np, 1);
labels_prev = ones(np, 1);
% initialize distances between data points and cluster centroids
D = zeros(np,k);
it = 1;
% iterate to find correct centers
while ~isequal(labels_cur,labels_prev) % until labels don't change
    labels_prev=labels_cur; % update labels from previous step
    % calculate squared distances (metric) between data points and
    % trial cluster centers
    for i = 1:k
        D(:,i) = (data(:,1)-xc(i,1)).^2 + (data(:,2)-xc(i,2)).^2;
    end
    % find indices (labels) of points which have the minimum distances
    % to each of the cluster
    % in order to label data points which cluster they belong to
    [~,labels_cur] = min(D,[],2);
    % update trial cluster centers with the actual centeroids of the
    % corresponding set of points
    for i_clust = 1:k
        xc(i_clust,:) = mean(data(labels_cur==i_clust,:),1);
    end
    it = it+1 
    disp(it)
    % plot trial clusters marked with different color
    figure(1),clf
    hold on
    plot(data(:,1),data(:,2),'k.', 'MarkerSize',12)
    plot(data(labels_cur==1,1),data(labels_cur==1,2),'r.', 'MarkerSize',12)
    plot(data(labels_cur==2,1),data(labels_cur==2,2),'g.', 'MarkerSize',12)
    plot(data(labels_cur==3,1),data(labels_cur==3,2),'b.', 'MarkerSize',12)
    plot(xc(:,1),xc(:,2),'ko', 'MarkerSize',12,'LineWidth',2)
    plot(xc0(:,1),xc0(:,2),'k+', 'MarkerSize',12,'LineWidth',2)
    drawnow
    title(['iteration = ',num2str(it)])
    pause(0.5)
end


