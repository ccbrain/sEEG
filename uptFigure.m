
% user defined clusters
preCluster=who('CLUSTER*');

% user defined no cluster voxels
noCluster=who('NOCLUSTER*');

if ~isempty(preCluster) || ~isempty(noCluster)
    tmpNoCluster=[];
    for i=1:length(noCluster)
        tmpNoCluster=[tmpNoCluster; eval(['NOCLUSTER' num2str(i)])];
    end
    for i=1:length(preCluster)
        tmpNoCluster=[tmpNoCluster; eval(['CLUSTER' num2str(i)])];
    end
    % new corrdinates for clustering, exclusing no-cluster and pre-defined
    % voxels
    newCor = setdiff([x_seed(validIndx) y_seed(validIndx) z_seed(validIndx)],tmpNoCluster,'rows');
    [GS_seed,clusters_seed]=  clustering(nEle-length(preCluster),newCor,ones(size(newCor,1),1));  % use the func from iElectrodes toolbox to give fixed solution
    
    % Current Cluster
    x_clustered=newCor(:,1);
    y_clustered=newCor(:,2);
    z_clustered=newCor(:,3);
    clustered=clusters_seed;
    
    % append user defined voxels
    for i=1:length(preCluster)
        tmpCor=eval(['CLUSTER' num2str(i)]);
        x_clustered=[x_clustered;tmpCor(:,1)];
        y_clustered=[y_clustered;tmpCor(:,2)];
        z_clustered=[z_clustered;tmpCor(:,3)];
        clustered=[clustered;repmat(nEle-length(preCluster)+i,size(tmpCor,1),1)];
        GS_seed=[GS_seed;mean(tmpCor)];
    end
else
    [GS_seed,clusters_seed]=  clustering (nEle,[x_seed(validIndx) y_seed(validIndx) z_seed(validIndx)],ones(length(x_seed(validIndx)),1));
    % Current Cluster
    x_clustered=x_seed(validIndx);
    y_clustered=y_seed(validIndx);
    z_clustered=z_seed(validIndx);
    clustered=clusters_seed;
end

% update figure
[az, el]=view;
XYZlim=[get(gca,'Xlim');get(gca,'Ylim');get(gca,'Zlim')];
cla;
S.seed=scatter3(x_clustered, y_clustered, z_clustered,36,clr(clustered,:),'Marker','.');
S.GS=scatter3(GS_seed(:,1),GS_seed(:,2),GS_seed(:,3),100,clr(1:nEle,:),'Marker','o','LineWidth',3);
S.all=scatter3(x_all, y_all, z_all,1,'Marker','o');

view(az,el);axis vis3d, box on;rotate3d on;
xlabel('x'),ylabel('y'),zlabel('z');
set(gca,'Xlim',XYZlim(1,:));
set(gca,'Ylim',XYZlim(2,:));
set(gca,'Zlim',XYZlim(3,:));
