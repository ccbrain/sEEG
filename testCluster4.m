% 1.Change from ANALYZE to nifti format for both CT2 and T1 in
% matlab using SPM. I do not know if FSL reslices the imaging during
% transformation so it is safter to do it in SPM
%  	
%  	f='T1'  % should change accordingly;
%  	V=spm_vol([f '.img']);
%  	ima=spm_read_vols(V);
%  	ima=ima+32768;       % the T1 has strange negatie values here, this will make it impossible to be processed. DO NOT need to to this for CT!
%  	V.fname=[f '.nii'];
%   spm_write_vol(V,ima);
% 
%   f='CT-2'  % should change accordingly;
%  	V=spm_vol([f '.img']);
%  	ima=spm_read_vols(V);
%  	ima(ima<600)=0;       % Get rid of scanner bed, wires, etc.
%  	V.fname=[f '.nii'];
%   spm_write_vol(V,ima);

% 2. Remove the .hdr and .img files to avoid confusing FSL

% 3. Coreg CT1.nii to CT2.nii in SPM, gives rCT1.nii
%       3.1 Take the difference map between CT2 and CT1 (CT2-CT1), can do it in
%           SPM imcalc, or in Matlab

% 4. Coreg T1.nii to CT2.nii in SPM, gives rT1.nii
%       4.1 Segmentation rT1.nii for future use
%       4.2. Brain extration of T1
%           bet T1.nii T1_brain.nii -R -m
%           fslmaths T1_brain_mask -ero T1_brain_mask_ero
%           gzip -d T1_brain_mask_ero.nii.gz   # FSL give compressed nii
%           file. Unzip it for SPM
%       4.3 apply the brain mask to CT2-CT1 nii file
% 5. Solution Below

% cd /cubric/scratch/sapjz/sEEG002;

clear all
close all

addpath(genpath('~/tools/iElectrodes'));
addpath(genpath('/cubric/software/spm'));
addpath(genpath('/cubric/collab/seeg/scripts'));

basedir = '/cubric/collab/seeg/analysis/055/Imaging/'; %change per patient
nEle = 013;  % num of electrodes from PatientInfo.xlsx - change per patient

cd(basedir);

firstThreshold = 650;  % the first-pass intensity threshold in CT (high enough to get rid of brain tissue residuals, low enough to avoid masking electodes)
CT_masked = fullfile(basedir,'brain_masked_CTdiff.nii'); % the CT nii file, take the difference between CT2 and CT1
brain_mask = fullfile(basedir,'brain_rT1F_mask_ero.nii'); % the binary brain mask from co-reged T1. "bet T1C.nii T1C_brain -R -m", and then erode the mask with default option "fslmaths rT1C_brain_mask -ero rT1C_brain_mask_ero"

%% imaging CT data (masked within brain)
V=spm_vol(CT_masked);
ima=spm_read_vols(V);

%% Brain mask
V_mask=spm_vol(brain_mask);
ima_mask=spm_read_vols(V_mask);

%% erode the brain mask to get the seed mask close to cortical surface
[xx,yy,zz] = ndgrid(-5:5);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 4.0;
ima_mask_erode = imerode(ima_mask,nhood);

[xx,yy,zz] = ndgrid(-5:5);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 4.0;
ima_mask_erode2 = imerode(ima_mask_erode,nhood);
ima_mask_seed=ima_mask_erode-ima_mask_erode2;  % 2nd erode

%% binarize CT data with a given threshold
ima_bw=(ima>firstThreshold);
ind_seed = find (ima_mask_seed & ima_bw);
[x_seed, y_seed, z_seed] = ind2sub(size(ima), ind_seed);
[x, y, z] = ind2sub(size(ima), find(ima_mask>0));

[x_all, y_all, z_all] = ind2sub(size(ima), find(ima_bw & ima_mask_erode));



% Keep only the electode handles
validIndx = find(x_seed<(min(x)+(max(x)-min(x))*0.25) | x_seed>(max(x)-(max(x)-min(x))*0.25) | z_seed>(min(z)+(max(z)-min(z))*0.5));

% save some memory
clear ima_mask_erode ima_mask_seed;

if ~isempty(dir(fullfile(basedir,'initCL.mat')))
    load(fullfile(basedir,'initCL.mat'));
else
    %% check initial seed clusters, correct if problem occurs here.
    % [clusters,GS] = kmeans([x_seed x_seed x_seed],nEle,'Start','sample','Replicates',500,'Display','Iter','Distance','cityblock');
    [GS_seed,clusters_seed]=  clustering (nEle,[x_seed(validIndx) y_seed(validIndx) z_seed(validIndx)],ones(length(x_seed(validIndx)),1));  % use the func from iElectrodes toolbox to give fixed solution
    % Current Cluster
    x_clustered=x_seed(validIndx);
    y_clustered=y_seed(validIndx);
    z_clustered=z_seed(validIndx);
    clustered=clusters_seed;
    
    % Quality check, do we have all electrodes in separate clusteres?
    %clr = lines(nEle);
    f=figure; hold on;
    set(f,'Position',[800 200 1100 800]);
    clr = lines(nEle);
    S.seed=scatter3(x_clustered, y_clustered, z_clustered,36,clr(clustered,:),'Marker','.');
    S.GS=scatter3(GS_seed(:,1),GS_seed(:,2),GS_seed(:,3),100,clr(1:nEle,:),'Marker','o','LineWidth',3);
    S.all=scatter3(x_all, y_all, z_all,1,'Marker','o');
    view(3);axis vis3d, box on;rotate3d on;
    xlabel('x'),ylabel('y'),zlabel('z');
    title('If needed, get CLUSTER and NOCLUSTER variables using the brush tool');
        p = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.85 .05 .1 .08],...
        'String','Refresh',...
        'CallBack','uptFigure');
    p2 = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.05 .85 .1 .08],...
        'String','Vox on/off',...
        'CallBack','uptFigure2');
        p3 = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.05 .75 .1 .08],...
        'String','Seed on/off',...
        'CallBack','uptFigure3');
    p4 = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.05 .65 .1 .08],...
        'String','GS on/off',...
        'CallBack','uptFigure4');
    q = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.05 .05 .1 .08],...
        'String','Close',...
        'CallBack','close(gcbf)');
    q2 = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.85 .85 .1 .08],...
        'String','Restore',...
        'CallBack','clearCLUSTERs');
    waitfor(f);

    % save the initial clustering results
    save(fullfile(basedir,'initCL.mat'),'clusters_seed','x_clustered','y_clustered','z_clustered','GS_seed','clustered');
end


%% for each seed cluster, set the seed in results
ima_clustered=zeros(size(ima));  % image matrix for clustered voxels
% set clusters from seeds
for i=1:nEle
    ima_clustered(sub2ind(size(ima_clustered),x_clustered(clustered==i),y_clustered(clustered==i),z_clustered(clustered==i))) = i;
end
% the rest voxels to be clustered within Erode2 Masks
ind_rest = find (ima_mask_erode2 & ima_bw);
[x, y, z] = ind2sub(size(ima), ind_rest);
newCor = setdiff([x y z],[x_clustered y_clustered z_clustered],'rows');
x=newCor(:,1);
y=newCor(:,2);
z=newCor(:,3);

% fit SVD lines per electrode
for i=1:nEle
    % fit a line
    D=[x_clustered(clustered==i) y_clustered(clustered==i) z_clustered(clustered==i)];
    [p0(i,:),d(:,i)]=svdfit(D);
    maxDist(i) = 20;
end

% % % plot the current SVD line
% clr = lines(nEle);
% FF=figure; hold on;
% scatter3(x_clustered, y_clustered, z_clustered,36,clr(clustered,:),'Marker','.');
% % scatter3(GS(:,1),GS(:,2),GS(:,3),100,clr,'Marker','o','LineWidth',3);
% view(3);axis vis3d, box on;rotate3d on;
% xlabel('x'),ylabel('y'),zlabel('z');
% t=-100:100;
% for k=1:nEle
%     for i=1:length(t)
%         P(i,:) = p0(k,:) + d(:,k)'.*t(i);
%     end
%     scatter3(P(:,1),P(:,2),P(:,3),'.');
%     clear P;
% end

%% IMPORTANT STEP:Plot brain surface and electode handles
% Based on this plot, we can assign electrode letters to clusteres
% Put result in a text file "EleMap" in a format (%ElectodeName %contactNumber %clusterID)
% P 14 1
% Y' 16 2
% Y 16 3
% P' 12 4
% X 16 5
% O' 12 6
% ...
clr = lines(nEle);
h=plotBrainandElec(ima_mask_erode2,ima_clustered,x_clustered,y_clustered,z_clustered,x,y,z,clustered,GS_seed,clr,0);

%% initial tracking from seeds
allSVDdist=Inf;
while (1)
    clear allDist;
    for i=1:nEle
        allDist(:,i)=svdfitDist(p0(i,:),d(:,i),[x_all y_all z_all]);
    end
    [~, tempCluster]=min(allDist,[],2);
    allSVDdist_cur=max(allDist(sub2ind(size(allDist),1:size(allDist,1),tempCluster')));
    
    % stop tracking if the SVD distance is stable (<2 voxel change)
    if (allSVDdist - allSVDdist_cur ) < 2
        break;
    % Otherwise, continue
    else
        allSVDdist=allSVDdist_cur;
        disp(['Current max of SVD distance: ' num2str(allSVDdist)]);
        clear allDist;
        steps = 100;
        k=0;
        while ~isempty(x) && k<20
            disp([num2str(length(x)) ' voxels remaining ...']);
            stepClusters=min(length(x),steps);
            clear allDist allIndx allDistfromGS;
            for i=1:nEle
                % Distance from GM of seeds
                allDistfromGS(:,i)=sqrt(sum(([x y z]-repmat(GS_seed(i,:),[size(x),1])).^2,2));
                % Refit SVD lines
                D=[x_clustered(clustered==i) y_clustered(clustered==i) z_clustered(clustered==i)];
                [p0(i,:),d(:,i)]=svdfit(D);
                allDist(:,i)=svdfitDist(p0(i,:),d(:,i),[x y z]);
            end
            
            % find voxels closest to seeds
            CostFunc=allDistfromGS;
            [tempVal,tempIndx]=sort(CostFunc(:));
            [tempVox, tempCluster]=ind2sub(size(CostFunc), tempIndx(1:stepClusters));
            tempDist = allDist(sub2ind(size(allDist),tempVox,tempCluster));
            
            strangIndx=tempDist>maxDist(tempCluster)';
            
            tempVox(strangIndx)=[];
            tempCluster(strangIndx)=[];
            
            if length(unique(tempVox))~=stepClusters
                warning('overlapping clusters');
            end
            
            for i=1:length(tempVox)
                ima_clustered(sub2ind(size(ima_clustered),x(tempVox(i)),y(tempVox(i)),z(tempVox(i)))) = tempCluster(i);
            end
            clustered = [clustered; tempCluster];
            
            x_clustered = [x_clustered;x(tempVox)];
            y_clustered = [y_clustered;y(tempVox)];
            z_clustered = [z_clustered;z(tempVox)];
            
            x(tempVox)=[];
            y(tempVox)=[];
            z(tempVox)=[];
            
            x(strangIndx)=[];
            y(strangIndx)=[];
            z(strangIndx)=[];
            k=k+1;
        end
    end
end

%% Further tracking for the rest of the voxels
% fit the remaining voxels along the SVD line
clear allDist;
for i=1:nEle
    allDist(:,i)=svdfitDist(p0(i,:),d(:,i),[x y z]);
end
[~, tempCluster]=min(allDist,[],2);
for i=1:length(x)
    ima_clustered(sub2ind(size(ima_clustered),x,y,z)) = tempCluster;
end
clustered = [clustered; tempCluster];
x_clustered=[x_clustered;x];
y_clustered=[y_clustered;y];
z_clustered=[z_clustered;z];

% Refit the SVD line
for i=1:nEle
    D=[x_clustered(clustered==i) y_clustered(clustered==i) z_clustered(clustered==i)];
    [p0(i,:),d(:,i)]=svdfit(D);
end

%% Recluster based on the final SVD line
clear allDist;
clear allDistfromMean;
clear allDistfromGS;
clear x y z;
%% Erode the brain mask to get the seed mask close to cortical surface
[xx,yy,zz] = ndgrid(-5:5);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 3.0;
ima_mask_erode3 = imerode(ima_mask,nhood);
[x, y, z] = ind2sub(size(ima), find(ima_bw & ima_mask_erode3));

for i=1:nEle
    % Distance to each SVD line
    allDist(:,i)=svdfitDist(p0(i,:),d(:,i),[x y z]);
    
    % Euclidean distance to the center of the line
    allDistfromMean(:,i)=sqrt(sum((repmat(p0(i,:),length(x),1)-[x y z]).^2,2));
    allDistfromGS(:,i)=sqrt(sum(([x y z]-repmat(GS_seed(i,:),[size(x),1])).^2,2));
end

% Cost function = combined euclidean distance to seeds, SVD distance and SVD value along the line
[~, tempCluster]=min(allDist.*20+allDistfromMean+allDistfromGS,[],2);

clustered=tempCluster;
x_clustered=x;
y_clustered=y;
z_clustered=z;

% Refit the SVD line
for i=1:nEle
    D=[x_clustered(tempCluster==i) y_clustered(tempCluster==i) z_clustered(tempCluster==i)];
    [p0(i,:),d(:,i)]=svdfit(D);
end

ima_clustered=zeros(size(ima));
for i=1:length(x_clustered)
    ima_clustered(sub2ind(size(ima_clustered),x_clustered,y_clustered,z_clustered)) = tempCluster;
end

% recheck if all electrodes are clustered correctly
h=plotBrainandElec(ima_mask_erode2,ima_clustered,x,y,z,[],[],[],clustered,GS_seed,clr,1);
t=-100:100;
for k=1:nEle
    for i=1:length(t)
        P(i,:) = p0(k,:) + d(:,k)'.*t(i);
    end
    scatter3(h.ax2,P(:,1),P(:,2),P(:,3),5,repmat(clr(k,:),size(P,1),1),'.');
    scatter3(h.ax1,P(:,1),P(:,2),P(:,3),5,repmat(clr(k,:),size(P,1),1),'.');
    clear P;
end


%% Localise each contact point
% Read contact file in a format (%ElectodeName %contactNumber %clusterID)
% Can also pass the 4th and 5th columns, which limit the range [4th column, 5th column] of the SVD
% line, to use for exclusively masking heads and tails.
fid = fopen(fullfile(basedir,'EleMap.txt'));
FC = textscan(fid, '%s%f%f%f%f%f');
fclose(fid);
Contacts.contName=FC{1};
Contacts.contNum=FC{2};
Contacts.cluster=FC{3};
Contacts.range = [FC{4} FC{5}];
Contacts.p0=p0;
Contacts.d=d;

% Results for a giving cluster
for i=1:nEle
    pEle=i;
    figure; hold on;
    
    allP=[x_clustered(clustered==pEle),y_clustered(clustered==pEle),z_clustered(clustered==pEle)];
    % only consider voxels with >700 intensity
    P=allP(ima(sub2ind(size(ima),allP(:,1),allP(:,2),allP(:,3)))>700,:);
    
    bbb=ima(sub2ind(size(ima),P(:,1),P(:,2),P(:,3)));
    linesp=svd3D_1D(p0(pEle,:), d(:,pEle)',P);
    linedist=svdfitDist(p0(pEle,:),d(:,pEle),P);
    
    scatter(linesp(linedist<20),bbb(linedist<20),1,repmat(linedist(linedist<20)./max(linedist),1,3));
    
    % clustering to get individual contact
    if sum(isnan(Contacts.range(pEle,:)))==2
        indx=linedist<2; 
    elseif sum(isnan(Contacts.range(pEle,:)))==0
        indx=linedist<2 & linesp<=max(Contacts.range(pEle,:)) & linesp>=min(Contacts.range(pEle,:));
    else
        error(['The ranage of electrode ' Contacts.contName{pEle} ' is wrong!']);
    end
    xx=sub2ind(size(ima),P(indx,1),P(indx,2),P(indx,3));
    [GS2,clusters2]=  clustering (Contacts.contNum(pEle),linesp(indx),ima(xx));
    
    linespGS = svd3D_1D(p0(pEle,:), d(:,pEle)',GS_seed(pEle,:)) ;
    plot([linespGS linespGS],[0 2000]);
    scatter(GS2(:,1),zeros(Contacts.contNum(pEle),1)+2000);
        
    % coordinate for each contact, the first index always the one
    % furthest to the handle
    GScontacts = sort(GS2(:,1),'ascend');
    if dist(linespGS,GScontacts(1))<dist(linespGS,GScontacts(end))
        GScontacts = sort(GS2(:,1),'descend');
    end
    for j=1:Contacts.contNum(pEle)
        Contacts.contcoor{pEle}(j).GS=GScontacts(j);
        tempcoor=round(svd1D_3d(p0(pEle,:), d(:,pEle)',GScontacts(j)));
        % find the max intensity and its coordinates in a 3x3x3 cube
        bnd = 1;
        [x y z] = meshgrid((tempcoor(1)-bnd):(tempcoor(1)+bnd),(tempcoor(2)-bnd):(tempcoor(2)+bnd),(tempcoor(3)-bnd):(tempcoor(3)+bnd));
        [~, b]=max(ima(sub2ind(size(ima),x(:), y(:) ,z(:))));
        Contacts.contcoor{pEle}(j).coor=[x(b), y(b) ,z(b)];
    end
end
save(fullfile(basedir,'Contacts'),'Contacts');


% project 1D point to 3D coordinate
%  round(svd1D_3d(p0(13,:), d(:,13)',64.3))