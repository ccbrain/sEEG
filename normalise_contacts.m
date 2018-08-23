function Elecs = normalise_contacts(sub, T1, contacts, saveoption)
%% normalise_contacts
%
%Performs segmentation and normalisation of patient T1 image using FSL,
%then applies warping matrix to native electrode contact coordinates to
%bring into standard space
%!!!segmentation/tissue probabilities per coordinate still not saving correctly - needs
%some debugging!!!
%
%Usage: Elecs = normalise_contacts(sub, T1, contactfile)
%where sub is a char array containing patient ID, T1 is the
%co-registered MRI, contact file is the output from contact localisation
%and saveoption determines whether output is saved to .mat file (1 = save,
%0= do not save).
%
%e.g. Elecs = normalise_contacts('002','nameofmyMRI.nii','mycontactfile.mat', 1);
%
%The default filename for T1 is 'rT1F.nii' and for contacts is 'Contacts.mat',
%so these can be left empty if naming convention holds.
%
%e.g. Elecs = normalise_contacts('002',[],[],1);
%
%Variable output is a structure containing electrode names, native voxel coordinates
%for each contact, tissue probabilities for each contact in native space
%and converted MNI coordinates for each contact. The function also
%generates output images using FSL -see line 34 of function for details.
%
%N.B. This is the function wrapper for ContactNormFSL.m

basedir = '/cubric/collab/seeg/analysis/';

cd(basedir);

m = sprintf('Processing subject %s', sub);
disp(m);

%specify subject directory according to 'sub' input
subjdir = sprintf('%s%s/Imaging/',basedir,sub);
cd(subjdir);

%Are all inputs specified or should I use defaults?
if ~exist('T1', 'var') || isempty(T1)
    T1 = fullfile(subjdir,'rT1F.nii');
end

if ~exist('contacts','var') || isempty(contacts)
    contacts = fullfile(subjdir, 'Contacts.mat');
end

%% Correct orientation/origin in T1 image
%these two are required inputs/defaults from above
contactfile = contacts;
regMRI = T1;

%these files will be generated below, if they do not already exist
corrMRI = fullfile(subjdir,'rT1FPA.nii'); %corrects posterior-anterior axis error in coregistered T1
betMRI = fullfile(subjdir, 'rT1FPA_bet.nii.gz'); %product of brain extraction
segMRI = fullfile(subjdir, 'rT1FPA_bet_seg.nii.gz'); %segmented T1 (individual images of tissue compartments are also produced)
MRIwarp = fullfile(subjdir, 'str2standard_warp.nii.gz'); %warping matrix from non-linear normalisation (linear and non-linear normalised images are also produced)

if exist((corrMRI), 'file')
    fprintf('Corrected MRI exists, skipping to next step \n');
elseif exist((regMRI), 'file')
    fprintf('.....Correcting T1 image orientation \n');
    !cp rT1F.nii rT1F_upt.nii
    !fslorient -deleteorient rT1F_upt.nii
    !fslswapdim rT1F_upt.nii -x -y z rT1FPA.nii
    !gzip -d rT1FPA.nii.gz
    newMRI = fullfile(subjdir, 'rT1FPA.nii');
    setorigin_center(newMRI)
else
    error('Warning: file does not exist: rT1F.nii - cannot swap dimensions \n')
end

%% Normalisation using FSL

%%run bet
if exist((betMRI), 'file')
    fprintf('Brain extracted MRI exists, skipping to next step \n');
else
    fprintf('.....Running bet \n');
    system(sprintf('bet %s/rT1FPA %s/rT1FPA_bet -R -m',subjdir, subjdir));
end

%%segmentation - output is seg_0 (CSF), seg_1 (grey matter) and seg_2
%%(white matter)
if exist((segMRI), 'file')
    fprintf('Segmented MRI exists, skipping to next step \n');
else
    fprintf('.....Running segmentation \n');
    system(sprintf('fast -g %s/rT1FPA_bet', subjdir));
    !gunzip -d rT1FPA_bet_pve_0.nii.gz
    !gunzip -d rT1FPA_bet_pve_1.nii.gz
    !gunzip -d rT1FPA_bet_pve_2.nii.gz
end

%%run linear and non-linear normalisation
if exist((MRIwarp), 'file')
    fprintf('Normalisation warping matrix exists, skipping to next step \n');
else
    fprintf('.....Running T1 normalisation \n');
    system(sprintf('flirt -in %srT1FPA_bet -ref /cubric/software/fsl/data/standard/MNI152_T1_1mm_brain -omat %sstr2standard.mat -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -cost mutualinfo',subjdir, subjdir));
    system(sprintf('flirt -in %srT1FPA_bet -ref /cubric/software/fsl/data/standard/MNI152_T1_1mm_brain -out %srT1FPA_2_standard_linear -init  %sstr2standard.mat -applyxfm -interp nearestneighbour', subjdir, subjdir, subjdir));
    system(sprintf('fnirt --in=%srT1FPA --ref=/cubric/software/fsl/data/standard/MNI152_T1_1mm --aff=%sstr2standard.mat --cout=%sstr2standard_warp --interp=spline', subjdir, subjdir, subjdir));
    system(sprintf('applywarp --in=%srT1FPA --ref=/cubric/software/fsl/data/standard/MNI152_T1_1mm_brain --warp=%sstr2standard_warp --out=%srT1FPA_2_standard --interp=spline',subjdir, subjdir, subjdir));
end

%% Load native electrode file, extract and normalise contact points

%%Specify MR files
ntractMask = fullfile([subjdir 'rT1FPA']);
refVo = fullfile([subjdir 'rT1FPA_2_standard']);
warpFile = fullfile([subjdir 'str2standard_warp']);

%%Define tissue compartments
tmp1 = spm_vol(fullfile(subjdir,'rT1FPA_bet_pve_1.nii'));
greymatter = spm_read_vols(tmp1);
tmp2 = spm_vol(fullfile(subjdir,'rT1FPA_bet_pve_2.nii'));
whitematter = spm_read_vols(tmp2);
tmp3 = spm_vol(fullfile(subjdir,'rT1FPA_bet_pve_0.nii'));
csf = spm_read_vols(tmp3);
% tmp1 = spm_vol(fullfile(subjdir,'c1rT1FPA.nii'));
% greymatter = spm_read_vols(tmp1);
% tmp2 = spm_vol(fullfile(subjdir,'c2rT1FPA.nii'));
% whitematter = spm_read_vols(tmp2);
% tmp3 = spm_vol(fullfile(subjdir,'c3rT1FPA.nii'));
% csf = spm_read_vols(tmp3);

clear tmp1 tmp2 tmp3;

%%load electrodes and extract coords
load(contactfile);

%%Load native electrodes into data structure and apply PA flip
TempE = struct;

for i = 1:size(Contacts.contcoor,2)
    TempE(i).contact = Contacts.contcoor{1,i};
end

Elecs = struct;

for j = 1:size(TempE,2)
    m = sprintf('Processing electrode %d of %d',j, size(TempE,2));
    disp(m);
    count=1;
    Elecs(j).name = Contacts.contName(j);
    for k = 1: size(TempE(j).contact,2)
        Elecs(j).contact(k).vox_native = TempE(j).contact(k).coor;
        Elecs(j).contact(k).vox_PAflip_native(1) = 513-(Elecs(j).contact(k).vox_native(1));
        Elecs(j).contact(k).vox_PAflip_native(2) = 513-(Elecs(j).contact(k).vox_native(2));
        Elecs(j).contact(k).vox_PAflip_native(3) = Elecs(j).contact(k).vox_native(3); %set native voxels and PA flipped coordinates
        %%assign tissue probabilities to native contacts
        currcont = Elecs(j).contact(k).vox_PAflip_native;
        Elecs(j).tissues(k).greymatter = greymatter(currcont);
        Elecs(j).tissues(k).whitematter = whitematter(currcont);
        Elecs(j).tissues(k).csf = csf(currcont);
        Conlist(k,:,:,:) = Elecs(j).contact(k).vox_PAflip_native;
        %%Now write a text file containing elec contacts to feed to img2imgcoord
        fID = fopen(sprintf('voxcoords_E%d.txt',j),'w');
        fprintf(fID,'%5d %5d %5d\n',Conlist);
        fclose(fID);
        count = count+1;
    end
    
    %%get normalised contact coordinate location
    [~,g] = system(sprintf('img2imgcoord -src %s -dest %s -warp %s -vox voxcoords_E%d.txt',ntractMask,refVo,warpFile,j));
    NormCoords = str2num(g(46:end));
    Elecs(j).normcoords = round(NormCoords);
    
end

fprintf('Cleaning up temporary files \n');
!rm voxcoords_E*.txt

if saveoption == 1
    fprintf('Finished, saving \n');
    save(fullfile(subjdir,'Elecs_MNI.mat'),'Elecs');
elseif saveoption == 0
    fprintf('Finished, output not saved \n');
end

end
