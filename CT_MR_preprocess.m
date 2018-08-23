%% Preprocessing of CT and MRI images in preparation for automatic sEEG electrode detection.
% Script assumes that scans are named T1F, CTpre and CTpost, are in ANALYZE
% format, and are organised by subject (according to 'ImagingList.txt' e.g. 001, 002,) in
% subdirectory /path/<subjectno>/Imaging/
% Preprocessing steps include: Conversion from ANALYZE to nifti, co-registration of
% pre-surgical CT and T1 MRI to post-surgical CT, skull-stripping of MRI
% and creation of brain-masked CT difference image (CT2-CT1).
% Written by Beth Routley and Jiaxiang Zhang- November 2017

clear all
close all

basedir = '/cubric/collab/seeg/analysis/';

cd(basedir);

temp = regexp(fileread('ImagingList.txt'), '\n', 'split'); %text file containing list of subject numbers for whom we have complete imaging
subs = vertcat(temp{:});
clear temp;

for i=1:1%length(subs)
    %% Navigate to subject folder
    
    subjdir = sprintf('%s%s/Imaging/',basedir,subs(i,:));
    cd(subjdir)
    
    m = sprintf('Processing subject %d of %d',i, length(subs));
    disp(m);
    
%     %% Convert CT1, CT2 and T1 from ANALYZE to nifti
%     mr = 'T1F';
%     V = spm_vol([mr '.img']);
%     ima = spm_read_vols(V);
%     if min(ima) < -3000; %Strange intensely negative values on background of T1 image, reset if necessary
%         ima = ima+32768;
%     end
%     V.fname = [mr '.nii'];
%     spm_write_vol(V,ima);
%     
%     ct1 = 'CTpre';
%     V = spm_vol([ct1 '.img']);
%     ima = spm_read_vols(V);
%     ima(ima<350)=0;       % Threshold image to get rid of scanner bed, wires, etc.
%     V.fname = [ct1 '.nii'];
%     spm_write_vol(V,ima);
%     
%     ct2 = 'CTpost';
%     V = spm_vol([ct2 '.img']);
%     ima = spm_read_vols(V);
%     ima(ima<350)=0;
%     V.fname = [ct2 '.nii'];
%     spm_write_vol(V,ima);
    
%     %% Remove the .hdr and .img files to avoid confusing FSL
%     if exist((sprintf('%sT1F.nii', subjdir)), 'file')
%         !rm T1F.hdr
%         !rm T1F.img
%     else
%         error('Warning: file does not exist: T1F.nii - nifti conversion unsuccessful')
%     end
%     
%     if exist((sprintf('%sCTpre.nii', subjdir)), 'file')
%         !rm CTpre.hdr
%         !rm CTpre.img
%     else
%         error('Warning: file does not exist: CTpre.nii - nifti conversion unsuccessful')
%     end
%     
%     if exist((sprintf('%sCTpost.nii', subjdir)), 'file')
%         !rm CTpost.hdr
%         !rm CTpost.img
%     else
%         error('Warning: file does not exist: CTpost.nii - nifti conversion unsuccessful')
%     end
%     
    
    %% Coregister CT1 to CT2
    fileCT1 = [subjdir 'CTpre.nii'];
    fileCT2 = [subjdir 'CTpost.nii'];
    VolCT2 = spm_vol(fileCT2); % reference image
    VolCT1 = spm_vol(fileCT1); % source (moved) image
    
    %estimate options (all spm defaults)
    eoptions.cost_fun = 'nmi'; %normalized mutual information
    eoptions.sep = [4 2];
    eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    eoptions.fwhm = [7 7];
    
    % estimate and save warping matrix
    matCTpre2CTpost = spm_coreg(VolCT2,VolCT1,eoptions);
    MatCTpre2CTpost = spm_matrix(matCTpre2CTpost);
    save ([subjdir '/MatCTpre2CTpost.mat'],'MatCTpre2CTpost');
    
    % apply transformation to CT1
    Mi  = inv(MatCTpre2CTpost);
    MM = VolCT1.mat;
    spm_get_space(fileCT1,Mi*MM); % change(set), and write to the input file
    
    %set vars and options for reslicing
    P = str2mat(fileCT2, fileCT1); % filenames for reslicing
    roptions.interp = 4;
    roptions.wrap = [0 0 0];
    roptions.mask = 0;
    roptions.mean = 0;
    roptions.which = 1;
    roptions.prefix = 'r';
    
    %reslice CT1
    spm_reslice(P,roptions);
    
    %     %check coreg worked - N.B. has graphics so do not use if running full loop/on cluster
    filerCT1 = [subjdir '/rCTpre.nii'];
    %     p = strvcat(fileCT2, filerCT1);
    %     spm_check_registration(p,{'CTpost ','rCTpre'});
    
    %% Create CT difference image
    Vi = strvcat(fileCT2, filerCT1);
    Vo = [subjdir 'CTdiff.nii'];
    spm_imcalc(Vi, Vo, 'i1-i2');
    
    %% Coregister T1 to CT2
    fileT1F = [subjdir 'T1F.nii'];
    fileCT2 = [subjdir 'CTpost.nii'];
    VolCT2 = spm_vol(fileCT2); % reference image
    VolT1F = spm_vol(fileT1F); % source (moved) image
    
    %estimate options (all spm defaults)
    eoptions.cost_fun = 'nmi'; %normalized mutual information
    eoptions.sep = [4 2];
    eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    eoptions.fwhm = [7 7];
    
    % estimate and save warping matrix
    matT1F2CTpost = spm_coreg(VolCT2,VolT1F,eoptions);
    MatT1F2CTpost   = spm_matrix(matT1F2CTpost);
    save ([subjdir '/MatT1F2CTpost.mat'],'MatT1F2CTpost');
    
    % apply transformation to T1
    Mi  = inv(MatT1F2CTpost);
    MM = VolT1F.mat;
    spm_get_space(fileT1F,Mi*MM); % change(set), and write to the input file
    
    %set vars and options for reslicing
    P = str2mat(fileCT2, fileT1F); % filenames for reslicing
    roptions.interp = 4;
    roptions.wrap = [0 0 0];
    roptions.mask = 0;
    roptions.mean = 0;
    roptions.which = 1;
    roptions.prefix = 'r';
    
    %reslice T1F
    spm_reslice(P,roptions);
    
%     %check coreg worked - N.B. has graphics so do not use if running full loop/on cluster
    filerT1F = [subjdir 'rT1F.nii'];
%     p = strvcat(fileCT2, filerT1F);
%     spm_check_registration(p,{'CTpost ','T1F'});
    
    %% Extract brain mask from coregistered T1 image
    
    fileBrainT1 = [subjdir 'brain_rT1F.nii'];
    eval (['!bet ' filerT1F ' ' fileBrainT1 ' -R -m']);
    %     eval (['!gunzip -d ' fileBrainT1]);
    %     eval (['!gunzip -d brain_rT1F_mask.nii']);
    fileBrainT1Mask = [subjdir 'brain_rT1F_mask'];
    fileT1MaskEro = [subjdir 'brain_rT1F_mask_ero.nii.gz'];
    eval (['!fslmaths ' fileBrainT1Mask ' -ero ' fileT1MaskEro]);
    eval (['!gunzip -d ' fileT1MaskEro]);
    fileT1MaskEro = [subjdir 'brain_rT1F_mask_ero.nii'];
    
    %% Apply the brain mask to CT difference image
    
    fileCTdiff = [subjdir 'CTdiff.nii'];
    fileMaskedCT = [subjdir 'brain_masked_CTdiff.nii.gz'];
    eval (['!fslmaths ' fileCTdiff ' -mas ' fileT1MaskEro ' ' fileMaskedCT]);
    eval (['!gunzip -d ' fileMaskedCT]);
    
    m = sprintf('Finished subject %d of %d',i, length(subs));
    disp(m);
    
    clearvars -except basedir subs
end

%Follow this with electrode localisation
