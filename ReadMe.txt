This dataset is made up of data from 100+ patients processed through an epilepsy centre in Guangzhou, China. 

For each patient (where dataset is complete), we have pre-op and post-op CT scans, a T1 MRI, sEEG recordings saved as .edf files and associated event markers. 

At the top level of the repository, see spreadsheet 'PatientInfo' for data, processing and demographic info/notes. The spreadsheet 'chinese_number_conversion.xlsx' give the correspondence between the original chinese character format and the assigned numeric IDs. Both files give indications about incomplete datasets. 

The archive is organised as follows:
1. Scripts - containing, unsurprisingly (!), scripts
2. Data - a complete record of the raw data for each individual along with surgical planning and clinical info
3. Analysis - again organised by individual, then containing separate folders for preprocessing of imaging and electrophys data
4. Outputs - containing drafts of manuscripts and conference abstracts. 

Image pre-processing steps have included: Co-registration of CT1 (pre-op) and MRI to CT2 (post-op), semi-automated localisation of electrode contacts, segmentation of T1 and identification of tissue compartment for each contact, normalisation of both T1 image and electrode contact coordinates. As of April 2018, the processing pipeline has used the following scripts: bulkrename.sh, CT_MR_preprocess.m, cluster_electrodes.m, normalise_contacts.m. All scripts are commented but see also precis below:  

1. bulkrename.sh -- should not need to be run again but converts data folders from original chinese names/characters to numeric IDs (see above)

2.  CT_MR_preprocess.m -- again, should not need to be run again for this dataset. This script completes all necessary preprocessing of CT and MRI images, including conversion to nifti, coregistration of pre-op CT and MRI to post-op CT, skull stripping of MRI and creation of brain-masked CT difference image for use in electrode clustering. 

3. cluster_electrodes.m -- this script performs clustering of sEEG electrodes. Starting from outer edges of mask for initial clustering and tracking inwards using k-means clustering and singular value decomposition line fitting, each voxel in the masked CT image is allocated automatically to an electrode in the array. Manual interaction via a GUI allows correction at the initial clustering stage. 

4. normalise_contacts.m -- this function performs segmentation and normalisation of the T1 image for a patient using FSL, and applied the resultant warping matrix to the previously generated electrode contact coordinates for that patient. The output is a .mat file containing a data structure which holds electrode names, native voxel coordinates, tissue probabilities for each contact in native space and MNI converted coordinates for each contact. 

sEEG data have not yet been preprocessed (May 2018). 
