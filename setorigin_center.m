function setorigin_center(files)
% set origin of image files to the center of xyz dimensions using spm
% functions
% Fumio Yamashita 2014.1.20
%% check arguments
if nargin == 0
    files = spm_select(Inf,'image','Select image files');
end
%% main loop
for i=1:size(files,1)
    file = deblank(files(i,:));
    st.vol = spm_vol(file);
    vs = st.vol.mat\eye(4);
    vs(1:3,4) = (st.vol.dim+1)/2;
    spm_get_space(st.vol.fname,inv(vs));
end