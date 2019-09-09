clear;

%% define path
path_read='..\..\BOLD_data\step1_4D_nii(rotate and mask)\3Dmasked\';
path_write='..\..\BOLD_data\step2_motion_corrected_3Dmasked_nii\';
path_write_translateMat=[path_write 'translateMat\'];
path_write_trajectory=[path_write 'trajectory\'];
Temp=load('..\..\Data\filename.mat','BOLD_filename_array');
filename_array=Temp.BOLD_filename_array;
%% define SPM job and run motion correction
t=2000;
matlabbatch=createSPM_reg2img1_job(path_read,filename_array(:,1:3)',t);
spm('defaults', 'FMRI');
spm;
spm_jobman('run', matlabbatch);

%% move files
N=size(matlabbatch,2);
status=zeros(3,N);
mkdir(path_other);
mkdir(path_write);
for i=1:3
    status(1,i)=movefile([path_read '3Dmasked_' char(filename_array(:,i)) '.mat'],[path_other 'correction_matrix\']);
    status(2,i)=movefile([path_read 'rp_3Dmasked_' char(filename_array(:,i)) '.txt'],[path_other 'trajectory\']);
    status(3,i)=movefile([path_read 'r3Dmasked_' char(filename_array(:,i)) '.nii'],path_write);
end
status

