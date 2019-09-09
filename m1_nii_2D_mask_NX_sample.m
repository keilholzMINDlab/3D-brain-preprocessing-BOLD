clear;
close all
clc
%% define path
addpath('C:\Users\nxu25\Dropbox (GaTech)\Mouse Data from Kai-hsiang’s Group\Code\Preprocessing_Xiaodi\tools');
path_read='../../Data/';
path_write_4D='..\..\BOLD_data\step1_4D_nii(rotate and mask)\';
path_mask='..\..\BOLD_data\Mask\';
path_quality='..\..\BOLD_data\Quality Assurance\Masking\';
mkdir(path_quality);

File=load([path_read 'filename.mat']);
% Temp=load('filename.mat','BOLD_filename_array','date_array','scan_number_array');
% filename_array=Temp.BOLD_filename_array;
% date_array=Temp.date_array;
% scan_number_array=Temp.scan_number_array;

N=size(File.filename_array,1);
for time_id=1
    filename=char(File.filename_array(time_id));
    nii = load_untouch_nii([path_read filename '.nii']);
    [x,y,z,t]=size(nii.img);
    images_rot_masked=zeros(y,x,z,t);
    mask_final=zeros(y,x,z);
    for iz=1:z   
        I=zeros(y,x,t);
        for it=1:t
            I_cur=squeeze(nii.img(:, :, iz, it));
            I(:,:,it)=rot90(I_cur,3);            
        end
        [mask,maskedImage] = segmentImage(I(:,:,1));    
        
        [cen_x,cen_y]=find(mask);
        cen_x=(max(cen_x)-min(cen_x))/2+min(cen_x);
        cen_y=(max(cen_y)-min(cen_y))/2+min(cen_y);  
        
        mask1=circshift(mask,-(round(cen_y-x/2)),2);
        mask1=circshift(mask1,-(floor(cen_x-y/2)),1);
        maskedImage1=circshift(maskedImage,-(round(cen_y-x/2)),2);
        maskedImage1=circshift(maskedImage1,-(floor(cen_x-y/2)),1);

   
        % fine segment        
        II=circshift(I,-(round(cen_y-x/2)),2);
        II=circshift(II,-(floor(cen_x-y/2)),1);
        [mask2,maskedImage2] = segmentImage2(II(:,:,1),mask1); %% mask3 is mask2 eroded by 2 pixels
        
        %% visualization
        h1=figure;clf;
        set(h1,'units','normalized','outerposition',[0 0 1 1]);
        tp_sample=randi(z);
        imagesc([II(:,:,tp_sample),  maskedImage1, II(:,:,tp_sample).*mask2]);
        axis equal
       
        %% quality assurance
        saveas(h1,[path_quality '\' filename 'zplane' num2str(iz) '_tps' num2str(tp_sample) '.png']);
        
        %% masking
        maskedImages_z=II.*repmat(mask2,[1,1,t]);
        images_rot_masked(:,:,iz,:)=permute(maskedImages_z,[1 2 4 3]);
%         images_rot_masked(:,:,iz,:)=squeeze(maskedImages_z);        

        mask_final(:,:,iz)=mask2;

    end 
        nii_rot_mask=make_nii(images_rot_masked);
        nii_rot_mask.untouch=0;

        %% save
%         save_nii(nii2,[path_write_3D '3D_' filename '.nii']);
%         save_nii(nii3,[path_write_4D '4D_' filename '.nii']);
%         save([path_mask 'mask_' filename '.mat'],'mask');
        time_id

end