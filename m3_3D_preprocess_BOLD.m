% clear;

%% define path
path_raw='..\..\BOLD_data\step1_4D_nii(rotate and mask)\3Dmasked\';
path_read='..\..\BOLD_data\step2_motion_corrected_3Dmasked_nii\';
path_write_nii='..\..\BOLD_data\step3_motion_corrected_3D_nii\';
path_write_mat='..\..\BOLD_data\processed\';
path_mask='..\..\BOLD_data\mask\';
Temp=load('..\..\Data\filename.mat','BOLD_filename_array');
mkdir(path_write_nii);
mkdir(path_write_mat);
filename_array=Temp.BOLD_filename_array;
% LFP_filename_array=Temp.LFP_filename_array;
N=size(filename_array,1);
status=zeros(1,N);

for rats_id=1:N
    filename=char(filename_array(rats_id,:));
%     LFP_filename=char(LFP_filename_array(rats_id,:));
    fprintf(['start process ' filename '... ' num2str(rats_id) '\n']);
    if exist([path_read 'r4D_' filename '.nii'],'file')==2
        %% loading
        nii = load_untouch_nii([path_read 'r4D_' filename '.nii']);
        BOLD_corrected=nii.img;        
        BOLD_corrected(isnan(BOLD_corrected))=0;        
        [x, y, z, t]=size(BOLD_corrected);
        
        nii_raw = load_untouch_nii([path_raw '3Dmasked_' filename '.nii']);
        BOLD_raw=nii_raw.img;

        for iz=1:z
            BOLD_corrected(:,:,iz,:)=BOLD_corrected(:,:,iz,:)/max(max(max(BOLD_corrected(:,:,iz,:))));
            BOLD_raw(:,:,iz,:)=BOLD_raw(:,:,iz,:)/max(max(max(BOLD_raw(:,:,iz,:))));
        end


        %% masking
        TempMask=load([path_mask 'mask_' filename '.mat'],'mask3D_final');
        mask=TempMask.mask3D_final;
        
        radius=1;
        decomposition=0;
        se=strel('disk', radius, decomposition);
        mask=imerode(mask, se);
        BOLD_corrected=BOLD_corrected.*repmat(mask,[1,1,1,t]);
        
        %% recenter
%         [x,y]=meshgrid(-32:31,-32:31);
%         x_cen1=mean(x(find(mask)));
%         y_cen1=mean(y(find(mask)));
%         BOLD_raw=circshift(BOLD_raw,-round(y_cen1),1);
%         BOLD_raw=circshift(BOLD_raw,-round(x_cen1),2);
%         BOLD_corrected=circshift(BOLD_corrected,-round(y_cen1),1);
%         BOLD_corrected=circshift(BOLD_corrected,-round(x_cen1),2);
%         mask=circshift(mask,-round(y_cen1),1);
%         mask=circshift(mask,-round(x_cen1),2);
        
        %% save nii
        nii_corrected2=make_nii(BOLD_corrected);
        nii_corrected2.untouch=0;
        save_nii(nii_corrected2,[path_write_nii 'r3D_' filename '.nii']);
        
        %% clip
%         BOLD_raw=BOLD_raw(9:56,9:56,:);
%         BOLD_corrected=BOLD_corrected(9:56,9:56,:);
%         mask=mask(9:56,9:56,:);
        
        %% spatial smoothing
        sigma=1;        
        mu = [0 0];
        Sigma = [sigma 0; 0 sigma];
        d=3;
        BOLD_smoothed=zeros(x,y,z,t);
        BOLD_smoothed_matlab=zeros(x,y,z,t);
        for it=1:t
%             temp=conv2(BOLD_corrected(:,:,i),gauss_win);
%             BOLD_smoothed(:,:,i)=temp(d+1:end-d,d+1:end-d);
%             BOLD_smoothed_matlab(:,:,it)=imgaussfilt(BOLD_corrected(:,:,i),sigma,'FilterDomain','spatial','FilterSize',d);
            BOLD_smoothed_matlab(:,:,:,it)=imgaussfilt3(BOLD_corrected(:,:,:,it),[sigma,sigma,sigma],'FilterDomain','spatial','FilterSize',d);
        end
%         diff=(BOLD_smoothed_matlab(:,:,1)-BOLD_smoothed(:,:,1));
%         norm(diff(:))
        
%         BOLD_smoothed=BOLD_smoothed.*repmat(mask,[1,1,1,t]); %option for 2D smoothing         
        figure        
        for iz=1:z
            subplot(2,8,iz)
            imagesc(BOLD_smoothed_matlab(:,:,iz,1))
        end
        BOLD_smooth_array=reshape(BOLD_smoothed_matlab,x*y*z,t)';
        
        %% Global Signal Regression
        global_signal=mean(BOLD_smooth_array,2);
        global_signal=global_signal-mean(global_signal);
        b_array=zeros(x*y*z,3);
        for ipixel=1:x*y*z
            if sum(BOLD_smooth_array(:,ipixel))~=0
                linear_trend=(1:length(global_signal))'/length(global_signal);
                X=[ones(size(global_signal)),global_signal,linear_trend];
                Y=BOLD_smooth_array(:,ipixel);
                p=regress(Y,X);
                b_array(ipixel,:)=p;
            end
        end 
        BOLD_gs_regressed_array=BOLD_smooth_array-X(:,2:3)*b_array(:,2:3)';
        BOLD_gs_regressed=reshape(BOLD_gs_regressed_array',x,y,z,t);
              
                
        %% temporal filtering
        %Bandpass filter
        BOLD_2b_filter=BOLD_gs_regressed_array; %BOLD_smooth_array;
        BOLD_filtered_array=zeros(size(BOLD_2b_filter));
        
        TR=0.3; %(s)
        fs=1/TR;
        fpass=0.01; %(Hz)
        deltaT=500;
        Fc1=0.0083; Fc2=0.1; %(Hz)
        for ipixel=1:x*y*z
            if sum(BOLD_2b_filter(:,ipixel))~=0     
                signal=BOLD_2b_filter(:,ipixel);
                signal_2bfiltered=[zeros(deltaT,1); signal-mean(signal); zeros(deltaT,1)];
                signal_filtered_zeros=bandpass(signal_2bfiltered,[Fc1 Fc2],fs);                
                BOLD_filtered_array(:,ipixel) = signal_filtered_zeros(deltaT+1:deltaT+length(signal));
            end
        end 
        BOLD_filtered=reshape(BOLD_filtered_array',x,y,z,t);
                
%         Xiaodi's code
%         if contains(LFP_filename,'iso')
% %             fs=2;lo_cutoff=0.01;up_cutoff=0.1;
% %             [ BOLD_filtered_array ] = myFilter( BOLD_gs_regressed_array, fs, lo_cutoff, up_cutoff );
% %             BOLD_filtered=permute(reshape(BOLD_filtered_array,t,48,48),[2,3,1]);
% %             
%             Hd = myIIR1;% 0.01Hz-0.1Hz
%             b=Hd.Numerator;
%             a=Hd.Denominator;
%             A=BOLD_gs_regressed_array-repmat(mean(BOLD_gs_regressed_array,1),t,1);
%             A=[zeros(300,y*x*z);A;zeros(300,y*x*z)];
%             B=filtfilt(Hd.Numerator,Hd.Denominator,A);
%             BOLD_filtered_array=B(300+1:300+t,:);
%             BOLD_filtered=permute(reshape(BOLD_filtered_array,t,y,x,z),[2,3,4,1]);
%             fprintf('iso\n')
%         elseif contains(LFP_filename,'med')
% %             fs=2;lo_cutoff=0.01;up_cutoff=0.25;
% %             [ BOLD_filtered_array ] = myFilter( BOLD_gs_regressed_array, fs, lo_cutoff, up_cutoff );
% %             BOLD_filtered=permute(reshape(BOLD_filtered_array,t,48,48),[2,3,1]);
%             
%             Hd = myIIR2;% 0.01Hz-0.25Hz
%             b=Hd.Numerator;
%             a=Hd.Denominator;
%             A=BOLD_gs_regressed_array-repmat(mean(BOLD_gs_regressed_array,1),t,1);
%             A=[zeros(300,y*x*z);A;zeros(300,y*x*z)];
%             B=filtfilt(Hd.Numerator,Hd.Denominator,A);
%             BOLD_filtered_array=B(301:t+300,:);
%             BOLD_filtered=permute(reshape(BOLD_filtered_array,t,y,x,z),[2,3,4,1]);
%             fprintf('med\n')
%         else
%             fprintf('error!not iso or med')
%         end
        %% quality assurance                
        status(rats_id)=1;
        Y=15;X=40;iz_select=11;
        figure
        subplot(511)
        plot(squeeze(BOLD_raw(Y,X,iz_select,:)))
        title('raw')
        subplot(512)
        plot(squeeze(BOLD_corrected(Y,X,iz_select,:)))
        title('spm corrected')
        subplot(513)
        plot(squeeze(BOLD_smoothed_matlab(Y,X,iz_select,:)))
        title('spatially smoothed')
        subplot(514)
        plot(squeeze(BOLD_gs_regressed(Y,X,iz_select,:)))
        title('gs regressed')
        subplot(515)
        plot(squeeze(BOLD_filtered(Y,X,iz_select,:)))
        title('temporally filtered')

        save([path_write_mat filename '.mat'],'BOLD_raw','BOLD_corrected','BOLD_smoothed','BOLD_gs_regressed','BOLD_filtered');
        fprintf('done\n');
    else
        fprintf('file not found\n');
        status(rats_id)=0;
    end
end