% Prototype program for CR-SIM reconstruction 
% Author:Chen Nanguang, Shen Shuhao
% Copyright@Optical Bioimaging Lab, National University of Singapore

clear all;
close all;
%% read image file
p_num=3;% phase shift times for each pattern orientation
N=178; % Number of z slices

filepath='D:\Manuscript\Line-scan SIM\Image data\mouse_brain_slit_60x\';%replace with your file's path
filename='Mouse_60xdentrite_stack';% the name format should be filename+'_0/120/240d'+'_zslice number'
fileformat='tif';

%% saving file
save_flag=1; % save the image results if save_flag equals 1;

%% parameter of the detection system
lambda=520;% fluorescence emission wavelength (emission maximum). unit: nm
psize=39; % psize=pixel size/magnification power. unit: nm
NA=0.8; %effective numerical aperture

%% parameter for reconstruction
alpha=0.4; % Merged baseband coefficient
wiener_factor=0.05; % Wiener filter factor
mask_factor=0.3;%a high-pass mask (fmask) is utilized to estimate the modulation vector;

show_initial_result_flag=1; % Show initial result if show_initial_result_flag equals 1;

show_corrected_result_flag=1; % Show corrected result if show_corrected_result_flag equals 1;

for z=1:N
    fprintf('Current z slice: No.%d\n',z);
 %% Read raw images from files    
    for jj=1:p_num
        noiseimage_raw(:,:,jj)=...
            double(imread([filepath,filename,'_',num2str((jj-1)*120),'d_', num2str(z), '.',fileformat]));
    end
    [imv,imh] = size(noiseimage_raw(:,:,1));
    noiseimage_raw = noiseimage_raw - min(noiseimage_raw(:));
    noiseimage = [noiseimage_raw zeros(imv,imv-imh,size(noiseimage_raw,3))];
   

    [xsize,ysize]=size(noiseimage(:,:,1));
    [Y,X]=meshgrid(1:ysize,1:xsize);

    PSF_edge = fspecial('gaussian',5,40);

    for jj=1:p_num
        noiseimage(:,:,jj)=edgetaper(noiseimage(:,:,jj),PSF_edge);
    end

    xc=floor(xsize/2+1);% the x-coordinate of the center
    yc=floor(ysize/2+1);% the y-coordinate of the center
    yr=Y-yc;
    xr=X-xc;
    R=sqrt((xr).^2+(yr).^2);% distance between the point (x,y) and center (xc,yc)

    %% Generate the PSF
    pixelnum=xsize;
    rpixel=NA*pixelnum*psize/lambda;
    cutoff=round(2*rpixel);% cutoff frequency
    ctfde=ones(pixelnum,pixelnum).*(R<=rpixel);
    ctfdeSignificantPix=numel(find(abs(ctfde)>eps(class(ctfde))));
    ifftscalede=numel(ctfde)/ctfdeSignificantPix;
    apsfde=fftshift(ifft2(ifftshift(ctfde)));
    ipsfde=ifftscalede*abs(apsfde).^2;
    OTFde=abs(fftshift(fft2(ifftshift(ipsfde))));
    clear apsfde ctfde temp X Y

    %% filter/deconvolution before using noiseimage
    widefield=sum(noiseimage,3); % Obtain widefield image

    separated_FT=zeros(xsize,ysize,3);
    noiseimagef=zeros(size(noiseimage));% Obtain the frequency domain of the raw images
    noiseimagef(:,:,1)=fftshift(fft2(noiseimage(:,:,1)));
    noiseimagef(:,:,2)=fftshift(fft2(noiseimage(:,:,2)));
    noiseimagef(:,:,3)=fftshift(fft2(noiseimage(:,:,3)));

    re0_temp=zeros(xsize,ysize);
    rep_temp=zeros(xsize,ysize);
    rem_temp=zeros(xsize,ysize);
    modulation_matrix=[1,1/2*exp(-1i*(pi*0)),1/2*exp(1i*(pi*0));...
        1,1/2*exp(-1i*(pi*2/3)),1/2*exp(1i*(pi*2/3));...
        1,1/2*exp(-1i*(pi*4/3)),1/2*exp(1i*(pi*4/3))];
    matrix_inv=inv(modulation_matrix);

    for jj=1:p_num

        re0_temp=matrix_inv(1,jj)*noiseimagef(:,:,jj)+re0_temp;
        rep_temp=matrix_inv(2,jj)*noiseimagef(:,:,jj)+rep_temp;
        rem_temp=matrix_inv(3,jj)*noiseimagef(:,:,jj)+rem_temp;
    end

    separated_FT(:,:,1)=re0_temp;
    separated_FT(:,:,2)=rep_temp;
    separated_FT(:,:,3)=rem_temp;

    clear re0_temp rep_temp rem_temp

    fmask=double(sqrt(xr.^2+yr.^2)>cutoff*mask_factor);
    [shiftvalue,~]=frequency_est_tirf_v2(separated_FT,0.8,fmask,show_initial_result_flag,mask_factor*cutoff);
    clear separated_FT

    shiftvalue(2,:)=shiftvalue(2,:)-shiftvalue(1,:);
    shiftvalue(3,:)=shiftvalue(3,:)-shiftvalue(1,:);
    shiftvalue(1,1)=0;
    shiftvalue(1,2)=0;


    %% phase correction with inverse matrix based algorithm
    search_range=0.4;%the max radius in the local search algorithm

    %obtain a more precise estimation of the period and the directon of sinusodial pattern
    [ precise_shift,~] = precise_frequency_tirf(noiseimagef,shiftvalue,search_range);

    % estimation the phase of each pattern

    [inv_phase] = separation_matrix_correction_v3(noiseimagef,precise_shift,OTFde);
    %% cross-correlation based algorithm
    sigma=0.1;
    % [cc_phase]=crosscorrelation_phase_est_SIM(noiseimagef,precise_shift,sigma,OTFde);

    %% auto-correlation based algorithm
    auto_phase=zeros(1,p_num);

    for jj=1:p_num
        f_temp=exact_shift(noiseimagef(:,:,jj),...
            [-precise_shift(2,1),-precise_shift(2,2)],1);

        auto_phase(jj)=angle(sum(sum(conj(noiseimagef(:,:,jj)).*f_temp)));
    end


    my_phase_temp=mod(-inv_phase,2*pi);
    my_phase_auto=mod(-auto_phase,2*pi);


    % inv_phase=auto_phase;
    % reconstruct with phases determined by the auto-correlation method

    % inv_phase=-cc_phase;
    % reconstruct with phases determined by the cross-correlation method
    % end



    %% separate different frequency component
    % n_filt is a notch-filter
    n_filt = 1 - exp(-0.05*R.^1.1);
    separated_FT=zeros(xsize,ysize,3);% store different bands of frequency component

    re0_temp=zeros(xsize,ysize);
    rep_temp=zeros(xsize,ysize);
    rem_temp=zeros(xsize,ysize);
    mi=0.5;
    modulation_matrix=[1,mi*exp(-1i*(inv_phase(1))),mi*exp(1i*(inv_phase(1)));...
        1,mi*exp(-1i*(inv_phase(2))),mi*exp(1i*(inv_phase(2)));...
        1,mi*exp(-1i*(inv_phase(3))),mi*exp(1i*(inv_phase(3)))];

    matrix_inv=inv(modulation_matrix);
    for jj=1:p_num
        re0_temp=matrix_inv(1,jj)*noiseimagef(:,:,jj)+re0_temp;
        rep_temp=matrix_inv(2,jj)*noiseimagef(:,:,jj)+rep_temp;
        rem_temp=matrix_inv(3,jj)*noiseimagef(:,:,jj)+rem_temp;
    end


    separated_FT(:,:,1)=re0_temp;
    separated_FT(:,:,2)=rep_temp.*n_filt;
    separated_FT(:,:,3)=rem_temp.*n_filt;

    %Obtaining optical-sectioned baseband
    tempf = (noiseimagef(:,:,1)-noiseimagef(:,:,2)).*n_filt;

    refm1 = repmat([0:ysize-1],xsize,1);
    refm2 = repmat([0:xsize-1]',1,ysize);
    refm = cos(2*pi*(refm1*precise_shift(2,2)/ysize+refm2*precise_shift(2,1)/xsize)+inv_phase(1));
    refm = refm - cos(2*pi*(refm1*precise_shift(2,2)/ysize+refm2*precise_shift(2,1)/xsize)+inv_phase(2));


    temps = real(ifft2(fftshift(tempf)));
    tempsdm = refm.*temps;
    tempf1 = fftshift(fft2(tempsdm));
    separated_FT(:,:,1) = tempf1.*(1-n_filt)+alpha*separated_FT(:,:,1).*n_filt;% Modified baseband

    [~,noise_ratio]=frequency_est_tirf_v2(separated_FT,0.8,fmask,show_corrected_result_flag,mask_factor*cutoff);

    clear noiseimagef
    
    clear re0_temp rem_temp rep_temp R X Y xr yr OTF_temp
    %% Frequency shift
    OTF_n=zeros(size(OTFde));
    OTF_nb=OTF_n;% OTF of reconstruct image
    ft_true=zeros(size(separated_FT));
    OTF_de_temp=abs(fftshift(fft2(ipsfde)));
    OTFcirc=double(OTFde./(wiener_factor^2+OTFde));

    for jj=1:3
        if jj~=1
            ft_true(:,:,jj)=exact_shift(separated_FT(:,:,jj),...
                [precise_shift(jj,1),precise_shift(jj,2)],1);

        else
            ft_true(:,:,jj)=separated_FT(:,:,jj);
        end

        OTF_n=circshift(OTF_de_temp,[shiftvalue(jj,1),shiftvalue(jj,2)])+OTF_n;
        OTF_nb=circshift(OTFcirc,[shiftvalue(jj,1),shiftvalue(jj,2)])+OTF_nb;
    end

    clear separated_FT
    OTF_n=OTF_n./max(max(abs(OTF_n)));
    OTF_nb=OTF_nb./max(max(OTF_nb));

    psf_nb=abs(fftshift(fft2(OTF_nb))); %the efficient PSF for SIM with pre-deconvolution

    psf_n=fftshift(ifft2(ifftshift(OTF_n)));
    psf_n=abs(psf_n); %the efficient PSF for SIM without pre-deconvolution

    mod_depth_temp=zeros(1,3);
    reference=ft_true(:,:,1);

    for jj=1:3
        if jj==1          
            mask_temp=OTFcirc.*circshift(OTF_de_temp,[shiftvalue(2,1),shiftvalue(2,2)]);
            mask_temp=mask_temp./(OTF_de_temp+wiener_factor^2);
            %without pre-deconvolution
            mod_depth_temp(jj)=sum(sum(conj(reference).*ft_true(:,:,jj).*mask_temp));
        else
            mod_depth_temp(jj)=sum(sum(conj(reference).*ft_true(:,:,jj)));
        end
        temp=mod_depth_temp(jj)./abs(mod_depth_temp(jj));
        if jj~=1
            ft_true(:,:,jj)=ft_true(:,:,jj)./temp;
        end
    end


    %% estimate the modulation depth via cross-correlation
    abs_mod=abs(mod_depth_temp);
    illu_max=max(abs_mod,[],2);

    abs_mod=abs_mod./illu_max(1);

    % the following is the best set of parameter among 10 sets. (Empirical value)
    illu_max=illu_max./max(illu_max);
    modulation_depth=mean(abs_mod(:,2:3),2);% modulation depth for each pattern orientation
    noise_ratio=1./noise_ratio;
    noise_suppress=noise_ratio(:,2)./max(noise_ratio(:,2));


    for jj=2:3
        myweight_factor=max([min([1/modulation_depth,2.5]),1])...
            *noise_suppress;
        ft_true(:,:,jj)=ft_true(:,:,jj)*myweight_factor;
    end

    %% Reconstruction
    FT_extended=sum(ft_true,3);
%   FT_extended(:,[1:800 end-800:end])=0;%Further filtering high-freq
    sft = shiftvalue(2,2);
    n_filt1 = [ones(xsize,sft), n_filt(1:xsize,1:end-sft)];
    n_filt1 = n_filt1.*[n_filt(1:xsize,sft+1:end),ones(xsize,sft)];
    n_filt1 = n_filt1.*[n_filt(1:xsize,sft*2+1:end),ones(xsize,sft*2)];
    n_filt1 = n_filt1.*[ones(xsize,sft*2), n_filt(1:xsize,1:end-sft*2)];
    
    reconstructed_im=real(ifft2(ifftshift(FT_extended.*(0.1+0.9*n_filt1))));
    reconstructed_im=real(reconstructed_im).*(real(reconstructed_im>0));
    reconstructed_im = reconstructed_im(1:imv,1:imh);
    reconstructed_im = imresize(reconstructed_im, [imv/2 imh]);
    reconstructed_im = imfilter(reconstructed_im, fspecial('gaussian',10,1));
    reconstructed_im = reconstructed_im.*(reconstructed_im>0);
    reconstructed_im=deconvlucy(reconstructed_im,psf_n(950:970,950:970),4);
    crsim(:,:,z)=reconstructed_im;
    % figure;imshow(reconstructed_im,[]);title('SIM');
  
    % widefield=deconvlucy(widefield,ipsfde,3);
    widefield1=deconvlucy(widefield,ipsfde,3);
 
    widefield1 = widefield1(1:imv,1:imh);
     widefield1 = imresize(widefield1, [imv/2 imh]);
    lscm(:,:,z)=widefield1;
    % figure;imshow(widefield1,[]);title('wide-field');

    close all;
end

if save_flag==1

    crsim = uint16(crsim/max(crsim(:))*65535);
    lscm = uint16(lscm/max(lscm(:))*65535);

    for i=1:N
        disp(i);
        if i==1
            imwrite(crsim(:,:,i),[filepath,'CRSIM_' filename '_stack.tif']);
            imwrite(lscm(:,:,i),[filepath,'LSCM_' filename '_stack.tif']);
        else
            imwrite(crsim(:,:,i),[filepath,'CRSIM_' filename '_stack.tif'],'WriteMode', 'append');
            imwrite(lscm(:,:,i),[filepath,'LSCM_' filename '_stack.tif'],'tif','WriteMode', 'append');
        end
    end
end

