function [ shiftvalue,background ] = frequency_est_tirf_v2(ft_im,suppress_noise_factor,fmask,if_show_factor,cutoff)
% estimate the modulation vector k0 via phase-only correlation
% if_show_factor=1:show the corresponding cross-correlation figure
% ft_im: separated/reweighed frequency component
[xsize,ysize,p_num]=size(ft_im);
norm_ft=zeros(size(ft_im));
im=norm_ft;
re_f=im;
reference=zeros(xsize,ysize);

for jj=1:p_num
    ft_max=max(max(abs(ft_im(:,:,jj))));
    ft_im(:,:,jj)=ft_im(:,:,jj)./ft_max;
    norm_ft(:,:,jj)=ft_im(:,:,jj)./(suppress_noise_factor+abs(ft_im(:,:,jj)));
    %          norm_ft(:,:,ii,jj)=ft_im(:,:,ii,jj);
    im(:,:,jj)=fftshift(fft2(norm_ft(:,:,jj)));
end
reference=fftshift(fft2(ifftshift(norm_ft(:,:,1))));




for jj=1:3
    temp=conj(reference).*im(:,:,jj);
    re_f(:,:,jj)=ifft2(ifftshift(temp));% phase-only correlation result
end



shiftvalue=zeros(3,2);%the rough estimation of the modulation frequency k0
% phase mask is introduced to provide a better recognitio of the frequency k0

for jj=1:3
    if jj==1
        [shiftvalue(jj,1),shiftvalue(jj,2)]=find(abs(re_f(:,:,jj))==max(max(abs(re_f(:,:,jj)))));
    else
        [shiftvalue(jj,1),shiftvalue(jj,2)]=find(abs(re_f(:,:,jj)).*fmask==max(max(abs(re_f(:,:,jj)).*fmask)));
    end
end


if if_show_factor==1
    cutoff=round(cutoff);
    circ_x=-cutoff:cutoff;
    circ_y(1:2*cutoff+1)=sqrt(cutoff^2-circ_x.^2);
    circ_y(2*cutoff+2:4*cutoff+2)=-flip(circ_y(1:2*cutoff+1));
    circ_x=[circ_x,flip(circ_x)];
    circ_x(4*cutoff+3)=circ_x(1);
    circ_y(4*cutoff+3)=circ_y(1);
    
    circ_x=circ_x+shiftvalue(1,2);
    circ_y=circ_y+shiftvalue(1,1);
    box_x=[3,3,-3,-3,3];
    box_y=[3,-3,-3,3,3];

    figure;imagesc(abs(re_f(:,:,2)));
    hold on;
    plot(circ_x,circ_y,'--g',...
        box_x+shiftvalue(2,2),box_y+shiftvalue(2,1),'-r');
    % the component within the green dashed circle is excluded when estimating k0
    % the center of the red solid box denotes the estimated k0

end

mask=zeros(xsize,ysize);
background=zeros(1,3);
center_value=zeros(1,3);

for jj=1:3
    xx=shiftvalue(jj,1);
    yy=shiftvalue(jj,2);
    center_value(jj)=abs(re_f(xx,yy,jj));
    mask(xx-3:xx+3,yy-3:yy+3)=1;
    mask(xx-1:xx+1,yy-1:yy+1)=0;
    background(jj)=sum(sum(abs(re_f(:,:,jj)).*mask))/40;
    % mean of background. 40=sum(mask(:)).

    mask(xx-3:xx+3,yy-3:yy+3)=0;
end

background=background./center_value;% ratio of the background to the peak


for jj=2:3
    if background(jj)>0.5
        background(jj)=1.7*background(jj);
    end
end


end