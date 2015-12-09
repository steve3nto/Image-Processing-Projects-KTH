% Wavelets
clear
im = im2double(imread('peppers512x512.tif'));
s = im(56,:);
load haar
W = haar;  %prototype

%generate filters from prototype
% Lo_R = W/norm(W);   %reconstruction LPF
% Lo_D = wrev(Lo_R);  %decomposition LPF
% Hi_R = qmf(Lo_R);   %reconstruction HPF
% Hi_D = wrev(Hi_R);  %decomposition HPF
% 
% %extend signals to remove edge problems
% s_ext = wextend('1D','sym',s,size(W,2));
% 
% scale = conv(s_ext,Lo_D,'same');
% scale = scale(1+size(W,2) : end-size(W,2));
% scale = downsample(scale,2);
% wave = conv(s_ext,Hi_D,'same');
% wave = wave(1+size(W,2) : end-size(W,2));
% wave = downsample(wave,2);
% 
% scale_ups = upsample(scale,2);
% scale_ups_ext = wextend('1D','sym',scale_ups,size(W,2));
% rec_lp = conv(scale_ups_ext,Lo_R,'same');
% rec_lp = rec_lp(1+size(W,2) : end-size(W,2));
% wave_ups = upsample(wave,2);
% wave_ups_ext = wextend('1D','sym',wave_ups,size(W,2));
% rec_hp = conv(wave_ups_ext,Hi_R,'same');
% rec_hp = rec_hp(1+size(W,2) : end-size(W,2));
% % eliminate the borders due to the tails of the convolution
% % rec_lp = rec_lp(2*size(W,2):end-2*size(W,2));
% % rec_hp = rec_hp(2*size(W,2):end-2*size(W,2));
% 
% s_rec = rec_lp + rec_hp;
% plot(s);
% hold
% plot(s_rec,'r');

%symmetrically extend signal to prevent border effect
%s_ext = wextend('1D','sym',s,length(W)-1);

[scaling1, wavelet1] = DWT_Analysis(s,W);
[scaling2, wavelet2] = DWT_Analysis(scaling1,W);
[scaling3, wavelet3] = DWT_Analysis(scaling2,W);
[scaling4, wavelet4] = DWT_Analysis(scaling3,W);
[scaling5, wavelet5] = DWT_Analysis(scaling4,W);
[scaling6, wavelet6] = DWT_Analysis(scaling5,W);
[scaling7, wavelet7] = DWT_Analysis(scaling6,W);
[scaling8, wavelet8] = DWT_Analysis(scaling7,W);



rec7 = DWT_Synthesis(scaling8,wavelet8,W);
rec6 = DWT_Synthesis(rec7,wavelet7,W);
rec5 = DWT_Synthesis(rec6,wavelet6,W);
rec4 = DWT_Synthesis(rec5,wavelet5,W);
rec3 = DWT_Synthesis(rec4,wavelet4,W);
rec2 = DWT_Synthesis(rec3,wavelet3,W);
rec1 = DWT_Synthesis(rec2,wavelet2,W);
rec  = DWT_Synthesis(rec1,wavelet1,W);

%rec = iwt1D(scaling1,wavelet1,W);

% cut out the part due to extension of the input
%rec = rec(length(W)-1:end-length(W));
plot(s)
hold
plot(rec,'r')

%% 2D part
clear
im = im2double(imread('peppers512x512.tif'));
s = im(1,:);
load db4
wavelet = db4;  %prototype


% for n = 1:size(im,2)   %apply to single columns
%     [a_col(:,n) d_col(:,n)] = wt1D(im(:,n)',wavelet);
% end
% 
% for m = 1:size(a_col,1)   %apply to rows
%     [a(m,:) d_ver(m,:)] = wt1D(a_col(m,:),wavelet);   %get LP and vertical details
% end
% 
% for m = 1:size(d_col,1)   %apply to rows
%     [d_hor(m,:) d_diag(m,:)] = wt1D(d_col(m,:),wavelet);  %get horizontal and diagonal details
% end


% pack the coefficients together
DWT = FWT2(im,wavelet,1);
im_rec = IFWT2_single(DWT,wavelet);

% % Inverse DWT
% a = DWT( 1:size(DWT,1)*0.5 , 1:size(DWT,2)*0.5 );
% h_hor = DWT( 1:size(DWT,1)*0.5 , 1+size(DWT,2)*0.5:end );
% h_ver = DWT( 1+size(DWT,1)*0.5:end , 1:size(DWT,2)*0.5 );
% h_diag = DWT( 1+size(DWT,1)*0.5:end , 1+size(DWT,2)*0.5:end );
% 
% %reconstruct LP horizontally
% for m = 1:size(a,1)   %apply to row one-by-one
%      reclp_hor(m,:) = DWT_Synthesis(a(m,:),h_hor(m,:),wavelet);   
% end
% % reconstruct HP horizontally
% for m = 1:size(a,1)   %apply to row one-by-one
%      rechp_hor(m,:) = DWT_Synthesis(h_ver(m,:),h_diag(m,:),wavelet);   
% end
% %reconstruct final image
% for n = 1:size(reclp_hor,2)   %apply to colun one-by-one
%     im_rec(:,n) = DWT_Synthesis(reclp_hor(:,n),rechp_hor(:,n),wavelet);
% end

%FUUUUCK
[CA,CH,CV,CD] = dwt2(im,'db10');
rec = idwt2(CA,CH,CV,CD,'db10');

subplot(1,2,1)
imshow(rec);
title('Matlab function')

subplot(1,2,2)
imshow(im_rec)
title('My implementation')