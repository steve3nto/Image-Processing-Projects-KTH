%% 2.1 Histogram Equalization
im = imread('lena512.bmp');
pdf_original = histcounts(im(:),[0:256]);
figure(1)
imshow(im);
figure(2)
imhist(im);

% reduce the contrast of the image by 
% linearly mapping the range [0, 255] to [b, b+a*255]
a = 0.2;
b = 50;
% these values map [0,255] to [50,101]
for i=1:size(im,1)
    for j=1:size(im,2)
        im_lowc(i,j) = min(max(a*im(i,j) + b, 0), 255);
    end
end
figure(3)
imshow(im_lowc);
figure(4)
imhist(im_lowc);

% Perform histogram equalization to improve the low contrast image
pdf = histcounts(im_lowc(:),[0:256]);
cdf = cumsum(pdf);
im_eq = uint8(zeros(512,512));
for i=1:size(im,1)
    for j=1:size(im,2)
        im_eq(i,j) = 255*cdf(im_lowc(i,j)+1)/(size(im_lowc,1)*size(im_lowc,2));
    end
end

% im_eq = im2uint8(im_eq);
pdf_eq = histcounts(im_eq(:),[0:256]);
figure(5)
imshow(im_eq);
figure(6)
imhist(im_eq);
% The histogram is not flat because discrete equalization is a one to one
% mapping between values(ignoring rounding that can merge together some bins)
% thus the histogram has the same "shape" (peaks are identifiable),
% but it is stretched over the whole dynamic range

%% 2.2 Image denoising

im = 255*im2double(imread('lena512.bmp'));

% add gaussian noise
n_gauss = sqrt(64).*randn(size(im));  %zero mean variance 64 gaussian noise
im_gauss = im + n_gauss;

% add salt-pepper noise
im_saltp = im;
n = mynoisegen('saltpepper', 512, 512, .05, .05);
im_saltp(n==0) = 0;
im_saltp(n==1) = 255;

% 3x3 mean filter with two vector kernels (it's a separable filter)
h_horizontal = (1/3)*ones(1,3);
h_vertical = h_horizontal';
im_gauss_lpf = conv2(h_vertical,h_horizontal,im_gauss,'same');

im_saltp_lpf = conv2(h_vertical,h_horizontal,im_saltp,'same');

% median filter
im_gauss_med = medfilt2(im_gauss);
im_saltp_med = medfilt2(im_saltp);


% Plot all figures and compare
figure(1)
imshow(im,[0 255]);
title('original image')
figure(2)
subplot(2,2,1);
imshow(im_gauss,[0 255]);
title('gaussian noise')
subplot(2,2,2);
imshow(im_gauss_lpf,[0 255]);
title('result of lpf of gaussian noise')
subplot(2,2,3);
imshow(im_saltp,[0 255]);
title('saltp noise')
subplot(2,2,4);
imshow(im_saltp_lpf,[0 255]);
title('result of lpf saltp noise')

% Plot all histograms and compare
figure(3)
imhist(im2uint8(im./255));
axis([0 255 0 max(hist(im(:),255))+50])
title('original histogram')
figure(4)
subplot(2,2,1);
imhist(im2uint8(im_gauss./255));
title('gaussian noise')
subplot(2,2,2);
imhist(im2uint8(im_gauss_lpf./255));
title('result of lpf of gaussian noise')
subplot(2,2,3);
imhist(im2uint8(im_saltp./255));
axis([-5 260 0 max(hist(im_saltp(:),255))+50])
title('saltp noise')
subplot(2,2,4);
imhist(im2uint8(im_saltp_lpf./255));
title('result of lpf saltp noise')

%% Frequency Domain Filtering
im = double(imread('lena512.bmp'));
% generate the blur kernel
blur = myblurgen('gaussian',8);
im_blur = conv2(im,blur,'same');
figure
imshow(im_blur,[]);
%remove edges
% im_blur_2 = im_blur(10:size(im_blur,1)-10,10:size(im_blur,2)-10);
% figure(2)
% imshow(im_blur_2,[]);
% just removing the edges doesn't remove discontinuities
% in the periodic repetition of the image

% %Apply windowing
% [M, N] = size(im);
% w1 = cos(linspace(-pi/2, pi/2, M));
% w2 = cos(linspace(-pi/2, pi/2, N));
% w = w1' * w2;
% 
% w1 = hanning(M);
% w2 = hanning(N);
% w = w1*w2';
% 
% im = im.*w;
% im_blur = im_blur.*w;
% figure
% imshow(im_blur,[]);

%plot magnitude spectra of original image

im_fft = fftshift(fft2(im));
A = abs(im_fft);
L = log(A + 1);
figure
imshow(L,[]);
title('magnitude spectrum of original img')
%plot magnitude spectrums
im_fft_blur = fftshift(fft2(im_blur));
A_blur = abs(im_fft_blur);
L_blur = log(A_blur + 1);
figure
imshow(L_blur,[]);
title('magnitude spectrum of blurred img')

% reconstruction
im_out = ifft2(ifftshift(im_fft));
im_blur_out = ifft2(ifftshift(im_fft_blur));
figure
imshow(im_out,[]);
title('original image');
figure
imshow(im_blur_out,[]);
title('blurred image');

% weiner filtering on distorted (blurred + nosied) image
%copy this to separate function

im = double(imread('boats512_outoffocus.bmp'));
b = myblurgen('gaussian',8); %% blur function
%g = imfilter(im,b,'same'); %% degraded image (w/o noise)
g = im;

noise_var = 0.0833; %% noise variance
noise = sqrt(noise_var).*randn(size(g)); %n represents quanization error
%g = g + noise;

figure
imshow(g,[]);
title('blurred image');

%starting here:
%function wiener = weiner(g,n,var)

%fss = fft(autocorr(g(:)));
%fnn = fft(autocorr(noise(:)));
%h = fss./(fss+fnn);

%nope
%imshow(conv2(g,h,'same'),[]);

imvar = var(g(:));

PSF = fspecial('gaussian',60,10);
g = edgetaper(g,PSF);

wfiltered = deconvwnr(g,b,noise_var/imvar);

figure
imshow(wfiltered,[0,255]);  % not sure why are we getting these edge artifacts
title('restored');

%try of own implementation:
K = noise_var/imvar; %SNR ratio estimation (arbitrary)
H = fftshift(fft2(b, 1024, 1024)); %fft of blur (passed as argument)
H_mag = H.*conj(H);

%filter = (1./H) .* (H_mag./(H_mag+K)); %more computations
filter = conj(H) ./ (H_mag + K);

out_freq = filter .* fftshift(fft2(g,1024,1024));

out = ifft2(fftshift(out_freq));
out = out(1:512,1:512);

figure
imshow(out,[0,255]);
title('self restored');
%end

