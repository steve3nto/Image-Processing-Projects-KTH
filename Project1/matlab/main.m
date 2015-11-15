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
im_gauss_med = ones(size(im_gauss));
for i=2:size(im_gauss)-1
    for j=2:size(im_gauss)-1
        roi = im_gauss(i-1:i+1,j-1:j+1); %3x3 mask around the pixel
        im_gauss_med(i,j) = median(roi(:));
    end
end


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
title('saltp noise')
subplot(2,2,4);
imhist(im2uint8(im_saltp_lpf./255));
title('result of lpf saltp noise')