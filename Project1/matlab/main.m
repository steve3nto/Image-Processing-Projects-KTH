% 2.1 Histogram Equalization
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

% 2.2 Image denoising

