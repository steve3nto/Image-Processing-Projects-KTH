function g = deblur(f, var)
%deblur Uses Wiener Filtering to deblur a noisy and blurred image
%   f in the noisy and blurred image
%   var is the variance of the noise present in f
%   g is the deblurred image



end

h = 1/9.*ones(3);
H = fft2(h,512,512);
yoo = log(abs(fftshift(H))+1);
imshow(yoo,[]);