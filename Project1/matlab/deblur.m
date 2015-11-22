function [f f_med] = deblur(g, b, n_var)
%Uses Wiener Filtering to deblur a noisy and blurred image
%   g in the noisy and blurred image
%   b is the convolution kernel for image deblurring
%   var is the variance of the noise present in f
%   f is the deblurred image
%   f_med is an optional alternative with a 3-by-3 median filter post-processing
%   to remove residual noise after deblurring

im_var = var(g(:)); %used to estimate the NSR
K = n_var / im_var; %NSR ratio estimation assuming it is constant

% taper the edges of the image to reduce ringing
PSF = fspecial('gaussian',60,10);   %large gaussian window is used
g_taper = edgetaper(g,PSF);

% Pad with zeros to avoid artifacts and increase fft resolution
H = fft2(b, 1024, 1024); %fft of blur 
H_mag = H.*conj(H);

filter = conj(H) ./ (H_mag + K);

% use symmetric padding to have the same fft resolution as H (1024 fft bins)
g_padded = padarray(g_taper,[256 256], 'symmetric');

out_freq = filter .* fft2(g_padded);

f = ifft2(out_freq);
% Apply 3-by-3 median filter to remove residual noise
f_med = medfilt2(f,[3 3]);
% crop only one replica of the image
f = f(249:760,249:760);
f_med = f_med(249:760,249:760);
end
