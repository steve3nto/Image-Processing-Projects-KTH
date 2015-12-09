
%% Task 2 - DCT based Image compression
%
x = im2double(imread('peppers512x512.tif'));

%Function definitions
uniform_quantizer = @(x,ssize) round(x/ssize)*ssize;
mse = @(x,y) sum(sum((y-x).^2))/(size(y,1) * size(y,2));
PSNR = @(D) 10*log10(255^2/D);

%DCT-1 and Inverse DCT of 8x8 white image
white_8x8 = ones(8);

ct = dctz(white_8x8); %dctz and idctz work with 8x8 input matrix hardcoded
ict = idctz(ct);

%Process the loaded image block by block

%dct version
% block_ct = blockproc(im,[8 8],@dctz);
% block_ict = blockproc(block_ct,[8 8],@idctz);
%
% figure;
% imshow(block_ct);
% figure
% imshow(block_ict);

%DCT-2
block_dct = blockproc(x,[8 8],@dctz2);
block_idct = blockproc(block_dct,[8 8],@idctz2);

figure;
imshow(x);
title('Original Image peppers 512x512');
figure;
imshow(block_dct);
title('DCT coefficients');
figure
imshow(block_idct);
title('Reconstruction without quantization');

%reconstructed image error before applying quantization
reconstruction_MSE_without_quantization = mse(block_idct, x)
PSNR(reconstruction_MSE_without_quantization)

%Applying uniform quantizer with step size of 1
ss = 1; %step size
y_coeff = uniform_quantizer(block_dct,ss);

y = blockproc(y_coeff,[8 8], @idctz2);

%reconstructed image error after applying quantization
reconstruction_MSE_with_quantization = mse(y, x)
PSNR(reconstruction_MSE_with_quantization)

% due to quantization errors the output can contain values "<0" and ">1"
% therefore we floor/ceil output range to [0,1] since we are working
% with doubles
y_bound = y;
y_bound( y_bound>1 ) = 1;
y_bound( y_bound<0 ) = 0;
reconstruction_MSE_with_quantization_and_bound = mse(y_bound, x)
PSNR(reconstruction_MSE_with_quantization_and_bound)

figure;
imshow(y);
title('Reconstructed image after quantization with step size 1');

%mse of coeffictients (before, after quatization)
coeffs_mse = mse(block_dct, y_coeff)

% to calculate real PSNR we have to switch back from [0,1] to [0,255] scale
%    we are using y_bound as output (this has values clipped to range [0,1]
y = y_bound;

x255 = round(255*x);
y255 = round(255*y); 

%D = mse(y,x); % PSNR=72dB (for ss=1), where x,y are in range [0,1]
D = mse(y255,x255); % PSNR=24dB (for q=1) which is more realistic
PSNR(D)


%% Distortion and Bit-rate estimation
%
% TODO: commment why they are the same (orthonormal transform?)
%

x = im2double(imread('peppers512x512.tif'));
x = 255*x;

x_dct = blockproc(x,[8 8],@dctz2);

ssizes = [1 2 4 8 16 32 64 128 256 512]';
psnrs = zeros(size(ssizes));
rates = zeros(size(ssizes));

for i=1:size(ssizes)
    s=ssizes(i);
    y_coeff = uniform_quantizer(x_dct,s);
    reconstruction_MSE = mse(x_dct, y_coeff);
    reconstruction_PSNR = PSNR(reconstruction_MSE);
    psnrs(i) = reconstruction_PSNR;
    y = blockproc(y_coeff,[8 8],@idctz2);
    
    %todo process blocks of all images?
    %calculate rates
    
    %figure;
    %imshow(blockproc(y_coeff,[8 8],@idctz2),[]);
end

plot(psnrs,ssizes)
