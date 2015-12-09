
%% Task 2 - DCT based Image compression
%
x = im2double(imread('peppers512x512.tif'));

%Function definitions
uniform_quantizer = @(x,ssize) round(x/ssize)*ssize;
mse = @(x,y) sum(sum((y-x).^2))/(size(y,1) * size(y,2));
PSNR = @(D) 10*log10(255^2./D);

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
ss = 6; %step size
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

i=4; %for testing
for i=1:size(ssizes)
    s=ssizes(i);
    y_coeff = uniform_quantizer(x_dct,s);
    reconstruction_MSE = mse(x_dct, y_coeff);
    reconstruction_PSNR = PSNR(reconstruction_MSE);
    psnrs(i) = reconstruction_PSNR;
    y = blockproc(y_coeff,[8 8],@idctz2);
    
    %todo process blocks of all images?
    %calculate rates (hardcoded)
    coefs = zeros(8,8,64*64);
    for w=1:8
        for h=1:8
            index = 1;
            for k=0:63
                for l=0:63
                    %get only DC coeffs
                    coefs(w,h,index) = y_coeff((8*k)+w,(8*l)+h);
                    index = index + 1;
                end
            end
        end
    end
    
    %now calculate R
    H=zeros(8,8);
    for w=1:8
        for h=1:8
            vals = squeeze(coefs(w,h,:));
            p = hist(vals,min(vals):s:max(vals));
            p = p/sum(p);
            
            H(w,h) = -sum(p.*log2(p+eps)); %fix eps
        end
    end
    
    rates(i) = mean2(H);
    
    %figure;
    %imshow(blockproc(y_coeff,[8 8],@idctz2),[]);
end


plot(rates, psnrs, '+-', 'linewidth', 2);
grid on;

%% Task 3
% filters are generated inside FWT2 using the DWTAnalysis and DWTSynthesis
% functions

im = 255*im2double(imread('harbour512x512.tif'));
s = im(1,:);
load db4  
wavelet = db4;  %prototype for the 8-tap daubechies filters
 
DWT = FWT2(im,wavelet,4);
% Show scale 4 DWT coefficients
imshow(DWT/255);
 
%Uniform Quantization of wavelet transform
Lo_R = wavelet/norm(wavelet);   %reconstruction LPF
Lo_D = wrev(Lo_R);  %decomposition LPF
Hi_R = qmf(Lo_R);   %reconstruction HPF
Hi_D = wrev(Hi_R);  %decomposition HPF
 
[CA,CH,CV,CD] = dwt2(im,Lo_R,Hi_R);
rec = idwt2(CA,CH,CV,CD,Lo_D,Hi_D);
 
% step size for the quantizer, smaller step size is better
stepq = 2.^[0 1 2 3 4 5 6 7 8 9];
step_count = [1:length(stepq)];
 
% Uniformely quantize values
for k = step_count
        CAq(:,:,k) = uniform_quantizer(CA,stepq(k));
        CHq(:,:,k) = uniform_quantizer(CH,stepq(k));
        CVq(:,:,k) = uniform_quantizer(CV,stepq(k));
        CDq(:,:,k) = uniform_quantizer(CD,stepq(k));
end
% Reconstruct images
for k = step_count
        recq(:,:,k) = idwt2(CAq(:,:,k),CHq(:,:,k),CVq(:,:,k),CDq(:,:,k),Lo_D,Hi_D);
        % compute mse
        mserr(k) = mse(im,recq(:,:,k));
end
 
mserrdb_wav = 10*log10(mserr);
Psnr_wav = PSNR(mserr);