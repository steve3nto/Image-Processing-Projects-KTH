% Computes bits given a frame of quantized block DCT coefficients
% using information entropy to estimate the lower bound of the
% bits needed
% 
% INPUT
% F - 8x8 blocked frame of quantized DCT coefficients 
% step_q - integer step used to quantize the DCT coefficients in F 
%
% OUTPUT
% rate - estimated bits needed to encode 1 dct coefficient
%
%
function rate = EntropyRate( F,step_q )
    
video_height = size(F,1);
video_width = size(F,2);

coefs = zeros(8,8,(video_width/8)*(video_height/8));
index = 1;
ww = [1:8];
hh = [1:8];
%stack 8x8 blocks
for w=1:(video_width/8)
    for h=1:(video_height/8)
        coefs(ww,hh,index) = F(8*(h-1)+hh,8*(w-1)+ww);
        index = index+1;
    end
end
%now calculate R
H=zeros(8,8);
    for w=1:8
        for h=1:8
            vals = squeeze(coefs(h,w,:));  %coefs of frequency h,w
            p = hist(vals,min(vals):(step_q/32):max(vals));
            p = p/sum(p);
            % Entropy matrix
            H(h,w) = -sum(p.*log2(p+eps)); %eps added for log of 0 vals
        end
end

%average over all coefficients
rate = mean2(H);

end

