function [rec, lp_rec, hp_rec] = DWT_Synthesis(a,d,p)
%iwt1D reconstructs a signal from the wavelet decomposed and downsampled
%scaling and wavelet functions generated with the prototype vector p
%
%INPUT
%a is the vector containing the approximation coeff or scaling function
%d is the vector containing the detail coefficients or wavelet function
%p is a row vector containing the prototype proportional to the
%scaling reconstruction filter
%
%OUTPUT
%rec is the reconstructed approximation signal
%lp_rec is the optional output, containing only the LP approximation part
%hp_rec is the optional output containing only the HP detail part

Lo_R = p/norm(p);   %reconstruction LPF
Lo_D = wrev(Lo_R);  %decomposition LPF
Hi_R = qmf(Lo_R);   %reconstruction HPF
Hi_D = wrev(Hi_R);  %decomposition HPF

% Check that approximation and detail vectors have same sizes
% and remove one sample if needed
if length(a) > length(d)
    a = a(1:end-1);   
elseif length(a) < length(d)
    d = d(1:end-1);     
end

% upsample
a_up = upsample(a,2);
d_up = upsample(d,2);
% symmetrically extend signal to prevent border effect
a_ext = wextend('1D','ppd',a_up,length(p)-1);
d_ext = wextend('1D','ppd',d_up,length(p)-1);
% Filter 
lp_rec = conv(a_ext,Lo_R,'same');
hp_rec = conv(d_ext,Hi_R,'same');
% sum to get the approximation
rec_uncut = lp_rec + hp_rec;
% cut out the part due to extension of the input
rec = rec_uncut(length(p):end-length(p)+1);

end

