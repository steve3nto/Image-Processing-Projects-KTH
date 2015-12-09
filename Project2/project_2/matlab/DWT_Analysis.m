function [a, d] = DWT_Analysis(s,p)
%wt1D returns the LPF and HPF subband
%of a single pass of wavelet decomposition
%specified by a prototype vector
%
%INPUT
%s is the signal to decompose
%p is a row vector containing the 
%prototype proportional to the scaling reconstruction filter
%
%OUTPUT
%a is the LP part, or scaling function, or approximation coefficients
%d is the HP part, or wavelet function, or detail coefficients

Lo_R = p/norm(p);   %reconstruction LPF
Lo_D = wrev(Lo_R);  %decomposition LPF
Hi_R = qmf(Lo_R);   %reconstruction HPF
Hi_D = wrev(Hi_R);  %decomposition HPF

%symmetrically extend signal to prevent border effect
s_ext = wextend('1D','ppd',s,length(p)-1);

%LPF to get approximation coefficients
a1 = conv(s_ext,Lo_D,'same');
%HPF to get detail coefficients
d1 = conv(s_ext,Hi_D,'same');
% cut out the part due to extension of the input
a_cut = a1(length(p)-1:end-length(p));
d_cut = d1(length(p)-1:end-length(p));
%Downsample by 2
a = downsample(a_cut,2);
d = downsample(d_cut,2);

end

