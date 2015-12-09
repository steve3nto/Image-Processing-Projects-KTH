function im_rec = IFWT2_single(Dwt,w)
% Computes a single synthesis of 2D-DWT from the image DWT coefficients
%                        |
%       low-pass-approx  |  horizontal-detail
%                        |
%       -------------------------------------
%                        |
%       vertical-detail  |  diagonal-detail
%                        |
%   
% INPUT
% Dwt: the image containing the DWT coefficients as seen above
% w: the tap-filter specifying the wavelet type ( see 'wfilters' )
% 
% OUTPUT
% im_rec: the reconstructed image
%

% unpack coefficients
a = Dwt( 1:size(Dwt,1)*0.5 , 1:size(Dwt,2)*0.5 );
h_hor = Dwt( 1:size(Dwt,1)*0.5 , 1+size(Dwt,2)*0.5:end );
h_ver = Dwt( 1+size(Dwt,1)*0.5:end , 1:size(Dwt,2)*0.5 );
h_diag = Dwt( 1+size(Dwt,1)*0.5:end , 1+size(Dwt,2)*0.5:end );

% reconstruct LP horizontally
for m = 1:size(a,1)   %apply to row one-by-one
     reclp_hor(m,:) = DWT_Synthesis(a(m,:),h_hor(m,:),w);   
end
% reconstruct HP horizontally
for m = 1:size(a,1)   %apply to row one-by-one
     rechp_hor(m,:) = DWT_Synthesis(h_ver(m,:),h_diag(m,:),w);   
end
% reconstruct final image
for n = 1:size(reclp_hor,2)   %apply to colun one-by-one
    im_rec(:,n) = DWT_Synthesis(reclp_hor(:,n),rechp_hor(:,n),w);
end

end

