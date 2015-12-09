function [DWT,a,d_hor,d_ver,d_diag] = FWT2_single(im,w)
% Computes a single analysis DWT of a grayscale image
% INPUT
% im: the image to be processed
% w: the tap-filter specifying the wavelet type ( see 'wfilters' )
% 
% OUTPUT
% DWT: same size as im, containing the wavelet coefficients so arranged
% in quadrants of the same size and a quarter the size of im
%
%                        |
%       low-pass-approx  |  horizontal-detail
%                        |
%       -------------------------------------
%                        |
%       vertical-detail  |  diagonal-detail
%                        |
%   
% a,b_hor,d_ver,d_diag: optional output containing the downsampled 
% subbands separate from one another instead of tiled together
% 

for n = 1:size(im,2)   %apply to single columns
    [a_col(:,n) d_col(:,n)] = DWT_Analysis(im(:,n),w);
end

for m = 1:size(a_col,1)   %apply to rows
    [a(m,:) d_ver(m,:)] = DWT_Analysis(a_col(m,:),w);   %get LP and vertical details
end

for m = 1:size(d_col,1)   %apply to rows
    [d_hor(m,:) d_diag(m,:)] = DWT_Analysis(d_col(m,:),w);  %get horizontal and diagonal details
end

% concatenate coefficients in a single image
DWT = cat(1,[a d_hor],[d_ver d_diag]);

end

