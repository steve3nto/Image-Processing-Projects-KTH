function DWT = FWT2(im,w,scale)
% Computes the arbitrary scale DWT of a grayscale image
% INPUT
% im: the image to be processed
% w: the tap-filter specifying the wavelet type ( see 'wfilters' )
% scale: the number of times the DWT is iterated to reach smaller
% resolutions. Must be an integer greater or equal to 1.
% OUTPUT
% DWT: same size as im, containing the wavelet coefficients so arranged
% in quadrants of the same size and a quarter the size of im
%                        |
%       low-pass-approx  |  horizontal-detail
%                        |
%       -------------------------------------
%                        |
%       vertical-detail  |  diagonal-detail
%                        |
%   Low pass approx is further divided in the same manner depending on the
%   scale input


[DWT,app1] = FWT2_single(im,w);

if scale >= 2
    a{1,1} = app1;
    for k=2:scale
        [Dwt{1,k},a{1,k},d_hor{1,k},d_ver{1,k},d_diag{1,k}] = FWT2_single(a{1,k-1},w);
        DWT(1:size(DWT,1)*(0.5)^(k-1),1:size(DWT,2)*(0.5)^(k-1)) = Dwt{1,k};
    end
end

end

