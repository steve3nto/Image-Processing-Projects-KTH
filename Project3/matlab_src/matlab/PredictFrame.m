% Predicts next frame from given frame and block motion vectors
% Input
%   imgI : The given frame, must be padded according to size of shifts
%   motionVect : The motion vectors in a [2 x N_blocks] matrix
%   bSize : Size of the block
%   shiftSize: vector containing [dy_max, dx_max]
%
% Ouput
%   predFrame : The motion compensated predicted frame
%


function predFrame = PredictFrame(F, motionVect, bSize, shiftSize)

dy_max = shiftSize(1);
dx_max = shiftSize(2);

[row, col] = size(F);

bCount = 1;
for i = 1+dy_max:bSize:row-bSize+1
    for j = 1+dx_max:bSize:col-bSize+1
        
        % dy is row(vertical) index
        % dx is col(horizontal) index
        % we are reading blocks line by line
        
        dy = motionVect(1,bCount);
        dx = motionVect(2,bCount);
        i_ref = i + dy;
        j_ref = j + dx;
        predF(i-dy_max:i-dy_max+bSize-1,j-dx_max:j-dx_max+bSize-1) = ...
            F(i_ref:i_ref+bSize-1, j_ref:j_ref+bSize-1);
    
        bCount = bCount + 1;
    end
end

predFrame = predF;
