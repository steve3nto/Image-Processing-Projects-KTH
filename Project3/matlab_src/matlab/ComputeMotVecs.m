% Computes the best block-based motion vectors for a single frame I, 
% using squared blocks of size bSize and whole pixel shifts from
% -10 to 10 in both x and y direction
% The cost used to find the best motion vector for every block is the MSE
%
% INPUT
% F - The current frame, must be padded around according to size of shifts
% nextF - the next frame to be predicted
% bSize - the size of the blocks
% shifts - [2 x N_SHIFTS] matrix with all possible shifts over which to
% search. Columns are of the form [dy dx]'
%
% OUTPUT
% mVecs - matrix of motion vectors with dy,dx along the vertical dimension
% indexes - vector of indices to quickly compute the entropy of the shifts

function [mVecs,mVecsIndexes] = ComputeMotVecs(F, nextF, bSize, shifts)

[vid_height, vid_width] = size(nextF);
dy_max = max(shifts(1,:));
dx_max = max(shifts(2,:));

err = zeros(vid_height/bSize,vid_width/bSize,length(shifts));
mVecsIndexes = zeros(1,vid_height/bSize*vid_width/bSize);
mVecs = zeros(2,vid_height/bSize*vid_width/bSize);
index = 1;
% compute prediction error (MSE) for all blocks and motion vectors
% and keep the shift giving the minimum error
for i=1:vid_height/bSize    %blocks vertically
    for j=1:vid_width/bSize      %blocks horizontally
        for d=1:length(shifts)           %motion_v displacements
            
            % "inside-padding" rows and columns indices of the blocks
            rows = 1+dy_max+(i-1)*bSize : 1+dy_max+(i-1)*bSize + bSize-1;
            cols = 1+dx_max+(j-1)*bSize : 1+dx_max+(j-1)*bSize + bSize-1;
            
            shifted_rows = rows + shifts(1,d);  % rows + dy
            shifted_cols = cols + shifts(2,d);  % rows + dx
            
            ref_rows = rows-dy_max;
            ref_cols = cols-dx_max;
            
            %compute MSE for all of them
            Diff2 = (F(shifted_rows,shifted_cols)-nextF(ref_rows,ref_cols)).^2;
            err(i,j,d) = sum(Diff2(:))/numel(Diff2);
        end
        %Find minimum MSEs and keep only the associated motion vectors
        %Unwarap matrix of blocks reading the blocks horizontally
        temp = find(err(i,j,:) == min(err(i,j,:)));
        if length(temp) > 1  % if more than one vector has minimum MSE
            temp = temp(1,1);  % keep the first one
        end
        mVecsIndexes(1,index) = temp;
        mVecs(:,index) = shifts(:,temp);
        index = index+1;
    end
end

end

