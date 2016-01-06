% Project 3
clear

% Function definitions
uniform_quantizer = @(x,ssize) round(x/ssize)*ssize;
PSNR = @(D) 10*log10(255^2./D);
lagrangian = @(n,q,r) n+0.2*(q^2)*r;

% Set globals
video_width = 176;   %must be divisible by bSize
video_height = 144;  %must be divisible by bSize
FPS = 30;
Nframes = 50;  %to be imported
% Block Size for motion compensation
bSize = 16;
% Size of search space for the shifts. (same for negative and positive)
dy_max = 10;
dx_max = 10;


% Compute all the possible shifts. For +- 10pxls there are 
% 441 possible motion vectors - 8bits needed without entropy coding.
index = 1;
n_vecs = (2*dy_max+1)*(2*dx_max+1);
shifts = zeros(2,n_vecs);
for i=-dy_max:dy_max
    for j=-dx_max:dx_max
        shifts(:,index) = [i,j];
        index = index+1;
    end
end

% import 50 frames of video
V = yuv_import_y('foreman_qcif.yuv',[video_width video_height],Nframes);

% deconstruct cell array
Frames = zeros(video_height,video_width,Nframes);
FramesPadded = zeros(video_height+2*dy_max,...
                     video_width+2*dx_max,Nframes);
for f=1:Nframes
    Frames(:,:,f) = V{f,1};
    % Pad with zeros around to handle the borders in the motion_vec_search
    FramesPadded(:,:,f) = padarray(V{f,1},[dy_max dx_max]);
end



% Quantize with different step sizes
q_step = 2.^[3:6];
Framesq1 = zeros(video_height,video_width,Nframes,length(q_step));
Reconstructed1 = zeros(video_height,video_width,Nframes,length(q_step));
Err1 = zeros(Nframes,length(q_step));
PSNR1 = zeros(Nframes,length(q_step));
Rate1 = zeros(Nframes,length(q_step));
for f=1:Nframes
    %Intra frame only (mode 1)
    for q=1:length(q_step)
        Framesq1(:,:,f,q) = ...
            uniform_quantizer(blockproc(Frames(:,:,f),[8 8],@dctz2),q_step(q));
        Reconstructed1(:,:,f,q) = ...
            blockproc(Framesq1(:,:,f,q),[8 8],@idctz2);
        %Compute MSEs and PSNRs for mode 1
        Diff2 = (Reconstructed1(:,:,f,q) - Frames(:,:,f)).^2;
        Err1(f,q) = sum(Diff2(:))/numel(Diff2(:));
        PSNR1(f,q) = 10*log10( (255^2)/Err1(f,q) );
        %Compute Kbits/Frame for mode 1
        Rate1(f,q) = EntropyRate(Framesq1(:,:,f,q),q_step(q))*...
            video_width*video_height/1000;
    end
end



% Pre allocate memory
MotVecs = zeros(2,99,Nframes-1);
MotVecsIndices = zeros(1,99,Nframes-1);
% Use motion vectors to predict next frame and compute difference image
Predicted3 = zeros(video_height,video_width,Nframes); %mode 3
Residual = zeros(video_height,video_width,Nframes);
Residualq = zeros(video_height,video_width,Nframes,length(q_step));
ResidualDCT = zeros(video_height,video_width,Nframes,length(q_step));
MotionCompq = zeros(video_height,video_width,Nframes,length(q_step));
%First frame cannot be predicted from past
Predicted3(:,:,1) = Frames(:,:,1); 
for q=1:length(q_step)
    MotionCompq(:,:,1,q) = Frames(:,:,1);
end

Err2 = zeros(Nframes,length(q_step));
PSNR2 = zeros(Nframes,length(q_step));
Rate2 = zeros(Nframes,length(q_step));
Err3 = zeros(Nframes,length(q_step));
PSNR3 = zeros(Nframes,length(q_step));
Rate3 = zeros(Nframes,length(q_step));


for f=1:Nframes-1   %at every frame look at the next one to decide mode
    
    % Motion Compensation encoder (mode 3)
    % Compute motion vectors for all frames
    % Exaustive search implemented in ComputeMotVecs.m
    [MotVecs(:,:,f),MotVecsIndices(:,:,f)] = ...
        ComputeMotVecs(FramesPadded(:,:,f),Frames(:,:,f+1),bSize,shifts);
    Predicted3(:,:,f+1) = PredictFrame(FramesPadded(:,:,f),...
                                MotVecs(:,:,f),bSize,[dy_max dx_max]);
    Residual(:,:,f+1) =  Frames(:,:,f+1) - Predicted3(:,:,f+1);
    
    % Estimate average bitrate needed for all motion vectors
    p = hist(MotVecsIndices(:,:,f),n_vecs)./n_vecs;
    H = -sum(p.*log2(p+eps)); %bits/block
    rate_vecs = H*size(MotVecs,2)/1000;   %Kbits/frame
    
    %quantize residuals, with q_step in 4th dimension
    for q=1:length(q_step)
        ResidualDCT(:,:,f+1,q) = blockproc(Residual(:,:,f+1),[8 8],@dctz2);
        Residualq(:,:,f+1,q) = ...
            uniform_quantizer(ResidualDCT(:,:,f+1,q),q_step(q));
        % Reconstruct with motion vector and quantized residuals
        MotionCompq(:,:,f+1,q) = Predicted3(:,:,f+1) + Residualq(:,:,f+1,q);
    
        % Compute Error,PSNR for mode 2 (next - current quantized frame)
        Diff2 = (Frames(:,:,f+1)-Reconstructed1(:,:,f,q)).^2;
        Err2(f+1,q) = sum(Diff2(:))/numel(Diff2(:));
        PSNR2(f+1,q) = 10*log10( (255^2)/Err2(f+1,q) );
        % Rate2 for comparisons can be kept to zero, since the 2 bits to
        % signal the mode are needed for all 3 modes
        
        % Compute Error,PSNR for mode 3
        Diff2 = (MotionCompq(:,:,f+1,q)-Frames(:,:,f+1)).^2;
        Err3(f+1,q) = sum(Diff2(:))/numel(Diff2(:));
        PSNR3(f+1,q) = 10*log10( (255^2)/Err3(f+1,q) );
        
        % Compute Bitrate for mode 3
        rate_residual = EntropyRate(ResidualDCT(:,:,f+1,q),q_step(q))*...
                        video_height*video_width/1000;
        Rate3(f+1,q) = rate_vecs + rate_residual;  %Kbits/frame mode 3
    end
end


% Play Videos
%implay(uint8(Frames),FPS); 
%implay(uint8(Residual),FPS);  
%implay(uint8(Predicted3),FPS);

% Prepare Videos for writing
% Framesw(:,:,1,:) = Frames;
% Residualw(:,:,1,:) = Residual;
% Predictedw(:,:,1,:) = Predicted3;
% MotionCompqw(:,:,1,:,:) = MotionCompq;
% % Write Videos
% v = VideoWriter('Results\Quantized residual.avi','Grayscale AVI');
% open(v);
% writeVideo(v,uint8(MotionCompqw(:,:,:,:,4)));
% close(v);

%%%%%%=======================================================%%%%


% % Use motion vectors to predict next frame and compute difference image
% Predicted3 = zeros(video_height,video_width,Nframes);
% Predicted3(:,:,1) = Frames(:,:,1); %First frame cannot be predicted from past
% Residual = zeros(video_height,video_width,Nframes);
% for f=1:Nframes-1
%     Predicted3(:,:,f+1) = PredictFrame(FramesPadded(:,:,f),...
%                                 MotVecs(:,:,f),bSize,[dy_max dx_max]);
%     Residual(:,:,f+1) =  Frames(:,:,f+1) - Predicted3(:,:,f+1);
% end

% %loop over frames
% for k=1:length(V)
%     %simmetrically extend around to handle edges for motion compensation
%     V_pad{k} = padarray(V{k},[10 10],'symmetric');
% end
% 
% %construct all possible motion vectors
% % There are 441 possible motion vectors --> 8 bits per vector needed.
% index = 1;
% for i=-10:10
%     for j=-10:10
%         motionv(index,:) = [i,j];
%         index = index+1;
%     end
% end

% % compute prediction error (mse) for all frames, blocks and motion vectors
% for k=1:length(V_pad)-1 %frames
%     for i=1:video_height/16    %rows
%         for j=1:video_width/16      %columns
%             for d=1:length(motionv)           %motionv displacements
%                 
%                 err(k,i,j,d) = ...
%                     mse(V_pad{k}(11+(i-1)*16+motionv(d,1):11+(i-1)*16+motionv(d,1)+15,...
%                     11+(j-1)*16+motionv(d,2):11+(j-1)*16+motionv(d,2)+15),...
%                     V_pad{k+1}(11+(i-1)*16:11+(i-1)*16+15,...
%                     11+(j-1)*16:11+(j-1)*16+15));
%                 
%             end
%             
%         end
%     end
% end
% 
% %Find minimum mses and keep only the associated motion vector indexes
% for k=1:length(V_pad)-1 %frames
%     for i=1:video_height/16    %rows
%         for j=1:video_width/16 
%     
%             temp = find(err(k,i,j,:)==min(err(k,i,j,:)));
%             if length(temp) > 1  % if more than one vector has minimum MSE
%                temp = temp(1,1);  % keep the first one
%             end
%             indexes(k,i,j) = temp;
%         
%         end
%     end
% end
% 
% %unwrap the indexes to compute a list of motion vectors for every frame
% for k=1:length(V_pad)-1
%     index = 1;
%     for i=1:video_height/16
%         for j=1:video_width/16
%             tempmv{k,index} = motionv(indexes(k,i,j),:);
%             index = index+1;
%         end
%     end
% end
% 
% for k=1:length(V_pad)-1  % frames
%     temp = [];
%     for i=1:size(tempmv,2)   %all blocks
%         temp = [temp;tempmv{k,i}];
%     end
%     MotionVectors{k} = temp';
% end



% %split image into blocks
% for i=1:length(V)
%     V{i,1} = mat2cell(V{i,1},16*ones(1,size(i,1)/8),16*ones(1,size(i,2)/8),3);
% end
%
% % Quantize with different step sizes
% q_step = 2.^[3:6];
%
% for j=1:length(q_step)
%     s = q_step(j);
%     for i=1:length(V)
%         Vdctq{i,j} = uniform_quantizer(Vdct{i,1},s);
%     end
% end
%
% % Compute MSE for every frame (along rows) and q_step (along columns)
% for j=1:length(q_step)
%     for i=1:length(V)
%         err(i,j) = mse(Vdct{i,1},Vdctq{i,j});
%     end
% end
% % Reconstruct and play quantized videos
% for j=1:length(q_step)
%     for i=1:length(V)
%         V_rec{i,j} = blockproc(Vdctq{i,j},[8 8],@idctz2);
%     end
% end
% %Vid[1,2] are images
% %Vid[3] are frames
% %Vid[4] are different quantizer steps
% for j=1:length(q_step)
%     for i=1:length(V)
%         Vid(:,:,i,j) = V_rec{i,j};
%     end
% end
%
% % play the best one and the worst one
% implay(uint8(Vid(:,:,:,1)),FPS);
% implay(uint8(Vid(:,:,:,4)),FPS);
% %
%
% %calculate mean PSNRs
% meanPSNRquality = mean(PSNR(err));
%
% rates = zeros(size(q_step)); %bitrate placeholder array
%
% %calculate bitrates/s per each qunatizer step
% for j=1:length(q_step) %for each quantizer step
%     s = q_step(j);
%     coefs = zeros(8,8,(video_width/8)*(video_height/8)*length(Vdctq));
%
%     for i=1:length(Vdctq) %for each frame
%         index = 1;
%         for w=1:(video_width/8)
%             for h=1:(video_height/8)
%                 for ww=1:8
%                     for hh=1:8
%                         coefs(ww,hh,(i-1)*(video_width/8)*(video_height/8)+index) = Vdctq{i,j}(8*(h-1)+hh,8*(w-1)+ww);
%                     end
%                 end
%                 index = index + 1;
%             end
%         end
%     end
%
%     %now calculate R
%     H=zeros(8,8);
%     for w=1:8
%         for h=1:8
%             vals = squeeze(coefs(w,h,:));
%             p = hist(vals,min(vals):s:max(vals));
%             p = p/sum(p);
%
%             H(w,h) = -sum(p.*log2(p+eps)); %eps added for log of 0 vals
%         end
%     end
%
%     rates(j) = mean2(H);
% end
%
% %covert to kbit/s
% ratesKBPS = rates .* ((video_height*video_width*FPS)/1000);
%
% figure;
% plot(ratesKBPS, meanPSNRquality, '+-', 'linewidth', 2);
% title('Performance vs bitrate (DCT)');
% grid on;
% xlabel('[Kbps] rate');
% ylabel('[dB] PSNR');
