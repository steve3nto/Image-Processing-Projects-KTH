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
% Different step sizes for Quantization of Block-DCT coefficients
q_step = 2.^[0:9];


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

FramesDCTq = zeros(video_height,video_width,Nframes,length(q_step));
Reconstructed1 = zeros(video_height,video_width,Nframes,length(q_step));
Err1 = zeros(Nframes,length(q_step));
PSNR1 = zeros(Nframes,length(q_step));
Rate1 = zeros(Nframes,length(q_step));
for f=1:Nframes
    %Intra frame only (mode 1)
    for q=1:length(q_step)
        FramesDCTq(:,:,f,q) = ...
            uniform_quantizer(blockproc(Frames(:,:,f),[8 8],@dctz2),q_step(q));
        Reconstructed1(:,:,f,q) = ...
            blockproc(FramesDCTq(:,:,f,q),[8 8],@idctz2);
        %Compute MSEs and PSNRs for mode 1
        Diff2 = (Reconstructed1(:,:,f,q) - Frames(:,:,f)).^2;
        Err1(f,q) = sum(Diff2(:))/numel(Diff2(:));
        PSNR1(f,q) = 10*log10( (255^2)/Err1(f,q) );
        %Compute bits/coeff for mode 1
        Rate1(f,q) = EntropyRate(FramesDCTq(:,:,f,q),q_step(q));
    end
end



% Pre allocate memory
MotVecs = zeros(2,99,Nframes-1);
MotVecsIndices = zeros(1,99,Nframes-1);
% Use motion vectors to predict next frame and compute difference image
Predicted3 = zeros(video_height,video_width,Nframes); %mode 3
Residual = zeros(video_height,video_width,Nframes);
ResidualDCT = zeros(video_height,video_width,Nframes);
Residualq = zeros(video_height,video_width,Nframes,length(q_step));
ResidualDCTq = zeros(video_height,video_width,Nframes,length(q_step));
MotionCompq = zeros(video_height,video_width,Nframes,length(q_step));
%First frame cannot be predicted from past
Predicted3(:,:,1) = Frames(:,:,1); 
for q=1:length(q_step)
    MotionCompq(:,:,1,q) = Frames(:,:,1);
end

Err2 = zeros(Nframes,length(q_step));
PSNR2 = zeros(Nframes,length(q_step));
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
    ResidualDCT(:,:,f+1) = blockproc(Residual(:,:,f+1),[8 8],@dctz2);
    
    % Estimate average bitrate needed for all motion vectors
    % motion vectors have usually low entropy, good for compression!
    p = hist(MotVecsIndices(:,:,f),n_vecs)./n_vecs;
    H = -sum(p.*log2(p+eps)); %bits/block
    rate_vecs = H/(bSize^2);   %bits/pixel
    
    %quantize residuals, with q_step in 4th dimension
    for q=1:length(q_step)
        ResidualDCTq(:,:,f+1,q) = ...
            uniform_quantizer(ResidualDCT(:,:,f+1),q_step(q));
        % Reconstruct with motion vector and quantized residuals
        Residualq(:,:,f+1,q) = blockproc(ResidualDCTq(:,:,f,q),[8 8],@idctz2);
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
        
        % Compute Bitrate for mode 3 in bits/pixel
        rate_residual = EntropyRate(ResidualDCTq(:,:,f+1,q),q_step(q));
                        
        Rate3(f+1,q) = rate_vecs + rate_residual;  %bits/pxl mode 3
    end
end

%Now we have rates and errors for all frames. To simplify average all the
%rates and set the rate as a constant in the optimization of the Lagrangian
%function to choose mode
AvgRate1 = mean(Rate1,1);
AvgRate3 = mean(Rate3,1);
% Play Videos
%implay(uint8(Frames),FPS); 
%implay(uint8(Residual),FPS);  
%implay(uint8(Predicted3),FPS);

%%%%%=====================================================%%%%%

%Compute avges for mode 1 only and plot performance curve
AvgPSNR1 = mean(PSNR1,1);
AvgRateKbps1 = mean(Rate1,1)*video_height*video_width*30/1000;

figure;
plot(fliplr(AvgRateKbps1), fliplr(AvgPSNR1), '+-', 'linewidth', 1);
title('Performance vs bitrate (BlockDCT inter-mode only)');
grid on;
xlabel('Rate [Kbps]');
ylabel('PSNR [dB]');

% Prepare Videos for writing
% Framesw(:,:,1,:) = Frames;   %original
% Residualw(:,:,1,:) = Residual;  %residual of mode 3
% Predictedw(:,:,1,:) = Predicted3; %mode 3 without summing residuals
% MotionCompqw(:,:,1,:,:) = MotionCompq;   %mode3 only
% Reconstructed1w(:,:,1,:,:) = Reconstructed1; %mode1 only
% % Write Videos
% v = VideoWriter('Results\OnlyMode1.avi','Grayscale AVI');
% open(v);
% writeVideo(v,uint8(Reconstructed1w(:,:,:,:,4)));
% close(v);

%%%%%%=======================================================%%%%

