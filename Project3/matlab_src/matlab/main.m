% Project 3
clear

% Import luminance of the video
video_width = 176;
video_height = 144;
FPS = 30;
V = yuv_import_y('foreman_qcif.yuv',[video_width video_height],50);

% %convert to uint8
% for i=1:length(V)
%     V{i,1} = uint8(V{i,1});
% end

%Function definitions
uniform_quantizer = @(x,ssize) round(x/ssize)*ssize;
mse = @(x,y) sum(sum((y-x).^2))/(size(y,1) * size(y,2));
PSNR = @(D) 10*log10(255^2./D);
lagrangian = @(n,q,r) n+0.2*(q^2)*r;

%DCT-2
for i=1:length(V)
    Vdct{i,1} = blockproc(V{i,1},[8 8],@dctz2);
end

% Quantize with different step sizes
q_step = 2.^[3:6];

for j=1:length(q_step)
    s = q_step(j);
    for i=1:length(V)
        Vdctq{i,j} = uniform_quantizer(Vdct{i,1},s);
    end
end

% Compute MSE for every frame (along rows) and q_step (along columns)
for j=1:length(q_step)
    for i=1:length(V)
        err(i,j) = mse(Vdct{i,1},Vdctq{i,j});
    end
end
% Reconstruct and play quantized videos
for j=1:length(q_step)
    for i=1:length(V)
        V_rec{i,j} = blockproc(Vdctq{i,j},[8 8],@idctz2);
    end
end
%Vid[1,2] are images
%Vid[3] are frames
%Vid[4] are different quantizer steps
for j=1:length(q_step)
    for i=1:length(V)
        Vid(:,:,i,j) = V_rec{i,j};
    end
end

% play the best one and the worst one
implay(uint8(Vid(:,:,:,1)),30);
implay(uint8(Vid(:,:,:,4)),30);
% 

%calculate mean PSNRs
meanPSNRquality = mean(PSNR(err));

rates = zeros(size(q_step)); %bitrate placeholder array

%calculate bitrates/s per each qunatizer step
for j=1:length(q_step) %for each quantizer step
    s = q_step(j);
    coefs = zeros(8,8,(video_width/8)*(video_height/8)*length(Vdctq));
    
    for i=1:length(Vdctq) %for each frame
        index = 1;
        for w=1:(video_width/8)
            for h=1:(video_height/8)
                for ww=1:8
                    for hh=1:8
                        coefs(ww,hh,(i-1)*(video_width/8)*(video_height/8)+index) = Vdctq{i,j}(8*(h-1)+hh,8*(w-1)+ww);
                    end
                end
                index = index + 1;
            end
        end
    end
    
    %now calculate R
    H=zeros(8,8);
    for w=1:8
        for h=1:8
            vals = squeeze(coefs(w,h,:));
            p = hist(vals,min(vals):s:max(vals));
            p = p/sum(p);
            
            H(w,h) = -sum(p.*log2(p+eps)); %eps added for log of 0 vals
        end
    end
    
    rates(j) = mean2(H);
end

%covert to kbit/s
ratesKBPS = rates .* ((video_height*video_width*FPS)/1000);

figure;
plot(ratesKBPS, meanPSNRquality, '+-', 'linewidth', 2);
title('Performance vs bitrate (DCT)');
grid on;
xlabel('[Kbps] rate');
ylabel('[dB] PSNR');
