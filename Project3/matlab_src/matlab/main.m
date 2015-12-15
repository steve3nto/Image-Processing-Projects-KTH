% Project 3
clear

% Import luminance of the video
V = yuv_import_y('foreman_qcif.yuv',[176 144],50);

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

for j=1:length(q_step)
    for i=1:length(V)
        Vid(:,:,i,j) = V_rec{i,j};
    end
end
% play the best one
implay(uint8(Vid(:,:,:,1)),30)
implay(uint8(Vid(:,:,:,4)),30)
% 