function[Y] = dctz(X)
%Discrete Cosine Transform of 8x8 block
%   X input matrix
%   M block size
%   Y output matrix

% Making sure it works for both matrix inputs and structs
if(isstruct(X))
    M = X.blockSize(1);
    X = X.data;
else
    M = size(X);
    M = M(1);
end

if (M~=8)
    error('unexpected block size, should be 8x8');
    Y = []
    return
end

for k=0:M-1
    for i=0:M-1
        alpha = sqrt(((i~=0)+1)/M);
        c = cos(((2*k+1)*i*pi)/(2*M));
        A(i+1,k+1) = alpha * c; %stupid matlab non-zero indexing
    end
end

A = dctmtx(8);

Y = A*X*A';

end