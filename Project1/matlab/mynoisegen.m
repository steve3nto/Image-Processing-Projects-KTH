% Generates a noise m-by-n image
% USAGE: noise = mynoisegen(type, m, n, param1, param2)
% 
% type:         
% 'uniform'     uniform noise between param1 and param2
%               default param1 = -1, param2 = 1
% 'gaussian'    gaussian noise of mean param1 and variance param2
%               default param1 = 0, param2 = 1
% 'saltpepper'  salt & pepper noise)
%               p(0) = param1, p(1) = param2
function noise = mynoisegen(type, m, n, param1, param2)
switch lower(type)
case 'uniform'
    if nargin <= 3
        param1 = -1;
        param2 = 1;
    end
    noise = param1 + (param2-param1).*rand(m,n);
case 'gaussian'
    if nargin <= 3
        param1 = 0;
        param2 = 1;
    end
    noise = param1 + sqrt(param2).*randn(m,n);
case 'saltpepper'
    if nargin <= 3
        param1 = 0.5;
        param2 = 1;
    end
    noise = ones(m,n) .* .5;
    
    nn = rand(m,n);
    noise(nn<=param1) = 0;
    noise(nn>param1 & nn<=param1+param2) = 1;
otherwise
    error('Unknown noise type.');
end