function varargout = calculateZernikeComplexity(mc_i,threshold)
%   given an NxN matrix containing data on a roughly circular disc,
%   calcualate the maximum radial degree of Zernike Polynomials required to
%   faithfully reconstruct
%
%
%
%
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023
% 
% Inputs:
%       -mc_i: an NxN matrix defining the data of interest
%
%
%
%
%

maximum = 30;
%Initialize grid for computing Zernike Polynomials
L = size(mc_i,1);
X = -1:2/(L-1):1;
[x,y] = meshgrid(X);
x = x(:); y = y(:);
[theta,r] = cart2pol(x,y); 

fits = [];
reconstructions = {};

for j = 1:maximum
    N = []; M = [];
    for n = 0:j
        N = [N n*ones(1,n+1)];
        M = [M -n:2:n];
    end
    is_in_circle = (r <= 1 );
    [a,Z] = calculateZernikeCoefficients(mc_i,j);
    mc_i(~is_in_circle) = 0;
    reconstructed = zeros(size(mc_i));
    reconstructed(is_in_circle) = Z*a;
    reconstructions = [reconstructions; reconstructed];
    fit = corr2(mc_i,reconstructed);
    fits = [fits; fit];
end
intersect = find(fits > threshold,1);

if isempty(intersect)
    intersect = max_order;
end

varargout{1} = intersect;
varargout{2} = fits;



end