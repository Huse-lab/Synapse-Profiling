function varargout = calculateZernikeCoefficients(mc_i,max_degree)
% Decompose the mean curvature data on a 
% NxN gridded synapse on the Zernike Polynomial Basis.
%
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023
% 
% depends on zernfun.m from Paul Fricker, Mathworks
% 
% Inputs:
%       -mc_i: an NxN matrix defining the data of interest
%       -max_degree: integer, function will calculate coefficients based on
%       polynomial functions up to a specified radial degree (n). It will use
%       all polynomials up to and included each azimuthal order (m).
%
% Outputs:
%       -a: a 1xk array specifying modal coefficients that correspond to
%       the underlying data transformed onto Zernike basis. 
%       -Z: the zernike polynomials projected onto the NxN grid space used
%       to calculate a

if ~(mod(max_degree,1)==0)
    error('Error: max degree not properly defined as an integer')
end



% Initialize zernike functions on the NxN grid same size as mc_i
L = size(mc_i,1);
X = -1:2/(L-1):1;
[x,y] = meshgrid(X);
x = x(:); y = y(:);
[theta,r] = cart2pol(x,y); 
N = []; M = [];
for n = 0:max_degree
    N = [N n*ones(1,n+1)];
    M = [M -n:2:n];
end
circle_mask = ( r <= 1 );
Z = zernfun(N,M,r(circle_mask),theta(circle_mask));

% Compute modal coefficients with left-divisor
a = Z\mc_i(circle_mask);

varargout{1} = a;
varargout{2} = Z;


end