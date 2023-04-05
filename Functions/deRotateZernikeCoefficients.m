function varargout = deRotateZernikeCoefficients(a)
% Convert modal coefficients from calculateZernikeCoefficients into
% de-Rotated coefficients by taking the quadrature sum of coefficients
% correseponding to the same magnitude degree n, and order m, with flipped
% signs: (n,m) and (n,-m). See Methods for details.  
%
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023
%
% Inputs
%       a: modal coefficients computed by calculateZernikeCoefficients.m
% Outputs
%       C: 1d array of de-Rotated coefficients 
%       Nr: 1xlength(C) array specifying the radiial degree n of each
%       element in C (optional)
%       Mr: 1xlength(C) array specificying the azimuthal order m of each
%       element in C (optional)



%check that coefficient array length is valid (must be a triangular number)
T = length(a);
is_tri = @(x) floor(sqrt(8*x+1))==sqrt(8*x+1); 
test = is_tri(T);
if ~test
    error('error: length of coefficient array  is not triangular, unable to infer length')
end


% Infer radial degree from length of coefficients
% T = n(n+1)/2
Z_degree = (-1 + sqrt(1+8*T))/2;
N = []; M = [];
for n = 0:Z_degree-1
    N = [N n*ones(1,n+1)];
    M = [M -n:2:n];
end

C = [];
Mr = [];
Nr = [];

for j = 1:length(M)
    if M(j) > 0
        c1 = a(j);
        c2 = a(j-M(j));
        Cnew = sqrt(c1^2 + c2^2);       %quadrature of rotationally related coefficients
        C = [C Cnew];
        Mr = [Mr M(j)];
        Nr = [Nr N(j)];
    elseif M(j) == 0
        Cnew = a(j);
        C = [C Cnew];
        Mr = [Mr M(j)];
        Nr = [Nr N(j)];
    end
end

varargout{1} = C;
varargout{2} = Nr;
varargout{3} = Mr;




end

