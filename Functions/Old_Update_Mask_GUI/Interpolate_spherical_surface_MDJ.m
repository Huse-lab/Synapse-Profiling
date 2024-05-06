function[IPsignal] = Interpolate_spherical_surface(Qtheta,Qphi,Ktheta,Kphi,Ksignal,method,nneighbours)
% Find nearest neighbours and use nearest neighbour or linear interpolation
% to interpolate values on a spherical surface. This can be used, for
% example, to go from a equidistant to a rectangular grid

% Check inputs and otherwise use default values
if nargin < 6
    method = 'linear';
end

if nargin < 7 && strcmpi(method,'linear')
    nneighbours = 5;
end

% Check that inputs are column vectors
if size(Ktheta,1) == 1
    Ktheta = Ktheta';
    Kphi   = Kphi';
end
if size(Ksignal,1) == 1
    Ksignal = Ksignal';
end
% Check that query points are row vectors
if size(Qtheta,2) == 1
    Qtheta = Qtheta';
    Qphi   = Qphi';
end

% Calculate euclidian distances between points. This is much faster, 
nOGpoints  = numel(Ktheta);
nQpoints   = numel(Qtheta);
[X,Y,Z]    = sph2cart(Ktheta,Kphi,ones(nOGpoints,1));
[Qx,Qy,Qz] = sph2cart(Qtheta,Qphi,ones(1,nQpoints) );
Eudist = pdist2([X,Y,Z],[Qx;Qy;Qz]');

if strcmpi(method,'nearest')
       
    [~,minloc] = min(Eudist,[],2);   
    IPsignal   = Ksignal(minloc);
    
elseif strcmpi(method,'linear')

    [~,sortidx] = sort(Eudist);
    sel         = sortidx(1:nneighbours,:);
    pgcdist     = real(acos(sin(Kphi(sel)).*sin(Qphi)+cos(Kphi(sel)).*cos(Qphi).*cos(Ktheta(sel)-Qtheta)));
    weights     = 1./pgcdist;
    IPsignal    = nansum(weights.*(Ksignal(sel)))./nansum(weights);

else
    
    error('Method should be "nearest" or "linear"')
    
end
    
