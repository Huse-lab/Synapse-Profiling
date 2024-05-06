% Plot Zernike charts of arbitrary sizes
% MDJ: 2022-06-14

%% Ask user for how many Zernike functions to visualize

k_max = input("How many Zernike functions to plot? ",'s');
k_max = str2num(k_max);
k_max = uint16((k_max));
if ~isa(k_max,'integer'); error("Please input an integer. "); end
if size(k_max,1) == 0; error("Please input an integer. "); end

%% Create the Zernike array

n_array = [0]; m_array = [0];

n = 1;
idx = 1;
while idx <= k_max
    m = -n;
    while m ~= n && idx <= k_max
        idx = idx+1;
        n_array = [n_array n];
        m_array = [m_array m];
        m = m+2;
    end    
    if m == n && idx <= k_max
        idx = idx+1;
        n_array = [n_array n];
        m_array = [m_array m];
        n = n+1;
    end
end

%% Calculate and plot Zernike functions

plotmax = 64; % Set the maximum # of Zernike plots per graph.
fig_n = 1 + floor(double(k_max)/plotmax); 
rem = double(mod(k_max,plotmax)); % For setting rows and columns of figures.

x = -1:0.01:1;
        [X,Y] = meshgrid(x,x);
        [theta,r] = cart2pol(X,Y);
        j = r<=1;
        z = nan(size(X));
        
for i=1:fig_n
    figure
    % Set subplot as square
    if i<fig_n % For full 36-panel figures
        n_row = ceil(sqrt(plotmax));
        n_col = ceil(plotmax/n_row);
    elseif i==fig_n % For the remainders in the last figure
        n_row = ceil(sqrt(rem));
        n_col = ceil(rem/n_row);
    end
    
    for k=1:plotmax
        idx_n = n_array(k+plotmax*(i-1));
        idx_m = m_array(k+plotmax*(i-1));        
    z(j) = zernfun(idx_n,idx_m,r(j),theta(j));
    subplot(n_row,n_col,k)
    pcolor(x,x,z), shading interp
    pbaspect([1 1 1]); axis off;
    title(strcat(num2str(idx_n),',',num2str(idx_m)))
%   colorcet('L8')
%   colormap(brewermap([],'PRGn'))
    colormap(turquoisebrown)
    end
    
end