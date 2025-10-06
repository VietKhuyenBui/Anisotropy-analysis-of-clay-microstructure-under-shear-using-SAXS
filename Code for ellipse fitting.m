% This 2D fabric orientation visualization technique is designed for 
% various 2D datasets. Here, the raw data from 2D SAXS patterns is used as an example.
% For more detailed information, please refer to the following journal paper:
% Bui, V. K., et al., "Anisotropy analysis of clay microstructure under
% shear using SAXS."
% If you use this technique for academic purposes, please cite the paper mentioned above.
%
% The patch size is determined by 'patchw', which corresponds to 2*patchw. 
% The patch size can also be an odd number.
% For larger values of patchw, it is recommended to increase the value of N
% accordingly.

% Read the SAXS Image
filename = 'KAO - LI=3 10s-1.dat';
data = readmatrix(filename);
x = data(:, 1);
y = data(:, 2);
I = data(:, 3);
in_roi = (x >= -2 & x <= 2 & y >= -2 & y <= 2);  
x = x(in_roi);  y = y(in_roi);  I = I(in_roi);

% Symmetry line y = 0
pivot_y = 0;
% Build the combined point set: Replace that lower strip with a mirrored copy of the upper strip (|x|<1 & y>0)
upper_src   = (x > -0.1 & x <  0.1 & y > pivot_y);   
lower_strip = (x > -0.1 & x <  0.1 & y < pivot_y);   
keep_mask = ~lower_strip;
x_keep = x(keep_mask);
y_keep = y(keep_mask);
I_keep = I(keep_mask);

% Mirror the upper strip into the lower strip: (x,y)->(-x,-y)
x_mirror = -x(upper_src)+0.02;
y_mirror = -y(upper_src);
I_mirror =  I(upper_src);

% Combined point
x_combined = [x_keep;   x_mirror];
y_combined = [y_keep;   y_mirror];
I_combined = [I_keep;   I_mirror];

% Interpolate to a regular grid for display/processing
nx = 300; ny = 300;
xg = linspace(min(x_combined), max(x_combined), nx);
yg = linspace(min(y_combined), max(y_combined), ny);
[Xq, Yq] = meshgrid(xg, yg);

F = scatteredInterpolant(x_combined, y_combined, I_combined, 'natural', 'none');
V = F(Xq, Yq);
V(isnan(V)) = 0;              
V(V <= 0)   = eps;             
V_log = log(V);

% Plot the 2D SAXS map and the raw points
figure; 
imagesc('XData', xg, 'YData', yg, 'CData', V_log);
set(gca, 'YDir', 'normal'); axis image;
colorbar; xlabel('X'); ylabel('Y'); hold on;
plot(x_combined, y_combined, 'b.', 'MarkerSize', 0.1);
xlim([-2 2]); ylim([-2 2]);

% Q-tensor coefficients via Monte Carlo weights (matches 300x300 grid)
patchw = 150;                            
n_nbmaskQ = zeros(patchw*2);
n_maskQ   = zeros(patchw*2, patchw*2, 2, 2);

N = 10000;
Niter = floor(patchw*2 * patchw*2 * 200 / N);
disp(Niter);

for j = 1:Niter
    r = (rand(N, 2) - 0.5) * 2 * patchw;     
    K = zeros(N,2);
    den = sqrt(r(:,1).^2 + r(:,2).^2);
    den(den == 0) = 1;
    K(:,1) = r(:,1) ./ den;
    K(:,2) = r(:,2) ./ den;
    xpix = floor(r + patchw) + 1;             
    for i = 1:N
        rr = xpix(i,2); cc = xpix(i,1);
        if rr >= 1 && rr <= patchw*2 && cc >= 1 && cc <= patchw*2
            n_nbmaskQ(rr, cc) = n_nbmaskQ(rr, cc) + 1;
            n_maskQ(rr, cc, 1, 1) = n_maskQ(rr, cc, 1, 1) + K(i,1)*K(i,1);
            n_maskQ(rr, cc, 1, 2) = n_maskQ(rr, cc, 1, 2) + K(i,1)*K(i,2);
            n_maskQ(rr, cc, 2, 1) = n_maskQ(rr, cc, 2, 1) + K(i,2)*K(i,1);
            n_maskQ(rr, cc, 2, 2) = n_maskQ(rr, cc, 2, 2) + K(i,2)*K(i,2);
        end
    end
end
n_maskQ(:,:,1,1) = n_maskQ(:,:,1,1) ./ max(n_nbmaskQ, 1);
n_maskQ(:,:,1,2) = n_maskQ(:,:,1,2) ./ max(n_nbmaskQ, 1);
n_maskQ(:,:,2,1) = n_maskQ(:,:,2,1) ./ max(n_nbmaskQ, 1);
n_maskQ(:,:,2,2) = n_maskQ(:,:,2,2) ./ max(n_nbmaskQ, 1);

% Build Q from weights and intensity field V (300x300)
Q = zeros(patchw*2, patchw*2, 2, 2);
Q(:,:,1,1) = n_maskQ(:,:,1,1) .* V;
Q(:,:,1,2) = n_maskQ(:,:,1,2) .* V;
Q(:,:,2,1) = n_maskQ(:,:,2,1) .* V;
Q(:,:,2,2) = n_maskQ(:,:,2,2) .* V;

Q2 = permute(sum(sum(Q,1),2), [3,4,1,2]);   
Q2(1,1) = Q2(1,1) - 1/2;                    
Q2(2,2) = Q2(2,2) - 1/2;
Q22 = sqrt(2) * Q2;
Q22 = Q22 / sum(V(:));                       

% Ellipse from Q principal axes
[eigVecs, eigVals] = eig(Q22);
evals = diag(eigVals);
[~, imax] = max(evals);

majorAxisLength = sqrt(max(evals));
minorAxisLength = sqrt(min(evals));

theta0 = atan2(eigVecs(2, imax), eigVecs(1, imax));
theta  = linspace(0, 2*pi, 180);
ex = majorAxisLength * cos(theta);
ey = minorAxisLength * sin(theta);

R = [cos(theta0) -sin(theta0); sin(theta0) cos(theta0)];
E  = R * [ex; ey];

cx = mean(x_combined);  cy = mean(y_combined);
E(1,:) = E(1,:) + cx;   E(2,:) = E(2,:) + cy;

plot(E(1,:), E(2,:), 'r-', 'LineWidth', 2);  
axis equal; hold off;

fprintf('Major Axis Length: %f\n', majorAxisLength);
fprintf('Minor Axis Length: %f\n', minorAxisLength);
fprintf('Orientation (deg): %f\n', rad2deg(theta0));
