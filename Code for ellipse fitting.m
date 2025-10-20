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

%% Data import and preprocessing
filename = 'Kao_LI=3 _ 0.1s-1.dat';
data = readmatrix(filename);
x = data(:, 1);
y = data(:, 2);
intensity = data(:, 3);
filtered_indices = (x >= -2 & x <= 2 & y >= 0 & y <= 2);
x_filtered = x(filtered_indices);
y_filtered = y(filtered_indices);
intensity_filtered = intensity(filtered_indices);
lowest_y_value = min(y_filtered);

%% Build the combined point set
positive_y_indices = (y_filtered > lowest_y_value);
x_rotated = -x_filtered(positive_y_indices);
y_rotated = -y_filtered(positive_y_indices) + 2 * lowest_y_value; 
intensity_positive = intensity_filtered(positive_y_indices);
x_combined = [x_filtered; x_rotated];
y_combined = [y_filtered; y_rotated];
intensity_combined = [intensity_filtered; intensity_positive];

%% Grid interpolation (q-space)
[Xq, Yq] = meshgrid(linspace(-2, 2, 300), linspace(-2, 2, 300));
[Xq, Yq] = meshgrid(linspace(min(x_combined), max(x_combined), 300), linspace(min(y_combined), max(y_combined), 300));
F = scatteredInterpolant(x_combined, y_combined, intensity_combined);
V = F(Xq, Yq);
V(isnan(V)) = 0;
V(V <= 0) = eps;
V_log = log10(V);

%% Visualize 2D SAXS pattern (log scale)
figure;
imagesc('XData', linspace(min(x_combined), max(x_combined), 300), ...
        'YData', linspace(min(y_combined), max(y_combined), 300), ...
        'CData', V_log);
set(gca, 'YDir', 'normal'); 
colorbar;
xlabel('X');
ylabel('Y');
hold on;
plot(x_combined, y_combined, 'b.', 'MarkerSize', 0.1); 
xlim([-2 2]);
ylim([-2 2]);
%% Monte carlo sampling: compute local Q tensor masks
patchw = 150;
n_nbmaskQ = zeros(patchw*2);
n_maskQ = zeros(patchw*2, patchw*2, 2, 2);
N = 10000;
Niter = floor(patchw*2 * patchw*2 * 200 / N); disp(Niter)
for j = 1:Niter
    r = (rand(N, 2) - 0.5) * 2 * patchw;
    K(:,1) = r(:,1) ./ (sqrt(r(:,1).^2 + r(:,2).^2));
    K(:,2) = r(:,2) ./ (sqrt(r(:,1).^2 + r(:,2).^2));
    x = floor(r + patchw) + 1;
    for i = 1:N
        n_nbmaskQ(x(i,2), x(i,1)) = n_nbmaskQ(x(i,2), x(i,1)) + 1;
        n_maskQ(x(i,2), x(i,1), 1, 1) = n_maskQ(x(i,2), x(i,1), 1, 1) + K(i,1)*K(i,1);
        n_maskQ(x(i,2), x(i,1), 1, 2) = n_maskQ(x(i,2), x(i,1), 1, 2) + K(i,1)*K(i,2);
        n_maskQ(x(i,2), x(i,1), 2, 1) = n_maskQ(x(i,2), x(i,1), 2, 1) + K(i,2)*K(i,1);
        n_maskQ(x(i,2), x(i,1), 2, 2) = n_maskQ(x(i,2), x(i,1), 2, 2) + K(i,2)*K(i,2);
    end
end
%% Normalize Q tensor masks and compute global Q tensor
n_maskQ(:,:,1,1) = n_maskQ(:,:,1,1) ./ n_nbmaskQ;
n_maskQ(:,:,1,2) = n_maskQ(:,:,1,2) ./ n_nbmaskQ;
n_maskQ(:,:,2,1) = n_maskQ(:,:,2,1) ./ n_nbmaskQ;
n_maskQ(:,:,2,2) = n_maskQ(:,:,2,2) ./ n_nbmaskQ;

Q(:,:,1,1) = n_maskQ(:,:,1,1) .* V;
Q(:,:,1,2) = n_maskQ(:,:,1,2) .* V;
Q(:,:,2,1) = n_maskQ(:,:,2,1) .* V;
Q(:,:,2,2) = n_maskQ(:,:,2,2) .* V;

Q2 = permute(sum(sum(Q,1),2), [3,4,1,2]);  
Qbar = Q2 / max(sum(V(:)), eps);             
Q22 = sqrt(2) * (Qbar - 0.5 * eye(2));    
T = Q22;                                     
[eigVecs, eigVals] = eig(Q22);
evals = diag(eigVals);
[~, imax] = max(evals);

%% Eigen decomposition and ellipse parameters
evals_Qbar = eig(Qbar);
majorAxisLength = sqrt(max(evals_Qbar));
minorAxisLength = sqrt(min(evals_Qbar));
a = majorAxisLength / max(minorAxisLength, eps);
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
fprintf('a: %f\n', a);
