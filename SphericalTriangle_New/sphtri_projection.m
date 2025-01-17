% Define vertices of the spherical triangle
A = [1, 0, 0];
B = [0, 1, 0];
C = [0, 0, 1];

% Normalize to get centroid at the north pole
centroid = (A + B + C) / norm(A + B + C);

% Create a sphere
theta = linspace(0, pi, 50);
phi = linspace(0, 2*pi, 50);
[theta, phi] = meshgrid(theta, phi);
x = sin(theta) .* cos(phi);
y = sin(theta) .* sin(phi);
z = cos(theta);

% Plot the sphere
figure;
surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
colormap([1 1 0]); % Yellow color for the sphere

% Plot the spherical triangle
fill3([A(1) B(1) C(1)], [A(2) B(2) C(2)], [A(3) B(3) C(3)], 'r', 'FaceAlpha', 0.6);
fill3([A(1) centroid(1) B(1)], [A(2) centroid(2) B(2)], [A(3) centroid(3) B(3)], 'b', 'FaceAlpha', 0.6);
fill3([B(1) centroid(1) C(1)], [B(2) centroid(2) C(2)], [B(3) centroid(3) C(3)], 'm', 'FaceAlpha', 0.6);
fill3([C(1) centroid(1) A(1)], [C(2) centroid(2) A(2)], [C(3) centroid(3) A(3)], 'c', 'FaceAlpha', 0.6);

% Plot the z-axis
plot3([0 0], [0 0], [-1 1], 'k', 'LineWidth', 2);

% Project the triangle onto the xy-plane
proj_A = A; proj_A(3) = 0;
proj_B = B; proj_B(3) = 0;
proj_C = C; proj_C(3) = 0;
proj_centroid = centroid; proj_centroid(3) = 0;

% Plot the projection
fill3([proj_A(1) proj_B(1) proj_C(1)], [proj_A(2) proj_B(2) proj_C(2)], [proj_A(3) proj_B(3) proj_C(3)], 'r', 'FaceAlpha', 0.6);
fill3([proj_A(1) proj_centroid(1) proj_B(1)], [proj_A(2) proj_centroid(2) proj_B(2)], [proj_A(3) proj_centroid(3) proj_B(3)], 'b', 'FaceAlpha', 0.6);
fill3([proj_B(1) proj_centroid(1) proj_C(1)], [proj_B(2) proj_centroid(2) proj_C(2)], [proj_B(3) proj_centroid(3) proj_C(3)], 'm', 'FaceAlpha', 0.6);
fill3([proj_C(1) proj_centroid(1) proj_A(1)], [proj_C(2) proj_centroid(2) proj_A(2)], [proj_C(3) proj_centroid(3) proj_A(3)], 'c', 'FaceAlpha', 0.6);

% Adjust view
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
view([1,1,0.5]);
grid on;
hold off;
