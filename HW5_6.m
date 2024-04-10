clear all;
close all;

h = 0.1; %dx and dy
dt = 1; 
N = 4; %num vertices of polygon
eps = 0.1;
a = 1; %radius of circumscribed circle of the polygon
gamma = [1,-1,1,-1]; %strength of vortex
x = -11:h:11;
y = -11:h:11;
t = 0:dt:100;
%x0 = (a-.2)*cos(2*pi/N); y0 = (a-.2)*sin(2*pi/N);
x0 = -.5; y0 = .2;
z0 = [1+10i, -1+10i, -1-10i, 1-10i]; % Initial positions of vortices

for k = 1:length(t)
    U{k} = zeros(length(x), length(y));
    V{k} = zeros(length(x), length(y));
end
[pvx,pvy,vort_vel] = vortex_pos(N,gamma,z0,t,a);
for k = 1:length(t)
    for m = 1:N
        for i = 1:length(x)
            for j = 1:length(y)
                xm = pvx(:,k); ym = pvy(:,k);
                dx(1) = x(i)-xm(1); dy(1) = y(j)-ym(1);
                dx(2) = x(i)-xm(2); dy(2) = y(j)-ym(2);
                dx(3) = x(i)-xm(3); dy(3) = y(j)-ym(3);
                dx(4) = x(i)-xm(4); dy(4) = y(j)-ym(4);

                % Exclude points too close to the vortex center
                if any(sqrt(dx(:).^2 + dy(:).^2) > 0.01) 
                    % Compute velocity
                    U{k}(i,j) = U{k}(i,j) + (gamma(m)./(2*pi))*(ym(m)-y(j))./((x(i)-xm(m)).^2+(y(j)-ym(m)).^2);
                    V{k}(i,j) = V{k}(i,j) + (gamma(m)/(2*pi))*(xm(m)-x(j))./((x(i)-xm(m)).^2+(y(j)-ym(m)).^2);
                    % Check for excessively high velocity and set to zero
                    exceed_limit_indices = abs(U{k}.^2 + V{k}.^2) > 10;
                    U{k}(exceed_limit_indices) = 0;
                    V{k}(exceed_limit_indices) = 0;
                    
                    % Replace excessively high velocities with the average of surrounding 10 points
                    [rows, cols] = find(exceed_limit_indices);
                    for p = 1:length(rows)
                        i = rows(p);
                        j = cols(p);
                        
                        % Define region around the point
                        i_start = max(i-5, 1);
                        i_end = min(i+5, length(x));
                        j_start = max(j-5, 1);
                        j_end = min(j+5, length(y));
                        
                        % Exclude the point itself from averaging
                        point_index = (i - i_start) * (i_end - i_start + 1) + (j - j_start + 1);
                        % Exclude the point itself from averaging
                        surrounding_U(i - i_start + 1, j - j_start + 1) = NaN;
                        surrounding_V(i - i_start + 1, j - j_start + 1) = NaN;
                        
                        % Compute average of surrounding velocities
                        avg_U = nanmean(surrounding_U(:));
                        avg_V = nanmean(surrounding_V(:));
                        % Replace NaN values with the computed average
                        surrounding_U(isnan(surrounding_U)) = avg_U;
                        surrounding_V(isnan(surrounding_V)) = avg_V;
                        
                        % Replace excessively high velocity with average
                        U{k}(i, j) = avg_U;
                        V{k}(i, j) = avg_V;
                    end
                else
                    % Check for excessively high velocity and set to zero
                    exceed_limit_indices = abs(U{k}.^2 + V{k}.^2) > 10;
                    U{k}(exceed_limit_indices) = 0;
                    V{k}(exceed_limit_indices) = 0;
                    
                    % Replace excessively high velocities with the average of surrounding 10 points
                    [rows, cols] = find(exceed_limit_indices);
                    for p = 1:length(rows)
                        i = rows(p);
                        j = cols(p);
                        
                        % Define region around the point
                        i_start = max(i-5, 1);
                        i_end = min(i+5, length(x));
                        j_start = max(j-5, 1);
                        j_end = min(j+5, length(y));
                        
                        % Exclude the point itself from averaging
                        point_index = (i - i_start) * (i_end - i_start + 1) + (j - j_start + 1);
                        % Exclude the point itself from averaging
                        surrounding_U(i - i_start + 1, j - j_start + 1) = NaN;
                        surrounding_V(i - i_start + 1, j - j_start + 1) = NaN;
                        
                        % Compute average of surrounding velocities
                        avg_U = nanmean(surrounding_U(:));
                        avg_V = nanmean(surrounding_V(:));
                        % Replace NaN values with the computed average
                        surrounding_U(isnan(surrounding_U)) = avg_U;
                        surrounding_V(isnan(surrounding_V)) = avg_V;
                        % Replace excessively high velocity with average
                        U{k}(i, j) = avg_U;
                        V{k}(i, j) = avg_V;
                    end
                end
            end
        end
    end
end

for k = 1:length(t)
    %[x_particle, y_particle] = particle_plot(x0, y0, U{k}, V{k}, t, dt, x, y);
    
    % Create contour plot
    figure(k);
    contour(x, y, abs(U{k}.^2 + sqrt(-1)*V{k}.^2));
    hold on
    
    % Plot particle position
    %scatter(x_particle, y_particle, 'b', 'filled');
    
    % Plot fixed points
    % for m = 1:N
    %     scatter(a*cos(2*pi*m/N), a*sin(2*pi*m/N), 'ro', 'filled');
    % end
    
    hold off
    title(['t = ', num2str(t(k))])
    xlabel('X')
    ylabel('Y')
    zlabel('Magnitude of Velocity Field')
    filename = sprintf('magnitude_of_velocity_field_for_N_%d_a_%.2f_gamma_%.2f_t=_%.2f.png', N, a, gamma(1),t(k));
    saveas(gcf, filename);
    close; % Close figure to avoid displaying multiple figures
end

%Compute the vorticity
vorticity = cell(size(U));
angular_velocity = zeros(length(t), 1);

for k = 1:length(t)
    [dV_dx, dU_dy] = gradient(V{k}, h, h);
    vorticity{k} = dV_dx - dU_dy;

    % Compute the magnitude of the vorticity as the angular velocity
    angular_velocity(k) = mean(abs(vorticity{k}(:))); % Taking the average for all grid points
    angular_velocity_analytic(k) = 3*gamma(1)*(N-1)/(4*pi*a^2);
end

% Plot the angular velocity
figure;
hold on
plot(t, angular_velocity)
plot(t, angular_velocity_analytic, '--r')
hold off
xlabel('Time');
ylabel('Angular Velocity');
title('Angular Velocity vs Time');

% Annotate the plot with parameters
parameters_text = {sprintf('a = %.2f', a), sprintf('gamma(1) = %.2f', gamma(1)), sprintf('N = %d', N)};
text(0.7, 0.9, parameters_text, 'Units', 'normalized');

% Save the plot as a PNG file
filename = sprintf('angular_velocity_plot_a_%.2f_gamma_%.2f_N_%d.png', a, gamma(1), N);
saveas(gcf, filename);

function [x_particle,y_particle] = particle_plot(x0,y0,U,V,t,dt,x,y)
    x_particle = zeros(length(t), 1); % x-coordinate of particle
    y_particle = zeros(length(t), 1); % y-coordinate of particle

    for k = 1:length(t)
        % RK4 integration for particle position
        if k == 1
            x_particle(k) = x0;
            y_particle(k) = y0;
        else
            k1 = dt * interp2(x, y, U, x_particle(k-1), y_particle(k-1));
            l1 = dt * interp2(x, y, V, x_particle(k-1), y_particle(k-1));
            k2 = dt * interp2(x, y, U, x_particle(k-1) + k1/2, y_particle(k-1) + l1/2);
            l2 = dt * interp2(x, y, V, x_particle(k-1) + k1/2, y_particle(k-1) + l1/2);
            k3 = dt * interp2(x, y, U, x_particle(k-1) + k2/2, y_particle(k-1) + l2/2);
            l3 = dt * interp2(x, y, V, x_particle(k-1) + k2/2, y_particle(k-1) + l2/2);
            k4 = dt * interp2(x, y, U, x_particle(k-1) + k3, y_particle(k-1) + l3);
            l4 = dt * interp2(x, y, V, x_particle(k-1) + k3, y_particle(k-1) + l3);
            
            x_particle(k) = x_particle(k-1) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
            y_particle(k) = y_particle(k-1) + (1/6) * (l1 + 2*l2 + 2*l3 + l4);
        end
    end

end

function [pvx,pvy,v] = vortex_pos(N,gamma,z0,t,a)
    
    % Time parameters
    t0 = t(1); % Initial time
    tf = t(end); % Final time
    dt = t(2)-t(1); % Time step
    
    % Preallocate arrays to store positions and velocities
    z = zeros(N, tf/dt + 1);
    v = zeros(N, tf/dt + 1);
    
    % Set initial positions
    z(:, 1) = z0.';
    
    % Perform time integration using RK4
    for k = 1:(tf/dt)
        % Calculate velocities of each vortex pair
        for n = 1:N
            omega_n = 0;
            for m = 1:N
                if m ~= n
                    omega_n = a*(-1*(2*pi*n/N)*sin(2*pi*n*t(k)/N)+sqrt(-1)*(2*pi*n/N)*cos(2*pi*n*t(k)/N));
                end
            end
            v(n, k) = omega_n;
        end
        
        % Update positions using RK4
        k1 = dt * v(:, k);
        k2 = dt * (v(:, k) + 0.5 * k1);
        k3 = dt * (v(:, k) + 0.5 * k2);
        k4 = dt * (v(:, k) + k3);
        
        z(:, k+1) = z(:, k) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    pvx = real(z(:,:)); pvy = imag(z(:,:));
end