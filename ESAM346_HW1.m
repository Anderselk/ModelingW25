function rk_Anders(N, T, x, y, g, p, q)

    % Glossary
    % N: Nsteps
    % T: Final time (ends at t=T)
    % (x,y): Vortex positions 
    % g: Vortex strengths 
    % (p,q): Particle positions 
    
    dt = T / N;  
    num_vortices = length(x);  % # of vortices
    num_particles = length(p);  % # of particles

    % Arrays to store trajectories
    x_traj = zeros(num_vortices, N+1);
    y_traj = zeros(num_vortices, N+1);
    p_traj = zeros(num_particles, N+1);
    q_traj = zeros(num_particles, N+1);

    % Init Conditions
    x_traj(:, 1) = x;
    y_traj(:, 1) = y;
    p_traj(:, 1) = p;
    q_traj(:, 1) = q;

    % RK4 integration loop
    for k = 1:N
        % k1 values
        [Fx1, Fy1, Fp1, Fq1] = F(x_traj(:, k), y_traj(:, k), g, p_traj(:, k), q_traj(:, k));

        k1x = dt * Fx1;
        k1y = dt * Fy1;
        k1p = dt * Fp1;
        k1q = dt * Fq1;

        % k2 values
        [Fx2, Fy2, Fp2, Fq2] = F(x_traj(:, k) + 0.5 * k1x, y_traj(:, k) + 0.5 * k1y, g, ...
                                   p_traj(:, k) + 0.5 * k1p, q_traj(:, k) + 0.5 * k1q);

        k2x = dt * Fx2;
        k2y = dt * Fy2;
        k2p = dt * Fp2;
        k2q = dt * Fq2;

        % k3 values
        [Fx3, Fy3, Fp3, Fq3] = F(x_traj(:, k) + 0.5 * k2x, y_traj(:, k) + 0.5 * k2y, g, ...
                                   p_traj(:, k) + 0.5 * k2p, q_traj(:, k) + 0.5 * k2q);

        k3x = dt * Fx3;
        k3y = dt * Fy3;
        k3p = dt * Fp3;
        k3q = dt * Fq3;

        % k4 values
        [Fx4, Fy4, Fp4, Fq4] = F(x_traj(:, k) + k3x, y_traj(:, k) + k3y, g, ...
                                   p_traj(:, k) + k3p, q_traj(:, k) + k3q);

        k4x = dt * Fx4;
        k4y = dt * Fy4;
        k4p = dt * Fp4;
        k4q = dt * Fq4;

        % Update vortex & particle
        x_traj(:, k+1) = x_traj(:, k) + (k1x + 2*k2x + 2*k3x + k4x) / 6;
        y_traj(:, k+1) = y_traj(:, k) + (k1y + 2*k2y + 2*k3y + k4y) / 6;
        p_traj(:, k+1) = p_traj(:, k) + (k1p + 2*k2p + 2*k3p + k4p) / 6;
        q_traj(:, k+1) = q_traj(:, k) + (k1q + 2*k2q + 2*k3q + k4q) / 6;
    end

    % Plot results: (MUST BE vortices blue, particles red)
    figure;
    hold on;
    
    % Plot vortex trajectories
    for i = 1:num_vortices
        plot(x_traj(i, :), y_traj(i, :), 'b-', 'DisplayName', sprintf('Vortex %d', i));
    end

    % Plot particle trajectories
    for i = 1:num_particles
        plot(p_traj(i, :), q_traj(i, :), 'r-', 'DisplayName', sprintf('Particle %d', i));
    end
    
    % Add funny words to make it look good
    xlabel('x');
    ylabel('y');
    title('Trajectories of Vortices and Particles');
    legend;
    grid on;
    hold off;
end


% Calling the function:

% Parameters
Gamma = 1;                      %VortexStrength
x0 = 1; y0 = 0;                 %InitConditions
r0 = sqrt(x0^2 + y0^2);         %RadiusOfMotion
T = 4 * pi^2;                   %PeriodOfRevolution
nsteps = [50, 100, 200, 400]; 

% Preallocate for error storage
errors = zeros(size(nsteps));

figure;
hold on;

for i = 1:length(nsteps)
    % Time step size
    N = nsteps(i);
    dt = T / N;
    
    % Initialize arrays for Euler's method
    x = zeros(1, N+1);
    y = zeros(1, N+1);
    t = linspace(0, T, N+1);
    
    % Initial conditions
    x(1) = x0;
    y(1) = y0;
    
    % Euler's method loop
    for n = 1:N
        r_squared = x(n)^2 + y(n)^2;
        dxdt = -y(n) / (2 * pi * r_squared);
        dydt = x(n) / (2 * pi * r_squared);
        
        % Update positions
        x(n+1) = x(n) + dt * dxdt;
        y(n+1) = y(n) + dt * dydt;
    end
    
    % Plot trajectory
    plot(x, y, 'DisplayName', sprintf('N = %d', N));
    
    % Compute error at t = T
    errors(i) = sqrt((x0 - x(end))^2 + (y0 - y(end))^2);
end

% Formatting the plot
legend;
xlabel('x');
ylabel('y');
title('Trajectories for different N using Euler''s method');
grid on;
hold off;

% Display errors
disp('Errors at t = T:');
for i = 1:length(nsteps)
    fprintf('N = %d, Error = %.6f\n', nsteps(i), errors(i));
end

function ESAM346_HW1_P1c(Nsteps, T, x0, y0, Gamma)
    % Inputs:
    % Nsteps - Vector of time step counts, e.g., [50, 100, 200, 400]
    % T - Final time
    % x0, y0 - Initial position of the particle
    % Gamma - Strength of the vortex
    Gamma = 1;                      %VortexStrength
    x0 = 1; y0 = 0;                 %InitConditions
    T = 4 * pi^2;                   %PeriodOfRevolution
    Nsteps = [50, 100, 200, 400];
    % Initialize exact solution and errors
    errors = zeros(size(Nsteps));

    % Exact trajectory at t = T
    r0 = sqrt(x0^2 + y0^2); % Initial radius
    omega = Gamma / (2 * pi * r0^2); % Angular velocity
    x_exact = r0 * cos(omega * T);
    y_exact = r0 * sin(omega * T);
    
    % Plot setup
    figure;
    hold on;
    colors = lines(length(Nsteps));
    
    for i = 1:length(Nsteps)
        % Current step size and number of steps
        N = Nsteps(i);
        dt = T / N;
        
        % Initialize positions
        x = x0;
        y = y0;
        trajectory_x = zeros(N+1, 1);
        trajectory_y = zeros(N+1, 1);
        trajectory_x(1) = x;
        trajectory_y(1) = y;
        
        % RK4 method loop
        for n = 1:N
            % Step 1: Get k1
            [Fx_k1, Fy_k1] = F(x, y, Gamma, x, y);
            k1x = dt * Fx_k1;
            k1y = dt * Fy_k1;
            
            % Step 2: Get k2
            x_k2 = x + 0.5 * k1x;
            y_k2 = y + 0.5 * k1y;
            [Fx_k2, Fy_k2] = F(x_k2, y_k2, Gamma, x, y);
            k2x = dt * Fx_k2;
            k2y = dt * Fy_k2;
            
            % Step 3: Get k3
            x_k3 = x + 0.5 * k2x;
            y_k3 = y + 0.5 * k2y;
            [Fx_k3, Fy_k3] = F(x_k3, y_k3, Gamma, x, y);
            k3x = dt * Fx_k3;
            k3y = dt * Fy_k3;
            
            % Step 4: Get k4
            x_k4 = x + k3x;
            y_k4 = y + k3y;
            [Fx_k4, Fy_k4] = F(x_k4, y_k4, Gamma, x, y);
            k4x = dt * Fx_k4;
            k4y = dt * Fy_k4;
            
            % Update positions
            x = x + (k1x + 2*k2x + 2*k3x + k4x) / 6;
            y = y + (k1y + 2*k2y + 2*k3y + k4y) / 6;
            
            % Debugging: Print derivatives and updates
            fprintf('Step %d/%d: Fx = %.6f, Fy = %.6f, x = %.6f, y = %.6f\n', n, N, Fx_k1, Fy_k1, x, y);
            
            % Store trajectory
            trajectory_x(n+1) = x;
            trajectory_y(n+1) = y;
        end
        
        % Compute error at T
        errors(i) = sqrt((x - x_exact)^2 + (y - y_exact)^2);
        
        % Plot trajectory
        plot(trajectory_x, trajectory_y, 'Color', colors(i, :), 'DisplayName', sprintf('N = %d', N));
    end
    
    % Finalize plot
    xlabel('x');
    ylabel('y');
    title('Trajectories for RK4 with Different Time Steps (Debugged)');
    legend show;
    grid on;
    hold off;
    
    % Display errors
    fprintf('Errors at T for each N:\n');
    for i = 1:length(Nsteps)
        fprintf('N = %d: Error = %.6f\n', Nsteps(i), errors(i));
    end
end