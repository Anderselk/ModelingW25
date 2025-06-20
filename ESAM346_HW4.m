function ESAM346_HW4()
    % fiber_Anders function
    function [t, x, y, a] = fiber_Anders(N, T, b, x0, y0, a0)
        dt = T / N;
        t = linspace(0, T, N+1);
        x = zeros(1, N+1);
        y = zeros(1, N+1);
        a = zeros(1, N+1);
        x(1) = x0;
        y(1) = y0;
        a(1) = a0;

        for n = 1:N
            r = sqrt(x(n)^2 + y(n)^2);
            b_r = b(r);

            k1x = cos(a(n));
            k1y = sin(a(n));
            k1a = -b_r * ((x(n) * -sin(a(n)) + y(n) * cos(a(n))) / r);

            k2x = cos(a(n) + 0.5 * dt * k1a);
            k2y = sin(a(n) + 0.5 * dt * k1a);
            k2a = -b_r * (((x(n) + 0.5 * dt * k1x) * -sin(a(n) + 0.5 * dt * k1a) + ...
                           (y(n) + 0.5 * dt * k1y) * cos(a(n) + 0.5 * dt * k1a)) / r);

            k3x = cos(a(n) + 0.5 * dt * k2a);
            k3y = sin(a(n) + 0.5 * dt * k2a);
            k3a = -b_r * (((x(n) + 0.5 * dt * k2x) * -sin(a(n) + 0.5 * dt * k2a) + ...
                           (y(n) + 0.5 * dt * k2y) * cos(a(n) + 0.5 * dt * k2a)) / r);

            k4x = cos(a(n) + dt * k3a);
            k4y = sin(a(n) + dt * k3a);
            k4a = -b_r * (((x(n) + dt * k3x) * -sin(a(n) + dt * k3a) + ...
                           (y(n) + dt * k3y) * cos(a(n) + dt * k3a)) / r);

            x(n+1) = x(n) + (dt / 6) * (k1x + 2 * k2x + 2 * k3x + k4x);
            y(n+1) = y(n) + (dt / 6) * (k1y + 2 * k2y + 2 * k3y + k4y);
            a(n+1) = a(n) + (dt / 6) * (k1a + 2 * k2a + 2 * k3a + k4a);
        end
    end

    % fibersde_Anders function
    function [t, x, y, a] = fibersde_Anders(N, T, b, x0, y0, a0, A)
        dt = T / N;
        t = linspace(0, T, N+1);
        x = zeros(1, N+1);
        y = zeros(1, N+1);
        a = zeros(1, N+1);
        x(1) = x0;
        y(1) = y0;
        a(1) = a0;

        for n = 1:N
            r = sqrt(x(n)^2 + y(n)^2);
            b_r = b(r);
            dx_drift = cos(a(n));
            dy_drift = sin(a(n));
            da_drift = -b_r * ((x(n) * -sin(a(n)) + y(n) * cos(a(n))) / r);

            dW = sqrt(dt) * randn;  
            dx_diff = A * dW;
            dy_diff = A * dW;
            da_diff = A * dW;

            x(n+1) = x(n) + dx_drift * dt + dx_diff;
            y(n+1) = y(n) + dy_drift * dt + dy_diff;
            a(n+1) = a(n) + da_drift * dt + da_diff;
        end
    end

    % Run problems 1–4 (fiber_Anders)
    clc; clear; close all;

    T = 2 * pi;    
    a0 = pi / 2;   
    N_values = [100, 200, 400]; 
    fprintf('Problem 1: Circle with b(r) = r\n');
    x0 = 1; y0 = 0; b1 = @(r) r;  
    errors1 = zeros(1, length(N_values));

    figure; hold on;
    for i = 1:length(N_values)
        N = N_values(i);
        [t, x, y, a] = fiber_Anders(N, T, b1, x0, y0, a0);
        final_error = norm([x(end) - 1, y(end) - 0]);
        errors1(i) = final_error;
        fprintf('N = %d, Final Error = %.10f\n', N, final_error);
        plot(x, y, 'DisplayName', sprintf('Numerical (N = %d)', N));
    end
    theta = linspace(0, 2*pi, 1000);
    plot(cos(theta), sin(theta), 'k--', 'DisplayName', 'Exact Circle');
    xlabel('x'); ylabel('y'); title('Trajectory of Fiber Extrusion (Problem 1)');
    legend; axis equal; grid on; hold off;

    convergence_rate1 = diff(log(errors1)) ./ diff(log(N_values));
    fprintf('\nProblem 1 Convergence Rate (should be close to 4):\n');
    for i = 1:length(convergence_rate1)
        fprintf('Rate between N=%d and N=%d: %.4f\n', N_values(i), N_values(i+1), convergence_rate1(i));
    end

    fprintf('\nProblem 2: b(r) = 1/r for two initial radii\n');
    b2 = @(r) 1 / r;  
    initial_radii = [0.5, 1.5];

    for j = 1:length(initial_radii)
        x0 = initial_radii(j); y0 = 0;
        figure; hold on;
        for i = 1:length(N_values)
            N = N_values(i);
            [t, x, y, a] = fiber_Anders(N, T, b2, x0, y0, a0);
            final_error = norm([x(end) - x0, y(end) - 0]);
            fprintf('Radius = %.1f, N = %d, Final Error = %.10f\n', x0, N, final_error);
            plot(x, y, 'DisplayName', sprintf('Numerical (N = %d)', N));
        end
        theta = linspace(0, 2*pi, 1000);
        plot(x0 * cos(theta), x0 * sin(theta), 'k--', 'DisplayName', 'Exact Circle');
        xlabel('x'); ylabel('y'); title(sprintf('Problem 2, r = %.1f', x0));
        legend; axis equal; grid on; hold off;
    end

    fprintf('\nProblem 3: b(r) = 1/r, T = pi/2, x0 = 1, y0 = 0, a0 = 3*pi/4\n');
    T = pi / 2; x0 = 1; y0 = 0; a0 = 3*pi/4;
    b3 = @(r) 1 / r;
    N_values3 = [400, 800, 1600, 3200]; errors3 = zeros(1, length(N_values3));
    figure; hold on;
    for i = 1:length(N_values3)
        N = N_values3(i);
        [t, x, y, a] = fiber_Anders(N, T, b3, x0, y0, a0);
        final_error = norm([x(end) - x0, y(end) - 0]);
        errors3(i) = final_error;
        fprintf('N = %d, Final Error = %.10f\n', N, final_error);
        plot(x, y, 'DisplayName', sprintf('Numerical (N = %d)', N));
    end
    plot(1, 0, 'ro', 'DisplayName', 'Initial Position');
    xlabel('x'); ylabel('y'); title('Trajectory of Fiber Extrusion (Problem 3)');
    legend; axis equal; grid on; hold off;

    convergence_rate3 = diff(log(errors3)) ./ diff(log(N_values3));
    fprintf('\nProblem 3 Convergence Rate (should be close to 4):\n');
    for i = 1:length(convergence_rate3)
        fprintf('Rate between N=%d and N=%d: %.4f\n', N_values3(i), N_values3(i+1), convergence_rate3(i));
    end

    fprintf('\nProblem 4: Additional Input Values and Bending Functions\n');
    additional_cases = {
        {1600, 106 * pi, @(r) r, 1, 0, pi / 3, 'Trajectory for (N=1600, T=106π, b=@(r)r, x0=1, y0=0, a0=π/3)'}
        {1600, 91 * pi, @(r) r, 0.5, 0, pi / 2, 'Trajectory for (N=1600, T=91π, b=@(r)r, x0=0.5, y0=0, a0=π/2)'}
        {1600, 28 * pi, @(r) (1/r), 1, 0, pi / 2 + 0.01, 'Trajectory for (N=1600, T=28π, b=@(r)(1/r), x0=1, y0=0, a0=π/2+0.01)'}
    };

    for j = 1:length(additional_cases)
        params = additional_cases{j};
        N = params{1}; T = params{2}; b = params{3};
        x0 = params{4}; y0 = params{5}; a0 = params{6};
        title_text = params{7};
        [t, x, y, a] = fiber_Anders(N, T, b, x0, y0, a0);
        figure; hold on;
        plot(x, y, 'b-', 'DisplayName', sprintf('Numerical (N = %d)', N));
        plot(x0, y0, 'ro', 'DisplayName', 'Initial Position');
        xlabel('x'); ylabel('y'); title(title_text);
        legend; axis equal; grid on; hold off;
    end

    % Run problems 5–8 (fibersde_Anders)
    clc; clear A a0 b N N_values T x0 y0;

    a0 = pi / 2; b = @(r) r; x0 = 1; y0 = 0; T = 2 * pi;
    A = 0; N_values = [200, 400, 800]; errors = zeros(1, length(N_values));
    figure; hold on;
    fprintf('Convergence Study for Deterministic Part (A = 0)\n');
    for i = 1:length(N_values)
        N = N_values(i);
        [t, x, y, a] = fibersde_Anders(N, T, b, x0, y0, a0, A);
        final_error = norm([x(end) - 1, y(end) - 0]);
        errors(i) = final_error;
        fprintf('N = %d, Final Error = %.10e\n', N, final_error)
        plot(x, y, 'DisplayName', sprintf('Numerical (N = %d)', N));
    end
    theta = linspace(0, 2*pi, 1000);
    plot(cos(theta), sin(theta), 'k--', 'DisplayName', 'Exact Circle');
    xlabel('x'); ylabel('y'); title('Trajectory (Deterministic, A = 0)');
    legend; axis equal; grid on; hold off;

    convergence_rate = diff(log(errors)) ./ diff(log(N_values));
    fprintf('\nConvergence Rates (should be close to 1 for Euler-Maruyama):\n');
    for i = 1:length(convergence_rate)
        fprintf('Rate between N=%d and N=%d: %.4f\n', N_values(i), N_values(i+1), convergence_rate(i));
    end
    p = mean(convergence_rate); 
    N_est = ceil(N_values(end) * (errors(end) / 1e-3)^(1 / p));
    fprintf('\nEstimated N for error < 10^-3: %d\n', N_est);

    fprintf('\nStochastic Simulations with A = 0.1, 0.5, 1 (N = 3051)\n');
    A_values = [0.1, 0.5, 1]; N = 3051; num_runs = 10;
    for j = 1:length(A_values)
        A = A_values(j); figure; hold on;
        fprintf('\nNoise Level A = %.1f\n', A);
        for k = 1:num_runs
            [t, x, y, a] = fibersde_Anders(N, T, b, x0, y0, a0, A);
            plot(x, y, 'DisplayName', sprintf('Run %d', k));
        end
        plot(cos(theta), sin(theta), 'k--', 'DisplayName', 'Exact Circle');
        xlabel('x'); ylabel('y'); title(sprintf('Stochastic Trajectories (A = %.1f, N = %d)', A, N));
        legend; axis equal; grid on; hold off;
    end

    fprintf('\nSingle Sample Path for Problem 4 Cases (A = 0.1)\n');
    problem4_cases = {
        {1600, 106 * pi, @(r) r, 1, 0, pi / 3},
        {1600, 91 * pi, @(r) r, 0.5, 0, pi / 2},
        {1600, 28 * pi, @(r) 1 / r, 1, 0, pi / 2 + 0.01}
    };

    for k = 1:length(problem4_cases)
        params = problem4_cases{k};
        N = params{1}; T = params{2}; b = params{3};
        x0 = params{4}; y0 = params{5}; a0 = params{6}; A = 0.1;
        [t, x, y, a] = fibersde_Anders(N, T, b, x0, y0, a0, A);
        figure; plot(x, y, 'b', 'DisplayName', 'Sample Path'); hold on;
        plot(x0, y0, 'ro', 'DisplayName', 'Initial Position');
        xlabel('x'); ylabel('y');
        title(sprintf('Single Sample Path for Problem 4 Case %d (A = 0.1)', k));
        legend; axis equal; grid on; hold off;
    end

    fprintf('\nComparing Terminal Point Distributions for b(r) = r and b(r) = 1/r\n');
    N = 3051; T = 2 * pi; x0 = 1; y0 = 0; a0 = pi / 2; A = 0.1; num_runs = 2000;

    figure; hold on; title('Terminal Point Distribution (b(r) = r, A = 0.1)');
    xlabel('x'); ylabel('y'); axis equal; grid on;
    for k = 1:num_runs
        [t, x, y, a] = fibersde_Anders(N, T, @(r) r, x0, y0, a0, A);
        plot(x(end), y(end), 'k.', 'MarkerSize', 2);  
    end
    plot(1, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'No-Noise Terminal Point');
    hold off;

    figure; hold on; title('Terminal Point Distribution (b(r) = 1/r, A = 0.1)');
    xlabel('x'); ylabel('y'); axis equal; grid on;
    for k = 1:num_runs
        [t, x, y, a] = fibersde_Anders(N, T, @(r) 1 / r, x0, y0, a0, A);
        plot(x(end), y(end), 'b.', 'MarkerSize', 2);  
    end
    plot(1, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'No-Noise Terminal Point');
    hold off;

    fprintf('Terminal point distribution comparison complete.\n');
end
