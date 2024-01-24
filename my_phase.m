
% Parameters
l = 1.0;    % length of the pendulum
g = 10;    % acceleration due to gravity
dt = 0.01;  % time step for integration
num_steps = 1000;  % number of steps for integration

% Generate and display the phase portrait
generate_phase_portrait(l, g, dt, num_steps);

% Function to calculate the derivatives of the pendulum equations
function [theta_dot, omega_dot] = pendulum_derivatives(theta, omega, l, g)
    theta_dot = omega;
    omega_dot = -g / l * sin(theta);
end

% Function to perform Euler integration
function [theta, omega] = euler_integration(theta, omega, l, g, dt)
    [theta_dot, omega_dot] = pendulum_derivatives(theta, omega, l, g);
    theta = theta + theta_dot * dt;
    omega = omega + omega_dot * dt;
end

% Function to generate the phase portrait
function generate_phase_portrait(l, g, dt, num_steps)
    % Initial conditions
    theta_initial = linspace(-2 * pi, 2 * pi, 200);
    omega_initial = zeros(size(theta_initial));

    % Arrays to store results
    theta_values = zeros(num_steps, length(theta_initial));
    omega_values = zeros(num_steps, length(theta_initial));

    % Perform integration for each set of initial conditions
    for i = 1:length(theta_initial)
        theta = theta_initial(i);
        omega = omega_initial(i);

        for step = 1:num_steps
            [theta, omega] = euler_integration(theta, omega, l, g, dt);
            theta_values(step, i) = theta;
            omega_values(step, i) = omega;
        end
    end

    % Plot the phase portrait
    figure;
    hold on;

    % Plot arrows using quiver
    quiver(theta_values(1:end-1, :), omega_values(1:end-1, :), ...
        diff(theta_values), diff(omega_values), 0.5, 'Color', 'b');

    % Mark fixed points
    scatter([0, pi], [0, 0], 100, 'red', 'filled', 'o', 'LineWidth', 1.5, 'MarkerEdgeColor', 'black');

    % Mark attractors (for illustrative purposes, you can modify based on your analysis)
    %scatter([0, pi, 3, -3], [0, 0, -2, 2], 100, 'black', '*', 'LineWidth', 1.5);

    title('Pendulum Phase Portrait','FontSize',25);
    xlabel('Theta (radians)','FontSize',20);
    ylabel('Omega (angular velocity)','FontSize',20);
    legend('Arrows','Fixed Points','FontSize', 20);
    grid on;
    hold off;
end
