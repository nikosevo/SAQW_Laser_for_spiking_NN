% LIF Neuron Parameters
tau_m = 10;          % Membrane time constant (ms)
V_rest = -70;        % Resting membrane potential (mV)
R = 10;              % Membrane resistance (Mohm)
V_th = -55;          % Threshold potential (mV)
V_reset = V_rest;    % Reset potential (mV)

% Time parameters
dt = 0.1;            % Time step (ms)
t_max = 100;         % Maximum simulation time (ms)

% Initial conditions
V0 = V_rest;         % Initial membrane potential

% Arrays to store results
time = 0:dt:t_max;
V_values = zeros(size(time));

% Simulation loop
V = V0;
for i = 1:length(time)
    % LIF neuron dynamics
    dV_dt = -(V - V_rest) / tau_m;
    V = V + dV_dt * dt;
    
    % Spike generation (reset if threshold is crossed)
    if V >= V_th
        V = V_reset;
    end
    
    % Store membrane potential at each time step
    V_values(i) = V;
end

% Plotting the membrane potential trajectory with spikes (stair plot)
figure;
stairs(time, V_values, 'b', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('LIF Neuron Membrane Potential Trajectory with Spikes');
grid on;
  