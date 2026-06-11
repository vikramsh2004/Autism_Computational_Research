%% NMDA EPSC code
clear; clc

%% Simulation setup
N_num = 1;
tspan = 500; % ms
dt = 0.01; % time step for Euler method
loop = ceil(tspan/dt); % number of Euler iterations
time = (0:loop-1).*dt; % time vector in ms

%% Postsynaptic voltage used to calculate current
v_post_NMDA = -70; % mV

%% NMDA constants
Cmax_NMDA = 1;       % max neurotransmitter concentration
Cdur_NMDA = 1;       % neurotransmitter pulse duration, ms
alpha_NMDA = 0.072; % forward binding rate
beta_NMDA = 0.0066; % backward unbinding rate
E_NMDA = 0;         % reversal potential, mV
g_NMDA = 0.0002;    % maximum conductance
mg_NMDA = 1;        % external magnesium concentration, mM

%% Neurotransmitter pulse
nt_NMDA = zeros(N_num, loop);

NMDA_start = round(25/dt);
NMDA_stop = NMDA_start + round(Cdur_NMDA/dt);

nt_NMDA(:, NMDA_start:NMDA_stop) = Cmax_NMDA;

%% Initial receptor values
r_NMDA = zeros(N_num, loop);

%% Magnesium block and current arrays
B_NMDA = zeros(N_num, loop);
I_NMDA = zeros(N_num, loop);

%% Euler method
for i = 1:loop-1

    % NMDA receptor open fraction
    dNMDA_rdt(:,i) = ...
        alpha_NMDA .* nt_NMDA(:,i) .* (1 - r_NMDA(:,i)) ...
        - beta_NMDA .* r_NMDA(:,i);

    r_NMDA(:,i+1) = r_NMDA(:,i) + dt .* dNMDA_rdt(:,i);

    % Magnesium block
    B_NMDA(:,i) = ...
        1 ./ (1 + exp(0.062 .* (-v_post_NMDA)) .* (mg_NMDA ./ 3.57));

    % NMDA current
    I_NMDA(:,i) = ...
        g_NMDA .* B_NMDA(:,i) .* r_NMDA(:,i) .* (E_NMDA - v_post_NMDA);

end

% Fill final magnesium block and current values for plotting.
B_NMDA(:,end) = ...
    1 ./ (1 + exp(0.062 .* (-v_post_NMDA)) .* (mg_NMDA ./ 3.57));
I_NMDA(:,end) = ...
    g_NMDA .* B_NMDA(:,end) .* r_NMDA(:,end) .* (E_NMDA - v_post_NMDA);

%% Plot NMDA components in a new figure each run
figure
set(gcf, 'Name', 'NMDA EPSC', 'Color', 'w');

subplot(4,1,1)
plot(time, nt_NMDA, 'k', 'LineWidth', 1.5)
title('NMDA Neurotransmitter Pulse')
xlabel('Time (ms)')
ylabel('NT concentration')
grid on

subplot(4,1,2)
plot(time, r_NMDA, 'm', 'LineWidth', 1.5)
title('NMDA Receptor Open Fraction')
xlabel('Time (ms)')
ylabel('r_{NMDA}')
grid on

subplot(4,1,3)
plot(time, B_NMDA, 'Color', [0.3 0.1 0.7], 'LineWidth', 1.5)
title('NMDA Magnesium Block')
xlabel('Time (ms)')
ylabel('B_{NMDA}')
grid on

subplot(4,1,4)
plot(time, I_NMDA, 'b', 'LineWidth', 1.5)
yline(0,'k:')
title('NMDA EPSC')
xlabel('Time (ms)')
ylabel('Current')
grid on
