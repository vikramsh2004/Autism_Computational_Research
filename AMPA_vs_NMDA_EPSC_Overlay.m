%% AMPA vs NMDA EPSC — Overlaid Comparison
% Runs both the AMPA and NMDA synapse models on a shared time base and
% overlays their comparable traces (neurotransmitter pulse, receptor open
% fraction, and EPSC) on the same axes so they can be compared visually.
% The NMDA current here is computed WITHOUT the magnesium block.
clear; clc

%% Distinct colors so AMPA and NMDA are easy to tell apart in every panel
col_AMPA = [0.00 0.45 0.74];   % blue
col_NMDA = [0.85 0.33 0.10];   % orange-red

%% Simulation setup (shared by both models)
N_num = 1;
tspan = 500; % ms
dt = 0.01; % time step for Euler method
loop = ceil(tspan/dt); % number of Euler iterations
time = (0:loop-1).*dt; % time vector in ms

%% Postsynaptic voltage used to calculate current
v_post = -70; % mV (shared)

%% AMPA constants
Cmax_AMPA = 1;      % max neurotransmitter concentration
Cdur_AMPA = 1;      % neurotransmitter pulse duration, ms
alpha_AMPA = 1.1;   % forward binding rate
beta_AMPA = 0.5;    % backward unbinding rate
E_AMPA = 0;         % reversal potential, mV
g_AMPA = 0.0001;    % maximum conductance

%% NMDA constants
Cmax_NMDA = 1;      % max neurotransmitter concentration
Cdur_NMDA = 1;      % neurotransmitter pulse duration, ms
alpha_NMDA = 0.072; % forward binding rate
beta_NMDA = 0.015;  % backward unbinding rate
E_NMDA = 0;         % reversal potential, mV
g_NMDA = 0.0002;    % maximum conductance

%% Neurotransmitter pulses (identical timing for both)
nt_AMPA = zeros(N_num, loop);
nt_NMDA = zeros(N_num, loop);

AMPA_start = round(25/dt);
AMPA_stop = AMPA_start + round(Cdur_AMPA/dt);
nt_AMPA(:, AMPA_start:AMPA_stop) = Cmax_AMPA;

NMDA_start = round(25/dt);
NMDA_stop = NMDA_start + round(Cdur_NMDA/dt);
nt_NMDA(:, NMDA_start:NMDA_stop) = Cmax_NMDA;

%% Receptor / current arrays
r_AMPA = zeros(N_num, loop);
I_AMPA = zeros(N_num, loop);

r_NMDA = zeros(N_num, loop);
I_NMDA = zeros(N_num, loop);

%% Euler method (both models integrated together)
for i = 1:loop-1

    % AMPA receptor open fraction
    dAMPA_rdt(:,i) = ...
        alpha_AMPA .* nt_AMPA(:,i) .* (1 - r_AMPA(:,i)) ...
        - beta_AMPA .* r_AMPA(:,i);
    r_AMPA(:,i+1) = r_AMPA(:,i) + dt .* dAMPA_rdt(:,i);

    % AMPA current
    I_AMPA(:,i) = g_AMPA .* r_AMPA(:,i) .* (E_AMPA - v_post);

    % NMDA receptor open fraction
    dNMDA_rdt(:,i) = ...
        alpha_NMDA .* nt_NMDA(:,i) .* (1 - r_NMDA(:,i)) ...
        - beta_NMDA .* r_NMDA(:,i);
    r_NMDA(:,i+1) = r_NMDA(:,i) + dt .* dNMDA_rdt(:,i);

    % NMDA current (magnesium block removed: B = 1)
    I_NMDA(:,i) = ...
        g_NMDA .* r_NMDA(:,i) .* (E_NMDA - v_post);

end

% Fill final current value for plotting (no magnesium block).
I_NMDA(:,end) = ...
    g_NMDA .* r_NMDA(:,end) .* (E_NMDA - v_post);

%% Plot: AMPA and NMDA overlaid on shared axes
figure
set(gcf, 'Name', 'AMPA vs NMDA EPSC (Overlay)', 'Color', 'w');

% Panel 1: Neurotransmitter pulse
subplot(3,1,1)
plot(time, nt_AMPA, 'Color', col_AMPA, 'LineWidth', 1.5, 'DisplayName', 'AMPA NT')
hold on
plot(time, nt_NMDA, '--', 'Color', col_NMDA, 'LineWidth', 1.5, 'DisplayName', 'NMDA NT')
title('Neurotransmitter Pulse')
xlabel('Time (ms)')
ylabel('NT concentration')
legend('show')
grid on

% Panel 2: Receptor open fraction
subplot(3,1,2)
plot(time, r_AMPA, 'Color', col_AMPA, 'LineWidth', 1.5, 'DisplayName', 'r_{AMPA}')
hold on
plot(time, r_NMDA, 'Color', col_NMDA, 'LineWidth', 1.5, 'DisplayName', 'r_{NMDA}')
title('Receptor Open Fraction')
xlabel('Time (ms)')
ylabel('r')
legend('show')
grid on

% Panel 3: EPSC current (NMDA without magnesium block)
subplot(3,1,3)
plot(time, I_AMPA, 'Color', col_AMPA, 'LineWidth', 1.5, 'DisplayName', 'I_{AMPA}')
hold on
plot(time, I_NMDA, 'Color', col_NMDA, 'LineWidth', 1.5, 'DisplayName', 'I_{NMDA} (no Mg block)')
yline(0,'k:', 'HandleVisibility', 'off')
title('EPSC (AMPA vs NMDA, no Mg block)')
xlabel('Time (ms)')
ylabel('Current')
legend('show')
grid on
