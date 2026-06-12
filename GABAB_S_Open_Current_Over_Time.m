%% GABA_B S, open fraction, and current over time
clear; clc; % Omitted 'close all' so old windows stay open

%% Simulation setup
dt = 0.01; % ms
tspan = 500; % ms
time = 0:dt:tspan;
nt = numel(time);

%% GABA_B constants (Destexhe)
Cmax_GABAb = 1; % max neurotransmitter concentration
Cdur_GABAb = 1; % neurotransmitter pulse duration (ms)
K1_GABAb = 0.09; % forward binding rate
K2_GABAb = 0.008; % backward unbinding rate
K3_GABAb = 0.18; % S/G-protein production rate
K4_GABAb = 0.034; % S/G-protein decay rate
Kd_GABAb = 100; % dissociation constant of K+ channel
n_GABAb = 4; % number of binding sites
E_GABAb = -95; % reversal potential (mV)
g_GABAb = 0.0001; % maximum conductance

%% Postsynaptic voltage used to calculate current
v_post_GABA = -50; % mV

%% Neurotransmitter pulse
nt_GABAb = zeros(1, nt);
GABAb_start = 75; % ms
GABAb_stop = GABAb_start + Cdur_GABAb;
nt_GABAb(time >= GABAb_start & time <= GABAb_stop) = Cmax_GABAb;

%% State variables
r_GABAb = zeros(1, nt);
S_GABAb = zeros(1, nt);
GABAb_open = zeros(1, nt);
I_GABAb = zeros(1, nt);

%% Euler method
for i = 1:nt-1
    % GABA_B receptor activation
    drdt = K1_GABAb * nt_GABAb(i) * (1 - r_GABAb(i)) ...
        - K2_GABAb * r_GABAb(i);
    r_GABAb(i+1) = r_GABAb(i) + dt * drdt;

    % S/G-protein activation for GABA_B
    dSdt = K3_GABAb * r_GABAb(i) - K4_GABAb * S_GABAb(i);
    S_GABAb(i+1) = S_GABAb(i) + dt * dSdt;

    % GABA_B channel open fraction and current
    GABAb_open(i) = S_GABAb(i)^n_GABAb ...
        / (S_GABAb(i)^n_GABAb + Kd_GABAb);
    I_GABAb(i) = g_GABAb * GABAb_open(i) * (E_GABAb - v_post_GABA);
end

% Calculate final sample for plotting
GABAb_open(end) = S_GABAb(end)^n_GABAb ...
    / (S_GABAb(end)^n_GABAb + Kd_GABAb);
I_GABAb(end) = g_GABAb * GABAb_open(end) * (E_GABAb - v_post_GABA);

%% Plot S over time
figure('Name','GABA_B S over Time','Color','w');
plot(time, S_GABAb, 'Color', [0.9 0.5 0.1], 'LineWidth', 1.5)
title('GABA_B S over time')
xlabel('Time (ms)')
ylabel('S')
grid on

%% Plot S^n / (S^n + Kd) over time
figure('Name','GABA_B Open Fraction','Color','w');
plot(time, GABAb_open, 'Color', [0.3 0.1 0.7], 'LineWidth', 1.5)
title('GABA_B S^n / (S^n + K_d) over time')
xlabel('Time (ms)')
ylabel('S^n / (S^n + K_d)')
grid on

%% Plot GABA_B current over time
figure('Name','GABA_B Current','Color','w');
plot(time, I_GABAb, 'Color', [0.3 0.1 0.7], 'LineWidth', 1.5)
yline(0, 'k:')
title('GABA_B current over time')
xlabel('Time (ms)')
ylabel('Current')
grid on
