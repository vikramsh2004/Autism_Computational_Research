%% Plot GABA_B S^n / (S^n + Kd) over time
% This script reuses the GABA_B receptor and G-protein kinetics, then
% evaluates the Hill-style channel activation equation using S over time.
clear; clc; % Omitted 'close all' so old windows stay open

%% Simulation setup
N_num = 1;
tspan = 500; % ms
dt = 0.01; % time step for Euler method
loop = floor(tspan/dt) + 1; % include the final tspan sample
time = (0:loop-1).*dt; % time vector in ms

%% GABA_B constants (Destexhe)
Cmax_GABAb = 1; % max NT concentration
Cdur_GABAb = 1; % NT pulse duration
K1_GABAb = 0.09; % forward binding rate
K2_GABAb = 0.008; % backward unbinding rate
K3_GABAb = 0.18; % G-protein production rate
K4_GABAb = 0.034; % G-protein decay rate
Kd_GABAb = 100; % dissociation constant of K+ channel
n_GABAb = 4; % number of G-protein binding sites

%% Neurotransmitter pulse
nt_GABAb = zeros(N_num, loop);
GABAb_start = round(75/dt) + 1;
GABAb_stop = GABAb_start + round(Cdur_GABAb/dt);
GABAb_stop = min(GABAb_stop, loop);
nt_GABAb(:, GABAb_start:GABAb_stop) = Cmax_GABAb;

%% Initial receptor, G-protein, and equation values
r_GABAb = zeros(N_num, loop);
S_GABAb = zeros(N_num, loop);
S_hill_GABAb = zeros(N_num, loop);

%% Euler method
for i = 1:loop-1
    % GABA_B receptor activation
    dGABAb_rdt = (K1_GABAb.*nt_GABAb(:, i).*(1-r_GABAb(:, i))) ...
        - (K2_GABAb.*r_GABAb(:, i));
    r_GABAb(:, i+1) = r_GABAb(:, i) + dt.*dGABAb_rdt;

    % S is the G-protein activation variable driven by receptor activation.
    dGABAb_Sdt = (K3_GABAb.*r_GABAb(:, i)) - (K4_GABAb.*S_GABAb(:, i));
    S_GABAb(:, i+1) = S_GABAb(:, i) + dt.*dGABAb_Sdt;

    % Evaluate S^n / (S^n + Kd)
    S_hill_GABAb(:, i) = (S_GABAb(:, i).^n_GABAb) ...
        ./ (S_GABAb(:, i).^n_GABAb + Kd_GABAb);
end

% Calculate final sample for plotting
S_hill_GABAb(:, end) = (S_GABAb(:, end).^n_GABAb) ...
    ./ (S_GABAb(:, end).^n_GABAb + Kd_GABAb);

%% Plot S^n / (S^n + Kd) over time
figure('Name', 'GABA_B S Hill equation', 'Color', 'w');
plot(time, S_hill_GABAb, 'Color', [0.3 0.1 0.7], 'LineWidth', 1.5)
title('GABA_B S^n / (S^n + K_d) over time')
xlabel('Time (ms)')
ylabel('S^n / (S^n + K_d)')
grid on
