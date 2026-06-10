%% GABAa receptor dynamics only
% Standalone version of the GABAa synapse from the EPSC/IPSC script.
% Edit the parameters below to explore this one receptor in isolation.
clear; close all; clc

%% Simulation setup
N_num = 1;
tspan = 500; % ms
dt = 0.01; % time step for Euler method
loop = ceil(tspan/dt) + 1; % number of time samples
time = (0:loop-1).*dt; % time vector in ms

%% Postsynaptic voltage used to calculate current
v_post_GABA = -50; % mV

%% GABAa constants (Destexhe)
Cmax_GABAa = 1; % max neurotransmitter concentration
Cdur_GABAa = 1; % neurotransmitter pulse duration, ms
alpha_GABAa = 5; % forward binding rate
beta_GABAa = 0.18; % backward unbinding rate
E_GABAa = -80; % reversal potential, mV
g_GABAa = 0.0001; % maximum conductance

%% Neurotransmitter pulse
GABA_start_ms = 75; % pulse onset, ms
nt_GABAa = zeros(N_num, loop);

GABA_start = round(GABA_start_ms/dt) + 1;
GABA_stop = min(loop, GABA_start + round(Cdur_GABAa/dt) - 1);
nt_GABAa(:, GABA_start:GABA_stop) = Cmax_GABAa;

%% Initial receptor value
r_GABAa = zeros(N_num, loop);

%% Current array
I_GABAa = zeros(N_num, loop);
dGABAa_rdt = zeros(N_num, loop);

%% Euler method
for i = 1:loop-1
    % GABAa receptor open fraction
    dGABAa_rdt(:, i) = (alpha_GABAa.*nt_GABAa(:, i).*(1-r_GABAa(:, i))) ...
        - (beta_GABAa.*r_GABAa(:, i));
    r_GABAa(:, i+1) = r_GABAa(:, i) + dt.*dGABAa_rdt(:, i);

    % GABAa inhibitory postsynaptic current
    I_GABAa(:, i) = g_GABAa.*r_GABAa(:, i).*(E_GABAa-v_post_GABA);
end

% Fill the final current sample after the last receptor update.
I_GABAa(:, end) = g_GABAa.*r_GABAa(:, end).*(E_GABAa-v_post_GABA);

%% Summary values
[peak_current, peak_idx] = min(I_GABAa, [], 2);
peak_time = time(peak_idx);
charge_transfer = trapz(time, I_GABAa, 2);

fprintf('Peak GABAa current: %.6g at %.2f ms\n', peak_current, peak_time);
fprintf('Charge transfer: %.6g current*ms\n', charge_transfer);

%% Plot GABAa dynamics
figure('Name','GABAa Dynamics Only','Color','w');

subplot(4,1,1)
plot(time, nt_GABAa, 'Color', [0.1 0.5 0.1], 'LineWidth', 1.5)
title('GABAa neurotransmitter pulse')
xlabel('Time (ms)')
ylabel('[NT]')
grid on

subplot(4,1,2)
plot(time, r_GABAa, 'b', 'LineWidth', 1.5)
title('GABAa receptor open fraction')
xlabel('Time (ms)')
ylabel('r_{GABAa}')
ylim([0 1])
grid on

subplot(4,1,3)
plot(time, I_GABAa, 'b', 'LineWidth', 1.5)
yline(0,'k:')
title('GABAa IPSC')
xlabel('Time (ms)')
ylabel('Current')
grid on

subplot(4,1,4)
plot(time, cumtrapz(time, I_GABAa), 'k', 'LineWidth', 1.5)
yline(0,'k:')
title('Cumulative GABAa charge transfer')
xlabel('Time (ms)')
ylabel('Current*ms')
grid on
