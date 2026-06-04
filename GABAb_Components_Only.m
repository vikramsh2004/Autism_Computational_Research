%% GABA_B receptor components only
clear; close all; clc

%% Simulation setup
N_num = 1;
tspan = 500; % ms
dt = 0.01; % time step for Euler method
loop = ceil(tspan/dt); % number of Euler iterations
time = (0:loop-1).*dt; % time vector in ms

%% GABA_B constants (Destexhe)
Cmax_GABAb = 1;    % max NT concentration
Cdur_GABAb = 1;    % NT pulse duration
K1_GABAb = 0.09;    % forward binding rate
K2_GABAb = 0.0012;    % backward unbinding rate
K3_GABAb = 0.18;    % G-protein production rate
K4_GABAb = 0.034;    % G-protein decay rate
KD_GABAb = 100;    % dissociation constant of K+ channel
n_GABAb = 4;    % number of G-protein binding sites
E_GABAb = -95;    % reversal potential
g_GABAb = 0.0001;    % maximum conductance

%% Postsynaptic voltage used to calculate current
v_post_GABA = -50; % mV

%% Neurotransmitter pulse
nt_GABAb = zeros(N_num,loop);

GABAb_start = round(75/dt);
GABAb_stop = GABAb_start + round(Cdur_GABAb/dt);
nt_GABAb(:,GABAb_start:GABAb_stop) = Cmax_GABAb;

%% Initial receptor and G-protein values
r_GABAb = zeros(N_num,loop);
G_GABAb = zeros(N_num,loop);
GABAb_open = zeros(N_num,loop);

%% Current array
I_GABAb = zeros(N_num,loop);

%% Euler method
for i = 1:loop-1

% GABA_B receptor activation
dGABAb_rdt(:,i) = (K1_GABAb.*nt_GABAb(:,i).*(1-r_GABAb(:,i))) - (K2_GABAb.*r_GABAb(:,i));
r_GABAb(:,i+1) = r_GABAb(:,i) + dt.*dGABAb_rdt(:,i);

% G-protein activation for GABA_B
dGABAb_Gdt(:,i) = (K3_GABAb.*r_GABAb(:,i)) - (K4_GABAb.*G_GABAb(:,i));
G_GABAb(:,i+1) = G_GABAb(:,i) + dt.*dGABAb_Gdt(:,i);

% GABA_B channel open fraction and current
GABAb_open(:,i) = (G_GABAb(:,i).^n_GABAb)./(G_GABAb(:,i).^n_GABAb + KD_GABAb);
I_GABAb(:,i) = g_GABAb.*GABAb_open(:,i).*(E_GABAb-v_post_GABA);
end

% Calculate final sample for plotting
GABAb_open(:,end) = (G_GABAb(:,end).^n_GABAb)./(G_GABAb(:,end).^n_GABAb + KD_GABAb);
I_GABAb(:,end) = g_GABAb.*GABAb_open(:,end).*(E_GABAb-v_post_GABA);

%% Plot GABA_B components
figure('Name','GABA_B Components','Color','w');

subplot(5,1,1)
plot(time, nt_GABAb, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5)
title('GABA_B neurotransmitter pulse')
xlabel('Time (ms)')
ylabel('[GABA]')
grid on

subplot(5,1,2)
plot(time, r_GABAb, 'Color', [0.1 0.5 0.2], 'LineWidth', 1.5)
title('GABA_B receptor activation')
xlabel('Time (ms)')
ylabel('r_{GABA_B}')
grid on

subplot(5,1,3)
plot(time, G_GABAb, 'Color', [0.9 0.5 0.1], 'LineWidth', 1.5)
title('G-protein activation')
xlabel('Time (ms)')
ylabel('G')
grid on

subplot(5,1,4)
plot(time, GABAb_open, 'Color', [0.3 0.1 0.7], 'LineWidth', 1.5)
title('GABA_B K^+ channel open fraction')
xlabel('Time (ms)')
ylabel('Open fraction')
grid on

subplot(5,1,5)
plot(time, I_GABAb, 'Color', [0.3 0.1 0.7], 'LineWidth', 1.5)
yline(0,'k:')
title('GABA_B IPSC')
xlabel('Time (ms)')
ylabel('Current')
grid on
