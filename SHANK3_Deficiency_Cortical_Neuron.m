% =========================================================================
% Biophysical Simulation of SHANK3 Deficiency in a Cortical Neuron
% Model: Single-compartment Hodgkin-Huxley with intrinsic and synaptic currents
% =========================================================================

clear; clc; close all;

%% 1. Simulation Parameters
dt = 0.01;              % Time step (ms)
T_total = 350;          % Total simulation time (ms)
time = 0:dt:T_total;    % Time vector
N = length(time);       % Number of time steps

%% 2. Biophysical Parameters (WT vs SHANK3)
% Common Reversal Potentials (mV)
E_Na = 50; E_K = -77; E_L = -54.4; E_h = -43; E_Ca = 120;
E_exc = 0; E_inh = -70;

% Membrane Capacitance (uF/cm^2)
% SHANK3 spines are smaller, effectively reducing total local capacitance
Cm_WT = 1.0; 
Cm_SH = 0.8; 

% Leak Conductance (mS/cm^2)
g_L_WT = 0.3; g_L_SH = 0.3;

% -- INTRINSIC CURRENTS --
g_Na_WT = 120;  g_Na_SH = 120; 
g_K_WT = 36;    g_K_SH = 25;    % SHANK3: Reduced K+ for homeostatic excitability
g_h_WT = 0.2;   g_h_SH = 0.04;  % SHANK3: Loss of HCN channel scaffolding
g_T_WT = 0.05;  g_T_SH = 0.25;  % SHANK3: Upregulated T-type Calcium current

% -- SYNAPTIC CURRENTS --
g_AMPA_WT = 0.2;    g_AMPA_SH = 0.08;   % SHANK3: Reduced AMPA density
g_NMDA_WT = 0.1;    g_NMDA_SH = 0.04;   % SHANK3: Reduced NMDA density
tau_NMDA_WT = 50;   tau_NMDA_SH = 80;   % SHANK3: GluN2B shift (slower decay)
tau_AMPA = 2; 

%% 3. Stimulation Protocol
% We apply a hyperpolarizing pulse to observe input resistance and Ih,
% followed by a burst of excitatory synaptic inputs.
I_inj = zeros(1, N);
I_inj(time > 30 & time < 100) = -2.0; % Hyperpolarizing pulse

spike_times = [150, 160, 170]; % Excitatory synaptic input times (ms)

%% 4. Initialization Arrays
V_WT = -65 * ones(1, N); V_SH = -65 * ones(1, N);

% Gating variables (m, h, n for HH; r for Ih; u for IT)
m_WT=0; h_WT=1; n_WT=0; r_WT=0; u_WT=0;
m_SH=0; h_SH=1; n_SH=0; r_SH=0; u_SH=0;

% Synaptic gating variables
s_AMPA_WT=0; s_NMDA_WT=0;
s_AMPA_SH=0; s_NMDA_SH=0;

% Arrays to store currents for plotting
I_NMDA_WT_arr = zeros(1,N); I_NMDA_SH_arr = zeros(1,N);
I_T_WT_arr = zeros(1,N);    I_T_SH_arr = zeros(1,N);
I_h_WT_arr = zeros(1,N);    I_h_SH_arr = zeros(1,N);

%% 5. Main Euler Integration Loop
for i = 1:N-1
    
    % --- Current Voltage ---
    v_W = V_WT(i) + 1e-6; % Add eps to avoid divide by zero at specific voltages
    v_S = V_SH(i) + 1e-6;
    
    % --- Synaptic Activation ---
    % Add discrete spikes to synaptic variables
    if any(abs(time(i) - spike_times) < dt)
        s_AMPA_WT = s_AMPA_WT + 0.5; s_NMDA_WT = s_NMDA_WT + 0.5;
        s_AMPA_SH = s_AMPA_SH + 0.5; s_NMDA_SH = s_NMDA_SH + 0.5;
    end
    
    % Synaptic Decay
    s_AMPA_WT = s_AMPA_WT - (s_AMPA_WT / tau_AMPA) * dt;
    s_NMDA_WT = s_NMDA_WT - (s_NMDA_WT / tau_NMDA_WT) * dt;
    s_AMPA_SH = s_AMPA_SH - (s_AMPA_SH / tau_AMPA) * dt;
    s_NMDA_SH = s_NMDA_SH - (s_NMDA_SH / tau_NMDA_SH) * dt;

    % Mg2+ Block (Jahr & Stevens 1990) - Extracellular Mg = 1.0 mM
    B_WT = 1 / (1 + (1.0 / 3.57) * exp(-0.062 * v_W));
    B_SH = 1 / (1 + (1.0 / 3.57) * exp(-0.062 * v_S));
    
    % Synaptic Currents
    I_AMPA_WT = g_AMPA_WT * s_AMPA_WT * (v_W - E_exc);
    I_NMDA_WT = g_NMDA_WT * s_NMDA_WT * B_WT * (v_W - E_exc);
    I_AMPA_SH = g_AMPA_SH * s_AMPA_SH * (v_S - E_exc);
    I_NMDA_SH = g_NMDA_SH * s_NMDA_SH * B_SH * (v_S - E_exc);
    
    % --- WT Kinematics (Na, K, h, T) ---
    am = 0.1*(v_W+40)/(1-exp(-(v_W+40)/10)); bm = 4.0*exp(-(v_W+65)/18);
    ah = 0.07*exp(-(v_W+65)/20);             bh = 1.0/(1+exp(-(v_W+35)/10));
    an = 0.01*(v_W+55)/(1-exp(-(v_W+55)/10)); bn = 0.125*exp(-(v_W+65)/80);
    
    r_inf = 1/(1+exp((v_W+80)/8)); tau_r = 100; % Ih kinetics
    u_inf = 1/(1+exp((v_W+80)/5)); tau_u = 50;  % IT slow inactivation
    s_inf_WT = 1/(1+exp(-(v_W+60)/6));          % IT fast activation
    
    m_WT = m_WT + (am*(1-m_WT) - bm*m_WT)*dt;
    h_WT = h_WT + (ah*(1-h_WT) - bh*h_WT)*dt;
    n_WT = n_WT + (an*(1-n_WT) - bn*n_WT)*dt;
    r_WT = r_WT + ((r_inf - r_WT)/tau_r)*dt;
    u_WT = u_WT + ((u_inf - u_WT)/tau_u)*dt;
    
    % WT Currents
    I_Na_W = g_Na_WT * m_WT^3 * h_WT * (v_W - E_Na);
    I_K_W  = g_K_WT * n_WT^4 * (v_W - E_K);
    I_L_W  = g_L_WT * (v_W - E_L);
    I_h_W  = g_h_WT * r_WT * (v_W - E_h);
    I_T_W  = g_T_WT * s_inf_WT^2 * u_WT * (v_W - E_Ca);
    
    % --- SHANK3 Kinematics (Na, K, h, T) ---
    am_S = 0.1*(v_S+40)/(1-exp(-(v_S+40)/10)); bm_S = 4.0*exp(-(v_S+65)/18);
    ah_S = 0.07*exp(-(v_S+65)/20);             bh_S = 1.0/(1+exp(-(v_S+35)/10));
    an_S = 0.01*(v_S+55)/(1-exp(-(v_S+55)/10)); bn_S = 0.125*exp(-(v_S+65)/80);
    
    r_inf_S = 1/(1+exp((v_S+80)/8)); 
    u_inf_S = 1/(1+exp((v_S+80)/5)); 
    s_inf_SH = 1/(1+exp(-(v_S+60)/6));
    
    m_SH = m_SH + (am_S*(1-m_SH) - bm_S*m_SH)*dt;
    h_SH = h_SH + (ah_S*(1-h_SH) - bh_S*h_SH)*dt;
    n_SH = n_SH + (an_S*(1-n_SH) - bn_S*n_SH)*dt;
    r_SH = r_SH + ((r_inf_S - r_SH)/tau_r)*dt;
    u_SH = u_SH + ((u_inf_S - u_SH)/tau_u)*dt;
    
    % SHANK3 Currents
    I_Na_S = g_Na_SH * m_SH^3 * h_SH * (v_S - E_Na);
    I_K_S  = g_K_SH * n_SH^4 * (v_S - E_K);
    I_L_S  = g_L_SH * (v_S - E_L);
    I_h_S  = g_h_SH * r_SH * (v_S - E_h);
    I_T_S  = g_T_SH * s_inf_SH^2 * u_SH * (v_S - E_Ca);
    
    % --- Voltage Update (dV/dt) ---
    V_WT(i+1) = v_W + (dt/Cm_WT) * (I_inj(i) - (I_Na_W + I_K_W + I_L_W + I_h_W + I_T_W + I_AMPA_WT + I_NMDA_WT));
    V_SH(i+1) = v_S + (dt/Cm_SH) * (I_inj(i) - (I_Na_S + I_K_S + I_L_S + I_h_S + I_T_S + I_AMPA_SH + I_NMDA_SH));
    
    % Store for plotting
    I_NMDA_WT_arr(i) = I_NMDA_WT; I_NMDA_SH_arr(i) = I_NMDA_SH;
    I_T_WT_arr(i) = I_T_W;        I_T_SH_arr(i) = I_T_S;
    I_h_WT_arr(i) = I_h_W;        I_h_SH_arr(i) = I_h_S;
end

%% 6. Plotting the Results
% Wild-Type and SHANK3 KO are shown on SEPARATE plots (no overlay):
%   Left  column  -> Wild-Type
%   Right column  -> SHANK3 KO
% Each row shares the same y-axis limits so the two genotypes stay
% directly comparable side-by-side.
figure('Name', 'SHANK3 Biophysical Simulation', 'Position', [100, 100, 1100, 750]);

% ---- Shared y-limits per row (computed so both panels use one scale) ----
% pad_lim() adds a 5% margin (min 1 unit) so paired panels use one y-scale.
pad_lim = @(d) [min(d(:)), max(d(:))] + [-1, 1] * max(0.05*(max(d(:))-min(d(:))), 1);
V_lim      = pad_lim([V_WT, V_SH]);
I_int_lim  = pad_lim([I_h_WT_arr, I_h_SH_arr, I_T_WT_arr, I_T_SH_arr]);
I_nmda_lim = pad_lim([I_NMDA_WT_arr, I_NMDA_SH_arr]);

% 6.1 Membrane Voltage -- Wild-Type
subplot(3,2,1);
plot(time, V_WT, 'k', 'LineWidth', 1.5);
title('Membrane Voltage: Wild-Type');
ylabel('Voltage (mV)'); ylim(V_lim); grid on;
text(50, V_lim(1)+3, '\leftarrow Hyperpolarizing Pulse', 'FontSize', 9);
text(150, V_lim(2)-8, '\leftarrow Synaptic Burst', 'FontSize', 9);

% 6.1b Membrane Voltage -- SHANK3 KO
subplot(3,2,2);
plot(time, V_SH, 'r', 'LineWidth', 1.5);
title('Membrane Voltage: SHANK3 KO (Hyperexcitable)');
ylabel('Voltage (mV)'); ylim(V_lim); grid on;
text(50, V_lim(1)+3, '\leftarrow Hyperpolarizing Pulse', 'FontSize', 9);
text(150, V_lim(2)-8, '\leftarrow Synaptic Burst', 'FontSize', 9);

% 6.2 Intrinsic Currents (Ih and IT) -- Wild-Type
subplot(3,2,3);
plot(time, I_h_WT_arr, 'k', 'LineWidth', 1.2); hold on;
plot(time, I_T_WT_arr, 'b', 'LineWidth', 1.2);
title('Intrinsic Currents (Wild-Type)');
ylabel('Current (\muA/cm^2)'); ylim(I_int_lim); grid on;
legend('Ih (WT)', 'IT (WT)');

% 6.2b Intrinsic Currents (Ih and IT) -- SHANK3 KO
subplot(3,2,4);
plot(time, I_h_SH_arr, 'k', 'LineWidth', 1.2); hold on;
plot(time, I_T_SH_arr, 'b', 'LineWidth', 1.2);
title('Intrinsic Currents (SHANK3 KO)');
ylabel('Current (\muA/cm^2)'); ylim(I_int_lim); grid on;
legend('Ih (SHANK3) - Loss of Leak', 'IT (SHANK3) - Rebound Burst');

% 6.3 NMDA Current -- Wild-Type
subplot(3,2,5);
plot(time, I_NMDA_WT_arr, 'g', 'LineWidth', 1.5);
title('NMDA Current (Wild-Type)');
xlabel('Time (ms)'); ylabel('Current (\muA/cm^2)');
ylim(I_nmda_lim); xlim([140 250]); grid on;
legend('NMDA (WT)');

% 6.3b NMDA Current -- SHANK3 KO
subplot(3,2,6);
plot(time, I_NMDA_SH_arr, 'g', 'LineWidth', 1.5);
title('NMDA Current & Mg2+ Block Failure (SHANK3 KO)');
xlabel('Time (ms)'); ylabel('Current (\muA/cm^2)');
ylim(I_nmda_lim); xlim([140 250]); grid on;
legend('NMDA (SHANK3) - Starved of Voltage');

disp('Simulation Complete. Observing plots...');
