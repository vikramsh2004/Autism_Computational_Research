%% PFC Neuron: VTA Dopamine Pooling + AMPA/NMDA Upregulation -> Spiking
% ========================================================================
% WHAT THIS SCRIPT SIMULATES
% ------------------------------------------------------------------------
% A single prefrontal-cortex (PFC) pyramidal neuron receives:
%   (1) a steady background glutamatergic drive (AMPA + NMDA synapses), and
%   (2) a dopaminergic (DA) projection from the ventral tegmental area (VTA).
%
% The VTA fires TONICALLY at baseline and switches to a PHASIC BURST during
% a salient event (here 500-1000 ms, e.g. a social interaction). Each DA
% spike releases dopamine that POOLS in the extracellular space and then
% decays (time constant tau_DA). The pooled dopamine acts through a D1-like
% term (D1mod = 1 + D1) that UP-REGULATES the amplitude of both the AMPA and
% NMDA synaptic currents -- i.e. dopamine sets the synaptic GAIN of the PFC
% neuron.
%
% HOW TO READ THE OUTPUT  (Figure 1, top -> bottom, shared time axis)
%   1. VTA drive          : tonic vs phasic firing rate + DA spike raster.
%                           The shaded band marks the phasic burst window.
%   2. Extracellular [DA]  : dopamine pools during the burst, then decays.
%   3. Dopamine gain D1mod : the actual modulation factor applied to the
%                           synapses -- THE mechanism (peak % annotated).
%   4. DA-added synaptic   : the EXTRA AMPA & NMDA drive contributed purely
%      drive               : by dopamine (Delta g = g*r*[DA]). Zero at rest,
%                           grows during the burst -> "up-regulation".
%   5. Membrane voltage    : the PFC neuron's response.
%
%   Figure 2 quantifies the effect per epoch (Baseline / DA burst / Recovery):
%   the % synaptic up-regulation and the mean DA-added AMPA/NMDA drive.
%
% Tunable knobs to explore: r_phasic_VTA, tau_DA, B_D1 (DA release gain),
% g_AMPA, g_NMDA, and the burst window (burst_on / burst_off).
%
% NOTE ON EFFECT SIZE: with the default B_D1 = 0.05 the dopamine pool reaches
% only D1 ~ 0.006, i.e. a ~0.6% synaptic gain increase (the plots autoscale
% so this small, real effect is still clearly visible in its timing/shape).
% Increase B_D1 (e.g. 2-5) to move into the strong-modulation regime where
% the up-regulation becomes large; the figures rescale automatically.
% ========================================================================

clear; close all; clc;
rng(42);                   % reproducible Poisson input for stable figures

% 1. SIMULATION PARAMETERS
dt = 0.01;                 % Time step (ms)
t_max = 1500;              % Total simulation time (ms)
t = 0:dt:t_max;            
N = length(t);             

% Phasic burst window (VTA salient-event burst)
burst_on  = 500;           % ms
burst_off = 1000;          % ms

% 2. INPUT PARAMETERS (Glutamate & Dopamine)
r_pre_Glu = 30;            % Steady background Glutamate rate (Hz)
release_steps = round(1/dt); 
counter_Glu = 0;           

r_tonic_VTA = 4;           % Baseline VTA dopamine firing (Hz)
r_phasic_VTA = 40;         % Burst VTA firing (Hz) - The trigger for things like social interaction
tau_DA = 300;              % Dopamine decay constant (ms)
A_D1 = 1 / tau_DA;         % Dopamine decay rate
B_D1 = 0.05;               % Scaling coefficient for DA release

% 3. MEMBRANE & CHANNEL PARAMETERS
C = 1.0;                   
g_L = 0.1;   E_L = -65;    
g_Na = 120;  E_Na = 50;    
g_K = 36;    E_K = -77;    
E_syn = 0;   Mg = 1.0;                  

% Base conductances (Set intentionally low so it only fires when modulated)
g_AMPA = 0.10;  alpha_AMPA = 1.0;   beta_AMPA = 0.5;         
g_NMDA = 0.05;  alpha_NMDA = 0.1;   beta_NMDA = 0.015;      

% 4. INITIALIZE VARIABLES
V = zeros(1, N);
V(1) = E_L;                

m = zeros(1, N); h = zeros(1, N); n = zeros(1, N);
v_m = V(1) + 40; a_m = 0.1 * v_m / (1 - exp(-v_m/10)); b_m = 4 * exp(-(V(1)+65)/18);
a_h = 0.07 * exp(-(V(1)+65)/20); b_h = 1 / (1 + exp(-(V(1)+35)/10));
v_n = V(1) + 55; a_n = 0.01 * v_n / (1 - exp(-v_n/10)); b_n = 0.125 * exp(-(V(1)+65)/80);
m(1) = a_m / (a_m + b_m); h(1) = a_h / (a_h + b_h); n(1) = a_n / (a_n + b_n);

r_AMPA = zeros(1, N);
r_NMDA = zeros(1, N);
D1 = zeros(1, N);          
I_AMPA = zeros(1, N);
I_NMDA = zeros(1, N);

% Recording arrays for visualising the dopamine mechanism
VTA_rate = zeros(1, N);    % instantaneous VTA firing rate (Hz)
DA_event = false(1, N);    % logical: a DA spike occurred at this step

% 5. MAIN INTEGRATION LOOP
for i = 1:N-1
    
    % --- A. INPUT GENERATORS ---
    % 1. Constant background Glutamate
    if poissrnd(r_pre_Glu * (dt/1000)) > 0
        counter_Glu = release_steps;
    end
    T_Glu = double(counter_Glu > 0);
    
    % 2. VTA Dopamine Burst (phasic during the salient-event window)
    if t(i) >= burst_on && t(i) <= burst_off
        VTA_rate(i) = r_phasic_VTA;
    else
        VTA_rate(i) = r_tonic_VTA;
    end
    lambda_VTA = VTA_rate(i) * (dt/1000);
    T_D = double(poissrnd(lambda_VTA) > 0);
    DA_event(i) = T_D > 0;
    
    % --- B. SYNAPTIC & DOPAMINE KINETICS ---
    r_AMPA(i+1) = r_AMPA(i) + dt * ((alpha_AMPA * T_Glu * (1 - r_AMPA(i))) - (beta_AMPA * r_AMPA(i)));
    r_NMDA(i+1) = r_NMDA(i) + dt * ((alpha_NMDA * T_Glu * (1 - r_NMDA(i))) - (beta_NMDA * r_NMDA(i)));
    
    % Extracellular Dopamine Pooling & Modulation Upregulation
    dD1dt = -A_D1 * D1(i) + B_D1 * T_D;
    D1(i+1) = D1(i) + dt * dD1dt;
    D1mod = 1 + D1(i);
    
    % --- C. MODULATED SYNAPTIC CURRENTS ---
    B_V = 1 / (1 + (Mg / 3.57) * exp(-0.062 * V(i)));
    
    % D1mod directly scales up the amplitude of both AMPA and NMDA
    I_AMPA(i) = g_AMPA * r_AMPA(i) * (V(i) - E_syn) * D1mod;
    I_NMDA(i) = g_NMDA * r_NMDA(i) * B_V * (V(i) - E_syn) * D1mod; 
    I_Total_syn = I_AMPA(i) + I_NMDA(i);
    
    % --- D. HODGKIN-HUXLEY ACTIVE CHANNELS ---
    v_m = V(i) + 40;
    if abs(v_m) < 1e-6, a_m = 1; else, a_m = 0.1 * v_m / (1 - exp(-v_m/10)); end
    b_m = 4 * exp(-(V(i)+65)/18);
    a_h = 0.07 * exp(-(V(i)+65)/20); b_h = 1 / (1 + exp(-(V(i)+35)/10));
    v_n = V(i) + 55;
    if abs(v_n) < 1e-6, a_n = 0.1; else, a_n = 0.01 * v_n / (1 - exp(-v_n/10)); end
    b_n = 0.125 * exp(-(V(i)+65)/80);

    m(i+1) = m(i) + dt * (a_m * (1 - m(i)) - b_m * m(i));
    h(i+1) = h(i) + dt * (a_h * (1 - h(i)) - b_h * h(i));
    n(i+1) = n(i) + dt * (a_n * (1 - n(i)) - b_n * n(i));

    I_Na = g_Na * (m(i)^3) * h(i) * (V(i) - E_Na);
    I_K  = g_K  * (n(i)^4) * (V(i) - E_K);
    I_leak  = g_L  * (V(i) - E_L);

    % --- E. MEMBRANE VOLTAGE UPDATE ---
    dV = (-I_leak - I_Total_syn - I_Na - I_K) / C;
    V(i+1) = V(i) + dt * dV;
    
    counter_Glu = max(0, counter_Glu - 1);
end
I_AMPA(end) = I_AMPA(end-1); I_NMDA(end) = I_NMDA(end-1);
VTA_rate(end) = VTA_rate(end-1);

% ------------------------------------------------------------------------
% 6. DERIVED QUANTITIES FOR VISUALISING THE MODULATION
% ------------------------------------------------------------------------
% The dopamine gain that multiplies both synapses.
D1mod  = 1 + D1;                                  % modulation factor (>= 1)

% Effective synaptic drive g_syn = g * gating, WITHOUT and WITH dopamine.
% (Plotting the synaptic DRIVE isolates the dopamine modulation from the
%  fast, spike-driven swings in the (V - E_syn) term of the raw current.)
gA_base = g_AMPA .* r_AMPA;                        % AMPA drive, no DA
gA_mod  = gA_base .* D1mod;                        % AMPA drive, DA-modulated
gN_base = g_NMDA .* r_NMDA;                        % NMDA drive, no DA
gN_mod  = gN_base .* D1mod;                        % NMDA drive, DA-modulated

% The EXTRA synaptic drive contributed purely by dopamine (= drive * [DA]).
dgA = gA_mod - gA_base;                            % DA-added AMPA drive
dgN = gN_mod - gN_base;                            % DA-added NMDA drive

% Detect PFC spikes (upward crossing of 0 mV).
spk_idx   = find(V(1:end-1) <= 0 & V(2:end) > 0);
spk_times = t(spk_idx);

% Per-epoch summaries: Baseline | DA burst | Recovery.
epoch_edges = [0 burst_on burst_off t_max];
epoch_names = {'Baseline','DA burst','Recovery'};
mGain = zeros(1,3); mdgA = zeros(1,3); mdgN = zeros(1,3); nSpk = zeros(1,3);
for e = 1:3
    lo = epoch_edges(e); hi = epoch_edges(e+1);
    sel = t >= lo & t < hi;
    mGain(e) = mean(D1mod(sel));
    mdgA(e)  = mean(dgA(sel));
    mdgN(e)  = mean(dgN(sel));
    nSpk(e)  = sum(spk_times >= lo & spk_times < hi);
end

% Console summary.
peak_pct = (max(D1mod) - 1) * 100;
fprintf('\n%-10s | %8s | %11s | %11s | %6s\n', ...
        'Epoch','Gain(x)','+AMPA drive','+NMDA drive','Spikes');
fprintf('%s\n', repmat('-',1,58));
for e = 1:3
    fprintf('%-10s | %8.4f | %11.3e | %11.3e | %6d\n', ...
            epoch_names{e}, mGain(e), mdgA(e), mdgN(e), nSpk(e));
end
fprintf('%s\n', repmat('-',1,58));
fprintf('Peak dopamine gain: +%.2f%% of AMPA/NMDA amplitude (B_D1 = %.3g).\n', peak_pct, B_D1);
fprintf('Dopamine up-regulates AMPA/NMDA drive specifically during the burst.\n\n');

% ------------------------------------------------------------------------
% 7. FIGURE 1 : THE DOPAMINE MODULATION CASCADE
%    VTA burst -> [DA] pooling -> synaptic gain -> AMPA/NMDA up-regulation
% ------------------------------------------------------------------------
cDA   = [0.10 0.65 0.75];   % dopamine (teal)
cAMPA = [0.85 0.30 0.20];   % AMPA (red)
cNMDA = [0.55 0.25 0.70];   % NMDA (purple)
cV    = [0.15 0.30 0.75];   % voltage (blue)

figure('Name','Dopamine modulation cascade','Color','w','Position',[60 40 960 980]);
ax = zeros(1,5);

% --- Panel 1: VTA drive (tonic vs phasic) + DA spike raster ---
ax(1) = subplot(5,1,1);
plot(t, VTA_rate, 'Color',[0.2 0.2 0.2], 'LineWidth',1.4); hold on;
plot(t(DA_event), repmat(r_phasic_VTA*1.12, 1, nnz(DA_event)), '|', ...
     'Color',cDA, 'MarkerSize',7, 'LineWidth',1.0);
ylim([0 r_phasic_VTA*1.3]);
ylabel({'VTA rate','(Hz)'});
title('1) VTA dopamine drive: tonic baseline \rightarrow phasic burst  (ticks = DA spikes)');

% --- Panel 2: Extracellular dopamine pooling ---
ax(2) = subplot(5,1,2);
area(t, D1, 'FaceColor',cDA, 'EdgeColor',cDA, 'LineWidth',1.2); hold on;
try, set(get(gca,'Children'),'FaceAlpha',0.35); catch, end
ylabel({'[DA]','(a.u.)'});
title('2) Extracellular dopamine pools during the burst, then decays (\tau_{DA})');

% --- Panel 3: Dopamine modulation gain applied to the synapses ---
ax(3) = subplot(5,1,3);
plot(t, D1mod, 'Color',cDA, 'LineWidth',1.8); hold on;
ylabel('D1mod (\times)');
title(sprintf('3) Dopamine gain D1mod = 1 + [DA]  (multiplies AMPA & NMDA;  peak +%.2f%%)', peak_pct));

% --- Panel 4: DA-added synaptic drive (the up-regulation itself) ---
ax(4) = subplot(5,1,4);
hA = plot(t, dgA, 'Color',cAMPA, 'LineWidth',1.3); hold on;
hN = plot(t, dgN, 'Color',cNMDA, 'LineWidth',1.3);
ylabel({'DA-added','drive (a.u.)'});
legend([hA hN], {'AMPA','NMDA'}, 'Location','northeast');
title('4) Extra AMPA & NMDA drive contributed by dopamine  ( g\cdotr\cdot[DA] )');

% --- Panel 5: Membrane voltage ---
ax(5) = subplot(5,1,5);
plot(t, V, 'Color',cV, 'LineWidth',1.0); hold on;
if ~isempty(spk_times)
    plot(spk_times, repmat(28,1,numel(spk_times)), 'v', ...
         'MarkerFaceColor',cV, 'MarkerEdgeColor',cV, 'MarkerSize',4);
end
ylim([-90 40]);
ylabel({'V_m','(mV)'});
xlabel('Time (ms)');
title(sprintf('5) PFC membrane voltage  (%d spike(s))', numel(spk_times)));

% Shade the phasic burst window on every panel and link the time axes.
for k = 1:5
    yl = ylim(ax(k));
    p = patch('XData',[burst_on burst_off burst_off burst_on], ...
              'YData',[yl(1) yl(1) yl(2) yl(2)], ...
              'FaceColor',[1.0 0.85 0.45], 'EdgeColor','none', ...
              'HandleVisibility','off', 'Parent',ax(k));
    try, set(p,'FaceAlpha',0.18); catch, end
    try, uistack(p,'bottom'); catch, end
    ylim(ax(k), yl);
    grid(ax(k),'on');
    xlim(ax(k),[0 t_max]);
end
linkaxes(ax,'x');

% ------------------------------------------------------------------------
% 8. FIGURE 2 : QUANTIFYING THE MODULATION PER EPOCH
% ------------------------------------------------------------------------
figure('Name','Dopamine modulation summary','Color','w','Position',[120 120 940 380]);

subplot(1,2,1);
bar((mGain - 1) * 100, 0.6, 'FaceColor', cDA);
set(gca,'XTickLabel',epoch_names);
ylabel('Synaptic up-regulation (%)'); grid on;
title('Mean dopamine gain by epoch');

subplot(1,2,2);
hb = bar([mdgA(:), mdgN(:)], 'grouped');
try, hb(1).FaceColor = cAMPA; hb(2).FaceColor = cNMDA; catch, end
set(gca,'XTickLabel',epoch_names);
ylabel('Mean DA-added drive (a.u.)'); grid on;
legend({'AMPA','NMDA'},'Location','northwest');
title('Mean DA-added AMPA/NMDA drive by epoch');
