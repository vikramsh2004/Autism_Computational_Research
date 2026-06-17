%% TRN <-> Thalamic neuron microcircuit
% One thalamocortical (TC) relay neuron and one thalamic reticular nucleus
% (TRN) neuron are coupled with reciprocal synapses:
%   TC  -> TRN: glutamatergic excitation with AMPA and NMDA components
%   TRN -> TC : GABAergic inhibition with GABA_A and GABA_B components
%
% Synaptic kinetics are adapted from the Destexhe-style receptor equations
% in the reference script. Intrinsic membrane dynamics use a compact
% Hodgkin-Huxley model so each neuron can spike while transmitter release can
% be run as a staggered visualization protocol or as spike-triggered release.

clear; close all; clc

%% Simulation setup
dt = 0.01; % ms
t_end = 300; % ms
time = 0:dt:t_end;
nt = numel(time);

%% Intrinsic membrane parameters
Cm = 1.0; % uF/cm^2

% TC relay neuron
gNa_TC = 120; ENa_TC = 50;
gK_TC = 36; EK_TC = -77;
gL_TC = 0.3; EL_TC = -54.4;

% TRN neuron
gNa_TRN = 100; ENa_TRN = 50;
gK_TRN = 30; EK_TRN = -77;
gL_TRN = 0.25; EL_TRN = -58;

spike_threshold = 0; % mV, upward crossing used for spike-triggered release

%% External drive
% Staggered stimuli make the two directions of the loop easier to see:
% TC drive evokes TC -> TRN excitation first, and later TRN drive evokes
% TRN -> TC inhibition after the excitatory response has separated in time.
I_ext_TC = zeros(1, nt);
I_ext_TRN = zeros(1, nt);

TC_drive_early = time >= 20 & time < 38;
TC_drive_late = time >= 170 & time < 188;
TRN_drive_early = time >= 90 & time < 108;
TRN_drive_late = time >= 235 & time < 253;

I_ext_TC(TC_drive_early) = 10;
I_ext_TC(TC_drive_late) = 8;
I_ext_TRN(TRN_drive_early) = 9;
I_ext_TRN(TRN_drive_late) = 8;

% By default, the demo uses staggered presynaptic release times so the
% glutamatergic and GABAergic components are visually separated. Set this
% false to make transmitter release occur only from threshold-crossing
% presynaptic spikes.
use_staggered_release_protocol = true;
TC_release_times = [25 175]; % ms, TC -> TRN AMPA/NMDA release
TRN_release_times = [95 235]; % ms, TRN -> TC GABA_A/GABA_B release

%% Synaptic constants
% The receptor kinetics, neurotransmitter pulse durations, reversal
% potentials, and Mg2+ value below use the Destexhe values from the
% reference script. Synaptic conductances are scaled for this two-cell
% conductance-based demo so the postsynaptic responses are visible.

% TC -> TRN excitation
Cmax_AMPA = 1;
Cdur_AMPA = 1; % ms
alpha_AMPA = 1.1;
beta_AMPA = 0.19;
E_AMPA = 0;
g_AMPA_TC_TRN = 0.18; % mS/cm^2

Cmax_NMDA = 1;
Cdur_NMDA = 1; % ms
alpha_NMDA = 0.072;
beta_NMDA = 0.0066;
E_NMDA = 0;
g_NMDA_TC_TRN = 0.12; % mS/cm^2
mg_NMDA = 1; % mM

% TRN -> TC inhibition
Cmax_GABAa = 1;
Cdur_GABAa = 1; % ms
alpha_GABAa = 5;
beta_GABAa = 0.18;
E_GABAa = -80;
g_GABAa_TRN_TC = 0.20; % mS/cm^2

Cmax_GABAb = 1;
Cdur_GABAb = 1; % ms
K1_GABAb = 0.09;
K2_GABAb = 0.0012;
K3_GABAb = 0.18;
K4_GABAb = 0.034;
KD_GABAb = 0.08; % normalized G-protein half-activation for visibility
n_GABAb = 4;
E_GABAb = -95;
g_GABAb_TRN_TC = 0.04; % mS/cm^2

%% State variables
v_TC = -65 * ones(1, nt);
v_TRN = -65 * ones(1, nt);

m_TC = zeros(1, nt); h_TC = zeros(1, nt); n_TC = zeros(1, nt);
m_TRN = zeros(1, nt); h_TRN = zeros(1, nt); n_TRN = zeros(1, nt);

[m_TC(1), h_TC(1), n_TC(1)] = init_hh_gates(v_TC(1));
[m_TRN(1), h_TRN(1), n_TRN(1)] = init_hh_gates(v_TRN(1));

%% Synaptic transmitter, receptor, and current arrays
nt_TC_to_TRN_AMPA = zeros(1, nt);
nt_TC_to_TRN_NMDA = zeros(1, nt);
nt_TRN_to_TC_GABAa = zeros(1, nt);
nt_TRN_to_TC_GABAb = zeros(1, nt);

r_AMPA_TRN = zeros(1, nt);
r_NMDA_TRN = zeros(1, nt);
B_NMDA_TRN = zeros(1, nt);

r_GABAa_TC = zeros(1, nt);
r_GABAb_TC = zeros(1, nt);
G_GABAb_TC = zeros(1, nt);
GABAb_open_TC = zeros(1, nt);

I_AMPA_TC_TRN = zeros(1, nt);
I_NMDA_TC_TRN = zeros(1, nt);
I_GABAa_TRN_TC = zeros(1, nt);
I_GABAb_TRN_TC = zeros(1, nt);
I_syn_TRN = zeros(1, nt);
I_syn_TC = zeros(1, nt);

TC_spike_events = false(1, nt);
TRN_spike_events = false(1, nt);
TC_release_events = false(1, nt);
TRN_release_events = false(1, nt);
TC_scheduled_release_events = event_mask(time, TC_release_times);
TRN_scheduled_release_events = event_mask(time, TRN_release_times);

tc_release_AMPA_steps = 0;
tc_release_NMDA_steps = 0;
trn_release_GABAa_steps = 0;
trn_release_GABAb_steps = 0;

AMPA_pulse_steps = max(1, round(Cdur_AMPA / dt));
NMDA_pulse_steps = max(1, round(Cdur_NMDA / dt));
GABAa_pulse_steps = max(1, round(Cdur_GABAa / dt));
GABAb_pulse_steps = max(1, round(Cdur_GABAb / dt));

%% Euler method
for k = 1:nt-1
    % Track spikes separately from transmitter release. In the default
    % visualization protocol, release is scheduled so excitation and
    % inhibition are staggered; disabling the protocol restores purely
    % spike-triggered release.
    TC_spike = k > 1 && v_TC(k-1) < spike_threshold && v_TC(k) >= spike_threshold;
    TRN_spike = k > 1 && v_TRN(k-1) < spike_threshold && v_TRN(k) >= spike_threshold;

    TC_spike_events(k) = TC_spike;
    TRN_spike_events(k) = TRN_spike;

    TC_release_triggered = (use_staggered_release_protocol && ...
        TC_scheduled_release_events(k)) || ...
        (~use_staggered_release_protocol && TC_spike);
    TRN_release_triggered = (use_staggered_release_protocol && ...
        TRN_scheduled_release_events(k)) || ...
        (~use_staggered_release_protocol && TRN_spike);

    TC_release_events(k) = TC_release_triggered;
    TRN_release_events(k) = TRN_release_triggered;

    tc_release_AMPA_steps = max(tc_release_AMPA_steps, ...
        AMPA_pulse_steps * TC_release_triggered);
    tc_release_NMDA_steps = max(tc_release_NMDA_steps, ...
        NMDA_pulse_steps * TC_release_triggered);
    trn_release_GABAa_steps = max(trn_release_GABAa_steps, ...
        GABAa_pulse_steps * TRN_release_triggered);
    trn_release_GABAb_steps = max(trn_release_GABAb_steps, ...
        GABAb_pulse_steps * TRN_release_triggered);

    AMPA_pulse_active = tc_release_AMPA_steps > 0;
    NMDA_pulse_active = tc_release_NMDA_steps > 0;
    GABAa_pulse_active = trn_release_GABAa_steps > 0;
    GABAb_pulse_active = trn_release_GABAb_steps > 0;

    nt_TC_to_TRN_AMPA(k) = Cmax_AMPA * AMPA_pulse_active;
    nt_TC_to_TRN_NMDA(k) = Cmax_NMDA * NMDA_pulse_active;
    nt_TRN_to_TC_GABAa(k) = Cmax_GABAa * GABAa_pulse_active;
    nt_TRN_to_TC_GABAb(k) = Cmax_GABAb * GABAb_pulse_active;

    tc_release_AMPA_steps = max(tc_release_AMPA_steps - AMPA_pulse_active, 0);
    tc_release_NMDA_steps = max(tc_release_NMDA_steps - NMDA_pulse_active, 0);
    trn_release_GABAa_steps = max(trn_release_GABAa_steps - GABAa_pulse_active, 0);
    trn_release_GABAb_steps = max(trn_release_GABAb_steps - GABAb_pulse_active, 0);

    % TC -> TRN AMPA receptor open fraction.
    dAMPA_rdt = alpha_AMPA * nt_TC_to_TRN_AMPA(k) * (1 - r_AMPA_TRN(k)) ...
        - beta_AMPA * r_AMPA_TRN(k);
    r_AMPA_TRN(k+1) = clamp01(r_AMPA_TRN(k) + dt * dAMPA_rdt);

    % TC -> TRN NMDA receptor open fraction with voltage-dependent Mg block.
    dNMDA_rdt = alpha_NMDA * nt_TC_to_TRN_NMDA(k) * (1 - r_NMDA_TRN(k)) ...
        - beta_NMDA * r_NMDA_TRN(k);
    r_NMDA_TRN(k+1) = clamp01(r_NMDA_TRN(k) + dt * dNMDA_rdt);
    B_NMDA_TRN(k) = magnesium_block(v_TRN(k), mg_NMDA);

    % TRN -> TC GABA_A receptor open fraction.
    dGABAa_rdt = alpha_GABAa * nt_TRN_to_TC_GABAa(k) * (1 - r_GABAa_TC(k)) ...
        - beta_GABAa * r_GABAa_TC(k);
    r_GABAa_TC(k+1) = clamp01(r_GABAa_TC(k) + dt * dGABAa_rdt);

    % TRN -> TC GABA_B receptor and G-protein activation.
    dGABAb_rdt = K1_GABAb * nt_TRN_to_TC_GABAb(k) * (1 - r_GABAb_TC(k)) ...
        - K2_GABAb * r_GABAb_TC(k);
    r_GABAb_TC(k+1) = clamp01(r_GABAb_TC(k) + dt * dGABAb_rdt);

    dGABAb_Gdt = K3_GABAb * r_GABAb_TC(k) - K4_GABAb * G_GABAb_TC(k);
    G_GABAb_TC(k+1) = max(G_GABAb_TC(k) + dt * dGABAb_Gdt, 0);
    GABAb_open_TC(k) = gprotein_channel_open(G_GABAb_TC(k), n_GABAb, KD_GABAb);

    % Synaptic currents. Positive values depolarize the postsynaptic cell;
    % negative values hyperpolarize it.
    I_AMPA_TC_TRN(k) = g_AMPA_TC_TRN * r_AMPA_TRN(k) * (E_AMPA - v_TRN(k));
    I_NMDA_TC_TRN(k) = g_NMDA_TC_TRN * B_NMDA_TRN(k) * r_NMDA_TRN(k) * (E_NMDA - v_TRN(k));
    I_GABAa_TRN_TC(k) = g_GABAa_TRN_TC * r_GABAa_TC(k) * (E_GABAa - v_TC(k));
    I_GABAb_TRN_TC(k) = g_GABAb_TRN_TC * GABAb_open_TC(k) * (E_GABAb - v_TC(k));

    I_syn_TRN(k) = I_AMPA_TC_TRN(k) + I_NMDA_TC_TRN(k);
    I_syn_TC(k) = I_GABAa_TRN_TC(k) + I_GABAb_TRN_TC(k);

    % Intrinsic HH gating updates.
    [am, bm, ah, bh, an, bn] = hh_rates(v_TC(k));
    m_TC(k+1) = clamp01(m_TC(k) + dt * (am * (1 - m_TC(k)) - bm * m_TC(k)));
    h_TC(k+1) = clamp01(h_TC(k) + dt * (ah * (1 - h_TC(k)) - bh * h_TC(k)));
    n_TC(k+1) = clamp01(n_TC(k) + dt * (an * (1 - n_TC(k)) - bn * n_TC(k)));

    [am, bm, ah, bh, an, bn] = hh_rates(v_TRN(k));
    m_TRN(k+1) = clamp01(m_TRN(k) + dt * (am * (1 - m_TRN(k)) - bm * m_TRN(k)));
    h_TRN(k+1) = clamp01(h_TRN(k) + dt * (ah * (1 - h_TRN(k)) - bh * h_TRN(k)));
    n_TRN(k+1) = clamp01(n_TRN(k) + dt * (an * (1 - n_TRN(k)) - bn * n_TRN(k)));

    INa_TC = gNa_TC * m_TC(k)^3 * h_TC(k) * (v_TC(k) - ENa_TC);
    IK_TC = gK_TC * n_TC(k)^4 * (v_TC(k) - EK_TC);
    IL_TC = gL_TC * (v_TC(k) - EL_TC);

    INa_TRN = gNa_TRN * m_TRN(k)^3 * h_TRN(k) * (v_TRN(k) - ENa_TRN);
    IK_TRN = gK_TRN * n_TRN(k)^4 * (v_TRN(k) - EK_TRN);
    IL_TRN = gL_TRN * (v_TRN(k) - EL_TRN);

    % Membrane voltage updates.
    dv_TC = (I_ext_TC(k) + I_syn_TC(k) - INa_TC - IK_TC - IL_TC) / Cm;
    dv_TRN = (I_ext_TRN(k) + I_syn_TRN(k) - INa_TRN - IK_TRN - IL_TRN) / Cm;

    v_TC(k+1) = v_TC(k) + dt * dv_TC;
    v_TRN(k+1) = v_TRN(k) + dt * dv_TRN;
end

%% Final point calculations
B_NMDA_TRN(end) = magnesium_block(v_TRN(end), mg_NMDA);
GABAb_open_TC(end) = gprotein_channel_open(G_GABAb_TC(end), n_GABAb, KD_GABAb);

I_AMPA_TC_TRN(end) = g_AMPA_TC_TRN * r_AMPA_TRN(end) * (E_AMPA - v_TRN(end));
I_NMDA_TC_TRN(end) = g_NMDA_TC_TRN * B_NMDA_TRN(end) * r_NMDA_TRN(end) * (E_NMDA - v_TRN(end));
I_GABAa_TRN_TC(end) = g_GABAa_TRN_TC * r_GABAa_TC(end) * (E_GABAa - v_TC(end));
I_GABAb_TRN_TC(end) = g_GABAb_TRN_TC * GABAb_open_TC(end) * (E_GABAb - v_TC(end));
I_syn_TRN(end) = I_AMPA_TC_TRN(end) + I_NMDA_TC_TRN(end);
I_syn_TC(end) = I_GABAa_TRN_TC(end) + I_GABAb_TRN_TC(end);

%% Summary in command window
fprintf('TC spikes:  %d\n', sum(TC_spike_events));
fprintf('TRN spikes: %d\n', sum(TRN_spike_events));
fprintf('TC -> TRN release pulses:  %d\n', sum(TC_release_events));
fprintf('TRN -> TC release pulses: %d\n', sum(TRN_release_events));
fprintf('Peak TC -> TRN AMPA current:  %.3f\n', max(I_AMPA_TC_TRN));
fprintf('Peak TC -> TRN NMDA current:  %.3f\n', max(I_NMDA_TC_TRN));
fprintf('Peak TRN -> TC GABA_A current: %.3f\n', min(I_GABAa_TRN_TC));
fprintf('Peak TRN -> TC GABA_B current: %.3f\n', min(I_GABAb_TRN_TC));

%% Plot 1: membrane potentials and external drive
figure('Name', 'TRN-Thalamic Microcircuit: Voltages and Drive', 'Color', 'w');

subplot(3, 1, 1)
plot(time, I_ext_TC, 'k', 'LineWidth', 1.2, 'DisplayName', 'TC drive')
hold on
plot(time, I_ext_TRN, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'DisplayName', 'TRN drive')
title('External input currents')
xlabel('Time (ms)')
ylabel('I_{ext}')
legend('Location', 'best')
grid on

subplot(3, 1, 2)
plot(time, v_TC, 'r', 'LineWidth', 1.2, 'DisplayName', 'Thalamic TC neuron')
hold on
plot(time, v_TRN, 'b', 'LineWidth', 1.2, 'DisplayName', 'TRN neuron')
plot(time, spike_threshold * ones(size(time)), 'k:', 'DisplayName', 'Release threshold')
title('Membrane potentials')
xlabel('Time (ms)')
ylabel('V_m (mV)')
legend('Location', 'best')
grid on

subplot(3, 1, 3)
stem(time(TC_release_events), ones(1, sum(TC_release_events)), ...
    'r', 'Marker', 'none', 'DisplayName', 'TC -> TRN release')
hold on
stem(time(TRN_release_events), -ones(1, sum(TRN_release_events)), ...
    'b', 'Marker', 'none', 'DisplayName', 'TRN -> TC release')
title('Staggered transmitter release events')
xlabel('Time (ms)')
ylabel('Release event')
ylim([-1.5 1.5])
legend('Location', 'best')
grid on

%% Plot 2: neurotransmitter pulses
figure('Name', 'TRN-Thalamic Microcircuit: Neurotransmitter Pulses', 'Color', 'w');

subplot(2, 1, 1)
plot(time, nt_TC_to_TRN_AMPA, 'r', 'LineWidth', 1.1, 'DisplayName', 'AMPA pulse')
hold on
plot(time, nt_TC_to_TRN_NMDA, 'm--', 'LineWidth', 1.1, 'DisplayName', 'NMDA pulse')
title('TC -> TRN glutamate pulses')
xlabel('Time (ms)')
ylabel('[NT]')
legend('Location', 'best')
grid on

subplot(2, 1, 2)
plot(time, nt_TRN_to_TC_GABAa, 'b', 'LineWidth', 1.1, 'DisplayName', 'GABA_A pulse')
hold on
plot(time, nt_TRN_to_TC_GABAb, 'Color', [0.3 0.1 0.7], ...
    'LineStyle', '--', 'LineWidth', 1.1, 'DisplayName', 'GABA_B pulse')
title('TRN -> TC GABA pulses')
xlabel('Time (ms)')
ylabel('[NT]')
legend('Location', 'best')
grid on

%% Plot 3: receptor and channel states
figure('Name', 'TRN-Thalamic Microcircuit: Receptor States', 'Color', 'w');

subplot(2, 1, 1)
plot(time, r_AMPA_TRN, 'r', 'LineWidth', 1.2, 'DisplayName', 'TRN AMPA open fraction')
hold on
plot(time, r_NMDA_TRN, 'm', 'LineWidth', 1.2, 'DisplayName', 'TRN NMDA open fraction')
plot(time, B_NMDA_TRN, 'k--', 'LineWidth', 1.0, 'DisplayName', 'NMDA Mg block factor')
title('Excitatory receptor state on TRN neuron')
xlabel('Time (ms)')
ylabel('State')
legend('Location', 'best')
grid on

subplot(2, 1, 2)
plot(time, r_GABAa_TC, 'b', 'LineWidth', 1.2, 'DisplayName', 'TC GABA_A open fraction')
hold on
plot(time, r_GABAb_TC, 'Color', [0.3 0.1 0.7], ...
    'LineWidth', 1.2, 'DisplayName', 'TC GABA_B receptor activation')
plot(time, GABAb_open_TC, 'k--', 'LineWidth', 1.0, 'DisplayName', 'GABA_B K channel open')
title('Inhibitory receptor state on TC neuron')
xlabel('Time (ms)')
ylabel('State')
legend('Location', 'best')
grid on

%% Plot 4: synaptic currents
figure('Name', 'TRN-Thalamic Microcircuit: Synaptic Currents', 'Color', 'w');

subplot(3, 1, 1)
plot(time, I_AMPA_TC_TRN, 'r', 'LineWidth', 1.2, 'DisplayName', 'AMPA')
hold on
plot(time, I_NMDA_TC_TRN, 'm', 'LineWidth', 1.2, 'DisplayName', 'NMDA')
plot(time, I_syn_TRN, 'k', 'LineWidth', 1.3, 'DisplayName', 'Total TC -> TRN')
title('Excitatory current into TRN neuron')
xlabel('Time (ms)')
ylabel('Current')
legend('Location', 'best')
grid on

subplot(3, 1, 2)
plot(time, I_GABAa_TRN_TC, 'b', 'LineWidth', 1.2, 'DisplayName', 'GABA_A')
hold on
plot(time, I_GABAb_TRN_TC, 'Color', [0.3 0.1 0.7], ...
    'LineWidth', 1.2, 'DisplayName', 'GABA_B')
plot(time, I_syn_TC, 'k', 'LineWidth', 1.3, 'DisplayName', 'Total TRN -> TC')
plot(time, zeros(size(time)), 'k:')
title('Inhibitory current into TC neuron')
xlabel('Time (ms)')
ylabel('Current')
legend('Location', 'best')
grid on

subplot(3, 1, 3)
plot(time, I_syn_TRN, 'Color', [0.8 0.1 0.1], 'LineWidth', 1.2, ...
    'DisplayName', 'Net excitatory input to TRN')
hold on
plot(time, I_syn_TC, 'Color', [0.1 0.1 0.8], 'LineWidth', 1.2, ...
    'DisplayName', 'Net inhibitory input to TC')
plot(time, zeros(size(time)), 'k:')
title('Reciprocal microcircuit currents')
xlabel('Time (ms)')
ylabel('Current')
legend('Location', 'best')
grid on

%% Plot 5: phase-plane view of coupled voltages
figure('Name', 'TRN-Thalamic Microcircuit: Coupled Voltage State', 'Color', 'w');
plot(v_TC, v_TRN, 'k', 'LineWidth', 1.0)
hold on
plot(v_TC(1), v_TRN(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start')
plot(v_TC(end), v_TRN(end), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'End')
title('TC and TRN voltage trajectory')
xlabel('TC V_m (mV)')
ylabel('TRN V_m (mV)')
legend('Location', 'best')
grid on

%% Helper functions
function [m0, h0, n0] = init_hh_gates(V)
    [am, bm, ah, bh, an, bn] = hh_rates(V);
    m0 = am / (am + bm);
    h0 = ah / (ah + bh);
    n0 = an / (an + bn);
end

function [am, bm, ah, bh, an, bn] = hh_rates(V)
    am = alpha_m(V); bm = beta_m(V);
    ah = alpha_h(V); bh = beta_h(V);
    an = alpha_n(V); bn = beta_n(V);
end

function val = alpha_m(V)
    x = V + 40;
    if abs(x) < 1e-6
        val = 1.0;
    else
        val = 0.1 * x / (1 - exp(-x / 10));
    end
end

function val = beta_m(V)
    val = 4 * exp(-(V + 65) / 18);
end

function val = alpha_h(V)
    val = 0.07 * exp(-(V + 65) / 20);
end

function val = beta_h(V)
    val = 1 / (1 + exp(-(V + 35) / 10));
end

function val = alpha_n(V)
    x = V + 55;
    if abs(x) < 1e-6
        val = 0.1;
    else
        val = 0.01 * x / (1 - exp(-x / 10));
    end
end

function val = beta_n(V)
    val = 0.125 * exp(-(V + 65) / 80);
end

function val = clamp01(x)
    val = min(max(x, 0), 1);
end

function events = event_mask(time, event_times)
    events = false(size(time));

    for event_time = event_times
        if event_time >= time(1) && event_time <= time(end)
            [~, event_idx] = min(abs(time - event_time));
            events(event_idx) = true;
        end
    end
end

function val = magnesium_block(V, mg)
    val = 1 / (1 + exp(0.062 * (-V)) * (mg / 3.57));
end

function val = gprotein_channel_open(G, n_sites, KD)
    val = (G^n_sites) / (G^n_sites + KD);
end
