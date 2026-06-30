%% Corticothalamic 3-Neuron Microcircuit with T-type Ca2+ Channels on TC
% Demonstrates rebound spiking / anode-break excitation on the TC relay neuron
% after TRN-driven inhibition, using T-type calcium and H-current kinetics.
%
% Combines:
%   - 3-neuron corticothalamic microcircuit (Ctx, TC, TRN)
%   - T-type Ca2+ (I_T) and H-current (I_H) on the TC neuron
%
% Protocol: prolonged TRN bias drives TRN spiking -> GABAergic inhibition of TC;
%           TC rebound upon release is mediated by I_T (no post-inhibitory TC drive).

clear; clc

%% Simulation setup
dt = 0.01; % ms
t_end = 350; % ms
time = 0:dt:t_end;
nt = numel(time);

%% Intrinsic membrane parameters
Cm = 1.0; % uF/cm^2

% Cortical (Ctx) neuron
gNa_Ctx = 100; ENa_Ctx = 50;
gK_Ctx = 30;  EK_Ctx = -77;
gL_Ctx = 0.1;  EL_Ctx = -65;

% TC relay neuron (HH + T-type Ca2+ + H-current)
gNa_TC = 120; ENa_TC = 50;
gK_TC = 36;  EK_TC = -77;
gL_TC = 0.3;  EL_TC = -54.4;
gT_TC = 1.2;  ECa = 120;   % mS/cm^2, mV — T-type Ca2+ (slightly above isolated demo)
gH_TC = 0.05; EH = -43;    % mS/cm^2, mV — H-current supports rebound with I_T

% TRN neuron
gNa_TRN = 100; ENa_TRN = 50;
gK_TRN = 30;  EK_TRN = -77;
gL_TRN = 0.25; EL_TRN = -58;

spike_threshold = 0; % mV

%% External drive
I_ext_TC = zeros(1, nt);
I_ext_TRN = zeros(1, nt);
I_ext_Ctx = zeros(1, nt);

% Prolonged TRN drive (deep hyperpolarization of TC via GABA from TRN spikes)
I_ext_TRN(time >= 150 & time < 180) = 10;

% Rebound TC drive intentionally disabled — observe pure I_T-mediated rebound
% I_ext_TC(time >= 180 & time < 300) = 8;

use_staggered_release_protocol = false;
TC_release_times = []; TRN_release_times = []; Ctx_release_times = [];

%% Synaptic constants
Cmax_AMPA = 1; Cdur_AMPA = 1; alpha_AMPA = 1.1; beta_AMPA = 0.19; E_AMPA = 0;
Cmax_NMDA = 1; Cdur_NMDA = 1; alpha_NMDA = 0.072; beta_NMDA = 0.0066; E_NMDA = 0; mg_NMDA = 1;
Cmax_GABAa = 1; Cdur_GABAa = 1; alpha_GABAa = 5; beta_GABAa = 0.30; E_GABAa = -80;
% beta_GABAa increased vs standard (0.18) so IPSCs decay quickly enough for I_T rebound
Cmax_GABAb = 1; Cdur_GABAb = 1; K1_GABAb = 0.09; K2_GABAb = 0.0012; K3_GABAb = 0.18; K4_GABAb = 0.034; KD_GABAb = 0.08; n_GABAb = 4; E_GABAb = -95;

% Synaptic max conductances (mS/cm^2)
g_AMPA_TC_TRN  = 0.18; g_NMDA_TC_TRN  = 0.12;
g_AMPA_TC_Ctx  = 0.15; g_NMDA_TC_Ctx  = 0.10;
g_AMPA_Ctx_TC  = 0.16; g_NMDA_Ctx_TC  = 0.11;
g_AMPA_Ctx_TRN = 0.14; g_NMDA_Ctx_TRN = 0.09;
g_GABAa_TRN_TC = 0.35; g_GABAb_TRN_TC = 0.00;
% GABA_A tuned for deep TC hyperpolarization; GABA_B disabled (slow GABAB masks rebound)

%% State variables
v_Ctx = -65 * ones(1, nt); v_TC = -65 * ones(1, nt); v_TRN = -65 * ones(1, nt);
m_Ctx = zeros(1, nt); h_Ctx = zeros(1, nt); n_Ctx = zeros(1, nt);
m_TC  = zeros(1, nt); h_TC  = zeros(1, nt); n_TC  = zeros(1, nt);
m_TRN = zeros(1, nt); h_TRN = zeros(1, nt); n_TRN = zeros(1, nt);

% T-type inactivation (h_T) and H-current activation (r) on TC
ht_TC = zeros(1, nt);
r_TC = zeros(1, nt);

[m_Ctx(1), h_Ctx(1), n_Ctx(1)] = init_hh_gates(v_Ctx(1));
[m_TC(1), h_TC(1), n_TC(1)]   = init_hh_gates(v_TC(1));
[m_TRN(1), h_TRN(1), n_TRN(1)] = init_hh_gates(v_TRN(1));
ht_TC(1) = hInf_T(v_TC(1));
r_TC(1) = rInf_H(v_TC(1));

%% Synaptic arrays
nt_TC_AMPA = zeros(1, nt);  nt_TC_NMDA = zeros(1, nt);
nt_Ctx_AMPA = zeros(1, nt); nt_Ctx_NMDA = zeros(1, nt);
nt_TRN_GABAa = zeros(1, nt); nt_TRN_GABAb = zeros(1, nt);

r_AMPA_TRN_from_TC = zeros(1, nt); r_NMDA_TRN_from_TC = zeros(1, nt);
r_AMPA_Ctx_from_TC = zeros(1, nt); r_NMDA_Ctx_from_TC = zeros(1, nt);
r_AMPA_TC_from_Ctx = zeros(1, nt); r_NMDA_TC_from_Ctx = zeros(1, nt);
r_AMPA_TRN_from_Ctx = zeros(1, nt); r_NMDA_TRN_from_Ctx = zeros(1, nt);
r_GABAa_TC = zeros(1, nt); r_GABAb_TC = zeros(1, nt); G_GABAb_TC = zeros(1, nt);
B_NMDA_TRN = zeros(1, nt); B_NMDA_TC = zeros(1, nt); B_NMDA_Ctx = zeros(1, nt);

I_syn_Ctx = zeros(1, nt); I_syn_TC = zeros(1, nt); I_syn_TRN = zeros(1, nt);
I_T_TC = zeros(1, nt); I_H_TC = zeros(1, nt);
I_GABAa_TRN_TC = zeros(1, nt); I_GABAb_TRN_TC = zeros(1, nt);

Ctx_spike_events = zeros(1, nt); TC_spike_events = zeros(1, nt); TRN_spike_events = zeros(1, nt);
Ctx_release_events = zeros(1, nt); TC_release_events = zeros(1, nt); TRN_release_events = zeros(1, nt);
Ctx_sched_events = zeros(1, nt); TC_sched_events = zeros(1, nt); TRN_sched_events = zeros(1, nt);

tc_amp_steps = 0; tc_nmd_steps = 0; ctx_amp_steps = 0; ctx_nmd_steps = 0; trn_ga_steps = 0; trn_gb_steps = 0;
AMPA_pulse_steps = max(1, round(Cdur_AMPA / dt));
NMDA_pulse_steps = max(1, round(Cdur_NMDA / dt));
GABAa_pulse_steps = max(1, round(Cdur_GABAa / dt));
GABAb_pulse_steps = max(1, round(Cdur_GABAb / dt));

%% Simulation loop
for k = 1:nt-1
    Ctx_spike = k > 1 && v_Ctx(k-1) < spike_threshold && v_Ctx(k) >= spike_threshold;
    TC_spike  = k > 1 && v_TC(k-1)  < spike_threshold && v_TC(k)  >= spike_threshold;
    TRN_spike = k > 1 && v_TRN(k-1) < spike_threshold && v_TRN(k) >= spike_threshold;

    if Ctx_spike, Ctx_spike_events(k) = 1; end
    if TC_spike,  TC_spike_events(k)  = 1; end
    if TRN_spike, TRN_spike_events(k) = 1; end

    if use_staggered_release_protocol
        Ctx_rel = Ctx_sched_events(k) == 1;
        TC_rel = TC_sched_events(k) == 1;
        TRN_rel = TRN_sched_events(k) == 1;
    else
        Ctx_rel = Ctx_spike; TC_rel = TC_spike; TRN_rel = TRN_spike;
    end

    if Ctx_rel
        Ctx_release_events(k) = 1;
        ctx_amp_steps = AMPA_pulse_steps; ctx_nmd_steps = NMDA_pulse_steps;
    end
    if TC_rel
        TC_release_events(k) = 1;
        tc_amp_steps = AMPA_pulse_steps; tc_nmd_steps = NMDA_pulse_steps;
    end
    if TRN_rel
        TRN_release_events(k) = 1;
        trn_ga_steps = GABAa_pulse_steps; trn_gb_steps = GABAb_pulse_steps;
    end

    if ctx_amp_steps > 0, nt_Ctx_AMPA(k) = Cmax_AMPA; ctx_amp_steps = ctx_amp_steps - 1; end
    if ctx_nmd_steps > 0, nt_Ctx_NMDA(k) = Cmax_NMDA; ctx_nmd_steps = ctx_nmd_steps - 1; end
    if tc_amp_steps  > 0, nt_TC_AMPA(k)  = Cmax_AMPA; tc_amp_steps  = tc_amp_steps - 1;  end
    if tc_nmd_steps  > 0, nt_TC_NMDA(k)  = Cmax_NMDA; tc_nmd_steps  = tc_nmd_steps - 1;  end
    % Sustained GABA release while TRN is externally driven (high-rate firing proxy)
    if I_ext_TRN(k) > 0
        nt_TRN_GABAa(k) = Cmax_GABAa;
        nt_TRN_GABAb(k) = Cmax_GABAb;
    end

    if trn_ga_steps  > 0, nt_TRN_GABAa(k) = Cmax_GABAa; trn_ga_steps = trn_ga_steps - 1; end
    if trn_gb_steps  > 0, nt_TRN_GABAb(k) = Cmax_GABAb; trn_gb_steps = trn_gb_steps - 1; end

    %% Synaptic kinetics
    r_AMPA_TRN_from_TC(k+1) = clamp01(r_AMPA_TRN_from_TC(k) + dt * (alpha_AMPA * nt_TC_AMPA(k) * (1 - r_AMPA_TRN_from_TC(k)) - beta_AMPA * r_AMPA_TRN_from_TC(k)));
    r_NMDA_TRN_from_TC(k+1) = clamp01(r_NMDA_TRN_from_TC(k) + dt * (alpha_NMDA * nt_TC_NMDA(k) * (1 - r_NMDA_TRN_from_TC(k)) - beta_NMDA * r_NMDA_TRN_from_TC(k)));
    r_AMPA_Ctx_from_TC(k+1) = clamp01(r_AMPA_Ctx_from_TC(k) + dt * (alpha_AMPA * nt_TC_AMPA(k) * (1 - r_AMPA_Ctx_from_TC(k)) - beta_AMPA * r_AMPA_Ctx_from_TC(k)));
    r_NMDA_Ctx_from_TC(k+1) = clamp01(r_NMDA_Ctx_from_TC(k) + dt * (alpha_NMDA * nt_TC_NMDA(k) * (1 - r_NMDA_Ctx_from_TC(k)) - beta_NMDA * r_NMDA_Ctx_from_TC(k)));
    r_AMPA_TC_from_Ctx(k+1) = clamp01(r_AMPA_TC_from_Ctx(k) + dt * (alpha_AMPA * nt_Ctx_AMPA(k) * (1 - r_AMPA_TC_from_Ctx(k)) - beta_AMPA * r_AMPA_TC_from_Ctx(k)));
    r_NMDA_TC_from_Ctx(k+1) = clamp01(r_NMDA_TC_from_Ctx(k) + dt * (alpha_NMDA * nt_Ctx_NMDA(k) * (1 - r_NMDA_TC_from_Ctx(k)) - beta_NMDA * r_NMDA_TC_from_Ctx(k)));
    r_AMPA_TRN_from_Ctx(k+1) = clamp01(r_AMPA_TRN_from_Ctx(k) + dt * (alpha_AMPA * nt_Ctx_AMPA(k) * (1 - r_AMPA_TRN_from_Ctx(k)) - beta_AMPA * r_AMPA_TRN_from_Ctx(k)));
    r_NMDA_TRN_from_Ctx(k+1) = clamp01(r_NMDA_TRN_from_Ctx(k) + dt * (alpha_NMDA * nt_Ctx_NMDA(k) * (1 - r_NMDA_TRN_from_Ctx(k)) - beta_NMDA * r_NMDA_TRN_from_Ctx(k)));
    r_GABAa_TC(k+1) = clamp01(r_GABAa_TC(k) + dt * (alpha_GABAa * nt_TRN_GABAa(k) * (1 - r_GABAa_TC(k)) - beta_GABAa * r_GABAa_TC(k)));
    r_GABAb_TC(k+1) = clamp01(r_GABAb_TC(k) + dt * (K1_GABAb * nt_TRN_GABAb(k) * (1 - r_GABAb_TC(k)) - K2_GABAb * r_GABAb_TC(k)));
    G_GABAb_TC(k+1) = max(G_GABAb_TC(k) + dt * (K3_GABAb * r_GABAb_TC(k) - K4_GABAb * G_GABAb_TC(k)), 0);

    B_NMDA_TRN(k) = magnesium_block(v_TRN(k), mg_NMDA);
    B_NMDA_TC(k)  = magnesium_block(v_TC(k), mg_NMDA);
    B_NMDA_Ctx(k) = magnesium_block(v_Ctx(k), mg_NMDA);

    %% T-type and H-current gating (TC only)
    ht_TC(k+1) = ht_TC(k) + dt * ((hInf_T(v_TC(k)) - ht_TC(k)) / tauH_T(v_TC(k)));
    r_TC(k+1) = rInf_H(v_TC(k)) + (r_TC(k) - rInf_H(v_TC(k))) * exp(-dt / tauR_H(v_TC(k)));

    %% Currents
    I_AMPA_TC_Ctx = g_AMPA_TC_Ctx * r_AMPA_Ctx_from_TC(k) * (E_AMPA - v_Ctx(k));
    I_NMDA_TC_Ctx = g_NMDA_TC_Ctx * B_NMDA_Ctx(k) * r_NMDA_Ctx_from_TC(k) * (E_NMDA - v_Ctx(k));
    I_syn_Ctx(k)  = I_AMPA_TC_Ctx + I_NMDA_TC_Ctx;

    I_AMPA_Ctx_TC  = g_AMPA_Ctx_TC * r_AMPA_TC_from_Ctx(k) * (E_AMPA - v_TC(k));
    I_NMDA_Ctx_TC  = g_NMDA_Ctx_TC * B_NMDA_TC(k) * r_NMDA_TC_from_Ctx(k) * (E_NMDA - v_TC(k));
    I_GABAa_TRN_TC(k) = g_GABAa_TRN_TC * r_GABAa_TC(k) * (E_GABAa - v_TC(k));
    I_GABAb_TRN_TC(k) = g_GABAb_TRN_TC * gprotein_channel_open(G_GABAb_TC(k), n_GABAb, KD_GABAb) * (E_GABAb - v_TC(k));
    I_syn_TC(k) = I_AMPA_Ctx_TC + I_NMDA_Ctx_TC + I_GABAa_TRN_TC(k) + I_GABAb_TRN_TC(k);

    I_AMPA_TC_TRN  = g_AMPA_TC_TRN * r_AMPA_TRN_from_TC(k) * (E_AMPA - v_TRN(k));
    I_NMDA_TC_TRN  = g_NMDA_TC_TRN * B_NMDA_TRN(k) * r_NMDA_TRN_from_TC(k) * (E_NMDA - v_TRN(k));
    I_AMPA_Ctx_TRN = g_AMPA_Ctx_TRN * r_AMPA_TRN_from_Ctx(k) * (E_AMPA - v_TRN(k));
    I_NMDA_Ctx_TRN = g_NMDA_Ctx_TRN * B_NMDA_TRN(k) * r_NMDA_TRN_from_Ctx(k) * (E_NMDA - v_TRN(k));
    I_syn_TRN(k)   = I_AMPA_TC_TRN + I_NMDA_TC_TRN + I_AMPA_Ctx_TRN + I_NMDA_Ctx_TRN;

    I_T_TC(k) = gT_TC * mInf_T(v_TC(k))^2 * ht_TC(k) * (ECa - v_TC(k));
    I_H_TC(k) = gH_TC * r_TC(k) * (EH - v_TC(k));

    %% Hodgkin-Huxley updates
    [am, bm, ah, bh, an, bn] = hh_rates(v_Ctx(k));
    m_Ctx(k+1) = clamp01(m_Ctx(k) + dt * (am * (1 - m_Ctx(k)) - bm * m_Ctx(k)));
    h_Ctx(k+1) = clamp01(h_Ctx(k) + dt * (ah * (1 - h_Ctx(k)) - bh * h_Ctx(k)));
    n_Ctx(k+1) = clamp01(n_Ctx(k) + dt * (an * (1 - n_Ctx(k)) - bn * n_Ctx(k)));
    INa_Ctx = gNa_Ctx * m_Ctx(k)^3 * h_Ctx(k) * (v_Ctx(k) - ENa_Ctx);
    IK_Ctx  = gK_Ctx * n_Ctx(k)^4 * (v_Ctx(k) - EK_Ctx);
    IL_Ctx  = gL_Ctx * (v_Ctx(k) - EL_Ctx);

    [am, bm, ah, bh, an, bn] = hh_rates(v_TC(k));
    m_TC(k+1) = clamp01(m_TC(k) + dt * (am * (1 - m_TC(k)) - bm * m_TC(k)));
    h_TC(k+1) = clamp01(h_TC(k) + dt * (ah * (1 - h_TC(k)) - bh * h_TC(k)));
    n_TC(k+1) = clamp01(n_TC(k) + dt * (an * (1 - n_TC(k)) - bn * n_TC(k)));
    INa_TC = gNa_TC * m_TC(k)^3 * h_TC(k) * (v_TC(k) - ENa_TC);
    IK_TC  = gK_TC * n_TC(k)^4 * (v_TC(k) - EK_TC);
    IL_TC  = gL_TC * (v_TC(k) - EL_TC);

    [am, bm, ah, bh, an, bn] = hh_rates(v_TRN(k));
    m_TRN(k+1) = clamp01(m_TRN(k) + dt * (am * (1 - m_TRN(k)) - bm * m_TRN(k)));
    h_TRN(k+1) = clamp01(h_TRN(k) + dt * (ah * (1 - h_TRN(k)) - bh * h_TRN(k)));
    n_TRN(k+1) = clamp01(n_TRN(k) + dt * (an * (1 - n_TRN(k)) - bn * n_TRN(k)));
    INa_TRN = gNa_TRN * m_TRN(k)^3 * h_TRN(k) * (v_TRN(k) - ENa_TRN);
    IK_TRN  = gK_TRN * n_TRN(k)^4 * (v_TRN(k) - EK_TRN);
    IL_TRN  = gL_TRN * (v_TRN(k) - EL_TRN);

    %% Voltage integration (I_T and I_H add depolarizing drive to TC)
    v_Ctx(k+1) = v_Ctx(k) + dt * ((I_ext_Ctx(k) + I_syn_Ctx(k) - INa_Ctx - IK_Ctx - IL_Ctx) / Cm);
    v_TC(k+1)  = v_TC(k)  + dt * ((I_ext_TC(k)  + I_syn_TC(k)  + I_T_TC(k)  + I_H_TC(k)  - INa_TC  - IK_TC  - IL_TC)  / Cm);
    v_TRN(k+1) = v_TRN(k) + dt * ((I_ext_TRN(k) + I_syn_TRN(k) - INa_TRN - IK_TRN - IL_TRN) / Cm);
end

%% Final-point alignments
B_NMDA_TRN(end) = magnesium_block(v_TRN(end), mg_NMDA);
B_NMDA_TC(end)  = magnesium_block(v_TC(end), mg_NMDA);
B_NMDA_Ctx(end) = magnesium_block(v_Ctx(end), mg_NMDA);
I_syn_Ctx(end)  = g_AMPA_TC_Ctx * r_AMPA_Ctx_from_TC(end) * (E_AMPA - v_Ctx(end)) + g_NMDA_TC_Ctx * B_NMDA_Ctx(end) * r_NMDA_Ctx_from_TC(end) * (E_NMDA - v_Ctx(end));
I_syn_TC(end)   = g_AMPA_Ctx_TC * r_AMPA_TC_from_Ctx(end) * (E_AMPA - v_TC(end)) + g_NMDA_Ctx_TC * B_NMDA_TC(end) * r_NMDA_TC_from_Ctx(end) * (E_NMDA - v_TC(end)) + g_GABAa_TRN_TC * r_GABAa_TC(end) * (E_GABAa - v_TC(end)) + g_GABAb_TRN_TC * gprotein_channel_open(G_GABAb_TC(end), n_GABAb, KD_GABAb) * (E_GABAb - v_TC(end));
I_syn_TRN(end)  = g_AMPA_TC_TRN * r_AMPA_TRN_from_TC(end) * (E_AMPA - v_TRN(end)) + g_NMDA_TC_TRN * B_NMDA_TRN(end) * r_NMDA_TRN_from_TC(end) * (E_NMDA - v_TRN(end)) + g_AMPA_Ctx_TRN * r_AMPA_TRN_from_Ctx(end) * (E_AMPA - v_TRN(end)) + g_NMDA_Ctx_TRN * B_NMDA_TRN(end) * r_NMDA_TRN_from_Ctx(end) * (E_NMDA - v_TRN(end));
I_T_TC(end) = gT_TC * mInf_T(v_TC(end))^2 * ht_TC(end) * (ECa - v_TC(end));
I_H_TC(end) = gH_TC * r_TC(end) * (EH - v_TC(end));
I_GABAa_TRN_TC(end) = g_GABAa_TRN_TC * r_GABAa_TC(end) * (E_GABAa - v_TC(end));
I_GABAb_TRN_TC(end) = g_GABAb_TRN_TC * gprotein_channel_open(G_GABAb_TC(end), n_GABAb, KD_GABAb) * (E_GABAb - v_TC(end));

%% Metrics
post_inhib_mask = time >= 180;
tc_rebound_spikes = sum(TC_spike_events & post_inhib_mask);
[~, peak_it_idx] = max(I_T_TC(post_inhib_mask));
post_inhib_times = time(post_inhib_mask);

fprintf('--- Spike Metrics ---\n');
fprintf('Ctx spikes: %d | TC spikes: %d | TRN spikes: %d\n', sum(Ctx_spike_events), sum(TC_spike_events), sum(TRN_spike_events));
fprintf('TC rebound spikes (t >= 180 ms): %d\n', tc_rebound_spikes);
fprintf('Peak TC I_T after inhibition: %.3f uA/cm^2 at t = %.1f ms\n', max(I_T_TC(post_inhib_mask)), post_inhib_times(peak_it_idx));
fprintf('Min TC voltage during TRN drive: %.1f mV\n', min(v_TC(time >= 150 & time < 200)));

%% Microcircuit figure
figure('Name', '3-Neuron Microcircuit — T-type Rebound', 'Color', 'w', 'Position', [100 100 800 650]);

subplot(3, 1, 1)
plot(time, I_ext_Ctx, 'g', 'LineWidth', 1.2, 'DisplayName', 'Ctx drive'); hold on
plot(time, I_ext_TC, 'r', 'LineWidth', 1.2, 'DisplayName', 'TC drive')
plot(time, I_ext_TRN, 'b', 'LineWidth', 1.2, 'DisplayName', 'TRN bias')
title('External Injection Currents')
xlabel('Time (ms)'); ylabel('I_{ext} (\muA/cm^2)'); legend('Location', 'best'); grid on

subplot(3, 1, 2)
plot(time, v_Ctx, 'g', 'LineWidth', 1.2, 'DisplayName', 'Cortical (Ctx)'); hold on
plot(time, v_TC, 'r', 'LineWidth', 1.2, 'DisplayName', 'Thalamic (TC)')
plot(time, v_TRN, 'b', 'LineWidth', 1.2, 'DisplayName', 'Reticular (TRN)')
line([0 t_end], [spike_threshold spike_threshold], 'Color', 'k', 'LineStyle', ':', 'DisplayName', 'Threshold')
xline(150, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'HandleVisibility', 'off')
xline(180, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'HandleVisibility', 'off')
title('Membrane Potentials (TC rebound after TRN inhibition release)')
xlabel('Time (ms)'); ylabel('V_m (mV)'); legend('Location', 'best'); grid on

subplot(3, 1, 3)
stem(time(Ctx_release_events == 1), 1.2 * ones(1, sum(Ctx_release_events)), 'g', 'Marker', 'none', 'LineWidth', 1.5, 'DisplayName', 'Ctx Release')
hold on
stem(time(TC_release_events == 1), 1.0 * ones(1, sum(TC_release_events)), 'r', 'Marker', 'none', 'LineWidth', 1.5, 'DisplayName', 'TC Release')
stem(time(TRN_release_events == 1), -1.0 * ones(1, sum(TRN_release_events)), 'b', 'Marker', 'none', 'LineWidth', 1.5, 'DisplayName', 'TRN Release')
title('Neurotransmitter Release Events')
xlabel('Time (ms)'); ylabel('Event Pathway'); ylim([-1.5 1.5]); legend('Location', 'best'); grid on

%% T-type rebound mechanism figure (mirrors I_T demo panels for TC)
figure('Name', 'TC T-type Ca2+ Rebound Mechanism', 'Color', 'w', 'Position', [120 120 750 700]);

subplot(4, 1, 1)
plot(time, v_TC, 'k', 'LineWidth', 1.4)
ylabel('V (mV)')
title('TC membrane potential during TRN-driven inhibition and rebound')
grid on

subplot(4, 1, 2)
plot(time, I_T_TC, 'r', 'LineWidth', 1.4)
ylabel('I_T (\muA/cm^2)')
title('T-type calcium current on TC')
grid on

subplot(4, 1, 3)
plot(time, ht_TC, 'b', 'LineWidth', 1.4); hold on
plot(time, mInf_T(v_TC), 'r--', 'LineWidth', 1.2)
ylabel('Gate')
title('h_T (solid) and m_{\infty,T}(V) (dashed)')
legend('h_T', 'm_{\infty,T}', 'Location', 'best')
grid on

subplot(4, 1, 4)
plot(time, I_H_TC, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.2); hold on
plot(time, I_GABAa_TRN_TC + I_GABAb_TRN_TC, 'b', 'LineWidth', 1.0)
xlabel('Time (ms)')
ylabel('Current (\muA/cm^2)')
title('H-current (green) and net TRN GABA inhibition (blue)')
legend('I_H', 'I_{GABA,TRN\rightarrowTC}', 'Location', 'best')
grid on

%% Helper functions
function [m0, h0, n0] = init_hh_gates(V)
    [am, bm, ah, bh, an, bn] = hh_rates(V);
    m0 = am / (am + bm); h0 = ah / (ah + bh); n0 = an / (an + bn);
end

function [am, bm, ah, bh, an, bn] = hh_rates(V)
    am = alpha_m(V); bm = beta_m(V); ah = alpha_h(V); bh = beta_h(V); an = alpha_n(V); bn = beta_n(V);
end

function val = alpha_m(V), x = V + 40; if abs(x) < 1e-6, val = 1.0; else, val = 0.1 * x / (1 - exp(-x / 10)); end; end
function val = beta_m(V),  val = 4 * exp(-(V + 65) / 18); end
function val = alpha_h(V), val = 0.07 * exp(-(V + 65) / 20); end
function val = beta_h(V),  val = 1 / (1 + exp(-(V + 35) / 10)); end
function val = alpha_n(V), x = V + 55; if abs(x) < 1e-6, val = 0.1; else, val = 0.01 * x / (1 - exp(-x / 10)); end; end
function val = beta_n(V),  val = 0.125 * exp(-(V + 65) / 80); end
function val = clamp01(x), val = min(max(x, 0), 1); end
function val = magnesium_block(V, mg), val = 1 / (1 + exp(0.062 * (-V)) * (mg / 3.57)); end
function val = gprotein_channel_open(G, n_sites, KD), val = (G^n_sites) / (G^n_sites + KD); end

% T-type / H-current kinetics live in mInf_T.m, hInf_T.m, tauH_T.m, rInf_H.m, tauR_H.m
