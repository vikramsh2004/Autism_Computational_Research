%% Minimal HH Loop: TC <-> TRN and TC -> Cx (Longer AMPA EPSP)
% 1 neuron per region
% HH Na/K/leak intrinsic currents
% AMPA excitation and GABA_A inhibition
% AMPA kinetics are slowed to broaden EPSP shape.

%% Simulation setup
dt = 0.01; % ms
t_end = 400; % ms (extended from 200 ms)
time = 0:dt:t_end;
nt = numel(time);

%% Intrinsic membrane parameters
Cm = 1.0; % uF/cm^2

% TC
gNa_TC = 120; ENa_TC = 50;
gK_TC = 36; EK_TC = -77;
gL_TC = 0.3; EL_TC = -54.4;

% TRN
gNa_TRN = 120; ENa_TRN = 50;
gK_TRN = 36; EK_TRN = -77;
gL_TRN = 0.3; EL_TRN = -54.4;

% Cortex
gNa_Cx = 120; ENa_Cx = 50; %#ok<NASGU>
gK_Cx = 36; EK_Cx = -77; %#ok<NASGU>
gL_Cx = 0.3; EL_Cx = -54.4; %#ok<NASGU>

%% Synapses
% TC -> TRN excitatory AMPA
g_AMPA_TC_TRN = 0.20; % mS/cm^2
E_AMPA = 0; % mV

% TRN -> TC inhibitory GABA_A
g_GABA_TRN_TC = 0.40; % mS/cm^2
E_GABA = -80; % mV

% Synaptic kinetics
% Slower AMPA kinetics to produce a broader, longer-lasting EPSP.
tau_rise_AMPA = 1.0; % ms
tau_decay_AMPA = 6.0; % ms
alpha_AMPA = 1 / tau_rise_AMPA; % 1/ms
beta_AMPA = 1 / tau_decay_AMPA; % 1/ms

alpha_GABA = 5.0; % 1/ms
beta_GABA = 0.2; % 1/ms

spikeThresh = 0; % simple release threshold

%% State variables
% Membrane voltages
v_TC = -65 * ones(1, nt);
v_TRN = -65 * ones(1, nt);

% HH gating variables
m_TC = zeros(1, nt); h_TC = zeros(1, nt); n_TC = zeros(1, nt);
m_TRN = zeros(1, nt); h_TRN = zeros(1, nt); n_TRN = zeros(1, nt);

% Synaptic gating variables
s_AMPA_TC_TRN = zeros(1, nt);
s_GABA_TRN_TC = zeros(1, nt);

% Synaptic currents
I_AMPA_TC_TRN = zeros(1, nt);
I_GABA_TRN_TC = zeros(1, nt);

%% Initialize HH gates
[m_TC(1), h_TC(1), n_TC(1)] = init_gates_longampa(v_TC(1));
[m_TRN(1), h_TRN(1), n_TRN(1)] = init_gates_longampa(v_TRN(1));

%% External input
I_ext_TC = zeros(1, nt);
I_ext_TC(2000:5000) = 10; % brief ext pulse into TC

%% Main update loop
for k = 1:nt-1

   % Presynaptic release detection
   TC_spike = double(v_TC(k) > spikeThresh);
   TRN_spike = double(v_TRN(k) > spikeThresh);

   % Synapse updates
   % TC -> TRN AMPA
   ds = alpha_AMPA * TC_spike * (1 - s_AMPA_TC_TRN(k)) ...
     - beta_AMPA * s_AMPA_TC_TRN(k);
   s_AMPA_TC_TRN(k+1) = s_AMPA_TC_TRN(k) + dt * ds;
   s_AMPA_TC_TRN(k+1) = min(max(s_AMPA_TC_TRN(k+1), 0), 1);

   % TRN -> TC GABA_A
   ds = alpha_GABA * TRN_spike * (1 - s_GABA_TRN_TC(k)) ...
     - beta_GABA * s_GABA_TRN_TC(k);
   s_GABA_TRN_TC(k+1) = s_GABA_TRN_TC(k) + dt * ds;
   s_GABA_TRN_TC(k+1) = min(max(s_GABA_TRN_TC(k+1), 0), 1);

   % Synaptic currents
   I_AMPA_TC_TRN(k) = g_AMPA_TC_TRN * s_AMPA_TC_TRN(k) * (E_AMPA - v_TRN(k));
   I_GABA_TRN_TC(k) = g_GABA_TRN_TC * s_GABA_TRN_TC(k) * (E_GABA - v_TC(k));

   % Update HH gates: TC
   [am,bm,ah,bh,an,bn] = hh_rates_longampa(v_TC(k));
   m_TC(k+1) = m_TC(k) + dt * (am*(1-m_TC(k)) - bm*m_TC(k));
   h_TC(k+1) = h_TC(k) + dt * (ah*(1-h_TC(k)) - bh*h_TC(k));
   n_TC(k+1) = n_TC(k) + dt * (an*(1-n_TC(k)) - bn*n_TC(k));

   % Update HH gates: TRN
   [am,bm,ah,bh,an,bn] = hh_rates_longampa(v_TRN(k));
   m_TRN(k+1) = m_TRN(k) + dt * (am*(1-m_TRN(k)) - bm*m_TRN(k));
   h_TRN(k+1) = h_TRN(k) + dt * (ah*(1-h_TRN(k)) - bh*h_TRN(k));
   n_TRN(k+1) = n_TRN(k) + dt * (an*(1-n_TRN(k)) - bn*n_TRN(k));

   % Intrinsic currents
   INa_TC = gNa_TC * m_TC(k)^3 * h_TC(k) * (v_TC(k) - ENa_TC);
   IK_TC = gK_TC * n_TC(k)^4 * (v_TC(k) - EK_TC);
   IL_TC = gL_TC * (v_TC(k) - EL_TC);

   INa_TRN = gNa_TRN * m_TRN(k)^3 * h_TRN(k) * (v_TRN(k) - ENa_TRN);
   IK_TRN = gK_TRN * n_TRN(k)^4 * (v_TRN(k) - EK_TRN);
   IL_TRN = gL_TRN * (v_TRN(k) - EL_TRN);

   % Voltage updates
   dv_TC = (I_ext_TC(k) + I_GABA_TRN_TC(k) - INa_TC - IK_TC - IL_TC) / Cm;
   dv_TRN = (I_AMPA_TC_TRN(k) - INa_TRN - IK_TRN - IL_TRN) / Cm;
   v_TC(k+1) = v_TC(k) + dt * dv_TC;
   v_TRN(k+1) = v_TRN(k) + dt * dv_TRN;
end

% Current calculations
I_AMPA_TC_TRN(end) = g_AMPA_TC_TRN * s_AMPA_TC_TRN(end) * (E_AMPA - v_TRN(end));
I_GABA_TRN_TC(end) = g_GABA_TRN_TC * s_GABA_TRN_TC(end) * (E_GABA - v_TC(end));

%% Plotting
figure;
subplot(4,1,1)
plot(time, I_ext_TC, 'LineWidth', 1.2)
title('External input to TC')
xlabel('Time (ms)')
ylabel('I_{ext}')
grid on

subplot(4,1,2)
plot(time, v_TC, 'LineWidth', 1.3, 'DisplayName', 'TC')
hold on
plot(time, v_TRN, 'LineWidth', 1.2, 'DisplayName', 'TRN')
title('Membrane potentials')
xlabel('Time (ms)')
ylabel('v_m (mV)')
legend
grid on

subplot(4,1,3)
plot(time, I_AMPA_TC_TRN, 'LineWidth', 1.2, 'DisplayName', 'TC -> TRN AMPA')
hold on
plot(time, I_GABA_TRN_TC, 'LineWidth', 1.2, 'DisplayName', 'TRN -> TC GABA_A')
title('Synaptic currents')
xlabel('Time (ms)')
ylabel('Current')
legend
grid on

subplot(4,1,4)
plot(time, s_AMPA_TC_TRN, 'LineWidth', 1.2)
title('AMPA gating variable s_{AMPA} (broader kinetics)')
xlabel('Time (ms)')
ylabel('s_{AMPA}')
grid on

%% Helper functions
function [m0,h0,n0] = init_gates_longampa(V)
    [am,bm,ah,bh,an,bn] = hh_rates_longampa(V);
    m0 = am/(am+bm);
    h0 = ah/(ah+bh);
    n0 = an/(an+bn);
end

function [am,bm,ah,bh,an,bn] = hh_rates_longampa(V)
    am = alpha_m_longampa(V);
    bm = beta_m_longampa(V);
    ah = alpha_h_longampa(V);
    bh = beta_h_longampa(V);
    an = alpha_n_longampa(V);
    bn = beta_n_longampa(V);
end

function val = alpha_m_longampa(V)
    x = V + 40;
    if abs(x) < 1e-6
        val = 1.0;
    else
        val = 0.1*x/(1-exp(-x/10));
    end
end

function val = beta_m_longampa(V)
    val = 4*exp(-(V+65)/18);
end

function val = alpha_h_longampa(V)
    val = 0.07*exp(-(V+65)/20);
end

function val = beta_h_longampa(V)
    val = 1/(1+exp(-(V+35)/10));
end

function val = alpha_n_longampa(V)
    x = V + 55;
    if abs(x) < 1e-6
        val = 0.1;
    else
        val = 0.01*x/(1-exp(-x/10));
    end
end

function val = beta_n_longampa(V)
    val = 0.125*exp(-(V+65)/80);
end
