function out = L5_PFC_Pyramidal_Compartmental_Model(opts)
%% ========================================================================
%  L5 PFC PYRAMIDAL NEURON - 4-COMPARTMENT CONDUCTANCE-BASED MODEL
%  ========================================================================
%  Builds and simulates a Layer-5 prefrontal-cortex (PFC) pyramidal neuron
%  as a four-compartment cable, following the layout worked out in the
%  "L5 PFC Neuron Compartmental Layout" design document.
%
%  COMPARTMENTS (electrotonic tree, soma = central hub):
%
%        [ Basal Dendrites ]
%                |  g_c,BS
%            [ Soma ]
%                |  g_c,SP
%        [ Proximal Apical Trunk ]
%                |  g_c,PD
%        [ Distal Apical Tuft ]
%
%  PER-COMPARTMENT CHANNELS (each compartment carries ONLY its own set):
%
%    SOMA    : I_Na  (fast Na, spike generator)
%              I_Kdr (delayed rectifier K, repolarization)
%              I_M   (Kv7/KCNQ M-current, spike-frequency adaptation)
%              I_L   (leak)
%
%    BASAL   : I_Na  (low density, supports bAPs)
%              I_Kdr (low density)
%              I_A   (Kv4.2/4.3 A-type K, regulates EPSP summation)
%              I_L   (leak)
%              + AMPA / NMDA synapses  (mediodorsal-thalamic feedforward)
%
%    PROXIMAL: I_h    (HCN, sets dendritic Rin / resting potential)
%              I_CaL  (Cav1.2/1.3 L-type Ca, plateau potentials)
%              I_L    (leak)
%              + AMPA / NMDA synapses  (local / neuromodulatory site)
%              + intracellular Ca pool (fed by I_CaL)
%
%    DISTAL  : I_h    (HCN, HIGH density - increases distally)
%              I_CaHVA(Cav2.x P/Q/N/R-type HVA Ca, dendritic Ca spikes)
%              I_SK   (SK/KCa2 Ca-activated K, terminates Ca spikes / mAHP)
%              I_L    (leak)
%              + AMPA / NMDA synapses  (ACC top-down feedback)
%              + intracellular Ca pool (fed by I_CaHVA, gates I_SK)
%
%  Intrinsic channel gating kinetics (I_M, I_A, I_CaL, I_CaHVA, I_SK) use
%  the alpha/beta (and Hill, for SK) formulations developed one-by-one in
%  the design document. I_Na / I_Kdr use classic Hodgkin-Huxley rates and
%  I_h uses a Kole-type HCN activation.
%
%  NEUROMODULATION / PATHOLOGY HOOKS (see "mod" struct below) are exposed
%  as simple gain factors so VTA-dopamine (D1/D2) modulation and
%  SHANK3-ASD channelopathy (AMPA down-regulation, dendritic I_h loss,
%  mGluR5 uncoupling) can be layered on later without editing the core.
%
%  Units: V [mV], t [ms], g [mS/cm^2], I [uA/cm^2], C [uF/cm^2],
%         [Ca] [uM]. With g in mS/cm^2 and V in mV, g*V is in uA/cm^2,
%         so dV/dt = I/C is in mV/ms (self-consistent).
% =========================================================================

close all;
if nargin < 1 || isempty(opts), opts = struct(); end
rng(42);

% Compartment index constants
SOMA = 1; BAS = 2; PROX = 3; DIST = 4;
nComp = 4;
cName = {'Soma', 'Basal', 'Proximal', 'Distal'};

%% ---------------------- 1. SIMULATION PARAMETERS ------------------------
dt    = getdef(opts,'dt',    0.01);   % time step (ms)
t_max = getdef(opts,'t_max', 600);    % total time (ms)
t     = 0:dt:t_max;
N     = numel(t);
Cm    = 1.0;                  % membrane capacitance (uF/cm^2), all comps

%% ---------------------- 2. MORPHOLOGY & COUPLING ------------------------
% Effective-cylinder surface areas from the layout document (um^2).
area_um2 = [1570.8, 2199.1, 3769.9, 1885.0];   % [Soma Basal Prox Distal]
area_cm2 = area_um2 * 1e-8;                     % 1 um^2 = 1e-8 cm^2

% Inter-compartmental coupling conductances g_c (nS), from the axial-
% resistance calculation in the document (R_a = 150 ohm*cm):
%   Basal<->Soma      : 63.9 nS
%   Soma <->Proximal  : 23.5 nS
%   Proximal<->Distal :  8.77 nS
gc_nS = zeros(nComp);
gc_nS(SOMA,BAS)  = 63.9;  gc_nS(BAS,SOMA)  = 63.9;
gc_nS(SOMA,PROX) = 23.5;  gc_nS(PROX,SOMA) = 23.5;
gc_nS(PROX,DIST) = 8.77;  gc_nS(DIST,PROX) = 8.77;

% Convert to a per-compartment coupling DENSITY (mS/cm^2). The absolute
% axial current g_c*(Vj-Vi) is symmetric, but when expressed as a current
% density it is normalized by the *target* compartment's own area, so the
% density matrix is asymmetric (matches the table in the document, e.g.
% Basal<->Soma = 4.07 mS/cm^2 from the soma, 2.91 mS/cm^2 from the basal).
gc_den = zeros(nComp);
for i = 1:nComp
    for j = 1:nComp
        if gc_nS(i,j) > 0
            gc_den(i,j) = gc_nS(i,j) * 1e-6 / area_cm2(i);  % mS/cm^2
        end
    end
end

%% ---------------------- 3. REVERSAL POTENTIALS --------------------------
E_Na  =  50;    % mV
E_K   = -90;    % mV
E_Ca  = 120;    % mV
E_h   = -45;    % mV (HCN, mixed cation)
E_syn =   0;    % mV (AMPA/NMDA)

%% ---------------------- 4. MAXIMAL CONDUCTANCE DENSITIES ----------------
% Each row is a compartment; absent channels are simply left at 0.
g.Na   = zeros(1,nComp);  g.Kdr = zeros(1,nComp);  g.M   = zeros(1,nComp);
g.A    = zeros(1,nComp);  g.h   = zeros(1,nComp);  g.CaL = zeros(1,nComp);
g.CaHVA= zeros(1,nComp);  g.SK  = zeros(1,nComp);  g.L   = zeros(1,nComp);
E.L    = zeros(1,nComp);  % per-compartment leak reversal (sets rest)

% -- Soma: spike generator --
g.Na(SOMA)  = 120;  g.Kdr(SOMA) = 36;   g.M(SOMA) = 1.5;
g.L(SOMA)   = 0.10; E.L(SOMA)   = -70;

% -- Basal: low-density active + A-type --
g.Na(BAS)   = 20;   g.Kdr(BAS)  = 8;    g.A(BAS)  = 12;
g.L(BAS)    = 0.10; E.L(BAS)    = -70;

% -- Proximal apical trunk: HCN + L-type Ca --
g.h(PROX)   = 0.5;  g.CaL(PROX) = 0.5;
g.L(PROX)   = 0.15; E.L(PROX)   = -72;

% -- Distal apical tuft: high HCN (distal gradient) + HVA Ca + SK --
g.h(DIST)   = 1.0;  g.CaHVA(DIST) = getdef(opts,'gCaHVA_dist',2.0); g.SK(DIST) = 3.0;
g.L(DIST)   = 0.20; E.L(DIST)   = -72;

%% ---------------------- 5. SYNAPSE PARAMETERS ---------------------------
% AMPA/NMDA kinetics (ms^-1) from the synaptic script in the document.
syn.alpha_AMPA = 1.1;   syn.beta_AMPA = 0.5;
syn.alpha_NMDA = 0.072; syn.beta_NMDA = 0.015;
syn.Mg         = 1.0;   % extracellular [Mg2+] (mM) for NMDA block

% Synaptic conductance densities (mS/cm^2). Present only on dendrites.
g.AMPA = zeros(1,nComp);  g.NMDA = zeros(1,nComp);
g.AMPA(BAS)  = getdef(opts,'gAMPA_bas', 0.35); % MD-thalamic feedforward
g.NMDA(BAS)  = getdef(opts,'gNMDA_bas', 0.20);
g.AMPA(DIST) = getdef(opts,'gAMPA_dist',0.60); % ACC top-down feedback
g.NMDA(DIST) = getdef(opts,'gNMDA_dist',0.80);

% Presynaptic Poisson glutamate rates (Hz) driving each dendritic synapse.
rate_Glu    = zeros(1,nComp);
rate_Glu(BAS)  = getdef(opts,'rate_basal', 40);   % MD-thalamic feedforward
rate_Glu(DIST) = getdef(opts,'rate_distal',20);   % ACC top-down (tonic)
rel_steps   = round(1/dt);   % 1 ms glutamate pulse duration (samples)
glu_counter = zeros(1,nComp);

% Optional transient high-rate glutamate BURST onto the distal tuft, used to
% probe BAC firing / dendritic Ca spikes (set rate to 0 to disable).
burst_on   = getdef(opts,'burst_on',   300);   % ms
burst_off  = getdef(opts,'burst_off',  340);   % ms
burst_rate = getdef(opts,'burst_rate', 0);     % Hz onto distal tuft (0 = off)

%% ---------------------- 6. CALCIUM DYNAMICS -----------------------------
% Single-pool accumulation model (proximal <- I_CaL, distal <- I_CaHVA):
%   d[Ca]/dt = -phi_ca * I_Ca - ([Ca]-Ca_rest)/tau_ca
Ca_rest = 0.05;    % resting [Ca]i (uM)
phi_ca  = 0.005;   % influx scaling (uM per (uA/cm^2) per ms)
tau_ca  = 80;      % Ca extrusion time constant (ms)
% SK Ca-gating (Hill):
SK_EC50 = 0.5;     % half-activation [Ca]i (uM)
SK_nH   = 4.0;     % Hill coefficient
SK_tau  = 10.0;    % SK gating time constant (ms)

%% ---------------------- 7. NEUROMODULATION / PATHOLOGY HOOKS ------------
% Default = 1.0 (no modulation). These are the knobs for later VTA-DA and
% SHANK3-ASD experiments; they scale maximal conductances in place.
mod.D1_NMDA_gain    = getdef(opts,'D1_NMDA_gain', 1.00); % D1 (VTA) enhances NMDA  (e.g. 1.35)
mod.D2_AMPA_gain    = getdef(opts,'D2_AMPA_gain', 1.00); % D2 (VTA) suppresses AMPA (e.g. 0.85)
mod.D1_CaL_gain     = getdef(opts,'D1_CaL_gain',  1.00); % D1-PKA enhances L-type Ca
mod.SHANK3_AMPA     = getdef(opts,'SHANK3_AMPA',  1.00); % SHANK3 loss -> AMPA down-reg (<1)
mod.SHANK3_h        = getdef(opts,'SHANK3_h',     1.00); % SHANK3 loss -> dendritic I_h loss (<1)
mod.SHANK3_NMDA     = getdef(opts,'SHANK3_NMDA',  1.00); % PSD disruption -> NMDA change
% Apply hooks to the working conductance set.
g.AMPA = g.AMPA .* mod.D2_AMPA_gain .* mod.SHANK3_AMPA;
g.NMDA = g.NMDA .* mod.D1_NMDA_gain .* mod.SHANK3_NMDA;
g.CaL  = g.CaL  .* mod.D1_CaL_gain;
g.h    = g.h    .* mod.SHANK3_h;

%% ---------------------- 8. EXTERNAL CURRENT INJECTION -------------------
% Somatic depolarizing step to probe intrinsic firing + adaptation.
I_inj = zeros(nComp, N);
% (i) Somatic depolarizing step -> probes intrinsic firing + adaptation.
stim_on  = getdef(opts,'stim_on', 150);   % ms
stim_off = getdef(opts,'stim_off',350);   % ms
I_step   = getdef(opts,'I_step',  45.0);  % uA/cm^2 (somatic)
I_inj(SOMA, t >= stim_on & t < stim_off) = I_step;

% (ii) Distal apical current pulse -> probes the apical Ca/SK machinery
%      (dendritic Ca spike, like a dendritic patch-electrode injection).
dist_on  = getdef(opts,'dist_on', 430);   % ms
dist_off = getdef(opts,'dist_off',480);   % ms
I_dist   = getdef(opts,'I_dist',  40.0);  % uA/cm^2 (distal tuft)
I_inj(DIST, t >= dist_on & t < dist_off) = I_dist;

%% ---------------------- 9. STATE ALLOCATION ----------------------------
V  = E.L(:);                      % nComp x 1 initial voltages (leak rest)
Vrec = zeros(nComp, N);  Vrec(:,1) = V;

% Gating variables (scalars per compartment, steady-state initialized)
m_Na = zeros(1,nComp); h_Na = zeros(1,nComp); n_Kdr = zeros(1,nComp);
m_M  = zeros(1,nComp);
a_A  = zeros(1,nComp); b_A  = zeros(1,nComp);
m_h  = zeros(1,nComp);
s_CaL= zeros(1,nComp); f_CaL= zeros(1,nComp);
m_HVA= zeros(1,nComp); h_HVA= zeros(1,nComp);
z_SK = zeros(1,nComp);
r_AMPA = zeros(1,nComp); r_NMDA = zeros(1,nComp);

Ca = Ca_rest * ones(1,nComp);     % intracellular [Ca]i per compartment

for c = 1:nComp
    Vc = V(c);
    [am,bm,ah,bh,an,bn] = hh_rates(Vc);
    m_Na(c) = am/(am+bm);  h_Na(c) = ah/(ah+bh);  n_Kdr(c) = an/(an+bn);
    [aM,bM] = rates_M(Vc);              m_M(c)  = aM/(aM+bM);
    [aa,ba,ab,bb] = rates_A(Vc);        a_A(c)  = aa/(aa+ba); b_A(c) = ab/(ab+bb);
    m_h(c)  = hcn_inf(Vc);
    [as_,bs_,af_,bf_] = rates_CaL(Vc);  s_CaL(c)= as_/(as_+bs_); f_CaL(c)= af_/(af_+bf_);
    [amv,bmv,ahv,bhv] = rates_CaHVA(Vc);m_HVA(c)= amv/(amv+bmv); h_HVA(c)= ahv/(ahv+bhv);
    z_SK(c) = sk_inf(Ca(c), SK_EC50, SK_nH);
end

% Recording arrays for key traces
Ca_rec = zeros(nComp, N);  Ca_rec(:,1) = Ca';
I_M_rec    = zeros(1,N);   % somatic M-current
I_SK_rec   = zeros(1,N);   % distal SK current
I_CaHVA_rec= zeros(1,N);   % distal HVA Ca current
I_h_rec    = zeros(1,N);   % distal HCN current
I_syn_rec  = zeros(nComp,N); % total synaptic current per compartment

%% ---------------------- 10. MAIN INTEGRATION LOOP ----------------------
for k = 1:N-1
    Iion   = zeros(nComp,1);   % total intrinsic+synaptic current per comp

    % --- Presynaptic Poisson glutamate (Bernoulli approx. per step) ---
    for c = [BAS DIST]
        r_c = rate_Glu(c);
        if c == DIST && t(k) >= burst_on && t(k) < burst_off
            r_c = burst_rate;                 % transient distal tuft burst
        end
        if r_c > 0 && rand < r_c * dt/1000
            glu_counter(c) = rel_steps;
        end
    end
    T_Glu = double(glu_counter > 0);          % transmitter present (0/1)

    for c = 1:nComp
        Vc = V(c);
        Ic = 0;                                % accumulate g*(V-E) terms

        % ---- Fast Na + delayed-rectifier K (soma, basal) ----
        if g.Na(c) > 0
            [am,bm,ah,bh,an,bn] = hh_rates(Vc);
            m_Na(c)  = m_Na(c)  + dt*(am*(1-m_Na(c))  - bm*m_Na(c));
            h_Na(c)  = h_Na(c)  + dt*(ah*(1-h_Na(c))  - bh*h_Na(c));
            n_Kdr(c) = n_Kdr(c) + dt*(an*(1-n_Kdr(c)) - bn*n_Kdr(c));
            Ic = Ic + g.Na(c) *m_Na(c)^3*h_Na(c)*(Vc - E_Na);
            Ic = Ic + g.Kdr(c)*n_Kdr(c)^4      *(Vc - E_K);
        end

        % ---- M-current (soma) ----
        if g.M(c) > 0
            [aM,bM] = rates_M(Vc);
            m_M(c) = m_M(c) + dt*(aM*(1-m_M(c)) - bM*m_M(c));
            IM = g.M(c)*m_M(c)*(Vc - E_K);
            Ic = Ic + IM;
            if c == SOMA, I_M_rec(k) = IM; end
        end

        % ---- A-type K (basal) ----
        if g.A(c) > 0
            [aa,ba,ab,bb] = rates_A(Vc);
            a_A(c) = a_A(c) + dt*(aa*(1-a_A(c)) - ba*a_A(c));
            b_A(c) = b_A(c) + dt*(ab*(1-b_A(c)) - bb*b_A(c));
            Ic = Ic + g.A(c)*a_A(c)^3*b_A(c)*(Vc - E_K);
        end

        % ---- HCN / I_h (proximal, distal) ----
        if g.h(c) > 0
            minf = hcn_inf(Vc);  tauh = hcn_tau(Vc);
            m_h(c) = m_h(c) + dt*((minf - m_h(c))/tauh);
            Ih = g.h(c)*m_h(c)*(Vc - E_h);
            Ic = Ic + Ih;
            if c == DIST, I_h_rec(k) = Ih; end
        end

        % ---- L-type Ca (proximal) ----
        I_Ca_c = 0;
        if g.CaL(c) > 0
            [as_,bs_,af_,bf_] = rates_CaL(Vc);
            s_CaL(c) = s_CaL(c) + dt*(as_*(1-s_CaL(c)) - bs_*s_CaL(c));
            f_CaL(c) = f_CaL(c) + dt*(af_*(1-f_CaL(c)) - bf_*f_CaL(c));
            ICaL = g.CaL(c)*s_CaL(c)^2*f_CaL(c)*(Vc - E_Ca);
            Ic = Ic + ICaL;  I_Ca_c = I_Ca_c + ICaL;
        end

        % ---- HVA Ca (distal) ----
        if g.CaHVA(c) > 0
            [amv,bmv,ahv,bhv] = rates_CaHVA(Vc);
            m_HVA(c) = m_HVA(c) + dt*(amv*(1-m_HVA(c)) - bmv*m_HVA(c));
            h_HVA(c) = h_HVA(c) + dt*(ahv*(1-h_HVA(c)) - bhv*h_HVA(c));
            ICaHVA = g.CaHVA(c)*m_HVA(c)^2*h_HVA(c)*(Vc - E_Ca);
            Ic = Ic + ICaHVA;  I_Ca_c = I_Ca_c + ICaHVA;
            if c == DIST, I_CaHVA_rec(k) = ICaHVA; end
        end

        % ---- SK Ca-activated K (distal) ----
        if g.SK(c) > 0
            zinf = sk_inf(Ca(c), SK_EC50, SK_nH);
            z_SK(c) = z_SK(c) + dt*((zinf - z_SK(c))/SK_tau);
            ISK = g.SK(c)*z_SK(c)*(Vc - E_K);
            Ic = Ic + ISK;
            if c == DIST, I_SK_rec(k) = ISK; end
        end

        % ---- Leak ----
        Ic = Ic + g.L(c)*(Vc - E.L(c));

        % ---- AMPA / NMDA synapses (basal, distal) ----
        Isyn = 0;
        if g.AMPA(c) > 0 || g.NMDA(c) > 0
            r_AMPA(c) = r_AMPA(c) + dt*(syn.alpha_AMPA*T_Glu(c)*(1-r_AMPA(c)) - syn.beta_AMPA*r_AMPA(c));
            r_NMDA(c) = r_NMDA(c) + dt*(syn.alpha_NMDA*T_Glu(c)*(1-r_NMDA(c)) - syn.beta_NMDA*r_NMDA(c));
            Bv = 1/(1 + (syn.Mg/3.57)*exp(-0.062*Vc));   % Mg2+ block
            Iampa = g.AMPA(c)*r_AMPA(c)*(Vc - E_syn);
            Inmda = g.NMDA(c)*r_NMDA(c)*Bv*(Vc - E_syn);
            Isyn = Iampa + Inmda;
            Ic = Ic + Isyn;
        end
        I_syn_rec(c,k) = Isyn;

        % ---- Update [Ca]i for compartments with Ca channels ----
        if g.CaL(c) > 0 || g.CaHVA(c) > 0
            Ca(c) = Ca(c) + dt*(-phi_ca*I_Ca_c - (Ca(c)-Ca_rest)/tau_ca);
            Ca(c) = max(Ca(c), Ca_rest);
        end

        Iion(c) = Ic;
    end

    % --- Axial coupling + injection -> voltage update (all comps) ---
    for c = 1:nComp
        Icouple = 0;
        for j = 1:nComp
            if gc_den(c,j) > 0
                Icouple = Icouple + gc_den(c,j)*(V(j) - V(c));
            end
        end
        dV = (-Iion(c) + Icouple + I_inj(c,k)) / Cm;
        Vrec(c,k+1) = V(c) + dt*dV;
    end
    V = Vrec(:,k+1);

    % Decrement glutamate pulse counters
    glu_counter = max(0, glu_counter - 1);
    Ca_rec(:,k+1) = Ca';
end

% Backfill final samples of recorded currents
I_M_rec(end)=I_M_rec(end-1);   I_SK_rec(end)=I_SK_rec(end-1);
I_CaHVA_rec(end)=I_CaHVA_rec(end-1); I_h_rec(end)=I_h_rec(end-1);
I_syn_rec(:,end)=I_syn_rec(:,end-1);

%% ---------------------- 11. PLOTTING -----------------------------------
col = lines(nComp);

fig = figure('Name','L5 PFC 4-Compartment Pyramidal Neuron', ...
             'Position',[40 40 900 950]);

% (a) Membrane voltages of all four compartments
subplot(4,1,1); hold on;
for c = 1:nComp
    plot(t, Vrec(c,:), 'Color', col(c,:), 'LineWidth', 1.2, ...
         'DisplayName', cName{c});
end
title(sprintf(['Membrane Potential by Compartment  ' ...
    '(somatic step %g-%g ms; distal Ca-spike pulse %g-%g ms)'], ...
    stim_on, stim_off, dist_on, dist_off));
ylabel('V_m (mV)'); legend('Location','northeast'); grid on; xlim([0 t_max]);

% (b) Intracellular calcium (proximal L-type pool, distal HVA pool)
subplot(4,1,2); hold on;
plot(t, Ca_rec(PROX,:), 'Color', col(PROX,:), 'LineWidth', 1.3, ...
     'DisplayName','Proximal [Ca] (I_{CaL})');
plot(t, Ca_rec(DIST,:), 'Color', col(DIST,:), 'LineWidth', 1.3, ...
     'DisplayName','Distal [Ca] (I_{CaHVA})');
title('Dendritic Intracellular Calcium');
ylabel('[Ca^{2+}]_i (\muM)'); legend('Location','northeast'); grid on; xlim([0 t_max]);

% (c) Key intrinsic currents: somatic I_M, distal I_SK, I_CaHVA, I_h
subplot(4,1,3); hold on;
plot(t, I_M_rec,    'LineWidth', 1.1, 'DisplayName','Soma I_M');
plot(t, I_SK_rec,   'LineWidth', 1.1, 'DisplayName','Distal I_{SK}');
plot(t, I_CaHVA_rec,'LineWidth', 1.1, 'DisplayName','Distal I_{CaHVA}');
plot(t, I_h_rec,    'LineWidth', 1.1, 'DisplayName','Distal I_h');
title('Selected Intrinsic Currents');
ylabel('I (\muA/cm^2)'); legend('Location','northeast'); grid on; xlim([0 t_max]);

% (d) Dendritic synaptic currents (AMPA+NMDA)
subplot(4,1,4); hold on;
plot(t, I_syn_rec(BAS,:),  'Color', col(BAS,:),  'LineWidth', 1.0, ...
     'DisplayName','Basal AMPA+NMDA (MD thalamus)');
plot(t, I_syn_rec(DIST,:), 'Color', col(DIST,:), 'LineWidth', 1.0, ...
     'DisplayName','Distal AMPA+NMDA (ACC)');
title('Dendritic Synaptic Currents (Poisson glutamate drive)');
xlabel('Time (ms)'); ylabel('I_{syn} (\muA/cm^2)');
legend('Location','northeast'); grid on; xlim([0 t_max]);

% Optional output (for analysis / testing); ignored in normal script use.
if nargout > 0
    out = struct();
    out.t = t;  out.V = Vrec;  out.Ca = Ca_rec;
    out.I_M = I_M_rec;  out.I_SK = I_SK_rec;
    out.I_CaHVA = I_CaHVA_rec;  out.I_h = I_h_rec;
    out.I_syn = I_syn_rec;  out.gc_den = gc_den;  out.area_cm2 = area_cm2;
    out.idx = struct('SOMA',SOMA,'BAS',BAS,'PROX',PROX,'DIST',DIST);
end

end  % ===================== end main function ============================


%% ======================= LOCAL RATE FUNCTIONS ==========================

% ---- Fast Na & delayed-rectifier K (classic Hodgkin-Huxley, rest ~ -65) ----
function [am,bm,ah,bh,an,bn] = hh_rates(V)
    xm = V + 40;
    if abs(xm) < 1e-6, am = 1.0; else, am = 0.1*xm/(1-exp(-xm/10)); end
    bm = 4*exp(-(V+65)/18);
    ah = 0.07*exp(-(V+65)/20);
    bh = 1/(1+exp(-(V+35)/10));
    xn = V + 55;
    if abs(xn) < 1e-6, an = 0.1; else, an = 0.01*xn/(1-exp(-xn/10)); end
    bn = 0.125*exp(-(V+65)/80);
end

% ---- M-current I_M (Kv7): slow, non-inactivating, V_half = -35 mV ----
function [a,b] = rates_M(V)
    k = 0.003;
    a = k / (1 + exp(-(V+35)/10));
    b = (k * exp(-(V+35)/10)) / (1 + exp(-(V+35)/10));
end

% ---- A-type K I_A: fast activation (a), fast/moderate inactivation (b) ----
function [aa,ba,ab,bb] = rates_A(V)
    k_a = 0.40; k_b = 0.05;
    aa = k_a / (1 + exp(-(V+30)/10));
    ba = (k_a * exp(-(V+30)/10)) / (1 + exp(-(V+30)/10));
    ab = k_b / (1 + exp((V+60)/7));
    bb = (k_b * exp((V+60)/7)) / (1 + exp((V+60)/7));
end

% ---- HCN I_h activation (Kole-type): activates on hyperpolarization ----
function minf = hcn_inf(V)
    minf = 1 / (1 + exp((V+82)/9));
end
function tau = hcn_tau(V)
    tau = 1 / (exp(-14.59 - 0.086*V) + exp(-1.87 + 0.0701*V));
end

% ---- L-type Ca I_CaL: fast activation (s^2), slow VDI inactivation (f) ----
function [as_,bs_,af_,bf_] = rates_CaL(V)
    k_s = 0.40; k_f = 0.01;
    as_ = k_s / (1 + exp(-(V+15)/8));
    bs_ = (k_s * exp(-(V+15)/8)) / (1 + exp(-(V+15)/8));
    af_ = k_f / (1 + exp((V+35)/8));
    bf_ = (k_f * exp((V+35)/8)) / (1 + exp((V+35)/8));
end

% ---- Non-L HVA Ca I_CaHVA: fast activation (m^2), moderate inactivation (h) ----
function [am,bm,ah,bh] = rates_CaHVA(V)
    k_m = 0.50; k_h = 0.02;
    am = k_m / (1 + exp(-(V+10)/8));
    bm = (k_m * exp(-(V+10)/8)) / (1 + exp(-(V+10)/8));
    ah = k_h / (1 + exp((V+40)/10));
    bh = (k_h * exp((V+40)/10)) / (1 + exp((V+40)/10));
end

% ---- SK Ca-activated K: voltage-independent Hill activation by [Ca]i ----
function zinf = sk_inf(Ca, EC50, nH)
    zinf = (Ca^nH) / (Ca^nH + EC50^nH);
end

% ---- Option helper: return opts.field if present, else the default ----
function val = getdef(opts, field, default)
    if isfield(opts, field) && ~isempty(opts.(field))
        val = opts.(field);
    else
        val = default;
    end
end
