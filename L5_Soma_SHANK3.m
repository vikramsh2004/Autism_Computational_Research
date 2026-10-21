function out = L5_Soma_SHANK3(opts)
%% L5 PFC PYRAMIDAL NEURON - ISOLATED SOMATIC COMPARTMENT (SHANK3-KO ready)
% =========================================================================
% WHAT THIS SIMULATES
% -------------------------------------------------------------------------
% A single, isolated somatic compartment of a Layer-5 prefrontal-cortex
% (PFC) pyramidal neuron, built as the first "unit model" of the larger
% 4-compartment cable (soma / basal / proximal / distal). Per the model
% design document, the soma is the spike generator and the seat of
% spike-frequency adaptation (SFA). It contains exactly four conductances
% plus passive properties:
%
%   I_Na   fast Na+          (Nav1.6)       g_Na  = 120.0 mS/cm^2  E_Na = +50
%   I_Kdr  delayed rectifier (Kv3.1/3.2)    g_Kdr =  36.0 mS/cm^2  E_K  = -90
%   I_M    M-type K+         (Kv7/KCNQ)     g_M   =   1.5 mS/cm^2  E_K  = -90
%   I_L    passive leak                     g_L   =   0.10 mS/cm^2 E_L  = -70
%   Cm = 1.0 uF/cm^2,  A = 1570.8 um^2
%
% SHANK3-KO PATHOLOGY (soma)
% -------------------------------------------------------------------------
% SHANK3 is a post-synaptic-density scaffold confined to dendritic spines,
% so there are NO direct structural changes to somatic channels. The somatic
% phenotype is a MALADAPTIVE HOMEOSTATIC compensation for the chronic loss of
% dendritic excitatory drive. Two biologically grounded intrinsic changes
% (the only two supported by the SHANK3-ASD literature for an isolated soma):
%
%   (1) Somatic M-current down-regulation
%         g_M : 1.5 -> 1.0 mS/cm^2   (~33% reduction)
%       -> removes the slow K+ "brake" -> REDUCED spike-frequency adaptation,
%          steeper somatic F-I curve.
%
%   (2) Hyperpolarizing shift of fast-Na+ activation (AIS remodeling)
%         V_1/2,m : -40 -> -43 mV    (-3 mV shift)
%       -> LOWERS spike-initiation threshold (~-50 -> ~-53 mV)
%          -> lower rheobase, soma hyperexcitable to weak depolarizations.
%
% Everything else at the soma (g_Kdr spike repolarization kinetics, passive
% leak g_L / E_L, resting potential) is deliberately held FIXED so the unit
% test isolates the two intrinsic SHANK3 changes without confounds.
%
% USAGE
% -------------------------------------------------------------------------
%   out = L5_Soma_SHANK3();                       % Wild-type default
%   out = L5_Soma_SHANK3(struct('state','SHANK3'));% SHANK3-KO preset
%   out = L5_Soma_SHANK3(struct('I_step',20));     % custom step current
%
% Any field may be overridden via opts; 'state' only sets the two pathology
% knobs (gM, na_vshift) and may be further overridden by explicit fields.
%
% Units: V [mV], t [ms], g [mS/cm^2], I [uA/cm^2], C [uF/cm^2].
% =========================================================================

if nargin < 1 || isempty(opts), opts = struct(); end

%% 0. PATHOLOGY STATE PRESET ----------------------------------------------
state = getdef(opts, 'state', 'WT');
switch upper(state)
    case {'WT','CONTROL','WILDTYPE','WILD-TYPE'}
        gM_def       = 1.5;   % mS/cm^2  (full M-current)
        na_vshift_def = 0.0;  % mV       (standard Na activation)
    case {'SHANK3','SHANK3-KO','KO','ASD'}
        gM_def       = 1.0;   % mS/cm^2  (~33% down-regulation)
        na_vshift_def = -3.0; % mV       (-3 mV hyperpolarizing Na shift)
    otherwise
        error('L5_Soma_SHANK3:badState', 'Unknown state "%s".', state);
end

%% 1. SIMULATION PARAMETERS -----------------------------------------------
dt    = getdef(opts, 'dt',    0.01);   % time step (ms)
t_max = getdef(opts, 't_max', 500);    % total time (ms)
t     = 0:dt:t_max;
N     = numel(t);
Cm    = getdef(opts, 'Cm', 1.0);       % membrane capacitance (uF/cm^2)

%% 2. REVERSAL POTENTIALS --------------------------------------------------
E_Na = 50;    % mV
E_K  = -90;   % mV
E_L  = getdef(opts, 'EL', -70);   % mV (leak reversal sets rest)

%% 3. MAXIMAL CONDUCTANCE DENSITIES (mS/cm^2) -----------------------------
gNa  = getdef(opts, 'gNa',  120.0);
gKdr = getdef(opts, 'gKdr',  36.0);
gM   = getdef(opts, 'gM',    gM_def);          % SHANK3 target #1
gL   = getdef(opts, 'gL',     0.10);

% M-current kinetic scale: sets the activation time constant (tau ~= 1/mM_k
% ms) WITHOUT changing the steady-state activation curve m_inf (which only
% depends on V_1/2 = -35 mV). The document used mM_k = 0.003 (tau ~ 330 ms);
% physiological Kv7/M-current activation is faster, so the default is 0.02
% (tau ~ 50 ms) which lets spike-frequency adaptation express within a
% sustained step. Set opts.mM_k = 0.003 to reproduce the document exactly.
mM_k = getdef(opts, 'mM_k', 0.02);

% SHANK3 target #2: hyperpolarizing shift of Na+ activation (m-gate only).
% na_vshift < 0 shifts V_1/2,m in the hyperpolarizing direction.
na_vshift = getdef(opts, 'na_vshift', na_vshift_def);

%% 4. CURRENT INJECTION ----------------------------------------------------
stim_on  = getdef(opts, 'stim_on',  150);   % ms
stim_off = getdef(opts, 'stim_off', 350);   % ms (200 ms step per the doc)
I_step   = getdef(opts, 'I_step',    45.0); % uA/cm^2 somatic step

I_inj = zeros(1, N);
I_inj(t >= stim_on & t < stim_off) = I_step;

%% 5. STATE ALLOCATION & STEADY-STATE INITIALIZATION ----------------------
V = E_L;                       % start at leak rest
Vrec = zeros(1, N); Vrec(1) = V;

[am,bm,ah,bh,an,bn] = hh_rates(V, na_vshift);
m_Na = am/(am+bm);
h_Na = ah/(ah+bh);
n_Kdr = an/(an+bn);
[aM,bM] = rates_M(V, mM_k);
m_M = aM/(aM+bM);

%% 6. MAIN INTEGRATION LOOP (forward Euler) -------------------------------
for k = 1:N-1
    Vc = V;

    % Fast Na+ (m^3 h) with SHANK3 activation shift on the m-gate
    [am,bm,ah,bh,an,bn] = hh_rates(Vc, na_vshift);
    m_Na  = m_Na  + dt*(am*(1-m_Na)  - bm*m_Na);
    h_Na  = h_Na  + dt*(ah*(1-h_Na)  - bh*h_Na);
    n_Kdr = n_Kdr + dt*(an*(1-n_Kdr) - bn*n_Kdr);
    I_Na  = gNa  * m_Na^3 * h_Na * (Vc - E_Na);
    I_Kdr = gKdr * n_Kdr^4        * (Vc - E_K);

    % M-current (slow, non-inactivating K+) -> spike-frequency adaptation
    [aM,bM] = rates_M(Vc, mM_k);
    m_M = m_M + dt*(aM*(1-m_M) - bM*m_M);
    I_M = gM * m_M * (Vc - E_K);

    % Passive leak
    I_L = gL * (Vc - E_L);

    % Membrane update
    Iion = I_Na + I_Kdr + I_M + I_L;
    dV   = (-Iion + I_inj(k)) / Cm;
    V    = Vc + dt*dV;
    Vrec(k+1) = V;
end

%% 7. SPIKE DETECTION & ADAPTATION METRICS --------------------------------
Vth = getdef(opts, 'Vth', 0.0);  % mV upward-crossing threshold
above = Vrec > Vth;
cross = find(~above(1:end-1) & above(2:end));   % upward crossings
spike_t = t(cross + 1);

% Spikes that occur during the stimulus window
in_stim = spike_t >= stim_on & spike_t <= stim_off + 5;
spike_stim = spike_t(in_stim);
nSpikes = numel(spike_stim);

ISI = diff(spike_stim);                 % ms
if numel(ISI) >= 2
    ISI_ratio = ISI(end) / ISI(1);      % last/first  (>1 => adaptation)
else
    ISI_ratio = NaN;
end

% Resting potential: mean over the 50 ms before stimulus onset
pre = t < stim_on & t > (stim_on - 50);
Vrest = mean(Vrec(pre));
Vpeak = max(Vrec);

%% 8. OUTPUT ---------------------------------------------------------------
out = struct();
out.state     = upper(state);
out.t         = t;
out.V         = Vrec;
out.I_inj     = I_inj;
out.spike_t   = spike_t;
out.spike_stim = spike_stim;
out.nSpikes   = nSpikes;
out.ISI       = ISI;
out.ISI_ratio = ISI_ratio;
out.Vrest     = Vrest;
out.Vpeak     = Vpeak;
out.params    = struct('gNa',gNa,'gKdr',gKdr,'gM',gM,'gL',gL, ...
                       'na_vshift',na_vshift,'I_step',I_step, ...
                       'stim_on',stim_on,'stim_off',stim_off, ...
                       'E_Na',E_Na,'E_K',E_K,'E_L',E_L,'Cm',Cm);
end % ===================== end main function =============================


%% ------------------------- LOCAL RATE FUNCTIONS -------------------------

function [am,bm,ah,bh,an,bn] = hh_rates(V, na_vshift)
% Classic Hodgkin-Huxley fast-Na / delayed-rectifier-K rates (rest ~ -65).
% na_vshift (mV) shifts ONLY the Na+ activation (m) gate: the rates are
% evaluated at (V - na_vshift), so na_vshift < 0 moves V_1/2,m in the
% hyperpolarizing direction (SHANK3 AIS remodeling: -40 -> -43 mV).
if nargin < 2, na_vshift = 0; end
Vm = V - na_vshift;                 % m-gate "sees" a shifted voltage

xm = Vm + 40;
if abs(xm) < 1e-6, am = 1.0; else, am = 0.1*xm/(1-exp(-xm/10)); end
bm = 4*exp(-(Vm+65)/18);

% Inactivation (h) and K+ activation (n) use the true membrane voltage V.
ah = 0.07*exp(-(V+65)/20);
bh = 1/(1+exp(-(V+35)/10));

xn = V + 55;
if abs(xn) < 1e-6, an = 0.1; else, an = 0.01*xn/(1-exp(-xn/10)); end
bn = 0.125*exp(-(V+65)/80);
end

function [a,b] = rates_M(V, k)
% M-current I_M (Kv7/KCNQ): slow, non-inactivating, V_half = -35 mV.
% k scales the forward/backward rates and hence the activation time
% constant tau = 1/(a+b) = 1/k; the steady-state activation
% m_inf = a/(a+b) = 1/(1+exp(-(V+35)/10)) is independent of k.
if nargin < 2 || isempty(k), k = 0.02; end
a = k / (1 + exp(-(V+35)/10));
b = (k * exp(-(V+35)/10)) / (1 + exp(-(V+35)/10));
end

function val = getdef(opts, field, default)
% Return opts.(field) if present and non-empty, else the default.
if isfield(opts, field) && ~isempty(opts.(field))
    val = opts.(field);
else
    val = default;
end
end
