%% TRN T-type Ca2+ Channel Over-expression  (ASD / PTEN-like model)
% ========================================================================
% WHAT THIS SCRIPT SIMULATES
% ------------------------------------------------------------------------
% The thalamic reticular nucleus (TRN) is a shell of GABAergic neurons that
% inhibit thalamocortical (TC) relay cells and thereby GATE / FILTER the
% flow of sensory information to cortex. TRN neurons rely on low-threshold
% T-type Ca2+ channels (I_T) to fire rhythmic bursts; those bursts drive a
% strong, SUSTAINED ("long-lasting") GABA inhibition onto TC cells.
%
% Several genetic ASD models (e.g. PTEN mutations) show TRN neurons with
% OVER-EXPRESSED / MUTANT T-type channels. The channels are hyperactive
% (a huge, poorly-inactivating inward Ca2+ current) yet the neuron
% PARADOXICALLY fails to deliver long-lasting inhibition. This model shows
% why: the excess, non-inactivating T-type current depolarizes the TRN into
% DEPOLARIZATION BLOCK. The membrane is hyper-excited, but it can no longer
% fire the crisp rhythmic spikes needed to inhibit TC. With inhibition gone,
% the TC relay cell escapes and fires a full spike train to the same
% sensory drive -> a minimal in-silico model of SENSORY OVERLOAD.
%
% THE MUTANT IS MODELLED WITH TWO T-TYPE CHANGES
%   (1) OVER-EXPRESSION      -> larger maximal conductance   (gT_TRN up)
%   (2) LOSS OF INACTIVATION -> a floor on the inactivation gate h, so the
%                               channel cannot fully close    (hfloor up)
%
% MODEL (2-cell reciprocal loop, 1 neuron per region)
%   TC  cell : Hodgkin-Huxley Na/K/leak + T-type Ca2+ (rebound-capable)
%   TRN cell : Hodgkin-Huxley Na/K/leak + T-type Ca2+  (the ASD locus)
%   Synapses : TC --AMPA--> TRN (excitation)
%              TRN --GABA_A--> TC (inhibition; slow decay = long-lasting)
%   Drive    : a sustained "sensory" current to TC and a background
%              corticothalamic current to TRN during the stimulus.
%   Release  : event-based (per spike), so a cell in depolarization block
%              (no threshold crossings) releases almost no transmitter.
%
% HOW TO READ THE OUTPUT
%   Fig 1  Membrane potentials (Control vs Over-expression) for TRN and TC.
%          -> Control TRN bursts rhythmically and keeps TC quiet (gated).
%          -> Over-expression TRN locks into a depolarized plateau (block);
%             TC then fires a full train  = SENSORY OVERLOAD.
%   Fig 2  Mechanism: the (hyperactive) TRN T-type current, the TRN
%          inactivation gate, and the GABA conductance delivered to TC
%          (sustained for control, collapses for the mutant).
%   Fig 3  Dose-response: sweeping T-type over-expression shows the paradox
%          as an INVERTED-U -- TRN output first rises (hyperactive) then
%          COLLAPSES into block, while TC output jumps up (overload).
%
% Tunable knobs to explore: cond().gT_TRN, cond().hfloor, P.gGABA, P.bG
% (GABA decay = how "long-lasting" inhibition is), P.Iext_TC, P.Iext_TRN.
% ========================================================================

clear; close all; clc;

%% ------------------------- Time base -----------------------------------
P.dt    = 0.01;                 % ms
P.t_end = 400;                  % ms
P.time  = 0:P.dt:P.t_end;
P.nt    = numel(P.time);

%% ------------------------- Fixed biophysics ----------------------------
P.Cm  = 1.0;                    % uF/cm^2
P.gL  = 0.15;                   % leak conductance (high input resistance)
P.EL_TC  = -70;                 % TC leak reversal
P.EL_TRN = -80;                 % TRN rests hyperpolarized -> I_T de-inactivated
P.gT_TC  = 1.5;                 % TC T-type (lets TC do post-inhibitory rebound)
P.E_Ca   = 120;                 % T-type Ca2+ reversal

% Synapses
P.gAMPA = 0.30;                 % TC -> TRN excitation
P.gGABA = 0.60;                 % TRN -> TC inhibition
P.EGABA = -85;                  % GABA_A reversal
P.bG    = 0.05;                 % GABA decay rate (small = SLOW = long-lasting)

% Drive (a maintained sensory volley)
on  = round(40  / P.dt) + 1;
off = round(380 / P.dt) + 1;
P.Iext_TC  = zeros(1, P.nt);  P.Iext_TC(on:off)  = 7.0;  % sensory drive to TC
P.Iext_TRN = zeros(1, P.nt);  P.Iext_TRN(on:off) = 3.0;  % corticothalamic drive to TRN

%% ------------------------- Conditions ----------------------------------
% Control  : physiological T-type conductance, normal inactivation.
% ASD/PTEN : T-type OVER-EXPRESSION (gT up) + LOSS OF INACTIVATION (hfloor up).
cond(1) = struct('name','Control TRN',                'gT_TRN',3.0, 'hfloor',0.0, 'color',[0.10 0.45 0.80]);
cond(2) = struct('name','T-type Over-expression (ASD)','gT_TRN',8.0, 'hfloor',0.5, 'color',[0.85 0.20 0.20]);

%% ------------------------- Run the two conditions ----------------------
R = cell(1, numel(cond));
fprintf('\n%-32s | TC spikes | TRN spikes | GABA charge (a.u.)\n', 'Condition');
fprintf('%s\n', repmat('-',1,78));
for c = 1:numel(cond)
    R{c} = run_condition(P, cond(c));
    fprintf('%-32s | %8d  | %9d  | %10.1f\n', ...
        cond(c).name, R{c}.TC_spikes, R{c}.TRN_spikes, R{c}.GABA_charge);
end
fprintf('%s\n', repmat('-',1,78));
fprintf(['Paradox: over-expression makes the TRN HYPERACTIVE (huge T-type current)\n' ...
         'but it enters depolarization block, so it delivers LESS inhibition and the\n' ...
         'TC relay cell fires MORE  ->  sensory overload.\n\n']);

%% ------------------------- Dose-response sweep -------------------------
gT_scan   = 2:0.5:12;
TC_scan   = zeros(size(gT_scan));
rel_scan  = zeros(size(gT_scan));
for i = 1:numel(gT_scan)
    r = run_condition(P, struct('gT_TRN',gT_scan(i), 'hfloor',0.5, 'name','scan','color',[0 0 0]));
    TC_scan(i)  = r.TC_spikes;      % TC output (sensory overload index)
    rel_scan(i) = r.TRN_releases;   % effective TRN inhibitory spikes
end

%% =======================================================================
%% FIGURE 1 : Membrane potentials, Control vs Over-expression
%% =======================================================================
figure('Name','Membrane Potentials','Color','w','Position',[60 60 1050 640]);

subplot(2,2,1)
plot(P.time, R{1}.v_TRN, 'LineWidth',1.1, 'Color',cond(1).color); grid on
title('Control : TRN  (rhythmic bursting)'); ylabel('V_m (mV)'); ylim([-100 60]);

subplot(2,2,2)
plot(P.time, R{2}.v_TRN, 'LineWidth',1.1, 'Color',cond(2).color); grid on
title('Over-expression : TRN  (depolarization block)'); ylabel('V_m (mV)'); ylim([-100 60]);

subplot(2,2,3)
plot(P.time, R{1}.v_TC, 'LineWidth',1.1, 'Color',[0.15 0.15 0.15]); grid on
title(sprintf('Control : TC relay  (%d spikes = gated)', R{1}.TC_spikes));
xlabel('Time (ms)'); ylabel('V_m (mV)'); ylim([-100 60]);

subplot(2,2,4)
plot(P.time, R{2}.v_TC, 'LineWidth',1.1, 'Color',[0.15 0.15 0.15]); grid on
title(sprintf('Over-expression : TC relay  (%d spikes = OVERLOAD)', R{2}.TC_spikes));
xlabel('Time (ms)'); ylabel('V_m (mV)'); ylim([-100 60]);

sgtitle('Fig 1  |  Hyperactive T-type -> TRN depolarization block -> TC sensory overload');

%% =======================================================================
%% FIGURE 2 : Mechanism
%% =======================================================================
figure('Name','Mechanism','Color','w','Position',[100 100 1050 720]);

subplot(3,1,1)
plot(P.time, R{1}.I_T_TRN, 'LineWidth',1.3, 'Color',cond(1).color, 'DisplayName',cond(1).name); hold on
plot(P.time, R{2}.I_T_TRN, 'LineWidth',1.3, 'Color',cond(2).color, 'DisplayName',cond(2).name);
title('TRN T-type Ca^{2+} current  (mutant = huge, sustained = "hyperactive")');
ylabel('I_T (\muA/cm^2)'); legend('Location','best'); grid on

subplot(3,1,2)
plot(P.time, R{1}.hT_TRN, 'LineWidth',1.3, 'Color',cond(1).color, 'DisplayName',cond(1).name); hold on
plot(P.time, R{2}.hT_TRN, 'LineWidth',1.3, 'Color',cond(2).color, 'DisplayName',cond(2).name);
title('TRN T-type inactivation gate  (mutant cannot inactivate -> stays open)');
ylabel('h_{eff}'); ylim([0 1]); legend('Location','best'); grid on

subplot(3,1,3)
plot(P.time, R{1}.g_GABA, 'LineWidth',1.4, 'Color',cond(1).color, 'DisplayName',cond(1).name); hold on
plot(P.time, R{2}.g_GABA, 'LineWidth',1.4, 'Color',cond(2).color, 'DisplayName',cond(2).name);
title('GABA_A conductance delivered to TC  (control = sustained; mutant = collapses)');
xlabel('Time (ms)'); ylabel('g_{GABA} (mS/cm^2)'); legend('Location','best'); grid on

sgtitle('Fig 2  |  Paradox: an over-active T-type current that FAILS to sustain inhibition');

%% =======================================================================
%% FIGURE 3 : Dose-response (the inverted-U paradox) + summary
%% =======================================================================
figure('Name','Dose-response','Color','w','Position',[140 140 1000 420]);

subplot(1,2,1)
yyaxis left
plot(gT_scan, rel_scan, 'o-','LineWidth',1.5); ylabel('TRN effective spikes (inhibition)');
yyaxis right
plot(gT_scan, TC_scan, 's-','LineWidth',1.5); ylabel('TC spikes (sensory overload)');
xlabel('TRN T-type conductance  gT_{TRN}  (over-expression \rightarrow)');
title('Dose-response: TRN output collapses, TC output explodes');
grid on

subplot(1,2,2)
b = bar([R{1}.TC_spikes, R{2}.TC_spikes]); b.FaceColor='flat';
b.CData(1,:)=cond(1).color; b.CData(2,:)=cond(2).color;
set(gca,'XTickLabel',{'Control','Over-expr.'}); ylabel('TC spike count');
title('Sensory throughput to cortex'); grid on

sgtitle('Fig 3  |  More T-type channels -> less inhibition -> more TC firing (paradox)');

%% ======================================================================
%% ========================= LOCAL FUNCTIONS ============================
%% ======================================================================
function out = run_condition(P, C)
% Integrate the TC<->TRN loop for one condition (forward Euler).
    nt = P.nt; dt = P.dt;

    vTC  = -70;  vTRN = -80;
    [mTC,hTC,nTC]    = init_gates(vTC);
    [mTRN,hTRN,nTRN] = init_gates(vTRN);
    hT_TC  = hinf_T(vTC);        % T-type inactivation gates
    hT_TRN = hinf_T(vTRN);
    sA = 0; sG = 0;              % AMPA / GABA_A synaptic gates

    v_TC=zeros(1,nt); v_TRN=zeros(1,nt);
    I_T_TRN=zeros(1,nt); hT_TRN_rec=zeros(1,nt); g_GABA=zeros(1,nt);
    prev_vTC = vTC; prev_vTRN = vTRN; releases = 0;

    for k = 1:nt
        v_TC(k)=vTC; v_TRN(k)=vTRN; g_GABA(k)=P.gGABA*sG;

        % ---- event-based synaptic release (upward crossing of 0 mV) ----
        if prev_vTC  <= 0 && vTC  > 0, sA = min(sA+0.6, 1); end
        if prev_vTRN <= 0 && vTRN > 0, sG = min(sG+0.5, 1); releases = releases + 1; end
        prev_vTC = vTC; prev_vTRN = vTRN;

        % ---- synaptic gates decay (GABA is slow -> long-lasting) ----
        sA = sA + dt*(-1.0*sA);
        sG = sG + dt*(-P.bG*sG);
        I_AMPA = P.gAMPA * sA * (0      - vTRN);   % onto TRN
        I_GABA = P.gGABA * sG * (P.EGABA - vTC);   % onto TC

        % ---- Hodgkin-Huxley gate updates ----
        [am,bm,ah,bh,an,bn] = hh_rates(vTC);
        mTC=mTC+dt*(am*(1-mTC)-bm*mTC); hTC=hTC+dt*(ah*(1-hTC)-bh*hTC); nTC=nTC+dt*(an*(1-nTC)-bn*nTC);
        [am,bm,ah,bh,an,bn] = hh_rates(vTRN);
        mTRN=mTRN+dt*(am*(1-mTRN)-bm*mTRN); hTRN=hTRN+dt*(ah*(1-hTRN)-bh*hTRN); nTRN=nTRN+dt*(an*(1-nTRN)-bn*nTRN);

        % ---- T-type inactivation gates ----
        hT_TC  = hT_TC  + dt*(hinf_T(vTC)  - hT_TC ) / tauh_T(vTC);
        hT_TRN = hT_TRN + dt*(hinf_T(vTRN) - hT_TRN) / tauh_T(vTRN);
        % Mutant: a floor on inactivation -> channel cannot fully close.
        hEff_TRN = C.hfloor + (1 - C.hfloor)*hT_TRN;
        hT_TRN_rec(k) = hEff_TRN;

        % ---- intrinsic currents: TC ----
        INa = 120*mTC^3*hTC*(vTC-50);
        IK  = 36 *nTC^4    *(vTC+77);
        IL  = P.gL*(vTC - P.EL_TC);
        IT  = P.gT_TC*minf_T(vTC)^2*hT_TC*(vTC - P.E_Ca);
        dvTC = (P.Iext_TC(k) + I_GABA - INa - IK - IL - IT) / P.Cm;

        % ---- intrinsic currents: TRN (T-type = the ASD locus) ----
        INa2 = 120*mTRN^3*hTRN*(vTRN-50);
        IK2  = 36 *nTRN^4     *(vTRN+77);
        IL2  = P.gL*(vTRN - P.EL_TRN);
        IT2  = C.gT_TRN*minf_T(vTRN)^2*hEff_TRN*(vTRN - P.E_Ca);
        I_T_TRN(k) = IT2;
        dvTRN = (P.Iext_TRN(k) + I_AMPA - INa2 - IK2 - IL2 - IT2) / P.Cm;

        vTC  = vTC  + dt*dvTC;
        vTRN = vTRN + dt*dvTRN;
    end

    out.v_TC=v_TC; out.v_TRN=v_TRN;
    out.I_T_TRN=I_T_TRN; out.hT_TRN=hT_TRN_rec; out.g_GABA=g_GABA;
    out.TC_spikes  = count_spikes(v_TC);
    out.TRN_spikes = count_spikes(v_TRN);
    out.TRN_releases = releases;
    out.GABA_charge = sum(g_GABA)*dt;   % time-integrated inhibition delivered
end

function n = count_spikes(v)
    n = sum(v(1:end-1) <= 0 & v(2:end) > 0);
end

%% ---- T-type Ca2+ channel gating (Destexhe/Huguenard-style) ----
function m = minf_T(V), m = 1 ./ (1 + exp(-(V + 52) / 7.4)); end   % instantaneous activation
function h = hinf_T(V), h = 1 ./ (1 + exp( (V + 80) / 5 )); end    % steady-state inactivation
function tau = tauh_T(V)                                            % inactivation time constant (ms)
    tau = 30 + 220 ./ (1 + exp((V + 65) / 8));
    tau = max(tau, 1.0);
end

%% ---- Hodgkin-Huxley Na/K gating ----
function [m0,h0,n0] = init_gates(V)
    [am,bm,ah,bh,an,bn] = hh_rates(V);
    m0 = am/(am+bm);  h0 = ah/(ah+bh);  n0 = an/(an+bn);
end
function [am,bm,ah,bh,an,bn] = hh_rates(V)
    am = alpha_m(V); bm = beta_m(V);
    ah = alpha_h(V); bh = beta_h(V);
    an = alpha_n(V); bn = beta_n(V);
end
function val = alpha_m(V)
    x = V + 40;
    if abs(x) < 1e-6, val = 1.0; else, val = 0.1*x/(1-exp(-x/10)); end
end
function val = beta_m(V),  val = 4*exp(-(V+65)/18); end
function val = alpha_h(V), val = 0.07*exp(-(V+65)/20); end
function val = beta_h(V),  val = 1/(1+exp(-(V+35)/10)); end
function val = alpha_n(V)
    x = V + 55;
    if abs(x) < 1e-6, val = 0.1; else, val = 0.01*x/(1-exp(-x/10)); end
end
function val = beta_n(V),  val = 0.125*exp(-(V+65)/80); end
