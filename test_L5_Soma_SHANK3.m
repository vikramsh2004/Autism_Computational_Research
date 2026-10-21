%% DIAGNOSTIC TEST for the isolated somatic compartment (WT vs SHANK3-KO)
% =========================================================================
% Verifies that L5_Soma_SHANK3.m is "working correctly": that it produces a
% physiologically valid spike generator AND that the two SHANK3-KO knobs
% reproduce the three somatic pathology signatures described in the model
% design document:
%
%   (1) LOWER RHEOBASE   (hyperexcitability)      <- Na-activation -3 mV shift
%   (2) HIGHER SPIKE COUNT at a fixed step (steeper F-I)  <- both knobs
%   (3) REDUCED SPIKE-FREQUENCY ADAPTATION (lower ISI_last/ISI_first) <- gM down
%
% It also checks SPECIFICITY (passive leak / Kdr repolarization / resting
% potential are NOT changed, so the phenotype is purely intrinsic
% excitability) and MECHANISTIC ATTRIBUTION (each knob, applied alone,
% produces its expected effect).
%
% Run:   octave-cli --eval test_L5_Soma_SHANK3   (exits non-zero on failure)
% or in MATLAB/Octave:   test_L5_Soma_SHANK3
% =========================================================================
function test_L5_Soma_SHANK3()

more off;
fprintf('\n===== L5 SOMA (SHANK3-KO) DIAGNOSTIC TEST =====\n\n');

results = {};   % {name, passed(bool), detail}

%% ---- Baseline runs at the documented protocol (45 uA/cm^2, 200 ms) -----
wt = L5_Soma_SHANK3(struct('state','WT'));
ko = L5_Soma_SHANK3(struct('state','SHANK3'));

fprintf('Protocol: %g uA/cm^2 step, %g-%g ms.\n', ...
        wt.params.I_step, wt.params.stim_on, wt.params.stim_off);
fprintf('  WT     : Vrest=%6.2f mV  Vpeak=%5.1f mV  nSpikes=%2d  ISIratio=%.3f\n', ...
        wt.Vrest, wt.Vpeak, wt.nSpikes, wt.ISI_ratio);
fprintf('  SHANK3 : Vrest=%6.2f mV  Vpeak=%5.1f mV  nSpikes=%2d  ISIratio=%.3f\n\n', ...
        ko.Vrest, ko.Vpeak, ko.nSpikes, ko.ISI_ratio);

%% ---- A. PHYSIOLOGICAL SANITY (is the soma a valid spike generator?) ----
results(end+1,:) = check('WT resting potential physiological (-82..-66 mV)', ...
    wt.Vrest > -82 && wt.Vrest < -66, sprintf('Vrest=%.2f mV', wt.Vrest));

results(end+1,:) = check('WT action potentials overshoot (Vpeak > +10 mV)', ...
    wt.Vpeak > 10, sprintf('Vpeak=%.1f mV', wt.Vpeak));

results(end+1,:) = check('WT fires repetitively under step (nSpikes >= 5)', ...
    wt.nSpikes >= 5, sprintf('nSpikes=%d', wt.nSpikes));

results(end+1,:) = check('SHANK3 fires repetitively under step (nSpikes >= 5)', ...
    ko.nSpikes >= 5, sprintf('nSpikes=%d', ko.nSpikes));

results(end+1,:) = check('WT shows spike-frequency adaptation (ISIratio >= 1.0)', ...
    wt.ISI_ratio >= 1.0 - 1e-6, sprintf('ISIratio=%.3f', wt.ISI_ratio));

%% ---- B. SHANK3 PATHOLOGY SIGNATURES (the three documented readouts) ----
rb_wt = rheobase('WT');
rb_ko = rheobase('SHANK3');
fprintf('Rheobase:  WT=%.1f uA/cm^2   SHANK3=%.1f uA/cm^2\n\n', rb_wt, rb_ko);

% (1) Hyperexcitability: lower rheobase under SHANK3 loss.
results(end+1,:) = check('SIGNATURE 1  lower rheobase (SHANK3 < WT)', ...
    rb_ko < rb_wt - 1e-9, sprintf('WT=%.1f  SHANK3=%.1f', rb_wt, rb_ko));

% (2) Steeper F-I: more spikes for the same 45 uA/cm^2 step.
results(end+1,:) = check('SIGNATURE 2  higher spike count (SHANK3 > WT, +>=2)', ...
    ko.nSpikes >= wt.nSpikes + 2, sprintf('WT=%d  SHANK3=%d', wt.nSpikes, ko.nSpikes));

% (3) Reduced spike-frequency adaptation: lower ISI_last/ISI_first.
results(end+1,:) = check('SIGNATURE 3  reduced adaptation (SHANK3 ISIratio < WT, gap>=0.02)', ...
    (wt.ISI_ratio - ko.ISI_ratio) >= 0.02, ...
    sprintf('WT=%.3f  SHANK3=%.3f  gap=%.3f', wt.ISI_ratio, ko.ISI_ratio, wt.ISI_ratio-ko.ISI_ratio));

%% ---- C. PARAMETER WIRING (the two knobs are actually set) --------------
results(end+1,:) = check('KNOB 1  M-current down-regulated (gM: SHANK3 < WT, ~33%)', ...
    ko.params.gM < wt.params.gM && abs(ko.params.gM/wt.params.gM - 0.667) < 0.1, ...
    sprintf('gM WT=%.2f  SHANK3=%.2f  (%.0f%% of WT)', ...
            wt.params.gM, ko.params.gM, 100*ko.params.gM/wt.params.gM));

results(end+1,:) = check('KNOB 2  Na activation hyperpolarized (na_vshift < 0, ~-3 mV)', ...
    ko.params.na_vshift < 0 && abs(ko.params.na_vshift + 3) < 0.5, ...
    sprintf('na_vshift=%.1f mV', ko.params.na_vshift));

%% ---- D. SPECIFICITY (no confounding changes to passive / Kdr) ----------
results(end+1,:) = check('SPECIFICITY  passive/Kdr unchanged (gKdr, gL, E_L identical)', ...
    wt.params.gKdr==ko.params.gKdr && wt.params.gL==ko.params.gL && wt.params.E_L==ko.params.E_L, ...
    sprintf('gKdr=%.1f gL=%.2f E_L=%.0f', ko.params.gKdr, ko.params.gL, ko.params.E_L));

results(end+1,:) = check('SPECIFICITY  resting potential ~unchanged (|dVrest| < 2 mV)', ...
    abs(wt.Vrest - ko.Vrest) < 2.0, sprintf('|dVrest|=%.2f mV', abs(wt.Vrest-ko.Vrest)));

%% ---- E. MECHANISTIC ATTRIBUTION (each knob alone -> its own effect) ----
% Na shift ALONE (keep WT gM) must lower the rheobase.
rb_na = rheobase_opts(struct('state','WT','na_vshift',-3));
results(end+1,:) = check('MECHANISM  Na -3 mV shift alone lowers rheobase vs WT', ...
    rb_na < rb_wt - 1e-9, sprintf('WT=%.1f  Na-shift-only=%.1f', rb_wt, rb_na));

% gM reduction ALONE (no Na shift) must reduce adaptation & raise spike count.
gmonly = L5_Soma_SHANK3(struct('state','WT','gM',1.0));
results(end+1,:) = check('MECHANISM  gM 1.5->1.0 alone reduces adaptation vs WT', ...
    gmonly.ISI_ratio < wt.ISI_ratio - 1e-6, ...
    sprintf('WT ISIratio=%.3f  gM-only=%.3f', wt.ISI_ratio, gmonly.ISI_ratio));

results(end+1,:) = check('MECHANISM  gM 1.5->1.0 alone raises spike count vs WT', ...
    gmonly.nSpikes >= wt.nSpikes, ...
    sprintf('WT n=%d  gM-only n=%d', wt.nSpikes, gmonly.nSpikes));

%% ---- SUMMARY -----------------------------------------------------------
fprintf('\n----------------------------------------------------------------\n');
nPass = 0; nFail = 0;
for i = 1:size(results,1)
    name = results{i,1}; passed = results{i,2}; detail = results{i,3};
    if passed, tag = 'PASS'; nPass=nPass+1; else, tag = 'FAIL'; nFail=nFail+1; end
    fprintf('  [%s] %s   (%s)\n', tag, name, detail);
end
fprintf('----------------------------------------------------------------\n');
fprintf('  %d passed, %d failed, %d total\n\n', nPass, nFail, nPass+nFail);

if nFail > 0
    error('test_L5_Soma_SHANK3:FAILED', '%d diagnostic check(s) failed.', nFail);
else
    fprintf('ALL DIAGNOSTIC CHECKS PASSED: soma compartment reproduces SHANK3-KO pathology.\n\n');
end

end % ===================== end main test function =======================


%% ======================= LOCAL HELPER FUNCTIONS =========================

function row = check(name, passed, detail)
    if nargin < 3, detail = ''; end
    row = {name, logical(passed), detail};
end

function rb = rheobase(state)
    rb = rheobase_opts(struct('state', state));
end

function rb = rheobase_opts(baseopts)
    % Minimal step current (uA/cm^2) that evokes >=1 spike during a 300 ms
    % sustained step. Coarse 0.5 uA/cm^2 resolution is enough to separate
    % the WT and SHANK3 thresholds.
    rb = NaN;
    for I = 0.5:0.5:40
        o = L5_Soma_SHANK3(withfields(baseopts, 'I_step', I, ...
                'stim_on', 100, 'stim_off', 400, 't_max', 500));
        if numel(o.spike_stim) >= 1
            rb = I; return;
        end
    end
end

function s = withfields(s, varargin)
    % Set multiple fields at once: withfields(s,'a',1,'b',2).
    for i = 1:2:numel(varargin)
        s.(varargin{i}) = varargin{i+1};
    end
end
