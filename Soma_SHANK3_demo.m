function Soma_SHANK3_demo(outfile)
%% DEMO: Wild-type vs SHANK3-KO isolated somatic compartment
% =========================================================================
% Visualises the three SHANK3-KO somatic pathology signatures produced by
% L5_Soma_SHANK3.m:
%   (a) Membrane potential during a 45 uA/cm^2 step  (hyperexcitable, faster firing)
%   (b) F-I curve: spike count vs injected current   (lower rheobase, steeper slope)
%   (c) Instantaneous firing rate vs spike index     (reduced spike-freq. adaptation)
%
% Usage:  Soma_SHANK3_demo                 % shows the figure
%         Soma_SHANK3_demo('soma.png')     % also saves a PNG (headless OK)
% =========================================================================
if nargin < 1, outfile = ''; end

wt = L5_Soma_SHANK3(struct('state','WT'));
ko = L5_Soma_SHANK3(struct('state','SHANK3'));

% F-I sweep (light sims: only need up to stim_off for the spike count)
Isweep = 0:4:80;
nWT = zeros(size(Isweep)); nKO = zeros(size(Isweep));
for i = 1:numel(Isweep)
    a = L5_Soma_SHANK3(struct('state','WT',    'I_step',Isweep(i),'t_max',360));
    b = L5_Soma_SHANK3(struct('state','SHANK3','I_step',Isweep(i),'t_max',360));
    nWT(i) = a.nSpikes; nKO(i) = b.nSpikes;
end

cWT = [0.15 0.35 0.75];   % blue
cKO = [0.80 0.20 0.15];   % red

if ~isempty(outfile)
    % Render/print headlessly (no display / OpenGL) via the gnuplot toolkit,
    % bound directly to this figure handle.
    fig = figure('visible','off','Position',[40 40 950 950]);
    graphics_toolkit(fig, 'gnuplot');
else
    fig = figure('Position',[40 40 950 950]);
end

% (a) Voltage traces
subplot(3,1,1); hold on;
plot(wt.t, wt.V, 'Color', cWT, 'LineWidth', 1.0);
plot(ko.t, ko.V, 'Color', cKO, 'LineWidth', 1.0);
plot([wt.params.stim_on wt.params.stim_on], [-85 55], 'k:', 'HandleVisibility','off');
plot([wt.params.stim_off wt.params.stim_off], [-85 55], 'k:', 'HandleVisibility','off');
title('Somatic membrane potential: 45 \muA/cm^2 step (WT vs SHANK3-KO)');
xlabel('Time (ms)'); ylabel('V_m (mV)');
legend({sprintf('WT (n=%d, ISIr=%.2f)',wt.nSpikes,wt.ISI_ratio), ...
        sprintf('SHANK3 (n=%d, ISIr=%.2f)',ko.nSpikes,ko.ISI_ratio)}, ...
        'Location','northeast');
grid on; xlim([wt.params.stim_on-20 wt.params.stim_off+20]);

% (b) F-I curve
subplot(3,1,2); hold on;
plot(Isweep, nWT, '-o', 'Color', cWT, 'MarkerFaceColor', cWT, 'MarkerSize',3);
plot(Isweep, nKO, '-o', 'Color', cKO, 'MarkerFaceColor', cKO, 'MarkerSize',3);
title('F-I curve: lower rheobase & steeper slope under SHANK3-KO');
xlabel('Injected current (\muA/cm^2)'); ylabel('Spikes per 200 ms step');
legend({'WT','SHANK3-KO'}, 'Location','northwest'); grid on;

% (c) Adaptation: instantaneous rate vs spike index
subplot(3,1,3); hold on;
fWT = 1000 ./ wt.ISI; fKO = 1000 ./ ko.ISI;
plot(1:numel(fWT), fWT, '-o', 'Color', cWT, 'MarkerFaceColor', cWT, 'MarkerSize',3);
plot(1:numel(fKO), fKO, '-o', 'Color', cKO, 'MarkerFaceColor', cKO, 'MarkerSize',3);
title('Spike-frequency adaptation: WT slows (adapts), SHANK3-KO sustains firing');
xlabel('Inter-spike interval index'); ylabel('Instantaneous rate (Hz)');
legend({'WT','SHANK3-KO'}, 'Location','northeast'); grid on;

if ~isempty(outfile)
    print(fig, outfile, '-dpng', '-r110');
    fprintf('Saved figure to %s\n', outfile);
    close(fig);
end
end
