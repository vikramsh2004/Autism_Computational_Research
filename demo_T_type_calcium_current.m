%% Minimal T-type Calcium Current (I_T) Demonstration
% Runs only the I_T kinetics snippet in isolation.
% Two modes:
%   1) Voltage-clamp protocol (prescribed V) — no external current
%   2) Simple passive membrane + I_T + I_H — hyperpolarizing I_ext pulse

clear; clc

demo_mode = 'voltage_clamp'; % 'voltage_clamp' | 'passive_rebound'

%% Shared parameters
dt = 0.01;          % ms
t_end = 500;        % ms
time = 0:dt:t_end;
nt = numel(time);

gT = 1.0;           % mS/cm^2
ECa = 120;          % mV
gH = 0.05;          % mS/cm^2
EH = -43;           % mV
Cm = 1.0;           % uF/cm^2
gL = 0.3;           % mS/cm^2 (passive mode only)
EL = -54.4;         % mV

V = zeros(1, nt);
ht = zeros(1, nt);
r = zeros(1, nt);
I_T = zeros(1, nt);
I_H = zeros(1, nt);
I_ext = zeros(1, nt); % zero in voltage-clamp mode

switch demo_mode
    case 'voltage_clamp'
        % Prescribed voltage steps (no external current injected)
        V_hold = -70;   V_hyper = -90;   V_dep = -50;
        t_hyper_on = 100;  t_hyper_off = 300;
        t_dep_on = 300;    t_dep_off = 450;

        V(:) = V_hold;
        V(time >= t_hyper_on & time < t_hyper_off) = V_hyper;
        V(time >= t_dep_on & time < t_dep_off) = V_dep;

        fprintf('--- Stimulus protocol: voltage_clamp ---\n');
        fprintf('External current (I_ext): NONE (voltage is prescribed/clamped)\n');
        fprintf('  0–%.0f ms:     hold at %.0f mV\n', t_hyper_on, V_hold);
        fprintf('  %.0f–%.0f ms: hyperpolarize to %.0f mV (h_T de-inactivates)\n', t_hyper_on, t_hyper_off, V_hyper);
        fprintf('  %.0f–%.0f ms: depolarize to %.0f mV (I_T activation window)\n', t_dep_on, t_dep_off, V_dep);
        fprintf('  %.0f–%.0f ms: return to %.0f mV\n\n', t_dep_off, t_end, V_hold);

        ht(1) = hInf_T(V(1));
        r(1) = rInf_H(V(1));

        for i = 1:nt-1
            I_T(i) = gT * mInf_T(V(i))^2 * ht(i) * (ECa - V(i));
            I_H(i) = gH * r(i) * (EH - V(i));

            ht(i+1) = ht(i) + dt * ((hInf_T(V(i)) - ht(i)) / tauH_T(V(i)));
            r(i+1) = rInf_H(V(i)) + (r(i) - rInf_H(V(i))) * exp(-dt / tauR_H(V(i)));
        end

    case 'passive_rebound'
        I_ext_amp = -12; % uA/cm^2, hyperpolarizing
        t_ext_on = 150;
        t_ext_off = 250;

        V(1) = -65;
        ht(1) = hInf_T(V(1));
        r(1) = rInf_H(V(1));

        I_ext(time >= t_ext_on & time < t_ext_off) = I_ext_amp;

        fprintf('--- Stimulus protocol: passive_rebound ---\n');
        fprintf('External current (I_ext): %.1f uA/cm^2 from %.0f to %.0f ms\n', I_ext_amp, t_ext_on, t_ext_off);
        fprintf('  (hyperpolarizing pulse — removes T-type inactivation, then rebound when pulse ends)\n');
        fprintf('Voltage is computed from I_ext + I_T + I_H + leak (no clamp).\n\n');

        for i = 1:nt-1
            I_T(i) = gT * mInf_T(V(i))^2 * ht(i) * (ECa - V(i));
            I_H(i) = gH * r(i) * (EH - V(i));
            IL = gL * (EL - V(i));

            V(i+1) = V(i) + dt * (I_ext(i) + I_T(i) + I_H(i) + IL) / Cm;

            ht(i+1) = ht(i) + dt * ((hInf_T(V(i)) - ht(i)) / tauH_T(V(i)));
            r(i+1) = rInf_H(V(i)) + (r(i) - rInf_H(V(i))) * exp(-dt / tauR_H(V(i)));
        end
end

I_T(end) = gT * mInf_T(V(end))^2 * ht(end) * (ECa - V(end));
I_H(end) = gH * r(end) * (EH - V(end));

fprintf('Peak I_T: %.3f uA/cm^2 at t = %.1f ms\n', max(I_T), time(find(I_T == max(I_T), 1, 'first')));

%% Plot
figure('Name', 'T-type Ca2+ Current Demo', 'Color', 'w', 'Position', [100 100 750 780]);

if strcmp(demo_mode, 'voltage_clamp')
    subplot(5, 1, 1)
    plot(time, V, 'k', 'LineWidth', 1.4)
    hold on
    shade_epoch(t_hyper_on, t_hyper_off, 'Hyperpolarize')
    shade_epoch(t_dep_on, t_dep_off, 'Depolarize')
    ylabel('V (mV)')
    title('Prescribed voltage protocol (I_{ext} = 0 throughout)')
    legend('V_{cmd}', 'Location', 'best')
    grid on

    row_V = 2;
else
    subplot(5, 1, 1)
    plot(time, I_ext, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4)
    hold on
    shade_epoch(t_ext_on, t_ext_off, sprintf('I_{ext} = %.1f', I_ext_amp))
    ylabel('I_{ext} (\muA/cm^2)')
    title(sprintf('External current: %.1f uA/cm^2, %.0f–%.0f ms', I_ext_amp, t_ext_on, t_ext_off))
    grid on

    subplot(5, 1, 2)
    plot(time, V, 'k', 'LineWidth', 1.4)
    hold on
    shade_epoch(t_ext_on, t_ext_off, 'I_{ext} on')
    ylabel('V (mV)')
    title('Membrane voltage (free evolution)')
    grid on

    row_V = 3;
end

subplot(5, 1, row_V)
plot(time, I_T, 'r', 'LineWidth', 1.4)
ylabel('I_T (\muA/cm^2)')
title('T-type calcium current')
grid on

subplot(5, 1, row_V + 1)
plot(time, ht, 'b', 'LineWidth', 1.4); hold on
plot(time, mInf_T(V), 'r--', 'LineWidth', 1.2)
ylabel('Gate')
title('h_T (solid) and m_{\infty,T}(V) (dashed)')
legend('h_T', 'm_{\infty,T}', 'Location', 'best')
grid on

subplot(5, 1, row_V + 2)
plot(time, I_H, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.2)
xlabel('Time (ms)')
ylabel('I_H (\muA/cm^2)')
title('H-current')
grid on

sgtitle(sprintf('T-type Ca^{2+} demo — %s mode', demo_mode), 'FontWeight', 'bold')

%% Local functions
function shade_epoch(t_on, t_off, label_text)
    yl = ylim;
    patch([t_on t_off t_off t_on], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.9 0.95 1.0], 'FaceAlpha', 0.35, 'EdgeColor', 'none')
    text((t_on + t_off) / 2, yl(2) - 0.08 * (yl(2) - yl(1)), label_text, ...
        'HorizontalAlignment', 'center', 'FontSize', 8)
    uistack(findobj(gca, 'Type', 'line'), 'top')
end

function val = mInf_T(V)
    val = 1 ./ (1 + exp(-(V + 52) / 7.4));
end

function val = hInf_T(V)
    val = 1 ./ (1 + exp((V + 80) / 5));
end

function val = tauH_T(V)
    val = 28 + exp(-(V + 22) / 10.5);
    val(V < -80) = exp((V(V < -80) + 467) / 66.6);
end

function val = rInf_H(V)
    val = 1 ./ (1 + exp((V + 75) / 5.5));
end

function tau = tauR_H(V)
    tau = 20 + 200 ./ (exp((V + 75) / 10) + exp(-(V + 75) / 10));
end
