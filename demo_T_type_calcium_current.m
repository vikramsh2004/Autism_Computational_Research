%% Minimal T-type Calcium Current (I_T) Demonstration
% Runs only the I_T kinetics snippet in isolation.
% Two modes:
%   1) Voltage-clamp protocol (prescribed V) — clearest view of I_T kinetics
%   2) Simple passive membrane + I_T + I_H — shows rebound depolarization

clear; clc

demo_mode = 'voltage_clamp'; % 'voltage_clamp' | 'passive_rebound'

%% Shared parameters
dt = 0.01;          % ms
t_end = 500;        % ms
time = 0:dt:t_end;
nt = numel(time);

gT = 1.0;           % mS/cm^2
ECa = 120;          % mV
gH = 0.05;          % mS/cm^2 (optional; helps rebound in passive mode)
EH = -43;           % mV
Cm = 1.0;           % uF/cm^2
gL = 0.3;           % mS/cm^2 (passive mode only)
EL = -54.4;         % mV

V = zeros(1, nt);
ht = zeros(1, nt);
r = zeros(1, nt);
I_T = zeros(1, nt);
I_H = zeros(1, nt);

switch demo_mode
    case 'voltage_clamp'
        % Classic low-threshold Ca2+ window test:
        % hyperpolarize to remove h-inactivation, then step into activation range
        V(:) = -70;
        V(time >= 100 & time < 300) = -90;
        V(time >= 300 & time < 450) = -50;

        ht(1) = hInf_T(V(1));
        r(1) = rInf_H(V(1));

        for i = 1:nt-1
            I_T(i) = gT * mInf_T(V(i))^2 * ht(i) * (ECa - V(i));
            I_H(i) = gH * r(i) * (EH - V(i));

            ht(i+1) = ht(i) + dt * ((hInf_T(V(i)) - ht(i)) / tauH_T(V(i)));
            r(i+1) = rInf_H(V(i)) + (r(i) - rInf_H(V(i))) * exp(-dt / tauR_H(V(i)));
        end

    case 'passive_rebound'
        % Brief inhibitory-like hyperpolarization, then release into I_T rebound
        V(1) = -65;
        ht(1) = hInf_T(V(1));
        r(1) = rInf_H(V(1));

        I_ext = zeros(1, nt);
        I_ext(time >= 150 & time < 250) = -12; % hyperpolarizing pulse

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

fprintf('Peak I_T: %.3f uA/cm^2 at t = %.1f ms\n', max(I_T), time(I_T == max(I_T)));

%% Plot
figure('Name', 'T-type Ca2+ Current Demo', 'Color', 'w', 'Position', [100 100 750 700]);

subplot(4, 1, 1)
plot(time, V, 'k', 'LineWidth', 1.4)
ylabel('V (mV)')
title(sprintf('T-type Ca^{2+} demo (%s)', demo_mode))
grid on

subplot(4, 1, 2)
plot(time, I_T, 'r', 'LineWidth', 1.4)
ylabel('I_T (\muA/cm^2)')
title('T-type calcium current')
grid on

subplot(4, 1, 3)
plot(time, ht, 'b', 'LineWidth', 1.4); hold on
plot(time, mInf_T(V), 'r--', 'LineWidth', 1.2)
ylabel('Gate')
title('h_T (solid) and m_{\infty,T}(V) (dashed)')
legend('h_T', 'm_{\infty,T}', 'Location', 'best')
grid on

subplot(4, 1, 4)
plot(time, I_H, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.2)
xlabel('Time (ms)')
ylabel('I_H (\muA/cm^2)')
title('H-current (included from same snippet)')
grid on

%% Local functions (your snippet)
function val = mInf_T(V)
    val = 1 / (1 + exp(-(V + 52) / 7.4));
end

function val = hInf_T(V)
    val = 1 / (1 + exp((V + 80) / 5));
end

function val = tauH_T(V)
    if V < -80
        val = exp((V + 467) / 66.6);
    else
        val = 28 + exp(-(V + 22) / 10.5);
    end
end

function val = rInf_H(V)
    val = 1 / (1 + exp((V + 75) / 5.5));
end

function tau = tauR_H(V)
    tau = 20 + 200 / (exp((V + 75) / 10) + exp(-(V + 75) / 10));
end
