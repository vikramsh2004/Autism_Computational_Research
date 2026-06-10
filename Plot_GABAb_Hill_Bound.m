%% Plot bounded GABA_B Hill/open-fraction term
clear; clc; % Keep existing figure windows open

%% GABA_B open-fraction constants
KD_GABAb = 100; % dissociation constant
n_GABAb = 4;   % number of G-protein binding sites

%% S range
S_half = KD_GABAb^(1/n_GABAb);
S = linspace(0, 5*S_half, 1000);

%% Bounded Hill/open-fraction term
GABAb_open = (S.^n_GABAb) ./ (S.^n_GABAb + KD_GABAb);
GABAb_open = min(max(GABAb_open, 0), 1);
GABAb_percent_open = 100 * GABAb_open;

%% Plot S^n/(S^n+Kd)
figure('Name','Bounded GABA_B Hill Term','Color','w');

plot(S, GABAb_percent_open, 'Color', [0.3 0.1 0.7], 'LineWidth', 1.5)
hold on
yline(100, 'k:', 'Upper bound = 100%');
xline(S_half, 'k--', 'Half activation');
hold off

title('Bounded GABA_B open fraction: S^n/(S^n+K_d)')
xlabel('S')
ylabel('GABA_B channels open (%)')
ylim([0 105])
grid on
