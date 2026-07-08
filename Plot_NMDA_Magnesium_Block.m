%% NMDA receptor magnesium block from Susin and Destexhe (2023)
% This script plots only the voltage-dependent magnesium block factor B(V)
% used for the NMDA receptor in the paper:
%
%   B(V) = 1 / (1 + exp(-0.062*V) * ([Mg2+]_o / 3.57))
%
% where V is membrane voltage in mV and [Mg2+]_o is the external magnesium
% concentration. The paper uses [Mg2+]_o = 1 mM.

clear; clc; close all;

%% Paper parameters
Mg_o = 1.0;          % external magnesium concentration (mM)
V = -100:0.1:60;    % membrane voltage range (mV)

%% Magnesium block factor
B = 1 ./ (1 + exp(-0.062 .* V) .* (Mg_o ./ 3.57));

%% Plot
figure('Color', 'w');
plot(V, B, 'LineWidth', 2);
grid on;
box on;

xlabel('Membrane voltage, V (mV)');
ylabel('NMDA magnesium block factor, B(V)');
title('Voltage-dependent Mg^{2+} block of NMDA receptors');
ylim([0 1]);
xlim([min(V) max(V)]);
