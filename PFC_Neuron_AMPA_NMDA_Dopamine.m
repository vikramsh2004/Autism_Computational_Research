% PFC Neuron: AMPA/NMDA Synaptic Inputs (Passive Membrane with Dopamine Pooling)
% Re-running this script overlays the new traces on top of the existing
% figure instead of opening a fresh window each time.

% Use rng('shuffle') so each run produces a different Poisson realization;
% with a fixed seed every run would be identical and overlays would coincide.
rng('shuffle');

% 1. SIMULATION PARAMETERS
dt = 0.01;                 % Time step (ms)
t_max = 300;               % Total simulation time (ms)
t = 0:dt:t_max;            
N = length(t);             

% 2. INPUT PARAMETERS (Poisson Rates)
r_pre_Glu = 30;            % Background Glutamate rate (Hz)
r_pre_DA = 40;             % Background Dopamine firing rate (Hz)
release_steps = round(1/dt); % release duration
counter_Glu = 0;           
counter_DA = 0;

% 3. MEMBRANE & CHANNEL PARAMETERS
C = 1.0;                   
g_L = 0.1;   E_L = -65;    
g_Na = 0;    E_Na = 50;    % Set to 0 
g_K  = 0;    E_K = -77;    % Set to 0
E_syn = 0;   Mg = 1.0;                  

% Synaptic Base Conductances
g_AMPA = 0.30;  alpha_AMPA = 1.1;   beta_AMPA = 0.5;         
g_NMDA = 0.20;  alpha_NMDA = 0.072; beta_NMDA = 0.015;       

% Dopamine (D1) Parameters
B_D1 = 0.1;                % Scaling coefficient for DA release
tau = 50;                  % Decay time constant for D1 (ms)
A_D1 = 1 / tau;            % Decay rate

% 4. INITIALIZE VARIABLES
V = zeros(1, N);
V(1) = E_L;                
r_AMPA = zeros(1, N);
r_NMDA = zeros(1, N);
I_AMPA = zeros(1, N);
I_NMDA = zeros(1, N);
D1 = zeros(1, N);          % Dopamine pooled concentration (nM)

% 5. MAIN INTEGRATION LOOP
for i = 1:N-1
    % INPUT GENERATORS (Glutamate)
    if poissrnd(r_pre_Glu * (dt/1000)) > 0
        counter_Glu = release_steps;
    end
    T_Glu = double(counter_Glu > 0);
    
    % INPUT GENERATORS (Dopamine)
    if poissrnd(r_pre_DA * (dt/1000)) > 0
        counter_DA = release_steps;
    end
    T_D = double(counter_DA > 0);
    
    % DOPAMINE INPUT AND LEAKY INTEGRATION (POOLING)
    dD1dt = -A_D1 * D1(i) + B_D1 * T_D;
    D1(i+1) = D1(i) + dt * dD1dt; 
    D1mod = 1 + D1(i);
    
    % SYNAPTIC KINETICS
    r_AMPA(i+1) = r_AMPA(i) + dt * ((alpha_AMPA * T_Glu * (1 - r_AMPA(i))) - (beta_AMPA * r_AMPA(i)));
    r_NMDA(i+1) = r_NMDA(i) + dt * ((alpha_NMDA * T_Glu * (1 - r_NMDA(i))) - (beta_NMDA * r_NMDA(i)));
    
    % SYNAPTIC CURRENTS (Modulated by D1mod)
    B_V = 1 / (1 + (Mg / 3.57) * exp(-0.062 * V(i)));
    
    I_AMPA(i) = g_AMPA * r_AMPA(i) * (V(i) - E_syn) * D1mod;
    I_NMDA(i) = g_NMDA * r_NMDA(i) * B_V * (V(i) - E_syn) * D1mod; 
    I_Total_syn = I_AMPA(i) + I_NMDA(i);
    
    % MEMBRANE VOLTAGE UPDATE (Passive Leak + Synaptic Currents Only)
    I_leak = g_L * (V(i) - E_L);
    dV = (-I_leak - I_Total_syn) / C;
    V(i+1) = V(i) + dt * dV;
    
    % Decrement counters
    counter_Glu = max(0, counter_Glu - 1);
    counter_DA = max(0, counter_DA - 1);
end
I_AMPA(end) = I_AMPA(end-1); I_NMDA(end) = I_NMDA(end-1);
D1mod_plot = 1 + D1; 

% 6. PLOTTING
% Reuse a single figure across runs so successive runs overlay. The figure
% is located by name; if it does not exist yet it is created.
figName = 'Passive EPSCs & Dopamine Pooling';
fig = findobj('Type', 'figure', 'Name', figName);
if isempty(fig)
    fig = figure('Name', figName, 'Position', [100, 50, 750, 950]);
    setappdata(fig, 'runCount', 0);
else
    fig = fig(1);
    figure(fig);   % bring existing figure to the front / make current
end

% Track how many times the script has been run against this figure and pick
% a distinct color for each run so overlaid traces are distinguishable.
runCount = getappdata(fig, 'runCount');
if isempty(runCount)
    runCount = 0;
end
runCount = runCount + 1;
setappdata(fig, 'runCount', runCount);

colorOrder = lines(7);
runColor = colorOrder(mod(runCount - 1, size(colorOrder, 1)) + 1, :);
runLabel = sprintf('Run %d', runCount);

% Panel 1: Sub-threshold Membrane Potential (EPSPs)
subplot(4,1,1);
hold on;
plot(t, V, 'Color', runColor, 'LineWidth', 1.2, 'DisplayName', runLabel);
title('Membrane Potential (Sub-threshold EPSPs)');
ylabel('Voltage (mV)');
legend('show');
grid on;

% Panel 2: Dopamine Pooling Over Time
subplot(4,1,2);
hold on;
plot(t, D1, 'Color', runColor, 'LineWidth', 1.2, 'DisplayName', runLabel);
title('Dopamine Pooling Concentration (D1)');
ylabel('Concentration (nM)');
legend('show');
grid on;

% Panel 3: AMPA EPSC
subplot(4,1,3);
hold on;
plot(t, I_AMPA, 'Color', runColor, 'LineWidth', 1.2, 'DisplayName', runLabel);
title('AMPA EPSC');
ylabel('Current');
legend('show');
grid on;

% Panel 4: NMDA EPSC
subplot(4,1,4);
hold on;
plot(t, I_NMDA, 'Color', runColor, 'LineWidth', 1.2, 'DisplayName', runLabel);
title('NMDA EPSC');
xlabel('Time (ms)');
ylabel('Current');
legend('show');
grid on;
