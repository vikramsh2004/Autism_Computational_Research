function val = hInf_T(V)
% Steady-state T-type Ca2+ inactivation (h_inf,T).
val = 1 ./ (1 + exp((V + 80) / 5));
end
