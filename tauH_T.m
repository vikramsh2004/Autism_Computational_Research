function val = tauH_T(V)
% T-type Ca2+ inactivation time constant (ms).
val = 28 + exp(-(V + 22) / 10.5);
val(V < -80) = exp((V(V < -80) + 467) / 66.6);
end
