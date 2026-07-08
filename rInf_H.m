function val = rInf_H(V)
% Steady-state H-current activation (Destexhe 1993).
val = 1 ./ (1 + exp((V + 75) / 5.5));
end
