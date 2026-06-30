function tau = tauR_H(V)
% H-current activation time constant (ms, Destexhe 1993).
tau = 20 + 200 ./ (exp((V + 75) / 10) + exp(-(V + 75) / 10));
end
