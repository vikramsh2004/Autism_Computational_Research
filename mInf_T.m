function val = mInf_T(V)
% Steady-state T-type Ca2+ activation (m_inf,T).
val = 1 ./ (1 + exp(-(V + 52) / 7.4));
end
