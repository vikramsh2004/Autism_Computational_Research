# Autism_Computational_Research

Computational models of autism-associated channelopathies in cortical and
thalamic neurons (conductance-based, Hodgkin–Huxley style). Scripts are
written to run in **GNU Octave** (`octave-cli`) as well as MATLAB.

## L5 PFC pyramidal neuron — isolated somatic compartment (SHANK3-KO)

A single-compartment "unit model" of the Layer-5 prefrontal-cortex pyramidal
soma — the first building block of the larger 4-compartment cable. It is the
spike generator and the seat of spike-frequency adaptation (SFA), and it can
be switched between a wild-type and a **SHANK3 knock-out** state.

### Files

| File | Purpose |
|------|---------|
| `L5_Soma_SHANK3.m`      | The soma compartment model (`out = L5_Soma_SHANK3(opts)`). |
| `test_L5_Soma_SHANK3.m` | Diagnostic test: verifies the model and the SHANK3 pathology. |
| `Soma_SHANK3_demo.m`    | Generates a WT-vs-SHANK3 comparison figure. |

### Channels (per the model design document)

`I_Na` (Nav1.6, 120 mS/cm²), `I_Kdr` (Kv3.x, 36 mS/cm²),
`I_M` (Kv7/KCNQ, 1.5 mS/cm²), passive leak `I_L` (0.10 mS/cm²),
`Cm = 1.0 µF/cm²`, area `1570.8 µm²`.

### SHANK3-KO pathology (soma)

SHANK3 is a dendritic post-synaptic-density scaffold, so the soma has **no
direct** channelopathy. Its phenotype is a *maladaptive homeostatic*
compensation for lost dendritic drive, captured by two intrinsic knobs:

1. **M-current down-regulation** — `g_M`: 1.5 → 1.0 mS/cm² (~33 %) → reduced
   spike-frequency adaptation, steeper F–I curve.
2. **Hyperpolarizing Na⁺ activation shift** (AIS remodeling) — `V₁₍₂,m`:
   −40 → −43 mV → lower spike threshold / lower rheobase (hyperexcitable).

All passive and repolarization properties are held fixed, so the phenotype
is purely a change in intrinsic excitability.

### Usage

```matlab
wt = L5_Soma_SHANK3(struct('state','WT'));        % wild-type
ko = L5_Soma_SHANK3(struct('state','SHANK3'));    % SHANK3 knock-out
```

```bash
octave-cli --eval "test_L5_Soma_SHANK3"                     # run the diagnostic test
octave-cli --eval "Soma_SHANK3_demo('soma_shank3_demo.png')" # render the figure
```

### Diagnostic test

`test_L5_Soma_SHANK3.m` checks (and exits non-zero on any failure):

- **Physiological sanity** — physiological resting potential, overshooting
  action potentials, repetitive firing, WT exhibits SFA.
- **Three SHANK3 signatures** — lower rheobase, higher spike count at a fixed
  step, and reduced ISI adaptation ratio under SHANK3-KO.
- **Parameter wiring & specificity** — the two knobs are set, and the passive
  leak / delayed-rectifier / resting potential are unchanged.
- **Mechanistic attribution** — the Na⁺ shift alone lowers rheobase; the
  `g_M` reduction alone reduces adaptation and raises spike count.

The M-current activation time constant uses a physiological `τ ≈ 50 ms`
(`mM_k = 0.02`); set `opts.mM_k = 0.003` to reproduce the document's slower
value exactly (the steady-state activation curve is identical either way).
