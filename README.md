# Autism_Computational_Research

MATLAB scripts for simple thalamic and thalamic reticular nucleus (TRN)
microcircuit simulations.

## Scripts

- `Minimal_Passive_TRN_Loop_EPSP_Only.m` - original minimal TC-to-TRN
  EPSP-only loop with a passive TRN neuron.
- `TRN_Thalamic_Microcircuit.m` - reciprocal TC/TRN microcircuit where the
  thalamic relay neuron excites TRN through AMPA and NMDA synaptic components,
  and the TRN neuron inhibits the thalamic relay neuron through GABA_A and
  GABA_B synaptic components.

## Running the reciprocal microcircuit

Open MATLAB and run:

```matlab
TRN_Thalamic_Microcircuit
```

By default, the script uses a staggered release protocol so TC-to-TRN
AMPA/NMDA pulses occur before TRN-to-TC GABA_A/GABA_B pulses. In the script,
set `use_staggered_release_protocol = false` to make transmitter release occur
only from presynaptic threshold-crossing spikes.

The script prints a short spike/current summary and opens several figures:

1. External drive, membrane voltages, and staggered release events.
2. TC-to-TRN glutamate and TRN-to-TC GABA neurotransmitter pulses.
3. AMPA/NMDA and GABA_A/GABA_B receptor or channel state variables.
4. Excitatory and inhibitory synaptic current components.
5. A phase-plane view of the coupled TC/TRN voltage trajectory.