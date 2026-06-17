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
- `TRN_TC_Delayed_IPSP_Experiment.m` - compares a normal TRN-to-TC IPSP with a
  delayed IPSP to test whether TC activity still decreases after delayed
  inhibition arrives.

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

## Running the delayed IPSP experiment

Open MATLAB and run:

```matlab
TRN_TC_Delayed_IPSP_Experiment
```

Change `ipsp_delay_ms` near the top of the script to test a different delay
between TRN release and GABA_A/GABA_B arrival at the TC cell. The script
compares a 0 ms delay condition against the delayed condition, then plots:

1. TC voltage with no-delay versus delayed inhibition.
2. TRN release events and delayed IPSP arrival times.
3. GABA_A/GABA_B inhibitory current components at the TC cell.
4. Mean TC voltage and TC spike counts before versus after IPSP arrival.
5. Zoomed TC responses around each delayed IPSP.