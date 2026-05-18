"""Very basic Hodgkin-Huxley neuron firing simulation.

Reference:
https://neuronaldynamics.epfl.ch/online/Ch2.S2.html

This script simulates the classic 4D Hodgkin-Huxley model:
- Membrane potential V (mV)
- Sodium activation m
- Sodium inactivation h
- Potassium activation n
"""

from __future__ import annotations

import csv
import math


# Membrane and channel parameters (classic squid giant axon values)
CM = 1.0  # uF/cm^2
G_NA = 120.0  # mS/cm^2
G_K = 36.0  # mS/cm^2
G_L = 0.3  # mS/cm^2
E_NA = 50.0  # mV
E_K = -77.0  # mV
E_L = -54.387  # mV


def alpha_n(v: float) -> float:
    x = v + 55.0
    if abs(x) < 1e-6:
        return 0.1
    return 0.01 * x / (1.0 - math.exp(-x / 10.0))


def beta_n(v: float) -> float:
    return 0.125 * math.exp(-(v + 65.0) / 80.0)


def alpha_m(v: float) -> float:
    x = v + 40.0
    if abs(x) < 1e-6:
        return 1.0
    return 0.1 * x / (1.0 - math.exp(-x / 10.0))


def beta_m(v: float) -> float:
    return 4.0 * math.exp(-(v + 65.0) / 18.0)


def alpha_h(v: float) -> float:
    return 0.07 * math.exp(-(v + 65.0) / 20.0)


def beta_h(v: float) -> float:
    return 1.0 / (1.0 + math.exp(-(v + 35.0) / 10.0))


def injected_current(t: float) -> float:
    """Simple square pulse current (uA/cm^2)."""
    return 10.0 if 10.0 <= t <= 40.0 else 0.0


def steady_state_gate(v: float, alpha_fn, beta_fn) -> float:
    a = alpha_fn(v)
    b = beta_fn(v)
    return a / (a + b)


def run_simulation(t_stop: float = 50.0, dt: float = 0.01) -> dict[str, list[float]]:
    steps = int(t_stop / dt) + 1

    t_values: list[float] = []
    v_values: list[float] = []
    m_values: list[float] = []
    h_values: list[float] = []
    n_values: list[float] = []
    i_values: list[float] = []

    # Initial state near resting potential
    v = -65.0
    m = steady_state_gate(v, alpha_m, beta_m)
    h = steady_state_gate(v, alpha_h, beta_h)
    n = steady_state_gate(v, alpha_n, beta_n)

    for k in range(steps):
        t = k * dt
        i_ext = injected_current(t)

        i_na = G_NA * (m**3) * h * (v - E_NA)
        i_k = G_K * (n**4) * (v - E_K)
        i_l = G_L * (v - E_L)

        dv_dt = (i_ext - i_na - i_k - i_l) / CM
        dm_dt = alpha_m(v) * (1.0 - m) - beta_m(v) * m
        dh_dt = alpha_h(v) * (1.0 - h) - beta_h(v) * h
        dn_dt = alpha_n(v) * (1.0 - n) - beta_n(v) * n

        v += dt * dv_dt
        m += dt * dm_dt
        h += dt * dh_dt
        n += dt * dn_dt

        t_values.append(t)
        v_values.append(v)
        m_values.append(m)
        h_values.append(h)
        n_values.append(n)
        i_values.append(i_ext)

    return {
        "t": t_values,
        "V": v_values,
        "m": m_values,
        "h": h_values,
        "n": n_values,
        "I": i_values,
    }


def count_spikes(v_trace: list[float], threshold: float = 0.0) -> int:
    spikes = 0
    for prev_v, cur_v in zip(v_trace, v_trace[1:]):
        if prev_v < threshold <= cur_v:
            spikes += 1
    return spikes


def write_csv(data: dict[str, list[float]], output_file: str = "hh_trace.csv") -> None:
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["t_ms", "V_mV", "m", "h", "n", "I_uA_per_cm2"])
        for row in zip(data["t"], data["V"], data["m"], data["h"], data["n"], data["I"]):
            writer.writerow(row)


def try_plot(data: dict[str, list[float]]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping plot.")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    ax1.plot(data["t"], data["V"], color="navy", lw=1.4)
    ax1.set_ylabel("V (mV)")
    ax1.set_title("Basic Hodgkin-Huxley Simulation")

    ax2.plot(data["t"], data["I"], color="darkred", lw=1.2)
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("I_ext (uA/cm^2)")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    sim_data = run_simulation(t_stop=50.0, dt=0.01)
    spikes = count_spikes(sim_data["V"])
    peak_v = max(sim_data["V"])
    min_v = min(sim_data["V"])

    print(f"Simulated {len(sim_data['t'])} steps.")
    print(f"Estimated spike count (0 mV crossing): {spikes}")
    print(f"Peak V: {peak_v:.2f} mV, Min V: {min_v:.2f} mV")

    write_csv(sim_data, "hh_trace.csv")
    print("Saved trace to hh_trace.csv")

    try_plot(sim_data)
