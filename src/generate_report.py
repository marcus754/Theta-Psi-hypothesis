
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import sys
sys.path.append(str(Path(__file__).resolve().parents[1]))

from background.frw_integrator import Params, pack_state, integrate, unpack_state, theta_from_x, theta_dot_from_x

def run_demo_and_plot(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    H0 = 0.05
    p = Params(gamma=1.2, V0=0.2, H_of_t=lambda t: H0)
    y0 = pack_state(0.2, 0.02, 0.05, 0.0)
    traj = integrate(0.0, 25.0, y0, p, h0=5e-3)

    t_list, theta_list, psi_list = [], [], []
    for t, y in traj:
        x, xdot, psi, psidot = unpack_state(y, p)
        th = theta_from_x(x, p)
        t_list.append(t); theta_list.append(th); psi_list.append(psi)

    # Save CSV
    csv_path = outdir / "background_demo.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["t","theta","psi"])
        for t,th,ps in zip(t_list, theta_list, psi_list):
            w.writerow([t, th, ps])

    # Plot (single plot, no colors specified)
    plt.figure()
    plt.plot(t_list, theta_list, label="theta(t)")
    plt.plot(t_list, psi_list, label="Psi(t)")
    plt.xlabel("t")
    plt.ylabel("fields")
    plt.legend()
    fig_path = outdir / "background_fields.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    return csv_path, fig_path

if __name__ == "__main__":
    outdir = Path("results")
    paths = run_demo_and_plot(outdir)
    print("Generated:", paths)
