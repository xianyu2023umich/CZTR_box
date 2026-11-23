#!/usr/bin/env python3
"""The magnetic field's evolution in a CZTR_box run.

Usage: /usr/bin/python3 tools/mag_evolution.py OUT_DIR out.png

From log.txt: max |B|, E_mag/E_kin and rms(div B dx)/rms|B| against time;
from the profiles: rms B_z and rms B_h at several heights against time, and
their profiles against z at selected times (colour: time). Prints the
growth rate of E_mag/E_kin (e-folding time) over the run's halves.
"""
import glob
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    d, png = sys.argv[1], sys.argv[2]
    log = np.loadtxt(f'{d}/log.txt', comments='#')
    log = log[log[:, 2] > 0]                                # (the first line: dt = 0)
    t = log[:, 1]/86400
    bmax, er, divb = log[:, 20], log[:, 21], log[:, 22]
    prof = []
    for f in sorted(glob.glob(f'{d}/profile_[0-9][0-9][0-9][0-9].dat')):
        with open(f) as fh:
            tt = float(fh.readline().split('=')[1].split()[0])
        prof.append((tt/86400, np.loadtxt(f, comments='#')))
    z = prof[0][1][:, 2]/1e3
    n = len(t)
    for a, b in ((0, n//2), (n//2, n - 1)):
        if er[a] > 0 and er[b] > 0 and t[b] > t[a]:
            g = np.log(er[b]/er[a])/(t[b] - t[a])
            print(f't {t[a]:.1f}-{t[b]:.1f} d: E_mag/E_kin {er[a]:.3e} -> {er[b]:.3e}, rate {g:.3f}/d '
                  f'(e-folding {1/g if g != 0 else np.inf:.2f} d)')
    print(f'max |B| {bmax[0]:.1f} -> {bmax[-1]:.1f} G; rms(div B dx)/rms|B| last {divb[-1]:.2e}')

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    ax[0, 0].semilogy(t, er); ax[0, 0].set_ylabel('E_mag / E_kin'); ax[0, 0].set_title('magnetic / kinetic energy')
    ax[0, 1].semilogy(t, bmax); ax[0, 1].set_ylabel('max |B| [G]'); ax[0, 1].set_title('maximum field')
    heights = [z[-1] - 2.0, z[-1] - 15.0, -80.0, -130.0, -180.0]
    cm = plt.get_cmap('plasma')
    tp = np.array([p[0] for p in prof])
    for q, zz in enumerate(heights):
        i = int(np.argmin(abs(z - zz)))
        c = cm(q/len(heights))
        ax[0, 2].semilogy(tp, [max(p[1][i, 21], 1e-6) for p in prof], color=c, label=f'B_h, z = {z[i]:.0f} Mm')
        ax[0, 2].semilogy(tp, [max(p[1][i, 20], 1e-6) for p in prof], color=c, ls='--')
    ax[0, 2].set_ylabel('rms B [G] (solid B_h, dashed B_z)'); ax[0, 2].set_title('rms field by height')
    for a in ax[0]:
        a.set_xlabel('t [d]'); a.grid(alpha=0.3)
    ax[0, 2].legend(fontsize=7)
    sel = np.unique(np.linspace(0, len(prof) - 1, 8).astype(int))
    cm = plt.get_cmap('viridis')
    for k, s in enumerate(sel):
        tt, a = prof[s]
        c = cm(k/max(len(sel) - 1, 1))
        ax[1, 0].semilogx(np.maximum(a[:, 21], 1e-6), z, color=c, label=f'{tt:.0f} d')
        ax[1, 1].semilogx(np.maximum(a[:, 20], 1e-6), z, color=c, label=f'{tt:.0f} d')
        # the field's energy over the flow's (rho0 based, per height)
        eb = (a[:, 20]**2 + a[:, 21]**2)/(8*np.pi)
        ek = 0.5*a[:, 3]*(a[:, 8]**2 + a[:, 9]**2)
        ax[1, 2].semilogx(np.maximum(eb/np.maximum(ek, 1e-30), 1e-12), z, color=c, label=f'{tt:.0f} d')
    ax[1, 0].set_xlabel('rms B_h [G]'); ax[1, 1].set_xlabel('rms B_z [G]'); ax[1, 2].set_xlabel('E_mag / E_kin per height')
    for a in ax[1]:
        a.set_ylabel('z [Mm]'); a.grid(alpha=0.3); a.legend(fontsize=7)
    fig.suptitle(f'{d}: the magnetic field')
    fig.tight_layout()
    fig.savefig(png, dpi=100)
    print('wrote', png)


if __name__ == '__main__':
    main()
