#!/usr/bin/env python3
"""The evolution of one CZTR_box run (e.g. a cascade level) from its log and
profiles.

Usage: /usr/bin/python3 tools/level_evolution.py OUT_DIR out.png

Top row, time series: max |v_z| and max |v_h| (log.txt); rms v_z at several
heights; the total flux F/F_sun at several heights (radiation + MLT +
resolved enthalpy and kinetic + numerical diffusion + cooling layer).
Bottom row, profiles against z at selected times (colour: time): rms v_z,
the total flux and its parts (the last profile), <s1>.
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
    prof = []
    for f in sorted(glob.glob(f'{d}/profile_[0-9][0-9][0-9][0-9].dat')):
        with open(f) as fh:
            t = float(fh.readline().split('=')[1].split()[0])
        prof.append((t, np.loadtxt(f, comments='#')))
    tp = np.array([p[0] for p in prof])/86400
    z = prof[0][1][:, 2]/1e3                                # [Mm]
    ztop = z[-1]
    ftot = lambda a: a[:, 10] + a[:, 11] + a[:, 12] + a[:, 14] + a[:, 15]

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    # max speeds
    ax[0, 0].semilogy(log[:, 1]/86400, np.maximum(log[:, 4], 1e-3)/1e2, label='max |v_z|')
    ax[0, 0].semilogy(log[:, 1]/86400, np.maximum(log[:, 5], 1e-3)/1e2, label='max |v_h|')
    ax[0, 0].set_ylabel('[m/s]'); ax[0, 0].set_title('maximum speeds')
    # rms v_z and the total flux at heights
    heights = [ztop - 2.0, ztop - 10.0, ztop - 25.0, ztop - 60.0, -120.0, -180.0]
    cm = plt.get_cmap('plasma')
    for q, zz in enumerate(heights):
        i = int(np.argmin(abs(z - zz)))
        c = cm(q/len(heights))
        ax[0, 1].semilogy(tp, [max(p[1][i, 8], 1e-3)/1e2 for p in prof], color=c, label=f'z = {z[i]:.0f} Mm')
        ax[0, 2].plot(tp, [ftot(p[1])[i] for p in prof], color=c, label=f'z = {z[i]:.0f} Mm')
    ax[0, 1].set_ylabel('rms v_z [m/s]'); ax[0, 1].set_title('rms v_z by height')
    ax[0, 2].set_ylabel('F_total / F_sun'); ax[0, 2].set_title('total flux by height')
    ax[0, 2].axhline(1.0, color='k', lw=0.6, ls=':')
    for a in ax[0]:
        a.set_xlabel('t [d]'); a.grid(alpha=0.3); a.legend(fontsize=7)
    # profiles at selected times
    sel = np.unique(np.linspace(0, len(prof) - 1, 9).astype(int))
    cm = plt.get_cmap('viridis')
    for k, s in enumerate(sel):
        t, a = prof[s]
        c = cm(k/max(len(sel) - 1, 1))
        lab = f'{t/86400:.0f} d'
        ax[1, 0].semilogx(np.maximum(a[:, 8], 1e-3)/1e2, z, color=c, label=lab)
        ax[1, 2].plot(a[:, 5], z, color=c, label=lab)
    t, a = prof[-1]
    ax[1, 1].plot(ftot(a), z, 'k', lw=2, label='total')
    for col, lab in ((10, 'radiation'), (11, 'MLT'), (12, 'resolved'), (14, 'diffusion'), (15, 'cooling layer')):
        ax[1, 1].plot(a[:, col], z, label=lab)
    ax[1, 1].set_title(f'flux parts at {t/86400:.1f} d')
    ax[1, 0].set_xlabel('rms v_z [m/s]'); ax[1, 0].set_title('rms v_z profiles')
    ax[1, 1].set_xlabel('F / F_sun')
    ax[1, 2].set_xlabel('<s1> [erg/g/K]'); ax[1, 2].set_title('mean entropy perturbation')
    for a_ in ax[1]:
        a_.set_ylabel('z [Mm]'); a_.grid(alpha=0.3); a_.legend(fontsize=7)
    fig.suptitle(f'{d}: top at z = {ztop:.1f} Mm')
    fig.tight_layout()
    fig.savefig(png, dpi=100)
    print('wrote', png)


if __name__ == '__main__':
    main()
