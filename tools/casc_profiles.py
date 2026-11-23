#!/usr/bin/env python3
"""Profiles across a radial-cascade insertion (CZTR_box).

Usage: /usr/bin/python3 tools/casc_profiles.py LEVEL_DIR [NEXT_LEVEL_DIR] [out.png]

Reads out_dir/profile_NNNN.dat of a level (and of the next level, started
from its export) and shows, against z near the old top:
  - the total flux F/F_sun (radiation + MLT + resolved enthalpy and kinetic
    + numerical diffusion + cooling layer),
  - rms(v_z), rms(v_h),
  - <s1> and rms_h(s1),
for the last profiles of the level and the profiles of the next level
(colour: time), and time series at heights around the old top. Prints a
table of the same.
"""
import glob
import sys

import numpy as np


def load(d):
    out = []
    for f in sorted(glob.glob(f'{d}/profile_[0-9][0-9][0-9][0-9].dat')):
        with open(f) as fh:
            t = float(fh.readline().split('=')[1].split()[0])
        a = np.loadtxt(f, comments='#')
        out.append((t, a))
    return out


def ftot(a):
    return a[:, 10] + a[:, 11] + a[:, 12] + a[:, 14] + a[:, 15]


def main():
    d1 = sys.argv[1]
    d2 = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].endswith('.png') else None
    png = sys.argv[-1] if sys.argv[-1].endswith('.png') else None
    p1 = load(d1)
    p2 = load(d2) if d2 else []
    z1top = p1[0][1][-1, 2]                       # the level's top cell [km]
    print(f'{d1}: {len(p1)} profiles, top cell z = {z1top:.0f} km' + (f'; {d2}: {len(p2)} profiles' if d2 else ''))
    # heights for the time series: below the old top, at it, in the new layer
    zs = [z1top - 12000.0, z1top - 4000.0, z1top - 1000.0]
    if p2:
        ztop2 = p2[0][1][-1, 2]
        zs += [z1top + 2000.0, 0.5*(z1top + ztop2), ztop2 - 500.0]
    rows = []
    for lev, prof in ((1, p1), (2, p2)):
        for t, a in prof:
            z = a[:, 2]
            vals = []
            for zz in zs:
                if zz > z[-1] + 1.0:
                    vals.append((np.nan, np.nan, np.nan))
                    continue
                i = int(np.argmin(abs(z - zz)))
                vals.append((ftot(a)[i], a[i, 8], a[i, 5]))
            rows.append((lev, t, vals))
    print('lev  t[d]   ' + '  '.join(f'z={zz/1e3:6.1f}Mm: F  rms_vz[m/s]  <s1>' for zz in zs))
    for lev, t, vals in rows:
        print(f'{lev:3d} {t/86400:6.2f}  ' + '  '.join(f'{v[0]:7.3f} {v[1]/100:8.1f} {v[2]:9.2e}' for v in vals))
    if not png:
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    zlo = z1top - 40000.0
    zhi = (p2[0][1][-1, 2] if p2 else z1top) + 500.0
    sets = [(p1[-3:], 'Greys', 'level 1')] + ([(p2, 'viridis', 'level 2')] if p2 else [])
    for prof, cmap, lab in sets:
        cm = plt.get_cmap(cmap)
        for k, (t, a) in enumerate(prof):
            z = a[:, 2]
            m = (z > zlo) & (z < zhi)
            c = cm(0.35 + 0.65*k/max(len(prof) - 1, 1))
            kw = dict(color=c, lw=1, label=f'{lab} {t/86400:.2f} d' if k in (0, len(prof) - 1) else None)
            ax[0, 0].plot(z[m]/1e3, ftot(a)[m], **kw)
            ax[0, 1].plot(z[m]/1e3, a[m, 8]/1e2, **kw)
            ax[0, 2].plot(z[m]/1e3, a[m, 9]/1e2, **kw)
            ax[1, 0].plot(z[m]/1e3, a[m, 5], **kw)
            ax[1, 1].plot(z[m]/1e3, a[m, 6], **kw)
    for a_, yl in zip([ax[0, 0], ax[0, 1], ax[0, 2], ax[1, 0], ax[1, 1]],
                      ['F_total / F_sun', 'rms v_z [m/s]', 'rms v_h [m/s]', '<s1> [erg/g/K]', 'rms_h(s1)']):
        a_.axvline(z1top/1e3, color='r', ls=':', lw=0.8)
        a_.set_xlabel('z [Mm]')
        a_.set_ylabel(yl)
        a_.grid(alpha=0.3)
    ax[0, 0].legend(fontsize=7)
    # time series
    tt = np.array([r[1] for r in rows])/86400
    for q, zz in enumerate(zs):
        ax[1, 2].plot(tt, [r[2][q][1]/1e2 for r in rows], label=f'z = {zz/1e3:.1f} Mm')
    ax[1, 2].set_xlabel('t [d]')
    ax[1, 2].set_ylabel('rms v_z [m/s]')
    ax[1, 2].legend(fontsize=7)
    ax[1, 2].grid(alpha=0.3)
    fig.suptitle(f'cascade insertion: {d1}' + (f' -> {d2}' if d2 else '') + ' (red: the old top cell)')
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    print('wrote', png)


if __name__ == '__main__':
    main()
