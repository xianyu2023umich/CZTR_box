# CZTR_box: a local Cartesian box of CZTR_3D

Design notes, 2026-11-22. Implemented the same day: see section 10.

## 1. Purpose

A high-resolution local testbed for CZTR_3D that runs on the office Mac
Studio (32 GB) or the laptop, for questions the global 0.47 deg run raises
but is too expensive to iterate on:

- the step-2 crash (2026-11-22): evacuation and temperature runaway in the
  low corona above a compact flux concentration (lat 12.1, lon 150.5), with
  the B-dependent heating H_B as the prime suspect;
- heating forms (density-scaled field heating, later AWSoM/WTD w±);
- relaxing the seeded field so the horizontal Alfven cap (va_cap_horizontal)
  can be used;
- near-surface and transition-region dynamics at the global run's resolution
  and finer.

The global 0.47 deg run needs ~145 GB; a box of the size below needs a few GB.

## 2. Decisions so far (with the user)

- **A fully Cartesian box, not a spherical wedge** (and not a spherical
  wedge at large radius mimicking one). A wedge's side boundaries are hard
  (theta cannot be periodic in a wedge); a box is periodic horizontally, as
  MURaM and Bifrost, and the Cartesian operators are faster and simple to
  write.
- **The background is the spherical one, unchanged in z.** The static
  background is fixed by the entropy profile, vertical force balance and the
  gravity law; force balance is the same in both geometries, so the
  spherical profile in a box is an exact hydrostatic state with an
  artificial gravity g(z) (the spherical 1/r^2 law). Everything near the
  surface is then directly comparable with the global run. The transport
  (MLT, grey radiation) is computed from the local gradients, so in a box it
  carries a constant flux; in the nearly adiabatic CZ that only changes the
  tiny superadiabatic gradient and the convective velocities (see the flux
  item in section 9), not the background.
- **Top boundary as MURaM's** (Rempel 2017, ApJ 834:10, sec. 2.7): a
  potential-field extrapolation for B, semi-transparent for flows (ghost
  cells 50% and 25% of the top domain values: mass can cross, vertical
  flows strongly damped), zero vertical conductive flux.

## 3. Approach: a fully Cartesian code

CZTR_box is a copy of CZTR_3D whose operators are rewritten for Cartesian
coordinates (the user's choice, 2026-11-22: faster, since the r and
sin(theta) metric work disappears from every loop, and not complicated to
code). The physics (the equations, EOS, MLT/grey transport, the atmosphere's
conduction, losses and heating, the diffusion forms, the radial
semi-implicit scheme) stays the same; only the geometry changes.

- **Domain:** ARBOR's box (arbor_domain_box), periodic in x and y (ARBOR
  supports periodic maps for a box; no ARBOR change needed). No Yin-Yang.
- **z is ARBOR's first direction,** as r is in CZTR_3D: the column
  machinery (radial blocks, si_radial's whole columns per rank, the
  transport columns, the radial copies of ghosts) carries over unchanged
  with "radial" read as "vertical". The z cells are the global grid's (a
  tabulated, non-uniform map, as CZTR_3D's radial map).
- **Gravity** g(z) from the background file (the spherical law: the
  artificial gravity of section 2), constant horizontally.
- **Plane-parallel geometry:** face areas and cell volumes constant in z,
  no curvature terms, no 1/r or cot(theta) factors, derivatives in x and y
  plain d/dx, d/dy.

The price is a fork: physics changes made later in CZTR_3D (or here) have to
be ported by hand. Keep the Cartesian changes confined to the geometry
(operators, metric factors, boundaries) so the physics files can be
compared and synced.

## 4. Grid and size (user, 2026-11-22)

- **Vertical:** the CZTR_3D 560-cell radial grid (the whole CZ, as CZTR_3D),
  trimmed at the radial block top closest to 100 Mm above the surface.
  ARBOR's table map (arbor_map_table) takes the faces as is.
- **Horizontal:** uniform, periodic, tunable (namelist): e.g. 32 x 32 or
  128 x 128 columns, the width up to ~1 R_sun (~700 Mm). 1 R_sun at 32 cells
  is 22 Mm (about 1.875 deg), at 128 cells 5.4 Mm (about 0.47 deg).
- **ARBOR support (checked):** arbor_domain_box (COORD_CARTESIAN metric:
  Area_F = 1, Vol_I = dx, h factors 1), per-direction periodic maps, and
  table maps all exist and are tested (test_m2/m4/m5/m6/m7; the RSST
  convection example is a periodic box). No ARBOR change needed. The
  convection example puts z last; CZTR_box puts it first to reuse the
  column code (just the map order; ARBOR doesn't care).

## 5. Boundaries

- **Sides:** periodic (the ARBOR option above).
- **Bottom:** the existing wall with the luminosity input, but the flux
  given directly (a new option): the code now converts L_bottom through the
  sphere's area 4 pi r^2, which a box doesn't have. Set it to the flux the
  global run carries at that depth (L_sun / 4 pi r_bot^2). Later perhaps a MURaM-like
  open bottom (symmetric mass flux and field, mean pressure fixed,
  perturbations damped).
- **Top:** a new wall_top option 'box' as in section 2: potential field
  from the top layer (FFT in the periodic box), semi-transparent flows,
  zero conductive flux. The global run's 'outflow' (wind) top is not used.

## 6. Initial state (user, 2026-11-22: no export, the box runs its own simulation)

CZTR_box is a general testbed, not only for the December deliverable, so it
starts from its own initial state and evolves its own convection, as
CZTR_3D did:

- **The CZ:** the background plus small random perturbations (CZTR_3D's
  hydro start), run hydro until the convection is statistically steady.
- **The field:** a seed as CZTR_3D's (whatever the current version uses
  for the CZ), then the MHD stage.
- **The atmosphere:** the global init profile, pressure-matched per column
  to the CZ top (atm_match_base), as step 2.
- **The atmosphere's field:** the potential field of the join's B_z,
  computed by FFT in the periodic box (exactly force-free and periodic;
  replaces PFSS/DPF; also a better start for the horizontal Alfven cap).

## 7. Code changes (a copy of CZTR_3D at a29a414)

The places where CZTR_3D writes the spherical geometry out, to rewrite:

1. **Domain and grid** (cz3d_main): arbor_domain_box with periodic x, y and
   the tabulated z map; no Yin-Yang, no refinement patches needed at first.
2. **Background** (cz3d_bg_mod): z instead of r; face areas and cell
   volumes per unit horizontal area (constant); g from the file.
3. **RHS** (rhs, rhs_mhd): the mass flux without r^2 and sin(theta), the
   advection and pressure gradient without 1/r, 1/(r sin theta), the
   curvature terms dropped, no Coriolis.
4. **curl_block, div_block:** Cartesian forms.
5. **Diffusion:** the face loops use ARBOR's metric (box metrics come
   from ARBOR); the resistive face currents (bcell_derivs, bfaces) without
   the h factors; the top-layer blend as is.
6. **Radial semi-implicit solves and the B radial implicit term:** the
   r_f^2 face areas and r^2 volumes become constants.
7. **Transport columns** (MLT, grey, grey_atm): constant areas and volumes;
   the bottom flux given directly (section 5).
8. **Atmosphere** (cztr3d_atm_mod): the conduction's face distances and
   gradients without r and sin(theta); the losses and heating unchanged.
9. **Time step, Alfven cap, GLM:** dx, dy instead of r dtheta, r sin(theta) dphi.
10. **Boundaries:** the top 'box' boundary (section 5); the sides periodic
    (ARBOR); the bottom wall as now with the flux option.
11. **Field seed:** the FFT potential field of the join's B_r (replaces
    PFSS/DPF).
12. **Initial state:** the box's own start (section 6); no import tool.
13. **Global pieces off or adapted:** Yin-Yang overlap weights in the
    diagnostics (box: uniform weights), the angular-momentum correction,
    rotation, the synthetic images (a view along z, or slanted), and the
    normalizations to the full sphere (L/L_sun in the log, the floors
    report's fluxes per box area).

## 8. First uses

1. Step 2's settings (H_B on) at 0.47 deg-like resolution: does the
   evacuation/runaway above strong concentrations appear? Compare with
   dbgA's box output when it arrives.
2. The same with H_B = 0 and with a density-scaled field heating.
3. The horizontal Alfven cap from a force-free start.
4. Resolution: 0.47 deg against 2-4x finer.

## 9. Decided / open

Decided (user, 2026-11-22): box size and resolution tunable (e.g. 32^2 or
128^2, up to ~1 R_sun wide); the whole CZ; the CZTR_3D 560 grid trimmed at
~100 Mm; the box runs its own simulation (no export); the physics is
CZTR_3D's current version.

Open:
- **The energy flux with the whole CZ.** In the sphere the flux falls as
  1/r^2: at the CZ base (0.71 R) it is ~2x the surface flux. In a box with
  constant areas the flux is the same at every depth. Proposal: give the
  bottom the surface flux (L_sun / 4 pi R_sun^2), so the photosphere and
  atmosphere get the solar flux; the deep CZ then carries ~half its
  spherical flux, so its convective velocities are ~(1/2)^(1/3) ~ 0.8 of the
  global run's and its superadiabaticity lower. The background
  (spherical, hydrostatic) is unaffected.
- **Width vs curvature:** 1 R_sun flat is fine for a testbed; only the
  comparison with the global run's patches is approximate.

## 10. Implementation status (2026-11-22)

Implemented as a copy of CZTR_3D a29a414 (`cztr_box.x`; `make`, the same
`config.mk` as CZTR_3D). Variable slots keep CZTR_3D's names: v_r, v_theta,
v_phi hold v_z, v_x, v_y; B_r, B_theta, B_phi hold B_z, B_x, B_y.

- Domain: `arbor_domain_box`, z first (the 1D grid's r as the vertical
  coordinate), x and y uniform and periodic. Namelist: `box_lx_mm`,
  `box_ly_mm`, `n_x`, `n_y`, `nb_x`, `nb_y` (and `nb_r`, `refine_h`, the patch
  as `patch_x_mm`, `patch_y_mm`, `patch_radius_mm`).
- Operators: RHS, curl, div, the time step, the Alfven cap, GLM, the column
  solves, the radial B implicit term, the resistive face currents and the
  atmosphere's conduction in Cartesian form; ARBOR's box metric in the flux
  divergences.
- Energy flux: the bottom gets F_bottom = L_bottom_Lsun L_sun/(4 pi R^2) per
  unit box area. The 1D radial transport (grey, MLT, the cooling layer) uses
  the sphere's relative face areas (r_f/R)^2 (`area_f`, volumes dz): its
  physical fluxes keep the spherical L(r), so the spherical background stays
  in balance, and the box flux area_f F is F_sun at every depth. (Without
  this the background's radiative flux at the CZ base, 1.9 F_sun, would
  exceed the bottom flux.) The resolved convection is not scaled: in the deep
  CZ it carries (r/R)^2 of its spherical flux (velocities ~0.8x at the base).
  Logs and profiles report fluxes as F/F_sun (columns still named L/Lsun).
- Top: `wall_top = 'box'`: v_z ghosts 0.5, 0.25 x the top cell's
  (`box_top_vz`), v_h and s1 copied, rho/rho0 copied; B by `b_bc_top`:
  'potential' (MURaM: the ghosts the potential field of the top layer's
  B_z, a DFT once per step, collective; use this), 'open' (zero gradient) or
  'radial'/'conductor'. 'radial' with a potential field puts a current sheet
  at the top (B_h != 0 there): v -> 700 km/s in 20 s in the first test.
  b_seed_test also checks the potential ghosts on the first step: exact to
  round-off (1e-14 G). No conductive flux through the top.
- Seeds: `b_seed` 'uniform' (B_z), 'horizontal' (B_x w(z), the toroidal
  analogue), 'potential' (FFT potential field of the join's B_z above the
  join); `b_seed_test` checks 'potential' against the analytic mode-(1,1)
  field: exact to output precision (5e-7 of 10 G) at all heights.
- Diagnostics: area weights, x/y positions [Mm] (footpoints, max Mach, NaN,
  probe), `region_*` output (CZTR_3D's box output), layer dumps with x, y.
- Removed: Yin-Yang, PFSS/DPF/magnetofriction, rotation/Coriolis/L_z,
  the synthetic EUV and continuum images (a vertical view is easy to add
  later).
- Fixed on the way (also in CZTR_3D): the atmosphere's ARBOR timers were
  created lazily in atm_block, so a rank without atmosphere blocks hung in
  the collective timer report at the end; si_radial's level check read the
  tree's whole key array instead of its leaves.

Tests (laptop, 32 x 32 columns over 696 Mm):
- `runs_box/smoke_hd.nml` (CZ only, hydro, 100 steps, 8 ranks): clean,
  reproducible checksum; box flux 1.0 F_sun at the bottom and top.
- `runs_box/smoke_atm.nml` (336 cells to 73.7 Mm, atmosphere, MHD, the
  potential test seed, step 2's heating, si_radial, 40 steps): clean; max
  Mach 0.4 at the transition region, the top quiet. (With si_radial whole
  columns go to ranks: 4 columns at 32^2 / 16^2 blocks; use more columns
  than ranks.)

Timing (2026-11-22, laptop, 8 ranks, runs_box/timing128.nml: 128 x 128
columns over 696 Mm (5.4 Mm), 336 cells, 768 blocks of 28 x 16 x 16, 5.5M
cells, step 2's physics, the potential top, 20 steps from the start):
0.645 s/step (0.94 us per cell-step per rank), ~0.9 GB per
rank, dt 0.435 s (the CFL of the top cells, which si_radial leaves explicit,
and the wave cap at 2 Mm, 0.449 s). Per stage (share of the wall time): rhs
21% (EOS 6, curls 5, Lorentz/psi 4, si_coeff 3, advection 3), diffusion 17%
(faces B 4, EOS 4, radial B implicit 3), radial implicit solves 12% (the
column solves 11%: mostly the transposition), ghost fills 12%, the interior
EOS during the fills 8%, energy transport 8%, atmosphere 5% + losses during
the fill 5%, stable_dt 4%, the potential top 2%, outputs 2.5%.

Optimization (2026-11-22, bit-identical: the -ffp-contract=off builds,
tools/build_nofma.sh, give the same checksums as before for smoke_atm and
smoke_hd; laptop timings vary +-10% run to run):
- The column solves in row groups (fixed k, nb(2) columns) instead of all
  nb(2) nb(3) columns of a stack at once: the work arrays stay in cache.
  128^2: 1.41 -> 0.64 s per 20 steps (16^2 blocks), 2.40 -> 0.64 s (32^2
  blocks); the step 12.9 -> ~11.1 s per 20 steps. (Generic: also for CZTR_3D.)
- The RHS derivatives fused (inline stencils dz4, dh4 with the block routines'
  arithmetic; ivdep and a branch-free pressure loop so they vectorize):
  advection 0.39 -> 0.30, pressure 0.17 -> 0.13, curls/div 0.60 -> 0.55,
  Lorentz/psi 0.48 -> 0.43 s per 20 steps; mass kept with rho v stored once
  (fused it recomputes the products, 0.17 -> 0.22). About 2% of the step: the
  RHS was not memory bound at this block size.
- 32^2 blocks: the fills, the potential top and the atmosphere ~10% cheaper,
  the step about the same as 16^2 within the noise.

Radial cascade (2026-11-22; the plan of CZ_3D/notes/radial_cascade.md):
CZ-only runs with the top at a radial-block face (n_r_cells 112, 140, 168,
196 = 0.954, 0.991, 0.998, 1.0 R), each exported and imported into the next;
the blocks above the export's top start from the init file (imp_top). Deep
tops: face wall, cooling layer (cool_width_hp 1, F_sun per box area), MLT
off; one RSST through the cascade (the user, 2026-11-22): anchored 500 km
deep (rsst_depth_km, the anchor's c_s ~9 km/s from the whole profile's CZ,
z_cz/cs_cz in cz3d_bg_mod, so also above a deep top), xi 25 at the base;
with no cell within grey_top_km the radiation is the explicit diffusion
everywhere. runs_box/casc1..4.nml. Chain 1 -> 2 -> 3 -> 4 tested by export
and import (20 steps each, 128^2, laptop, 8 ranks, ~0.40 us per cell-step per
rank): dt 48.7, 16.0, 4.9, 2.1 s (each the top cell's CFL with c ~9 km/s;
level 1 at the base cell); wall 0.09, 0.115, 0.135, 0.17 s per step. With the
plan's lengths (1 yr, 2 weeks, 3 days, 2 days) ~16 + 2.4 + 2 + 4 h on the
laptop at 128^2. (A fixed c_target 2 km/s gave 220, 72, 22 s for levels 1-3
but cannot be used at the photosphere.) Small flux transients where a new
layer meets the old top. Open: starting a new layer from a relaxed run's
export instead of the 1D init profile; the old top's cooled layer adjusts
after each insertion.

Potential top at deep walls and f-plane rotation (2026-11-22):
- b_bc_top = 'potential' with any top wall (the deep face walls of the
  cascade too). A cascade import with the field fills the new layer with the
  potential extension of the export's top B_z (cz_extend_potential), the
  field the potential top assumed: B continuous at the old top (tested:
  rms(div B dx)/rms|B| ~1e-9 before and after the insertion). (Chosen over a
  conductor top, which traps horizontal flux under the wall and blocks the
  magnetic helicity flux, and over the radial top, which forces B_h = 0.)
- f-plane: omega_rot, fplane_lat_deg (default 30), coriolis 'split'/'rk';
  Omega = omega_rot (sin(lat) e_z + cos(lat) e_y) constant over the box (y
  north; the full Coriolis term, horizontal component included). An eastward
  flow is deflected south (to its right) and upward. Without rotation the box
  has only a small-scale dynamo; with it a large-scale one is possible.
- runs_box/casc1..4: the potential top, the f-plane at 30 deg, 414 nHz.

Next: a hydro CZ run to a convective steady state
(the deep CZ's thermal adjustment to the box flux), then the atmosphere and
field; a vertical-view synthetic image.

Unresolved near-surface turbulence and the field (2026-11-23):
- Why: the 0.47 deg CZTR_3D step-2 crashes all sit at grid-scale kG surface
  field (mostly horizontal, ~1% of the area) that the Alfven cap lets pile up
  at the join. MLT carries the heat of the unresolved granulation but does
  nothing to B: no turbulent diffusion or pumping keeps the surface field weak
  and near radial, as surface flux transport needs (250-600 km^2/s; downward
  pumping in flux-transport dynamos, Karak & Cameron 2016).
- Turbulent closure, local per cell from the cell's MLT state (the same
  column evaluation as the heat flux; f_mlt v_MLT, v_MLT of the full MLT flux
  F/f, ell = alpha_mlt H_p): turb_c_eta (eta_t = c f v ell/3 added to the
  resistive eta: radial implicit, horizontal explicit, Ohmic heat),
  turb_c_nu (nu_t on v: horizontal explicit, radial implicit with the
  semi-implicit column solver, the kinetic energy lost into s1),
  turb_c_gamma (downward pumping gamma_z = -c f v in the induction only,
  E = (v + gamma_z e_z) x B, 0 beyond the walls). All default 0; need
  use_convection, si_radial, global steps (eta: b_diff_form = 'resistive').
  out_dir/turb_NNNN.dat: the horizontal means. Local, not a profile, so that
  heat flux, eta, nu, gamma come from one state (an LES-style subgrid
  closure) and the magnetic suppression acts where the field is.
- mlt_gough_tayler: MLT sees grad - grad_ad - delta, delta = B_z^2/(B_z^2 +
  4 pi Gamma_1 P) at each face (Gough & Tayler 1966; magnetic stellar models).
  Vertical field suppresses convection: umbral darkening; the closure
  follows through v_MLT. Kept (the user).
- turb_quench (option, default off): the coefficients times q = 1/(1 +
  B_h^2/B_eq^2), B_eq^2 = 4 pi rho v_MLT^2. A heuristic for the horizontal
  field (mean-field quenching is usually isotropic in |B|; Gough & Tayler
  covers B_z, so B_h only to avoid double counting). Near the surface B_eq
  ~ 90 G: strong horizontal flux is not the turbulence's to remove but
  rises by its own buoyancy, which needs the real Lorentz force (the cap).
- Tests (laptop, 8 ranks): defaults bit-identical (no-FMA smoke_hd 16225155724550,
  smoke_atm 24447136486753); level 4 at 32^2: <f v_MLT> 90 m/s at -28 Mm to
  1.8 km/s at -400 km, 0 in the top ~300 km (no MLT flux there), eta_t ~
  1.3-3.7e12 cm^2/s with c = 1; signs on a static seed: pumping moves B_h
  down, eta_t spreads it; nu_t stable on a convecting level-1 restart (+13%
  cost at 64^2); Gough & Tayler with a uniform 300 G B_z: MLT flux 1.00 ->
  0.80 at -2 Mm, 0.29 -> 0.04 at -0.4 Mm, top luminosity 0.99 -> 0.92 within
  6 min (heat collecting below). Open: c_gamma (1 gives 1-2 km/s near the
  top; flux transport suggests ~0.02-0.1); a convecting level-4 run with
  field, closure on/off.
