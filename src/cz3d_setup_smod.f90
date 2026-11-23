submodule (cz3d_app_mod) cz3d_setup_smod
  ! setup: initialisation, field seeds (uniform, horizontal, the box's potential field, smoothing), state import/export and matching, the column distribution, restart files
  ! (part of cz3d_app_mod: the module's data and routines by host association)
  implicit none

contains

  ! The number of state variables (5 hydro; 9 with the magnetic field)
  module integer function cz_n_var()
    cz_n_var = merge(nVarB, nVarCz, use_B)
  end function cz_n_var

  ! Refinement to refine_h levels in x and y everywhere and to patch_level
  ! near the patch centre (periodic distance; arbor_regrid)
  module subroutine cz_flag_block(self, iBlock, flags)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlock
    integer, intent(out) :: flags(3)
    real(dp) :: lx, ly, dx, dy
    integer :: j, k
    flags = 0
    if (Key_IB(KL+2, iBlock) < refine_h) flags(2:3) = 1
    if (Key_IB(KL+2, iBlock) >= patch_level) return
    lx = box_lx_mm*1.0e8_dp;  ly = box_ly_mm*1.0e8_dp
    do k = 1, nb(3)
      do j = 1, nb(2)
        dx = Metric_DB(2, iBlock)%X_I(j) - patch_x_mm*1.0e8_dp
        dy = Metric_DB(3, iBlock)%X_I(k) - patch_y_mm*1.0e8_dp
        dx = dx - lx*anint(dx/lx);  dy = dy - ly*anint(dy/ly)
        if (dx**2 + dy**2 <= (patch_radius_mm*1.0e8_dp)**2) then
          flags(2:3) = 1
          return
        end if
      end do
    end do
  end subroutine cz_flag_block

  ! After arbor_block_sync or a regrid (blocks may have changed): the
  ! per-block arrays (zeroed: the convective flux restarts at its targets),
  ! the subgrid factor, the scratch arrays (once)
  module subroutine cz_setup()
    integer :: i, lev
    real(dp) :: Hp, R, dxl
    if (mod(n_r, nb(1)) /= 0) call arbor_abort('cz3d: the radial block size must divide the radial cells')
    if (nb(1) < n_r .and. dt_budget > 0.0_dp) call arbor_abort('cz3d: the energy budget needs whole columns per block')
    if (mr_levels > 1 .and. (dt1 > 0.0_dp .or. dt_budget > 0.0_dp)) &
        call arbor_abort('cz3d: multirate needs dt1 = 0 and no energy budget')
    call cz_unbind()
    rsst_h2012 = rsst_form == 'hotta2012'
    glm_cfl = glm_ch == 'cfl'
    grey_i0_g = 1                                 ! first cell of the implicit grey solve (whole column)
    if (grey_top_km > 0.0_dp) then
      ! (a deep top, the cascade: no cell within grey_top_km of the photosphere;
      ! explicit radiative diffusion in the whole column)
      grey_i0_g = n_r + 1
      do i = 1, n_r
        if (-z_c(i) < grey_top_km*1.0e5_dp) then
          grey_i0_g = i
          exit
        end if
      end do
    end if
    if (use_atmosphere .and. mod(k_join, nb(1)) /= 0) call arbor_abort('cztr3d: the join cell must end a radial block '// &
        '(k_join a multiple of the radial block size)')
    if (nb(1) < merge(k_join, n_r, use_atmosphere) .and. grey_i0_g <= merge(k_join, n_r, use_atmosphere) - nb(1)) &
        call arbor_abort('cz3d: with radial blocks the implicit grey layer '// &
        '(grey_top_km) must lie within the top block')
    grey_i0 = grey_i0_g
    rsst_h2015 = rsst_form == 'hotta2015'
    if (use_atmosphere) call atm_setup()
    if (allocated(frz_B)) deallocate(frz_B, frz_have)
    if (freeze_side /= 'none') then
      allocate(frz_B(nb(1), nb(2), nb(3), cz_n_var(), nBlockMax), frz_have(nBlockMax))
      frz_have = .false.
    end if
    if (allocated(qrad_FB)) deallocate(qrad_FB, qdif_FB)
    if (allocated(qconv_FB)) deallocate(qconv_FB)
    if (associated(f_mlt_g)) deallocate(f_mlt_g)
    if (.not. qconv_face) then
      allocate(qconv_FB(0:nb(1), nb(2), nb(3), nBlockMax))
      qconv_FB = 0.0_dp
    end if
    allocate(qrad_FB(0:nb(1), nb(2), nb(3), nBlockMax))
    allocate(qdif_FB(0:nb(1), nb(2), nb(3), nBlockMax))
    qrad_FB = 0.0_dp
    qdif_FB = 0.0_dp
    ! the bottom face's flux: the solar surface flux (L_bottom_Lsun L_sun over
    ! the photosphere's 4 pi R^2), the same at every depth in the box
    F_bottom = L_bottom_Lsun*Lsun/(4.0_dp*pi*r_top_cz**2)
    ! the subgrid factor per horizontal level (the cell size halves per level)
    allocate(f_mlt_g(0:n_r, 0:max(refine_h, patch_level)))
    f_mlt_g = 1.0_dp
    do lev = 0, ubound(f_mlt_g, 2)
      dxl = min(box_lx_mm/real(n_x, dp), box_ly_mm/real(n_y, dp))*1.0e8_dp/2.0_dp**lev
      if (.not. mlt_switch) cycle
      do i = 0, n_r
        ! background pressure scale height at the face per horizontal cell size
        Hp = 0.5_dp*(p0(i)/(rho0(i)*g_c(i)) + p0(i+1)/(rho0(i+1)*g_c(i+1)))
        R = Hp/dxl
        f_mlt_g(i, lev) = 1.0_dp/(1.0_dp + (R/mlt_R0)**mlt_p)
      end do
    end do
    lev_view = refine_h
    f_mlt(0:n_r) => f_mlt_g(0:n_r, refine_h)
    if (.not. allocated(tloc)) then
      allocate(tloc(1:n_r))
      tloc = 0.0_dp
    end if
    if (si_radial .and. si_cols) call si_setup()     ! (after cz_columns_to_ranks)
    if (allocated(rho_G)) return
    allocate(tsc_g(1:n_r), rate_r_g(1:n_r))
    tsc => tsc_g
    rate_r => rate_r_g
    tsc = 1.0_dp
    rate_r = 0.0_dp
    allocate(rq_g(-1:n_r+2))
    rq => rq_g
    rq = 1.0_dp
    if (rsst_h2015) rq = 1.0_dp/xi**2
    call background_h_phi()
    allocate(rho0_w(-1:n_r+2), rq_w(-1:n_r+2), s0_w(-1:n_r+2), r_w(-1:n_r+2))
    rho0_w = rho0(-1:n_r+2)
    s0_w = s0(-1:n_r+2)
    r_w = r_c(-1:n_r+2)
    rq_w = rq(-1:n_r+2)
    allocate(rho_G(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng))
    allocate(p1_G, cs_G, W_G, lim_G, T_G, hP_G, mold=rho_G)
    nv = cz_n_var()
    ie_ = nv + 1
    allocate(D_C(nb(1), nb(2), nb(3)), R_CV(nb(1), nb(2), nb(3), nv))
    allocate(Dr_C, mold=D_C)
    ! (variable ie_ = nv + 1: the diffusion's energy flux)
    allocate(F1_FV(0:nb(1), nb(2), nb(3), ie_), F2_FV(nb(1), 0:nb(2), nb(3), ie_), &
        F3_FV(nb(1), nb(2), 0:nb(3), ie_), Div_CV(nb(1), nb(2), nb(3), ie_))
    allocate(C1_F(0:nb(1), nb(2), nb(3)), C2_F(nb(1), 0:nb(2), nb(3)), C3_F(nb(1), nb(2), 0:nb(3)))
    allocate(B_GV(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, 3))       ! (also the mass flux's scratch: rhs)
    if (use_B) allocate(C_CV(nb(1), nb(2), nb(3), 3), &
        J_CV(nb(1), nb(2), nb(3), 3), divB_C(nb(1), nb(2), nb(3)), vA_G(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng))
  end subroutine cz_setup

  ! Initial state: the 1D state in every column, plus an optional random s1
  ! (from each cell's global index, so independent of the number of ranks)
  module subroutine cz_init_block(self, iBlock)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlock
    integer :: i, j, k
    call cz_bind(iBlock)
    State_GVB(:, :, :, :, iBlock) = 0.0_dp
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          State_GVB(i, j, k, rho1_, iBlock) = rho1_init(i)/rq(i)
          State_GVB(i, j, k, vr_, iBlock) = v_init(i)
          State_GVB(i, j, k, s1_, iBlock) = s1_init(i)
          if (ds1_random > 0.0_dp .and. (.not. use_atmosphere .or. r_off + i <= k_join)) &
              State_GVB(i, j, k, s1_, iBlock) = s1_init(i) + ds1_random &
              *urand(Key_IB(KP, iBlock), r_off + i, Key_IB(KI+2, iBlock)*nb(2) + j, Key_IB(KI+3, iBlock)*nb(3) + k)
        end do
      end do
    end do
    if (len_trim(import_dir) > 0) call import_block(iBlock)
    if (use_B) call seed_field(iBlock)
    call cz_unbind()
  end subroutine cz_init_block

  ! The seed field added to the block's field: b_seed 'uniform' (B_z =
  ! b_seed_G) or 'horizontal' (B_x = b_seed_G w(z), w = 1 below b_seed_z1_km,
  ! 0 above b_seed_z2_km, cos^2 between); both divergence free
  subroutine seed_field(iBlock)
    integer, intent(in) :: iBlock
    real(dp) :: w, z1, z2
    integer :: i, j, k
    select case (b_seed)
    case ('potential')
      if (.not. b_seed_test) return
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            if (r_off + i > k_join) cycle
            State_GVB(i, j, k, br_, iBlock) = b_seed_G*cos(2.0_dp*pi*Metric_DB(2, iBlock)%X_I(j)/(box_lx_mm*1.0e8_dp)) &
                *cos(2.0_dp*pi*Metric_DB(3, iBlock)%X_I(k)/(box_ly_mm*1.0e8_dp))
          end do
        end do
      end do
    case ('uniform')
      State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_, iBlock) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_, iBlock) + b_seed_G
    case ('horizontal')
      z1 = b_seed_z1_km*1.0e5_dp
      z2 = b_seed_z2_km*1.0e5_dp
      do i = 1, nb(1)
        if (z_c(i) <= z1) then
          w = 1.0_dp
        else if (z_c(i) >= z2) then
          w = 0.0_dp
        else
          w = cos(0.5_dp*pi*(z_c(i) - z1)/(z2 - z1))**2
        end if
        State_GVB(i, 1:nb(2), 1:nb(3), bt_, iBlock) = State_GVB(i, 1:nb(2), 1:nb(3), bt_, iBlock) + b_seed_G*w
      end do
    end select
  end subroutine seed_field

  ! A hydrodynamic restart read into an MHD run (arbor_restart_read with
  ! n_var_file < nv): B and psi came back zero; add the seed field to every
  ! block (init_block is not called on a restart)
  module subroutine cz_field_from_hydro()
    integer :: jb
    character(len=200) :: msg
    if (.not. use_B) return
    do jb = 1, nBlock
      call cz_bind(iBlockList_I(jb))
      State_GVB(:, :, :, br_:psi_, iBlockList_I(jb)) = 0.0_dp
      call seed_field(iBlockList_I(jb))
    end do
    call cz_unbind()
    write(msg, '(A,A,A,F8.2,A)') '  magnetic field added to a hydrodynamic restart: seed ', trim(b_seed), ' (', b_seed_G, ' G)'
    call arbor_print(msg)
  end subroutine cz_field_from_hydro

  ! The box's potential fields: the DFT tables of the (n_x, n_y) map at the
  ! level refine_h (cell centres x_i = (i - 1/2) lx/nx; modes -n/2 .. n/2 - 1)
  subroutine pf_init()
    integer :: nx, ny, m, i
    real(dp) :: lx, ly
    if (allocated(pf_ex)) return
    nx = n_x*2**refine_h;  ny = n_y*2**refine_h
    lx = box_lx_mm*1.0e8_dp;  ly = box_ly_mm*1.0e8_dp
    allocate(pf_kx(nx), pf_ky(ny), pf_ex(nx, nx), pf_ey(ny, ny))
    do m = 1, nx
      pf_kx(m) = 2.0_dp*pi*real(m - 1 - nx*merge(1, 0, m - 1 >= nx/2), dp)/lx
      do i = 1, nx
        pf_ex(m, i) = exp(cmplx(0.0_dp, pf_kx(m)*(real(i, dp) - 0.5_dp)*lx/nx, dp))
      end do
    end do
    do m = 1, ny
      pf_ky(m) = 2.0_dp*pi*real(m - 1 - ny*merge(1, 0, m - 1 >= ny/2), dp)/ly
      do i = 1, ny
        pf_ey(m, i) = exp(cmplx(0.0_dp, pf_ky(m)*(real(i, dp) - 0.5_dp)*ly/ny, dp))
      end do
    end do
  end subroutine pf_init

  ! The Fourier coefficients bk (nx, ny) of B_z in the global cell layer ig
  ! (collective: the layer's map summed over the ranks, then the DFT as two
  ! matrix products, the same on every rank)
  subroutine pf_coeffs(ig, bk)
    integer, intent(in) :: ig
    complex(dp), intent(out) :: bk(:,:)
    real(dp), allocatable :: bz(:,:)
    integer :: jb, ib, i, j, k
    call pf_init()
    allocate(bz(size(pf_ex, 1), size(pf_ey, 1)))
    bz = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (Key_IB(KL+2, ib) /= refine_h .or. Key_IB(KL+3, ib) /= refine_h) &
          call arbor_abort('cztr_box: the potential fields need every block at the horizontal level refine_h')
      i = ig - Key_IB(KI+1, ib)*nb(1)
      if (i < 1 .or. i > nb(1)) cycle
      do k = 1, nb(3)
        do j = 1, nb(2)
          bz(Key_IB(KI+2, ib)*nb(2) + j, Key_IB(KI+3, ib)*nb(3) + k) = State_GVB(i, j, k, br_, ib)
        end do
      end do
    end do
    call arbor_sum_inplace(bz)
    bk = matmul(matmul(conjg(pf_ex), cmplx(bz, 0.0_dp, dp)), transpose(conjg(pf_ey)))/real(size(bz), dp)
  end subroutine pf_coeffs

  ! The potential field of the coefficients bk at height dz above their
  ! layer, on the cells of block ib: B(j, k, 1:3) = (B_z, B_x, B_y)
  subroutine pf_eval(ib, bk, dz, B)
    integer, intent(in) :: ib
    complex(dp), intent(in) :: bk(:,:)
    real(dp), intent(in) :: dz
    real(dp), intent(out) :: B(:,:,:)
    complex(dp), allocatable :: ak(:,:,:), px(:,:)
    complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
    real(dp) :: kk, dec
    integer :: mx, my, c, j0, k0
    allocate(ak(size(bk, 1), size(bk, 2), 3), px(nb(2), size(bk, 2)))
    do my = 1, size(bk, 2)
      do mx = 1, size(bk, 1)
        kk = sqrt(pf_kx(mx)**2 + pf_ky(my)**2)
        dec = exp(-kk*dz)
        if (kk == 0.0_dp) then
          ak(mx, my, :) = [bk(mx, my), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
        else if (dec < 1.0e-14_dp) then
          ak(mx, my, :) = (0.0_dp, 0.0_dp)
        else
          ak(mx, my, 1) = bk(mx, my)*dec
          ak(mx, my, 2) = -ii*pf_kx(mx)/kk*bk(mx, my)*dec
          ak(mx, my, 3) = -ii*pf_ky(my)/kk*bk(mx, my)*dec
        end if
      end do
    end do
    ! sum over mx at the block's x, then over my at its y
    j0 = Key_IB(KI+2, ib)*nb(2);  k0 = Key_IB(KI+3, ib)*nb(3)
    do c = 1, 3
      px = matmul(transpose(pf_ex(:, j0+1:j0+nb(2))), ak(:, :, c))
      B(:, :, c) = real(matmul(px, pf_ey(:, k0+1:k0+nb(3))))
    end do
  end subroutine pf_eval

  ! b_seed = 'potential': the field above the join (cells > k_join) replaced
  ! by the potential field of the join cell's B_z in the periodic box
  ! (cz_extend_potential)
  module subroutine cz_seed_potential()
    call cz_extend_potential(k_join)
  end subroutine cz_seed_potential

  ! The field above global cell ig0 (cells > ig0) replaced by the potential
  ! field of cell ig0's B_z (pf_coeffs, pf_eval: each mode decays as
  ! e^{-|k| dz} above the cell's centre; the mean b_0 continues as a uniform
  ! B_z); psi = 0 there. The seed (ig0 = k_join) and the radial cascade (ig0
  ! the export's top: the new layer's field is the potential extension the
  ! potential top assumed, B_z continuous)
  module subroutine cz_extend_potential(ig0)
    integer, intent(in) :: ig0
    complex(dp), allocatable :: bk(:,:)
    real(dp), allocatable :: Bc(:,:,:)
    real(dp) :: z0, t0, bmax
    integer :: jb, ib, i
    character(len=200) :: msg
    t0 = arbor_wtime()
    call pf_init()
    allocate(bk(size(pf_ex, 1), size(pf_ey, 1)), Bc(nb(2), nb(3), 3))
    call pf_coeffs(ig0, bk)
    call cz_unbind()
    z0 = z_c(ig0)
    bmax = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if ((Key_IB(KI+1, ib) + 1)*nb(1) <= ig0) cycle
      call cz_bind(ib)
      do i = 1, nb(1)
        if (r_off + i <= ig0) cycle
        call pf_eval(ib, bk, z_c(i) - z0, Bc)
        State_GVB(i, 1:nb(2), 1:nb(3), br_:bp_, ib) = Bc
        State_GVB(i, :, :, psi_, ib) = 0.0_dp
        bmax = max(bmax, maxval(abs(Bc)))
      end do
    end do
    call cz_unbind()
    bmax = arbor_max(bmax)
    write(msg, '(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A,F8.2,A)') '  potential field above cell ', ig0, ': ', size(bk, 1), &
        ' x ', size(bk, 2), ' modes, mean B_z ', real(bk(1, 1)), ' G, max |B| above ', bmax, ' G; ', &
        arbor_wtime() - t0, ' s'
    call arbor_print(msg)
  end subroutine cz_extend_potential

  ! b_bc_top = 'potential': the top ghosts' B (btop_GB) from the potential
  ! field of the top cell layer's B_z (collective; once per step)
  module subroutine cz_top_potential()
    complex(dp), allocatable :: bk(:,:)
    integer :: jb, ib, g
    call pf_init()
    if (.not. allocated(btop_GB)) allocate(btop_GB(ng, nb(2), nb(3), 3, nBlockMax))
    if (size(btop_GB, 5) /= nBlockMax) then
      deallocate(btop_GB)
      allocate(btop_GB(ng, nb(2), nb(3), 3, nBlockMax))
    end if
    allocate(bk(size(pf_ex, 1), size(pf_ey, 1)))
    call pf_coeffs(n_r, bk)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (.not. arbor_physical_side(ib, 1, 2)) cycle
      call cz_bind(ib)
      do g = 1, ng
        call pf_eval(ib, bk, z_c(nb(1) + g) - z_c(nb(1)), btop_GB(g, :, :, :, ib))
      end do
    end do
    call cz_unbind()
    if (b_seed_test .and. .not. btop_ready) call top_test()
    btop_ready = .true.
  contains
    ! b_seed_test, the first step: the ghosts against the analytic mode-(1,1)
    ! potential field of the seed (the top layer still the seed's)
    subroutine top_test()
      real(dp) :: k1x, k1y, kk, e, x, y, err, bmx, ref(3)
      integer :: j, k
      character(len=200) :: msg
      k1x = 2.0_dp*pi/(box_lx_mm*1.0e8_dp);  k1y = 2.0_dp*pi/(box_ly_mm*1.0e8_dp);  kk = hypot(k1x, k1y)
      err = 0.0_dp;  bmx = 0.0_dp
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        if (.not. arbor_physical_side(ib, 1, 2)) cycle
        call cz_bind(ib)
        do g = 1, ng
          e = b_seed_G*exp(-kk*(z_c(nb(1) + g) - z_c(k_join - r_off)))
          do k = 1, nb(3)
            do j = 1, nb(2)
              x = Metric_DB(2, ib)%X_I(j);  y = Metric_DB(3, ib)%X_I(k)
              ref = e*[cos(k1x*x)*cos(k1y*y), k1x/kk*sin(k1x*x)*cos(k1y*y), k1y/kk*cos(k1x*x)*sin(k1y*y)]
              err = max(err, maxval(abs(btop_GB(g, j, k, :, ib) - ref)))
              bmx = max(bmx, maxval(abs(ref)))
            end do
          end do
        end do
      end do
      call cz_unbind()
      err = arbor_max(err);  bmx = arbor_max(bmx)
      write(msg, '(A,ES10.3,A,ES10.3,A)') '  b_seed_test: the potential top ghosts vs analytic: max error ', err, &
          ' G (|B| up to ', bmx, ' G)'
      call arbor_print(msg)
    end subroutine top_test
  end subroutine cz_top_potential

  ! Every local block's interior state (all nv variables), one file per rank:
  ! out_dir/state_export_rNNNN.bin, stream access: n_block, nb(1:3), nv; then
  ! per block its key (KEY_LEN integers) and State(1:nb(1), 1:nb(2), 1:nb(3), 1:nv).
  ! The root also writes out_dir/state_export_meta.bin: t, the step, and the
  ! restart numbers (profile counter and schedule, the shell clocks: as
  ! cz_restart_save), so an import continues the run.
  module subroutine export_state_write(t, iStep)
    real(dp), intent(in) :: t
    integer, intent(in) :: iStep
    integer :: u, jb, ib
    character(len=300) :: name
    real(dp), allocatable :: Data_I(:)
    if (.not. export_state) return
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/state_export_r', my_rank, '.bin'
    open(newunit=u, file=trim(name), access='stream', form='unformatted', status='replace', action='write')
    write(u) nBlock, nb, nv
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      write(u) Key_IB(:, ib), State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:nv, ib)
    end do
    close(u)
    if (.not. is_root) return
    Data_I = [real(i_prof, dp), t_next_prof, t_prof0, tloc]
    open(newunit=u, file=trim(out_dir)//'/state_export_meta.bin', access='stream', form='unformatted', &
        status='replace', action='write')
    write(u) t, iStep, size(Data_I), Data_I
    close(u)
  end subroutine export_state_write

  ! Index of the export files in import_dir (every rank: the keys only, read
  ! by seeking): file, position and key of every exported block
  module subroutine import_open()
    integer :: u, n, nbf(3), nvf, f, k, m
    integer(i8) :: rec, pos
    character(len=300) :: name
    logical :: exists
    integer, allocatable :: key_K(:,:), file_K(:)
    integer(i8), allocatable :: pos_K(:)
    if (allocated(imp_key)) return
    m = 0
    allocate(key_K(KEY_LEN, 0), file_K(0), pos_K(0))
    do f = 0, 99999
      write(name, '(A,A,I4.4,A)') trim(import_dir), '/state_export_r', f, '.bin'
      inquire(file=trim(name), exist=exists)
      if (.not. exists) exit
      open(newunit=u, file=trim(name), access='stream', form='unformatted', status='old', action='read')
      read(u) n, nbf, nvf
      if (f == 0) then
        imp_nb = nbf
        imp_nv = nvf
      else if (any(nbf /= imp_nb) .or. nvf /= imp_nv) then
        call arbor_abort('cz3d: import: the export files differ in block size or variables')
      end if
      rec = 4_i8*KEY_LEN + 8_i8*product(int(nbf, i8))*nvf
      key_K = reshape([key_K, [(0, k = 1, KEY_LEN*n)]], [KEY_LEN, m + n])
      file_K = [file_K, [(f, k = 1, n)]]
      pos_K = [pos_K, [(0_i8, k = 1, n)]]
      do k = 1, n
        pos = 21_i8 + int(k - 1, i8)*rec                ! (stream positions count from 1; the header is 5 integers)
        read(u, pos=pos) key_K(:, m + k)
        pos_K(m + k) = pos
      end do
      m = m + n
      close(u)
    end do
    if (m == 0) call arbor_abort('cz3d: import: no state_export files in '//trim(import_dir))
    if (any(imp_nb(2:3) /= nb(2:3))) call arbor_abort('cz3d: import: the horizontal block size differs')
    ! the export's top (cells): a taller run (the cascade) starts its blocks
    ! above it from the init file (import_block)
    imp_top = maxval(key_K(KI+1, 1:m) + 1)*imp_nb(1)
    if (mod(imp_nb(1), nb(1)) /= 0) call arbor_abort('cz3d: import: the radial block size must divide the exported one')
    call move_alloc(key_K, imp_key)
    call move_alloc(file_K, imp_file)
    call move_alloc(pos_K, imp_pos)
  end subroutine import_open

  ! The finest horizontal level of the exported blocks
  module integer function import_max_level() result(l)
    l = maxval(imp_key(KL+2, :))
  end function import_max_level

  ! Refinement criterion for the import tree (arbor_tree_refine): refine a
  ! block until an exported block has its horizontal position and level
  module subroutine import_criterion(key, x_lo, x_hi, flags)
    integer, intent(in) :: key(KEY_LEN)
    real(dp), intent(in) :: x_lo(3), x_hi(3)
    integer, intent(out) :: flags(3)
    flags = 0
    if (import_record(key, .false.) == 0) flags(2:3) = 1
  end subroutine import_criterion

  ! The exported block holding block key (horizontal position and level
  ! equal; with radial: also the radial block that contains it), 0 if none
  integer function import_record(key, radial) result(r)
    integer, intent(in) :: key(KEY_LEN)
    logical, intent(in) :: radial
    integer :: m
    r = 0
    do m = 1, size(imp_key, 2)
      if (imp_key(KP, m) /= key(KP)) cycle
      if (any(imp_key(KL+2:KL+3, m) /= key(KL+2:KL+3))) cycle
      if (any(imp_key(KI+2:KI+3, m) /= key(KI+2:KI+3))) cycle
      if (radial) then
        if (imp_key(KI+1, m) /= (key(KI+1)*nb(1))/imp_nb(1)) cycle
      end if
      r = m
      return
    end do
  end function import_record

  ! The block's state from the export files in import_dir: its radial part of
  ! the exported block at the same place (the export may hold whole columns
  ! or larger radial blocks); the variables both runs have are copied, the
  ! others stay as set (B, psi: 0). With the atmosphere, a block above the
  ! join without an exported block keeps its initial state (a CZ-only export,
  ! e.g. CZ_3D's, under the atmosphere of the init file)
  subroutine import_block(iBlock)
    integer, intent(in) :: iBlock
    real(dp), allocatable :: buf(:,:,:,:)
    integer :: m, u, i0
    character(len=300) :: name
    call import_open()
    m = import_record(Key_IB(:, iBlock), .true.)
    if (m == 0 .and. use_atmosphere .and. Key_IB(KI+1, iBlock)*nb(1) >= k_join) return
    ! the cascade: a block above the export's top keeps the init file's state
    if (m == 0 .and. Key_IB(KI+1, iBlock)*nb(1) >= imp_top) return
    if (m == 0) call arbor_abort('cz3d: import: no exported block at this place (a different grid?)')
    allocate(buf(imp_nb(1), imp_nb(2), imp_nb(3), imp_nv))
    write(name, '(A,A,I4.4,A)') trim(import_dir), '/state_export_r', imp_file(m), '.bin'
    open(newunit=u, file=trim(name), access='stream', form='unformatted', status='old', action='read')
    read(u, pos=imp_pos(m) + 4_i8*KEY_LEN) buf
    close(u)
    i0 = Key_IB(KI+1, iBlock)*nb(1) - imp_key(KI+1, m)*imp_nb(1)
    State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:min(nv, imp_nv), iBlock) = buf(i0+1:i0+nb(1), :, :, 1:min(nv, imp_nv))
  end subroutine import_block

  ! import_cz_on_restart: the CZ blocks (below the join) from the export
  module subroutine cz_import_cz()
    integer :: jb, ib, n
    character(len=200) :: msg
    n = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (Key_IB(KI+1, ib)*nb(1) + nb(1) > k_join) cycle
      call import_block(ib)
      n = n + 1
    end do
    n = arbor_sum(n)
    write(msg, '(A,I0,A,A)') '  import_cz_on_restart: ', n, ' CZ blocks from ', trim(import_dir)
    call arbor_print(msg)
  end subroutine cz_import_cz

  ! atm_match_base: each column's lower atmosphere pressure-matched to its
  ! own top CZ cell (see cz3d_params_mod; after an import, before the field
  ! seed). f per column from the blocks holding cell k_join, gathered in a
  ! (1, x, y) map at the run's horizontal level; the atmosphere
  ! cells below atm_match_z_km(2) then get init_file's density times
  ! f_eff(z) at init_file's T (s from the EOS), v unchanged
  module subroutine cz_match_base()
    real(dp), allocatable :: fm(:,:,:)
    real(dp) :: out(n_field), lnp_ref, rho, rp, Tp, sp, rn, sn, zk, w, fe, fmin, fmax, qmin
    integer :: nth, nph, jb, ib, i, j, k, jg, kg, ig, n
    character(len=200) :: msg
    nth = n_x*2**refine_h
    nph = n_y*2**refine_h
    allocate(fm(1, nth, nph))
    fm = 0.0_dp
    call cz_unbind()
    call eos_eval(tab, rho0(k_join) + rho1_init(k_join), s0(k_join) + s1_init(k_join), out)
    lnp_ref = out(F_LNP)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (Key_IB(KL+2, ib) /= refine_h .or. Key_IB(KL+3, ib) /= refine_h) &
          call arbor_abort('cz3d: atm_match_base needs every block at the horizontal level refine_h')
      i = k_join - Key_IB(KI+1, ib)*nb(1)
      if (i < 1 .or. i > nb(1)) cycle
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
          call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
          fm(Key_IB(KP, ib), Key_IB(KI+2, ib)*nb(2) + j, Key_IB(KI+3, ib)*nb(3) + k) = exp(out(F_LNP) - lnp_ref)
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(fm)
    fmin = minval(fm);  fmax = maxval(fm)
    if (fmin <= 0.0_dp) call arbor_abort('cz3d: atm_match_base: a column without its top CZ cell')
    qmin = huge(1.0_dp);  n = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do i = 1, nb(1)
        ig = Key_IB(KI+1, ib)*nb(1) + i
        zk = z_c(i)/1.0e5_dp
        if (ig <= k_join .or. zk >= atm_match_z_km(2)) cycle
        w = 1.0_dp
        if (zk > atm_match_z_km(1)) w = cos(0.5_dp*pi*(zk - atm_match_z_km(1))/(atm_match_z_km(2) - atm_match_z_km(1)))**2
        rp = rho0(i) + rho1_init(i)
        sp = s0(i) + s1_init(i)
        call eos_eval(tab, rp, sp, out)
        Tp = exp(out(F_LNT))
        do k = 1, nb(3)
          do j = 1, nb(2)
            jg = Key_IB(KI+2, ib)*nb(2) + j
            kg = Key_IB(KI+3, ib)*nb(3) + k
            fe = 1.0_dp + (fm(Key_IB(KP, ib), jg, kg) - 1.0_dp)*w
            rn = fe*rp
            sn = eos_s_from_rho_T(tab, rn, Tp, s_guess=sp)
            State_GVB(i, j, k, rho1_, ib) = (rn - rho0(i))/rq(i)
            State_GVB(i, j, k, s1_, ib) = sn - s0(i)
            qmin = min(qmin, rn/rho0(i))
            n = n + 1
          end do
        end do
      end do
    end do
    call cz_unbind()
    qmin = arbor_min(qmin);  n = arbor_sum(n)
    write(msg, '(A,F7.4,A,F7.4,A,I0,A,F6.0,A,F6.0,A,F7.4)') '  atm_match_base: f = p(top CZ cell)/p(init_file) in ', fmin, &
        ' .. ', fmax, '; ', n, ' atmosphere cells rescaled (full below ', atm_match_z_km(1), ' km, none above ', &
        atm_match_z_km(2), ' km); min rho/rho0 there ', qmin
    call arbor_print(msg)
    if (qmin < 0.9_dp) call arbor_print('  atm_match_base: WARNING rho < 0.9 rho0 in the lower atmosphere (the background rule rho0 <= rho)')
  end subroutine cz_match_base

  ! The time, step and restart numbers of the export (state_export_meta.bin),
  ! if there: t and iStep set, the profile counter and the shell clocks
  ! restored; found false (nothing changed) without the file
  module subroutine import_meta(t, iStep, found)
    real(dp), intent(inout) :: t
    integer, intent(inout) :: iStep
    logical, intent(out) :: found
    integer :: u, nd
    real(dp), allocatable :: Data_I(:)
    inquire(file=trim(import_dir)//'/state_export_meta.bin', exist=found)
    if (.not. found) return
    open(newunit=u, file=trim(import_dir)//'/state_export_meta.bin', access='stream', form='unformatted', &
        status='old', action='read')
    read(u) t, iStep, nd
    allocate(Data_I(nd))
    read(u) Data_I
    close(u)
    i_prof = nint(Data_I(1))
    t_next_prof = Data_I(2)
    t_prof0 = Data_I(3)
    if (nd == 3 + n_r) then
      tloc = Data_I(4:)
    else
      tloc = t                                            ! (a different radial grid: all clocks at the run's time)
    end if
  end subroutine import_meta

  ! Uniform in [-1, 1] from (patch, i, j, k): Park-Miller steps (no overflow)
  real(dp) function urand(p, i, j, k)
    integer, intent(in) :: p, i, j, k
    integer(i8), parameter :: a = 48271_i8, m = 2147483647_i8
    integer(i8) :: x
    integer :: n
    x = 12345_i8 + 999331_i8*int(p, i8)
    do n = 1, 3
      x = mod(a*mod(x + int(i, i8), m), m)
      x = mod(a*mod(x + 7919_i8*int(j, i8), m), m)
      x = mod(a*mod(x + 104729_i8*int(k, i8), m), m)
    end do
    urand = 2.0_dp*real(x, dp)/real(m, dp) - 1.0_dp
  end function urand

  ! si_radial: the columns (patch, theta block, phi block) to ranks, each
  ! column's radial blocks on one rank (the column solves need no messages).
  ! The columns in the order their first blocks appear on the tree's Hilbert
  ! curve, cut into n_ranks contiguous stretches; the blocks move with their
  ! data (arbor_redistribute). All blocks at one horizontal level.
  module subroutine cz_columns_to_ranks(tree)
    type(arbor_tree_t), intent(in) :: tree
    integer, allocatable :: cid(:), pos(:), own(:), lev(:), per(:)
    integer :: j, nth_b, nph_b, ncol, n
    character(len=200) :: msg
    associate (key => tree%key(:, 1:tree%n_leaf))
    if (any(key(KL+2, :) /= key(KL+2, 1)) .or. any(key(KL+3, :) /= key(KL+2, 1))) &
        call arbor_abort('cz3d: si_radial needs every block at one horizontal level')
    nth_b = maxval(key(KI+2, :)) + 1
    nph_b = maxval(key(KI+3, :)) + 1
    ncol = maxval(key(KP, :))*nth_b*nph_b
    end associate
    allocate(cid(tree%n_leaf), pos(0:ncol-1), own(tree%n_leaf), lev(tree%n_leaf), per(0:n_ranks-1))
    pos = -1
    n = 0
    do j = 1, tree%n_leaf
      cid(j) = ((tree%key(KP, j) - 1)*nth_b + tree%key(KI+2, j))*nph_b + tree%key(KI+3, j)
      if (pos(cid(j)) < 0) then
        pos(cid(j)) = n
        n = n + 1
      end if
    end do
    per = 0
    do j = 0, ncol - 1
      if (pos(j) >= 0) per(int((int(pos(j), 8)*n_ranks)/n)) = per(int((int(pos(j), 8)*n_ranks)/n)) + 1
    end do
    do j = 1, tree%n_leaf
      own(j) = int((int(pos(cid(j)), 8)*n_ranks)/n)
    end do
    lev = 1
    call arbor_redistribute(tree, own, lev)
    call arbor_ghost_rebuild(tree)
    si_cols = .true.
    write(msg, '(A,I0,A,I0,A,I0,A,I0,A)') '  si_radial: ', n, ' columns to ', n_ranks, ' ranks (', minval(per), &
        ' .. ', maxval(per), ' per rank)'
    call arbor_print(msg)
  end subroutine cz_columns_to_ranks

  ! The implicit range (whole radial blocks covering si_z_km), the local
  ! columns and the per-block arrays (after every (re)distribution)
  subroutine si_setup()
    integer :: g0, g1, i, jb, ib, c, q, nk
    character(len=200) :: msg
    call cz_unbind()
    g0 = 1;  g1 = n_r
    do i = 1, n_r
      if (z_c(i) >= si_z_km(1)*1.0e5_dp) then
        g0 = i
        exit
      end if
    end do
    do i = n_r, 1, -1
      if (z_c(i) <= si_z_km(2)*1.0e5_dp) then
        g1 = i
        exit
      end if
    end do
    if (g1 < g0) call arbor_abort('cz3d: si_z_km holds no cell')
    si_kb0 = (g0 - 1)/nb(1)
    si_kb1 = (g1 - 1)/nb(1)
    si_i0 = si_kb0*nb(1) + 1
    si_i1 = (si_kb1 + 1)*nb(1)
    nk = si_kb1 - si_kb0 + 1
    if (allocated(si_cb)) deallocate(si_cb)
    allocate(si_cb(nk, nBlock))
    si_cb = 0
    si_ncol = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (.not. si_in(ib)) cycle
      c = 0
      do q = 1, si_ncol
        ib_c: block
          integer :: ib1
          ib1 = maxval(si_cb(:, q))
          if (all(Key_IB([KP, KI+2, KI+3], ib1) == Key_IB([KP, KI+2, KI+3], ib))) c = q
        end block ib_c
        if (c > 0) exit
      end do
      if (c == 0) then
        si_ncol = si_ncol + 1
        c = si_ncol
      end if
      si_cb(Key_IB(KI+1, ib) - si_kb0 + 1, c) = ib
    end do
    if (any(si_cb(:, 1:si_ncol) == 0)) &
        call arbor_abort('cz3d: si_radial: a column is not whole on its rank (cz_columns_to_ranks not called?)')
    if (allocated(siS_B)) deallocate(siS_B, siR_B)
    allocate(siS_B(nb(1), nb(2), nb(3), nBlockMax), siR_B(nb(1), nb(2), nb(3), nv, nBlockMax))
    siS_B = 0.0_dp
    if (allocated(siSa_B)) deallocate(siSa_B)
    if (si_alfven) then
      allocate(siSa_B(nb(1), nb(2), nb(3), nBlockMax))
      siSa_B = 0.0_dp
    end if
    if (allocated(siV_B)) deallocate(siV_B, siX_B)
    if (si_diff_precond) allocate(siV_B(nb(1), nb(2), nb(3), nBlockMax), siX_B(nb(1), nb(2), nb(3), nBlockMax))
    if (allocated(EtaB_GB)) deallocate(EtaB_GB)
    if (use_B .and. use_diffusion .and. b_diff_form == 'resistive') then
      if (si_i0 /= 1 .or. si_i1 /= n_r) &
          call arbor_abort("cz3d: b_diff_form = 'resistive' needs the whole column implicit (si_z_km default)")
      allocate(EtaB_GB(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, 1, size(State_GVB, 5)))
      EtaB_GB = 0.0_dp
    end if
    write(msg, '(A,I0,A,I0,A,F9.1,A,F10.1,A)') '  si_radial: implicit cells ', si_i0, ' .. ', si_i1, ' (z ', &
        z_c(si_i0)/1.0e5_dp, ' .. ', z_c(si_i1)/1.0e5_dp, ' km; the two cells at each end explicit)'
    call arbor_print(msg)
  end subroutine si_setup

  ! b_smooth_cells > 0 (see cz3d_params_mod): the CZ's field smoothed before
  ! the potential seed, each component of B on each layer by the same
  ! explicit horizontal diffusion (independent of z, so it commutes with the
  ! vertical part of div B): per pass
  !   B += D (d_x^2 B + d_y^2 B),  D = 0.1 min(dx, dy)^2
  ! n = b_smooth_cells^2/0.2 passes (diffusion length b_smooth_cells cells).
  ! The CZ cells only (with the atmosphere: up to the join); psi set to 0.
  module subroutine cz_smooth_b()
    integer :: n, it, jb, ib, i, j, k, c
    real(dp) :: acc(4), t0, w, dth, dph, st, sp, dd
    real(dp), allocatable :: B0(:,:,:,:,:), X(:,:,:)
    character(len=200) :: msg
    if (b_smooth_cells <= 0.0_dp) return
    t0 = arbor_wtime()
    n = max(1, ceiling(b_smooth_cells**2/0.2_dp))
    allocate(B0(nb(1), nb(2), nb(3), br_:bp_, size(State_GVB, 5)), X(nb(1), nb(2), nb(3)))
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      B0(:, :, :, :, ib) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib)
    end do
    do it = 1, n
      call arbor_fill_ghosts(State_GVB, iVarMin=br_, iVarMax=bp_, n_ghost=1, faces_only=.true.)
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        if (use_atmosphere .and. Key_IB(KI+1, ib)*nb(1) >= k_join) cycle
        associate (dth_ => Metric_DB(2, ib)%Dx_I, dph_ => Metric_DB(3, ib)%Dx_I)
        do c = br_, bp_
          do k = 1, nb(3)
            do j = 1, nb(2)
              dth = dth_(j)
              dph = dph_(k)
              dd = 0.1_dp*min(dth, dph)**2
              X(:, j, k) = State_GVB(1:nb(1), j, k, c, ib) + dd*( &
                  (State_GVB(1:nb(1), j+1, k, c, ib) - 2.0_dp*State_GVB(1:nb(1), j, k, c, ib) &
                  + State_GVB(1:nb(1), j-1, k, c, ib))/dth**2 &
                  + (State_GVB(1:nb(1), j, k+1, c, ib) - 2.0_dp*State_GVB(1:nb(1), j, k, c, ib) &
                  + State_GVB(1:nb(1), j, k-1, c, ib))/dph**2)
            end do
          end do
          do i = 1, nb(1)
            if (use_atmosphere .and. Key_IB(KI+1, ib)*nb(1) + i > k_join) cycle
            State_GVB(i, 1:nb(2), 1:nb(3), c, ib) = X(i, :, :)
          end do
        end do
        end associate
      end do
    end do
    ! report: rms |B - B0| / rms |B| over the CZ's top 3 Mm and over the whole CZ (volume weighted)
    acc = 0.0_dp
    call cz_unbind()
    w = r_f(merge(k_join, n_r, use_atmosphere))
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      State_GVB(:, :, :, psi_, ib) = 0.0_dp
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            if (use_atmosphere .and. r_off + i > k_join) cycle
            dd = (r_f(i) - r_f(i-1))*Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
            st = sum((State_GVB(i, j, k, br_:bp_, ib) - B0(i, j, k, :, ib))**2)
            sp = sum(B0(i, j, k, :, ib)**2)
            acc(3:4) = acc(3:4) + dd*[st, sp]
            if (w - r_c(i) < 3.0e8_dp) acc(1:2) = acc(1:2) + dd*[st, sp]
          end do
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(acc)
    write(msg, '(A,F6.2,A,I0,A,ES10.3,A,ES10.3,A,F8.2,A)') '  b_smooth: ', b_smooth_cells, ' cells, ', n, &
        ' passes: rms|dB|/rms|B| ', sqrt(acc(1)/max(acc(2), tiny(1.0_dp))), ' (top 3 Mm of the CZ), ', &
        sqrt(acc(3)/max(acc(4), tiny(1.0_dp))), ' (CZ); ', arbor_wtime() - t0, ' s'
    call arbor_print(msg)
  end subroutine cz_smooth_b

  ! Background enthalpy h0(r) (EOS core) and potential Phi(r) (Phi(1) = 0)
  subroutine background_h_phi()
    type(eos_state) :: st
    real(dp) :: out(n_field)
    integer :: i
    allocate(h0_r_g(-1:n_r+2), phi_r_g(-1:n_r+2))
    h0_r => h0_r_g
    phi_r => phi_r_g
    do i = -1, n_r + 2
      call eos_eval(tab, rho0(i), s0(i), out)
      call eos_get_from_rho_T(rho0(i), exp(out(F_LNT)), st, tab%comp, tab%coulomb)
      h0_r(i) = st%e + st%p/st%rho
    end do
    phi_r(1) = 0.0_dp
    do i = 1, n_r + 1
      phi_r(i+1) = phi_r(i) + 0.5_dp*(g_c(i) + g_c(i+1))*(r_c(i+1) - r_c(i))
    end do
    do i = 0, -1, -1
      phi_r(i) = phi_r(i+1) - 0.5_dp*(g_c(i) + g_c(i+1))*(r_c(i+1) - r_c(i))
    end do
  end subroutine background_h_phi

  ! The profile counter and schedule go into restart files
  module subroutine cz_restart_save(self, Data_I)
    class(cz3d_app_t), intent(inout) :: self
    real(dp), allocatable, intent(out) :: Data_I(:)
    Data_I = [real(i_prof, dp), t_next_prof, t_prof0, tloc]
  end subroutine cz_restart_save

  module subroutine cz_restart_load(self, Data_I)
    class(cz3d_app_t), intent(inout) :: self
    real(dp), intent(in) :: Data_I(:)
    if (size(Data_I) /= 2 .and. size(Data_I) /= 3 .and. size(Data_I) /= 3 + n_r) &
        call arbor_abort('cz3d: restart file from another application')
    i_prof = nint(Data_I(1))
    t_next_prof = Data_I(2)
    t_prof0 = 0.0_dp
    if (size(Data_I) >= 3) t_prof0 = Data_I(3)           ! (files before 2026-10-11 hold two numbers)
    if (.not. allocated(tloc)) allocate(tloc(1:n_r))
    tloc = -1.0_dp                                        ! (older files: set to the run's time by the caller)
    if (size(Data_I) == 3 + n_r) tloc = Data_I(4:)
  end subroutine cz_restart_load

end submodule cz3d_setup_smod
