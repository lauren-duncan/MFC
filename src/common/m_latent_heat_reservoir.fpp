!> Latent heat reservoir method for phase change
!> Store the energy exceeding internal energy at saturation state to a
!> reservoir \Lambda; when the stored energy is greater than the latent
!> heat of vaporization, phase change occurs
module m_latent_heat_reservoir

#ifndef MFC_POST_PROCESS

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use ieee_arithmetic

    ! ==========================================================================

    implicit none

    private; 
    public :: s_initialize_reservoir_module, &
              s_latent_heat_reservoir, &
              s_finalize_reservoir_solver_module


    !> @name Parameters required for latent heat reservoir method
    !> @{
    integer, parameter :: max_iter = 1e8        !< max # of iterations
    real(kind(0d0)), parameter :: pCr = 4.94d7   !< Critical water pressure
    real(kind(0d0)), parameter :: TCr = 385.05 + 273.15  !< Critical water temperature
    real(kind(0d0)), parameter :: mixM = 1.0d-8 !< threshold for 'mixture cell'.
    integer, parameter :: lp = 1    !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2    !< index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    real(kind(0d0)) :: A, B, C, D
    !> @}

    !$acc declare create(max_iter,pCr,TCr,mixM,lp,vp,A,B,C,D)


contains

    !>  The purpose of this subroutine is to initialize the reservoir module

    subroutine s_initialize_reservoir_module
        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) &
             + qvps(vp) - qvps(lp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

        D = ((gs_min(lp) - 1.0d0)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

    end subroutine s_initialize_reservoir_module

    !>  This subroutine is created to activate the latent heat reservoir
        !!  method for phase transition.
    subroutine s_latent_heat_reservoir(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0.0d0)) rM, TvF, MCT !< total reacting mass and volume
        real(kind(0.0d0)) rMb, TvFb  !< total reacting mass and volume after phase change
        real(kind(0.0d0)) pS, alphal, alphav, rhol, rhov, Tl, Tv, el, ev, Lambda
                          !< vapor and liquid states
        real(kind(0.0d0)) TSat, LSat, el_Sat  !< saturation states
        !$acc declare create(rM, TvF, MCT, rMb, TvFb, pS, alphal, alphav, rhol, rhov, Tl, Tv, el, ev, Lambda, TSat, LSat, sl_Sat)


        !< Generic loop iterators
        integer :: i, j, k, l

        ! starting reservoir solver
        !! ! compute pressure by converting conservative to primitive variables
        !! call s_convert_conservative_to_primitive_variables(q_cons_vf, q_prim_vf)

        !$acc parallel loop collapse(3) gang vector default(present) private(p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok,pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, dynE, rhos, rho, rM, m1, m2, MCT, TvF)
        do j = 0, m
            do k = 0, n
                do l = 0, p

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process.
                    rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! calculating the total reacting volume fraction.
                    ! This helps for update volume fraction after phase change
                    TvF = q_cons_vf(lp + advxb - 1)%sf(j, k, l) + q_cons_vf(vp + advxb - 1)%sf(j, k, l)

                    ! correcting negative (recating) mass fraction values in case they happen
                    call s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)


                    ! calculations of the states of the reacting liquid
                    ! and vapor
                    pS = q_prim_vf(E_idx)%sf(j, k, l) 
                    alphal = q_cons_vf(lp + adv_idx%beg - 1)%sf(j, k, l)
                    alphav = q_cons_vf(vp + adv_idx%beg - 1)%sf(j, k, l)

                    if (alphal > 0.0d0) then

                        rhol = q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) / &
                             q_cons_vf(lp + adv_idx%beg - 1)%sf(j, k, l)
                        rhov = q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) / &
                             q_cons_vf(vp + adv_idx%beg - 1)%sf(j, k, l)
                        Tl = (pS + ps_inf(lp))*gammas(lp)/cvs(lp)/rhol
                        Tv = (pS + ps_inf(vp))*gammas(vp)/cvs(vp)/rhov
                        el = (pS * gammas(lp) + pi_infs(lp))/rhol + qvs(lp)
                        ev = (pS * gammas(vp) + pi_infs(vp))/rhov + qvs(vp)
                        Lambda = q_cons_vf(lam_idx)%sf(j, k, l)/alphal/rhol

                        if((pCr < pS) .and. (TCr < Tl)) then
                            print *, "Warning: Exceed critical point!"
                        end if 

                    !<  ! Mixture density
                    !<  rho = 0.0d0
                    !<  !$acc loop seq
                    !<  do i = 1, num_fluids
                    !<      rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)
                    !<  end do
                    !<  ! kinetic energy
                    !<  dynE = 0.0d0;
                    !<  !$acc loop seq
                    !<  do i = momxb, momxe
                    !<      dynE = dynE + 5.0d-1*q_cons_vf(i)%sf(j, k, l)**2/rho     
                    !<  end do
                    !<  ! calculate the internal energy as the total energy
                    !<  ! minus the kinetic energy
                    !<  rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                        ! calculations of saturation temperature and latent
                        ! heat for vaporization under current pressure
                        call s_TSat(pS, TSat, Tl)
                        call s_LSat(TSat, LSat)

                        ! calculate the internal energy for liquid under the
                        ! vaporization temperature
                        el_Sat = ((pS + gs_min(lp)*ps_inf(lp))/(pS + ps_inf(lp)))*  &
                             cvs(lp)*TSat + qvs(lp)

                        ! compare current temperature
                        if (Tl < TSat)  then
                                 if(Lambda > 0.0d0) then
                                ! move energy from reservoir to increase
                                ! the temperature
                                el = el + MIN(el_Sat - el, Lambda)
                                Lambda = Lambda - MIN(el_Sat - el, Lambda)
                                end if
                        else
                                 ! move additional energy to the reservoir
                                 el = el_Sat
                                 Lambda = Lambda + (el - el_Sat) 
                        end if


                        ! compare reservoir with the latent heat of
                        ! vaporization
                        if (Lambda >= LSat) then
                                ! phase change occurs
                                ! release the energy and update the state variables 
                                rhov = (alphav * rhov + alphal *rhol)/(alphav + alphal)
                                ev = (alphav * ev + alphal * (el + Lambda))/(alphav + alphal)
                                rhol = 0.0d0
                                el = 0.0d0
                                alphav = alphav + alphal
                                alphal = 0.0d0
                                Lambda = 0.0d0
                        end if

                        !Update the conservative variables
                        q_cons_vf(lp + adv_idx%beg - 1)%sf(j, k, l) = alphal
                        q_cons_vf(vp + adv_idx%beg - 1)%sf(j, k, l) = alphav
                        q_cons_vf(lp + cont_idx%beg - 1)%sf(j, k, l) = alphal * rhol
                        q_cons_vf(vp + cont_idx%beg - 1)%sf(j, k, l) = alphav * rhov
                        q_cons_vf(lp + intxb - 1)%sf(j, k, l) = q_cons_vf(lp + contxb - 1)%sf(j, k, l)*el
                        q_cons_vf(vp + intxb - 1)%sf(j, k, l) = q_cons_vf(vp + contxb - 1)%sf(j, k, l)*ev
                        q_cons_vf(lam_idx)%sf(j, k, l) = alphal * rhol * Lambda


                        !Check total volume/mass/energy conservation
                        !By hypothesis, these should not change after
                        !phase change process

                        TvFb = q_cons_vf(lp + advxb - 1)%sf(j, k, l) + q_cons_vf(vp + advxb - 1)%sf(j, k, l)
                        if (abs(TvF - TvFb) > 1e-3) then
                                print *, "Warning: Total volume fraction of the &
                                reaacting fluids are not conservative after phase change"
                        end if

                        rMb = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)
                        if (abs(rM - rMb) > 1e-3) then
                                print *, "Warning: Total mass of the &
                                reaacting fluids are not conservative after phase change"
                        end if



                    end if

                end do
            end do
        end do

    end subroutine s_latent_heat_reservoir


    !>  This auxiliary subroutine corrects the partial densities of the REACTING fluids in case one of them is negative
        !!      but their sum is positive. Inert phases are not corrected at this moment
        !!  @param MCT partial density correction parameter
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rM sum of the reacting masses
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
    subroutine s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)
        !$acc routine seq

        !> @name variables for the correction of the reacting partial densities
        !> @{
        real(kind(0.0d0)), intent(out) :: MCT
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(kind(0.0d0)), intent(inout) :: rM
        integer, intent(in) :: j, k, l
        !> @}
        if (rM < 0.0d0) then

            if ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) >= -1.0d0*mixM) .and. &
                (q_cons_vf(vp + contxb - 1)%sf(j, k, l) >= -1.0d0*mixM)) then

                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = 0.0d0

                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = 0.0d0

                rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

            end if

        end if

        ! Defining the correction in terms of an absolute value might not be the best practice.
        ! Maybe a good way to do this is to partition the partial densities, giving a small percentage of the total reacting density
        MCT = 2*mixM

        ! correcting the partial densities of the reacting fluids. What to do for the nonreacting ones?
        if (q_cons_vf(lp + contxb - 1)%sf(j, k, l) < 0.0d0) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = MCT*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

        elseif (q_cons_vf(vp + contxb - 1)%sf(j, k, l) < 0.0d0) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = MCT*rM

        end if
    end subroutine s_correct_partial_densities

    !>  This auxiliary subroutine finds the Saturation temperature for a given
        !!      saturation pressure through a newton solver
        !!  @param pSat Saturation Pressure
        !!  @param TSat Saturation Temperature
        !!  @param TSIn equilibrium Temperature
    subroutine s_TSat(pSat, TSat, TSIn)
        !$acc routine seq

        real(kind(0.0d0)), intent(in) :: pSat
        real(kind(0.0d0)), intent(out) :: TSat
        real(kind(0.0d0)), intent(in) :: TSIn

        real(kind(0.0d0)) :: dFdT, FT, Om !< auxiliary variables

        ! Generic loop iterators
        integer :: ns

        if ((pSat == 0.0d0) .and. (TSIn == 0.0d0)) then

            ! assigning Saturation temperature
            TSat = 0.0d0

        else

            ! calculating initial estimate for temperature in the TSat procedure. I will also use this variable to
            ! iterate over the Newton's solver
            TSat = TSIn

            ! iteration counter
            ns = 0

            ! underrelaxation factor
            Om = 1.0d-3
            do while ((DABS(FT) > ptgalpha_eps) .or. (ns == 0))
                ! increasing counter
                ns = ns + 1

                ! calculating residual
                FT = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                           *(1 - DLOG(TSat)) - (qvps(lp) - qvps(vp)) &
                           + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

                ! calculating the jacobian
                dFdT = &
                    -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*DLOG(TSat) &
                    - (qvps(lp) - qvps(vp)) &
                    + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                    - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))

                ! updating saturation temperature
                TSat = TSat - Om*FT/dFdT

            end do

        end if

    end subroutine s_TSat


    !>  This auxiliary subroutine finds the latent heat of vaporization
        !!  for a given saturation temperature
        !!  @param TSat Saturation Temperature
        !!  @param LSat Latent heat of vaporization
    subroutine s_LSat(TSat, LSat)
        !$acc routine seq

        real(kind(0.0d0)), intent(in) :: TSat
        real(kind(0.0d0)), intent(out) :: LSat

        real(kind(0.0d0)) :: hl, hv
 
        hl = gs_min(lp)*cvs(lp)*TSat + qvs(lp)
        hv = gs_min(vp)*cvs(vp)*TSat + qvs(vp)

        LSat = hv - hl;

    end subroutine s_LSat



    !>  This subroutine finalizes the latent heat reservoir module
    subroutine s_finalize_reservoir_solver_module
    end subroutine

#endif

end module m_latent_heat_reservoir
