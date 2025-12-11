! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! START: ICE ALLOCATION
! ----------------------------------------------------------------------
! File: SIA2_env_ice_allocation.inc.f90
! Purpose: maps ice to new boundaries
!
! INPUT VARIALBES:
!
! ice - current physical ice vars
! bio - current biology vars
! int_z - current internal ice depth
! th_new - new thicknesses
! id_new - new ice depths of layers
! int_z_new - the new number of internal layers
! new_layer_top - where among the old layers that the new ice column should
!            begin (used during sublimation)
! c_gl - congelation ice growth layer (m), to be added to bottom internal layer
!         If c_gl is greater than zero, the salinity of the newly grown ice is also required
!         as the variable:
! c_gl_sal - salinity of congelation growth layer
! flooded - if part of the ice column has been flooded, this is the depth of
!         flooding in meters. If flooded is greater than zero, then these
!         additional variables must be defined for assignment to flooded
!         portions of layers:
! T_flood - temp (degC)
! S_flood - salinity (psu)
! d_flood - density of bulk ice (g/m^3)
! bs_flood - brine salinity (psu)
! bd_flood - density of brine (g/m^3)
! bv_flood - brine volume (ppt)
! heat_flood - heat required to bring ice from T_ref to T_flood
!
! OUTPUT VARIALBES:
!
! ice(mi)th
! ice(mi)id
! ice(mi)t
! ice(mi)s
! ice(mi)d
! ice(mi)bs
! ice(mi)bd
! ice(mi)bv
! ice(mi)heat
! ice(mi)smalg
! ice(mi)NO3
! ice(mi)NH4
! ice(mi)PO4
! ice(mi)SiOH4
!
! INTERNAL VARS REQUIRED:
!
! integer :: ii,jj
! double precision :: interim_top,interim_bot,old_layer_top,new_layer_bot
! double precision :: old_layer_bot,z1,z2,dz,heat_total,dz_total
! double precision :: t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean
! double precision, dimension(z_max) :: t_new,s_new,d_new,bs_new,bd_new
! double precision, dimension(z_max) :: bv_new,smalg_new
! double precision, dimension(z_max+1) :: NO3_new,NH4_new,PO4_new,SiOH4_new
!
! ----------------------------------------------------------------------

    ! setup vars to track new and old layer heights
    old_layer_top = 0.
    new_layer_bot = 0.
    old_layer_bot = 0.

    c_gl_copy = c_gl

    ! record flooded into var that we can change during boundary adjustment
    f_depth = flooded*ice(sc,ic,mi)%snow%d(1)/max(IceD/2.d0,ice(sc,ic,mi)%snow%d(1))
    mf_depth = melt_flood

    ! determine if we are to move algae with growing ice
    if ((c_gl_copy .lt. 0.) .or. (alg_mig .eq. 0) .or. (flooded .gt. 0.)) then
        alg_fixed = .true.
    else
        alg_fixed = .false.
    endif

    ! don't grow if depth will be maxed out by growth
    if (.NOT. maxed_depth) then

        ! need to know when re-mapping of old ice column begins - set var to tell when this happens
        remapping = .false.


        debug_z = ice(sc,ic,mi)%no3(1:z_max)

!            if(ice(sc,ic,mi)%s(1) .le. 0.12 .and. ice(sc,ic,mi)%no3(1) .gt. 10.) then
!                    testvar = -1
!                    testvar = sqrt(testvar)
!            endif

        ! calculate algal migration, and store in an array for latter use in
        ! re-mapping
        tmp1 = dhdt_cm_s*86400.
        print *,'dhdt_cm_d',tmp1
        print *,'c_gl_copy',c_gl_copy
            ! short circuit algal migration - move down with any speed of ice growth -- step 1/2
        if ((alg_mig .gt. 0) .and. (c_gl_copy .ne. 0.)) then

            if (alg_mig .eq. 1) then

              ! Algae can move at speed of ice bottom growth/melt up to alg_mig_crit speed
              tmp2 = abs(tmp1)
              if (tmp2 .gt. alg_mig_crit) then
                  tmp2 = c_gl_copy * alg_mig_crit/tmp2
              else
                  tmp2 = c_gl_copy
              endif
              call siesta_remap_tracer(ice(sc,ic,mi)%th,th_new,ice(sc,ic,mi)%smalg, &
                                       tmp2,int_z,int_z_new,.FALSE.,abs(s_gl),smalg_pre_map)


            elseif (alg_mig .eq. 2 .and. dhdt_cm_s .gt. 0.) then
!            elseif (alg_mig .eq. 2) then
              ! *************** Lavioe algal migration, in effect ****************
              ! algae move at linearly decreasing speed up to alg_mig_crit
              tmp2 = abs(tmp1)
              if (tmp2 .gt. alg_mig_crit) then
                  tmp2 = 0.d0  !
              else
                  tmp2 = c_gl_copy*(1.-tmp2/alg_mig_crit)
              endif
              ! *******************************************************************

              call siesta_remap_tracer(ice(sc,ic,mi)%th,th_new,ice(sc,ic,mi)%smalg, &
                                       tmp2,int_z,int_z_new,.FALSE.,abs(s_gl),smalg_pre_map)
            else

                ! default for when algae is not migrating
              call siesta_remap_tracer(ice(sc,ic,mi)%th,th_new,ice(sc,ic,mi)%smalg, &
                                       0.d0,int_z,int_z_new,.FALSE.,abs(s_gl),smalg_pre_map)


            endif

        else
                ! default for when algae is not migrating
              call siesta_remap_tracer(ice(sc,ic,mi)%th,th_new,ice(sc,ic,mi)%smalg, &
                                       0.d0,int_z,int_z_new,.FALSE.,abs(s_gl),smalg_pre_map)

        endif



        do ii=1,int_z_new+z_sk

            ! These vars keep track of new ice physics while layers are mixing
            dz_total = 0.
            d_new(ii) = 0.
            s_new(ii) = 0.
            NO3_new(ii) = 0.
            NH4_new(ii) = 0.
            PO4_new(ii) = 0.
            SiOH4_new(ii) = 0.
            smalg_new(ii) = 0.
            poc_new(ii) = 0.
            heat_total = 0.

            ! flood top layers first, then integrate old ice into new
            if (f_depth .gt. 0.) then
                f_depth = f_depth - th_new(ii)
                if (f_depth .lt. 0.) then
                    dz = f_depth + th_new(ii)
                else
                    dz = th_new(ii)
                endif

                ! record tracers
                tmp1 = S_flood*dz
                tmp2 = tmp1/f(mi)%s
                s_new(ii) = s_new(ii) + tmp1
                d_new(ii) = d_new(ii) + d_flood*dz
                heat_total = heat_total + heat_flood*dz
                NO3_new(ii) = NO3_new(ii) + f(mi)%no3*tmp2
                NH4_new(ii) = NH4_new(ii) + f(mi)%nh4*tmp2
                PO4_new(ii) = PO4_new(ii) + f(mi)%po4*tmp2
                SiOH4_new(ii) = SiOH4_new(ii) + f(mi)%sioh4*tmp2
                poc_new(ii) = poc_new(ii) + f(mi)%poc*tmp2
                smalg_new(ii) = smalg_new(ii) + alg_wc*dz

                ! keeping track of emerging new layer thickness
                dz_total = dz_total + dz

                ! assume layers that flood have a kick-ass brine channel network for flooding/desal
                ice(sc,ic,mi)%bced(ii) = 1

            elseif (mf_depth .gt. 0.) then
                mf_depth = mf_depth - th_new(ii)
                if (mf_depth .lt. 0.) then
                    dz = mf_depth + th_new(ii)
                else
                    dz = th_new(ii)
                endif

                ! record tracers
                z_snow = ice(sc,ic,mi)%snow%z
                tmp1 = 1. - ice(sc,ic,mi)%snow%d(z_snow)/IceD

                s_new(ii) = s_new(ii) + S_flood*dz
                d_new(ii) = d_new(ii) + d_flood*dz
                heat_total = heat_total + heat_flood*dz
                NO3_new(ii) = NO3_new(ii) + no3_pond_new*dz*tmp1
                NH4_new(ii) = NH4_new(ii) + nh4_pond_new*dz*tmp1
                PO4_new(ii) = PO4_new(ii) + po4_pond_new*dz*tmp1
                SiOH4_new(ii) = SiOH4_new(ii) + sioh4_pond_new*dz*tmp1
                poc_new(ii) = poc_new(ii) + poc_pond_new*dz*tmp1
                smalg_new(ii) = smalg_new(ii) + smalg_pond_new*dz*tmp1

                ! keeping track of emerging new layer thickness
                dz_total = dz_total + dz

            endif

            if (dz_total .lt. th_new(ii)) then

                ! find heights of new layer with regard to old thicknesses
                ! new_layer_top should be zero even if flooding has occured,
                ! because it is really a reference to the old ice pack.
                if (.NOT. remapping) then
                    remapping = .TRUE.
                    if (new_layer_top .lt. 0.) then
                        new_layer_top = 0.
                    endif
                        new_layer_bot = new_layer_top + th_new(ii) - dz_total
                else
                    new_layer_top = new_layer_bot
                    new_layer_bot = new_layer_bot + th_new(ii)
                endif

                ! find 1st old layer that contains part of new layer, going down
                jj = 0
                old_layer_bot = 0.
                do while(new_layer_top .ge. old_layer_bot .and. jj .le. sk_z)
                    jj=jj+1
                    if (jj .eq. 43) then
!                                            print *,'jj exceeded max layers in ice remapping'
!                                                 print *,'mi: ',mi,'layer: ',ii,'in sia2_env_ice_remap line 228'
                    endif
                    ! old thickness...
                    old_layer_bot = old_layer_bot + ice(sc,ic,mi)%th(jj)
                enddo
                old_layer_top = old_layer_bot - ice(sc,ic,mi)%th(jj)
                ! now jj is OLD layer where NEW layer ii starts...

                ! find total heat/salt from multiple layers/partial layers that make up
                ! the new layer
                do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. sk_z))


                    ! ----- NEW LAYER GEOMETRIES ------------------------------
                    interim_top = max(old_layer_top,new_layer_top)
                    interim_bot = min(old_layer_bot,new_layer_bot)

                    z1 = interim_top - old_layer_top
                    z2 = interim_bot - old_layer_top
                    dz = z2-z1

                    if (ii .le. int_z_new) then
                        if (jj .le. int_z .and. ii .le. int_z_new) then
                        ! redistribution of congelation ice to congelation ice (basal growth)
                            tmp1 = dz
                            s_mean = ice(sc,ic,mi)%s(jj)
                            d_mean = ice(sc,ic,mi)%d(jj)
                            heat_mean = ice(sc,ic,mi)%heat(jj)
                        else
                        ! redistribution of skeletal ice into congelation ice (basal growth)
                            tmp1 = min(1.,c_gl_sal/ice(sc,ic,mi)%s(jj))*dz ! remove brine volume
                            t_mean = f(mi)%t
                            s_mean = c_gl_sal
#include "sia2_env_brine_calc.inc.f90"
                        endif
                    else
                        ! redistribution of congelation ice to skeletal ice (basal melt)
                        tmp1 = dz
                        t_mean = f(mi)%t
                        s_mean = f(mi)%s/2.d0
#include "sia2_env_brine_calc.inc.f90"
                    endif

                    ! record tracers into new layer
                    s_new(ii) = s_new(ii) + s_mean*dz
                    d_new(ii) = d_new(ii) + d_mean*dz
                    heat_total = heat_total + heat_mean*dz
                    no3_new(ii) = no3_new(ii) + ice(sc,ic,mi)%no3(jj)*tmp1
                    nh4_new(ii) = nh4_new(ii) + ice(sc,ic,mi)%nh4(jj)*tmp1
                    po4_new(ii) = po4_new(ii) + ice(sc,ic,mi)%po4(jj)*tmp1
                    sioh4_new(ii) = sioh4_new(ii) + ice(sc,ic,mi)%sioh4(jj)*tmp1
                    poc_new(ii) = poc_new(ii) + ice(sc,ic,mi)%poc(jj)*tmp1
                    ! assuming algae doesn't dilute(is not lost) during skeletal->congelation freeze
                    !smalg_new(ii) = smalg_new(ii) + smalg_pre_map(jj)*dz
                    smalg_new(ii) = smalg_new(ii) + ice(sc,ic,mi)%smalg(jj)*dz

                    ! setup variable for next partial layer addition
                    dz_total = dz_total + dz

                    jj=jj+1
                    if (jj .le. sk_z) then
                        ! find boundaries of next old layer
                        old_layer_top = old_layer_bot
                        old_layer_bot = old_layer_bot + ice(sc,ic,mi)%th(jj)
                    endif

                enddo

                if (dz_total .lt. th_new(ii) .and. (jj .gt. sk_z)) then
                ! making new ice out of seawater

                    if (ii .le. int_z_new .and. c_gl_copy .gt. 0.d0) then ! these tests should alway be true together
                    ! new congelation ice out of seawater - fast growth!!
                        !dz = min(c_gl_copy,th_new(ii)-dz_total)
                        dz = th_new(ii)-dz_total
                        t_mean = f(mi)%t
                        s_mean = c_gl_sal
                    else
                        ! new skeletal ice out of seawater
                        dz = th_new(ii) - dz_total
                        s_mean = f(mi)%s*c_5
                        t_mean = min(f(mi)%t,f(mi)%s*mu)
                    endif
#include "sia2_env_brine_calc.inc.f90"
                    ! fraction seawater concentrations by new brine volume
                    tmp1 = dz*bv_mean*c_001

                    ! record tracers into new layer
                    s_new(ii) = s_new(ii) + s_mean*dz
                    d_new(ii) = d_new(ii) + d_mean*dz
                    heat_total = heat_total + heat_mean*dz
!                    no3_new(ii) = no3_new(ii) + f(mi)%no3*tmp1 ! extending into mixed layer
!                    nh4_new(ii) = nh4_new(ii) + f(mi)%nh4*tmp1 ! extending into mixed layer
!                    po4_new(ii) = po4_new(ii) + f(mi)%po4*tmp1 ! extending into mixed layer
!                    sioh4_new(ii) = sioh4_new(ii) + f(mi)%sioh4*tmp1 ! extending into mixed layer
                    no3_new(ii) = no3_new(ii) + ice(sc,ic,mi)%no3(ml_z)*tmp1 ! extending into mixed layer
                    nh4_new(ii) = nh4_new(ii) + ice(sc,ic,mi)%nh4(ml_z)*tmp1 ! extending into mixed layer
                    po4_new(ii) = po4_new(ii) + ice(sc,ic,mi)%po4(ml_z)*tmp1 ! extending into mixed layer
                    sioh4_new(ii) = sioh4_new(ii) + ice(sc,ic,mi)%sioh4(ml_z)*tmp1 ! extending into mixed layer
                    smalg_new(ii) = smalg_new(ii) + alg_wc*tmp1
                    poc_new(ii) = poc_new(ii) + f(mi)%poc*tmp1

                    ! new method for algal migration - using smalg_pre_map
                    smalg_pre_map(ii) = smalg_pre_map(ii) + alg_wc*tmp1/th_new(ii) ! convert

                    ! update boundaries for next layer
                    c_gl_copy = c_gl_copy - dz
                    dz_total = th_new(ii) ! round out layer - this is the last way to add mass

                endif

                ! round out dz_total whether we are shrinking or growing,
                ! otherwise division below can create/take away tracers
                dz_total = th_new(ii)

            endif     ! end of "if dz_total .lt. th_new(ii)

            if (ii .le. int_z_new) then
                ! solve for new layer temp. from heat/salt/density
                dz_mean = dz_total
                if (dz_total .eq. 0.) then
                    print *,mi,ic,i_ease,j_ease,'**********************************************'
                    print *,'t---',ice(sc,ic,mi)%t
                    print *,'t_prev---',t_prev
                    print *,'t_last---',t_last
                    print *,'t_next---',t_next
                    print *,'ki---',ki
                endif
                d_new(ii) = d_new(ii)/dz_total
                s_new(ii) = s_new(ii)/dz_total
                d_mean = d_new(ii)
                s_mean = s_new(ii)
                heat_mean = heat_total ! division by dz is included in calc below
                if (s_mean .le. 0.) then
                    print *,'Salinity (smean) is less than zero.'
                    print *,'mi: ',mi,'layer: ',iii,'in sia2_env_ice_remap line 363'
                endif

#include "sia2_env_temp_from_heat.inc.f90"

                ! check to make sure we don't get warmer than the ocean
                !if (ice(sc,ic,mi)%id(ii) .le. ice(sc,ic,mi)%fbh) then
                !        t_mean = min(t_mean,S_mean*mu)
                !else
                !        t_mean = min(t_mean,f(mi)%t)
                !endif

                t_new(ii) = t_mean

            else
                ! we know skeletal temp, save computations
                t_new(ii) = f(mi)%t
                s_new(ii) = s_new(ii)/dz_total
                d_new(ii) = d_new(ii)/dz_total
            endif

        enddo ! end of new grid portioning

        ! update layer indices
        int_z = int_z_new
        sk_1 = int_z + 1
        sk_z = int_z + z_sk
        ice(sc,ic,mi)%z = sk_z

        ! update ice structure
        do ii=1,sk_z

            ! assign new layer boundaries
            ice(sc,ic,mi)%th(ii) = th_new(ii)
            ice(sc,ic,mi)%id(ii) = id_new(ii)

            ! assign new ice physics
            ice(sc,ic,mi)%t(ii) = t_new(ii)
            ice(sc,ic,mi)%s(ii) = s_new(ii)
            t_mean = ice(sc,ic,mi)%t(ii)
            s_mean = max(min_sal,ice(sc,ic,mi)%s(ii))

#include "sia2_env_brine_calc.inc.f90"

            ice(sc,ic,mi)%d(ii) = d_mean
            ice(sc,ic,mi)%bs(ii) = bs_mean
            ice(sc,ic,mi)%bd(ii) = bd_mean
            ice(sc,ic,mi)%bv(ii) = bv_mean
            ice(sc,ic,mi)%heat(ii) = heat_mean


            ! update tracers
            ice(sc,ic,mi)%no3(ii) = NO3_new(ii)/th_new(ii)
            ice(sc,ic,mi)%nh4(ii) = NH4_new(ii)/th_new(ii)
            ice(sc,ic,mi)%po4(ii) = PO4_new(ii)/th_new(ii)
            ice(sc,ic,mi)%sioh4(ii) = SiOH4_new(ii)/th_new(ii)
            ice(sc,ic,mi)%poc(ii) = poc_new(ii)/th_new(ii)
            ice(sc,ic,mi)%smalg(ii) = smalg_new(ii)/th_new(ii) !

        enddo


        ! new method for algal migration!
        !print *,ice(sc,ic,mi)%smalg
        !print *, smalg_pre_map
        ice(sc,ic,mi)%smalg = smalg_pre_map

        dz_sk = 0.
#include "sia2_env_update_skeletal.inc.f90"

        ! adjust surface temp if there was no snow, but flooding occured (to prevent numerical instability)
        ! uh, this can't happen...
        !if (flooded .gt. 0 .and. m(mi)%z_snow .eq. 0) then
        !         ice(sc,ic,mi)%ts = ice(sc,ic,mi)%t(1)    ! make equal to 1st layer temp
        !endif

    endif     ! end of "don't grow if depth will be maxed out by growth"

!             if(ice(sc,ic,mi)%s(1) .le. 0.12 .and. ice(sc,ic,mi)%no3(1) .gt. 10.) then
!                     testvar = -1
!                     testvar = sqrt(testvar)
!             endif


            ! platelet crap
            ! -------------------------------------------------------------------
!             if (c_gl_copy .gt. 0. .and. use_pl .eq. 1) then
!                    ! move platelet layer down
!                    tmp1 = c_gl_copy/pl_th
!                    tmp2 = 1. - tmp1
!                    do ii = 1,pl_max-1
!                            pl%no3(ii) = pl%no3(ii)*tmp2 + pl%no3(ii+1)*tmp1
!                            pl%nh4(ii) = pl%nh4(ii)*tmp2 + pl%nh4(ii+1)*tmp1
!                            pl%po4(ii) = pl%po4(ii)*tmp2 + pl%po4(ii+1)*tmp1
!                            pl%sioh4(ii) = pl%sioh4(ii)*tmp2 + pl%sioh4(ii+1)*tmp1
!                            pl%smalg(ii) = pl%smalg(ii)*tmp2 + pl%smalg(ii+1)*tmp1
!                            pl%poc(ii) = pl%poc(ii)*tmp2 + pl%poc(ii+1)*tmp1
!                    enddo
!         pl%no3(pl_max) = pl%no3(pl_max)*tmp2 + f(mi)%no3*tmp1
!         pl%nh4(pl_max) = pl%nh4(pl_max)*tmp2 + f(mi)%nh4*tmp1
!         pl%po4(pl_max) = pl%po4(pl_max)*tmp2 + f(mi)%po4*tmp1
!         pl%sioh4(pl_max) = pl%sioh4(pl_max)*tmp2 + f(mi)%sioh4*tmp1
!         pl%smalg(pl_max) = pl%smalg(pl_max)*tmp2 + min_alg*tmp1
!         pl%poc(pl_max) = pl%poc(pl_max)*tmp2 + f(mi)%poc*tmp1
!             endif

! ----------------------------------------------------------------------
! End: Adjust Boundaries
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
