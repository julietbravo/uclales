!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999, 2001, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module srfc

use lsmdata

  real    :: zrough  = 0.01
  real    :: ubmin   = 0.20
  real    :: dthcon  = 100.0
  real    :: drtcon  = 0.0

  logical :: lhomflx = .false.
  real    :: rh_srf  = 1.
  real    :: drag    = -1.

  ! BvS: use local surface temp in radiation calculation
  ! instead of global average (SST)
  logical :: lrad_ts_loc = .false.  

  ! --------------------------
  ! some real-time statistics
  real    :: Hg,LEg,Gg,oblg,ustarg

contains
  !
  ! --------------------------------------------------------------------------
  ! SURFACE: Calculates surface fluxes using an algorithm chosen by ISFCTYP
  ! and fills the appropriate 2D arrays
  !
  !   default: specified thermo-fluxes (drtcon, dthcon)
  !   isfctyp=1: specified surface layer gradients (drtcon, dthcon)
  !   isfctyp=2: fixed lower boundary of water at certain sst
  !   isfctyp=3: bulk aerodynamic law with coefficients (drtcon, dthcon)
  !   isfctyp=4: regulate surface temperature to yield a constant surface buoyancy flux
  !   isfctyp=5: surface temperature and humidity determined using LSM (van Heerwaarden)
  !
  subroutine surface(sst,time_in)

    use defs, only: vonk, p00, rcp, cpr, g, cp, alvl, ep2
    use grid, only: nzp, nxp, nyp, a_up, a_vp, a_theta, vapor, zt, psrf,   &
         th00, umean, vmean, dn0, level, a_ustar, a_tstar, a_rstar,        &
         uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc, nstep, a_tskin, a_qskin,  &
         isfctyp, a_phiw, a_tsoil, a_Wl, obl, &
         a_lflxu, a_lflxd, a_sflxu, a_sflxd, rkalpha, rkbeta, a_Qnet, dt, a_G0, wspd, pi0, pi1
    use thrm, only: rslf
    use stat, only: sfc_stat, sflg
    use util, only : get_avg3
    use mpi_interface, only : nypg, nxpg, double_array_par_sum, double_scalar_par_sum, myid

    implicit none

    real, optional, intent (inout) :: sst, time_in
    integer :: i, j, k, iterate 
    real :: zs, bflx, ffact, sst1, bflx1, Vbulk, Vzt, usum, tskinsum_loc
    real (kind=8) :: bfl(2), bfg(2)

    real :: dtdz(nxp,nyp), drdz(nxp,nyp), usfc(nxp,nyp), vsfc(nxp,nyp) &
            ,bfct(nxp,nyp), mnflx(5), flxarr(5,nxp,nyp)
            !,wspd(nxp,nyp), bfct(nxp,nyp), mnflx(5), flxarr(5,nxp,nyp)

    ! BvS: for isfctyp=55 
    !---------------------------------------
    real     :: tsbar ! Average surface T
    real     :: rhcpa,gkm,gkp,soiltend,temp(20),zeff,exnera,exners

    drdz(:,:)   = 0.

    select case(isfctyp)
    !
    ! use prescribed surface gradients dthcon, drton from NAMELIST
    ! use then similarity theory to compute the fluxes 
    !
    case(1)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       do j=3,nyp-2
         do i=3,nxp-2
           dtdz(i,j) = dthcon
           if(level>0) drdz(i,j)=drtcon
         end do
       end do
       zs = zrough
  
       if(level>0) then
         call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,a_ustar,a_tstar,obl,drt=drdz,rstar=a_rstar)
         call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,  &
                       uw_sfc,vw_sfc,wt_sfc,ww_sfc,wq_sfc,a_rstar)
       else
         call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,a_ustar,a_tstar,obl)
         call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,  &
                       uw_sfc,vw_sfc,wt_sfc,ww_sfc)
       end if

    !
    ! ----------------------------------------------------------------------
    ! use prescribed SST and assume qsurf=qsat (i.e. ocean) to compute
    ! gradients. Then use similarity theory to predict the fluxes. 
    !
    case(2)

       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       usum = 0.
       do j=3,nyp-2
         do i=3,nxp-2
           dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
           if(level>0) drdz(i,j) = vapor(2,i,j) - rslf(psrf,sst)
           bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
           usum = usum + a_ustar(i,j)
         end do
       end do
       usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
       zs = zrough
       if (zrough <= 0.) zs = max(0.0001,(0.016/g)*usum**2)

       if(level>0) then
         call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,a_ustar,a_tstar,obl,drt=drdz,rstar=a_rstar)
         call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,  &
                       uw_sfc,vw_sfc,wt_sfc,ww_sfc,wq_sfc,a_rstar)
       else
         call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,a_ustar,a_tstar,obl)
         call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,  &
                       uw_sfc,vw_sfc,wt_sfc,ww_sfc)
       end if

    !
    ! ----------------------------------------------------------------------
    ! drtcon (wq=Ch*u*dth, Garrat p.55) and using prescribed sst 
    ! and qsurf=qsat (ocean); note that here zrough is not the roughness
    !
    case(3)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       do j=3,nyp-2
         do i=3,nxp-2
           dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
           if(level>0) drdz(i,j) = vapor(2,i,j) - rh_srf*rslf(psrf,sst)

           if (ubmin > 0.) then
             a_ustar(i,j) = sqrt(zrough)* wspd(i,j)
           else
             a_ustar(i,j) = abs(ubmin)
           end if

           if (drag > 0) then
             a_tstar(i,j) =  drag*dtdz(i,j)/a_ustar(i,j)
             if(level>0) a_rstar(i,j) =  drag*drdz(i,j)/a_ustar(i,j)
           else
             a_tstar(i,j) =  dthcon * wspd(i,j)*dtdz(i,j)/a_ustar(i,j)
             if(level>0) a_rstar(i,j) =  drtcon * wspd(i,j)*drdz(i,j)/a_ustar(i,j)
           endif

           bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
         end do
       end do


       if(level>0) then
         call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,  &
                       uw_sfc,vw_sfc,wt_sfc,ww_sfc,wq_sfc,a_rstar)
       else
         call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,  &
                       uw_sfc,vw_sfc,wt_sfc,ww_sfc)
       end if

       if (lhomflx) then
         flxarr(1,:,:) = uw_sfc
         flxarr(2,:,:) = vw_sfc
         flxarr(3,:,:) = ww_sfc
         flxarr(4,:,:) = wt_sfc
         if(level>0) then
           flxarr(5,:,:) = wq_sfc
           call get_avg3(5,nxp,nyp,flxarr,mnflx)
         else
           call get_avg3(4,nxp,nyp,flxarr,mnflx)
         end if 
         uw_sfc = mnflx(1)
         vw_sfc = mnflx(2)
         ww_sfc = mnflx(3)
         wt_sfc = mnflx(4)
         if(level>0) wq_sfc = mnflx(5)
       end if

    !
    ! ----------------------------------------------------------------------
    ! fix surface temperature to yield a constant surface buoyancy flux
    ! dthcon
    !
    case(4)
       Vzt   = 10.* (log(zt(2)/zrough)/log(10./zrough))
       Vbulk = Vzt * (vonk/log(zt(2)/zrough))**2

       bfl(:) = 0.
       do j=3,nyp-2
         do i=3,nxp-2
           bfl(1) = bfl(1)+a_theta(2,i,j)
           bfl(2) = bfl(2)+vapor(2,i,j)
         end do
       end do

       call double_array_par_sum(bfl,bfg,2)

       bfg(2) = bfg(2)/real((nxpg-4)*(nypg-4))
       bfg(1) = bfg(1)/real((nxpg-4)*(nypg-4))

       do iterate=1,5
         bflx  = ((sst -bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst) -bfg(2))) &
               * 0.5*(dn0(1)+dn0(2))*cp*Vbulk
         sst1 = sst + 0.1
         bflx1 = ((sst1-bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst1)-bfg(2))) &
               * 0.5*(dn0(1)+dn0(2))*cp*Vbulk
         sst  = sst + 0.1* (dthcon - bflx) / (bflx1-bflx)
       end do

       do j=3,nyp-2
         do i=3,nxp-2
           wt_sfc(i,j) = Vbulk * (sst -a_theta(2,i,j))
           wq_sfc(i,j) = Vbulk * (rslf(psrf,sst) - vapor(2,i,j))
           wspd(i,j)    = max(0.1,                                    &
                          sqrt((a_up(2,i,j)+umean)**2+(a_vp(2,i,j)+vmean)**2))
           bflx         = wt_sfc(i,j)*g/bfg(1) + g*ep2*wq_sfc(i,j)
           a_ustar(i,j) = diag_ustar(zt(2),zrough,bflx,wspd(i,j))
           uw_sfc(i,j)  = -a_ustar(i,j)*a_ustar(i,j)                  &
                *(a_up(2,i,j)+umean)/wspd(i,j)
           vw_sfc(i,j)  = -a_ustar(i,j)*a_ustar(i,j)                  &
                *(a_vp(2,i,j)+vmean)/wspd(i,j)
           ww_sfc(i,j)  = 0.
           a_rstar(i,j) = wq_sfc(i,j)/a_ustar(i,j)
           a_tstar(i,j) = wt_sfc(i,j)/a_ustar(i,j)
          end do
       end do  

    !
    ! ----------------------------------------------------------------------
    ! Malte: Get surface fluxes using a land surface model (van Heerwaarden, DALES)
    !
    !case(5)
    !   !Initialize Land Surface -> BvS called from init.f90 
    !   if (init_lsm) then
    !     call initlsm(sst,time_in)
    !     init_lsm = .false.
    !   end if

    !   !Local or filtered variables for flux calculation
    !   if (local) then
    !     !Here we are assuming that MO theory is also valid locally
    !     u0bar(:,:,:) = a_up(:,:,:)
    !     v0bar(:,:,:) = a_vp(:,:,:)
    !     thetaav(:,:) = a_theta(2,:,:)
    !     vaporav(:,:) = vapor(2,:,:)
    !     tskinav(:,:) = a_tskin(:,:)
    !     qskinav(:,:) = a_qskin(:,:)
    !   else
    !     stop('filtering not supported in LSM')
    !   end if

    !   !Calculate surface wind for flux calculation
    !   call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,u0bar,v0bar,umean,vmean)

    !   ! a) Calculate Monin Obuhkov Length from surface scalars
    !   do j=3,nyp-2
    !     do i=3,nxp-2
    !       dtdz(i,j) = thetaav(i,j) - tskinav(i,j)
    !       drdz(i,j) = vaporav(i,j) - qskinav(i,j)
    !     end do
    !   end do

    !   tskinavg = sum(a_tskin(3:(nxp-2),3:(nyp-2)))/(nxp-4)/(nyp-4)
    !   sst      = tskinavg*(psrf/p00)**(rcp)
    !   call srfcscls(nxp,nyp,zt(2),zrough,tskinavg,wspd,dtdz,a_ustar,a_tstar,obl,drt=drdz,rstar=a_rstar)

    !   !Calculate the drag coefficients and aerodynamic resistance
    !   do j=3,nyp-2
    !     do i=3,nxp-2
    !       cm(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
    !                  obl(i,j)) + psim(z0m(i,j) / obl(i,j))) ** 2.
    !       cs(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
    !                  obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zt(2) / &
    !                  z0h(i,j)) - psih(zt(2) / obl(i,j)) + psih(z0h(i,j) / &
    !                  obl(i,j)))
    !       ra(i,j) =  1. / (cs(i,j)* wspd(i,j))
    !     end do
    !   end do

    !   !Get skin temperature and humidity from land surface model (van Heerwaarden)
    !   call lsm

    !   !Update tskin and qskin (local or filtered) for flux calculation
    !   if (local) then
    !     tskinav(:,:) = a_tskin(:,:)
    !     qskinav(:,:) = a_qskin(:,:)
    !   else
    !     !Insert filter method here
    !   end if

    !   !Calculate the surface fluxes with bulk law (Fairall, 2003)
    !   !Fluxes in kinematic form with units [K m/s] and [kg/kg m/s]
    !   do j=3, nyp-2
    !     do i=3, nxp-2
    !       wt_sfc(i,j)  = - (thetaav(i,j) - tskinav(i,j)) / ra(i,j) 
    !       wq_sfc(i,j)  = - (vaporav(i,j) - qskinav(i,j)) / ra(i,j)
    !       uw_sfc(i,j)  = - a_ustar(i,j)*a_ustar(i,j)                  &
    !                        *(a_up(2,i,j)+umean)/wspd(i,j)
    !       vw_sfc(i,j)  = - a_ustar(i,j)*a_ustar(i,j)                  &
    !                        *(a_vp(2,i,j)+vmean)/wspd(i,j)

    !       a_rstar(i,j) = - wq_sfc(i,j)/a_ustar(i,j)
    !       a_tstar(i,j) = - wt_sfc(i,j)/a_ustar(i,j)
    !     end do
    !   end do

    !
    ! ----------------------------------------------------------------------
    ! BvS; Malte's LSM, with some changes. getobl() works better than the sfcscls() 
    !      method for stable conditions and the transition from NBL-CBL-NBL
    !
    case(5)
       if (init_lsm) then
         call initlsm(sst,time_in)
         init_lsm = .false.
       end if

       ! Temp, for getobl(). To-do: fix :)
       u0bar(:,:,:) = a_up(:,:,:)
       v0bar(:,:,:) = a_vp(:,:,:)
       thetaav(:,:) = a_theta(2,:,:)
       vaporav(:,:) = vapor(2,:,:)
       tskinav(:,:) = a_tskin(:,:)
       qskinav(:,:) = a_qskin(:,:)

       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       call getobl(wspd)

       ! Get drag coefficients and aerodynamic resistance
       do j=3,nyp-2
         do i=3,nxp-2
           cm(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
                      obl(i,j)) + psim(z0m(i,j) / obl(i,j))) ** 2.
           cs(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
                      obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zt(2) / &
                      z0h(i,j)) - psih(zt(2) / obl(i,j)) + psih(z0h(i,j) / obl(i,j)))
           ra(i,j) =  1. / (cs(i,j)* wspd(i,j))
         end do
       end do

       !tskinavg = sum(a_tskin(3:(nxp-2),3:(nyp-2)))/(nxp-4)/(nyp-4)
       ! BvS: calculate global averaged skin temperature
       tskinsum_loc = sum(a_tskin(3:(nxp-2),3:(nyp-2)))
       call double_scalar_par_sum(tskinsum_loc, tskinavg)
       tskinavg  = tskinavg / ((nxpg-4)*(nypg-4)) 

       ! Set SST to domain averaged skin temperature for the radiation code
       sst      = tskinavg*(psrf/p00)**(rcp)

       ! Call the LSM
       call lsm

       ! Calculate surface fluxes
       do j=3, nyp-2
         do i=3, nxp-2
           a_ustar(i,j) = vonk * wspd(i,j)  / (log(zt(2) / z0m(i,j)) - psim(zt(2) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))
           wt_sfc(i,j)  = - (a_theta(2,i,j) - a_tskin(i,j)) / ra(i,j)
           wq_sfc(i,j)  = - (vaporav(i,j) - qskinav(i,j))   / ra(i,j)
           uw_sfc(i,j)  = - a_ustar(i,j) * a_ustar(i,j) * (a_up(2,i,j)+umean)/wspd(i,j)
           vw_sfc(i,j)  = - a_ustar(i,j) * a_ustar(i,j) * (a_vp(2,i,j)+vmean)/wspd(i,j)
           a_tstar(i,j) = - wt_sfc(i,j)/a_ustar(i,j)
           a_rstar(i,j) = - wq_sfc(i,j)/a_ustar(i,j)
         end do
       end do

    !
    ! ----------------------------------------------------------------------
    ! Bart: simplified LSM for case without moisture
    !
    case(55)
       if(init_lsm) call initlsm_simple   ! shouldn't be necessary

       exnera = ((p00*((pi0(2)+pi1(2))/cp)**cpr)/p00)**rcp
       exners = (psrf/p00)**rcp

       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)

       ! Option to filter wind speed
       select case(imostloc)
         case(1)   ! Use local filtered wind
           call filterwind(nxp,nyp,wspd,1)
         case(2)   ! Use global average of wind
           call filterwind(nxp,nyp,wspd,2)
         case default
           ! do nothing
       end select

       call getobl(wspd)

       ! Get drag coefficients and aerodynamic resistance
       do j=3,nyp-2
         do i=3,nxp-2
           cm(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
                      obl(i,j)) + psim(z0m(i,j) / obl(i,j))) ** 2.
           cs(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
                      obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zt(2) / &
                      z0h(i,j)) - psih(zt(2) / obl(i,j)) + psih(z0h(i,j) / obl(i,j)))
           ra(i,j) =  1. / (cs(i,j)* wspd(i,j))
         end do
       end do

       if (dolsm) then ! --------------------------------------------------
         ! Solve surface temperature from SEB
         do j=3,nyp-2
           do i=3,nxp-2
             a_Qnet(i,j)  = a_lflxd(2,i,j) - a_lflxu(2,i,j) + a_sflxd(2,i,j) - a_sflxu(2,i,j)
             rhcpa        = dn0(2)*cp/ra(i,j)
             a_tskin(i,j) = (a_Qnet(i,j) + rhcpa * a_theta(2,i,j)*exnera + labsk * a_tsoil(1,i,j)) / (rhcpa + labsk)
           end do
         end do

         ! Integrate soil heat diffusion
         tsoilm(:,:,:) = a_tsoil(:,:,:)  
         do j=3, nyp-2
           do i=3, nxp-2
             do k=1,3 
               if(k==1) then
                 gkm = labsk * (a_tskin(i,j) - tsoilm(k,i,j))
               else
                 gkm = -lambdab * ((tsoilm(k,i,j) - tsoilm(k-1,i,j)) / (zsoilc(k) - zsoilc(k-1)))
               end if
               gkp   = -lambdab * ((tsoilm(k+1,i,j) - tsoilm(k,i,j))   / (zsoilc(k+1) - zsoilc(k)))
               soiltend = -(gkp - gkm) / (zsoil(k+1) - zsoil(k)) 

               if(nstep == 2) then
                 a_tsoil(k,i,j) = a_tsoil(k,i,j) + rkbeta(nstep)  * (dt/rhoCs) * soiltend &
                                                 + rkalpha(nstep) * (dt/rhoCs) * soiltendm(k,i,j)
               else
                 a_tsoil(k,i,j) = a_tsoil(k,i,j) + rkalpha(nstep) * (dt/rhoCs) * soiltend &
                                                 + rkbeta(nstep)  * (dt/rhoCs) * soiltendm(k,i,j)
               end if

               if(nstep==3) then
                 soiltendm(k,i,j) = 0.
               else                              
                 soiltendm(k,i,j) = soiltend
               end if
             end do
           end do
         end do
       else ! --------------------------------------------------
         a_tskin(:,:) = sst * (p00/psrf)**rcp
       end if

       a_tskin = a_tskin / exners ! to potential for flux and obl calculation in next t-step
 
       ! Calculate surface fluxes
       do j=3, nyp-2
         do i=3, nxp-2
           a_ustar(i,j) = vonk * wspd(i,j)  / (log(zt(2) / z0m(i,j)) - psim(zt(2) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))
           wt_sfc(i,j)  = - (a_theta(2,i,j) - a_tskin(i,j)) / ra(i,j)
           uw_sfc(i,j)  = - a_ustar(i,j)*a_ustar(i,j) * (a_up(2,i,j)+umean)/wspd(i,j)
           vw_sfc(i,j)  = - a_ustar(i,j)*a_ustar(i,j) * (a_vp(2,i,j)+vmean)/wspd(i,j)
           a_tstar(i,j) = - wt_sfc(i,j)/a_ustar(i,j)
           if(dolsm) a_G0(i,j) = labsk * (a_tskin(i,j)*exners - a_tsoil(1,i,j))
         end do
       end do

    !
    ! ----------------------------------------------------------------------
    ! fix thermodynamic fluxes at surface given values in energetic 
    ! units and calculate  momentum fluxes from similarity theory
    !
    case default
       ffact = 1.
       wt_sfc(1,1) = ffact* dthcon / (0.5*(dn0(1)+dn0(2)) * cp)
       if(level>0) wq_sfc(1,1)  = ffact* drtcon/(0.5*(dn0(1)+dn0(2))*alvl)

       if (zrough <= 0.) then
         usum = 0.
         do j=3,nyp-2
           do i=3,nxp-2
             usum = usum + a_ustar(i,j)

           end do
         end do
         usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
         zs = max(0.0001,(0.016/g)*usum**2) !Charnock for flow over sea
       else
         zs = zrough
       end if

       do j=3,nyp-2
         do i=3,nxp-2
           wt_sfc(i,j) = wt_sfc(1,1) 
           if(level>0) wq_sfc(i,j) = wq_sfc(1,1)

           wspd(i,j)    = max(0.1,                                    &
                          sqrt((a_up(2,i,j)+umean)**2+(a_vp(2,i,j)+vmean)**2))

           if (ubmin > 0.) then
             bflx = g*wt_sfc(1,1)/th00
           if (level >= 2) bflx = bflx + g*ep2*wq_sfc(i,j)
             a_ustar(i,j) = diag_ustar(zt(2),zs,bflx,wspd(i,j))
           else
             a_ustar(i,j) = abs(ubmin)
           end if

           ffact = a_ustar(i,j)*a_ustar(i,j)/wspd(i,j)
           uw_sfc(i,j)  = -ffact*(a_up(2,i,j)+umean)
           vw_sfc(i,j)  = -ffact*(a_vp(2,i,j)+vmean)
           ww_sfc(i,j)  = 0.
           if(level>0) a_rstar(i,j) = -wq_sfc(i,j)/a_ustar(i,j)
           a_tstar(i,j) = -wt_sfc(i,j)/a_ustar(i,j)
         end do
       end do

    end select

    if (sflg) then
      if(level==0) then
        call sfc_stat(nxp, nyp, wt_sfc, a_ustar, sst)
      else
        call sfc_stat(nxp, nyp, wt_sfc, a_ustar, sst, wq_sfc)
      end if
    end if

    if(isfctyp >= 5) call srfcstat

    return

  end subroutine surface

  !
  ! -------------------------------------------------------------------
  ! BvS: srfcstat: some global statistics of surface
  !
  subroutine srfcstat
    use mpi_interface, only : myid,double_scalar_par_sum,nxpg,nypg
    use grid, only          : wt_sfc, wq_sfc, obl, a_ustar, a_G0, nxp, nyp, dn0, level
    use defs, only          : cp, alvl
    implicit none

    real    :: Hl,Gl,LEl,obll,ustarl
    integer :: i,j

    Hl     = sum(wt_sfc (3:nxp-2,3:nyp-2))
    if(level>0) LEl = sum(wq_sfc(3:nxp-2,3:nyp-2))
    Gl     = sum(a_G0   (3:nxp-2,3:nyp-2))
    obll   = sum(obl    (3:nxp-2,3:nyp-2))
    ustarl = sum(a_ustar(3:nxp-2,3:nyp-2))
    
    call double_scalar_par_sum(Hl,Hg)
    if(level>0) call double_scalar_par_sum(LEl,LEg)
    call double_scalar_par_sum(Gl,Gg)
    call double_scalar_par_sum(obll,oblg)
    call double_scalar_par_sum(ustarl,ustarg)
    Hg     = (Hg     / ((nxpg-4)*(nypg-4))) * cp * (dn0(1)+dn0(2))/2.
    Gg     = (Gg     / ((nxpg-4)*(nypg-4))) 
    oblg   = (oblg   / ((nxpg-4)*(nypg-4))) 
    ustarg = (ustarg / ((nxpg-4)*(nypg-4)))

    if(level>0) then
      LEg    = (LEg    / ((nxpg-4)*(nypg-4))) * alvl * (dn0(1)+dn0(2))/2.
    else
      LEg    = 0.
    end if

  end subroutine srfcstat

  !
  ! -------------------------------------------------------------------
  ! BvS: filterwind: calculate local or global average of wind
  !
  subroutine filterwind(n2,n3,wind,ftype)
    use mpi_interface, only : myid,double_scalar_par_sum,nxpg,nypg
    implicit none

    integer, intent(in) :: n2,n3
    real, intent(inout) :: wind(n2,n3)
    integer, intent(in) :: ftype
    real                :: wind0(n2,n3)

    real                :: uloc,uglob
    integer             :: i,j
    integer, parameter  :: nfc = 1   ! BvS: can only be one because of nghost!!
    real, parameter     :: fac = 1. / (1+nfc+nfc)**2.

    select case(ftype)
      case(1)   ! Filter wind locally (top hat)
        wind0 = wind
        do j=3,n3-2
           do i=3,n2-2
             wind(i,j) = sum(wind0(i-nfc:i+nfc,j-nfc:j+nfc)) * fac
           end do
        end do 

      case(2)   ! Average wind globally
        uloc     = sum(wind(3:n2-2,3:n3-2))
        call double_scalar_par_sum(uloc,uglob)
        uglob    = (uglob / ((nxpg-4)*(nypg-4))) 
        
        do j=3,n3-2
           do i=3,n2-2
             wind(i,j) = uglob
           end do
        end do 

      case default
        stop 'filter type not supported'
    end select

  end subroutine filterwind

  !
  ! -------------------------------------------------------------------
  ! GET_SWNDS: returns surface winds valid at cell centers
  !
  subroutine get_swnds(n1,n2,n3,usfc,vsfc,wspd,up,vp,umean,vmean)

    implicit none

    integer, intent (in) :: n1, n2, n3
    real, intent (in)    :: up(n1,n2,n3), vp(n1,n2,n3), umean, vmean
    real, intent (out)   :: usfc(n2,n3), vsfc(n2,n3), wspd(n2,n3)

    integer :: i, j, ii, jj

    ! BvS: include one ghost cell
    do j=2,n3-1
       jj = j-1
       do i=2,n2-1
          ii = i-1
          usfc(i,j) = (up(2,i,j)+up(2,ii,j))*0.5+umean
          vsfc(i,j) = (vp(2,i,j)+vp(2,i,jj))*0.5+vmean
          wspd(i,j) = max(abs(ubmin),sqrt(usfc(i,j)**2+vsfc(i,j)**2))
       enddo
    enddo

  end subroutine get_swnds

  !
  ! ----------------------------------------------------------------------
  ! FUNCTION GET_USTAR:  returns value of ustar using the below
  ! similarity functions and a specified buoyancy flux (bflx) given in
  ! kinematic units
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !
  real function diag_ustar(z,z0,bflx,wnd)

    use defs, only: vonk

    implicit none

    real, parameter      :: am   =  4.7   !   "          "         "
    real, parameter      :: bm   = 16.0   !   "          "         "
    real, parameter      :: eps  = 1.e-10 ! non-zero, small number

    real, intent (in)    :: z             ! height where u locates
    real, intent (in)    :: z0            ! momentum roughness height
    real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
    real, intent (in)    :: wnd           ! wind speed at z

    integer :: iterate
    real    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

    lnz   = log(z/z0)
    klnz  = vonk/lnz
    c1    = 3.14159/2. - 3.*log(2.)

    ustar =  wnd*klnz
    if (bflx /= 0.0) then
       do iterate=1,4
          lmo   = -(ustar**3)/(bflx*vonk + eps)
          zeta  = z/lmo
          if (zeta > 0.) then
             ustar =  vonk*wnd  /(lnz + am*zeta)
          else
             x     = sqrt( sqrt( 1.0 - bm*zeta ) )
             psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
             ustar = wnd*vonk/(lnz - psi1)
          end if
       end do
    end if

    diag_ustar = ustar

    return
  end function diag_ustar

  !
  ! ----------------------------------------------------------------------
  ! Subroutine srfcscls:  returns scale values based on Businger/Dye
  ! similarity functions.
  !
  ! phi_h (zeta > 0) =  Pr * (1 + ah * zeta)
  ! phi_h (zeta < 0) =  Pr * (1 - bh * zeta)^(-1/2)
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !

  subroutine srfcscls(n2,n3,z,z0,th00,u,dth,ustar,tstar,obl,drt,rstar)
    use defs, only : vonk, g, ep2
    use grid ,only: nstep, runtype, level, a_theta, vapor, isfctyp
    use mpi_interface, only : myid

    implicit none
    real, parameter     :: pr   =  0.74  ! prandlt number
    real, parameter     :: eps  = 1.e-10 ! non-zero, small number

    integer, intent(in) :: n2,n3         ! span of indicies covering plane
    real, intent(in)    :: z             ! height where u & T locate
    real, intent(in)    :: z0            ! momentum roughness height
    real, intent(in)    :: th00          ! reference temperature
    real, intent(in)    :: u(n2,n3)      ! velocities at z
    real, intent(in)    :: dth(n2,n3)    ! theta (th(z) - th(z0))
    real, intent(inout) :: ustar(n2,n3)  ! scale velocity
    real, intent(inout) :: tstar(n2,n3)  ! scale temperature
    real, intent(inout) :: obl(n2,n3)    ! Obukhov Length
    real, intent(in),optional     :: drt(n2,n3)    ! qt(z) - qt(z0)
    real, intent(inout), optional :: rstar(n2,n3)  ! scale value of qt

    logical, save       :: first_call=.True.
    integer             :: i,j,iterate,iter
    real                :: lnz, klnz, betg
    real                :: zeta, lmo, dtv
    real                :: thv,Rib,Lstart,Lend,Lold,fx,fxdif,Ldif,zeff,cs0
    logical             :: exititer

    lnz   = log(z/z0)
    klnz  = vonk/lnz
    betg  = th00/g

    do j=3,n3-2
       do i=3,n2-2
          if(level>0) then
            dtv = dth(i,j) + ep2*th00*drt(i,j)
          else
            dtv = dth(i,j)
          end if

          !
          ! Neutral case
          ! 
          if (dtv == 0.) then
            ustar(i,j) =  vonk*u(i,j)/lnz
            tstar(i,j) =  vonk*dtv/(pr*lnz)
            lmo        = -1.e10

          !
          ! start iterations from values at previous tstep, 
          ! unless the sign has changed or if it is the first call, then 
          ! use neutral values.
          !
          else

             if ((runtype=='INITIAL' .and. first_call) .or.( tstar(i,j)*dtv <= 0.)) then
               ustar(i,j) =  vonk*u(i,j)/lnz
               tstar(i,j) =  vonk*dtv/(pr*lnz)
               lmo        = -1.e10
             end if

             if(ustar(i,j) == 0) ustar(i,j) = 0.1

             Lold  = 1e9
             Ldif  = 1e9
             iter  = 0
             exititer = .false.
             !do iterate = 1,100
             do while(abs(Ldif)>0.1)
               lmo        = betg*ustar(i,j)**2/(vonk*tstar(i,j))
               Ldif       = lmo - Lold
               Lold       = lmo

               if ((dtv < 0) .and. (lmo > -0.001)) lmo = -0.001999
               if ((dtv > 0) .and. (lmo < +0.001)) lmo = +0.001777

               ! BvS : Following ECMWF, limit z/L for very stable conditions
               !if(z/lmo > 5.) then
               !  zeff = lmo * 5.
               !  exititer = .true.
               !else
               !  zeff = z
               !end if

               !zeta       = zeff/lmo
               !ustar(i,j) = u(i,j)*vonk/(log(zeff/z0) - psim(zeta))
               !if(ustar(i,j)<0.) ustar(i,j) = 0.1
               !tstar(i,j) = (dtv*vonk/pr)/(log(zeff/z0) - psih(zeta))

               !if(exititer) then
               !  lmo        = zeff/5. 
               !  exit
               !end if

               if (isfctyp==5) then  ! Use local values for z0m, z0h
                 cs0 = vonk**2. / (log(z/z0m(i,j)) - psim(z/lmo) + psim(z0m(i,j)/lmo)) / &
                                  (log(z/z0h(i,j)) - psih(z/lmo) + psih(z0h(i,j)/lmo))
                 ustar(i,j) = vonk*u(i,j) / (log(z/z0m(i,j)) - psim(z/lmo) + psim(z0m(i,j)/lmo))
               else                  ! Use global values for z0
                 cs0 = vonk**2. / (log(z/z0) - psim(z/lmo) + psim(z0/lmo)) / &
                                  (log(z/z0) - psih(z/lmo) + psih(z0/lmo))
                 ustar(i,j) = vonk*u(i,j) / (log(z/z0) - psim(z/lmo) + psim(z0/lmo))
               end if
               tstar(i,j) = dtv * cs0 * u(i,j) / ustar(i,j)    

               iter = iter + 1

               ! Limit L for day/night transitions
               if(lmo > 1e6)  lmo = 1e6
               if(lmo < -1e6) lmo = -1e6 

               if(iter>10000) then
                 print*,'Obukh. length not converged, myid=',myid,'i,j=',i,j
                 stop
               end if
             end do
          end if
 
          obl(i,j) = lmo
          if(present(rstar)) rstar(i,j) = tstar(i,j)*drt(i,j)/(dtv + eps)
          tstar(i,j) = tstar(i,j)*dth(i,j)/(dtv + eps)

       end do
    end do

    first_call = .False.

    return
  end subroutine srfcscls


  !
  ! ----------------------------------------------------------------------
  ! subroutine: sfcflxs:  this routine returns the surface fluxes based
  ! on manton-cotton algebraic surface layer equations.
  !
  subroutine sfcflxs(n2,n3,vk,ubar,u,v,xx,us,ts,uw,vw,tw,ww,rw,rstar)
    implicit none
    real, parameter             :: cc=4.7,eps=1.e-20
    integer, intent(in)         :: n2,n3
    real, intent(in)            :: ubar(n2,n3),u(n2,n3),v(n2,n3),xx(n2,n3),vk
    real, intent(in)            :: us(n2,n3),ts(n2,n3)
    real, intent(out)           :: uw(n2,n3),vw(n2,n3),tw(n2,n3),ww(n2,n3)
    real, intent(in), optional  :: rstar(n2,n3)
    real, intent(out), optional :: rw(n2,n3)

    real    :: x(n2,n3),y(n2,n3)
    integer :: i,j

    do j=3,n3-2
       do i=3,n2-2
          uw(i,j)=-(u(i,j)/(ubar(i,j)+eps))*us(i,j)**2
          vw(i,j)=-(v(i,j)/(ubar(i,j)+eps))*us(i,j)**2
          tw(i,j)=-ts(i,j)*us(i,j)
          if(present(rw)) rw(i,j)=-rstar(i,j)*us(i,j)

          x(i,j) = xx(i,j)*vk*ts(i,j)*(ubar(i,j)/us(i,j))**2
          x(i,j) = x(i,j)*sqrt(sqrt(1.-15.*min(0.,x(i,j)))) &
               /(1.0+cc*max(0.,x(i,j)))
          y(i,j) =sqrt((1.-2.86*x(i,j))/(1.+x(i,j)* &
               (-5.39+x(i,j)*6.998 )))
          ww(i,j)=(0.27*max(6.25*(1.-x(i,j))*y(i,j),eps)-&
               1.18*x(i,j)*y(i,j))*us(i,j)**2
       enddo
    enddo
    return

  end subroutine sfcflxs

  !
  ! ----------------------------------------------------------------------
  ! subroutine: getobl - Calculates the Obuhkov length iteratively. (van Heerwarden)
  !
  subroutine getobl(wspd)
    use defs, only: g,ep2
    use grid, only: nzp,nxp,nyp,a_up,a_vp,a_theta,umean,vmean,vapor,zt,&
                    u0,v0,obl,level,a_theta,a_tskin,a_up,a_vp
    implicit none
    real, intent(in) :: wspd(nxp,nyp)
    integer          :: i,j,iter
    !real             :: upcu, vpcv
    real             :: L, Lstart, Lend, Lold
    real             :: Rib, fx, fxdif !, wspd2
    real             :: thetavbar, tvskinbar
     
    do j=3,nyp-2
      do i=3,nxp-2
        if(level==0) then
          thetavbar   = a_theta(2,i,j)
          tvskinbar   = a_tskin(i,j)
        else
          !thetavbar   = thetaav(i,j) * (1. + ep2 * vaporav(i,j))
          !tvskinbar   = tskinav(i,j) * (1. + ep2 * qskinav(i,j))
          thetavbar   = thetaav(i,j) * (1. + ep2 * vaporav(i,j))
          tvskinbar   = tskinav(i,j) * (1. + ep2 * qskinav(i,j))
        end if

        !upcu    = 0.5 * (a_up(2,i,j) + a_up(2,i+1,j)) + umean
        !vpcv    = 0.5 * (a_vp(2,i,j) + a_vp(2,i,j+1)) + vmean
        !wspd2   = max(abs(ubmin), upcu**2. + vpcv**2.)
        Rib     = g / thetavbar * zt(2) * (thetavbar - tvskinbar) / wspd(i,j)**2.

        iter = 0
        if(obl(i,j) == 0.0) then
          L = 1e5
        else
          L = obl(i,j)
        end if

        if (Rib * L < 0. .or. abs(L) == 1e5) then
           if(Rib > 0) L = 0.01
           if(Rib < 0) L = -0.01
        end if

        do while (.true.)
           iter    = iter + 1
           Lold    = L
           fx      = Rib - zt(2) / L * (log(zt(2) / z0h(i,j)) - psih(zt(2)  &
                     / L) + psih(z0h(i,j) / L)) / (log(zt(2) / z0m(i,j)) -  &
                     psim(zt(2) / L) + psim(z0m(i,j) / L)) ** 2.

           Lstart  = L - 0.001*L
           Lend    = L + 0.001*L

           fxdif   = ( (- zt(2) / Lstart * (log(zt(2) / z0h(i,j)) - psih(zt(2) / Lstart) + psih(z0h(i,j) / Lstart)) / (log(zt(2) / z0m(i,j)) - psim(zt(2) / Lstart) + psim(z0m(i,j) / Lstart)) ** 2.) - (-zt(2) / Lend * (log(zt(2) / z0h(i,j)) - psih(zt(2) / Lend) + psih(z0h(i,j) / Lend)) / (log(zt(2) / z0m(i,j)) - psim(zt(2) / Lend) + psim(z0m(i,j) / Lend)) ** 2.) ) / (Lstart - Lend)

           L       = L - fx / fxdif

           if(Rib * L < 0. .or. abs(L) == 1e5) then
             if(Rib > 0) L = 0.01
             if(Rib < 0) L = -0.01
           end if

           if(abs(L) > 1e9) exit  ! BvS prevent runaway L in neutral limit
           if(iter > 990) print*,wspd(i,j),thetavbar,tvskinbar,zt(2),Rib,L  ! BvS debug
           if(abs(L - Lold) < 0.0001) exit
           if(iter > 1000) stop 'Obukhov length calculation does not converge!'
         end do

         obl(i,j) = L

      end do
    end do

    oblav = L

    return
  end subroutine getobl

  !
  ! ----------------------------------------------------------
  ! Malte: Integrated stability function for momentum
  !
  function psim(zeta)
    implicit none

    real             :: psim
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      x     = (1. - 16. * zeta) ** (0.25)
      !psim = 3.14159265/2. - atan(x) + 2.*log((1+x)/2.) + log((1+x*x)/2.)
      psim  = 3.14159265/2. - 2. *atan(x) + log((1.+x)** 2. * (1. + x**2.)/ 8.)
    else
      !psim = - 4.7 * zeta
      psim  = -2./3. * (zeta - 5./0.35) * exp(-0.35 * zeta) - zeta - (10./3.) / 0.35 
    end if

    return
  end function psim

  !
  ! ----------------------------------------------------------------------
  ! Malte: Integrated stability function for heat
  !
  function psih(zeta)

    implicit none

    real             :: psih
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      x     = (1. - 15. * zeta) ** (0.25)
      psih  = 2. * log( (1. + x ** 2.) / 2. )
    else
      !psih  = - 4.7 * zeta
      psih  = -2./3. * (zeta - 5./0.35) * exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    end if

    return
  end function psih

  !
  ! ----------------------------------------------------------------------
  ! Malte: Calculate the surface humidity assuming saturation. 
  ! DALES (vanHeerwarden)
  !
  subroutine qtsurf
    use defs, only : tmelt, R, Rm, rcp, p00
    use grid, only : nxp, nyp, vapor, psrf, a_tskin, a_qskin

    implicit none
    real       :: exner, tsurf, qsatsurf, surfwet, es, qtsl
    integer    :: i,j

    qtsl = 0.
    do j=3,nyp-2
      do i=3,nxp-2
        exner      = (psrf / p00)**rcp
        tsurf      = a_tskin(i,j) * exner
        es = 0.611e3 * exp(17.2694 * (tsurf - tmelt) / (tsurf - 35.86))

        qsatsurf   = R / Rm * es / psrf
        surfwet    = ra(i,j) / (ra(i,j) + rsurf(i,j))
        a_qskin(i,j) = surfwet * qsatsurf + (1. - surfwet) * vapor(2,i,j)
      end do
    end do

    return
  end subroutine qtsurf

  !
  ! ----------------------------------------------------------
  ! Malte: LSM to calculate temperature and moisture at the surface [DALES (van Heerwaarden)]
  ! PLEASE NOTE: LSM variables are defined and initialized in lsmdata.f90
  !
  subroutine lsm

    use defs, only: p00, stefan, rcp, cp, R, alvl, rowt, g
    use grid, only: nzp, nxp, nyp, a_up, a_vp, a_theta, vapor, liquid, zt, &
                    psrf, th00, umean, vmean, dn0, iradtyp, dt, &
                    a_lflxu, a_lflxd, a_sflxu, a_sflxd, nstep, press, &
                    a_lflxu_avn, a_lflxd_avn, a_sflxu_avn, a_sflxd_avn, &
                    a_tsoil, a_phiw, a_tskin, a_Wl, a_qskin, a_Qnet, a_G0
    use mpi_interface, only: myid

    integer  :: i, j, k
    real     :: f1, f2, f3, f4, fsoil !Correction functions for Jarvis-Stewart
    real     :: lflxu_av, lflxd_av, sflxu_av, sflxd_av
    real     :: tsurfm, Tatm, qskinn, exnerair, exner, pair, thetaavg, rhosurf
    real     :: e, esat, qsat, desatdT, dqsatdT, Acoef, Bcoef
    real     :: fH, fLE, fLEveg, fLEsoil, fLEliq, LEveg, LEsoil, LEliq
    real     :: Wlmx, rk3coef=0

    thetaavg  = sum(a_theta(2,3:(nxp-2),3:(nyp-2)))/(nxp-4)/(nyp-4)
    pair      = psrf*exp(-1.*g*(zt(2)-zt(1))/2/R/thetaavg)
    exner     = (psrf/p00)**rcp
    exnerair  = (pair/p00)**rcp   
    rhosurf   = (dn0(1)+dn0(2)) / 2.

    !"1.0 - Compute water content per layer
    do j = 3,nyp-2
      do i = 3,nxp-2

        phitot(i,j) = 0.0
        do k = 1, ksoilmax
          phitot(i,j) = phitot(i,j) + a_phiw(k,i,j) * dzsoil(k)
        end do
        phitot(i,j) = phitot(i,j) / zsoil(ksoilmax)

        do k = 1, ksoilmax
          phifrac(k,i,j) = a_phiw(k,i,j)*dzsoil(k) / zsoil(ksoilmax) / phitot(i,j)
        end do

      end do
    end do

    do j = 3, nyp-2
      do i = 3, nxp-2

        ! 1.1 - Calculate net radiation (average over nradtime for full radiation)
        if (iradtyp .eq. 4) then
          if (nstep == 1) then

            a_sflxd_avn(2:nradtime,i,j) = a_sflxd_avn(1:nradtime-1,i,j)
            a_sflxu_avn(2:nradtime,i,j) = a_sflxu_avn(1:nradtime-1,i,j)
            a_lflxd_avn(2:nradtime,i,j) = a_lflxd_avn(1:nradtime-1,i,j)
            a_lflxu_avn(2:nradtime,i,j) = a_lflxu_avn(1:nradtime-1,i,j)

            a_sflxd_avn(1,i,j) = a_sflxd(2,i,j)
            a_sflxu_avn(1,i,j) = a_sflxu(2,i,j)
            a_lflxd_avn(1,i,j) = a_lflxd(2,i,j)
            a_lflxu_avn(1,i,j) = a_lflxu(2,i,j)

          end if

          !Surface radiation averaged over nradtime but depends on location
          sflxd_av = sum(a_sflxd_avn(:,i,j))/nradtime
          sflxu_av = sum(a_sflxu_avn(:,i,j))/nradtime
          lflxd_av = sum(a_lflxd_avn(:,i,j))/nradtime
          lflxu_av = sum(a_lflxu_avn(:,i,j))/nradtime

          Qnetn(i,j) = (sflxd_av - sflxu_av + lflxd_av - lflxu_av)
          a_Qnet(i,j) = Qnetn(i,j)

        else if (iradtyp .eq. 5) then ! BvS: simple surface radiation doesn't require averaging

          sflxd_av = a_sflxd(2,i,j)
          sflxu_av = a_sflxu(2,i,j)
          lflxd_av = a_lflxd(2,i,j)
          lflxu_av = a_lflxu(2,i,j)

          Qnetn(i,j) = (sflxd_av - sflxu_av + lflxd_av - lflxu_av)
          a_Qnet(i,j) = Qnetn(i,j)

        else ! Not using full radiation: use average radiation from Namelist
          a_Qnet(i,j) = Qnetav
        end if

        !" 2.1 - Calculate the surface resistance with vegetation
        !" a) Stomatal opening as a function of incoming short wave radiation
        if ((iradtyp .eq. 4) .or. (iradtyp .eq. 5)) then
          f1 = 1. /min(1., (0.004 * max(0.,sflxd_av) + 0.05) &
                  / (0.81 * (0.004 * max(0.,sflxd_av) + 1.)) )
        else
          f1 = 1.
        end if

        !" b) Soil moisture availability
        f2   = (phifc - phiwp) / (phitot(i,j) - phiwp)
        !" Put upper bound f2 in case of very dry soils and prevent less than 1
        f2   = max(f2, 1.)
        f2   = min(1.e8, f2)

        !" c) Response of stomata to vapor deficit of atmosphere
        Tatm = a_theta(2,i,j)*exnerair
        esat = 0.611e3 * exp(17.2694 * (Tatm - 273.16) / (Tatm - 35.86))
        e    = vapor(2,i,j) * psrf / (0.622 + vapor(2,i,j))
        f3   = 1. / exp(-gD(i,j) * (esat - e) / 100.)

        !" d) Response to temperature
        f4   = 1./ (1. - 0.0016 * (298.0 - Tatm) ** 2.)
        rsveg(i,j) = rsmin(i,j) / LAI(i,j) * f1 * f2 * f3 * f4

        !" 2.2 - Calculate soil resistance based on ECMWF method
        fsoil = (phifc - phiwp) / (a_phiw(1,i,j) - phiwp)
        fsoil = max(fsoil, 1.)
        fsoil = min(1.e8, fsoil)
        rssoil(i,j) = rssoilmin(i,j) * fsoil

        !" 2.3 - Calculate the heat transport properties of the soil.
        !" Put in init function, as we don't have prognostic soil moisture atm.

        !"Save temperature and liquid water from previous timestep 
        if (nstep == 1) then
          tskinm(i,j) = a_tskin(i,j)
          Wlm(i,j)    = a_Wl(i,j)
        end if

        !"Solve the surface temperature implicitly including variations in LWout
        tsurfm  = tskinm(i,j) * exner
        esat    = 0.611e3 * exp(17.2694 * (tsurfm-273.16) / (tsurfm-35.86))
        qsat    = 0.622 * esat / (psrf-esat)
        desatdT = esat * (17.2694 / (tsurfm-35.86) - 17.2694 *  &
                  (tsurfm-273.16) / (tsurfm-35.86)**2.)
        dqsatdT = 0.622 * desatdT / (psrf-desatdT)

        !" Remove LWup from Qnet calculation (if running without radiation)
        a_Qnet(i,j) = a_Qnet(i,j) + stefan * (tskinm(i,j)*exner)**4.

        !" Allow for dew fall and calculate dew water on leaves
        if(qsat - vapor(2,i,j) < 0.) then
          rsveg(i,j)  = 0.
          rssoil(i,j) = 0.
        end if

        Wlmx      = LAI(i,j) * Wmax
        a_Wl(i,j) = min(a_Wl(i,j), Wlmx)
        cliq(i,j) = a_Wl(i,j) / Wlmx

        !" Calculate coefficients for surface fluxes
        fH      = rhosurf * cp / ra(i,j) 
        fLEveg  = (1. - cliq(i,j)) * cveg(i,j)  * rhosurf * alvl / (ra(i,j) + rsveg(i,j))
        fLEsoil = (1. - cveg(i,j))              * rhosurf * alvl / (ra(i,j) + rssoil(i,j))
        fLEliq  = (cliq(i,j)       * cveg(i,j)) * rhosurf * alvl / (ra(i,j))
        fLE     = fLEveg + fLEsoil + fLEliq

        !" Weighted timestep in runge-kutta scheme
        rk3coef = dt / (4 - float(nstep))

        !" Compute skin temperature from linarized surface energy balance
        Acoef   = a_Qnet(i,j) - stefan * (tskinm(i,j)*exner) **4. &
                  + 4. * stefan * (tskinm(i,j)*exner)**4. + fH*Tatm &
                  + fLE * (dqsatdT* (tskinm(i,j)*exner) - qsat + vapor(2,i,j)) &
                  + lambdaskin(i,j) * a_tsoil(1,i,j)
        Bcoef   = 4. * stefan * (tskinm(i,j)* exner) ** 3. + fH &
                  + fLE * dqsatdT + lambdaskin(i,j)

        if (Cskin(i,j) == 0.) then
          a_tskin(i,j) = Acoef * Bcoef ** (-1.) / exner
        else
          a_tskin(i,j) = (1. + rk3coef/Cskin(i,j) * Bcoef) ** (-1.) / exner * ((tskinm(i,j)*exner) + rk3coef/Cskin(i,j) * Acoef) 
        end if

        a_Qnet(i,j) = a_Qnet(i,j) - (stefan* (tskinm(i,j)*exner)**4.  &
                    + 4.*stefan * (tskinm(i,j)*exner)**3. *(a_tskin(i,j)*exner - &
                    tskinm(i,j)*exner))

        a_G0(i,j)   = lambdaskin(i,j) * ( a_tskin(i,j) * exner - a_tsoil(1,i,j) )

        qskinn    = (dqsatdT * (a_tskin(i,j)*exner - tskinm(i,j)*exner) + qsat)

        LE(i,j)   = - fLE     * ( vapor(2,i,j) - qskinn)
        LEveg     = - fLEveg  * ( vapor(2,i,j) - qskinn)
        LEsoil    = - fLEsoil * ( vapor(2,i,j) - qskinn)
        LEliq     = - fLEliq  * ( vapor(2,i,j) - qskinn)

        !if(i==5 .and. j==5) print*,LE(i,j),LEveg,LEsoil,LEliq

        if(LE(i,j) == 0.) then
          rsurf(i,j) = 1.e8           ! "WHY???
        else
          rsurf(i,j) = -rhosurf * alvl * (vapor(2,i,j) - qskinn) / LE(i,j) - ra(i,j) 
        end if

        H(i,j)    = - fH  * ( a_theta(2,i,j)*exner - a_tskin(i,j)*exner )
        tndskin(i,j) = Cskin(i,j)*(a_tskin(i,j) - tskinm(i,j)) * exner / rk3coef

        !" In case of dew formation, allow all water to enter skin reservoir a_Wl
        if(qsat - vapor(2,i,j) < 0.) then
          a_Wl(i,j) =  Wlm(i,j) - rk3coef * ((LEliq + LEveg) / (rowt * alvl))
        else
          a_Wl(i,j) =  Wlm(i,j) - rk3coef * ( LEliq          / (rowt * alvl))
        end if

        !" Save temperature and liquid water from previous timestep 
        if(nstep == 1) then
          tsoilm(:,i,j) = a_tsoil(:,i,j)
          phiwm(:,i,j)  = a_phiw(:,i,j)
        end if

        !" Calculate soil heat capacity and conductivity(based on water content)
        do k = 1, ksoilmax
          pCs(k,i,j)    = (1. - phi) * pCm + a_phiw(k,i,j) * pCw
          Ke            = log10(max(0.1,(a_phiw(k,i,j)/phi))) + 1.
          lambda(k,i,j) = Ke * (lambdasat - lambdadry) + lambdadry
        end do
     
        !" Calculate soil heat conductivity at half levels
        do k = 1, ksoilmax-1
          lambdah(k,i,j) = 0.5 * (lambda(k,i,j) * dzsoil(k) &
                           + lambda(k+1,i,j) * dzsoil(k+1)) / dzsoilh(k)
        end do
        lambdah(ksoilmax,i,j) = lambda(ksoilmax,i,j)

        !" Calculate soil moisture conductivity and difusivity 
        do k = 1, ksoilmax
          gammas(k,i,j)  = gammasat * (a_phiw(k,i,j) / phi) ** (2. * bc + 3.)
          lambdas(k,i,j) = bc * gammasat * (-1.) * psisat / phi &
                           * (a_phiw(k,i,j) / phi) ** (bc + 2.)
        end do







        !" Calculate soil moisture conductivity and difusivity at half levels
        do k = 1, ksoilmax-1
          lambdash(k,i,j) = 0.5 * (lambdas(k,i,j) * dzsoil(k) &
                            + lambdas(k+1,i,j) * dzsoil(k+1)) / dzsoilh(k)
          gammash(k,i,j)  = 0.5 * (gammas(k,i,j)  * dzsoil(k) &
                            + gammas(k+1,i,j)  * dzsoil(k+1)) / dzsoilh(k)
        end do
        lambdash(ksoilmax,i,j) = lambdas(ksoilmax,i,j)

        !" Solve the diffusion equation for the heat transport
        a_tsoil(1,i,j) = tsoilm(1,i,j) + rk3coef * ( lambdah(1,i,j) * (a_tsoil(2,i,j) - a_tsoil(1,i,j)) / dzsoilh(1) + a_G0(i,j) ) / dzsoil(1) / pCs(1,i,j)

        do k = 2, ksoilmax-1
          a_tsoil(k,i,j) = tsoilm(k,i,j) + rk3coef / pCs(k,i,j) * ( lambdah(k,i,j) * (a_tsoil(k+1,i,j) - a_tsoil(k,i,j)) / dzsoilh(k) - lambdah(k-1,i,j) * (a_tsoil(k,i,j) - a_tsoil(k-1,i,j)) / dzsoilh(k-1) ) / dzsoil(k)
        end do

        a_tsoil(ksoilmax,i,j) = tsoilm(ksoilmax,i,j) + rk3coef / pCs(ksoilmax,i,j) * ( lambda(ksoilmax,i,j) * (tsoildeep(i,j) - a_tsoil(ksoilmax,i,j)) / dzsoil(ksoilmax) - lambdah(ksoilmax-1,i,j) * (a_tsoil(ksoilmax,i,j) - a_tsoil(ksoilmax-1,i,j)) / dzsoil(ksoilmax-1) ) / dzsoil(ksoilmax)

        !" Solve the diffusion equation for the moisture transport (closed bottom for now)
        a_phiw(1,i,j) = phiwm(1,i,j) + rk3coef * ( lambdash(1,i,j) * (a_phiw(2,i,j) - a_phiw(1,i,j)) / dzsoilh(1) - gammash(1,i,j) - (phifrac(1,i,j) * LEveg + LEsoil) / (rowt*alvl)) / dzsoil(1)

        do k = 2, ksoilmax-1
          a_phiw(k,i,j) = phiwm(k,i,j) + rk3coef * ( lambdash(k,i,j) * (a_phiw(k+1,i,j) - a_phiw(k,i,j)) / dzsoilh(k) - gammash(k,i,j) - lambdash(k-1,i,j) * (a_phiw(k,i,j) - a_phiw(k-1,i,j)) / dzsoilh(k-1) + gammash(k-1,i,j) - (phifrac(k,i,j) * LEveg) / (rowt*alvl)) / dzsoil(k)
        end do

        a_phiw(ksoilmax,i,j) = phiwm(ksoilmax,i,j) + rk3coef * (- lambdash(ksoilmax-1,i,j) * (a_phiw(ksoilmax,i,j) - a_phiw(ksoilmax-1,i,j)) / dzsoil(ksoilmax-1) + gammash(ksoilmax-1,i,j) - (phifrac(ksoilmax,i,j) * LEveg) / (rowt*alvl) ) / dzsoil(ksoilmax)

      end do
    end do
     
    call qtsurf


  end subroutine lsm

end module srfc

