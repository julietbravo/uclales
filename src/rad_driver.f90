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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module radiation

  use defs, only          : cp, rcp, cpr, rowt, roice, p00, pi, nv1, nv, SolarConstant
  use fuliou, only        : rad, minSolarZenithCosForVis
  use mpi_interface, only : myid, pecount, double_array_par_sum_root, broadcast_array, nxpg, nypg

  implicit none
  character (len=10), parameter :: background = 'backrad_in'

  logical, save :: first_time = .True.
  real, allocatable, save :: &
                   pp(:),   pt(:),   ph(:),   po(:),   pre(:), pde(:), &
                   plwc(:), piwc(:), prwc(:), pgwc(:), fds(:), fus(:), fdir(:), fuir(:)

  integer :: k, i, j, npts
  real    :: ee, u0, day, time, alat, zz
  logical :: fixed_sun = .false.
  logical :: radMcICA = .true.

  ! BvS
  logical :: heating_atmosphere = .true. ! flag to disable radiative heating in atmosphere
  logical :: single_column_rad = .false. ! use domain mean profiles for (single column) radiation, apply to 3d field
  real    :: rad_eff_radius = 1.

  ! BvS: for simple surface radiation
  logical :: fixed_lwin = .false.
  real    :: flwin      = 300.
  ! end BvS

  contains
    subroutine d4stream(n1,n2,n3,alat,time,sknt,sfc_albedo,CCN,dn0, &
                        pi0,pi1,dzi_m,pip,th,rv,rc,tt,rflx,sflx,lflxu,&
                        lflxd,sflxu,sflxd,albedo,lflxu_toa,lflxd_toa,&
                        sflxu_toa,sflxd_toa,rr,ice,nice,grp)

      integer, intent (in) :: n1, n2, n3
      real, intent (in)    :: alat, time, sknt, sfc_albedo, CCN
      real, dimension (n1), intent (in)                 :: dn0, pi0, pi1, dzi_m
      real, dimension (n1,n2,n3), intent (in)           :: pip, th, rv, rc
      real, optional, dimension (n1,n2,n3), intent (in) :: rr,ice,nice,grp
      real, dimension (n1,n2,n3), intent (inout)        :: tt, rflx, sflx, lflxu, lflxd, sflxu, sflxd
      real, dimension (n2,n3), intent (out),optional    :: albedo, lflxu_toa, lflxd_toa, sflxu_toa, sflxd_toa

      integer :: kk
      real    :: ijtoti ! (1/((nxp-4)*(nyp-4))
      real    :: xfact, prw, pri, p0(n1), exner(n1), pres(n1)
      real    :: th_suml(n1), th_p(n1), rv_suml(n1), rv_p(n1), ttend(n1)

      if (first_time) then
        p0(n1) = (p00*(pi0(n1)/cp)**cpr) / 100.
        p0(n1-1) = (p00*(pi0(n1-1)/cp)**cpr) / 100.
        call setup(background,n1,npts,nv1,nv,p0,pi0)
        first_time = .False.
        if (allocated(pre))   pre(:) = 0.
        if (allocated(pde))   pde(:) = 0.
        if (allocated(piwc)) piwc(:) = 0.
        if (allocated(prwc)) prwc(:) = 0.
        if (allocated(plwc)) plwc(:) = 0.
        if (allocated(pgwc)) pgwc(:) = 0.
      end if
      !
      ! initialize surface albedo, emissivity and skin temperature.
      !
      ee = 1.0
      !
      ! determine the solar geometery, as measured by u0, the cosine of the
      ! solar zenith angle
      !
      if (.not. fixed_sun) then
        u0 = zenith(alat,time)
      end if

      !
      ! call the radiation
      !
      prw = (4./3.)*pi*rowt
      pri = (3.*sqrt(3.)/8.)*roice
 
      ! 
      ! Calculate domain mean profiles and use them for radiation 
      !
      if(single_column_rad) then
    
        ! Calculate domain mean temperature and moisture profiles
        do k=1,n1
          ! Sum per MPI process
          th_suml(k) = sum(th(k,3:n2-2,3:n3-2))
          rv_suml(k) = sum(rv(k,3:n2-2,3:n3-2))
          ! Calculate pressure 
          exner(k) = (pi0(k)+pi1(k))/cp ! neglect pip(k,i,j) for now
          pres(k)  = p00 * (exner(k))**cpr   
        end do

        call double_array_par_sum_root(th_suml, th_p, n1)
        call double_array_par_sum_root(rv_suml, rv_p, n1)

        ! Only myid==0 does the radiation
        if(myid == 0) then
          ijtoti = 1./((nxpg-4)*(nypg-4)) 
          do k=1,n1
            ! Domain mean profiles:
            th_p(k) = th_p(k) * ijtoti 
            rv_p(k) = rv_p(k) * ijtoti 
          end do

          ! Put input radiation reversed into arrays:
          pp(nv1) = 0.5*(pres(1)+pres(2)) / 100.
          do k=2,n1
            kk = nv-(k-2) ! reversed index

            pt(kk)   = th_p(k) * exner(k)
            ph(kk)   = max(0., rv_p(k))

            plwc(kk) = 0. ! for now no liquid etc.
            pre(kk)  = 0.
            prwc(kk) = 0.
            piwc(kk) = 0.
            pde(kk)  = 0.
            pgwc(kk) = 0.

            if (k < n1) pp(kk) = 0.5*(pres(k)+pres(k+1)) / 100.

          end do 
 
          pp(nv-n1+2) = pres(n1)/100. - 0.5*(pres(n1-1)-pres(n1)) / 100.

          ! Radiation solver
          call rad(sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                   fds, fus, fdir, fuir, plwc=plwc, pre=pre, useMcICA=radMcICA)
        end if
        
        ! Broadcast radiation profiles to all MPI processes 
        call broadcast_array(fds,  0, nv1)
        call broadcast_array(fus,  0, nv1)
        call broadcast_array(fdir, 0, nv1)
        call broadcast_array(fuir, 0, nv1)

        do j=3,n3-2
          do i=3,n2-2
            do k=1,n1
              kk = nv1 - (k-1)
              sflx(k,i,j)  = fus(kk)  - fds(kk)
              sflxu(k,i,j) = fus(kk)
              sflxd(k,i,j) = fds(kk)
              lflxu(k,i,j) = fuir(kk)
              lflxd(k,i,j) = fdir(kk)
              rflx(k,i,j)  = sflx(k,i,j) + fuir(kk) - fdir(kk)
            end do

            if (present(albedo)) then
              if (u0 > minSolarZenithCosForVis) then
                albedo(i,j) = fus(1)/fds(1)
              else
                albedo(i,j) = -999.
              end if
            end if

            if (present(sflxu_toa)) then
              if (u0 > minSolarZenithCosForVis) then
                sflxu_toa(i,j) = fus(1)
              else
                sflxu_toa(i,j) = -999.
              end if
            end if
            if (present(sflxd_toa)) then
              if (u0 > minSolarZenithCosForVis) then
                sflxd_toa(i,j) = fds(1)
              else
                sflxd_toa(i,j) = -999.
              end if
            end if
            if (present(lflxu_toa)) then
              lflxu_toa(i,j) = fuir(1)
            end if
            if (present(lflxd_toa)) then
              lflxd_toa(i,j) = fdir(1)
            end if

            if(heating_atmosphere) then
              do k=2,n1-3
                xfact  = dzi_m(k)/(cp*dn0(k)*exner(k))
                tt(k,i,j) = tt(k,i,j) - (rflx(k,i,j) - rflx(k-1,i,j))*xfact
              end do
            end if
          end do
        end do
      ! 
      ! Radiation for each column 
      !
      else 
        do j=3,n3-2
          do i=3,n2-2
            do k=1,n1
              exner(k) = (pi0(k)+pi1(k)+pip(k,i,j))/cp
              pres(k) = p00 * (exner(k))**cpr
            end do
            pp(nv1) = 0.5*(pres(1)+pres(2)) / 100.
            do k=2,n1
              kk = nv-(k-2)
              pt(kk) = th(k,i,j)*exner(k)
              ph(kk) = max(0.,rv(k,i,j))
              plwc(kk) = 1000.*dn0(k)*max(0.,rc(k,i,j))
              pre(kk) = rad_eff_radius*1.e6*(plwc(kk)/(1000.*prw*CCN*dn0(k)))**(1./3.)
              pre(kk) = min(max(pre(kk),4.18),31.23)
              if (plwc(kk).le.0.) pre(kk) = 0.
              if (present(rr)) then
                prwc(kk) = 1000.*dn0(k)*rr(k,i,j)
              else
                prwc(kk) = 0.
              end if
              if (present(ice)) then
                piwc(kk) = 1000.*dn0(k)*ice(k,i,j)
                if (nice(k,i,j).gt.0.0) then
                  pde(kk) = 1.e6*(piwc(kk)/(1000.*pri*nice(k,i,j)*dn0(k)))**(1./3.)
                  pde(kk) = min(max(pde(kk),20.),180.)
                else
                  pde(kk) = 0.0
                endif
              else
                piwc(kk) = 0.
                pde(kk) = 0.0
              end if
              if (present(grp)) then
                pgwc(kk) = 1000.*dn0(k)*grp(k,i,j)
              else
                pgwc(kk) = 0.
              end if
              if (k < n1) pp(kk) = 0.5*(pres(k)+pres(k+1)) / 100.
            end do
            pp(nv-n1+2) = pres(n1)/100. - 0.5*(pres(n1-1)-pres(n1)) / 100.

            if (present(ice).and.present(grp)) then
              call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                   fds, fus, fdir, fuir, plwc=plwc, pre=pre, piwc=piwc, pde=pde, pgwc=pgwc, useMcICA=radMcICA)
            else
              call rad( sfc_albedo, u0, SolarConstant, sknt, ee, pp, pt, ph, po,&
                   fds, fus, fdir, fuir, plwc=plwc, pre=pre, useMcICA=radMcICA)
            end if

            do k=1,n1
              kk = nv1 - (k-1)
              sflx(k,i,j)  = fus(kk)  - fds(kk)
              sflxu(k,i,j) = fus(kk)
              sflxd(k,i,j) = fds(kk)
              lflxu(k,i,j) = fuir(kk)
              lflxd(k,i,j) = fdir(kk)
              rflx(k,i,j)  = sflx(k,i,j) + fuir(kk) - fdir(kk)
            end do

            if (present(albedo)) then
              if (u0 > minSolarZenithCosForVis) then
                albedo(i,j) = fus(1)/fds(1)
              else
                albedo(i,j) = -999.
              end if
            end if

            if (present(sflxu_toa)) then
              if (u0 > minSolarZenithCosForVis) then
                sflxu_toa(i,j) = fus(1)
              else
                sflxu_toa(i,j) = -999.
              end if
            end if
            if (present(sflxd_toa)) then
              if (u0 > minSolarZenithCosForVis) then
                sflxd_toa(i,j) = fds(1)
              else
                sflxd_toa(i,j) = -999.
              end if
            end if
            if (present(lflxu_toa)) then
              lflxu_toa(i,j) = fuir(1)
            end if
            if (present(lflxd_toa)) then
              lflxd_toa(i,j) = fdir(1)
            end if

            if(heating_atmosphere) then
              do k=2,n1-3
                xfact  = dzi_m(k)/(cp*dn0(k)*exner(k))
                tt(k,i,j) = tt(k,i,j) - (rflx(k,i,j) - rflx(k-1,i,j))*xfact
              end do
            end if
          end do
        end do
      end if

    end subroutine d4stream

  !
  ! ---------------------------------------------------------------------------
  ! BvS: Simple parameterized surface radiation for LSM
  !
  subroutine surfacerad(alat,time)
    use grid, only   : sfc_albedo,a_theta,a_lflxu,a_lflxd,a_sflxu,a_sflxd,a_tskin,nxp,nyp,a_pexnr,pi0,pi1
    use defs, only   : stefan, cp
    real, intent(in) :: time, alat
    integer          :: i,j
    real             :: tr, exner
    real, parameter  :: epss = 1.  
    real, parameter  :: epsa = 0.8   

    ! Assumes longitude = 0.
    u0 = max(0.,zenith(alat,time))
    tr = (0.6 + 0.2 * u0)

    do j=3,nyp-2
      do i=3,nxp-2
        exner          = (pi0(2)+pi1(2)+a_pexnr(2,i,j))/cp
        a_sflxd(2,i,j) = SolarConstant * tr * u0
        a_sflxu(2,i,j) = sfc_albedo * a_sflxd(2,i,j)
        if(fixed_lwin) then
          a_lflxd(2,i,j) = flwin
        else
          a_lflxd(2,i,j) = epsa * stefan * (a_theta(2,i,j)*exner)**4.
        end if
        a_lflxu(2,i,j) = epss * stefan * a_tskin(i,j)**4.
      end do
    end do

    !print*,'SEB:',a_sflxd(2,10,10),a_sflxu(2,10,10),a_lflxd(2,10,10),a_lflxu(2,10,10)
    !print*,'NET:',a_sflxd(2,10,10)-a_sflxu(2,10,10)+a_lflxd(2,10,10)-a_lflxu(2,10,10)

  end subroutine surfacerad

  ! ---------------------------------------------------------------------------
  ! sets up the input data to extend through an atmopshere of appreiciable
  ! depth using a background sounding specified as a parameter, match this to
  ! the original sounding using p0 as this does not depend on time and thus
  ! allows us to recompute the same background matching after a history start
  !
  subroutine setup(background,n1,npts,nv1,nv,zp,pi0)

    character (len=10), intent (in) :: background
    integer, intent (in) :: n1
    integer, intent (out):: npts,nv1,nv
    real, intent (in)    :: zp(n1), pi0(n1)

    real, allocatable  :: sp(:), st(:), sh(:), so(:), sl(:)

    integer :: k, ns, norig, index
    logical :: blend
    real    :: pa, pb, ptop, ptest, test, dp1, dp2, dp3, Tsurf, pp2

    open ( unit = 08, file = background, status = 'old' )
    print *, 'Reading Background Sounding: ',background
    read (08,*) Tsurf, ns
    allocate ( sp(ns), st(ns), sh(ns), so(ns), sl(ns))
    do k=1,ns
       read ( 08, *) sp(k), st(k), sh(k), so(k), sl(k)
    enddo
    close (08)

    !
    ! identify what part, if any, of background sounding to use
    !
    ptop = zp(n1)
    if (sp(2) < ptop) then
       pa = sp(1)
       pb = sp(2)
       k = 3
       do while (sp(k) < ptop)
          pa = pb
          pb = sp(k)
          k  = k+1
       end do
       k=k-1           ! identify first level above top of input
       blend = .True.
    else
       blend = .False.
    end if
    !
    ! if blend is true then the free atmosphere above the sounding will be
    ! specified based on the specified background climatology, here the
    ! pressure levels for this part of the sounding are determined
    !
    if (blend) then
       dp1 = pb-pa
       dp2 = ptop - pb
       dp3 = zp(n1-1) - zp(n1)
       if (dp1 > 2.*dp2) k = k-1 ! first level is too close, blend from prev
       npts  = k
       norig = k
       ptest = sp(k)
       test = ptop-ptest
       do while (test > 2*dp3)
          ptest = (ptest+ptop)*0.5
          test  = ptop-ptest
          npts  = npts + 1
       end do
       nv1 = npts + n1
    else
       nv1 = n1
    end if
    nv = nv1-1
    !
    ! allocate the arrays for the sounding data to be used in the radiation
    ! profile and then fill them first with the sounding data, by afill, then
    ! by interpolating the background profile at pressures less than the
    ! pressure at the top fo the sounding
    !
    allocate (pp(nv1),fds(nv1),fus(nv1),fdir(nv1),fuir(nv1))
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),prwc(nv),piwc(nv),pgwc(nv))

    if (blend) then
       pp(1:norig) = sp(1:norig)
       pt(1:norig) = st(1:norig)
       ph(1:norig) = sh(1:norig)
       po(1:norig) = so(1:norig)

       do k=norig+1,npts
          pp(k) = (ptop + pp(k-1))*0.5
          index = getindex(sp,ns,pp(k))
          pt(k) =  intrpl(sp(index),st(index),sp(index+1),st(index+1),pp(k))
          ph(k) =  intrpl(sp(index),sh(index),sp(index+1),sh(index+1),pp(k))
          po(k) =  intrpl(sp(index),so(index),sp(index+1),so(index+1),pp(k))
       end do
      ! !
      ! ! set the ozone constant below the reference profile
      ! !
      ! do k=npts+1,nv
      !    po(k) =  po(npts)
      ! end do

      ! interpolate ozone profile
       do k=npts+1,nv
            pp2 = (p00*(pi0(nv-k+2)/cp)**cpr) / 100.
            index  = getindex(sp,ns,pp2)
            po(k) =  intrpl(sp(index),so(index),sp(index+1),so(index+1),pp2)
         end do

    end if
  
  end subroutine setup
  ! ---------------------------------------------------------------------------
  !  locate the index closest to a value
  !
  integer function getindex(x,n,xval)

    integer, intent (in)  :: n
    real,    intent (in)  :: x(n),xval

    integer :: ia, ib

    ia=1
    ib=n
    if (xval < x(1)) then
       getindex = 1
    elseif (xval > x(n)) then
       getindex = n-1
    else
       getindex = (ia+ib)/2
       do while (getindex /= ia .or. ib /= getindex+1)
          getindex = (ia+ib)/2
          if ((xval-x(getindex)) >= 0.0) then
             ia = getindex
          else
             ib = getindex
          end if
       end do
    endif

  end function getindex

  ! ---------------------------------------------------------------------------
  ! linear interpolation between two points,
  !
  real function intrpl(x1,y1,x2,y2,x)

    real, intent (in)  :: x1,y1,x2,y2,x

    real :: slope

    slope  = (y2-y1)/(x2 - x1 + epsilon(1.))
    intrpl = y1+slope*(x-x1)

  end function intrpl

  ! ---------------------------------------------------------------------------
  ! Return the cosine of the solar zenith angle give the decimal day and
  ! the latitude
  !
  real function zenith(alat,time)

    real, intent (in)  :: alat, time

    real :: lamda, d, sig, del, h, day

    day    = floor(time)
    lamda  = alat*pi/180.
    d      = 2.*pi*int(time)/365.
    sig    = d + pi/180.*(279.9340 + 1.914827*sin(d) - 0.7952*cos(d) &
         &                      + 0.019938*sin(2.*d) - 0.00162*cos(2.*d))
    del    = asin(sin(23.4439*pi/180.)*sin(sig))
    h      = 2.*pi*((time-day)-0.5)
    zenith = sin(lamda)*sin(del) + cos(lamda)*cos(del)*cos(h)

  end function zenith

end module radiation
