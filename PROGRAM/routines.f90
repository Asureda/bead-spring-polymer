module subroutines
  use READ_DATA
  IMPLICIT NONE
  contains
!Routines to perform the simulation, 4 different potentials
  subroutine init_positions(r)
    IMPLICIT NONE
    integer :: i
    real(8) :: r(:,:)

    DO i = 1,N
      r(i,:) =i*(/1.d0,0.d0,0.d0/)
    END DO

  end subroutine init_positions

  subroutine integrator(r,vel,f)
    IMPLICIT NONE
    integer :: steps, i, j
    real(8) :: r(:,:),vel(:,:),f(:,:),v(N,3), new_r(3), factor,rnd_noise(3)
    factor = sqrt(2.d0*dt)
      F = 0.d0
      !call spring_interaction(r,F)
      !call wca_interaction(r,F)
      call fene_interaction(r,F)
      call wca_interaction2(r,F)
      DO i = 1,N
        CALL RANDOM_NOISE(rnd_noise)
        new_r = r(i,:) + dt*F(i,:)/T + factor*rnd_noise
        r(i,:) = new_r
        !r(i,:) = r(i,:)-int(2d0*r(i,:)/L)*L
      END DO
  end subroutine integrator

  SUBROUTINE RANDOM_NOISE(rnd)
      IMPLICIT NONE
      real(8) :: rnd(:),rn(4)
      REAL :: r1279,rnd1,rnd2,rnd3,rnd4,pi,rndmax
      pi = 4.d0*atan(1.d0)
      rndmax=200.d0

      50 rn=(/r1279(),r1279(),r1279(),r1279()/)
      rnd(1)=sqrt(-2.*log(rn(1)))*cos(2.*pi*rn(2))
      rnd(2)=sqrt(-2.*log(rn(1)))*sin(2.*pi*rn(2))
      rnd(3)=sqrt(-2.*log(rn(3)))*cos(2.*pi*rn(4))
      if (ANY(abs(rnd)>rndmax)) go to 50
      rnd(1)=sqrt(-2.*log(rn(1)))*cos(2.*pi*rn(2))
      rnd(2)=sqrt(-2.*log(rn(1)))*sin(2.*pi*rn(2))
      rnd(3)=sqrt(-2.*log(rn(3)))*cos(2.*pi*rn(4))

  END SUBROUTINE

  function PBC(x,L)
    IMPLICIT NONE
    real(8) :: x,L,PBC
    PBC = x-int(2d0*x/L)*L
    return
  end function PBC

  subroutine WCA(d,ff,cutoff)
    IMPLICIT NONE
    real(8) :: d,cutoff,ff,pot,r2i,r6i,virij,eps,sigma,cutoff2
    eps=1.d0
    sigma=1.d0
    cutoff =1.12246205 !2**(1.0d0/6.0d0)
    cutoff2= cutoff*cutoff
    if (d<cutoff2) then
        r2i = sigma**2/d
		    r6i = r2i**3
		    !pot = 4*eps*(r6i*r6i-r6i+1/4.d0)
        !ff=(48d0/d**13d0)-(24d0/d**7d0)
        virij = 48*(r6i*r6i-0.5D0*r6i)
        ff = virij*r2i
    end if
    return
  end subroutine WCA

    subroutine spring_interaction(r,F)
      IMPLICIT NONE
        integer :: i
        real(8) :: r(:,:),f(:,:),k,dx,dy,dz,dist,force,Dd
        k = 30.d0*T/(b*b)
        Dd = 1.d0
        do i=1,N-1
          dx=r(i,1)-r(i+1,1)
          dy=r(i,2)-r(i+1,2)
          dz=r(i,3)-r(i+1,3)
          dist=(dx**2d0+dy**2d0+dz**2d0)**0.5
          force=-K*(dist-Dd)*(1.d0/dist)
          F(i,1)=F(i,1)+force*dx
          F(i,2)=F(i,2)+force*dy
          F(i,3)=F(i,3)+force*dz
          F(i+1,1)=F(i+1,1)-force*dx
          F(i+1,2)=F(i+1,2)-force*dy
          F(i+1,3)=F(i+1,3)-force*dz
        end do
      end subroutine

    subroutine wca_interaction(r,F)
      IMPLICIT NONE
      integer :: i,j
      real(8) :: cutoff, pot
      real(8) :: dx,dy,dz,d,ff,d2
      real(8), dimension(:,:) :: r,F
      DO i = 1,N-2
        DO j = i+2, N
          dx=r(i,1)-r(j,1)
          dy=r(i,2)-r(j,2)
          dz=r(i,3)-r(j,3)
          d=(dx**2d0+dy**2d0+dz**2d0)**0.5
          d2 =d*d
          CALL WCA(d2,ff,cutoff)
          F(i,1)=F(i,1)+ff*dx
          F(i,2)=F(i,2)+ff*dy
          F(i,3)=F(i,3)+ff*dz
          F(j,1)=F(j,1)-ff*dx
          F(j,2)=F(j,2)-ff*dy
          F(j,3)=F(j,3)-ff*dz
        END DO
      END DO
    end subroutine
    subroutine wca_interaction2(r,F)
      IMPLICIT NONE
      integer :: i,j
      real(8) :: cutoff, pot,eps,sigma,cutoff2
      real(8) :: dx,dy,dz,d,ff,d1,d2,r2i,r6i,virij
      real(8), dimension(:,:) :: r,F
      eps=1.d0
      sigma=1.d0
      cutoff =1.12246205 !2**(1.0d0/6.0d0)
      cutoff2= cutoff*cutoff
      DO i = 1,N-1
        DO j = i+1, N
          dx=r(i,1)-r(j,1)
          if (abs(dx)<cutoff) then
            dy=r(i,2)-r(j,2)
            d1 = dx**2.d0+dy**2.d0
            if (d1<cutoff2) then
              dz=r(i,3)-r(j,3)
              d2 = d1 + dz**2.d0
              if (d2<cutoff2) then
                r2i = sigma**2/d2
                r6i = r2i**3
                virij = 48*(r6i*r6i-0.5D0*r6i)
                ff = virij*r2i
                F(i,1)=F(i,1)+ff*dx
                F(i,2)=F(i,2)+ff*dy
                F(i,3)=F(i,3)+ff*dz
                F(j,1)=F(j,1)-ff*dx
                F(j,2)=F(j,2)-ff*dy
                F(j,3)=F(j,3)-ff*dz
              end if
            end if
          end if
        END DO
      END DO
    end subroutine

    subroutine fene_interaction(r,F)
      IMPLICIT NONE
        integer :: i
        real(8) :: r(:,:),f(:,:),k,dx,dy,dz,dist,force,Dd,Ro,eps,sigma
        eps=1.d0
        sigma=1.d0
        k = 3.d0*T/(b*b)
        k = 3.d0*eps/sigma**2
        Dd = 1.d0
        Ro = 2.d0
        do i=1,N-1
          dx=r(i,1)-r(i+1,1)
          dy=r(i,2)-r(i+1,2)
          dz=r(i,3)-r(i+1,3)
          dist=(dx**2d0+dy**2d0+dz**2d0)**0.5
          if (abs(dist)<Ro) then
          force=-K*(1/dist)*(dist-Dd)/(1.d0-(dist/Ro)**2.d0)
          end if
          F(i,1)=F(i,1)+force*dx
          F(i,2)=F(i,2)+force*dy
          F(i,3)=F(i,3)+force*dz
          F(i+1,1)=F(i+1,1)-force*dx
          F(i+1,2)=F(i+1,2)-force*dy
          F(i+1,3)=F(i+1,3)-force*dz
        end do
      end subroutine

    function center_of_mass(r)
      IMPLICIT NONE
      integer :: i
      real(8) :: center_of_mass, r(:,:),d,center_of_massx,center_of_massy,center_of_massz
      center_of_massx = sum(r(:,1))/dble(N)
      center_of_massy = sum(r(:,2))/dble(N)
      center_of_massz = sum(r(:,3))/dble(N)
      center_of_mass = sqrt(center_of_massx**2d0+center_of_massy**2d0 + center_of_massz**2d0)
      !center_of_mass = (/center_of_massx,center_of_massy,center_of_massz/)
    end function center_of_mass

     function radius_gyration(r)
       IMPLICIT NONE
       integer :: i
       real(8) ::  r(:,:),d,radius,radius_gyration,cmx,cmy,cmz
       real(8) :: radius_gyrationx,radius_gyrationy,radius_gyrationz

        cmx = sum(r(:,1))/dble(N)
        cmy = sum(r(:,2))/dble(N)
        cmz = sum(r(:,3))/dble(N)
       do i =1,N
         radius_gyrationx =radius_gyrationx + (r(i,1)-cmx)**2
         radius_gyrationy =radius_gyrationy + (r(i,2)-cmy)**2
         radius_gyrationz =radius_gyrationz + (r(i,3)-cmz)**2
       end do
       radius_gyrationx =radius_gyrationx/dble(N)
       radius_gyrationy =radius_gyrationy/dble(N)
       radius_gyrationz =radius_gyrationz/dble(N)
       radius_gyration = sqrt(radius_gyrationx**2d0+radius_gyrationy**2d0 + radius_gyrationz**2d0)

     end function radius_gyration

    subroutine sample(iter,rp)
      IMPLICIT NONE
      integer :: i,iter
      real(8) :: rp(:,:),tmp,sum2,end2end,dist,pi,end1
      pi = 4.d0*atan(1.d0)
      DO i = 1,N-1
        dist = (rp(i,1)-rp(i+1,1))**2.+(rp(i,2)-rp(i+1,2))**2.+(rp(i,3)-rp(i+1,3))**2.
      END DO

      end1 = ((rp(1,1)-rp(N,1))**2.+(rp(1,2)-rp(N,2))**2.+(rp(1,3)-rp(N,3))**2.)**0.5
      end2end=end1**2.d0
      WRITE(13,*)iter,T,end2end, end1
    end subroutine sample
end module subroutines
