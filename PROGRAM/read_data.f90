MODULE READ_DATA
    IMPLICIT NONE
    include 'r1279block.h'
    INTEGER :: N
    REAL :: b , T,L,dt,t_eq,t_prod,t_meas
    !
    INTEGER :: steps_eq, steps_prod, steps_measure
    !
    INTEGER*16, DIMENSION(:), ALLOCATABLE :: histo_y
    INTEGER :: Kmax
    LOGICAL :: traj_print

    CONTAINS
    SUBROUTINE READ_FROM_FILES()
        IMPLICIT NONE
        INTEGER i
        OPEN (11,FILE='parameters.dat',STATUS='old')
        READ(11,*) N
        READ(11,*) b
        READ(11,*) Kmax
        READ(11,*) L
        READ(11,*) T
        READ(11,*) dt
        READ(11,*) t_eq
        READ(11,*) t_prod
        READ(11,*) t_meas
        READ(11,*) traj_print
        close(11)
        steps_eq=int(t_eq/dt)
        steps_prod=int(t_prod/dt)
        steps_measure=int(t_meas/dt)
        print*,steps_measure,steps_prod,steps_eq

  END SUBROUTINE READ_FROM_FILES
  END MODULE READ_DATA
