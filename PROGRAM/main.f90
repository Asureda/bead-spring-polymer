PROGRAM POLYMER
    !---------------------------------------------
    !      IMPORTING ALL MODULES
    !---------------------------------------------
    use READ_DATA
    use subroutines

    IMPLICIT NONE
    INTEGER :: i,j,seed,k
    REAL*8 :: norm
    CALL READ_FROM_FILES()
	  CALL DIMENSION_VARS()

    OPEN(12, FILE='positions.xyz')
	  OPEN(13,FILE='thermo.dat')
    OPEN(14,FILE='gyration.dat')
!The xyz is write in this way to better visualization in vmd
    DO k=1,Kmax
      print*,'k', k
    CALL SETSEED(seed)
    CALL SETR1279(seed)
    CALL init_positions(positions)

    WRITE(12,*)N
    WRITE(12,*)
    WRITE(12,*)'H', positions(1,:)
    DO i=2, N-1
        WRITE(12,*)'O', positions(i,:)
    END DO
    WRITE(12,*)'S', positions(N,:)

    DO i=1,steps_eq
        CALL integrator(positions,velos,forces)
        IF ((mod(i,10*steps_measure).eq.0).and.(traj_print.eqv..TRUE.))THEN
            WRITE(12,*)N
            WRITE(12,*)
            WRITE(12,*)'H', positions(1,:)
            DO j=2, N-1
                WRITE(12,*)'O', positions(j,:)
            END DO
            WRITE(12,*)'S', positions(N,:)

        END IF
    END DO

    WRITE(12,*)N
    WRITE(12,*)
    WRITE(12,*)'H', positions(1,:)
    DO i=2, N-1
        WRITE(12,*)'O', positions(i,:)
    END DO
    WRITE(12,*)'S', positions(N,:)

    DO i=1,steps_prod
        CALL integrator(positions,velos,forces)
        IF (mod(i,steps_measure).eq.0)THEN
            CALL SAMPLE(i,positions)
            WRITE(14,*) i, center_of_mass(positions), radius_gyration(positions),radius_gyration(positions)**0.5
        END IF

        IF ((mod(i,10*steps_measure).eq.0).and.(traj_print.eqv..TRUE.))THEN
            WRITE(12,*)N
            WRITE(12,*)
            WRITE(12,*)'H', positions(1,:)
            DO j=2, N-1
                WRITE(12,*)'O', positions(j,:)
            END DO
            WRITE(12,*)'S', positions(N,:)

        END IF



    END DO

    WRITE(12,*)N
    WRITE(12,*)
    WRITE(12,*)'H', positions(1,:)
    DO j=2, N-1
        WRITE(12,*)'O', positions(j,:)
    END DO
    WRITE(12,*)'S', positions(N,:)

  END DO
print*,'End program'
    END PROGRAM POLYMER
