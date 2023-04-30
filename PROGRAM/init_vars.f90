MODULE INIT_VARS
    use READ_DATA
    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: positions, velos, forces

    CONTAINS
    SUBROUTINE DIMENSION_VARS()
        IMPLICIT NONE
        ALLOCATE(positions(N,3),velos(N,3),forces(N,3))
        positions=0.
        velos=0.
        forces=0.

    END SUBROUTINE DIMENSION_VARS
END MODULE INIT_VARS
