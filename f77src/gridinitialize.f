      subroutine leastsquares(ictxt, nprow, npcol, myrow, mycol, m, n, 
     $           rowblock, colblock, submat, maxLLD, locLLD, LOCc, 
     $           locLLD_V, icpu, IA, JA, U, S, V, lwork, work, info)
      integer ictxt, nprow, npcol, myrow, mycol, m, n, rowblock
      integer colblock, info, locLLD, icpu, rowindex, lwork
      integer maxLLD, k, j, locLLD_V
c
      double precision submat(maxLLD, LOCc), V(LOCc, locLLD_V)
c     double precision S(n), U(maxLLD, LOCc), work(lwork+2)
      double precision S(n), U(locLLD, LOCc), work(lwork+2)
      integer desca(9), descu(9), descvt(9)
c     
      EXTERNAL BLACS_EXIT, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $         DESCINIT, PDGESVD, SL_INIT
c			initialize an NPROW x NPCOL process grid, retrieve a default system context
      CALL SL_INIT(ictxt, nprow, npcol)
c			obtain process row and column index
      CALL BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)
c      WRITE( 6, FMT = 10)icpu, myrow, mycol 
c
c10    FORMAT( 'cpu ', I3, ' running on cpu grid ', I3, ' by ', I3 )
c     
c			initialize descriptor vector for A, U, and V^T
      CALL DESCINIT(desca, m, n, rowblock, colblock, 0, 0, 
     $              ictxt, locLLD, info)
c
c      WRITE( 6, FMT = 11)icpu, info
c
c11    FORMAT( 'cpu ', I3, ' info for desca: ', I3  )
c
      CALL DESCINIT(descu, m, n, rowblock, colblock, 0, 0, 
     $              ictxt, locLLD, info)
c 
c      WRITE( 6, FMT = 12)icpu, info
c
c12    FORMAT( 'cpu ', I3, ' info for descu: ', I3  )

      CALL DESCINIT(descvt, n, n, colblock, colblock, 0, 0, ictxt, 
     $              locLLD_V, info)
c    
c      WRITE( 6, FMT = 20)icpu, lwork 
c
c20    FORMAT( 'cpu ', I3, ' info for lwork: ', I10  )
c
      CALL PDGESVD('V', 'V', m, n, submat, 1, 1, desca, S, U, 
     $              1, 1, descu, V, 1, 1, descvt, work, lwork, 
     $              info)
c
      IF (info.ne.0) THEN
      WRITE( 6, FMT = 30)icpu, info 
c
30    FORMAT( 'cpu ', I3, ' failed in svd with info: ', I3  )
c
      end if
      return
      end
c
