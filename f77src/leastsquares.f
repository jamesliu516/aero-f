      subroutine leastsquares(m, n, nrhs, subMatA, ia, ja, desc_a, 
     $    subMatB, ib, jb, desc_b, work, lwork, subMatLLD, icpu,info)
      integer m, n, nrhs, ia, ja, ib, jb, subMatLLD, info, icpu,lwork
      integer desc_a(9), desc_b(9)
      double precision subMatA(subMatLLD, n), subMatB(subMatLLD, nrhs)
      double precision work(lwork)
c     
      EXTERNAL PDGELS
c
      WRITE( 6, FMT = 20)icpu, lwork
c
20    FORMAT( 'cpu ', I3, ' lwork: ', I8)
c      
      CALL PDGELS('N', m, n, nrhs, subMatA, ia, ja, desc_a,
     $            subMatB, ib, jb, desc_b, work, lwork, info )
c
c     ================================
c     indicate if error
c     ================================
c
      IF (info.ne.0) THEN
      WRITE( 6, FMT = 30)icpu, info 
c
30    FORMAT( 'cpu ', I3, ' failed in least squares with info: ', I3  )
c
      end if
      return
      end
c
