      subroutine lworksizels( desc_a, desc_b, ia, ja, m, n, nrhs,lwork)
c 
c     =======================================
c			PURPOSE: compute LWORK >= LTAU + MAX( LWF, LWS ) which is the size of the work array required for the current least squares problem(from pdgels.f comments)
c			INPUTS: desc_a, desc_b, ia, ja, m, n, nrhs
c     OUTPUTS: lwork
c     =======================================
c 
      integer desc_a (9), desc_b (9), ia, ja, m, n, nrhs, lwork 
     $        m_a, n_a, mb_a, nb_a, rsrc_a, csrc_a, mb_b, nb_b, rsrc_b
     $        csrc_b, iroffa, icoffa, iarow, iacol, mpa0, nqa0, iroffb
     $        ibrow, ibcol, mpb0, npb0, nrhsqb0, lwf, lws, ltau
c     
c     external subroutines
      EXTERNAL BLACS_GRIDINFO
c
c     external functions
      integer NUMROC, INDXG2P
      EXTERNAL NUMROC, INDXG2P
c 
c     =======================================
c     extract needed information from global arrays desc_a and desc_b 
c     =======================================
c 
      ictxt = desc_a(2)
      m_a = desc_a(3)
      n_a = desc_a(4)
      mb_a = desc_a(5)
      nb_a = desc_a(6)
      rsrc_a = desc_a(7)
      csrc_a = desc_a(8)
      mb_b = desc_b(5)
      nb_b = desc_b(6)
      rsrc_b = desc_b(7)
      csrc_b = desc_b(8)
c 
c     input: ictxt, output: nprow, npcol, myrow, mycol
c 
      CALL BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)
c     
c     =======================================
c     compute the size of work (follows pdgels.f comments)
c     =======================================
c
      iroffa = mod( ia-1, mb_a )
      icoffa = mod( ja-1, nb_a )
      iarow = INDXG2P( ia, mb_a, myrow, rsrc_a, nprow )
      iacol = INDXG2P( ja, nb_a, mycol, csrc_a, npcol )
      mpa0 = NUMROC( m+iroffa, mb_a, myrow, iarow, nprow )
      nqa0 = NUMROC( n+icoffa, nb_a, mycol, iacol, npcol )
c
      iroffb = mod( ib-1, mb_b )
      icoffb = mod( jb-1, nb_b )
      ibrow = INDXG2P( ib, mb_b, myrow, rsrc_b, nprow )
      ibcol = INDXG2P( jb, nb_b, mycol, csrc_b, npcol )
      mpb0 = NUMROC( m+iroffb, mb_b, myrow, ibrow, nprow )
      npb0 = NUMROC( n+iroffb, mb_b, myrow, ibrow, nprow )
      nrhsqb0 = NUMROC( nrhs+icoffb, nb_b, mycol, ibcol, npcol )
c
      lwf  = nb_a * ( mpa0 + nqa0 + nb_a )
      lws  = max((nb_a*(nb_a-1))/2,(nrhsqb0 + mpb0)*nb_a )+nb_a * nb_a
      ltau = NUMROC( ja+min(m,n)-1, nb_a, mycol, csrc_a, npcol )
      lwork = ltau + max(lwf, lws) + 1
c
      return
      end
c
