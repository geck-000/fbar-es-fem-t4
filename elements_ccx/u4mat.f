!
!     SPAX: mixed displacement/pressure tetrahedron for near-incompressible
!     phases -- the element CalculiX does not have and Abaqus calls C3D4H.
!
!     FORMULATION: MINI, P1(+bubble)/P1.
!
!     Why not the literal P1/P0 that Abaqus documents for C3D4H: on a linear
!     tet the strain is constant over the element, so an element-constant
!     pressure represents the divergence exactly, static condensation of that
!     pressure is algebraically the mean-dilatation B-bar operator, and B-bar
!     on a one-point tet is the identity. A faithful P1/P0 element in ccx is
!     therefore bit-for-bit plain C3D4 -- measured, in
!     patches/0002-bbar-mean-dilatation.patch. It cannot close any gap.
!
!     P1/P1 without enrichment violates the inf-sup condition (checkerboard
!     pressure). MINI restores it parameter-free by enriching the displacement
!     with a cubic bubble, which is element-internal and so condenses away
!     locally, leaving 4 nodes x 4 dof.
!
!     DOF LAYOUT: node-major, 4 dof per node -- 1,2,3 displacement, 4 pressure.
!     Declare in the deck with
!
!         *USER ELEMENT,TYPE=U4,NODES=4,INTEGRATIONPOINTS=1,MAXDOF=4
!
!     MAXDOF=4 raises mi(2) in allocation.f, which is what gives mastruct.c
!     room for the fourth nodal dof.
!
!     THE SYSTEM IS SYMMETRIC INDEFINITE (saddle point). Incomplete-Cholesky
!     PCG cannot be used on it; solve with SOLVER=SPOOLES.
!
!
!     u4mat: the whole element operator in one place. e_c3d_u4 turns it into a
!     stiffness matrix and resultsmech_u4 turns it into internal forces, so the
!     two cannot drift apart -- the failure mode that leaves
!     0002-bbar-mean-dilatation.patch with an equilibrium gap of 1.0.
!
!     Returns the bubble-condensed blocks aull (12x12), bl (4x12), cc (4x4),
!     and the pieces abi/alb/bb needed to reconstruct the bubble amplitude.
!
      subroutine u4mat(xl,um,xk,aull,bl,cc,abi,alb,bb,nelem)
!
!     Closed-form element operator. No quadrature loop.
!
!     Every integral needed is a monomial in the barycentric coordinates over a
!     straight tetrahedron, for which
!
!         int L1^a L2^b L3^c L4^d dV = 6V a! b! c! d! / (a+b+c+d+3)!
!
!     With g_i = grad(L_i) constant and sum_i g_i = 0, and the cubic bubble
!     b = 256 L1 L2 L3 L4:
!
!         int L_i dV        = V/4
!         int L_i L_j dV    = V/10 (i=j), V/20 (i/=j)
!         int L_i db/dx_d   = -(256 V/840) g_i[d]
!         int db/dx_c db/dx_d = (256^2 V/15120) sum_i g_i[c] g_i[d]
!
!     A_lb IS EXACTLY ZERO. The bubble vanishes on every face, so
!     int grad(b) dV = closed surface integral of b n dS = 0, and every
!     linear-bubble coupling term carries a factor of it. Verified numerically
!     against the previous 15-point implementation: max|A_lb|/max|A_ll| = 7e-16.
!
!     That collapses the condensation to one line -- the bubble's ONLY effect is
!     a stabilisation added to the pressure block:
!
!         A_ll' = A_ll ,  B_l' = B_l ,  C' = C + B_b A_bb^-1 B_b^T
!
!     so U4 is P1/P1 with a parameter-free stabilisation, which is the standard
!     characterisation of MINI. Two consequences worth stating. The element can
!     declare INTEGRATIONPOINTS=1, which matters because ccx sizes sti, eme,
!     xstiff and stx at mi(1) for EVERY element in the model -- 1.74 GB against
!     0.12 GB on a 345k-element cell. And because C' is invertible, the pressure
!     admits a global Schur elimination A + B^T C'^-1 B, which is symmetric
!     POSITIVE DEFINITE; see ../elements_ccx/README.md on solvers.
!
      implicit none
!
      integer i,j,c,d,ii,jj,nelem
      real*8 xl(3,4),um,xk,aull(12,12),bl(4,12),cc(4,4),abi(3,3),
     &     alb(12,3),bb(4,3),abb(3,3),shp(4,4),xsj,vol,g(3,4),
     &     gm(3,3),gbb,det,fac,cbb,cbl,tmp2(3,4),trc,trs,theta,
     &     stab(4,3),strow(4),xlam
      character*132 spaxenv
      integer istab,ivol
      save istab,ivol
      data istab /-1/
      data ivol /-1/
!
!     gradients and volume: constant over a straight tet, so one evaluation at
!     the centroid is exact. iflag=3 -- iflag=2 returns before the inverse
!     Jacobian is applied and yields LOCAL derivatives.
!
      call shape4tet(0.25d0,0.25d0,0.25d0,xl,xsj,shp,3)
      vol=xsj/6.d0
      do i=1,4
        do j=1,3
          g(j,i)=shp(j,i)
        enddo
      enddo
!
!     deviatoric linear block
!
      do i=1,4
        do c=1,3
          ii=3*(i-1)+c
          do j=1,4
            do d=1,3
              jj=3*(j-1)+d
              fac=0.d0
              if(c.eq.d) fac=g(1,i)*g(1,j)+g(2,i)*g(2,j)+g(3,i)*g(3,j)
              aull(ii,jj)=vol*(um*(fac+g(d,i)*g(c,j))
     &             -2.d0*um/3.d0*g(c,i)*g(d,j))
            enddo
          enddo
        enddo
      enddo
!
!     A_lb is identically zero; kept in the interface so the recovery routine
!     can use the same condensation algebra without special-casing.
!
      do i=1,12
        do j=1,3
          alb(i,j)=0.d0
        enddo
      enddo
!
!     bubble block. gm(c,d) = sum_i g_i[c] g_i[d]
!
      cbb=256.d0*256.d0*vol/15120.d0
      do c=1,3
        do d=1,3
          gm(c,d)=g(c,1)*g(d,1)+g(c,2)*g(d,2)+g(c,3)*g(d,3)
     &           +g(c,4)*g(d,4)
        enddo
      enddo
      gbb=cbb*(gm(1,1)+gm(2,2)+gm(3,3))
      do c=1,3
        do d=1,3
          abb(c,d)=um/3.d0*cbb*gm(c,d)
        enddo
        abb(c,c)=abb(c,c)+um*gbb
      enddo
!
!     THE BUBBLE'S VOLUMETRIC ENERGY.
!
!     abb above is deviatoric only, which leaves the bubble able to dilate
!     inside the element at a cost of order mu. The P1 pressure resists only
!     the projection int q div(b); everything orthogonal to that four-
!     dimensional space is unresisted, so with K/mu = 5000 the element acquires
!     a dilatation channel 5000x softer than the physics it bypasses. That is
!     the 440x excess compliance measured on one element, and it is why
!     unscaled MINI reads within 4.5% of a DRAINED cell.
!
!     Scaling the condensed stabilisation cannot fix this -- CCX_U4_STAB=CAPPED
!     tried, and bought under-stabilisation instead (pressure oscillation 1.034
!     against 0.754 for MINI, above every displacement control on the mesh).
!     The missing term is not a constant, it is energy.
!
!     STIFFB adds the bubble's own volumetric stiffness, lambda (div b)^2 with
!     lambda = K - 2mu/3, so abb scales with K rather than mu and the condensed
!     term S = B_b abb^-1 B_b^T drops to the order of C. It slightly double
!     counts the P1-visible part of the bubble divergence, which errs stiff --
!     the safe direction, and small because most of div(b) is orthogonal to P1.
!
!     int db_c/dx_c db_d/dx_d over the tet is cbb*gm(c,d), the same integral
!     already formed above, so the volumetric term costs nothing extra.
!
      if(ivol.lt.0) then
        spaxenv=' '
        call getenv('CCX_U4_STAB',spaxenv)
        if(spaxenv(1:6).eq.'STIFFB') then
          ivol=1
        else
          ivol=0
        endif
      endif
      if(ivol.eq.1) then
        xlam=xk-2.d0*um/3.d0
        do c=1,3
          do d=1,3
            abb(c,d)=abb(c,d)+xlam*cbb*gm(c,d)
          enddo
        enddo
      endif
!
!     couplings and compressibility
!
      cbl=256.d0*vol/840.d0
      do i=1,4
        do j=1,4
          do d=1,3
            bl(i,3*(j-1)+d)=vol/4.d0*g(d,j)
          enddo
        enddo
        do d=1,3
          bb(i,d)=-cbl*g(d,i)
        enddo
        do j=1,4
          cc(i,j)=vol/20.d0/xk
        enddo
        cc(i,i)=vol/10.d0/xk
      enddo
!
!     abi = abb^-1
!
      det=abb(1,1)*(abb(2,2)*abb(3,3)-abb(2,3)*abb(3,2))
     &   -abb(1,2)*(abb(2,1)*abb(3,3)-abb(2,3)*abb(3,1))
     &   +abb(1,3)*(abb(2,1)*abb(3,2)-abb(2,2)*abb(3,1))
      if(dabs(det).lt.1.d-30) then
        write(*,*) '*ERROR in u4mat: singular bubble block, element',
     &       nelem
        write(*,*) '       mu,K,vol =',um,xk,vol
        call exit(201)
      endif
      abi(1,1)=(abb(2,2)*abb(3,3)-abb(2,3)*abb(3,2))/det
      abi(1,2)=(abb(1,3)*abb(3,2)-abb(1,2)*abb(3,3))/det
      abi(1,3)=(abb(1,2)*abb(2,3)-abb(1,3)*abb(2,2))/det
      abi(2,1)=(abb(2,3)*abb(3,1)-abb(2,1)*abb(3,3))/det
      abi(2,2)=(abb(1,1)*abb(3,3)-abb(1,3)*abb(3,1))/det
      abi(2,3)=(abb(1,3)*abb(2,1)-abb(1,1)*abb(2,3))/det
      abi(3,1)=(abb(2,1)*abb(3,2)-abb(2,2)*abb(3,1))/det
      abi(3,2)=(abb(1,2)*abb(3,1)-abb(1,1)*abb(3,2))/det
      abi(3,3)=(abb(1,1)*abb(2,2)-abb(1,2)*abb(2,1))/det
!
!     the only bubble effect: C' = C + theta * B_b abb^-1 B_b^T
!
!     SCALING THE STABILISATION FOR NEAR-INCOMPRESSIBLE ELASTICITY.
!
!     MINI's bubble stabilisation is derived for Stokes, where the physical
!     compressibility term is absent. Its size relative to that term is
!
!         S/C  ~  0.05 K/mu        -- independent of h
!
!     so for the brine, K/mu = 5000 and S is a few hundred times the physics.
!     Measured on one element: (C+S)/C = 439. The consequence is not subtle --
!     on a layered cell the undrained answer (K = 2.2 GPa) lands within 4.5% of
!     a genuinely DRAINED one (K = 2.2 MPa). The brine's bulk stiffness, a
!     factor of a thousand, is simply gone.
!
!     The stabilisation is still correct in the sense that matters structurally
!     -- it annihilates uniform pressure exactly (verified to 1e-15), so it
!     penalises only pressure gradients. It is the magnitude that is wrong for
!     this material contrast.
!
!     CCX_U4_STAB selects the remedy:
!
!       MINI   (default) unscaled, theta = 1. Correct for Stokes, and what the
!              textbook element is. Keep as the reference.
!       CAPPED theta = tr(C)/(tr(C)+tr(S)), so the stabilisation can never
!              exceed the physical compressibility it is perturbing. When the
!              material is compressible (S << C) theta -> 1 and this IS MINI;
!              when S >> C it caps at S_eff ~ C. Parameter-free, and it relies
!              on the pressure mass term -- which is nonzero here because the
!              brine is nu = 0.49993, not exactly 1/2 -- supplying the rest of
!              the stability itself.
!
!     CAPPED is a deliberate departure from textbook MINI and is validated
!     against Abaqus C3D4H, not asserted. See ../elements_ccx/README.md.
!
      if(istab.lt.0) then
        spaxenv=' '
        call getenv('CCX_U4_STAB',spaxenv)
        if(spaxenv(1:6).eq.'CAPPED') then
          istab=1
        else
          istab=0
        endif
      endif
!
      do i=1,3
        do j=1,4
          tmp2(i,j)=abi(i,1)*bb(j,1)+abi(i,2)*bb(j,2)+abi(i,3)*bb(j,3)
        enddo
      enddo
      trc=0.d0
      trs=0.d0
      do i=1,4
        trc=trc+cc(i,i)
        do j=1,3
          stab(i,j)=bb(i,j)
        enddo
      enddo
      do i=1,4
        strow(i)=bb(i,1)*tmp2(1,i)+bb(i,2)*tmp2(2,i)+bb(i,3)*tmp2(3,i)
        trs=trs+strow(i)
      enddo
      theta=1.d0
      if((istab.eq.1).and.(trc+trs.gt.0.d0)) theta=trc/(trc+trs)
      do i=1,4
        do j=1,4
          cc(i,j)=cc(i,j)+theta*(bb(i,1)*tmp2(1,j)+bb(i,2)*tmp2(2,j)
     &         +bb(i,3)*tmp2(3,j))
        enddo
      enddo
!
      return
      end
