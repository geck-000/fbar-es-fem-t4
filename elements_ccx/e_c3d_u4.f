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
      subroutine e_c3d_u4(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &     ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &     alzero,ielmat,ielorien,norien,orab,ntmat_,
     &     t0,t1,ithermal,vold,iperturb,nelemload,
     &     sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &     nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &     matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &     ttime,time,istep,iinc,coriolis,xloadold,reltime,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &     ne0,ipkon,thicke,integerglob,doubleglob,tieset,istartset,
     &     iendset,ialset,ntie,nasym,ielprop,prop)
!
      implicit none
!
      integer mass(*),stiffness,buckling,rhsi,coriolis
!
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*),amat
      character*81 tieset(3,*)
!
      integer konl(4),nelemload(2,*),nbody,nelem,mi(*),kon(*),
     &     ielprop(*),mattyp,ithermal(*),iperturb(*),nload,
     &     idist,i,j,c,d,ir,ic,nmethod,kk,nelcon(2,*),nrhcon(*),
     &     nalcon(2,*),ielmat(mi(3),*),ielorien(mi(3),*),ipkon(*),
     &     indexe,ntmat_,norien,ihyper,iexpl,kode,imat,iorien,
     &     istiff,ncmat_,intscheme,istep,iinc,ipompc(*),nodempc(3,*),
     &     nmpc,ikmpc(*),ilmpc(*),ne0,istartset(*),iendset(*),
     &     ialset(*),ntie,integerglob(*),nasym,nplicon(0:ntmat_,*),
     &     nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),xl(3,4),veold(0:mi(2),*),rho,s(60,60),bodyfx(3),
     &     ff(60),elconloc(ncmat_),coords(3),p1(3),p2(3),eth(6),
     &     rhcon(0:1,ntmat_,*),reltime,prop(*),alcon(0:6,ntmat_,*),
     &     alzero(*),orab(7,*),t0(*),t1(*),xloadold(2,*),
     &     vold(0:mi(2),*),xload(2,*),omx,e,un,um,xk,t0l,t1l,
     &     coefmpc(*),sm(60,60),sti(6,mi(1),*),stx(6,mi(1),*),
     &     thicke(mi(3),*),doubleglob(*),
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),plconloc(802),dtime,ttime,time,
     &     elcon(0:ncmat_,ntmat_,*),stiff(21)
!
      real*8 aull(12,12),bl(4,12),cc(4,4),abi(3,3),alb(12,3),bb(4,3)
!
      indexe=ipkon(nelem)
      do i=1,4
        konl(i)=kon(indexe+i)
        do j=1,3
          xl(j,i)=co(j,konl(i))
        enddo
      enddo
!
      do i=1,60
        ff(i)=0.d0
        do j=1,60
          s(i,j)=0.d0
          sm(i,j)=0.d0
        enddo
      enddo
      if(stiffness.eq.0) return
!
!     isotropic linear elasticity only. The mixed form needs mu and K; lambda
!     never appears, which is the point -- it is the term that blows up as
!     nu approaches 1/2.
!
      imat=ielmat(1,nelem)
!
!     E and nu straight from the material card. materialdata_me caches into
!     xstiff during assembly and does not re-derive the constants in the
!     results pass, so calling it there returns zeros; reading elcon keeps the
!     stiffness and the recovery on identical constants with no hidden state.
!     Isothermal isotropic elasticity only -- elcon(0,.,.) is the temperature,
!     elcon(1,.,.) is E and elcon(2,.,.) is nu.
!
      if(nelcon(1,imat).ne.2) then
        write(*,*) '*ERROR in U4: element',nelem,' needs an isotropic'
        write(*,*) '       *ELASTIC card (2 constants)'
        call exit(201)
      endif
      e=elcon(1,1,imat)
      un=elcon(2,1,imat)
      um=e/(2.d0*(1.d0+un))
      xk=e/(3.d0*(1.d0-2.d0*un))
!
      call u4mat(xl,um,xk,aull,bl,cc,abi,alb,bb,nelem)
!
!     scatter into the node-major 4-dof-per-node element matrix:
!     dof 1-3 displacement, dof 4 pressure. The pressure block is negative,
!     so the assembled system is symmetric INDEFINITE -- SPOOLES only.
!
      do i=1,4
        do c=1,3
          ir=4*(i-1)+c
          do j=1,4
            do d=1,3
              s(ir,4*(j-1)+d)=aull(3*(i-1)+c,3*(j-1)+d)
            enddo
            s(ir,4*(j-1)+4)=bl(j,3*(i-1)+c)
          enddo
        enddo
        ir=4*(i-1)+4
        do j=1,4
          do d=1,3
            s(ir,4*(j-1)+d)=bl(i,3*(j-1)+d)
          enddo
          s(ir,4*(j-1)+4)=-cc(i,j)
        enddo
      enddo
!
      return
      end
