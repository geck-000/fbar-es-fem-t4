!
!     SPAX: U5 -- linear tetrahedron carrying the DEVIATORIC stiffness only.
!
!     Half of the nodal-averaged B-bar pair. The volumetric response is removed
!     from the element and rebuilt from nodal averages outside it, so that the
!     bulk modulus never enters an element-local operator and cannot lock:
!
!        K = K_dev  (this element, per tet)
!          + K_vol  (per node: theta_a tied to the surrounding displacements by
!                    an *EQUATION, given energy by a grounded SPRING1 of
!                    stiffness K*V_a -- see nodalbbar.py)
!
!     Why this escapes what killed U4. There the stabilisation and the spurious
!     compliance were the same term, S = B_b A_bb^-1 B_b^T, so shrinking one
!     shrank the other and no setting was both stable and stiff enough. Here
!     stability comes from averaging over the node patch -- a geometric
!     construction with no 1/mu scaling anywhere -- and there is no pressure
!     unknown at all, so the assembled system stays symmetric POSITIVE
!     DEFINITE and incomplete-Cholesky PCG applies.
!
!         *USER ELEMENT,TYPE=U5,NODES=4,INTEGRATIONPOINTS=1,MAXDOF=3
!
!     STRESS OUTPUT IS DELIBERATELY NOT THE ELEMENT'S OWN CONSTITUTIVE
!     RESPONSE. resultsmech_u5 reports sigma = 2 mu dev(eps) + K div(u) using
!     the ELEMENT's divergence, not the nodal average, even though the
!     stiffness uses neither. That is correct for the homogenisation, and the
!     reason is an identity: theta_a is the V-weighted mean of the surrounding
!     element divergences, so
!
!         sum_a V_a theta_a = sum_a sum_(e at a) (V_e/4) div(u)|_e
!                           = sum_e V_e div(u)|_e
!
!     -- every element feeds four nodes. The volume-INTEGRATED volumetric
!     stress is therefore identical computed element-wise or nodally, so the
!     volume-averaged stress the homogenisation reads is exact even though the
!     pointwise split is not. Internal forces stay deviatoric, matching this
!     element's stiffness; the springs supply the rest of the reaction, so
!     equilibrium_gap remains a real check.
!
      subroutine e_c3d_u5(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
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
      character*80 matname(*)
      character*81 tieset(3,*)
!
      integer konl(4),nelemload(2,*),nbody,nelem,mi(*),kon(*),
     &     ielprop(*),ithermal(*),iperturb(*),nload,idist,i,j,c,d,
     &     ir,nmethod,nelcon(2,*),nrhcon(*),nalcon(2,*),
     &     ielmat(mi(3),*),ielorien(mi(3),*),ipkon(*),indexe,ntmat_,
     &     norien,iexpl,imat,ncmat_,intscheme,istep,iinc,ipompc(*),
     &     nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),ne0,istartset(*),
     &     iendset(*),ialset(*),ntie,integerglob(*),nasym,
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),xl(3,4),veold(0:mi(2),*),s(150,150),bodyfx(3),
     &     ff(150),p1(3),p2(3),rhcon(0:1,ntmat_,*),reltime,prop(*),
     &     alcon(0:6,ntmat_,*),alzero(*),orab(7,*),t0(*),t1(*),
     &     xloadold(2,*),vold(0:mi(2),*),xload(2,*),omx,e,un,um,
     &     coefmpc(*),sm(150,150),sti(6,mi(1),*),stx(6,mi(1),*),
     &     thicke(mi(3),*),doubleglob(*),
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),dtime,ttime,time,
     &     elcon(0:ncmat_,ntmat_,*)
!
      real*8 shp(4,4),xsj,vol,g(3,4),fac,gh
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
      imat=ielmat(1,nelem)
      if(nelcon(1,imat).ne.2) then
        write(*,*) '*ERROR in e_c3d_u5: element',nelem,' needs an'
        write(*,*) '       isotropic *ELASTIC card (2 constants)'
        call exit(201)
      endif
      e=elcon(1,1,imat)
      un=elcon(2,1,imat)
      um=e/(2.d0*(1.d0+un))
!
!     gradients and volume are constant over a straight tet; iflag=3, because
!     iflag=2 returns before the inverse Jacobian is applied.
!
      call shape4tet(0.25d0,0.25d0,0.25d0,xl,xsj,shp,3)
      vol=xsj/6.d0
      do i=1,4
        do j=1,3
          g(j,i)=shp(j,i)
        enddo
      enddo
!
!     2 mu dev(eps):dev(eps). For scalar shape functions with gradients g_i,
!     g_j and components c,d:
!         eps:eps = (delta_cd (g_i.g_j) + g_i[d] g_j[c])/2
!         div div = g_i[c] g_j[d]
!
      do i=1,4
        do c=1,3
          ir=3*(i-1)+c
          do j=1,4
            do d=1,3
              fac=0.d0
              if(c.eq.d) fac=g(1,i)*g(1,j)+g(2,i)*g(2,j)+g(3,i)*g(3,j)
              s(ir,3*(j-1)+d)=vol*(um*(fac+g(d,i)*g(c,j))
     &             -2.d0*um/3.d0*g(c,i)*g(d,j))
            enddo
          enddo
        enddo
      enddo
!
      return
      end
