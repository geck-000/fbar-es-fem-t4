!
!     SPAX: stress and internal-force recovery for the U4 mixed tetrahedron.
!
!     Deliberately re-derives the element operator through u4mat rather than
!     recomputing forces from a standard B. That is the mistake that leaves
!     0002-bbar-mean-dilatation.patch reporting an equilibrium gap of 1.0:
!     patch the stiffness but recover the forces with the unpatched operator
!     and the reaction cross-check -- the only reference-free convergence
!     evidence this repository has -- stops meaning anything. Here the same
!     matrix produces both, so equilibrium_gap stays a real check.
!
      subroutine resultsmech_u4(co,kon,ipkon,lakon,v,stx,elcon,nelcon,
     &     rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,
     &     ntmat_,t0,t1,ithermal,iperturb,fn,iout,qa,vold,nmethod,
     &     dtime,plicon,nplicon,plkcon,nplkcon,xstiff,npmat_,matname,
     &     mi,ncmat_,calcul_fn,calcul_qa,nal,nelem)
!
      implicit none
!
      character*8 lakon(*)
      character*80 matname(*),amat
!
      integer kon(*),ipkon(*),konl(4),mi(*),ielmat(mi(3),*),
     &     ielorien(mi(3),*),norien,ntmat_,ithermal(*),iperturb(*),
     &     iout,nmethod,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),
     &     npmat_,ncmat_,calcul_fn,calcul_qa,nal,nelem,nelcon(2,*),
     &     nrhcon(*),nalcon(2,*),indexe,i,j,c,imat,kk
!
      real*8 co(3,*),v(0:mi(2),*),stx(6,mi(1),*),
     &     elcon(0:ncmat_,ntmat_,*),
     &     rhcon(0:1,ntmat_,*),alcon(0:6,ntmat_,*),alzero(*),orab(7,*),
     &     t0(*),t1(*),fn(0:mi(2),*),qa(*),vold(0:mi(2),*),dtime,
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),stiff(21),elconloc(ncmat_),eth(6),
     &     plconloc(802),coords(3),rho,t0l,t1l
!
      real*8 xl(3,4),um,xk,e,un,aull(12,12),bl(4,12),cc(4,4),abi(3,3),
     &     alb(12,3),bb(4,3),ul(12),pl(4),ub(3),fu(12),fp(4),
     &     shp(4,4),xsj,gl(3,4),eps(3,3),
     &     tr,pg,sig(3,3)
!
!
      indexe=ipkon(nelem)
      do i=1,4
        konl(i)=kon(indexe+i)
        do j=1,3
          xl(j,i)=co(j,konl(i))
        enddo
      enddo
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
!     element solution vector
!
      do i=1,4
        do c=1,3
          ul(3*(i-1)+c)=v(c,konl(i))
        enddo
        pl(i)=v(4,konl(i))
      enddo
!
!     bubble amplitude, from the same condensation the stiffness used
!
      do i=1,3
        ub(i)=0.d0
        do j=1,12
          ub(i)=ub(i)-abi(i,1)*alb(j,1)*ul(j)-abi(i,2)*alb(j,2)*ul(j)
     &         -abi(i,3)*alb(j,3)*ul(j)
        enddo
        do j=1,4
          ub(i)=ub(i)-abi(i,1)*bb(j,1)*pl(j)-abi(i,2)*bb(j,2)*pl(j)
     &         -abi(i,3)*bb(j,3)*pl(j)
        enddo
      enddo
!
!     internal forces: exactly the condensed operator times the solution
!
      do i=1,12
        fu(i)=0.d0
        do j=1,12
          fu(i)=fu(i)+aull(i,j)*ul(j)
        enddo
        do j=1,4
          fu(i)=fu(i)+bl(j,i)*pl(j)
        enddo
      enddo
      do i=1,4
        fp(i)=0.d0
        do j=1,12
          fp(i)=fp(i)+bl(i,j)*ul(j)
        enddo
        do j=1,4
          fp(i)=fp(i)-cc(i,j)*pl(j)
        enddo
      enddo
!
      if(calcul_fn.eq.1) then
        do i=1,4
          do c=1,3
            fn(c,konl(i))=fn(c,konl(i))+fu(3*(i-1)+c)
          enddo
          fn(4,konl(i))=fn(4,konl(i))+fp(i)
        enddo
      endif
!
!     Stress. With closed-form integration there is one output point and it
!     carries the EXACT element volume average, which is the quantity the
!     homogenisation reads. That is exact rather than approximate because the
!     bubble contributes nothing to the mean strain -- int grad(b) dV = 0 --
!     so the average strain is the linear part alone, and the average pressure
!     is (p1+p2+p3+p4)/4 since int L_i dV = V/4 for each.
!
!     This also removes a trap. ccx's .dat reader collapses integration points
!     by ARITHMETIC mean, which equals the volume average only for equal
!     quadrature weights. The 15-point tet rule does not have equal weights, so
!     the previous version was feeding the homogenisation a slightly wrong
!     element average. One exact point cannot be wrong that way.
!
      call shape4tet(0.25d0,0.25d0,0.25d0,xl,xsj,shp,3)
      do i=1,4
        do j=1,3
          gl(j,i)=shp(j,i)
        enddo
      enddo
!
      do i=1,3
        do j=1,3
          eps(i,j)=0.d0
        enddo
      enddo
      do i=1,4
        do c=1,3
          do j=1,3
            eps(c,j)=eps(c,j)+0.5d0*ul(3*(i-1)+c)*gl(j,i)
            eps(j,c)=eps(j,c)+0.5d0*ul(3*(i-1)+c)*gl(j,i)
          enddo
        enddo
      enddo
!
      pg=(pl(1)+pl(2)+pl(3)+pl(4))/4.d0
      tr=(eps(1,1)+eps(2,2)+eps(3,3))/3.d0
      do i=1,3
        do j=1,3
          sig(i,j)=2.d0*um*eps(i,j)
        enddo
        sig(i,i)=sig(i,i)-2.d0*um*tr+pg
      enddo
!
      stx(1,1,nelem)=sig(1,1)
      stx(2,1,nelem)=sig(2,2)
      stx(3,1,nelem)=sig(3,3)
      stx(4,1,nelem)=sig(1,2)
      stx(5,1,nelem)=sig(1,3)
      stx(6,1,nelem)=sig(2,3)
!
      nal=nal+4
!
      return
      end
