!
!     SPAX: stress and internal-force recovery for the U5 deviatoric tet.
!
!     The two outputs deliberately use DIFFERENT operators, and the reason is
!     an identity rather than an oversight:
!
!     fn  = deviatoric only, matching this element's stiffness exactly. The
!           volumetric reaction is carried by the nodal springs and their
!           *EQUATIONs, so the reference-point reaction is still complete and
!           equilibrium_gap remains a real check.
!
!     stx = 2 mu dev(eps) + K div(u) with the ELEMENT's own divergence. The
!           volume-INTEGRATED volumetric stress is the same whether summed over
!           elements or over nodal patches, because theta_a is the V-weighted
!           mean of the surrounding element divergences and each element feeds
!           exactly four nodes:
!
!               sum_a V_a theta_a = sum_e V_e div(u)|_e
!
!           so the volume-averaged stress the homogenisation reads is exact,
!           even though the pointwise split between the two is not.
!
      subroutine resultsmech_u5(co,kon,ipkon,lakon,v,stx,elcon,nelcon,
     &     rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,
     &     ntmat_,t0,t1,ithermal,iperturb,fn,iout,qa,vold,nmethod,
     &     dtime,plicon,nplicon,plkcon,nplkcon,xstiff,npmat_,matname,
     &     mi,ncmat_,calcul_fn,calcul_qa,nal,nelem)
!
      implicit none
!
      character*8 lakon(*)
      character*80 matname(*)
!
      integer kon(*),ipkon(*),konl(4),mi(*),ielmat(mi(3),*),
     &     ielorien(mi(3),*),norien,ntmat_,ithermal(*),iperturb(*),
     &     iout,nmethod,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),
     &     npmat_,ncmat_,calcul_fn,calcul_qa,nal,nelem,nelcon(2,*),
     &     nrhcon(*),nalcon(2,*),indexe,i,j,c,d,imat
!
      real*8 co(3,*),v(0:mi(2),*),stx(6,mi(1),*),
     &     elcon(0:ncmat_,ntmat_,*),
     &     rhcon(0:1,ntmat_,*),alcon(0:6,ntmat_,*),alzero(*),orab(7,*),
     &     t0(*),t1(*),fn(0:mi(2),*),qa(*),vold(0:mi(2),*),dtime,
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*)
!
      real*8 xl(3,4),um,xk,e,un,shp(4,4),xsj,vol,g(3,4),ul(12),
     &     eps(3,3),tr,sig(3,3),fu(12),fac,gh,divu
!
!$omp critical(u5cnt)
!$omp end critical(u5cnt)
      indexe=ipkon(nelem)
      do i=1,4
        konl(i)=kon(indexe+i)
        do j=1,3
          xl(j,i)=co(j,konl(i))
        enddo
      enddo
!
      imat=ielmat(1,nelem)
      e=elcon(1,1,imat)
      un=elcon(2,1,imat)
      um=e/(2.d0*(1.d0+un))
      xk=e/(3.d0*(1.d0-2.d0*un))
!
      call shape4tet(0.25d0,0.25d0,0.25d0,xl,xsj,shp,3)
      vol=xsj/6.d0
      do i=1,4
        do j=1,3
          g(j,i)=shp(j,i)
        enddo
        do c=1,3
          ul(3*(i-1)+c)=v(c,konl(i))
        enddo
      enddo
!
!     strain, divergence
!
      do i=1,3
        do j=1,3
          eps(i,j)=0.d0
        enddo
      enddo
      do i=1,4
        do c=1,3
          do j=1,3
            eps(c,j)=eps(c,j)+0.5d0*ul(3*(i-1)+c)*g(j,i)
            eps(j,c)=eps(j,c)+0.5d0*ul(3*(i-1)+c)*g(j,i)
          enddo
        enddo
      enddo
      divu=eps(1,1)+eps(2,2)+eps(3,3)
      tr=divu/3.d0
!
!     internal forces: the deviatoric operator, matching the stiffness
!
      do i=1,4
        do c=1,3
          fu(3*(i-1)+c)=0.d0
          do j=1,4
            do d=1,3
              fac=0.d0
              if(c.eq.d) fac=g(1,i)*g(1,j)+g(2,i)*g(2,j)+g(3,i)*g(3,j)
              fu(3*(i-1)+c)=fu(3*(i-1)+c)
     &             +vol*(um*(fac+g(d,i)*g(c,j))
     &             -2.d0*um/3.d0*g(c,i)*g(d,j))*ul(3*(j-1)+d)
            enddo
          enddo
        enddo
      enddo
      if(calcul_fn.eq.1) then
        do i=1,4
          do c=1,3
            fn(c,konl(i))=fn(c,konl(i))+fu(3*(i-1)+c)
          enddo
        enddo
      endif
!
!     stress: deviatoric plus the element's own volumetric part -- see header
!
      do i=1,3
        do j=1,3
          sig(i,j)=2.d0*um*eps(i,j)
        enddo
        sig(i,i)=sig(i,i)-2.d0*um*tr+xk*divu
      enddo
      stx(1,1,nelem)=sig(1,1)
      stx(2,1,nelem)=sig(2,2)
      stx(3,1,nelem)=sig(3,3)
      stx(4,1,nelem)=sig(1,2)
      stx(5,1,nelem)=sig(1,3)
      stx(6,1,nelem)=sig(2,3)
!
      nal=nal+3*4
!
      return
      end
