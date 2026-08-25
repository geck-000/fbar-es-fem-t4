!
!     SPAX: U6 -- the nodal-averaged B-bar volumetric patch.
!
!     One element per node of the U5 phase. Its connectivity is that node's
!     1-ring, with THE PATCH CENTRE FIRST, and its stiffness is the rank-one
!
!         K_vol^a = K V_a b^a (b^a)^T ,   theta_a = b^a . u
!
!         V_a   = sum over U5 elements at a of V_e/4
!         b_n^a = (1/V_a) sum over those elements of (V_e/4) grad(L_n)|_e
!
!     WHY THIS IS AN ELEMENT AND NOT *EQUATION + SPRING1. The first attempt
!     carried theta as an extra dof tied by an *EQUATION and given energy by a
!     grounded spring. Every piece of that was verified exact in isolation --
!     including a two-node case, theta = 3u with k = 1 and F = 1, giving
!     u = 1/9 to every digit -- and a single patch on a real mesh reproduced an
!     independent Python assembly exactly. But adjacent patches overlap by
!     construction, so each mesh dof appears as an independent term in ~24
!     equations, and with the full set ccx returned a confined-compression
!     reaction 1.55x the closed-form answer. Delivered as an element the MPC
!     machinery is bypassed entirely.
!
!     The patch geometry is not in the deck: b^a depends on which subsets of the
!     ring form tets, which a node list does not determine. So the element finds
!     its own tets through a node->element map over the U5 elements, built once
!     and cached. The build is wrapped in a critical section because mafillsm
!     assembles under OpenMP.
!
!         *USER ELEMENT,TYPE=U6,NODES=<ring size>,INTEGRATIONPOINTS=1,MAXDOF=3
!
!     Keep U6 out of any *EL PRINT set: it has no material volume, and
!     printoutelem.f would try to integrate it with a shape function chosen
!     from its node count.
!
      subroutine e_c3d_u6(co,kon,lakonl,s,ff,nelem,elcon,nelcon,
     &     ielmat,mi,ncmat_,ntmat_,ipkon,lakon,ne,stiffness)
!
      implicit none
!
      character*8 lakonl,lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ielmat(mi(3),*),ncmat_,ntmat_,
     &     nelcon(2,*),nelem,ne,stiffness,i,j,c,d,ii,jj,nope,indexe,
     &     acen,imat,konl(255)
      real*8 co(3,*),s(150,150),ff(150),elcon(0:ncmat_,ntmat_,*),
     &     e,un,xk,va,kva,bb(3,255)
!
      nope=ichar(lakonl(8:8))
!
!     SPAX: hard capacity guard.  The patch spans the centre node's whole
!     1-ring, so nope follows mesh valence.  s is dimensioned (150,150)
!     in mafillsm and along the whole user-element path; overrunning it
!     used to corrupt the assembly silently.  Stop loudly instead.
!
      if(3*nope.gt.150) then
        write(*,*) '*ERROR in e_c3d_u6: patch element ',nelem
        write(*,*) '       spans ',nope,' nodes (',3*nope,' DOF).'
        write(*,*) '       The element matrix holds 150 DOF (50 nodes).'
        write(*,*) '       Raise nduser in mafillsm.f and the s/sm/ff'
        write(*,*) '       dimensions in the e_c3d_u* family together.'
        call exit(201)
      endif
      indexe=ipkon(nelem)
      do i=1,nope
        konl(i)=kon(indexe+i)
      enddo
      acen=konl(1)
!
      do i=1,60
        ff(i)=0.d0
        do j=1,60
          s(i,j)=0.d0
        enddo
      enddo
      if(stiffness.eq.0) return
!
      call u6patch(co,kon,ipkon,lakon,ne,konl,nope,acen,va,kva,bb,
     &     nelem,ielmat,elcon,nelcon,mi,ncmat_,ntmat_,
     &     ielmat(1,nelem))
      if(va.le.0.d0) return
!
      do i=1,nope
        do c=1,3
          ii=3*(i-1)+c
          do j=1,nope
            do d=1,3
              jj=3*(j-1)+d
              s(ii,jj)=kva*bb(c,i)*bb(d,j)
            enddo
          enddo
        enddo
      enddo
!
      return
      end
