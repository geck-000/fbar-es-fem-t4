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
      subroutine u6patch(co,kon,ipkon,lakon,ne,konl,nope,acen,
     &     va,kva,bb,nelem,ielmat,elcon,nelcon,mi,ncmat_,ntmat_,
     &     imatf)
!
!     Shared by e_c3d_u6 (stiffness) and resultsmech_u6 (internal forces), so
!     the two cannot drift apart -- the failure that leaves the B-bar patch
!     reporting an equilibrium_gap of 1.0.
!
      implicit none
!
      character*8 lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ne,konl(*),nope,acen,nelem,
     &     i,j,k,c,ie,ipos,n4(4),ielmat(mi(3),*),nelcon(2,*),
     &     ncmat_,ntmat_,imat,imatf
      real*8 co(3,*),va,kva,bb(3,*),xl(3,4),shp(4,4),xsj,vol,w,
     &     elcon(0:ncmat_,ntmat_,*),ek,eun,xke
!
      integer mapdone,maxn,nlen
      integer, allocatable, save :: nstart(:),nlist(:)
      save mapdone,maxn
      data mapdone /0/
!
!
!     The map is built once and then read by every thread.  The test MUST
!     stay inside the critical region: with the cheap double-checked
!     `if(mapdone.eq.0)` outside it, a thread could see mapdone=1 before
!     the nstart/nlist writes were visible to it and read a half-built
!     map, which silently corrupted va/kva/bb for that patch -- in the
!     stiffness and in the force recovery alike.  Measured on the sphere
!     cell: equilibrium_gap 8.7e-08 at 1-2 threads, 7.3e-01 at 8.
!     One lock acquire per patch element is nothing next to that.
!
      call u6lock()
      if(mapdone.eq.0) then
          maxn=0
          nlen=0
          do i=1,ne
            if(ipkon(i).lt.0) cycle
            if(lakon(i)(1:2).ne.'U5') cycle
            do j=1,4
              k=kon(ipkon(i)+j)
              if(k.gt.maxn) maxn=k
              nlen=nlen+1
            enddo
          enddo
          allocate(nstart(maxn+2))
          allocate(nlist(max(nlen,1)))
          do i=1,maxn+2
            nstart(i)=0
          enddo
          do i=1,ne
            if(ipkon(i).lt.0) cycle
            if(lakon(i)(1:2).ne.'U5') cycle
            do j=1,4
              nstart(kon(ipkon(i)+j)+1)=nstart(kon(ipkon(i)+j)+1)+1
            enddo
          enddo
          do i=2,maxn+2
            nstart(i)=nstart(i)+nstart(i-1)
          enddo
          do i=1,ne
            if(ipkon(i).lt.0) cycle
            if(lakon(i)(1:2).ne.'U5') cycle
            do j=1,4
              k=kon(ipkon(i)+j)
              nstart(k)=nstart(k)+1
              nlist(nstart(k))=i
            enddo
          enddo
          do i=maxn+1,2,-1
            nstart(i)=nstart(i-1)
          enddo
          nstart(1)=0
        mapdone=1
      endif
      call u6unlock()
!
      va=0.d0
      kva=0.d0
      do i=1,nope
        do c=1,3
          bb(c,i)=0.d0
        enddo
      enddo
      if(acen.gt.maxn) then
        write(*,*) '*ERROR in u6patch: element',nelem,' centre node',
     &       acen,' is in no U5 element'
        call exit(201)
      endif
!
      do ie=nstart(acen)+1,nstart(acen+1)
        i=nlist(ie)
!
!       PER-MATERIAL PATCHES. Only elements of the patch's own material
!       contribute. Averaging theta across a phase boundary breaks the identity
!       the stress recovery depends on -- sum_a V_a theta_a = sum_e V_e div(u)|_e
!       holds for ONE K, and across the brine/ice interface it mixes
!       divergences over a 1000x modulus contrast. Measured: equilibrium_gap
!       2.5e-3 at L_mesh=0.0240 and 2.6e-1 at 0.0120 with a shared patch.
!       Splitting by material restores it exactly, and still gives the ice its
!       own averaging -- which it needs, because locking follows isochoric
!       DEFORMATION, not the material's own nu, and ice next to an
!       incompressible slab is forced to deform isochorically.
!
        if(ielmat(1,i).ne.imatf) cycle
        do j=1,4
          n4(j)=kon(ipkon(i)+j)
          do k=1,3
            xl(k,j)=co(k,n4(j))
          enddo
        enddo
        call shape4tet(0.25d0,0.25d0,0.25d0,xl,xsj,shp,3)
        vol=dabs(xsj)/6.d0
        w=vol/4.d0
        va=va+w
!
!       K-weighted volume. Replacing the element divergence by the nodal one in
!       sum_e (1/2) K_e V_e theta_e^2 gives (1/2) theta_a^2 sum_e K_e V_e/4, so
!       a patch spanning two phases carries sum K_e V_e/4, NOT K times V_a.
!       That matters here: an interface patch mixes K = 2.2 GPa with 9.25 GPa.
!
        imat=ielmat(1,i)
        ek=elcon(1,1,imat)
        eun=elcon(2,1,imat)
        xke=ek/(3.d0*(1.d0-2.d0*eun))
        kva=kva+xke*w
        do j=1,4
          ipos=0
          do k=1,nope
            if(konl(k).eq.n4(j)) then
              ipos=k
              exit
            endif
          enddo
          if(ipos.eq.0) then
            write(*,*) '*ERROR in u6patch: element',nelem,' node',
     &           n4(j),' not in its own connectivity'
            call exit(201)
          endif
          do c=1,3
            bb(c,ipos)=bb(c,ipos)+w*shp(c,j)
          enddo
        enddo
      enddo
!
      if(va.gt.0.d0) then
        do i=1,nope
          do c=1,3
            bb(c,i)=bb(c,i)/va
          enddo
        enddo
      endif
!
      return
      end
