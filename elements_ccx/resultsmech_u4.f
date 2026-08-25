!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2025 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine resultsmech_u4(fn,calcul_fn,calcul_qa,mi,nelem)
!
!     SPAX: U4 -- the null base tetrahedron of the F-barES-FEM-T4 element.
!
!     U4 contributes no internal force, so its nodal-force block is zeroed;
!     no stress is written.
!
      implicit none
!
      integer calcul_fn,calcul_qa,mi(*),nelem,i
!
      real*8 fn(0:mi(2),*)
!
      if(calcul_fn.eq.1) then
         do i=0,mi(2)
            fn(i,nelem)=0.d0
         enddo
      endif
!
      return
      end
