c**********************************************************
      subroutine elem308 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     Material model by Sevan Goenezen, written by Sevan Goenezen and optimized 
c     by Jean-Francois.           
c     3D, linear tetrahedra, finite elasticity, incompressible
c     has deviatoric stress component
c     The strain energy density function is a modified Veronda-Westman
c      as in the paper by Sevan Goenezen and Assad Oberai (2010?)
c      W=\mu/(2.\gamma).(exp(\gamma.(J^(-2/3).I_1-3))-1)
c      S=-p.J.C^-1+\mu.J^(-2/3).(I_{dentity}-(I_1.C^-1)/3).exp(\gamma.J^(-2/3).I_1-3)
c      with this choice: trace(sigma) = -3 p
c     cite maniatty's paper on stabilized FEM
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem, s, ii
      parameter (mnode_elem = 4,mquad_elem=30)
      integer ielem, itask,k,ireg,iset
      integer l, ninte, iinte, inode, jnode,knode
      integer ievab, jevab, jdofn, idofn, i, j, q, r, t
      integer kk, tt
      integer tmp6, tmp7
      double precision xjaco, wtjac, temp, temp_der
      double precision deno, powe, alpha(2), beta(2), Treg(2)! variables for the regularization
      double precision h, tmp1, tmp2, tmp3, tmp4, tmp5, tauMult! variables for the stabilization
      double precision gamm, mu, Cdet, Inv1, Inv2, K1, K2(3,3), Fdet
      double precision shap(4,mnode_elem), pres(4)
      double precision sg(mquad_elem),tg(mquad_elem),zg(mquad_elem)
      double precision Ftmp(3,3), Ftmp_det, Sigma(3,3)
      double precision SecPK_grad(3,3), wg(mquad_elem)
      double precision Fdef(3,3), Ctens(3,3), dWdC(3,3), dJdC(3,3)
      double precision Cinv(3,3), SecPK(3,3), ident(3,3), Finv(3,3)
      double precision Ctang(3,3,3,3),d2JdC(3,3,3,3)
c JFD for debug      double precision Ctangt(3,3,3,3),d2JdCt(3,3,3,3)
      double precision Igeo(9,9),Lmat(3,3,3,3), Dgeomat(9,9)
      double precision Tp(4), FS(3,3), bb(mnode_elem,9,3)
      double precision dWdC_grad(3,3),udiff(3)
      double precision temp_dual(4*mnode_elem),prop_grad(2,4)
      double precision var_shape_prop(2,4)
      double precision dCinvdC(3,3,3,3)
      double precision dCinvdCt
c     Localized information
      integer nnode
      double precision pelem(*)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*)
      double precision uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      integer ndim! = ndime (local copy)
c-------------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask  
      
c     --------------------------------------------------------
c     Initialize the elemental parameters
c     --------------------------------------------------------
 1    elemvec_ndofn(1:nnode) = 4
      buildSymmetricMatrices=.false.
      return


c     --------------------------------------------------------
c     Read and store the material properties
c     --------------------------------------------------------
c     l          : no. of gauss pts/direction
c     ireg       : regularization type (0/1/2/3/4:none/H1/TVD/newTV/power)
c     alpha(1)   : regularization weight for the non-linear parameter
c     alpha(2)   : regularization weight for the shear modulus
c     beta(1)    : extra parameter for TVD - JFD ask Sevan for a better name
c     beta(2)    : extra parameter for TVD - JFD ask Sevan for a better name
c     Treg(1)    : second extra parameter
c     Treg(2)    : second extra parameter
c     tauMult    : stabilization factor
c     s          : stabilization type (0/1:off/on) - JFD use stype instead of s?
 2    read(iread,*) l,ireg,alpha(1),alpha(2),beta(1),beta(2),
     +     Treg(1),Treg(2),tauMult,s
      write(iwrit,205) l,ireg,alpha(1),alpha(2),beta(1),beta(2),
     +     Treg(1),Treg(2),tauMult,s
 205  format(/' finite elasticity (3D,4DOF) '/
     +       ' gauss pts/dir .........................',i12,/
     +       ' reg type(0/1/2/4:none/H1/TVD/log, etc).',i12,/
     +       ' regularization parameter 1 ............',1p,e16.4,/
     +       ' regularization parameter 2 ............',1p,e16.4,/
     +       ' extra parameter (TVD) .................',1p,e16.4,/! modify names here JFD
     +       ' extra parameter (TVD) .................',1p,e16.4,/! idem
     +       ' second extra parameter (TVD) ..........',1p,e16.4,/! modify names here JFD
     +       ' second extra parameter (TVD) ..........',1p,e16.4,/! idem
     +       ' stabilization factor ..................',1p,e16.4,/
     +       ' stabilization type (0/1:off/on) .......',i12)
c      nprop = 10
c     check for error in input file

      if (alpha(1).lt.0.0d0) then
        print*,'elem305.f: alpha(1) is less than 0.0d0: exiting'
        stop
      elseif (alpha(2).lt.0.0d0) then
        print*,'elem305.f: alpha(2) is less than 0.0d0: exiting'
        stop
      elseif (beta(1).le.0.0d0) then
        print*,'elem305.f: beta(1) is less or equal to 0.0d0: exiting'
        stop
      elseif (beta(2).le.0.0d0) then
        print*,'elem305.f: beta(2) is less or equal to 0.0d0: exiting'
        stop
      elseif (Treg(1).le.0.0d0) then
        print*,'elem305.f: Treg(1) is less or equal to 0.0d0: exiting'
        stop
      elseif (Treg(2).le.0.0d0) then
        print*,'elem305.f: Treg(2) is less or equal to 0.0d0: exiting'
        stop
      endif
     
c     note: ireg is checked in itask.eq.7 for keeping things simple
c     note: s is checked in itask.eq.3 for keeping things simple
      pelem(1) = dble(l)
      pelem(2) = dble(ireg)
      pelem(3) = alpha(1)
      pelem(4) = alpha(2)
      pelem(5) = beta(1)
      pelem(6) = beta(2)
      pelem(7) = Treg(1)
      pelem(8) = Treg(2)
      pelem(9) = tauMult
      pelem(10) = dble(s)
      return


c     --------------------------------------------------------
c     Build the elemental consistent tangent stiffness matrix (estif)
c      and the elemental RHS/residual (eforc)
c     --------------------------------------------------------
 3    l      = int(pelem(1))
      tauMult= pelem(9)
      s      = int(pelem(10))
      ndim = ndime

      
c     determine the characteristic length h of the element
c       (for tetrahedra: h is the length of the longest edge)
      tmp1=(xelem(1,1)-xelem(1,2))**2.0d0
      tmp2=(xelem(1,1)-xelem(1,3))**2.0d0
      tmp3=(xelem(1,1)-xelem(1,4))**2.0d0
      tmp4=(xelem(1,2)-xelem(1,3))**2.0d0
      tmp5=(xelem(1,2)-xelem(1,4))**2.0d0
      h=(xelem(1,3)-xelem(1,4))**2.0d0
      do i=2,ndim! 1, ndim
        tmp1 = tmp1 + (xelem(i,1)-xelem(i,2))**2.0d0
        tmp2 = tmp2 + (xelem(i,1)-xelem(i,3))**2.0d0
        tmp3 = tmp3 + (xelem(i,1)-xelem(i,4))**2.0d0
        tmp4 = tmp4 + (xelem(i,2)-xelem(i,3))**2.0d0
        tmp5 = tmp5 + (xelem(i,2)-xelem(i,4))**2.0d0
        h = h + (xelem(i,3)-xelem(i,4))**2.0d0
      enddo
      if (tmp1.gt.h) then
        h=tmp1
      endif
      if (tmp2.gt.h)then
        h=tmp2
      endif
      if (tmp3.gt.h) then
        h=tmp3
      endif
      if (tmp4.gt.h) then
        h=tmp4
      endif
      if (tmp5.gt.h) then
        h=tmp5
      endif
      h=sqrt(h)

c     initialize variables
      estif(:,:) = 0.0d0 ! elemental consistent tangent stiffness
      eforc(:) = 0.0d0! elemental RHS/residual
      ident(:,:) = 0.0d0! identity matrix
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0
      
      call gausstet(l,ninte,sg,tg,zg,wg)

c!$OMP PARALLEL DO
c!$OMP& PRIVATE(iinte,shap,xjaco,wtjac,mu,gamm,inode,temp,Fdef,j,
c!$OMP&  Fdet,Finv,pres,Ctens,r,Cdet,Cinv,Inv1,dJdC,tmp5,K1,tmp2,K2,tmp4,
c!$OMP&  tmp3,dWdC,SecPK,dCinvdCt,Ctang,d2JdC,Igeo,Lmat,q,k,t,Dgeomat,
c!$OMP&  bb,ievab,idofn,jevab,jnode,jdofn,Tp,FS)
c!$OMP& SHARED (ninte,ielem,sg,tg,zg,nnode,ndim,elnods,xelem,wg,
c!$OMP&  elemdata_nodal,s,tauMult,h,ident,uelem)
c!$OMP& REDUCTION(+:estif,eforc)! reduction on an array requires openmp2.5, reduction on allocatable openmp3.0
c!$acc region
      do iinte = 1,ninte! for all Gauss integration points
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,
     $               xjaco,.false.,nnode,ndim,elnods,xelem)
        wtjac = wg(iinte)*xjaco

        mu = 0.0d0! value of mu at the current Gauss point
        gamm = 0.0d0! value of gamma at the current Gauss point
        do inode = 1,nnode
          mu = mu + shap(4,inode)*elemdata_nodal(2,inode)
          gamm = gamm + shap(4,inode)*elemdata_nodal(1,inode)
        enddo

        if (s.eq.0) then! JFD find a better name than s? ask Sevan
          temp  = 0.0d0! JFD find a better name than temp? ask Sevan
        else if (s.eq.1) then
          temp = (0.5d0*tauMult*(h**2.0d0))/mu! find a better name than h? JFD ask sevan
        else
          write(iwrit,200) 
 200      format(4x,'elem305.f: Stabilization property must be 0 or 1')
          stop
        endif

c       compute the deformation gradient at Gauss Point
        Ftmp(:,:) = ident(:,:)
        do inode = 1,nnode
          do j = 1,ndim
            do i = 1,ndim
              Ftmp(i,j)=Ftmp(i,j)-uelem(i,inode)*shap(j,inode)
            end do
          end do
        end do
 
        Ftmp_det = Ftmp(1,1)*Ftmp(2,2)*Ftmp(3,3) +
     $             Ftmp(1,2)*Ftmp(2,3)*Ftmp(3,1) +
     $             Ftmp(1,3)*Ftmp(2,1)*Ftmp(3,2) -
     $             Ftmp(1,3)*Ftmp(2,2)*Ftmp(3,1) -
     $             Ftmp(1,2)*Ftmp(2,1)*Ftmp(3,3) -
     $             Ftmp(1,1)*Ftmp(2,3)*Ftmp(3,2)

        Fdef(1,1) = Ftmp(2,2)*Ftmp(3,3)-Ftmp(2,3)*Ftmp(3,2)
        Fdef(1,2) = Ftmp(1,3)*Ftmp(3,2)-Ftmp(1,2)*Ftmp(3,3)
        Fdef(1,3) = Ftmp(1,2)*Ftmp(2,3)-Ftmp(1,3)*Ftmp(2,2)
        Fdef(2,1) = Ftmp(2,3)*Ftmp(3,1)-Ftmp(2,1)*Ftmp(3,3)
        Fdef(2,2) = Ftmp(1,1)*Ftmp(3,3)-Ftmp(1,3)*Ftmp(3,1)
        Fdef(2,3) = Ftmp(1,3)*Ftmp(2,1)-Ftmp(1,1)*Ftmp(2,3)
        Fdef(3,1) = Ftmp(2,1)*Ftmp(3,2)-Ftmp(2,2)*Ftmp(3,1)
        Fdef(3,2) = Ftmp(1,2)*Ftmp(3,1)-Ftmp(1,1)*Ftmp(3,2)
        Fdef(3,3) = Ftmp(1,1)*Ftmp(2,2)-Ftmp(1,2)*Ftmp(2,1)
        Fdef(:,:) = (1.0d0/Ftmp_det)*Fdef(:,:)

c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)*Fdef(3,3) +
     $         Fdef(1,2)*Fdef(2,3)*Fdef(3,1) +
     $         Fdef(1,3)*Fdef(2,1)*Fdef(3,2) -
     $         Fdef(1,3)*Fdef(2,2)*Fdef(3,1) -
     $         Fdef(1,2)*Fdef(2,1)*Fdef(3,3) -
     $         Fdef(1,1)*Fdef(2,3)*Fdef(3,2)

        if (Fdet .lt. 0.0d0) then
          Print*,"elem305.f: Jacobian is negative at element", ielem,
     $           "  Jacob ", Fdet," ... exiting"
c          stop
        negJac = .FALSE.
        endif

c       Finv is inverse of deformation gradient at Gauss Point
        Finv(1,1) = Fdef(2,2)*Fdef(3,3)-Fdef(2,3)*Fdef(3,2)
        Finv(1,2) = Fdef(1,3)*Fdef(3,2)-Fdef(1,2)*Fdef(3,3)
        Finv(1,3) = Fdef(1,2)*Fdef(2,3)-Fdef(1,3)*Fdef(2,2)
        Finv(2,1) = Fdef(2,3)*Fdef(3,1)-Fdef(2,1)*Fdef(3,3)
        Finv(2,2) = Fdef(1,1)*Fdef(3,3)-Fdef(1,3)*Fdef(3,1)
        Finv(2,3) = Fdef(1,3)*Fdef(2,1)-Fdef(1,1)*Fdef(2,3)
        Finv(3,1) = Fdef(2,1)*Fdef(3,2)-Fdef(2,2)*Fdef(3,1)
        Finv(3,2) = Fdef(1,2)*Fdef(3,1)-Fdef(1,1)*Fdef(3,2)
        Finv(3,3) = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)
        Finv(:,:) = (1.0d0/Fdet)*Finv(:,:)

c       compute the spatial derivatives of the pressure and the pressure at the Gauss Point
c        (shap contain the derivative of the shape functions and its value)
        pres(1:4) = 0.0d0
        do i = 1, 4
          do inode = 1, nnode
            pres(i) = pres(i) + uelem(4,inode)*shap(i,inode)
          end do
        end do

c       compute the Cauchy tensor and its inverse
        Ctens(1:3,1:3) = 0.0d0
        do j = 1,ndim
          do i = j,ndim! Ctens is symmetric
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do
c       enforce the symmetries
        Ctens(1,2)=Ctens(2,1)
        Ctens(1,3)=Ctens(3,1)
        Ctens(2,3)=Ctens(3,2)

        Cdet = Ctens(1,1)*Ctens(2,2)*Ctens(3,3) +
     $         Ctens(1,2)*Ctens(2,3)*Ctens(3,1) +
     $         Ctens(1,3)*Ctens(2,1)*Ctens(3,2) -
     $         Ctens(1,3)*Ctens(2,2)*Ctens(3,1) -
     $         Ctens(1,2)*Ctens(2,1)*Ctens(3,3) -
     $         Ctens(1,1)*Ctens(2,3)*Ctens(3,2)
c JFD: note Cdet=Fdet*Fdet...
        Cinv(1,1) = Ctens(2,2)*Ctens(3,3)-Ctens(2,3)*Ctens(3,2)
        Cinv(1,2) = Ctens(1,3)*Ctens(3,2)-Ctens(1,2)*Ctens(3,3)
        Cinv(1,3) = Ctens(1,2)*Ctens(2,3)-Ctens(1,3)*Ctens(2,2)
        Cinv(2,1) = Cinv(1,2)! = Ctens(2,3)*Ctens(3,1)-Ctens(2,1)*Ctens(3,3)
        Cinv(2,2) = Ctens(1,1)*Ctens(3,3)-Ctens(1,3)*Ctens(3,1)
        Cinv(2,3) = Ctens(1,3)*Ctens(2,1)-Ctens(1,1)*Ctens(2,3)
        Cinv(3,1) = Cinv(1,3)! = Ctens(2,1)*Ctens(3,2)-Ctens(2,2)*Ctens(3,1)
        Cinv(3,2) = Cinv(2,3)! = Ctens(1,2)*Ctens(3,1)-Ctens(1,1)*Ctens(3,2)
        Cinv(3,3) = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)
        Cinv(:,:) = (1.0d0/Cdet)*Cinv(:,:)


c       principal invariant Inv1
        Inv1 = Ctens(1,1)+Ctens(2,2)+Ctens(3,3)
         
c       dJdC is the derivative of the Jacobian with respect to 
c       the Cauchy Green tensor
        dJdC(1:3,1:3) = 0.5d0*Fdet*Cinv(1:3,1:3)

c       compute the second Piola Kirchhoff stress SecPK 
c       and the material tangent Ctang
        tmp5 = Fdet**(-2.0d0/3.0d0)! tmp5 is re-used here to avoid computing the power multiple times
        K1 = tmp5*Inv1-3.0d0
        tmp2 = 1.0d0/3.0d0 ! tmp2 is reused to avoid computing the fraction multiple times
        K2(1:3,1:3) = (ident(1:3,1:3) - 
     $                tmp2*Cinv(1:3,1:3)*Inv1)*tmp5
        tmp5=tmp5/3.0d0
        tmp4=0.5d0*mu*exp(K1*gamm)! tmp4 is re-used here to avoid computing the exp multiple times
        tmp3=0.25d0*Fdet! tmp3 is reused to avoid computing the product multiple times

        dWdC(1:3,1:3) = tmp4*K2(1:3,1:3)

        SecPK(1:3,1:3) = 2.0d0*(dWdC(1:3,1:3)-pres(4)*dJdC(1:3,1:3))

c       compute Cauchy stress tensor
        Sigma(1:3,1:3) = 0.0d0
        do i = 1, ndim
         do r = 1, ndim
          do q = 1, ndim
           do j = 1, ndim
            Sigma(i,r) = Sigma(i,r) +
     $                   Fdef(i,j)*SecPK(j,q)*Fdef(r,q)/Fdet
           end do
          end do
         end do
        end do

c        if (ielem==90.AND.iinte==1) then
c          Print*,Fdef(1,:)
c          Print*,Fdef(2,:)
c          Print*,Fdef(3,:)
c          Print*
c          Print*,SecPK(1,:)
c          Print*,SecPK(2,:)
c          Print*,SecPK(3,:)
c        endif

c        if (.false.) then
c        do i = 1, ndime
c          do j = 1, ndime
c            do q = 1, ndime
c              do r = 1, ndime
c                dCinvdC(i,j,q,r) = -0.5d0*(Cinv(i,q)*Cinv(j,r)+
c     $                                   Cinv(i,r)*Cinv(j,q)) 
c                Ctang(i,j,q,r) = tmp4*((gamm*K2(i,j)*K2(q,r))-
c     $                          (1.0d0/3.0d0)*K2(i,j)*Cinv(q,r)-
c     $                          tmp5*(dCinvdC(i,j,q,r)*Inv1
c     $                          +Cinv(i,j)*ident(q,r)))
c                d2JdC(i,j,q,r) = 0.25d0*Fdet*(Cinv(q,r)*Cinv(i,j)
c     $                           +2.0d0*dCinvdC(i,j,q,r))
c              end do
c            end do
c          end do
c        end do
c        else
cc        Ctangt(:,:,:,:)=Ctang(:,:,:,:)
cc        d2JdCt(:,:,:,:)=d2JdC(:,:,:,:)

c        Ctang and d2JdC have minor and major symmetries
c        that is each one has 21 independent components
c       i=1,j=1,q=1,r=1
        dCinvdCt=-Cinv(1,1)*Cinv(1,1)
        Ctang(1,1,1,1) = tmp4*(gamm*K2(1,1)*K2(1,1) - tmp2*
     $            K2(1,1)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(1,1)))
        d2JdC(1,1,1,1)=tmp3*dCinvdCt
c       i=2,j=1,q=1,r=1
        dCinvdCt=-Cinv(2,1)*Cinv(1,1)
        Ctang(2,1,1,1) = tmp4*(gamm*K2(2,1)*K2(1,1) - tmp2*
     $            K2(2,1)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(2,1)))
        d2JdC(2,1,1,1)=tmp3*dCinvdCt
        Ctang(1,2,1,1)=Ctang(2,1,1,1)
        Ctang(1,1,2,1)=Ctang(2,1,1,1)
        d2JdC(1,2,1,1)=d2JdC(2,1,1,1)
        d2JdC(1,1,2,1)=d2JdC(2,1,1,1)
c       i=3,j=1,q=1,r=1
        dCinvdCt=-Cinv(3,1)*Cinv(1,1)
        Ctang(3,1,1,1) = tmp4*(gamm*K2(3,1)*K2(1,1) - tmp2*
     $            K2(3,1)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(3,1)))
        d2JdC(3,1,1,1)=tmp3*dCinvdCt
        Ctang(1,3,1,1)=Ctang(3,1,1,1)
        Ctang(1,1,3,1)=Ctang(3,1,1,1)
        d2JdC(1,3,1,1)=d2JdC(3,1,1,1)
        d2JdC(1,1,3,1)=d2JdC(3,1,1,1)
c       i=2,j=2,q=1,r=1
        dCinvdCt=-Cinv(2,1)*Cinv(2,1)
        Ctang(2,2,1,1) = tmp4*(gamm*K2(2,2)*K2(1,1) - tmp2*
     $            K2(2,2)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(2,2)))
        d2JdC(2,2,1,1)=tmp3*(Cinv(1,1)*Cinv(2,2)+2.0d0*dCinvdCt)
        Ctang(1,1,2,2)=Ctang(2,2,1,1)
        d2JdC(1,1,2,2)=d2JdC(2,2,1,1)
c       i=3,j=2,q=1,r=1
        dCinvdCt=-Cinv(3,1)*Cinv(2,1)
        Ctang(3,2,1,1) = tmp4*(gamm*K2(3,2)*K2(1,1) - tmp2*
     $            K2(3,2)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(3,2)))
        d2JdC(3,2,1,1)=tmp3*(Cinv(1,1)*Cinv(3,2)+2.0d0*dCinvdCt)
        Ctang(2,3,1,1)=Ctang(3,2,1,1)
        Ctang(1,1,3,2)=Ctang(3,2,1,1)
        d2JdC(2,3,1,1)=d2JdC(3,2,1,1)
        d2JdC(1,1,3,2)=d2JdC(3,2,1,1)
c       i=3,j=3,q=1,r=1
        dCinvdCt=-Cinv(3,1)*Cinv(3,1)
        Ctang(3,3,1,1) = tmp4*(gamm*K2(3,3)*K2(1,1) - tmp2*
     $            K2(3,3)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(3,3)))
        d2JdC(3,3,1,1)=tmp3*(Cinv(1,1)*Cinv(3,3)+2.0d0*dCinvdCt)
        Ctang(1,1,3,3)=Ctang(3,3,1,1)
        d2JdC(1,1,3,3)=d2JdC(3,3,1,1)
c       i=2,j=1,q=2,r=1
        dCinvdCt=-0.5d0*(Cinv(2,2)*Cinv(1,1)+Cinv(2,1)*Cinv(1,2))
        Ctang(2,1,2,1) = tmp4*(gamm*K2(2,1)*K2(2,1) - tmp2*
     $            K2(2,1)*Cinv(2,1) - tmp5*dCinvdCt*Inv1)
        d2JdC(2,1,2,1)=-tmp3*Cinv(2,2)*Cinv(1,1)
        Ctang(1,2,2,1)=Ctang(2,1,2,1)
        d2JdC(1,2,2,1)=d2JdC(2,1,2,1)
c       i=3,j=1,q=2,r=1
        dCinvdCt=-0.5d0*(Cinv(3,2)*Cinv(1,1)+Cinv(3,1)*Cinv(1,2))
        Ctang(3,1,2,1) = tmp4*(gamm*K2(3,1)*K2(2,1) - tmp2*
     $            K2(3,1)*Cinv(2,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,1,2,1)=-tmp3*Cinv(1,1)*Cinv(3,2)
        Ctang(1,3,2,1)=Ctang(3,1,2,1)
        Ctang(2,1,3,1)=Ctang(3,1,2,1)
        Ctang(1,2,3,1)=Ctang(3,1,2,1)
        d2JdC(1,3,2,1)=d2JdC(3,1,2,1)
        d2JdC(2,1,3,1)=d2JdC(3,1,2,1)
        d2JdC(1,2,3,1)=d2JdC(3,1,2,1)
c       i=2,j=2,q=2,r=1
        dCinvdCt=-Cinv(2,2)*Cinv(2,1)
        Ctang(2,2,2,1) = tmp4*(gamm*K2(2,2)*K2(2,1) - tmp2*
     $            K2(2,2)*Cinv(2,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(2,2,2,1)=tmp3*dCinvdCt
        Ctang(2,1,2,2)=Ctang(2,2,2,1)
        Ctang(1,2,2,2)=Ctang(2,2,2,1)
        d2JdC(2,1,2,2)=d2JdC(2,2,2,1)
        d2JdC(1,2,2,2)=d2JdC(2,2,2,1)
c       i=3,j=2,q=2,r=1
        dCinvdCt=-0.5d0*(Cinv(3,2)*Cinv(2,1)+Cinv(3,1)*Cinv(2,2))
        Ctang(3,2,2,1) = tmp4*(gamm*K2(3,2)*K2(2,1) - tmp2*
     $            K2(3,2)*Cinv(2,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,2,2,1)=-tmp3*Cinv(3,1)*Cinv(2,2)
        Ctang(2,3,2,1)=Ctang(3,2,2,1)
        Ctang(2,1,3,2)=Ctang(3,2,2,1)
        Ctang(1,2,3,2)=Ctang(3,2,2,1)
        d2JdC(2,3,2,1)=d2JdC(3,2,2,1)
        d2JdC(2,1,3,2)=d2JdC(3,2,2,1)
        d2JdC(1,2,3,2)=d2JdC(3,2,2,1)
c       i=3,j=3,q=2,r=1
        dCinvdCt=-Cinv(3,2)*Cinv(3,1)
        Ctang(3,3,2,1) = tmp4*(gamm*K2(3,3)*K2(2,1) - tmp2*
     $            K2(3,3)*Cinv(2,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,3,2,1)=tmp3*(Cinv(2,1)*Cinv(3,3)+2.0d0*dCinvdCt)
        Ctang(2,1,3,3)=Ctang(3,3,2,1)
        Ctang(1,2,3,3)=Ctang(3,3,2,1)
        d2JdC(2,1,3,3)=d2JdC(3,3,2,1)
        d2JdC(1,2,3,3)=d2JdC(3,3,2,1)
c       i=3,j=1,q=3,r=1
        dCinvdCt=-0.5d0*(Cinv(3,3)*Cinv(1,1)+Cinv(3,1)*Cinv(1,3))
        Ctang(3,1,3,1) = tmp4*(gamm*K2(3,1)*K2(3,1) - tmp2*
     $            K2(3,1)*Cinv(3,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,1,3,1)=-tmp3*Cinv(3,3)*Cinv(1,1)
        Ctang(1,3,3,1)=Ctang(3,1,3,1)
        d2JdC(1,3,3,1)=d2JdC(3,1,3,1)
c       i=2,j=2,q=3,r=1
        dCinvdCt=-Cinv(2,1)*Cinv(2,3)
        Ctang(2,2,3,1) = tmp4*(gamm*K2(2,2)*K2(3,1) - tmp2*
     $            K2(2,2)*Cinv(3,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(2,2,3,1)=tmp3*(Cinv(3,1)*Cinv(2,2)+2.0d0*dCinvdCt)
        Ctang(3,1,2,2)=Ctang(2,2,3,1)
        Ctang(1,3,2,2)=Ctang(2,2,3,1)
        d2JdC(3,1,2,2)=d2JdC(2,2,3,1)
        d2JdC(1,3,2,2)=d2JdC(2,2,3,1)
c       i=3,j=2,q=3,r=1
        dCinvdCt=-0.5d0*(Cinv(3,3)*Cinv(2,1)+Cinv(3,1)*Cinv(2,3))
        Ctang(3,2,3,1) = tmp4*(gamm*K2(3,2)*K2(3,1) - tmp2*
     $            K2(3,2)*Cinv(3,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,2,3,1)=-tmp3*Cinv(3,3)*Cinv(2,1)
        Ctang(2,3,3,1)=Ctang(3,2,3,1)
        Ctang(3,1,3,2)=Ctang(3,2,3,1)
        Ctang(1,3,3,2)=Ctang(3,2,3,1)
        d2JdC(2,3,3,1)=d2JdC(3,2,3,1)
        d2JdC(3,1,3,2)=d2JdC(3,2,3,1)
        d2JdC(1,3,3,2)=d2JdC(3,2,3,1)
c       i=3,j=3,q=3,r=1
        dCinvdCt=-Cinv(3,3)*Cinv(3,1)
        Ctang(3,3,3,1) = tmp4*(gamm*K2(3,3)*K2(3,1) - tmp2*
     $            K2(3,3)*Cinv(3,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,3,3,1)=tmp3*dCinvdCt
        Ctang(3,1,3,3)=Ctang(3,3,3,1)
        Ctang(1,3,3,3)=Ctang(3,3,3,1)
        d2JdC(3,1,3,3)=d2JdC(3,3,3,1)
        d2JdC(1,3,3,3)=d2JdC(3,3,3,1)
c       i=2,j=2,q=2,r=2
        dCinvdCt=-Cinv(2,2)*Cinv(2,2)
        Ctang(2,2,2,2) = tmp4*(gamm*K2(2,2)*K2(2,2) - tmp2*
     $            K2(2,2)*Cinv(2,2)-tmp5*(dCinvdCt*Inv1+Cinv(2,2)))
        d2JdC(2,2,2,2)=tmp3*dCinvdCt
c       i=3,j=2,q=2,r=2
        dCinvdCt=-Cinv(3,2)*Cinv(2,2)
        Ctang(3,2,2,2) = tmp4*(gamm*K2(3,2)*K2(2,2) - tmp2*
     $            K2(3,2)*Cinv(2,2)-tmp5*(dCinvdCt*Inv1+Cinv(3,2)))
        d2JdC(3,2,2,2)=tmp3*dCinvdCt
        Ctang(2,3,2,2)=Ctang(3,2,2,2)
        Ctang(2,2,3,2)=Ctang(3,2,2,2)
        d2JdC(2,3,2,2)=d2JdC(3,2,2,2)
        d2JdC(2,2,3,2)=d2JdC(3,2,2,2)
c       i=3,j=3,q=2,r=2
        dCinvdCt=-Cinv(3,2)*Cinv(3,2)
        Ctang(3,3,2,2) = tmp4*(gamm*K2(3,3)*K2(2,2) - tmp2*
     $            K2(3,3)*Cinv(2,2)-tmp5*(dCinvdCt*Inv1+Cinv(3,3)))
        d2JdC(3,3,2,2)=tmp3*(Cinv(2,2)*Cinv(3,3)+2.0d0*dCinvdCt)
        Ctang(2,2,3,3)=Ctang(3,3,2,2)
        d2JdC(2,2,3,3)=d2JdC(3,3,2,2)
c       i=3,j=2,q=3,r=2
        dCinvdCt=-0.5d0*(Cinv(3,3)*Cinv(2,2)+Cinv(3,2)*Cinv(3,2))
        Ctang(3,2,3,2) = tmp4*(gamm*K2(3,2)*K2(3,2) - tmp2*
     $            K2(3,2)*Cinv(3,2)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,2,3,2)=-tmp3*Cinv(3,3)*Cinv(2,2)
        Ctang(2,3,3,2)=Ctang(3,2,3,2)
        d2JdC(2,3,3,2)=d2JdC(3,2,3,2)
c       i=3,j=3,q=3,r=2
        dCinvdCt=-Cinv(3,3)*Cinv(3,2)
        Ctang(3,3,3,2) = tmp4*(gamm*K2(3,3)*K2(3,2) - tmp2*
     $            K2(3,3)*Cinv(3,2)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,3,3,2)=tmp3*dCinvdCt
        Ctang(3,2,3,3)=Ctang(3,3,3,2)
        Ctang(2,3,3,3)=Ctang(3,3,3,2)
        d2JdC(3,2,3,3)=d2JdC(3,3,3,2)
        d2JdC(2,3,3,3)=d2JdC(3,3,3,2)
c       i=3,j=3,q=3,r=3
        dCinvdCt=-Cinv(3,3)*Cinv(3,3)
        Ctang(3,3,3,3) = tmp4*(gamm*K2(3,3)*K2(3,3) - tmp2*
     $            K2(3,3)*Cinv(3,3)-tmp5*(dCinvdCt*Inv1+Cinv(3,3)))
        d2JdC(3,3,3,3)=tmp3*dCinvdCt
c       enforce the symmetries
        Ctang(:,:,1,2)=Ctang(:,:,2,1)
        Ctang(:,:,1,3)=Ctang(:,:,3,1)
        Ctang(:,:,2,3)=Ctang(:,:,3,2)
        d2JdC(:,:,1,2)=d2JdC(:,:,2,1)
        d2JdC(:,:,1,3)=d2JdC(:,:,3,1)
        d2JdC(:,:,2,3)=d2JdC(:,:,3,2)

c        d2JdCt(:,:,:,:)=d2JdCt(:,:,:,:)-d2JdC(:,:,:,:)
c        Ctangt(:,:,:,:)=Ctangt(:,:,:,:)-Ctang(:,:,:,:)
c        if (ielem.eq.2) then
c        if (iinte.eq.2) then
c        Print*,"1,1:",Ctangt(:,:,1,1)
c        Print*,"1,1:",Ctang(:,:,1,1)
c        Print*,"2,1:",Ctangt(:,:,2,1)
c        Print*,"2,1:",Ctang(:,:,2,1)
c        Print*,"3,1:",Ctangt(:,:,3,1)
c        Print*,"3,1:",Ctang(:,:,3,1)
c        Print*,"1,2:",Ctangt(:,:,1,2)
c        Print*,"1,2:",Ctang(:,:,1,2)
c        Print*,"2,2:",Ctangt(:,:,2,2)
c        Print*,"2,2:",Ctang(:,:,2,2)
c        Print*,"3,2:",Ctangt(:,:,3,2)
c        Print*,"3,2:",Ctang(:,:,3,2)
c        Print*,"1,3:",Ctangt(:,:,1,3)
c        Print*,"1,3:",Ctang(:,:,1,3)
c        Print*,"2,3:",Ctangt(:,:,2,3)
c        Print*,"2,3:",Ctang(:,:,2,3)
c        Print*,"3,3:",Ctangt(:,:,3,3)
c        Print*,"3,3:",Ctang(:,:,3,3)
c        endif
c        endif
c        if (.false.) then!if (ielem.eq.2) then
c        if (iinte.eq.2) then
c        Print*,"1,1:",d2JdCt(:,:,1,1)
c        Print*,"1,1:",d2JdC(:,:,1,1)
c        Print*,"2,1:",d2JdCt(:,:,2,1)
c        Print*,"2,1:",d2JdC(:,:,2,1)
c        Print*,"3,1:",d2JdCt(:,:,3,1)
c        Print*,"3,1:",d2JdC(:,:,3,1)
c        Print*,"1,2:",d2JdCt(:,:,1,2)
c        Print*,"1,2:",d2JdC(:,:,1,2)
c        Print*,"2,2:",d2JdCt(:,:,2,2)
c        Print*,"2,2:",d2JdC(:,:,2,2)
c        Print*,"3,2:",d2JdCt(:,:,3,2)
c        Print*,"3,2:",d2JdC(:,:,3,2)
c        Print*,"1,3:",d2JdCt(:,:,1,3)
c        Print*,"1,3:",d2JdC(:,:,1,3)
c        Print*,"2,3:",d2JdCt(:,:,2,3)
c        Print*,"2,3:",d2JdC(:,:,2,3)
c        Print*,"3,3:",d2JdCt(:,:,3,3)
c        Print*,"3,3:",d2JdC(:,:,3,3)
c        endif
c        endif
c        endif

c       construct the element stiffness matrix 
        Igeo(1,1) = Sigma(1,1)
        Igeo(1,2) = 0.0d0
        Igeo(2,2) = Sigma(2,2)
        Igeo(1,3) = 0.0d0
        Igeo(2,3) = 0.0d0
        Igeo(3,3) = Sigma(3,3)
        Igeo(1,4) = 0.5d0*Sigma(1,2)
        Igeo(2,4) = Igeo(1,4)! = 0.5d0*Sigma(2,1)
        Igeo(3,4) = 0.0d0
        Igeo(4,4) = 0.25d0*(Sigma(2,2)+Sigma(1,1))
        Igeo(1,5) = 0.5d0*Sigma(1,3)
        Igeo(2,5) = 0.0d0
        Igeo(3,5) = Igeo(1,5)! = 0.5d0*Sigma(3,1)
        Igeo(4,5) = 0.25d0*Sigma(2,3)
        Igeo(5,5) = 0.25d0*(Sigma(3,3)+Sigma(1,1))
        Igeo(1,6) = 0.0d0
        Igeo(2,6) = 0.5d0*Sigma(2,3)
        Igeo(3,6) = Igeo(2,6)! = 0.5d0*Sigma(3,2)
        Igeo(4,6) = 0.25d0*Sigma(1,3)
        Igeo(5,6) = 0.25d0*Sigma(1,2)
        Igeo(6,6) = 0.25d0*(Sigma(3,3)+Sigma(2,2))
        Igeo(1,7) = 0.5d0*Sigma(1,2)
        Igeo(2,7) = -0.5d0*Sigma(2,1)
        Igeo(3,7) = 0.0d0
        Igeo(4,7) = 0.25d0*(Sigma(2,2)-Sigma(1,1))
        Igeo(5,7) = 0.25d0*Sigma(3,2)
        Igeo(6,7) = -0.25d0*Sigma(3,1)
        Igeo(7,7) = 0.25d0*(Sigma(2,2)+Sigma(1,1))
        Igeo(1,8) = Igeo(1,5)! = 0.5d0*Sigma(1,3)
        Igeo(2,8) = 0.0d0
        Igeo(3,8) = -0.5d0*Sigma(3,1)
        Igeo(4,8) = Igeo(5,7)! = 0.25d0*Sigma(2,3)
        Igeo(5,8) = 0.25d0*(Sigma(3,3)-Sigma(1,1))
        Igeo(6,8) = -0.25d0*Sigma(2,1)
        Igeo(7,8) = Igeo(4,8)! = 0.25d0*Sigma(2,3)
        Igeo(8,8) = 0.25d0*(Sigma(3,3)+Sigma(1,1))
        Igeo(1,9) = 0.0d0
        Igeo(2,9) = Igeo(2,6)! = 0.5d0*Sigma(2,3)
        Igeo(3,9) = -0.5d0*Sigma(3,2)
        Igeo(4,9) = 0.25d0*Sigma(1,3)
        Igeo(5,9) = -0.25d0*Sigma(1,2)
        Igeo(6,9) = 0.25d0*(Sigma(3,3)-Sigma(2,2))
        Igeo(7,9) = -0.25d0*Sigma(1,3)
        Igeo(8,9) = Igeo(5,6)! = 0.25d0*Sigma(1,2)
        Igeo(9,9) = Igeo(6,6)! = 0.25d0*(Sigma(2,2)+Sigma(3,3))
c       enforce the symmetries
        do i = 1,8
          do j = (i+1),9
            Igeo(j,i) = Igeo(i,j)
          end do
        end do
c        if (ielem.eq.2) then
c        if (iinte.eq.2) then
c        Print*,"Igeo:",Igeo
c        endif
c        endif
 
        ! build Ctang=Ctang-pres(4)*d2JdC now (note d2JdC is reused after)
        Ctang(:,:,:,:)=Ctang(:,:,:,:)-pres(4)*d2JdC(:,:,:,:)
c       Lmat shows major symmetries [(i,j)<->(q,r)] JFD how can this be implemented fast?
        Lmat(:,:,:,:) = 0.0d0
        do r = 1, ndim
          do q = 1, ndim
            tmp6=q+(r-1)*ndim
            do j = 1, r! major symmetry
              do i = 1, ndim
                tmp7=i+(j-1)*ndim
                if (tmp7.le.tmp6) then! use more major symmetries
                  do k = 1, ndim
                    do t = 1, ndim
                     do kk = 1, ndim
                      do tt = 1, ndim
                        Lmat(i,j,q,r) = Lmat(i,j,q,r) +
     $                  4.0d0*Fdef(i,k)*Fdef(q,t)*Fdef(j,kk)*Fdef(r,tt)*
     $                  Ctang(kk,k,tt,t)/Fdet
                      end do
                     end do
                    end do
                  end do
                  Lmat(q,r,i,j)=Lmat(i,j,q,r)! major symmetry
                endif
              end do
            end do
          end do
        end do

        Dgeomat(1,1) = Lmat(1,1,1,1)
        Dgeomat(1,2) = Lmat(1,1,2,2)
        Dgeomat(2,2) = Lmat(2,2,2,2)
        Dgeomat(1,3) = Lmat(1,1,3,3)
        Dgeomat(2,3) = Lmat(2,2,3,3) 
        Dgeomat(3,3) = Lmat(3,3,3,3)
        Dgeomat(1,4) = 0.5d0*(Lmat(1,1,1,2)+Lmat(1,1,2,1))
        Dgeomat(2,4) = 0.5d0*(Lmat(2,2,1,2)+Lmat(2,2,2,1))
        Dgeomat(3,4) = 0.5d0*(Lmat(3,3,1,2)+Lmat(3,3,2,1))
        Dgeomat(4,4) = 0.25d0*(Lmat(1,2,1,2)+Lmat(2,1,2,1))
     $                   +0.5d0*Lmat(1,2,2,1)
        Dgeomat(1,5) = 0.5d0*(Lmat(1,1,1,3)+Lmat(1,1,3,1))
        Dgeomat(2,5) = 0.5d0*(Lmat(2,2,1,3)+Lmat(2,2,3,1))
        Dgeomat(3,5) = 0.5d0*(Lmat(3,3,1,3)+Lmat(3,3,3,1))
        Dgeomat(4,5) = 0.25d0*(Lmat(1,2,1,3)+Lmat(1,2,3,1))
     $                    + 0.25d0*(Lmat(2,1,1,3)+Lmat(2,1,3,1))
        Dgeomat(5,5) = 0.25d0*(Lmat(1,3,1,3)+Lmat(3,1,3,1))
     $                     +0.5d0*Lmat(1,3,3,1)
        Dgeomat(1,6) = 0.5d0*(Lmat(1,1,2,3)+Lmat(1,1,3,2))
        Dgeomat(2,6) = 0.5d0*(Lmat(2,2,2,3)+Lmat(2,2,3,2))
        Dgeomat(3,6) = 0.5d0*(Lmat(3,3,2,3)+Lmat(3,3,3,2))
        Dgeomat(4,6) = 0.25d0*(Lmat(1,2,2,3)+Lmat(1,2,3,2))
     $                     + 0.25d0*(Lmat(2,1,2,3)+Lmat(2,1,3,2))
        Dgeomat(5,6) = 0.25d0*(Lmat(1,3,2,3)+Lmat(1,3,3,2))
     $                     + 0.25d0*(Lmat(3,1,2,3)+Lmat(3,1,3,2))
        Dgeomat(6,6) = 0.25d0*(Lmat(2,3,2,3)+Lmat(3,2,3,2))
     $                    +0.5d0*Lmat(2,3,3,2)
        Dgeomat(1,7) = 0.5d0*(Lmat(1,1,1,2)-Lmat(1,1,2,1))
        Dgeomat(2,7) = 0.5d0*(Lmat(2,2,1,2)-Lmat(2,2,2,1))
        Dgeomat(3,7) = 0.5d0*(Lmat(3,3,1,2)-Lmat(3,3,2,1))
        Dgeomat(4,7) = 0.25d0*(Lmat(1,2,1,2)-Lmat(2,1,2,1))
        Dgeomat(5,7) = 0.25d0*(Lmat(1,3,1,2)-Lmat(1,3,2,1))
     $                     + 0.25d0*(Lmat(3,1,1,2)-Lmat(3,1,2,1))
        Dgeomat(6,7) = 0.25d0*(Lmat(2,3,1,2)-Lmat(2,3,2,1))
     $                     + 0.25d0*(Lmat(3,2,1,2)-Lmat(3,2,2,1))
        Dgeomat(7,7) = 0.25d0*(Lmat(1,2,1,2)+Lmat(2,1,2,1))
     $                   -0.5d0*Lmat(1,2,2,1)
        Dgeomat(1,8) = 0.5d0*(Lmat(1,1,1,3)-Lmat(1,1,3,1))
        Dgeomat(2,8) = 0.5d0*(Lmat(2,2,1,3)-Lmat(2,2,3,1))
        Dgeomat(3,8) = 0.5d0*(Lmat(3,3,1,3)-Lmat(3,3,3,1))
        Dgeomat(4,8) = 0.25d0*(Lmat(1,2,1,3)-Lmat(1,2,3,1))
     $                     + 0.25d0*(Lmat(2,1,1,3)-Lmat(2,1,3,1))
        Dgeomat(5,8) = 0.25d0*(Lmat(1,3,1,3)-Lmat(3,1,3,1))
        Dgeomat(6,8) = 0.25d0*(Lmat(2,3,1,3)-Lmat(2,3,3,1))
     $                     + 0.25d0*(Lmat(3,2,1,3)-Lmat(3,2,3,1))
        Dgeomat(7,8) = 0.25d0*(Lmat(1,2,1,3)-Lmat(1,2,3,1))
     $                     - 0.25d0*(Lmat(2,1,1,3)-Lmat(2,1,3,1))
        Dgeomat(8,8) = 0.25d0*(Lmat(1,3,1,3)+Lmat(3,1,3,1))
     $                   - 0.5d0*Lmat(1,3,3,1)
        Dgeomat(1,9) = 0.5d0*(Lmat(1,1,2,3)-Lmat(1,1,3,2))
        Dgeomat(2,9) = 0.5d0*(Lmat(2,2,2,3)-Lmat(2,2,3,2))
        Dgeomat(3,9) = 0.5d0*(Lmat(3,3,2,3)-Lmat(3,3,3,2))
        Dgeomat(4,9) = 0.25d0*(Lmat(1,2,2,3)-Lmat(1,2,3,2))
     $                     + 0.25d0*(Lmat(2,1,2,3)-Lmat(2,1,3,2))
        Dgeomat(5,9) = 0.25d0*(Lmat(1,3,2,3)-Lmat(1,3,3,2))
     $                     + 0.25d0*(Lmat(3,1,2,3)-Lmat(3,1,3,2))
        Dgeomat(6,9) = 0.25d0*(Lmat(2,3,2,3)-Lmat(3,2,3,2))
        Dgeomat(7,9) = 0.25d0*(Lmat(1,2,2,3)-Lmat(1,2,3,2))
     $                     - 0.25d0*(Lmat(2,1,2,3)-Lmat(2,1,3,2))
        Dgeomat(8,9) = 0.25d0*(Lmat(1,3,2,3)-Lmat(1,3,3,2))
     $                     - 0.25d0*(Lmat(3,1,2,3)-Lmat(3,1,3,2))
        Dgeomat(9,9) = 0.25d0*(Lmat(2,3,2,3)+Lmat(3,2,3,2))
     $                   - 0.5d0*Lmat(2,3,3,2)
c       enforce the symmetries
        do i = 1,8
          do j = (i+1),9
            Dgeomat(j,i) = Dgeomat(i,j)
          end do 
        end do

        Dgeomat(1:9,1:9) = Dgeomat(1:9,1:9)+Igeo(1:9,1:9)

c       create the b-matrix
        bb(:,:,:) = 0.0d0
        do inode = 1, nnode
          bb(inode,1,1) = shap(1,inode)
          bb(inode,2,2) = shap(2,inode)
          bb(inode,3,3) = shap(3,inode)
          bb(inode,4,1) = shap(2,inode)
          bb(inode,4,2) = shap(1,inode)
          bb(inode,5,1) = shap(3,inode)
          bb(inode,5,3) = shap(1,inode)
          bb(inode,6,2) = shap(3,inode)
          bb(inode,6,3) = shap(2,inode)
          bb(inode,7,1) = shap(2,inode)
          bb(inode,7,2) = -shap(1,inode)
          bb(inode,8,1) = shap(3,inode)
          bb(inode,8,3) = -shap(1,inode)
          bb(inode,9,2) = shap(3,inode)
          bb(inode,9,3) = -shap(2,inode)
        end do  

        ievab = 0
        do inode = 1, nnode
          do idofn = 1,4! 4 = elemvec_ndofn(inode)
            ievab = ievab+1
            jevab = 0
            do jnode = 1, nnode
              do jdofn = 1,4! 4 = elemvec_ndofn(jnode)
                jevab = jevab+1
                if (idofn.le.3  .AND. jdofn.eq.4) then
                  estif(ievab,jevab)=estif(ievab,jevab)-
     $                 shap(idofn,inode)*shap(4,jnode)*wtjac
                elseif (idofn.eq.4 .AND. jdofn.le.3) then
                  estif(ievab,jevab)=estif(ievab,jevab)+Fdet*
     $               (shap(1,jnode)*Finv(1,jdofn)+
     $                shap(2,jnode)*Finv(2,jdofn)+
     $                shap(3,jnode)*Finv(3,jdofn))
     $                                  *shap(4,inode)*wtjac! symmetric part
                elseif (idofn.eq.4  .AND. jdofn.eq.4) then
                  do i = 1,ndim
                    estif(ievab,jevab)=estif(ievab,jevab)+
     $                 temp*shap(i,jnode)*shap(i,inode)
     $                 *wtjac
                  end do
                else! (idofn.lt.4  .AND. jdofn.lt.4) then
                 if (ievab.le.jevab) then! use symmetries
                  do i = 1, 9
                    do j = 1, 9
                      estif(ievab,jevab)=estif(ievab,jevab)+
     $                         bb(inode,i,idofn)*Dgeomat(i,j)*
     $                         bb(jnode,j,jdofn)*wtjac
                    end do
                  end do
                 endif! the symmetric part is filled after the end of the Gauss integration
                end if
              end do! jdofn
            end do! jnode
          end do! idofn
        end do! inode

c       create element residual for right hand side
        Tp(1:4) = 0.0d0
        do i = 1, ndim
          Tp(i) = Tp(i) + temp*pres(i)
        end do

        Tp(4) = Fdet-1.0d0 

        FS(:,:) = 0.0d0
        do i = 1, ndim
          do j = 1, ndim
            do r = 1, ndim
               FS(i,j) = FS(i,j)+Fdef(i,r)*SecPK(r,j)
            end do
          end do
        end do
c       JFD: extend FS to account for Tp and avoid the if loop below?
        ievab = 0
        do inode = 1, nnode
          do i = 1,4! 4 = elemvec_ndofn(inode)
            ievab = ievab+1
            if (i.eq.4) then!4 = elemvec_ndofn(inode)
              do j = 1,4
                eforc(ievab)=eforc(ievab)+
     $                        shap(j,inode)*Tp(j)*wtjac 
              end do
            else
              do j = 1, 3
                eforc(ievab)=eforc(ievab)+
     $                        shap(j,inode)*Sigma(i,j)*wtjac
              end do
            end if
          end do
        end do! inode

      end do! iinte
c!$acc end region
c!$OMP END PARALLEL DO

      estif(2,1)=estif(1,2)
      estif(3,1)=estif(1,3)
      estif(5,1)=estif(1,5)
      estif(6,1)=estif(1,6)
      estif(7,1)=estif(1,7)
      estif(9,1)=estif(1,9)
      estif(10,1)=estif(1,10)
      estif(11,1)=estif(1,11)
      estif(13,1)=estif(1,13)
      estif(14,1)=estif(1,14)
      estif(15,1)=estif(1,15)
      estif(3,2)=estif(2,3)
      estif(5,2)=estif(2,5)
      estif(6,2)=estif(2,6)
      estif(7,2)=estif(2,7)
      estif(9,2)=estif(2,9)
      estif(10,2)=estif(2,10)
      estif(11,2)=estif(2,11)
      estif(13,2)=estif(2,13)
      estif(14,2)=estif(2,14)
      estif(15,2)=estif(2,15)
      estif(5,3)=estif(3,5)
      estif(6,3)=estif(3,6)
      estif(7,3)=estif(3,7)
      estif(9,3)=estif(3,9)
      estif(10,3)=estif(3,10)
      estif(11,3)=estif(3,11)
      estif(13,3)=estif(3,13)
      estif(14,3)=estif(3,14)
      estif(15,3)=estif(3,15)
      estif(6,5)=estif(5,6)
      estif(7,5)=estif(5,7)
      estif(9,5)=estif(5,9)
      estif(10,5)=estif(5,10)
      estif(11,5)=estif(5,11)
      estif(13,5)=estif(5,13)
      estif(14,5)=estif(5,14)
      estif(15,5)=estif(5,15)
      estif(7,6)=estif(6,7)
      estif(9,6)=estif(6,9)
      estif(10,6)=estif(6,10)
      estif(11,6)=estif(6,11)
      estif(13,6)=estif(6,13)
      estif(14,6)=estif(6,14)
      estif(15,6)=estif(6,15)
      estif(9,7)=estif(7,9)
      estif(10,7)=estif(7,10)
      estif(11,7)=estif(7,11)
      estif(13,7)=estif(7,13)
      estif(14,7)=estif(7,14)
      estif(15,7)=estif(7,15)
      estif(10,9)=estif(9,10)
      estif(11,9)=estif(9,11)
      estif(13,9)=estif(9,13)
      estif(14,9)=estif(9,14)
      estif(15,9)=estif(9,15)
      estif(11,10)=estif(10,11)
      estif(13,10)=estif(10,13)
      estif(14,10)=estif(10,14)
      estif(15,10)=estif(10,15)
      estif(13,11)=estif(11,13)
      estif(14,11)=estif(11,14)
      estif(15,11)=estif(11,15)
      estif(14,13)=estif(13,14)
      estif(15,13)=estif(13,15)
      estif(15,14)=estif(14,15)

      return


c     --------------------------------------------------------
c     Evaluate the elemantal RHS/residual from a source term
c     --------------------------------------------------------
 5    continue
      eforc(:) = 0.0d0! elemental RHS/residual
c     turned off for now
      return


c     --------------------------------------------------------
c     Build the elemental RHS/residual (eforc) for the dual problem
c     --------------------------------------------------------
 6    l  = int(pelem(1))
      ndim = ndime

      eforc(:) = 0.0d0! elemental RHS/residual

      call gausstet(l,ninte,sg,tg,zg,wg)

      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,xjaco,
     $              .false.,nnode,ndim,elnods,xelem)
        wtjac = wg(iinte)*xjaco

        ievab = 0
        do inode = 1, nnode
          do idofn = 1,4! 4 = elemvec_ndofn(inode)
            ievab = ievab+1
            do jnode = 1,nnode
c             do jdofn = 1,4! 4 = elemvec_ndofn(jnode)
c               if (idofn.eq.jdofn) then
c                 eforc(ievab) = eforc(ievab)-
c     $                 shap(4,inode)*shap(4,jnode)
c     $                *uelem_diff(jdofn,jnode)*wtjac
c                endif
c              enddo! jdofn
              jdofn=idofn
              eforc(ievab) = eforc(ievab)-
     $                 shap(4,inode)*shap(4,jnode)
     $                *uelem_diff(jdofn,jnode)*wtjac! twice jdofn ... ask Sevan JFD
            enddo! jnode
          enddo! idofn
        enddo! inode
      enddo! iinte

      return


c     --------------------------------------------------------
c     Compute the objective function (dataMatch+regularization)
c      and its gradient (egrad) on the element
c     --------------------------------------------------------
 7    l       = int(pelem(1))
      ireg    = int(pelem(2))
      alpha(1)= pelem(3)
      alpha(2)= pelem(4)
      beta(1) = pelem(5)
      beta(2) = pelem(6)
      Treg(1) = pelem(7)
      Treg(2) = pelem(8)
      tauMult = pelem(9)
      s       = int(pelem(10))

      elemDataMatch(ielem) = 0.0d0
      elemRegul(ielem) = 0.0d0
      l2grad1(ielem) = 0.0d0
      l2grad2(ielem) = 0.0d0
      egrad(:,:) = 0.0d0! initialize egrad

c     identity matrix ident
      ident(:,:) = 0.0d0
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0
      ndim = ndime

      powe = 0.5d0

      temp_dual(:) = 0.0d0
      ii=0

      do inode = 1,nnode
        do idofn = 1,4! 4 = elemvec_ndofn(inode)
          ii = ii+1
          temp_dual(ii)   = uelem_dual(idofn,inode)
        enddo
      enddo

c     determine the characteristic element length h
c       (for tetrahedra: h is the length of the longest edge)
      tmp1=(xelem(1,1)-xelem(1,2))**2.0d0
      tmp2=(xelem(1,1)-xelem(1,3))**2.0d0
      tmp3=(xelem(1,1)-xelem(1,4))**2.0d0
      tmp4=(xelem(1,2)-xelem(1,3))**2.0d0
      tmp5=(xelem(1,2)-xelem(1,4))**2.0d0
      h=(xelem(1,3)-xelem(1,4))**2.0d0
      do i=2,ndim! 1, ndim
        tmp1 = tmp1 + (xelem(i,1)-xelem(i,2))**2.0d0
        tmp2 = tmp2 + (xelem(i,1)-xelem(i,3))**2.0d0
        tmp3 = tmp3 + (xelem(i,1)-xelem(i,4))**2.0d0
        tmp4 = tmp4 + (xelem(i,2)-xelem(i,3))**2.0d0
        tmp5 = tmp5 + (xelem(i,2)-xelem(i,4))**2.0d0
        h = h + (xelem(i,3)-xelem(i,4))**2.0d0
      enddo
      if (tmp1.gt.h) then
        h=tmp1
      endif
      if (tmp2.gt.h)then
        h=tmp2
      endif
      if (tmp3.gt.h) then
        h=tmp3
      endif
      if (tmp4.gt.h) then
        h=tmp4
      endif
      if (tmp5.gt.h) then
        h=tmp5
      endif
      h=sqrt(h)
      
      call gausstet(l,ninte,sg,tg,zg,wg)

      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,xjaco,
     $              .false.,nnode,ndim,elnods,xelem)
        wtjac = wg(iinte)*xjaco

c       create prop_grad (contains prop gradient and value) 
        prop_grad(:,:)  = 0.0d0
        do iset = 1,2
          do idofn = 1,4      !three gradients+the value
            do inode = 1,nnode
               prop_grad(iset,idofn) = prop_grad(iset,idofn)+
     $                            shap(idofn,inode)*
     $                      elemdata_nodal(iset,inode)
            enddo
          enddo
        end do

c       compute the gradient and the objective function

        udiff(:)      = 0.0d0

c       compute the deformation gradient at Gauss Point
        Fdef(:,:) = ident(:,:)
        do i = 1,ndim
          do j = 1,ndim
            do inode = 1,nnode
              Fdef(i,j)=Fdef(i,j)+uelem(i,inode)*shap(j,inode)
            end do
          end do
        end do
c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)*Fdef(3,3) +
     $            Fdef(1,2)*Fdef(2,3)*Fdef(3,1) +
     $            Fdef(1,3)*Fdef(2,1)*Fdef(3,2) -
     $            Fdef(1,3)*Fdef(2,2)*Fdef(3,1) -
     $            Fdef(1,2)*Fdef(2,1)*Fdef(3,3) -
     $            Fdef(1,1)*Fdef(2,3)*Fdef(3,2)


c       Finv is inverse of deformation gradient at Gauss Point
        Finv(1,1) = Fdef(2,2)*Fdef(3,3)-Fdef(2,3)*Fdef(3,2)
        Finv(1,2) = Fdef(1,3)*Fdef(3,2)-Fdef(1,2)*Fdef(3,3)
        Finv(1,3) = Fdef(1,2)*Fdef(2,3)-Fdef(1,3)*Fdef(2,2)
        Finv(2,1) = Fdef(2,3)*Fdef(3,1)-Fdef(2,1)*Fdef(3,3)
        Finv(2,2) = Fdef(1,1)*Fdef(3,3)-Fdef(1,3)*Fdef(3,1)
        Finv(2,3) = Fdef(1,3)*Fdef(2,1)-Fdef(1,1)*Fdef(2,3)
        Finv(3,1) = Fdef(2,1)*Fdef(3,2)-Fdef(2,2)*Fdef(3,1)
        Finv(3,2) = Fdef(1,2)*Fdef(3,1)-Fdef(1,1)*Fdef(3,2)
        Finv(3,3) = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)
        Finv(:,:) = (1.0d0/Fdet)*Finv(:,:)

c       compute pressure and spatial derivatives at Gauss Point
        pres(1:4) = 0.0d0
        do i = 1, 4
          do inode = 1, nnode
            pres(i) = pres(i) + uelem(4,inode)*shap(i,inode)
          end do
        end do

c       compute the Cauchy tensor and its inverse
        Ctens(1:3,1:3) = 0.0d0
        do j = 1,ndim
          do i = j,ndim
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do
c       enforce the symmetries
        Ctens(1,2)=Ctens(2,1)
        Ctens(1,3)=Ctens(3,1)
        Ctens(2,3)=Ctens(3,2)

        Cdet = Ctens(1,1)*Ctens(2,2)*Ctens(3,3) +
     $            Ctens(1,2)*Ctens(2,3)*Ctens(3,1) +
     $            Ctens(1,3)*Ctens(2,1)*Ctens(3,2) -
     $            Ctens(1,3)*Ctens(2,2)*Ctens(3,1) -
     $            Ctens(1,2)*Ctens(2,1)*Ctens(3,3) -
     $            Ctens(1,1)*Ctens(2,3)*Ctens(3,2)


        Cinv(1,1) = Ctens(2,2)*Ctens(3,3)-Ctens(2,3)*Ctens(3,2)
        Cinv(1,2) = Ctens(1,3)*Ctens(3,2)-Ctens(1,2)*Ctens(3,3)
        Cinv(1,3) = Ctens(1,2)*Ctens(2,3)-Ctens(1,3)*Ctens(2,2)
        Cinv(2,1) = Cinv(1,2)! = Ctens(2,3)*Ctens(3,1)-Ctens(2,1)*Ctens(3,3)
        Cinv(2,2) = Ctens(1,1)*Ctens(3,3)-Ctens(1,3)*Ctens(3,1)
        Cinv(2,3) = Ctens(1,3)*Ctens(2,1)-Ctens(1,1)*Ctens(2,3)
        Cinv(3,1) = Cinv(1,3)! = Ctens(2,1)*Ctens(3,2)-Ctens(2,2)*Ctens(3,1)
        Cinv(3,2) = Cinv(2,3)! = Ctens(1,2)*Ctens(3,1)-Ctens(1,1)*Ctens(3,2)
        Cinv(3,3) = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)
        Cinv(:,:) = (1.0d0/Cdet)*Cinv(:,:)

c       principal invariant Inv1
        Inv1 = Ctens(1,1)+Ctens(2,2)+Ctens(3,3)

c       dJdC is the derivative of the Jacobian with respect to 
c       the Cauchy Green tensor
        dJdC(1:3,1:3) = 0.5d0*Fdet*Cinv(1:3,1:3)

        do iset = 1,2
          do knode = 1,nnode

c           set variations in the shape function for the property 
            var_shape_prop(:,:) = 0.0d0
            if (ielem_nod_grad(iset,knode).eq.1) then
              do idofn = 1,4
                var_shape_prop(iset,idofn) = shap(idofn,knode)
              enddo
            endif

            if (s.eq.0) then! JFD what is s????
              temp  = 0.0d0
              temp_der  = 0.0d0
            end if 
            if (s.eq.1) then! JFD make it a elseif
              temp = 0.5d0*(tauMult*(h**2.0d0))/prop_grad(2,4)
              temp_der = 0.0d0
            end if
            if (s.eq.1 .AND. iset.eq.2) then! JFD include in s.eq.1
              temp_der = -0.5d0*var_shape_prop(2,4)*
     $                (tauMult*(h**2.0d0))/(prop_grad(2,4)**2.0d0)
            end if

            Tp(1:4) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                Tp(i) = Tp(i) + 2.0d0*temp_der*dJdC(j,i)*pres(j)
              end do
            end do

c           Tp(4) = 0.0d0 is zero here

c           compute derivative of second Piola Kirchhoff stress SecPK 
c           with respect to the material parameters
         
            K1 = (Fdet**(-2.0d0/3.0d0))*Inv1-3.0d0
            K2(1:3,1:3) = (ident(1:3,1:3) - 
     $                   (1.0d0/3.0d0)*Cinv(1:3,1:3)*Inv1)
     $                   *(Fdet**(-2.0d0/3.0d0))

            if (iset.eq.1) then
              dWdC_grad(1:3,1:3) = 0.5d0*prop_grad(2,4)*
     $                         K1*var_shape_prop(iset,4)*
     $                         exp(K1*prop_grad(iset,4))*K2(1:3,1:3)
            else if (iset.eq.2) then
              dWdC_grad(1:3,1:3) = 0.5d0*var_shape_prop(iset,4)
     $                         *exp(K1*prop_grad(1,4))*K2(1:3,1:3)
            end if

            SecPK_grad(1:3,1:3) = 2.0d0*dWdC_grad(1:3,1:3)

            FS(:,:) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                do r = 1, ndim
                  FS(i,j) = FS(i,j)+Fdef(i,r)*SecPK_grad(r,j)
                enddo
              enddo
            enddo
         
            ievab = 0
            eforc(:) = 0.0d0
            do inode = 1,nnode
              do i = 1,4! 4 = elemvec_ndofn(inode)
                ievab = ievab+1
                if (i.eq.4) then! 4 = elemvec_ndofn(inode)
                  do j = 1,4
                    eforc(ievab)=eforc(ievab)+
     $                          shap(j,inode)*Tp(j)*wtjac
                  enddo
                else
                  do j = 1,3
                    eforc(ievab)=eforc(ievab)+
     $                          shap(j,inode)*FS(i,j)*wtjac
                  enddo
                endif
              enddo
            enddo! inode
 
            do ii = 1, 16
              egrad(iset,knode) = egrad(iset,knode) +
     $                            eforc(ii)*temp_dual(ii)
            enddo

c           account for the regularization term
c JFD            if (iset.eq.1) then
              if (ireg.eq.0) then
c               there is not supposed to be a regularization term: do nothing
              elseif (ireg.eq.1) then! H1 or Tikhonov reg.
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))
              elseif ((ireg.eq.2).or.(ireg.eq.21)) then! TVD reg.
                deno = sqrt(beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
     $              +prop_grad(iset,3)*prop_grad(iset,3))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $              var_shape_prop(iset,2)*prop_grad(iset,2)+
     $              var_shape_prop(iset,3)*prop_grad(iset,3))/deno
              elseif (ireg.eq.3) then! power reg.
                deno = sqrt(beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
     $              +prop_grad(iset,3)*prop_grad(iset,3))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              ((1+deno-beta(iset))**(powe-1.0d0))*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))/deno
              elseif (ireg.eq.31) then! power reg. (proposed by Paul Barbone)
                deno = beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
     $             +prop_grad(iset,3)*prop_grad(iset,3)
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))*
     $              (deno**(1-0.5d0*powe)-(deno-beta(iset)*beta(iset))*
     $            (1-0.5d0*powe)*(deno**(-0.5d0*powe)))/(deno**(2-powe))
              elseif (ireg.eq.4) then! logarithm reg.
                deno = sqrt(beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
     $              +prop_grad(iset,3)*prop_grad(iset,3))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))/(deno*
     $              Treg(iset)*(1+(deno-beta(iset))/Treg(iset)))
              elseif (ireg.eq.41) then! logarithm reg. (second implementation)
                deno = Treg(iset)*Treg(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2)+
     $              prop_grad(iset,3)*prop_grad(iset,3)
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))/deno
              elseif (ireg.eq.5) then! exponential reg. (contrast preserving)
                deno = sqrt(beta(iset)*beta(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2)+
     $              prop_grad(iset,3)*prop_grad(iset,3))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))*
     $              exp((beta(iset)-deno)/Treg(iset))/(Treg(iset)*deno)
              elseif (ireg.eq.6) then! fraction reg. (contrast preserving)
                deno = 1 + (prop_grad(iset,1)*prop_grad(iset,1)+
     $                      prop_grad(iset,2)*prop_grad(iset,2)+
     $                      prop_grad(iset,3)*prop_grad(iset,3))/
     $                     (Treg(iset)*Treg(iset))
                deno = Treg(iset)*Treg(iset)*deno*deno
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2)+
     $               var_shape_prop(iset,3)*prop_grad(iset,3))/deno
              else
                Print*,"elem305.f: ireg=",ireg,
     $                 " is not implemented: exiting"
                stop
              endif! ireg
c            elseif (iset.eq.2) then
c              if (ireg.eq.1) then
c                egrad(iset,knode) = egrad(iset,knode) + alpha(iset)*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2)+    
c     $              var_shape_prop(iset,3)*prop_grad(iset,3))*
c     $              wtjac
c              elseif (ireg.eq.2) then
c                deno = sqrt(
c     $              beta(iset)*beta(iset)+
c     $              prop_grad(2,1)*prop_grad(2,1)+
c     $              prop_grad(2,2)*prop_grad(2,2)+
c     $              prop_grad(2,3)*prop_grad(2,3))
c                egrad(iset,knode) = egrad(iset,knode) + 
c     $              alpha(iset)*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2)+
c     $              var_shape_prop(iset,3)*prop_grad(iset,3))*
c     $              wtjac/deno
c              endif! ireg
c            endif! iset
          enddo! knode
        enddo! iset
      
        do j=1,ndim
          do inode=1,nnode
            udiff(j) = udiff(j)+shap(4,inode)*uelem_diff(j,inode)
          enddo
        enddo
      
        elemDataMatch(ielem) = elemDataMatch(ielem) + 0.5d0*wtjac*(
     $        udiff(1)*udiff(1)+udiff(2)*udiff(2)+udiff(3)*udiff(3))

        if (ireg.eq.1) then! H1 or Tikhonov reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(2)*(prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)+prop_grad(2,3)*prop_grad(2,3))+
     $      alpha(1)*(prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3)))
        elseif (ireg.eq.2) then! TVD reg. with offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)+prop_grad(2,3)*prop_grad(2,3))+
     $      alpha(1)*sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3)))
        elseif (ireg.eq.21) then! TVD reg. without offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*(sqrt(beta(2)*beta(2)+prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2)+
     $           prop_grad(2,3)*prop_grad(2,3))-beta(2))+
     $      alpha(1)*(sqrt(beta(1)*beta(1)+prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))-beta(1)))
        elseif (ireg.eq.3) then! power reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)/powe*((1+sqrt(beta(2)*beta(2)+prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2)+
     $           prop_grad(2,3)*prop_grad(2,3))-beta(2))**(powe)-1.0d0)+
     $      alpha(1)/powe*((1+sqrt(beta(1)*beta(1)+prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))-beta(1))**(powe)-1.0d0))
        elseif (ireg.eq.31) then! power regularization (proposed by Paul Barbone)
c JFD: for powe.eq.1, this is another implementation of TVD
          deno = prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)+prop_grad(2,3)+prop_grad(2,3)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(2)*deno/(beta(2)*beta(2)+deno)**(1-0.5d0*powe))
          deno = prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(1)*deno/(beta(1)*beta(1)+deno)**(1-0.5d0*powe))
        elseif (ireg.eq.4) then! logarithm reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $       (alpha(2)*log(1+(sqrt(beta(2)*beta(2)+prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2)+
     $           prop_grad(2,3)*prop_grad(2,3))-beta(2))/Treg(2))+
     $       alpha(1)*log(1+(sqrt(beta(1)*beta(1)+prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))-beta(1))/Treg(1)))
        elseif (ireg.eq.41) then! logarithm reg. (second implementation to avoid using a sqrt operator and the associated beta that has a significant influence on the convergence speed)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $       (alpha(2)*log(1+(prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2)+prop_grad(2,3)*
     $           prop_grad(2,3))/(Treg(2)*Treg(2)))+
     $       alpha(1)*log(1+(prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)+prop_grad(1,3)*
     $           prop_grad(1,3))/(Treg(1)*Treg(1))))
        elseif (ireg.eq.5) then! exponential reg. (contrast preserving: no penalty for large jumps)
c JFD: proposed by Paul Barbone
c JFD: this regularization is very unstable
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $       (alpha(2)*(1.0d0 - exp((beta(2)-sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $         prop_grad(2,2)+prop_grad(2,3)*prop_grad(2,3)))/Treg(2)))+
     $       alpha(1)*(1.0d0-exp((beta(1)-sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $         prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3)))/Treg(1))))
        elseif (ireg.eq.6) then! fraction reg. (contrast preserving: no penalty for large jumps)
c JFD: this regularization is very unstable
          deno = (prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)+prop_grad(2,3)*prop_grad(2,3))/
     $           (Treg(2)*Treg(2))
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $                     (alpha(2)*deno/(1+deno))
          deno = (prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)+prop_grad(1,3)*prop_grad(1,3))/
     $           (Treg(1)*Treg(1))
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $                     (alpha(1)*deno/(1+deno))
        endif
        l2grad1(ielem) = l2grad1(ielem) + sqrt(prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))
        l2grad2(ielem) = l2grad2(ielem) + sqrt(prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2)+
     $           prop_grad(2,3)*prop_grad(2,3))
      enddo! iinte
      dataMatch=dataMatch+elemDataMatch(ielem)
      regularization = regularization + elemRegul(ielem)
      l2grad1(ielem) = l2grad1(ielem)/dble(ninte)
      l2grad2(ielem) = l2grad2(ielem)/dble(ninte)
      suml2grad1 = suml2grad1 + l2grad1(ielem)*l2grad1(ielem)
      suml2grad2 = suml2grad2 + l2grad2(ielem)*l2grad2(ielem)

      return


 4    return
 8    return
 9    return
 10   return
 11   return
 12   return
 13   return
 14   return
 15   return
 16   return
 17   return

      end

