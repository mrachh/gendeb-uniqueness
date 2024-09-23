      implicit real *8 (a-h,o-z)
      complex *16 zeps(2),zmu(2),zepssq(2),zmusq(2),zkp,zkm
      real *8 domeg
      complex *16 zdet,zdet2,zvals(4)
      complex *16 ima,z,z1,z2,z3,z4,z0,zroot1,zroot2,tuse,zroot1_pde
      real *8 err_est
      data ima/(0.0d0,1.0d0)/
      call prini(6,13)

      done = 1
      pi = atan(done)*4

      

      nn = 5
      rfac = 1.1838d0

      nmax = 1000
      tmin = 0.0d0
      vmin = 1.0d7
 1211 format(4(2x,e11.5))
      do i=1,nmax
        tt = 10.0d0*(i+0.0d0)/(nmax+0.0d0)
        tuse = tt
        call get_det_plasmon(tuse,rfac,nn,zdet,zdet2,err_est)
        write(111,1211) tt,real(zdet2),imag(zdet2),abs(zdet2)
        if(abs(zdet2).lt.vmin) then
          vmin = abs(zdet2)
          tmin = tt
        endif
      enddo
      call prin2('tmin=*',tmin,1)

      z0 = tmin
      z1 = tmin - ima*1.0d-2
      call prinf('nn=*',nn,1)
      call prin2('rfac=*',rfac,1)
      call find_plasmon_root(z0,z1,rfac,nn,zroot1)
      call prin2('zroot1=*',zroot1,2)
      call prin2_long('zroot1=*',zroot1,2)
      call get_det_plasmon(zroot1,rfac,nn,zdet,zdet2,err_est)
      call prin2('zdet at root=*',zdet,2)
      call prin2('zdet2 at root=*',zdet2,2)
      call prin2('error=*',err_est,1)
      call prin2('rfac=*',rfac,1)

      zroot1 = real(zroot1)

      call test_pde_resonance(zroot1,rfac,nn,zvals)
      call prin2('zvals=*',zvals,8)
      stop

      iind = 1
      z1 = tmin - ima*1.0d-2
      call find_pde_root(z0,z1,rfac,nn,iind,zroot1_pde)
      call prin2('zroot1_pde=*',zroot1_pde,2)
      stop

      z0 = tmin
      z1 = tmin - ima*1.0d-2
      call find_plasmon_root(z0,z1,rfac,nn,zroot2)
      call prin2('zroot2=*',zroot2,2)

      stop
      end
c
c
c
c
      subroutine find_pde_root(z0,z1,rfac,nn,iind,z)
      implicit real *8 (a-h,o-z)
      complex *16 z0,z1,z2,f0,f1,f2,ztmp,z,zvals(4)

      call test_pde_resonance(z0,rfac,nn,zvals)
      f0 = zvals(iind)

      call test_pde_resonance(z1,rfac,nn,zvals)
      f1 = zvals(iind)
      maxit = 100
      thresh = 1.0d-12
      ifdone = 0
      do iter = 1,maxit
        z2 = z1 - f1*(z1-z0)/(f1-f0)
        z = z2
        call test_pde_resonance(z2,rfac,nn,zvals)
        f2 = zvals(iind)
cc        call prinf('iter=*',iter,1)
cc        call prin2('z2=*',z2,2)
cc        call prin2('f2=*',f2,2)
        
        if(abs(f2).lt.thresh) then
          ifdone = ifdone + 1
        endif
        if(ifdone.gt.3) then
          call prin2('fun val at root=*',abs(f2),1)
          return
        endif
        z0 = z1
        f0 = f1
        z1 = z2
        f1 = f2
      enddo
      call prin2('error max iteration reached*',i,0)
      call prin2('last fun val=*',abs(f2),1)





      return
      end

c
c
c
c
c
c
c
      subroutine find_plasmon_root(z0,z1,rfac,nn,z)
      implicit real *8 (a-h,o-z)
      complex *16 z0,z1,z2,f0,f1,f2,ztmp,z

      call get_det_plasmon(z0,rfac,nn,ztmp,f0,err_est)
      call get_det_plasmon(z1,rfac,nn,ztmp,f1,err_est)
      maxit = 100
      thresh = 1.0d-12
      ifdone = 0
      do iter = 1,maxit
        if(abs(z1-z0).le.1.0d-14) return
        z2 = z1 - f1*(z1-z0)/(f1-f0)
        z = z2
        call get_det_plasmon(z2,rfac,nn,ztmp,f2,err_est)
        call prinf('iter=*',iter,1)
        call prin2('z2=*',z2,2)
        call prin2('f2=*',f2,2)
        
        if(abs(f2).lt.thresh) then
          ifdone = ifdone + 1
        endif
        if(ifdone.gt.3) then
          call prin2('fun val at root=*',abs(f2),1)
          return
        endif
        z0 = z1
        f0 = f1
        z1 = z2
        f1 = f2
      enddo
      call prin2('error max iteration reached*',i,0)
      call prin2('last fun val=*',abs(f2),1)





      return
      end

c
c
c
c
c

      subroutine get_det_plasmon(tt,rfac,nn,zdet,zdet2,err_est)
c
c   This subroutine computes the determinant of the 4 x4 matrix
c   corresponding to discretization of dielectric problem
c   using generalized debye representation for spherical
c   harmonic mode nn
c
c   The material parameters (ep_{+},\mu_{+},ep_{-},\mu_{-})
c   are  (t,t,i*Sqrt(rfac)*t, t)
c
c
      implicit real *8 (a-h,o-z)
      complex *16 zeps(2),zmu(2),zepssq(2),zmusq(2),zkp,zkm
      complex *16 zdet,zdet2,tt
      complex *16 ima
      real *8 domeg
      data ima/(0.0d0,1.0d0)/


      zkp = tt
      zkm = ima*sqrt(rfac)*tt
cc      zkm = tt
cc      zkp = ima*sqrt(rfac)*tt
      zmu(1) = zkp
      zmu(2) = zkp
      zeps(1) = zkm**2/zkp
      zeps(2) = zkp

      zmusq(1) = sqrt(zkp)
      zmusq(2) = sqrt(zkp)

      zepssq(1) = zkm/sqrt(zkp)
      zepssq(2) = sqrt(zkp)
      domeg = 1.0d0

      call get_det(zeps,zepssq,zmu,zmusq,
     1   domeg,nn,zdet,zdet2,err_est)

      print *, zeps(1) + zeps(2)
      print *, zmu(1) + zmu(2)

      return
      end
c
c
c
c

      subroutine test_pde_resonance(tt,rfac,nn,zvals)
c
c   This subroutine computes the denominators of mie scattering
c    on the unit sphere
c
c   The material parameters (ep_{+},\mu_{+},ep_{-},\mu_{-})
c   are  (t,t,i*Sqrt(rfac)*t, t)
c
c
      implicit real *8 (a-h,o-z)
      complex *16 zeps(2),zmu(2),zepssq(2),zmusq(2),zkp,zkm
      complex *16 zdet,zdet2,tt,zvals(4)
      complex *16, allocatable :: fjvals(:,:),fjder(:,:)
      complex *16, allocatable :: fhvals(:,:),fhder(:,:)
      complex *16 zk(2)
      complex *16 ima
      real *8 domeg
      data ima/(0.0d0,1.0d0)/


      zkp = tt
      zkm = ima*sqrt(rfac)*tt
cc      zkm = tt
cc      zkp = ima*sqrt(rfac)*tt
      zmu(1) = zkp
      zmu(2) = zkp
      zeps(1) = zkm**2/zkp
      zeps(2) = zkp

      zmusq(1) = sqrt(zkp)
      zmusq(2) = sqrt(zkp)

      zepssq(1) = zkm/sqrt(zkp)
      zepssq(2) = sqrt(zkp)
      domeg = 1.0d0


      ifder = 1
      rscale = 1.0d0
      
      njh = nn + 10
      allocate(fjvals(0:njh,2),fjder(0:njh,2))
      allocate(fhvals(0:njh,2),fhder(0:njh,2))

      do ii=1,2
        if(ii.eq.1) i = 2
        if(ii.eq.2) i = 1
        zk(i) = domeg*zepssq(i)*zmusq(i)

        call besseljs3d(njh,zk(i),rscale,fjvals(0,i),ifder,fjder(0,i))
        call h3dall(njh,zk(i),rscale,fhvals(0,i),ifder,fhder(0,i))
      enddo

      zvals(1) = zmu(1)*(fhvals(nn,2)+zk(2)*fhder(nn,2))*fjvals(nn,1)-
     1    zmu(2)*(fjvals(nn,1) + zk(1)*fjder(nn,1))*fhvals(nn,2) 
      zvals(3) = zmu(2)*(fhvals(nn,1)+zk(1)*fhder(nn,1))*fjvals(nn,2)-
     1    zmu(1)*(fjvals(nn,2) + zk(2)*fjder(nn,2))*fhvals(nn,1) 
      zvals(2) = zeps(1)*(fhvals(nn,2)+zk(2)*fhder(nn,2))*fjvals(nn,1)-
     1    zeps(2)*(fjvals(nn,1) + zk(1)*fjder(nn,1))*fhvals(nn,2) 
      zvals(4) = zeps(2)*(fhvals(nn,1)+zk(1)*fhder(nn,1))*fjvals(nn,2)-
     1    zeps(1)*(fjvals(nn,2) + zk(2)*fjder(nn,2))*fhvals(nn,1) 



      return
      end
c
c
c
c
c
c
c
c
      subroutine test_solver_get_lamsing(zeps,zmu,domeg,nnmax,
     1   zlamall,singall,errall)

      implicit real *8 (a-h,o-z)
      complex *16, allocatable :: fjvals(:,:),fhvals(:,:)
      complex *16, allocatable :: fjder(:,:),fhder(:,:)
      complex *16 zeps(2),zmu(2),zk(2)
      complex *16 zalpha(2),zbeta(2),zg1(2),zg2(2),zeta(2),za1(2),zd(2)
      complex *16 za2(2)
      complex *16 ima
      complex *16 zmat(4,4),zrhs(4),zsoln(4)
      complex *16 zmat2(4,4),zmat0(4,4),zmat0inv(4,4)
      complex *16 zeu(2),zex(2),zey(2)
      complex *16 zeu0(2),zex0(2),zey0(2)
      complex *16 zbu(2),zbx(2),zby(2)
      complex *16 zbu0(2),zbx0(2),zby0(2)
      complex *16 zsingl(4,4),zsingr(4,4)
      complex *16 zeigl(4,4),zeigr(4,4),zlam(4)
      complex *16 zlamall(4,nnmax)
      real *8 singall(4,nnmax),errall(nnmax)
      real *8 sing(4)

      data ima/(0.0d0,1.0d0)/
      call prini(6,13)
c
c
c  1 - interior
c  2 - exterior
c


      ifder = 1
      rscale = 1.0d0
      
      njh = nnmax + 10
      allocate(fjvals(0:njh,2),fjder(0:njh,2))
      allocate(fhvals(0:njh,2),fhder(0:njh,2))

      do i=1,2
        zk(i) = domeg*sqrt(zeps(i))*sqrt(zmu(i))

        call besseljs3d(njh,zk(i),rscale,fjvals(0,i),ifder,fjder(0,i))
        call h3dall(njh,zk(i),rscale,fhvals(0,i),ifder,fhder(0,i))
      enddo


      do nn=1,nnmax
        do i=1,2
          zalpha(i) = -ima*((nn+1)*fjvals(nn,i)*fhvals(nn-1,i)+
     1       nn*fjvals(nn+1,i)*fhvals(nn,i)-
     2       zk(i)*fjvals(nn+1,i)*fhvals(nn-1,i))
          zbeta(i) = fjvals(nn,i)*fhvals(nn,i)*zk(i)/ima
          za1(i) = ima*zk(i)**2*fjvals(nn,i)*fhder(nn,i)
          zg1(i) = ima*((fhvals(nn,i)+
     1       zk(i)*fhder(nn,i))*fjvals(nn,i)*zk(i)) 
          za2(i) = ima*zk(i)**2*fhvals(nn,i)*fjder(nn,i)
          zg2(i) = ima*((fjvals(nn,i)+
     1       zk(i)*fjder(nn,i))*fhvals(nn,i)*zk(i)) 
          zd(i) = -ima*zk(i)*fjvals(nn,i)*fhvals(nn,i)*
     1       sqrt(nn*(nn+1.0d0)) 
          zeta(i) =ima*sqrt(nn*(nn+1.0d0))*(fjvals(nn,i)*fhvals(nn-1,i)-
     1     fjvals(nn+1,i)*fhvals(nn,i))
        enddo
c
c  Now assemble the matrix
c
c  unknowns ordering r_{-} r_{+} q_{-} q_{+}
c  Equations ordering 
c   1. G_{0} \nabla_{s} (E_{+} - E^{-})_{tan} = G_{0} \nabla_{s} (E^{in})  
c   2. G_{0} \nabla_{s} (B_{+} - B^{-})_{tan} = G_{0} \nabla_{s} (B^{in})  
c   3. \eps_{+} E^{+}.n - \eps_{-} E^{-} \cdot n = \eps_{+} E^{in} \cdot n
c   4. \mu_{+} B^{+}.n - \mu_{-} B^{-} \cdot n = \mu_{+} B^{in} \cdot n
c
c
c   Note sign flipped on the right hand side in the above equation
c   by design since we are doing an analytical test and setting
c   E^{in} = E^{+} for some known E^{+} and similarly for B^{+}
c
        zmat(1,1) = -ima*zk(1)*sqrt(zmu(1))*ima*zk(1)*zalpha(1) - 
     1     sqrt(zmu(1))*zd(1)*sqrt(nn*(nn+1.0d0)) 
        zmat(1,2) = ima*zk(2)*sqrt(zmu(2))*ima*zk(2)*zalpha(2) + 
     1     sqrt(zmu(2))*zd(2)*sqrt(nn*(nn+1.0d0))
        zmat(1,3) = -zg1(2)*ima*domeg*zmu(1)*sqrt(zeps(1))
        zmat(1,4) = -zg2(1)*ima*domeg*zmu(2)*sqrt(zeps(2))
        zmat(1,1:4) = zmat(1,1:4)/(2*nn+1.0d0)
c
c
c
c
        zmat(2,1) = -zg1(2)*ima*domeg*zeps(1)*sqrt(zmu(1))
        zmat(2,2) = -zg2(1)*ima*domeg*zeps(2)*sqrt(zmu(2))
        zmat(2,3) = ima*zk(1)*sqrt(zeps(1))*ima*zk(1)*zalpha(1) + 
     1     sqrt(zeps(1))*zd(1)*sqrt(nn*(nn+1.0d0))
        zmat(2,4) = -ima*zk(2)*sqrt(zeps(2))*ima*zk(2)*zalpha(2) - 
     1     sqrt(zeps(2))*zd(2)*sqrt(nn*(nn+1.0d0))
        zmat(2,1:4) = zmat(2,1:4)/(2*nn+1.0d0)
c
c
c
c
        zmat(3,1) = zeps(1)*sqrt(zmu(1))*
     1      ((ima*zk(1))**2/sqrt(nn*(nn+1.0d0))*zeta(1)+za2(1))
        zmat(3,2) = -zeps(2)*sqrt(zmu(2))*
     1     ((ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*zeta(2)+za1(2))
        zmat(3,3) = zeps(2)*ima*domeg*zmu(1)*sqrt(zeps(1))*zbeta(2)
        zmat(3,4) = zeps(1)*ima*domeg*zmu(2)*sqrt(zeps(2))*zbeta(1)
c
c
c
c
        zmat(4,1) = zmu(2)*ima*domeg*zeps(1)*sqrt(zmu(1))*zbeta(2)
        zmat(4,2) = -zmu(1)*ima*domeg*zeps(2)*sqrt(zmu(2))*zbeta(1)
        zmat(4,3) = -zmu(1)*sqrt(zeps(1))*
     1      ((ima*zk(1))**2/sqrt(nn*(nn+1.0d0))*zeta(1)+za2(1))
        zmat(4,4) = zmu(2)*sqrt(zeps(2))*
     1     ((ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*zeta(2)+za1(2))


        
        zrhs(1) = 0
        zrhs(2) =-sqrt(nn*(nn+1.0d0))*zg1(2)/(2*nn+1.0d0)*sqrt(zeps(2)) 
        zrhs(3) = 0
        zrhs(4) = sqrt(nn*(nn+1.0d0))*zbeta(2)*zmu(2)*sqrt(zeps(2))

        zsoln(1) = 0
        zsoln(2) = 0
        zsoln(3) = 0
        zsoln(4) = 0

        call zgausselim(4,zmat,zrhs,info,zsoln,dcond)
c
c
c  Test exterior electric field
c

        rp = 0
        zeu0(2) = 0
        zex0(2) = -zbeta(2)*ima*zk(2)*sqrt(zmu(2))
        zey0(2) = 0

        zeu(2) = 0
        zex(2) = -ima*zk(2)*sqrt(zmu(2))*ima*domeg*zeps(1)*sqrt(zmu(1))/
     1     sqrt(zeps(2))*zsoln(1)/sqrt(nn*(nn+1.0d0))*zbeta(2) + 
     2     sqrt(zmu(2))*ima*zk(2)*zsoln(4)/sqrt(nn*(nn+1.0d0))*zg2(2)
        zey(2) = 0

        errep = abs(zex(2)-zex0(2))
        rp = abs(zex0(2))
        

c
c
c  test exterior magnetic field
c

        zbx0(2) = 0
        zby0(2) = sqrt(zeps(2))*sqrt(nn*(nn+1.0d0))*zbeta(2)
        zbu0(2) = -sqrt(zeps(2))*zg1(2)

        zbx(2) = 0
        zby(2) = sqrt(zeps(2))*(ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*
     1     zsoln(4)*zeta(2) + ima*domeg*zeps(1)*sqrt(zmu(1))*
     2     zsoln(1)*zbeta(2) + za1(2)*sqrt(zeps(2))*zsoln(4)
        zbu(2) = -sqrt(zeps(2))*(ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*
     1     zsoln(4)*zalpha(2) - ima*domeg*zeps(1)*sqrt(zmu(1))*
     2     zsoln(1)*zg1(2)/sqrt(nn*(nn+1.0d0)) - 
     3     zsoln(4)*sqrt(zeps(2))*zd(2)


        errbp = sqrt(abs(zby0(2)-zby(2))**2 + abs(zbu0(2)-zbu(2))**2)
        rbp = sqrt(abs(zby0(2))**2+abs(zbu0(2))**2)
        errall(nn) = sqrt(errbp**2 + errep**2)/sqrt(rbp**2 + rp**2) 
c
c  diagonal precondition
c
c

        do i=1,4
          do j=1,4
            zmat0(j,i) = 0
          enddo
        enddo
        zmat0(1,1) = sqrt(zmu(1))/4
        zmat0(1,2) = -sqrt(zmu(2))/4

        zmat0(2,3) = -sqrt(zeps(1))/4
        zmat0(2,4) = sqrt(zeps(2))/4
      
        zmat0(3,1) = zeps(1)*sqrt(zmu(1))/2
        zmat0(3,2) = zeps(2)*sqrt(zmu(2))/2
      
        zmat0(4,3) = -zmu(1)*sqrt(zeps(1))/2
        zmat0(4,4) = -zmu(2)*sqrt(zeps(2))/2
        call zinverse(4,zmat0,info,zmat0inv)

        call zmatmat(4,4,zmat0inv,4,zmat,zmat2)
c
c
c  compute singular values and eigenvalues of diagonally
c  preconditioned matrix
c


        call zeigs(4,zmat2,info,zeigl,zlam,zeigr)
        call zsvd(4,4,zmat2,zsingl,sing,zsingr)
        do i=1,4
          zlamall(i,nn) = zlam(i)
          singall(i,nn) = sing(i)
        enddo
      enddo



      

      return
      end






      subroutine test_solver_get_lamsing2(zeps,zepssq,zmu,zmusq,
     1   domeg,nnmax,zlamall,singall,errall)

      implicit real *8 (a-h,o-z)
      complex *16, allocatable :: fjvals(:,:),fhvals(:,:)
      complex *16, allocatable :: fjder(:,:),fhder(:,:)
      complex *16 zeps(2),zmu(2),zk(2),zepssq(2),zmusq(2)
      complex *16 zalpha(2),zbeta(2),zg1(2),zg2(2),zeta(2),za1(2),zd(2)
      complex *16 za2(2)
      complex *16 ima
      complex *16 zmat(4,4),zrhs(4),zsoln(4)
      complex *16 zmat2(4,4),zmat0(4,4),zmat0inv(4,4)
      complex *16 zeu(2),zex(2),zey(2)
      complex *16 zeu0(2),zex0(2),zey0(2)
      complex *16 zbu(2),zbx(2),zby(2)
      complex *16 zbu0(2),zbx0(2),zby0(2)
      complex *16 zsingl(4,4),zsingr(4,4)
      complex *16 zeigl(4,4),zeigr(4,4),zlam(4)
      complex *16 zlamall(4,nnmax)
      real *8 singall(4,nnmax),errall(nnmax)
      real *8 sing(4)

      data ima/(0.0d0,1.0d0)/
      call prini(6,13)
c
c
c  1 - interior
c  2 - exterior
c


      ifder = 1
      rscale = 1.0d0
      
      njh = nnmax + 10
      allocate(fjvals(0:njh,2),fjder(0:njh,2))
      allocate(fhvals(0:njh,2),fhder(0:njh,2))

      do i=1,2
        zk(i) = domeg*zepssq(i)*zmusq(i)

        call besseljs3d(njh,zk(i),rscale,fjvals(0,i),ifder,fjder(0,i))
        call h3dall(njh,zk(i),rscale,fhvals(0,i),ifder,fhder(0,i))
      enddo


      do nn=1,nnmax
        do i=1,2
          zalpha(i) = -ima*((nn+1)*fjvals(nn,i)*fhvals(nn-1,i)+
     1       nn*fjvals(nn+1,i)*fhvals(nn,i)-
     2       zk(i)*fjvals(nn+1,i)*fhvals(nn-1,i))
          zbeta(i) = fjvals(nn,i)*fhvals(nn,i)*zk(i)/ima
          za1(i) = ima*zk(i)**2*fjvals(nn,i)*fhder(nn,i)
          zg1(i) = ima*((fhvals(nn,i)+
     1       zk(i)*fhder(nn,i))*fjvals(nn,i)*zk(i)) 
          za2(i) = ima*zk(i)**2*fhvals(nn,i)*fjder(nn,i)
          zg2(i) = ima*((fjvals(nn,i)+
     1       zk(i)*fjder(nn,i))*fhvals(nn,i)*zk(i)) 
          zd(i) = -ima*zk(i)*fjvals(nn,i)*fhvals(nn,i)*
     1       sqrt(nn*(nn+1.0d0)) 
          zeta(i) =ima*sqrt(nn*(nn+1.0d0))*(fjvals(nn,i)*fhvals(nn-1,i)-
     1     fjvals(nn+1,i)*fhvals(nn,i))
        enddo
c
c  Now assemble the matrix
c
c  unknowns ordering r_{-} r_{+} q_{-} q_{+}
c  Equations ordering 
c   1. G_{0} \nabla_{s} (E_{+} - E^{-})_{tan} = G_{0} \nabla_{s} (E^{in})  
c   2. G_{0} \nabla_{s} (B_{+} - B^{-})_{tan} = G_{0} \nabla_{s} (B^{in})  
c   3. \eps_{+} E^{+}.n - \eps_{-} E^{-} \cdot n = \eps_{+} E^{in} \cdot n
c   4. \mu_{+} B^{+}.n - \mu_{-} B^{-} \cdot n = \mu_{+} B^{in} \cdot n
c
c
c   Note sign flipped on the right hand side in the above equation
c   by design since we are doing an analytical test and setting
c   E^{in} = E^{+} for some known E^{+} and similarly for B^{+}
c
        zmat(1,1) = -ima*zk(1)*zmusq(1)*ima*zk(1)*zalpha(1) - 
     1     zmusq(1)*zd(1)*sqrt(nn*(nn+1.0d0)) 
        zmat(1,2) = ima*zk(2)*zmusq(2)*ima*zk(2)*zalpha(2) + 
     1     zmusq(2)*zd(2)*sqrt(nn*(nn+1.0d0))
        zmat(1,3) = -zg1(2)*ima*domeg*zmu(1)*zepssq(1)
        zmat(1,4) = -zg2(1)*ima*domeg*zmu(2)*zepssq(2)
        zmat(1,1:4) = zmat(1,1:4)/(2*nn+1.0d0)
c
c
c
c
        zmat(2,1) = -zg1(2)*ima*domeg*zeps(1)*zmusq(1)
        zmat(2,2) = -zg2(1)*ima*domeg*zeps(2)*zmusq(2)
        zmat(2,3) = ima*zk(1)*zepssq(1)*ima*zk(1)*zalpha(1) + 
     1     zepssq(1)*zd(1)*sqrt(nn*(nn+1.0d0))
        zmat(2,4) = -ima*zk(2)*zepssq(2)*ima*zk(2)*zalpha(2) - 
     1     zepssq(2)*zd(2)*sqrt(nn*(nn+1.0d0))
        zmat(2,1:4) = zmat(2,1:4)/(2*nn+1.0d0)
c
c
c
c
        zmat(3,1) = zeps(1)*zmusq(1)*
     1      ((ima*zk(1))**2/sqrt(nn*(nn+1.0d0))*zeta(1)+za2(1))
        zmat(3,2) = -zeps(2)*zmusq(2)*
     1     ((ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*zeta(2)+za1(2))
        zmat(3,3) = zeps(2)*ima*domeg*zmu(1)*zepssq(1)*zbeta(2)
        zmat(3,4) = zeps(1)*ima*domeg*zmu(2)*zepssq(2)*zbeta(1)
c
c
c
c
        zmat(4,1) = zmu(2)*ima*domeg*zeps(1)*zmusq(1)*zbeta(2)
        zmat(4,2) = -zmu(1)*ima*domeg*zeps(2)*zmusq(2)*zbeta(1)
        zmat(4,3) = -zmu(1)*zepssq(1)*
     1      ((ima*zk(1))**2/sqrt(nn*(nn+1.0d0))*zeta(1)+za2(1))
        zmat(4,4) = zmu(2)*zepssq(2)*
     1     ((ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*zeta(2)+za1(2))


        
        zrhs(1) = 0
        zrhs(2) =-sqrt(nn*(nn+1.0d0))*zg1(2)/(2*nn+1.0d0)*zepssq(2) 
        zrhs(3) = 0
        zrhs(4) = sqrt(nn*(nn+1.0d0))*zbeta(2)*zmu(2)*zepssq(2)

        zsoln(1) = 0
        zsoln(2) = 0
        zsoln(3) = 0
        zsoln(4) = 0

        call zgausselim(4,zmat,zrhs,info,zsoln,dcond)
c
c
c  Test exterior electric field
c


        zeu0(2) = 0
        zex0(2) = -zbeta(2)*ima*zk(2)*zmusq(2)
        zey0(2) = 0

        zeu(2) = 0
        zex(2) = -ima*zk(2)*zmusq(2)*ima*domeg*zeps(1)*zmusq(1)/
     1     zepssq(2)*zsoln(1)/sqrt(nn*(nn+1.0d0))*zbeta(2) + 
     2     zmusq(2)*ima*zk(2)*zsoln(4)/sqrt(nn*(nn+1.0d0))*zg2(2)
        zey(2) = 0

        errep = abs(zex(2)-zex0(2))/abs(zex(2))

c
c
c  test exterior magnetic field
c

        zbx0(2) = 0
        zby0(2) = zepssq(2)*sqrt(nn*(nn+1.0d0))*zbeta(2)
        zbu0(2) = -zepssq(2)*zg1(2)

        zbx(2) = 0
        zby(2) = zepssq(2)*(ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*
     1     zsoln(4)*zeta(2) + ima*domeg*zeps(1)*zmusq(1)*
     2     zsoln(1)*zbeta(2) + za1(2)*zepssq(2)*zsoln(4)
        zbu(2) = -zepssq(2)*(ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*
     1     zsoln(4)*zalpha(2) - ima*domeg*zeps(1)*zmusq(1)*
     2     zsoln(1)*zg1(2)/sqrt(nn*(nn+1.0d0)) - 
     3     zsoln(4)*zepssq(2)*zd(2)


        errbp = abs(zby0(2)-zby(2)) + abs(zbu0(2)-zbu(2))
        rbp = abs(zby0(2))+abs(zbu0(2))
        errbp = errbp/rbp
        errall(nn) = max(errbp,errep)
c
c  diagonal precondition
c
c

        do i=1,4
          do j=1,4
            zmat0(j,i) = 0
          enddo
        enddo
        zmat0(1,1) = zmusq(1)/4
        zmat0(1,2) = -zmusq(2)/4

        zmat0(2,3) = -zepssq(1)/4
        zmat0(2,4) = zepssq(2)/4
      
        zmat0(3,1) = zeps(1)*zmusq(1)/2
        zmat0(3,2) = zeps(2)*zmusq(2)/2
      
        zmat0(4,3) = -zmu(1)*zepssq(1)/2
        zmat0(4,4) = -zmu(2)*zepssq(2)/2
        call zinverse(4,zmat0,info,zmat0inv)

        call zmatmat(4,4,zmat0inv,4,zmat,zmat2)
c
c
c  compute singular values and eigenvalues of diagonally
c  preconditioned matrix
c


        call zeigs(4,zmat2,info,zeigl,zlam,zeigr)
        call zsvd(4,4,zmat2,zsingl,sing,zsingr)
        do i=1,4
          zlamall(i,nn) = zlam(i)
          singall(i,nn) = sing(i)
        enddo
      enddo



      

      return
      end
c
c
c
c
c
c
      
      subroutine test_isin(rkp,rkm,isin)
      implicit real *8 (a-h,o-z)
      isin = 0
      done = 1
      pi = atan(done)*4
      if(abs(rkp-rkm).le.pi/2) then
        isin = 1
      endif
      if(abs(rkm-pi/2).le.1.0d-16.and.abs(rkp).le.1.0d-16) isin = 0
      if(abs(rkm-pi/2).le.1.0d-16.and.abs(rkp-pi).le.1.0d-16) isin = 0
      
      if(abs(rkp).le.1.0d-16.and.rkm.gt.pi/2) isin = 1
      if(abs(rkp-pi).le.1.0d-16.and.rkm.lt.pi/2) isin = 1

      return
      end
c
c
c
c
c


      subroutine comp_thet(rkp,rkm,thet)
      implicit real *8 (a-h,o-z)
      done = 1
      pi = atan(done)*4

      alpha = rkm - rkp
      if(alpha.ge.0.and.alpha.le.pi/2.and.rkp.gt.0.and.
     1    rkp.le.pi/2) thet = rkp
      if(alpha.ge.-pi/2.and.alpha.lt.0.and.rkm.gt.0
     1        .and.rkm.le.pi/2) thet = rkp+ 2*alpha
      if(rkm.gt.pi/2.and.rkm.lt.pi.and.rkp.gt.pi/2.and.rkp.lt.pi)
     1  thet = pi-rkp    

      if(abs(rkp).le.1.0d-16.and.rkm.lt.pi/2) thet=0
      if(abs(rkp).le.1.0d-16.and.rkm.gt.pi/2) thet=pi
      if(abs(rkm).le.1.0d-16) thet = -rkp
      if(abs(rkm-pi).le.1.0d-16) thet = pi-rkp
      if(abs(rkp-pi).le.1.0d-16.and.rkm.lt.pi/2) thet=pi
      if(abs(rkp-pi).le.1.0d-16.and.rkm.gt.pi/2) thet=0


      return
      end
c
c
c
c
c
      subroutine get_det(zeps,zepssq,zmu,zmusq,
     1   domeg,n,zdet,zdet2,erra)

      implicit real *8 (a-h,o-z)
      complex *16, allocatable :: fjvals(:,:),fhvals(:,:)
      complex *16, allocatable :: fjder(:,:),fhder(:,:)
      complex *16 zeps(2),zmu(2),zk(2),zepssq(2),zmusq(2)
      complex *16 zalpha(2),zbeta(2),zg1(2),zg2(2),zeta(2),za1(2),zd(2)
      complex *16 za2(2)
      complex *16 ima
      complex *16 zmat(4,4),zrhs(4),zsoln(4)
      complex *16 zmat2(4,4),zmat0(4,4),zmat0inv(4,4)
      complex *16 zeu(2),zex(2),zey(2)
      complex *16 zeu0(2),zex0(2),zey0(2)
      complex *16 zbu(2),zbx(2),zby(2)
      complex *16 zbu0(2),zbx0(2),zby0(2)
      complex *16 zsingl(4,4),zsingr(4,4)
      complex *16 zeigl(4,4),zeigr(4,4),zlam(4)
      complex *16 zdet,zdet2
      real *8 erra

      data ima/(0.0d0,1.0d0)/
      call prini(6,13)
c
c
c  1 - interior
c  2 - exterior
c
      ifder = 1
      rscale = 1.0d0
      
      njh = n + 10
      allocate(fjvals(0:njh,2),fjder(0:njh,2))
      allocate(fhvals(0:njh,2),fhder(0:njh,2))

      do i=1,2
        zk(i) = domeg*zepssq(i)*zmusq(i)

        call besseljs3d(njh,zk(i),rscale,fjvals(0,i),ifder,fjder(0,i))
        call h3dall(njh,zk(i),rscale,fhvals(0,i),ifder,fhder(0,i))
      enddo

      nn = n
      do i=1,2
        zalpha(i) = -ima*((nn+1)*fjvals(nn,i)*fhvals(nn-1,i)+
     1     nn*fjvals(nn+1,i)*fhvals(nn,i)-
     2     zk(i)*fjvals(nn+1,i)*fhvals(nn-1,i))
        zbeta(i) = fjvals(nn,i)*fhvals(nn,i)*zk(i)/ima
        za1(i) = ima*zk(i)**2*fjvals(nn,i)*fhder(nn,i)
        zg1(i) = ima*((fhvals(nn,i)+
     1       zk(i)*fhder(nn,i))*fjvals(nn,i)*zk(i)) 
        za2(i) = ima*zk(i)**2*fhvals(nn,i)*fjder(nn,i)
        zg2(i) = ima*((fjvals(nn,i)+
     1     zk(i)*fjder(nn,i))*fhvals(nn,i)*zk(i)) 
        zd(i) = -ima*zk(i)*fjvals(nn,i)*fhvals(nn,i)*
     1       sqrt(nn*(nn+1.0d0)) 
        zeta(i) =ima*sqrt(nn*(nn+1.0d0))*(fjvals(nn,i)*fhvals(nn-1,i)-
     1     fjvals(nn+1,i)*fhvals(nn,i))
      enddo
c
c  Now assemble the matrix
c
c  unknowns ordering r_{-} r_{+} q_{-} q_{+}
c  Equations ordering 
c   1. G_{0} \nabla_{s} (E_{+} - E^{-})_{tan} = G_{0} \nabla_{s} (E^{in})  
c   2. G_{0} \nabla_{s} (B_{+} - B^{-})_{tan} = G_{0} \nabla_{s} (B^{in})  
c   3. \eps_{+} E^{+}.n - \eps_{-} E^{-} \cdot n = \eps_{+} E^{in} \cdot n
c   4. \mu_{+} B^{+}.n - \mu_{-} B^{-} \cdot n = \mu_{+} B^{in} \cdot n
c
c
c   Note sign flipped on the right hand side in the above equation
c   by design since we are doing an analytical test and setting
c   E^{in} = E^{+} for some known E^{+} and similarly for B^{+}
c
      zmat(1,1) = -ima*zk(1)*zmusq(1)*ima*zk(1)*zalpha(1) - 
     1   zmusq(1)*zd(1)*sqrt(nn*(nn+1.0d0)) 
      zmat(1,2) = ima*zk(2)*zmusq(2)*ima*zk(2)*zalpha(2) + 
     1   zmusq(2)*zd(2)*sqrt(nn*(nn+1.0d0))
      zmat(1,3) = -zg1(2)*ima*domeg*zmu(1)*zepssq(1)
      zmat(1,4) = -zg2(1)*ima*domeg*zmu(2)*zepssq(2)
      zmat(1,1:4) = zmat(1,1:4)/(2*nn+1.0d0)
c
c
c
c
      zmat(2,1) = -zg1(2)*ima*domeg*zeps(1)*zmusq(1)
      zmat(2,2) = -zg2(1)*ima*domeg*zeps(2)*zmusq(2)
      zmat(2,3) = ima*zk(1)*zepssq(1)*ima*zk(1)*zalpha(1) + 
     1   zepssq(1)*zd(1)*sqrt(nn*(nn+1.0d0))
      zmat(2,4) = -ima*zk(2)*zepssq(2)*ima*zk(2)*zalpha(2) - 
     1     zepssq(2)*zd(2)*sqrt(nn*(nn+1.0d0))
      zmat(2,1:4) = zmat(2,1:4)/(2*nn+1.0d0)
c
c
c
c
      zmat(3,1) = zeps(1)*zmusq(1)*
     1   ((ima*zk(1))**2/sqrt(nn*(nn+1.0d0))*zeta(1)+za2(1))
      zmat(3,2) = -zeps(2)*zmusq(2)*
     1   ((ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*zeta(2)+za1(2))
      zmat(3,3) = zeps(2)*ima*domeg*zmu(1)*zepssq(1)*zbeta(2)
      zmat(3,4) = zeps(1)*ima*domeg*zmu(2)*zepssq(2)*zbeta(1)
c
c
c
c
      zmat(4,1) = zmu(2)*ima*domeg*zeps(1)*zmusq(1)*zbeta(2)
      zmat(4,2) = -zmu(1)*ima*domeg*zeps(2)*zmusq(2)*zbeta(1)
      zmat(4,3) = -zmu(1)*zepssq(1)*
     1    ((ima*zk(1))**2/sqrt(nn*(nn+1.0d0))*zeta(1)+za2(1))
      zmat(4,4) = zmu(2)*zepssq(2)*
     1   ((ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*zeta(2)+za1(2))


        
      zrhs(1) = 0
      zrhs(2) =-sqrt(nn*(nn+1.0d0))*zg1(2)/(2*nn+1.0d0)*zepssq(2) 
      zrhs(3) = 0
      zrhs(4) = sqrt(nn*(nn+1.0d0))*zbeta(2)*zmu(2)*zepssq(2)

      zsoln(1) = 0
      zsoln(2) = 0
      zsoln(3) = 0
      zsoln(4) = 0

      call zgausselim(4,zmat,zrhs,info,zsoln,dcond)
c
c
c  Test exterior electric field
c


      zeu0(2) = 0
      zex0(2) = -zbeta(2)*ima*zk(2)*zmusq(2)
      zey0(2) = 0

      zeu(2) = 0
      zex(2) = -ima*zk(2)*zmusq(2)*ima*domeg*zeps(1)*zmusq(1)/
     1   zepssq(2)*zsoln(1)/sqrt(nn*(nn+1.0d0))*zbeta(2) + 
     2     zmusq(2)*ima*zk(2)*zsoln(4)/sqrt(nn*(nn+1.0d0))*zg2(2)
      zey(2) = 0

      errep = abs(zex(2)-zex0(2))/abs(zex(2))

c
c
c  test exterior magnetic field
c

      zbx0(2) = 0
      zby0(2) = zepssq(2)*sqrt(nn*(nn+1.0d0))*zbeta(2)
      zbu0(2) = -zepssq(2)*zg1(2)

      zbx(2) = 0
      zby(2) = zepssq(2)*(ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*
     1   zsoln(4)*zeta(2) + ima*domeg*zeps(1)*zmusq(1)*
     2   zsoln(1)*zbeta(2) + za1(2)*zepssq(2)*zsoln(4)
      zbu(2) = -zepssq(2)*(ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*
     1   zsoln(4)*zalpha(2) - ima*domeg*zeps(1)*zmusq(1)*
     2   zsoln(1)*zg1(2)/sqrt(nn*(nn+1.0d0)) - 
     3   zsoln(4)*zepssq(2)*zd(2)


      errbp = abs(zby0(2)-zby(2)) + abs(zbu0(2)-zbu(2))
      rbp = abs(zby0(2))+abs(zbu0(2))
      errbp = errbp/rbp
      erra = max(errbp,errep)
c
c  diagonal precondition
c
c
      do i=1,4
        do j=1,4
          zmat0(j,i) = 0
        enddo
      enddo
      zmat0(1,1) = zmusq(1)/4
      zmat0(1,2) = -zmusq(2)/4

      zmat0(2,3) = -zepssq(1)/4
      zmat0(2,4) = zepssq(2)/4
      
      zmat0(3,1) = zeps(1)*zmusq(1)/2
      zmat0(3,2) = zeps(2)*zmusq(2)/2
      
      zmat0(4,3) = -zmu(1)*zepssq(1)/2
      zmat0(4,4) = -zmu(2)*zepssq(2)/2
      call zinverse(4,zmat0,info,zmat0inv)

      call zmatmat(4,4,zmat0inv,4,zmat,zmat2)
c
c
c  compute singular values and eigenvalues of diagonally
c  preconditioned matrix
c


      call zeigs(4,zmat2,info,zeigl,zlam,zeigr)
      zdet2 = 1
      do i=1,4
        zdet2 = zdet2*zlam(i)
      enddo

      call zeigs(4,zmat,info,zeigl,zlam,zeigr)
      zdet = 1
      do i=1,4
        zdet = zdet*zlam(i)
      enddo



      

      return
      end
