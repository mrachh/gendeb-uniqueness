      implicit real *8 (a-h,o-z)
      complex *16 zeps(2),zmu(2),zepssq(2),zmusq(2)
      real *8 domeg
      complex *16, allocatable :: zlamall(:,:)
      real *8, allocatable :: singall(:,:),errall(:)
      real *8, allocatable :: argepsm(:,:),argmum(:,:)
      integer, allocatable :: isin(:,:)
      real *8, allocatable :: rcond(:,:),rerr(:,:)
      complex *16 ima,z,z1,z2,z3,z4
      data ima/(0.0d0,1.0d0)/
      call prini(6,13)

      done = 1
      pi = atan(done)*4

      

      nnmax = 30
      allocate(zlamall(4,nnmax),singall(4,nnmax),errall(nnmax))


      nlat = 151
      allocate(argepsm(nlat,nlat),argmum(nlat,nlat),isin(nlat,nlat))
      allocate(rcond(nlat,nlat),rerr(nlat,nlat))

      rkpmin = 0 
      rkpmax = pi
      
      rkmmin = 0
      rkmmax = pi
      domeg = 1.0d0
 1111 format(2(2x,e11.5),2x,i2,2(2x,e11.5))      
      do i=1,nlat
        call prinf('i=*',i,1)
        do j=1,nlat
          
          rkp = rkpmax + (rkpmin-rkpmax)*(j-1.0d0)/(nlat-1.0d0)
          rkm = rkmmin + (rkmmax-rkmmin)*(i-1.0d0)/(nlat-1.0d0)


          call test_isin(rkp,rkm,isin(j,i))
          zeps(1) = exp(ima*(2*rkm-rkp))
          zeps(2) = 2*exp(ima*rkp)
          zmu(1) = exp(ima*rkp)
          zmu(2) = exp(ima*rkp)


          if(isin(j,i).eq.1) then
            call comp_thet(rkp,rkm,thet)

            if(imag(z1).lt.-1.0d-12) then
              call prinf('poor choice of thet1=*',i,0)
              call prin2('zeps(1)=*',zeps(1),2)
              call prin2('z=*',z1,2)
              call prin2('argkp=*',rkp,1)
              call prin2('argkm=*',rkm,1)
            else if(imag(z1).lt.0.and.imag(z1).ge.-1.0d-12) then
              z1 = real(z1)
            endif
            
            z2 = exp(-ima*thet)*zeps(2)
            if(imag(z2).lt.-1.0d-12) then
              call prinf('poor choice of thet2=*',i,0)
              call prin2('zeps(2)=*',zeps(2),2)
              call prin2('z=*',z2,2)
            else if(imag(z2).lt.0.and.imag(z2).ge.-1.0d-12) then
              z2 = real(z2)
            endif
            
            z3 = exp(ima*thet)*zmu(1)
            if(imag(z3).lt.-1.0d-12) then
              call prinf('poor choice of thet3=*',i,0)
              call prin2('zmu(1)=*',zmu(1),2)
              call prin2('z=*',z1,2)
            else if(imag(z3).lt.0.and.imag(z3).ge.-1.0d-12) then
              z3 = real(z3)
            endif
            
            z4 = exp(ima*thet)*zmu(2)
            if(imag(z1).lt.-1.0d-12) then
              call prinf('poor choice of thet4=*',i,0)
              call prin2('zmu(2)=*',zmu(2),2)
              call prin2('z=*',z4,2)
            else if(imag(z4).lt.0.and.imag(z4).ge.-1.0d-12) then
              z4 = real(z4)
            endif

            rphi = atan2(imag(z1),real(z1))
            zepssq(1) = exp(ima*thet/2)*exp(ima*rphi/2)
            
            rphi = atan2(imag(z2),real(z2))
            zepssq(2) = sqrt(2.0d0)*exp(ima*thet/2)*exp(ima*rphi/2)
            
            rphi = atan2(imag(z3),real(z3))
            zmusq(1) = exp(-ima*thet/2)*exp(ima*rphi/2)
            
            rphi = atan2(imag(z4),real(z4))
            zmusq(2) = exp(-ima*thet/2)*exp(ima*rphi/2)

            z = zepssq(1)*zmusq(1)
            if(imag(z).lt.-1.0d-12) then
              call prinf('poor choice of thet5=*',i,0)
              call prin2('z=*',z,2)
              call prin2('rkp=*',rkp,1)
              call prin2('rkm=*',rkm,1)
            endif
            z = zepssq(2)*zmusq(2)
            if(imag(z).lt.-1.0d-12) then
              call prinf('poor choice of thet6=*',i,0)
              call prin2('z=*',z,2)
              call prin2('rkp=*',rkp,1)
              call prin2('rkm=*',rkm,1)
            endif


          else
            zepssq(1) = sqrt(zeps(1))
            zepssq(2) = sqrt(zeps(2))
            zmusq(1) = sqrt(zmu(1))
            zmusq(2) = sqrt(zmu(2))
          endif
          do l=1,nnmax
            zlamall(1,l)=0
            zlamall(2,l)=0
            zlamall(3,l)=0
            zlamall(4,l)=0

            singall(1,l) = 0
            singall(2,l) = 0
            singall(3,l) = 0
            singall(4,l) = 0
            errall(l) = 0

          enddo
          call test_solver_get_lamsing2(zeps,zepssq,zmu,zmusq,domeg,
     1      nnmax,zlamall,singall,errall)
          rmax = singall(1,1)
          rmin = singall(1,1)
          err1 = errall(1)

          do l=1,nnmax
            do m=1,4
              if(singall(m,l).gt.rmax) rmax = singall(m,l)
              if(singall(m,l).lt.rmin) rmin = singall(m,l)
            enddo
            if(errall(l).gt.err1) err1 = errall(l)
          enddo
          rcond(j,i) = rmax/rmin
          write(84,1111) rkm,rkp,isin(j,i),
     1       rcond(j,i),err1

          if(abs(rkp).le.1.0d-16.and.rkm.gt.pi/2) then
            write(85,*) rkm,rkp,isin(j,i),rcond(j,i),err1
          endif

          if(abs(rkp-pi).le.1.0d-16.and.rkm.lt.pi/2) then
            write(86,*) rkm,rkp,isin(j,i),rcond(j,i),err1
          endif
        enddo
      enddo
      


      stop
      end





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




