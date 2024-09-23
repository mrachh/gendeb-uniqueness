      implicit real *8 (a-h,o-z)
      complex *16, allocatable :: fjvals(:,:),fhvals(:,:)
      complex *16, allocatable :: fjder(:,:),fhder(:,:)
      complex *16 zeps(2),zmu(2),zk(2)
      complex *16 zalpha(2),zbeta(2),zg1(2),zg2(2),zeta(2),za1(2),zd(2)
      complex *16 za2(2)
      complex *16 ima
      complex *16 zmat(4,4),zrhs(4),zsoln(4)
      complex *16 zeu(2),zex(2),zey(2)
      complex *16 zeu0(2),zex0(2),zey0(2)
      complex *16 zbu(2),zbx(2),zby(2)
      complex *16 zbu0(2),zbx0(2),zby0(2)

      data ima/(0.0d0,1.0d0)/
      call prini(6,13)
c
c
c  1 - interior
c  2 - exterior
c


      ifder = 1
      rscale = 1.0d0
      
      zeps(1) = 1.0d0
      zeps(2) = 1.2d0
      
      zmu(1) = 1.0d0
      zmu(2) = 1.2d0
      nn = 3
      njh = nn + 10
      allocate(fjvals(0:njh,2),fjder(0:njh,2))
      allocate(fhvals(0:njh,2),fhder(0:njh,2))

      domeg = 1.0d0
      do i=1,2
        zk(i) = domeg*sqrt(zeps(i))*sqrt(zmu(i))

        call besseljs3d(njh,zk(i),rscale,fjvals(0,i),ifder,fjder(0,i))
        call h3dall(njh,zk(i),rscale,fhvals(0,i),ifder,fhder(0,i))
        zalpha(i) = -ima*((nn+1)*fjvals(nn,i)*fhvals(nn-1,i)+
     1     nn*fjvals(nn+1,i)*fhvals(nn,i)-
     2     zk(i)*fjvals(nn+1,i)*fhvals(nn-1,i))
        zbeta(i) = fjvals(nn,i)*fhvals(nn,i)*zk(i)/ima
        za1(i) = ima*zk(i)**2*fjvals(nn,i)*fhder(nn,i)
        zg1(i) = ima*((fhvals(nn,i)+
     1     zk(i)*fhder(nn,i))*fjvals(nn,i)*zk(i)) 
        za2(i) = ima*zk(i)**2*fhvals(nn,i)*fjder(nn,i)
        zg2(i) = ima*((fjvals(nn,i)+
     1     zk(i)*fjder(nn,i))*fhvals(nn,i)*zk(i)) 
        zd(i) = -ima*zk(i)*fjvals(nn,i)*fhvals(nn,i)*
     1     sqrt(nn*(nn+1.0d0)) 
        zeta(i) = ima*sqrt(nn*(nn+1.0d0))*(fjvals(nn,i)*fhvals(nn-1,i)-
     1   fjvals(nn+1,i)*fhvals(nn,i))
      enddo
      call prin2('zbeta=*',zbeta,4)
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
     1   sqrt(zmu(1))*zd(1)*sqrt(nn*(nn+1.0d0)) 
      zmat(1,2) = ima*zk(2)*sqrt(zmu(2))*ima*zk(2)*zalpha(2) + 
     1   sqrt(zmu(2))*zd(2)*sqrt(nn*(nn+1.0d0))
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
     1   sqrt(zeps(1))*zd(1)*sqrt(nn*(nn+1.0d0))
      zmat(2,4) = -ima*zk(2)*sqrt(zeps(2))*ima*zk(2)*zalpha(2) - 
     1   sqrt(zeps(2))*zd(2)*sqrt(nn*(nn+1.0d0))
      zmat(2,1:4) = zmat(2,1:4)/(2*nn+1.0d0)
c
c
c
c
      zmat(3,1) = zeps(1)*sqrt(zmu(1))*
     1    ((ima*zk(1))**2/sqrt(nn*(nn+1.0d0))*zeta(1)+za2(1))
      zmat(3,2) = -zeps(2)*sqrt(zmu(2))*
     1   ((ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*zeta(2)+za1(2))
      zmat(3,3) = zeps(2)*ima*domeg*zmu(1)*sqrt(zeps(1))*zbeta(2)
      zmat(3,4) = zeps(1)*ima*domeg*zmu(2)*sqrt(zeps(2))*zbeta(1)
c
c
c
c
      zmat(4,1) = zmu(2)*ima*domeg*zeps(1)*sqrt(zmu(1))*zbeta(2)
      zmat(4,2) = -zmu(1)*ima*domeg*zeps(2)*sqrt(zmu(2))*zbeta(1)
      zmat(4,3) = -zmu(1)*sqrt(zeps(1))*
     1    ((ima*zk(1))**2/sqrt(nn*(nn+1.0d0))*zeta(1)+za2(1))
      zmat(4,4) = zmu(2)*sqrt(zeps(2))*
     1   ((ima*zk(2))**2/sqrt(nn*(nn+1.0d0))*zeta(2)+za1(2))


        
      zrhs(1) = 0
      zrhs(2) = -sqrt(nn*(nn+1.0d0))*zg1(2)/(2*nn+1.0d0)*sqrt(zeps(2)) 
      zrhs(3) = 0
      zrhs(4) = sqrt(nn*(nn+1.0d0))*zbeta(2)*sqrt(zeps(2))*zmu(2)

      call zgausselim(4,zmat,zrhs,info,zsoln,dcond)
      call prin2('zmat=*',zmat,32)
      call prin2('zrhs=*',zrhs,8)
      call prin2('zsoln=*',zsoln,8)
      call prin2('dcond=*',dcond,1)
c
c
c  Test exterior electric field
c


      zeu0(2) = 0
      zex0(2) = -zbeta(2)*ima*zk(2)*sqrt(zmu(2))
      zey0(2) = 0

      zeu(2) = 0
      zex(2) = -ima*zk(2)*sqrt(zmu(2))*ima*domeg*zeps(1)*sqrt(zmu(1))/
     1   sqrt(zeps(2))*zsoln(1)/sqrt(nn*(nn+1.0d0))*zbeta(2) + 
     2   sqrt(zmu(2))*ima*zk(2)*zsoln(4)/sqrt(nn*(nn+1.0d0))*zg2(2)
      zey(2) = 0

      errep = abs(zex(2)-zex0(2))/abs(zex(2))
      call prin2('error in e+=*',errep,1)

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

      call prin2('zby=*',zby(2),2)
      call prin2('zby0=*',zby0(2),2)
      call prin2('zbu=*',zbu(2),2)
      call prin2('zbu0=*',zbu0(2),2)

      errbp = abs(zby0(2)-zby(2)) + abs(zbu0(2)-zbu(2))
      rbp = abs(zby0(2))+abs(zbu0(2))
      errbp = errbp/rbp
      call prin2('errbp=*',errbp,1)




      

      stop
      end
