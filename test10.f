      implicit real *8 (a-h,o-z) 

      real *8, allocatable :: targs(:,:)
      complex *16, allocatable :: ynm_targ(:)
      complex *16, allocatable :: vynm_targ(:,:)
      complex *16, allocatable :: phinm_targ(:,:)
      complex *16, allocatable :: psinm_targ(:,:)

      complex *16 zalpha,zbeta,zgamma,zdelta,zeta,zteta,zk,ztetap
      complex *16 ztetam
      complex *16 fjvalst(0:100),fhvalst(0:100)
      complex *16 fjdert(0:100),fhdert(0:100)
      complex *16 z1,z2,z3,z4
      complex *16 zvec1(3),zvec2(3),zvec3(3)
      complex *16 dfdx, dfdy, dfdz, f(100), fvec(3,100)
      complex *16 dfvec(3,3)
      real *8 dvec1(3),dvec2(3),dvec3(3)

      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname


      logical isout0,isout1

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4
c
c  set the value of dzk
c
      zk = 0.19d0


      njh = 5
      ifder = 1
      rscale = 1.0d0
      ntargout = 1


c      
c     test r \cross \grad = \phi_nm upto scaling 
c
      ntarg = 7
      allocate(targs(3,ntarg))
      
      r0 = 1.3d0
      thet0 = hkrand(0)*pi
      phi0 = hkrand(0)*2*pi

      targs(1,1) = r0*sin(thet0)*cos(phi0)
      targs(2,1) = r0*sin(thet0)*sin(phi0)
      targs(3,1) = r0*cos(thet0)

      do i=2,7
        targs(1,i) = targs(1,1)
        targs(2,i) = targs(2,1)
        targs(3,i) = targs(3,1)

      enddo
      h = 1.0d-3
      targs(1,2) = targs(1,1) + h
      targs(1,3) = targs(1,1) - h

      targs(2,4) = targs(2,1) + h
      targs(2,5) = targs(2,1) - h
      
      targs(3,6) = targs(3,1) + h
      targs(3,7) = targs(3,7) - h

      allocate(ynm_targ(ntarg))
      allocate(vynm_targ(3,ntarg),phinm_targ(3,ntarg))
      allocate(psinm_targ(3,ntarg))

      nn = 3
      mm = 1


      call l3getsph_scalar_vec(mm,nn,3,ntarg,targs,ynm_targ,
     1   vynm_targ,psinm_targ,phinm_targ)

      
      dfdx = (ynm_targ(2) - ynm_targ(3))/2/h
      dfdy = (ynm_targ(4) - ynm_targ(5))/2/h
      dfdz = (ynm_targ(6) - ynm_targ(7))/2/h

      erra = abs(dfdx*r0 - psinm_targ(1,1))
      erra = erra + abs(dfdy*r0 - psinm_targ(2,1))
      erra = erra + abs(dfdz*r0 - psinm_targ(3,1))

      ra = abs(psinm_targ(1,1)) + abs(psinm_targ(2,1))
      ra = ra + abs(psinm_targ(3,1))
      
      erra =  erra/ra
      call prin2('error in r \grad =*', erra, 1)

      
      x = targs(1,1)
      y = targs(2,1)
      z = targs(3,1)
      
      zvec1(1) = (y*dfdz - z*dfdy)
      zvec1(2) = (z*dfdx - x*dfdz)
      zvec1(3) = (x*dfdy - y*dfdx)

      erra = abs(zvec1(1) - phinm_targ(1,1))
      erra = erra + abs(zvec1(2) - phinm_targ(2,1))
      erra = erra + abs(zvec1(3) - phinm_targ(3,1))

      ra = abs(phinm_targ(1,1)) + abs(phinm_targ(2,1))
      ra = ra + abs(phinm_targ(3,1))
      
      call prin2('error in r \times \curl=*',erra/ra,1)

c
c
c
c  now test Mnm= \curl r j_n(kr) Y_{nm} is what we expect it to be
c
      zk = 0.19d0
      ifder = 0
      do i=1,7
        x = targs(1,i)
        y = targs(2,i)
        z = targs(3,i)
        r = sqrt(x**2 + y**2 + z**2)
        z1 = zk*r
        call besseljs3d(nn, z1, rscale, fjvalst, ifder, fjdert)
        f(i) = fjvalst(nn)
        fvec(1,i) = x*f(i)*ynm_targ(i)
        fvec(2,i) = y*f(i)*ynm_targ(i)
        fvec(3,i) = z*f(i)*ynm_targ(i)
      enddo

      do i=1,3 
        dfdx = (fvec(i,2) - fvec(i,3))/2/h
        dfdy = (fvec(i,4) - fvec(i,5))/2/h
        dfdz = (fvec(i,6) - fvec(i,7))/2/h

        dfvec(1,i) = dfdx
        dfvec(2,i) = dfdy
        dfvec(3,i) = dfdz
      enddo

      zvec1(1) = dfvec(2,3) - dfvec(3,2)
      zvec1(2) = dfvec(3,1) - dfvec(1,3)
      zvec1(3) = dfvec(1,2) - dfvec(2,1)

      zvec2(1) = -f(1)*phinm_targ(1,1)
      zvec2(2) = -f(1)*phinm_targ(2,1)
      zvec2(3) = -f(1)*phinm_targ(3,1)

      erra = abs(zvec1(1) - zvec2(1)) +
     1       abs(zvec1(2) - zvec2(2)) +
     2       abs(zvec1(3) - zvec2(3)) 
      ra = abs(zvec1(1)) + abs(zvec1(2)) + abs(zvec1(3))
      
      erra = erra/ra
      call prin2('error in Mnm=*',erra,1)

c 
c  Now test Nnm = \nabla \times (-jn(kr) Phi_{nm})/k 
c 
c 
      do i=1,7
        fvec(1:3,i) = -f(i)*phinm_targ(1:3,i)
      enddo

      do i=1,3 
        dfdx = (fvec(i,2) - fvec(i,3))/2/h
        dfdy = (fvec(i,4) - fvec(i,5))/2/h
        dfdz = (fvec(i,6) - fvec(i,7))/2/h

        dfvec(1,i) = dfdx
        dfvec(2,i) = dfdy
        dfvec(3,i) = dfdz
      enddo

      zvec1(1) = (dfvec(2,3) - dfvec(3,2))/zk
      zvec1(2) = (dfvec(3,1) - dfvec(1,3))/zk
      zvec1(3) = (dfvec(1,2) - dfvec(2,1))/zk

      x = targs(1,1)
      y = targs(2,1)
      z = targs(3,1)
      r = sqrt(x**2 + y**2 + z**2)
      z1 = zk*r
      ifder = 1
      call besseljs3d(nn, z1, rscale, fjvalst, ifder, fjdert)
      z1 = nn*(nn + 1.0d0)*f(1)/(zk*r0)
      z2 = fjdert(nn) + f(1)/(zk*r0)
      zvec2(1:3) = z1*vynm_targ(1:3,1) + z2*psinm_targ(1:3,1) 

      erra = abs(zvec1(1) - zvec2(1)) +
     1       abs(zvec1(2) - zvec2(2)) +
     2       abs(zvec1(3) - zvec2(3)) 
      ra = abs(zvec1(1)) + abs(zvec1(2)) + abs(zvec1(3))

      call prin2('zvec1=*',zvec1,6)
      call prin2('zvec2=*',zvec2,6)
      
      erra = erra/ra
      call prin2('error in Nnm=*',erra,1)

      stop
      end



c
c
c
c
c
      subroutine l3getsph_scalar_vec(mm,nn,ndx,npts,xyzs,ynms,
     1    vynm,psinm,phinm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ima
      complex *16 ynms(npts)
      complex *16 vynm(3,npts),psinm(3,npts),phinm(3,npts)
      real *8 vtmp(3)
      real *8, allocatable :: wlege(:),ynm(:,:),ynmd(:,:)
      complex *16 zr,zt,zp
      data ima/(0.0d0,1.0d0)/

      nmax = nn+1
  
      nlege = nmax + 10
      lw = (nlege+1)**2*4
      allocate(wlege(lw),ynm(0:nmax,0:nmax),ynmd(0:nmax,0:nmax))
      call ylgndrfwini(nlege,wlege,lw,lused)
      
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        call cart2polar(xyzs(1,i),r,thet,phi)
        ctheta = cos(thet)
        rx = sin(thet)*cos(phi)
        ry = sin(thet)*sin(phi)
        rz = cos(thet)

        thetx = cos(thet)*cos(phi)
        thety = cos(thet)*sin(phi)
        thetz = -sin(thet)

        phix = -sin(phi)
        phiy = cos(phi)
        phiz = 0

        call ylgndr2sfw(nmax,ctheta,ynm,ynmd,wlege,nlege)

        vtmp(1) = x/r
        vtmp(2) = y/r
        vtmp(3) = z/r

        if(mm.eq.0) then

          ynms(i) = ynm(nn,0)
          vynm(1,i) = ynm(nn,0)*rx 
          vynm(2,i) = ynm(nn,0)*ry 
          vynm(3,i) = ynm(nn,0)*rz

          psinm(1,i) = -sin(thet)*ynmd(nn,0)*thetx
          psinm(2,i) = -sin(thet)*ynmd(nn,0)*thety
          psinm(3,i) = -sin(thet)*ynmd(nn,0)*thetz
        else
          zr = ynm(nn,abs(mm))*sin(thet)*exp(ima*mm*phi)
          ynms(i) = zr
          vynm(1,i) = zr*rx
          vynm(2,i) = zr*ry 
          vynm(3,i) = zr*rz

          zt = -ynmd(nn,abs(mm))*exp(ima*mm*phi)
          zp = ima*mm*ynm(nn,abs(mm))*exp(ima*mm*phi)

          psinm(1,i) = zt*thetx + zp*phix
          psinm(2,i) = zt*thety + zp*phiy
          psinm(3,i) = zt*thetz + zp*phiz
        endif
        call dzcross_prod3d(vtmp,psinm(1,i),phinm(1,i))
      enddo
       
      return
      end


