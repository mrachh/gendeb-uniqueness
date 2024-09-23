      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)
      real *8 dpars(2)

      integer, allocatable :: norders(:), ixyzs(:), iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: ynm(:,:,:)
      complex *16, allocatable :: psinm(:,:,:,:), phinm(:,:,:,:)
      complex *16, allocatable :: vynm(:,:,:,:)
      complex *16, allocatable :: einc(:,:), eincex(:,:)

      complex *16 zk

      complex *16 fjvals(0:100),fhvals(0:100),fjder(0:100),fhder(0:100)
      complex *16 frjvals(0:100),frhvals(0:100),frjder(0:100)
      complex *16 frhder(0:100)
      complex *16 z1, z2, z3, z4, zfac
      complex *16 zvec1(3),zvec2(3),zvec3(3)
      real *8 dvec1(3),dvec2(3),dvec3(3)

      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname

      real *8, allocatable :: w(:,:)
      real *8 c0(3)

      logical isout0,isout1

      complex *16 ztmp,ima
      procedure (), pointer :: fker
      external h3d_sgradx, h3d_sgrady, h3d_sgradz

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      norder = 8
      a = 1.0d0
      na = 4
      c0(1:3) = 0
      iptype0 = 1
      call get_sphere_npat_mem(a, na, c0, norder, iptype0, npatches,
     1   npts)
      
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      call get_sphere_npat(a, na, c0, norder, iptype0, npatches, npts,
     1  norders, ixyzs, iptype, srccoefs, srcvals)

      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


c
c       define rhs to be one of the ynm's
c
      nmax = 20
      allocate(ynm(0:nmax,-nmax:nmax,npts))
      allocate(psinm(3,0:nmax,-nmax:nmax,npts))
      allocate(phinm(3,0:nmax,-nmax:nmax,npts))
      allocate(vynm(3,0:nmax,-nmax:nmax,npts))

      call l3getsph_vec_all(nmax, 12, npts, srcvals, ynm, vynm, 
     1   psinm, phinm)

c
c  set the value of dzk
c
      zk = 1.0d0
      call jhfuns(nmax, zk, fjvals, fjder, fhvals, fhder, frjvals, 
     1   frjder, frhvals, frhder)


      allocate(einc(3,npts), eincex(3,npts))
      erra = 0
      ra = 0
      do i=1,npts
        eincex(1,i) = exp(ima*zk*srcvals(3,i))
        eincex(2,i) = 0
        eincex(3,i) = 0
        
        einc(1:3,i) = 0
        do n = 1,nmax
          z1 = n*(n+1.0d0)*fjvals(n)/zk
          z2 = fjder(n) + fjvals(n)/zk  
          zfac = ima**(n)*sqrt(2*n+1.0d0)/sqrt(n*(n+1.0d0))
          einc(1:3,i) = einc(1:3,i) + 
     1       zfac*(
     2       fjvals(n)*imag(phinm(1:3,n,1,i)) + 
     3       ima*(z1*real(vynm(1:3,n,1,i)) + 
     4            z2*real(psinm(1:3,n,1,i))))

        enddo
        erra = erra + abs(einc(1,i) - eincex(1,i))**2
        erra = erra + abs(einc(2,i) - eincex(2,i))**2
        erra = erra + abs(einc(3,i) - eincex(3,i))**2
        ra = ra + abs(eincex(1,i))**2
      enddo
      erra = sqrt(erra/ra)

      call prin2('einc=*', einc, 6)
      call prin2('eincex=*', eincex, 6)
      call prin2('error in pw expansion=*',erra,1)


      stop
      end

c
c      
c     
c   
c
c

      subroutine jhfuns(njh, zk, fjvals, fjder, fhvals, fhder, frjvals, 
     1   frjder, frhvals, frhder)
      implicit real *8 (a-h,o-z)
      integer njh
      complex *16 zk
      complex *16 fjvals(0:njh), fjder(0:njh)
      complex *16 frjvals(0:njh), frjder(0:njh)
      complex *16 fhvals(0:njh), fhder(0:njh)
      complex *16 frhvals(0:njh), frhder(0:njh)
      integer ifder
      real *8 rscale

      ifder = 1
      rscale = 1.0d0
      

      call besseljs3d(njh, zk, rscale, fjvals, ifder, fjder)
      call h3dall(njh, zk, rscale, fhvals, ifder, fhder)
      do i=0,njh
        frjvals(i) = zk*fjvals(i)
        frhvals(i) = zk*fhvals(i)
        frjder(i) = fjvals(i) + zk*fjder(i)
        frhder(i) = fhvals(i) + zk*fhder(i)
      enddo

      return
      end
c
c
c
c
c      
c
      subroutine l3getsph_vec_all(nmax, ndx, npts, xyzs, ynms, 
     1   vynm, psinm, phinm)
      implicit real *8 (a-h,o-z)
      integer nmax
      real *8 :: xyzs(ndx,npts)
      complex *16 ima
      complex *16 ynms(0:nmax,-nmax:nmax,npts)
      complex *16 vynm(3,0:nmax,-nmax:nmax,npts)
      complex *16 psinm(3,0:nmax,-nmax:nmax,npts)
      complex *16 phinm(3,0:nmax,-nmax:nmax,npts)
      real *8 vtmp(3)
      real *8, allocatable :: wlege(:),ynm(:,:),ynmd(:,:)
      complex *16 zr,zt,zp
      data ima/(0.0d0,1.0d0)/

  
      nlege = nmax + 10
      lw = (nlege+1)**2*4
      allocate(wlege(lw), ynm(0:nmax,0:nmax), ynmd(0:nmax,0:nmax))
      call ylgndrfwini(nlege,wlege,lw,lused)

      ynms = 0
      vynm = 0
      psinm = 0
      phinm = 0
  
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

        do nn = 0,nmax
          do mm = -nn,nn
            

            if(mm.eq.0) then
              ynms(nn,mm,i) = ynm(nn,0)
              vynm(1,nn,mm,i) = ynm(nn,0)*rx 
              vynm(2,nn,mm,i) = ynm(nn,0)*ry 
              vynm(3,nn,mm,i) = ynm(nn,0)*rz

              psinm(1,nn,mm,i) = -sin(thet)*ynmd(nn,0)*thetx
              psinm(2,nn,mm,i) = -sin(thet)*ynmd(nn,0)*thety
              psinm(3,nn,mm,i) = -sin(thet)*ynmd(nn,0)*thetz
            else
              zr = ynm(nn,abs(mm))*sin(thet)*exp(ima*mm*phi)
              ynms(nn,mm,i) = zr
              vynm(1,nn,mm,i) = zr*rx
              vynm(2,nn,mm,i) = zr*ry 
              vynm(3,nn,mm,i) = zr*rz

              zt = -ynmd(nn,abs(mm))*exp(ima*mm*phi)
              zp = ima*mm*ynm(nn,abs(mm))*exp(ima*mm*phi)

              psinm(1,nn,mm,i) = zt*thetx + zp*phix
              psinm(2,nn,mm,i) = zt*thety + zp*phiy
              psinm(3,nn,mm,i) = zt*thetz + zp*phiz
            endif

            if(mm.lt.0) then
              ynms(nn,mm,i) = (-1)**mm*ynms(nn,mm,i)
              vynm(1:3,nn,mm,i) = (-1)**mm*vynm(1:3,nn,mm,i)
              psinm(1:3,nn,mm,i) = (-1)**mm*psinm(1:3,nn,mm,i)
            endif
            call dzcross_prod3d(vtmp, psinm(1,nn,mm,i), 
     1        phinm(1,nn,mm,i))
          enddo
        enddo
      enddo
       
      return
      end



