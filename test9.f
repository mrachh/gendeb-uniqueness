      implicit real *8 (a-h,o-z)
      complex *16 zkm,zkp,ima,f1,f2,frjd,frhd,z0,z1,zroot,froot
      complex *16 fjsp(0:100),fjdsp(0:100),fjsm(0:100),fjdsm(0:100)
      complex *16 fhsp(0:100),fhdsp(0:100),fhsm(0:100),fhdsm(0:100)
      complex *16 frh,frj

      complex *16 znull(2),zmat(2,2)

      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      
      amin = -10.0d0
      amax = 10.0d0

      bmin = -10.0d0
      bmax = 10.0d0
c      amin = -0.55686d0 - 1.0d-3
c      amax = -0.55686d0 + 1.0d-3
      
c      bmin = 1.0572d0 - 1.0d-3 
c      bmax = 1.0572d0 + 1.0d-3 

      n = 2
      njh = n + 3
      zkm = 1.23d0 + 0.5d0*ima 
      rscale = 1.0d0
      ifder = 1
      call besseljs3d(njh,zkm,rscale,fjsm,ifder,fjdsm)
      frj = fjsm(n)*zkm
      frjd = fjsm(n) + zkm*fjdsm(n)
      print *, frjd
      call h3dall(njh,zkm,rscale,fhsm,ifder,fhdsm)


      goto 1111
      print *, "njh=",njh
      open(unit=33,file='test9.dat')

      
      nlat = 300
      do i=1,nlat
        if(mod(i,10).eq.1) print *, i
        do j=1,nlat
          a = amin + (amax-amin)*(i-1)/(nlat-1)
          b = bmax + (bmin-bmax)*(j-1)/(nlat-1)
          zkp = a + ima*b
          call besseljs3d(njh,zkp,rscale,fjsp,ifder,fjdsp)
          call h3dall(njh,zkp,rscale,fhsp,ifder,fhdsp)

          frhd = fhsp(n) + zkp*fhdsp(n)
          f1 = frhd*fjsm(n) - frjd*fhsp(n)
          f2 = zkm**2*frhd*fjsm(n) - zkp**2*frjd*fhsp(n)

          write(33,'(6(2x,e11.5))') a,b,real(f1),imag(f1),
     1       real(f2),imag(f2)
          
        enddo
      enddo

      close(33)
      stop
 1111 continue

      npt = 10000
      open(unit=33,file='revplasmonic_pde.dat')
      do i=1,npt
        tt = (i+0.0d0)/(npt+0.0d0)*10.0d0
        zkm = tt
        call besseljs3d(njh,zkm,rscale,fjsm,ifder,fjdsm)
        frj = fjsm(n)*zkm
        frjd = fjsm(n) + zkm*fjdsm(n)
        zkp = tt*ima*sqrt(1.1838d0)
        z0 = zkp+0.01
        z1 = zkp-0.01
      
        call find_root_pde(z0,z1,n,zkm,frjd,fjsm(n),zroot,froot)
        call prin2_long('zroot=*',zroot,2)
        call prin2_long('froot=*',froot,2)

        zkp = zroot
        call h3dall(njh,zkp,rscale,fhsp,ifder,fhdsp)
        frh = fhsp(n)*zkp
        frhd = fhsp(n) + zkp*fhdsp(n)

        znull(1) = frj
        znull(2) = -frh

        ra = sqrt(abs(znull(1))**2 + abs(znull(2))**2)
        znull(1) = znull(1)/ra
        znull(2) = znull(2)/ra

        zmat(1,1) = frhd/zkp
        zmat(1,2) = frjd/zkm
        zmat(2,1) = frh
        zmat(2,2) = frj

        rra = abs(zmat(1,1)*znull(1) + zmat(1,2)*znull(2))
        call prin2('error in null vector=*',rra,1)
        write(33,*) tt,real(zroot),imag(zroot),abs(froot),rra
      enddo



      stop
      end



      subroutine find_root_pde(z0,z1,n,zkm,frjd,fjm,zroot,froot)
      implicit real *8 (a-h,o-z)
      complex *16 z0,z1,frjd,fjm,zroot,froot,f0,f1,f2,z2,zkm

      call pde_fun(z0,n,zkm,frjd,fjm,f0)
      print *, "z0=",z0
      print *, "abs(f0)=",abs(f0)

      call pde_fun(z1,n,zkm,frjd,fjm,f1)
      print *, "z1=",z1
      print *, "abs(f1)=",abs(f1)
      maxit = 100
      thresh = 1.0d-12
      ifdone = 0
      do ii=1,maxit
        if(abs(z1-z0).lt.1.0d-12) return
        z2 = z1 - f1*(z1-z0)/(f1-f0)
        call pde_fun(z2,n,zkm,frjd,fjm,f2)
        zroot = z2
        froot = f2
        if(abs(f2).lt.thresh) ifdone = ifdone + 1
        if(ifdone.gt.3) then
          call prin2('fun val at root=*',abs(f2),1)
        endif
        z0 = z1
        f0 = f1
        z1 = z2
        f1 = f2
      enddo


      return
      end



      subroutine pde_fun(z,n,zkm,frjd,fjm,f0)
      implicit real *8 (a-h,o-z)
      complex *16 z,frjd,fjm,f0,zkm,frhd
      complex *16 fhs(0:100),fhds(0:100)

      rscale = 1.0d0
      ifder = 1
      njh = n + 2

      print *, "n=",n
      call h3dall(njh,z,rscale,fhs,ifder,fhds)
      frhd = fhs(n) + z*fhds(n)
      f0 = zkm**2*frhd*fjm - z**2*frjd*fhs(n)

      return
      end
