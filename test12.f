c
c  In this code, we test the signatures for the generalized debye representation
c  in both the interior and the exterior. The currents are computed analytically
c  from the charges skipping the Laplace beltrami solve in this test.
c  
c  Note that representations are given by:
c
c    E^{\pm} =  i \om \mu_{\pm} b3^{\pm} S_{k_\pm} [j^{\pm}] + 
c       a2^{\pm} \nabla S_{k_\pm} [r^{\pm}] + 
c       a3^{\pm} \nabla \times S_{k_\pm} [m^\pm]
c
c    H^{\pm} =  -i \om \ep_{\pm} a3^{\pm} S_{k_\pm} [m^{\pm}] + 
c       b2^{\pm} \nabla S_{k_\pm} [q^{\pm}] + 
c       b3^{\pm} \nabla \times S_{k_\pm} [j^\pm]
c

      implicit real *8 (a-h,o-z)
      real *8 om
      complex *16 epm(2), zmpm(2)
      complex *16 apms(2,2), bpms(2,2)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: wts(:)
      integer, allocatable :: norders(:), ixyzs(:), iptype(:)
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)

      complex *16, allocatable :: asigs(:,:,:,:)
      complex *16, allocatable :: bsigs(:,:,:,:)
      complex *16, allocatable :: dsigs(:,:,:,:)
      complex *16, allocatable :: zjsigs(:,:,:,:)
      complex *16, allocatable :: zmsigs(:,:,:,:)
      complex *16, allocatable :: esigs(:,:,:,:)
      complex *16, allocatable :: hsigs(:,:,:,:)

      complex *16, allocatable :: ynm(:,:,:)
      complex *16, allocatable :: vynm(:,:,:,:)
      complex *16, allocatable :: phinm(:,:,:,:)
      complex *16, allocatable :: psinm(:,:,:,:)

      real *8 c0(3)
      complex *16, allocatable :: zjpvals(:,:), zmpvals(:,:)
      complex *16, allocatable :: zrpvals(:), zqpvals(:)
      complex *16, allocatable :: zjmvals(:,:), zmmvals(:,:)
      complex *16, allocatable :: zrmvals(:), zqmvals(:)
      complex *16, allocatable :: zjuse(:,:), zmuse(:,:), zruse(:)
      complex *16, allocatable :: coefs(:,:,:)
      complex *16, allocatable :: ep(:,:), em(:,:)
      complex *16, allocatable :: hp(:,:), hm(:,:)
      complex *16, allocatable :: epex(:,:), emex(:,:)
      complex *16, allocatable :: hpex(:,:), hmex(:,:)
      complex *16 ima, cphi, cpsi, cy
      complex *16 zkp, zkm
      complex *16, allocatable :: wnearp(:,:), wnearm(:,:)
      data ima/(0.0d0, 1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4
c
c  get the geometry
c
      norder = 6
      a = 1.0d0
      na = 2
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
c  get the vector spherical harmonics
c

      nmax = 3
      allocate(ynm(0:nmax,-nmax:nmax,npts))
      allocate(psinm(3,0:nmax,-nmax:nmax,npts))
      allocate(phinm(3,0:nmax,-nmax:nmax,npts))
      allocate(vynm(3,0:nmax,-nmax:nmax,npts))

      call l3getsph_vec_all(nmax, 12, npts, srcvals, ynm, vynm, 
     1   psinm, phinm)

c
c
c  define the material parameters
c
      om = 1.1d0
      epm(1) = hkrand(0) + ima*hkrand(0) 
      epm(2) = hkrand(0) + ima*hkrand(0)

      epm(1) = 1.0d0

      zmpm(1) = hkrand(0) + ima*hkrand(0) 
      zmpm(2) = hkrand(0) + ima*hkrand(0)
      zmpm(1) = 1.0d0

c
c  define gendeb cosntants
c
      apms(1,1) = hkrand(0) + ima*hkrand(0)
      apms(1,2) = hkrand(0) + ima*hkrand(0)
      apms(2,1) = hkrand(0) + ima*hkrand(0)
      apms(2,2) = hkrand(0) + ima*hkrand(0)

      bpms(1,1) = hkrand(0) + ima*hkrand(0)
      bpms(1,2) = hkrand(0) + ima*hkrand(0)
      bpms(2,1) = hkrand(0) + ima*hkrand(0)
      bpms(2,2) = hkrand(0) + ima*hkrand(0)

c
c
c  get the signatures
c
      allocate(asigs(3,2,2,0:nmax))
      allocate(bsigs(3,2,2,0:nmax)) 
      allocate(dsigs(3,1,2,0:nmax)) 
      allocate(zjsigs(2,4,2,0:nmax)) 
      allocate(zmsigs(2,4,2,0:nmax)) 
      allocate(esigs(3,4,2,0:nmax)) 
      allocate(hsigs(3,4,2,0:nmax))

      print *, "building signatures"

      call prinf('nmax=*', nmax, 1)
      call prin2('om=*', om, 1)
      call prin2('epm=*', epm, 4)
      call prin2('zmpm=*', zmpm, 4)
      call prin2('apms=*', apms, 8)
      call prin2('bpms=*', bpms, 8)

      asigs = 0
      bsigs = 0
      dsigs = 0
      zjsigs = 0
      zmsigs = 0
      esigs = 0
      hsigs = 0
      
      
      call build_ops(nmax, om, epm, zmpm, apms, bpms, asigs, 
     1  bsigs, dsigs, zjsigs, zmsigs, esigs, hsigs)
      call prin2('asigs=*', asigs, 24*(nmax+1))
      call prin2('bsigs=*', bsigs, 24*(nmax+1))
      call prin2('dsigs=*', dsigs, 12*(nmax+1))
      call prin2('zjsigs=*', zjsigs, 32*(nmax+1))
      call prin2('zmsigs=*', zmsigs, 32*(nmax+1))
      call prin2('esigs=*', esigs, 48*(nmax+1))
      call prin2('hsigs=*', hsigs, 48*(nmax+1))
      

c
c  set up random coefs for spherical harmonic expansion for qpm, rpm
c    coefs(1,:,:) are coefs for rp
c    coefs(2,:,:) are coefs for rm
c    coefs(3,:,:) are coefs for qp
c    coefs(4,:,:) are coefs for qm
c  
      allocate(coefs(4,0:nmax,-nmax:nmax))
      coefs(1:4,0:nmax,-nmax:nmax) = 0
      do i=1,nmax
        do j=1,i
          coefs(1,i,j) = hkrand(0) + ima*hkrand(0)
          coefs(2,i,j) = hkrand(0) + ima*hkrand(0)
          coefs(3,i,j) = hkrand(0) + ima*hkrand(0)
          coefs(4,i,j) = hkrand(0) + ima*hkrand(0)
        enddo
      enddo

      allocate(zjpvals(3,npts), zmpvals(3,npts), zrpvals(npts))
      allocate(zqpvals(npts))

      allocate(zjmvals(3,npts), zmmvals(3,npts), zrmvals(npts))
      allocate(zqmvals(npts))


      do i=1,npts
        zrpvals(i) = 0
        zrmvals(i) = 0
        zqpvals(i) = 0
        zqmvals(i) = 0
        do j=1,nmax
          do k=1,j
            zrpvals(i) = zrpvals(i) + coefs(1,j,k)*ynm(j,k,i)
            zrmvals(i) = zrmvals(i) + coefs(2,j,k)*ynm(j,k,i)
            zqpvals(i) = zqpvals(i) + coefs(3,j,k)*ynm(j,k,i)
            zqmvals(i) = zqmvals(i) + coefs(4,j,k)*ynm(j,k,i)
          enddo
        enddo
c
c  now set up the currents based on the symbols
c
        zjpvals(1:3,i) = 0
        zjmvals(1:3,i) = 0
        zmpvals(1:3,i) = 0
        zmmvals(1:3,i) = 0
        do j=1,nmax
          do k=1,j
c
c  j^{+}
c

            cpsi = coefs(1,j,k)*zjsigs(1,1,1,j) + 
     1             coefs(2,j,k)*zjsigs(1,2,1,j) +       
     1             coefs(3,j,k)*zjsigs(1,3,1,j) +       
     1             coefs(4,j,k)*zjsigs(1,4,1,j) 

            cphi = coefs(1,j,k)*zjsigs(2,1,1,j) + 
     1             coefs(2,j,k)*zjsigs(2,2,1,j) +       
     1             coefs(3,j,k)*zjsigs(2,3,1,j) +       
     1             coefs(4,j,k)*zjsigs(2,4,1,j)

            zjpvals(1:3,i) = zjpvals(1:3,i) + cpsi*psinm(1:3,j,k,i) + 
     1                      cphi*phinm(1:3,j,k,i)
c
c  j^{-}
c
            cpsi = coefs(1,j,k)*zjsigs(1,1,2,j) + 
     1             coefs(2,j,k)*zjsigs(1,2,2,j) +       
     1             coefs(3,j,k)*zjsigs(1,3,2,j) +       
     1             coefs(4,j,k)*zjsigs(1,4,2,j) 

            cphi = coefs(1,j,k)*zjsigs(2,1,2,j) + 
     1             coefs(2,j,k)*zjsigs(2,2,2,j) +       
     1             coefs(3,j,k)*zjsigs(2,3,2,j) +       
     1             coefs(4,j,k)*zjsigs(2,4,2,j)

            zjmvals(1:3,i) = zjmvals(1:3,i) + cpsi*psinm(1:3,j,k,i) + 
     1                      cphi*phinm(1:3,j,k,i)

c
c  m^{+}
c
            cpsi = coefs(1,j,k)*zmsigs(1,1,1,j) + 
     1             coefs(2,j,k)*zmsigs(1,2,1,j) +       
     1             coefs(3,j,k)*zmsigs(1,3,1,j) +       
     1             coefs(4,j,k)*zmsigs(1,4,1,j) 

            cphi = coefs(1,j,k)*zmsigs(2,1,1,j) + 
     1             coefs(2,j,k)*zmsigs(2,2,1,j) +       
     1             coefs(3,j,k)*zmsigs(2,3,1,j) +       
     1             coefs(4,j,k)*zmsigs(2,4,1,j)

            zmpvals(1:3,i) = zmpvals(1:3,i) + cpsi*psinm(1:3,j,k,i) + 
     1                      cphi*phinm(1:3,j,k,i)

c
c  m^{-}
c
            cpsi = coefs(1,j,k)*zmsigs(1,1,2,j) + 
     1             coefs(2,j,k)*zmsigs(1,2,2,j) +       
     1             coefs(3,j,k)*zmsigs(1,3,2,j) +       
     1             coefs(4,j,k)*zmsigs(1,4,2,j) 

            cphi = coefs(1,j,k)*zmsigs(2,1,2,j) + 
     1             coefs(2,j,k)*zmsigs(2,2,2,j) +       
     1             coefs(3,j,k)*zmsigs(2,3,2,j) +       
     1             coefs(4,j,k)*zmsigs(2,4,2,j)

            zmmvals(1:3,i) = zmmvals(1:3,i) + cpsi*psinm(1:3,j,k,i) + 
     1                      cphi*phinm(1:3,j,k,i)

          enddo
        enddo
      enddo

      
      nnz = npts*npatches
      allocate(row_ptr(npts+1), col_ind(nnz), iquad(nnz+1))

      do i=1,npts
        row_ptr(i) = (i-1)*npatches+1
        do j=1,npatches
          col_ind(row_ptr(i)+j-1) = j
        enddo
      enddo
      row_ptr(npts+1) = npts*npatches + 1
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, 
     1   iquad)
      nquad = iquad(nnz+1) - 1
      allocate(wnearp(4,nquad), wnearm(4,nquad))

      zkp = om*sqrt(epm(1))*sqrt(zmpm(1))
      print *, "zkp=", zkp
      zkm = om*sqrt(epm(2))*sqrt(zmpm(2))
      
      print *, "starting exterior kernels"
      wnearp = 0
      wnearm = 0
      call get_quadrature_corrections(npatches, norders, ixyzs, iptype,
     1  npts, srccoefs, srcvals, zkp, nnz, row_ptr, col_ind, iquad, 
     2  nquad, wnearp) 
      call prin2('wnearp=*', wnearp, 24)

c      print *, "starting interior kernels"
c
c      call get_quadrature_corrections(npatches, norders, ixyzs, iptype,
c     1  npts, srccoefs, srcvals, zkm, nnz, row_ptr, col_ind, iquad, 
c     2  nquad, wnearm)
c      
      call prin2('wnearm=*', wnearm, 24)

      allocate(ep(3,npts), em(3,npts), hp(3,npts), hm(3,npts))
      allocate(epex(3,npts), emex(3,npts), hpex(3,npts), hmex(3,npts))
      allocate(zjuse(3,npts), zruse(npts), zmuse(3,npts))
c
c  Compute exterior electric field
c
      iout = 1

      print *, "computing exterior electric fields"
      do i=1,npts
c        zjuse(1:3,i) = ima*om*zmpm(1)*bpms(2,1)*zjpvals(1:3,i)
c        zruse(i) = apms(1,1)*zrpvals(i)
c        zmuse(1:3,i) = apms(2,1)*zmpvals(1:3,i)
         zjuse(1:3,i) = phinm(1:3,2,1,i)
         zruse(i) = 0
         zmuse(1:3,i) = 0
      enddo
      call prin2('zjuse=*', zjuse, 24)
      call prin2('zruse=*', zruse, 24)
      call prin2('zmuse=*', zmuse, 24)
      print *, "npts =", npts
      print *, "npatches =", npatches


     
      call compute_field(npatches, npts, ixyzs, srcvals, zjuse, zmuse, 
     1  zruse, nnz, row_ptr, col_ind, iquad, nquad, wnearp, iout, epex)
      
      call prin2('epex=*',epex,24)
c
c  test asigs by verifying epex
c
      call prin2('asigs=*',asigs(1,2,1,2), 6)
      erra = 0
      ra = 0
      do i=1,npts
        ep(1:3,i) = asigs(1,2,1,2)*psinm(1:3,2,1,i) + 
     1              asigs(2,2,1,2)*phinm(1:3,2,1,i) + 
     2              asigs(3,2,1,2)*vynm(1:3,2,1,i)
        erra = erra + abs(ep(1,i) - epex(1,i))**2
        erra = erra + abs(ep(2,i) - epex(2,i))**2
        erra = erra + abs(ep(3,i) - epex(3,i))**2

        ra = ra + abs(epex(1,i))**2
        ra = ra + abs(epex(2,i))**2
        ra = ra + abs(epex(3,i))**2
        
        if(i.le.5) then
          write(*,*) " "
          write(*,*) " "
          call prin2('ep=*',ep(1,i), 6)
          call prin2('epex=*',epex(1,i), 6)
        endif
      enddo

      erra = sqrt(erra/ra)
      print *, "error in s[psinm]=",erra
      stop


c
c  Compute exterior magnetic field
c
      do i=1,npts
        zjuse(1:3,i) = -ima*om*epm(1)*apms(2,1)*zmpvals(1:3,i)
        zruse(i) = bpms(1,1)*zqpvals(i)
        zmuse(1:3,i) = bpms(2,1)*zjpvals(1:3,i)
      enddo

      call compute_field(npatches, npts, ixyzs, srcvals, zjuse, zmuse, 
     1  zruse, nnz, row_ptr, col_ind, iquad, nquad, wnearp, iout, hpex)


      iout = -1
c
c  Compute interior electric field
c
      do i=1,npts
        zjuse(1:3,i) = ima*om*zmpm(2)*bpms(2,2)*zjmvals(1:3,i)
        zruse(i) = apms(1,2)*zrmvals(i)
        zmuse(1:3,i) = apms(2,2)*zmmvals(1:3,i)
      enddo

      call compute_field(npatches, npts, ixyzs, srcvals, zjuse, zmuse, 
     1  zruse, nnz, row_ptr, col_ind, iquad, nquad, wnearm, iout, emex)

c
c  Compute interior magnetic field
c
      do i=1,npts
        zjuse(1:3,i) = -ima*om*epm(2)*apms(2,2)*zmmvals(1:3,i)
        zruse(i) = bpms(1,2)*zqmvals(i)
        zmuse(1:3,i) = bpms(2,2)*zjmvals(1:3,i)
      enddo

      call compute_field(npatches, npts, ixyzs, srcvals, zjuse, zmuse, 
     1  zruse, nnz, row_ptr, col_ind, iquad, nquad, wnearm, iout, hmex)

c
c  Now compute the fields using the signatures
c

      do i=1,npts
        ep(1:3,i) = 0
        em(1:3,i) = 0
        hp(1:3,i) = 0
        hm(1:3,i) = 0
c
c  now set up the currents based on the symbols
c
        do j=1,nmax
          do k=1,j
c
c  e^{+}
c

            cpsi = coefs(1,j,k)*esigs(1,1,1,j) + 
     1             coefs(2,j,k)*esigs(1,2,1,j) +       
     1             coefs(3,j,k)*esigs(1,3,1,j) +       
     1             coefs(4,j,k)*esigs(1,4,1,j) 

            cphi = coefs(1,j,k)*esigs(2,1,1,j) + 
     1             coefs(2,j,k)*esigs(2,2,1,j) +       
     1             coefs(3,j,k)*esigs(2,3,1,j) +       
     1             coefs(4,j,k)*esigs(2,4,1,j)

            cy   = coefs(1,j,k)*esigs(3,1,1,j) + 
     1             coefs(2,j,k)*esigs(3,2,1,j) +       
     1             coefs(3,j,k)*esigs(3,3,1,j) +       
     1             coefs(4,j,k)*esigs(3,4,1,j)

            ep(1:3,i) = ep(1:3,i) + cpsi*psinm(1:3,j,k,i) + 
     1                  cphi*phinm(1:3,j,k,i) + 
     1                  cy*vynm(1:3,j,k,i)
c
c  e^{-}
c
            cpsi = coefs(1,j,k)*esigs(1,1,2,j) + 
     1             coefs(2,j,k)*esigs(1,2,2,j) +       
     1             coefs(3,j,k)*esigs(1,3,2,j) +       
     1             coefs(4,j,k)*esigs(1,4,2,j) 

            cphi = coefs(1,j,k)*esigs(2,1,2,j) + 
     1             coefs(2,j,k)*esigs(2,2,2,j) +       
     1             coefs(3,j,k)*esigs(2,3,2,j) +       
     1             coefs(4,j,k)*esigs(2,4,2,j)

            cy   = coefs(1,j,k)*esigs(3,1,2,j) + 
     1             coefs(2,j,k)*esigs(3,2,2,j) +       
     1             coefs(3,j,k)*esigs(3,3,2,j) +       
     1             coefs(4,j,k)*esigs(3,4,2,j)

            em(1:3,i) = em(1:3,i) + cpsi*psinm(1:3,j,k,i) + 
     1                  cphi*phinm(1:3,j,k,i) +
     1                  cy*vynm(1:3,j,k,i)

c
c  h^{+}
c
            cpsi = coefs(1,j,k)*hsigs(1,1,1,j) + 
     1             coefs(2,j,k)*hsigs(1,2,1,j) +       
     1             coefs(3,j,k)*hsigs(1,3,1,j) +       
     1             coefs(4,j,k)*hsigs(1,4,1,j) 

            cphi = coefs(1,j,k)*hsigs(2,1,1,j) + 
     1             coefs(2,j,k)*hsigs(2,2,1,j) +       
     1             coefs(3,j,k)*hsigs(2,3,1,j) +       
     1             coefs(4,j,k)*hsigs(2,4,1,j)

            cy   = coefs(1,j,k)*hsigs(3,1,1,j) + 
     1             coefs(2,j,k)*hsigs(3,2,1,j) +       
     1             coefs(3,j,k)*hsigs(3,3,1,j) +       
     1             coefs(4,j,k)*hsigs(3,4,1,j)

            hp(1:3,i) = hp(1:3,i) + cpsi*psinm(1:3,j,k,i) + 
     1                  cphi*phinm(1:3,j,k,i) + 
     1                  cy*vynm(1:3,j,k,i)

c
c  h^{-}
c
            cpsi = coefs(1,j,k)*hsigs(1,1,2,j) + 
     1             coefs(2,j,k)*hsigs(1,2,2,j) +       
     1             coefs(3,j,k)*hsigs(1,3,2,j) +       
     1             coefs(4,j,k)*hsigs(1,4,2,j) 

            cphi = coefs(1,j,k)*hsigs(2,1,2,j) + 
     1             coefs(2,j,k)*hsigs(2,2,2,j) +       
     1             coefs(3,j,k)*hsigs(2,3,2,j) +       
     1             coefs(4,j,k)*hsigs(2,4,2,j)

            cy   = coefs(1,j,k)*hsigs(3,1,2,j) + 
     1             coefs(2,j,k)*hsigs(3,2,2,j) +       
     1             coefs(3,j,k)*hsigs(3,3,2,j) +       
     1             coefs(4,j,k)*hsigs(3,4,2,j)

            hm(1:3,i) = hm(1:3,i) + cpsi*psinm(1:3,j,k,i) + 
     1                  cphi*phinm(1:3,j,k,i) + 
     1                  cy*vynm(1:3,j,k,i)

          enddo
        enddo
      enddo
      call prin2('ep=*',ep,24)

      errep = 0
      errem = 0
      errhp = 0
      errhm = 0
      rep = 0
      rem = 0
      rhp = 0
      rhm = 0

      do i=1,npts
        errep = errep + abs(ep(1,i) - epex(1,i))**2
        errep = errep + abs(ep(2,i) - epex(2,i))**2
        errep = errep + abs(ep(3,i) - epex(3,i))**2

        rep = rep + abs(epex(1,i))**2
        rep = rep + abs(epex(2,i))**2
        rep = rep + abs(epex(3,i))**2
        
        errem = errem + abs(em(1,i) - emex(1,i))**2
        errem = errem + abs(em(2,i) - emex(2,i))**2
        errem = errem + abs(em(3,i) - emex(3,i))**2

        rem = rem + abs(emex(1,i))**2
        rem = rem + abs(emex(2,i))**2
        rem = rem + abs(emex(3,i))**2
        
        errhp = errhp + abs(hp(1,i) - hpex(1,i))**2
        errhp = errhp + abs(hp(2,i) - hpex(2,i))**2
        errhp = errhp + abs(hp(3,i) - hpex(3,i))**2

        rhp = rhp + abs(hpex(1,i))**2
        rhp = rhp + abs(hpex(2,i))**2
        rhp = rhp + abs(hpex(3,i))**2
        
        errhm = errhm + abs(hm(1,i) - hmex(1,i))**2
        errhm = errhm + abs(hm(2,i) - hmex(2,i))**2
        errhm = errhm + abs(hm(3,i) - hmex(3,i))**2

        rhm = rhm + abs(hmex(1,i))**2
        rhm = rhm + abs(hmex(2,i))**2
        rhm = rhm + abs(hmex(3,i))**2
      enddo

      errep = sqrt(errep/rep)
      errem = sqrt(errem/rem)
      errhp = sqrt(errhp/rhp)
      errhm = sqrt(errhm/rhm)

      call prin2('errep = *', errep, 1)
      call prin2('errem = *', errem, 1)
      call prin2('errhp = *', errhp, 1)
      call prin2('errhm = *', errhm, 1)




      return
      end
c
c
c
c
      subroutine compute_field(npatches, npts, ixyzs, srcvals, zj, zm,
     1   zr, nnz, row_ptr, col_ind, iquad, nquad, wnear, iout, e)
c
c  This subroutine evaluates the field
c    e = S_{k}[zj] + \nabla S_{k}[zr] + \nabla \times S_{k}[zm]
c
c  If e is the electric field, then upto constant scaling factors, 
c  zj is the electric current, zm is the magnetic current, and zr
c  is the electric charge
c
c  with precomputed quadrature corrections for S_{k} and \nabla S_{k}
c
c  Input arguments:
c
c    - npatches: integer
c        number of patches
c    - npts: integer
c        total number of discretization points on the boundary
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in various
c        source arrays where information for patch i begins
c    - srcvals: double precision (12,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - zj: complex *16 (3,npts)
c        electric current
c    - zm: complex *16 (3,npts)
c        magnetic current
c    - zr: complex *16 (npts)
c        electric charge
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(npts+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer(nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear_ij array where quadrature for col_ind(i)
c        starts for the single kernel.
c    - nquad: integer
c        number of entries in wnear, i.e., number of near field entries
c        corresponding to each source target pair
c    - wnear: complex *16 (4,nquad)
c        Precomputed near field quadrature
c          * the first kernel is S_{k} 
c          * the second kernel is \partial_{x} S_{k} 
c          * the third kernel is \partial_{y} S_{k} 
c          * the fourth kernel is \partial_{z} S_{k}
c    - iout: integer
c        flag for interior or exterior problems
c          * iout = 1, corresponds to exterior problems
c          * iout = -1, corresponds to interior problems
c          * iout = 0, returns principal value part of the integrands
c    
c  Output arguments:
c    - e: complex *16 (3,npts)
c        electric field

      implicit real *8 (a-h,o-z)
      integer, intent(in) :: npatches, npts, ixyzs(npatches+1)
      real *8, intent(in) :: srcvals(12,npts)
      complex *16, intent(in) :: zj(3,npts), zr(npts), zm(3,npts)
      integer, intent(in) :: nnz, row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1), nquad, iout
      complex *16, intent(in) :: wnear(4,nquad)
      
      complex *16, intent(out) :: e(3,npts)
      complex *16 ztan, zn, zvec(3)

      zn = 0.5d0*iout
      ztan = -0.5d0*iout
      print *, "nnz=",nnz
      print *, "nquad=",nquad

      do i=1,npts
        e(1:3,i) = 0
        do j=row_ptr(i), row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            jpt = jstart+l-1
            jq = jquadstart+l-1
            e(1,i) = e(1,i) + wnear(3,jq)*zm(3,jpt) - 
     1             wnear(4,jq)*zm(2,jpt) + wnear(1,jq)*zj(1,jpt) + 
     2             wnear(2,jq)*zr(jpt)
            e(2,i) = e(2,i) + wnear(4,jq)*zm(1,jpt) - 
     1             wnear(2,jq)*zm(3,jpt) + wnear(1,jq)*zj(2,jpt) + 
     2             wnear(3,jq)*zr(jpt)
            e(3,i) = e(3,i) + wnear(2,jq)*zm(2,jpt) - 
     1             wnear(3,jq)*zm(1,jpt) + wnear(1,jq)*zj(3,jpt) + 
     2             wnear(4,jq)*zr(jpt)

          enddo
        enddo
        e(1:3,i) = e(1:3,i) + zn*srcvals(10:12,i)*zr(i)
        call dzcross_prod3d(srcvals(10:12,i), zm(1,i), zvec) 
        e(1:3,i) = e(1:3,i) + zvec(1:3)*ztan
      enddo


      return
      end

c
c
c
c
c
      subroutine get_quadrature_corrections(npatches, norders, ixyzs,
     1  iptype, npts, srccoefs, srcvals, zk, nnz, row_ptr, col_ind, 
     2  iquad, nquad, wnear)
c
c  This subroutine returns the quadrature corrections for s and grad s
c  on surface for the helmholtz greens function
c
c
c  Input arguments:
c
c    - npatches: integer
c        number of patches
c    - norders: integer(npatches)
c        order of discretization on each patch
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array where information for patch i begins
c    - iptype: integer(npatches)
c        type of patch
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        koornwinder expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$.
c    - srcvals: double precision (12,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - zk: complex *16
c        helmholtz wave number
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(npts+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer(nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear_ij array where quadrature for col_ind(i)
c        starts for the single kernel.
c    - nquad: integer
c        number of entries in wnear, i.e., number of near field entries
c        corresponding to each source target pair
c
c  Output arguments:
c    - wnear: complex *16(4,nquad)
c        Precomputed near field quadrature
c          * the first kernel is S_{k} 
c          * the second kernel is \partial_{x} S_{k} 
c          * the third kernel is \partial_{y} S_{k} 
c          * the fourth kernel is \partial_{z} S_{k} 

      implicit real *8 (a-h,o-z)
      integer, intent(in) :: npatches, norders(npatches)
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      integer, intent(in) :: npts
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      complex *16, intent(in) :: zk
      integer, intent(in) :: nnz, row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1), nquad
      complex *16, intent(out) :: wnear(4,nquad)

      complex *16, allocatable :: wneartmp(:)
      real *8 dpars
      integer ipars
      integer ipv, ndd, ndz, ndi
      real *8, allocatable :: uvs_src(:,:)
      integer, allocatable :: ipatch_id(:)
      
      external h3d_slp, h3d_sgradx, h3d_sgrady, h3d_sgradz
      
      allocate(uvs_src(2,npts), ipatch_id(npts))
      call prinf('iptype=*', iptype, npatches)
      print *, "npts=",npts
      print *, "nnz=",nnz
      print *, "nquad=",nquad
      print *, "iquad(end)=",iquad(nnz+1)
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, 
     1  ipatch_id, uvs_src)
      call prinf('ipatch_id=*',ipatch_id,npts)
      call prin2('uvs_src=*',uvs_src,2*npts)

      ndd = 0
      ndi = 0
      ndz = 1
      ipv = 1

      eps = 1.0d-8

      rfac0 = 1.25d0


      allocate(wneartmp(nquad))

      print *, "starting kernel 1"
      wneartmp(1:nquad) = 0
      nd = 12
      call zgetnearquad_ggq_guru(npatches, norders, ixyzs, iptype, 
     1  npts, srccoefs, srcvals, nd, npts, srcvals, ipatch_id, uvs_src,
     2  eps, ipv, h3d_slp, ndd, dpars, ndz, zk, ndi, ipars, nnz, 
     3  row_ptr, col_ind, iquad, rfac0, nquad, wneartmp)
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(1,i) = wneartmp(i)
      enddo
C$OMP END PARALLEL DO      

      print *, "starting kernel 2"
      wneartmp(1:nquad) = 0
      call zgetnearquad_ggq_guru(npatches, norders, ixyzs, iptype, 
     1  npts, srccoefs, srcvals, nd, npts, srcvals, ipatch_id, uvs_src,
     2  eps, ipv, h3d_sgradx, ndd, dpars, ndz, zk, ndi, ipars, nnz, 
     3  row_ptr, col_ind, iquad, rfac0, nquad, wneartmp)
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(2,i) = wneartmp(i)
      enddo
C$OMP END PARALLEL DO      
     
      print *, "Starting kernel 3"
      wneartmp(1:nquad) = 0
      call zgetnearquad_ggq_guru(npatches, norders, ixyzs, iptype, 
     1  npts, srccoefs, srcvals, nd, npts, srcvals, ipatch_id, uvs_src,
     2  eps, ipv, h3d_sgrady, ndd, dpars, ndz, zk, ndi, ipars, nnz, 
     3  row_ptr, col_ind, iquad, rfac0, nquad, wneartmp)
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(3,i) = wneartmp(i)
      enddo
C$OMP END PARALLEL DO      

      print *, "starting kernel 4"
      wneartmp(1:nquad) = 0
      call zgetnearquad_ggq_guru(npatches, norders, ixyzs, iptype, 
     1  npts, srccoefs, srcvals, nd, npts, srcvals, ipatch_id, uvs_src,
     2  eps, ipv, h3d_sgradz, ndd, dpars, ndz, zk, ndi, ipars, nnz, 
     3  row_ptr, col_ind, iquad, rfac0, nquad, wneartmp)
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(4,i) = wneartmp(i)
      enddo
C$OMP END PARALLEL DO      
     
      return
      end

c
c
c
c
c
c

      subroutine build_ops(nmax, om, epm, zmpm, apms, bpms, asigs, 
     1  bsigs, dsigs, zjsigs, zmsigs, esigs, hsigs)
c
c  This subroutine returns all the signatures in the spherical
c  harmonic basis for the generalized debye representation 
c  at non-zero frequency.       
c
c
c  Input arguments:
c    - nmax: integer
c        max spherical harmonic expansion order
c    - om: real *8
c        omega (wave number)
c    - epm: complex *16 (2)
c        permittivities, epm(1) is the exterior or + permittivity
c                        epm(2) is the interior or - permittivity
c    - zmpm: complex *16 (2)
c        permeabilities, zmpm(1) is the exterior or + permeabilities
c                        zmpm(2) is the interior or - permeabilities
c    - apms: complex *16(2,2)
c        The 'a_{2}' and 'a_{3}' coefficients in the paper
c        * apms(:,1) = (a_{2}^{+}, a_{3}^{+})
c        * apms(:,2) = (a_{2}^{-}, a_{3}^{-})
c    - bpms: complex *16(2,2)
c        The 'b_{2}' and 'b_{3}' coefficients in the paper
c        * bpms(:,1) = (b_{2}^{+}, b_{3}^{+})
c        * bpms(:,2) = (b_{2}^{-}, b_{3}^{-})
c
c  Output arguments:
c    - asigs: complex *16 (3,2,2,0:nmax)
c        signatures of S[] of a vector field
c        * asigs(:,:,1,n) is the mapping for v = S_{k+}[j_{n}], with
c          j_{n} expressed in the tangential vector spherical 
c          harmonics, and v expressed in the vector spherical harmonics
c        * asigs(:,:,2,n) is the corresponding mapping for 
c          v = S_{k-}[j_{n}]
c    - bsigs: complex *16 (3,2,2,0:nmax)
c        signatures of \nabla \times S[] of a vector field
c        * bsigs(:,:,1,n) is the mapping for 
c          v = \nabla \times S_{k+}[j_{n}], with
c          j_{n} expressed in the tangential vector spherical 
c          harmonics, and v expressed in the vector spherical harmonics
c        * bsigs(:,:,2,n) is the corresponding mapping for 
c          v = \nabla \times S_{k-}[j_{n}]
c    - dsigs: complex *16 (3,1,2,0:nmax)
c        signatures of \nabla S[] of a scalar function
c        * dsigs(:,:,1,n) is the mapping for 
c          v = \nabla S_{k+}[q_{n}], with
c          q_{n} being the scalar spherical harmonic, 
c          and v expressed in the vector spherical harmonics
c        * dsigs(:,:,2,n) is the corresponding mapping for 
c          v = \nabla S_{k-}[q_{n}]
c    - zjsigs: complex *16 (2,4,2,0:nmax)
c        coefficients that map strengths of charge densities 
c        expressed in scalar spherical harmonics, to 
c        coefficients of the electric currents in the tangetial 
c        vector spherical harmonic basis. 
c          * zjsigs(:,:,1,:) are the mappings for the exterior region
c          * zjsigs(:,:,2,:) are the mappings for the interior region
c    - zmsigs: complex *16 (2,4,2,0:nmax)
c        coefficients that map strengths of charge densities 
c        expressed in scalar spherical harmonics, to 
c        coefficients of the magnetic currents in the tangetial 
c        vector spherical harmonic basis. 
c          * zmsigs(:,:,1,:) are the mappings for the exterior region
c          * zmsigs(:,:,2,:) are the mappings for the interior region
c    - esigs: complex *16 (3,4,2,0:nmax)
c        coefficients that map strengths of charge densities 
c        expressed in scalar spherical harmonics, to 
c        coefficients of the electric field in the  
c        vector spherical harmonic basis
c          * esigs(:,:,1,:) are the mappings for the exterior region
c          * esigs(:,:,2,:) are the mappings for the interior region
c    - hsigs: complex *16 (3,4,2,0:nmax)
c        coefficients that map strengths of charge densities 
c        expressed in scalar spherical harmonics, to 
c        coefficients of the magnetic field in the  
c        vector spherical harmonic basis. 
c          * hsigs(:,:,1,:) are the mappings for the exterior region 
c          * hsigs(:,:,2,:) are the mappings for the interior region
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nmax
      real *8, intent(in) :: om
      complex *16, intent(in) :: epm(2), zmpm(2)
      complex *16, intent(in) :: apms(2,2), bpms(2,2)
      
      complex *16, intent(out) :: asigs(3,2,2,0:nmax) 
      complex *16, intent(out) :: bsigs(3,2,2,0:nmax) 
      complex *16, intent(out) :: dsigs(3,1,2,0:nmax) 
      complex *16, intent(out) :: zjsigs(2,4,2,0:nmax) 
      complex *16, intent(out) :: zmsigs(2,4,2,0:nmax) 
      complex *16, intent(out) :: esigs(3,4,2,0:nmax) 
      complex *16, intent(out) :: hsigs(3,4,2,0:nmax)

 
      complex *16, allocatable :: zjp(:), zjm(:)
      complex *16, allocatable :: zdjp(:), zdjm(:)
      complex *16, allocatable :: zjrp(:), zjrm(:)
      complex *16, allocatable :: zdjrp(:), zdjrm(:)

      complex *16, allocatable :: zhp(:), zhm(:)
      complex *16, allocatable :: zdhp(:), zdhm(:)
      complex *16, allocatable :: zhrp(:), zhrm(:)
      complex *16, allocatable :: zdhrp(:), zdhrm(:)

      complex *16 zkp, zkm
      complex *16 ima, z1
      complex *16 a3p, a2p, a3m, a2m
      complex *16 b3p, b2p, b3m, b2m
      data ima/(0.0d0,1.0d0)/


      allocate(zjp(0:nmax), zjm(0:nmax))
      allocate(zdjp(0:nmax), zdjm(0:nmax))
      allocate(zjrp(0:nmax), zjrm(0:nmax))
      allocate(zdjrp(0:nmax), zdjrm(0:nmax))

      allocate(zhp(0:nmax), zhm(0:nmax))
      allocate(zdhp(0:nmax), zdhm(0:nmax))
      allocate(zhrp(0:nmax), zhrm(0:nmax))
      allocate(zdhrp(0:nmax), zdhrm(0:nmax))

      zkp = om*sqrt(epm(1))*sqrt(zmpm(1))
      zkm = om*sqrt(epm(2))*sqrt(zmpm(2))
      call jhfuns(nmax, zkp, zjp, zdjp, zhp, zdhp, zjrp, zdjrp,
     1   zhrp, zdhrp)
      call prin2('zjp=*',zjp, 2*(nmax+1))
      call prin2('zdjp=*',zdjp, 2*(nmax+1))
      call prin2('zhp=*',zhp, 2*(nmax+1))
      call prin2('zdhp=*',zdhp, 2*(nmax+1))
      call prin2('zjrp=*',zjrp, 2*(nmax+1))
      call prin2('zdjrp=*',zdjrp, 2*(nmax+1))
      call prin2('zhrp=*',zhrp, 2*(nmax+1))
      call prin2('zdhrp=*',zdhrp,2*(nmax+1))
      call jhfuns(nmax, zkm, zjm, zdjm, zhm, zdhm, zjrm, zdjrm,
     1   zhrm, zdhrm)

      a2p = apms(1,1)
      a3p = apms(2,1)

      a2m = apms(1,2)
      a3m = apms(2,2)

      b2p = bpms(1,1)
      b3p = bpms(2,1)

      b2m = bpms(1,2)
      b3m = bpms(2,2)
      
      asigs(1:3,1:2,1:2,0:nmax) = 0
      bsigs(1:3,1:2,1:2,0:nmax) = 0
      dsigs(1:3,1,1:2,0:nmax) = 0
      zjsigs(1:2,1:4,1:2,0:nmax) = 0
      zmsigs(1:2,1:4,1:2,0:nmax) = 0
      esigs(1:3,1:4,1:2,0:nmax) = 0
      hsigs(1:3,1:4,1:2,0:nmax) = 0

      do n=1,nmax
c
c  compute asigs
c
        z1 = ima/zkp
        asigs(1,1,1,n) = z1*(zdjrp(n)*zdhrp(n) + 
     1      n*(n+1.0d0)*zjp(n)*zhp(n))
        asigs(3,1,1,n) = z1*n*(n+1)*(zjp(n)*zdhrp(n) + 
     1       zdjp(n)*zhrp(n))
        asigs(2,2,1,n) = z1*zjrp(n)*zhrp(n)

        z1 = ima/zkm
        asigs(1,1,2,n) = z1*(zdjrm(n)*zdhrm(n) + 
     1      n*(n+1.0d0)*zjm(n)*zhm(n))
        asigs(3,1,2,n) = z1*((n+1)*(zjm(n)*zdhrm(n) + zdjm(n)*zhrm(n)))
        asigs(2,2,2,n) = z1*zkm*zkm*zjm(n)*zhm(n)
c
c  compute bsigs
c
        bsigs(2,1,1,n) = -ima*zdhrp(n)*zjrp(n)
        bsigs(1,2,1,n) = -ima*zdjrp(n)*zhrp(n)
        bsigs(3,2,1,n) = -ima*n*(n+1)*zkp*zjp(n)*zhp(n)

        bsigs(2,1,2,n) = -ima*zdjrm(n)*zhrm(n)
        bsigs(1,2,2,n) = -ima*zdhrm(n)*zjrm(n)
        bsigs(3,2,2,n) = -ima*n*(n+1)*zkm*zjm(n)*zhm(n)
c
c  compute dsigs
c
        dsigs(1,1,1,n) = ima*zkp*zjp(n)*zhp(n)
        dsigs(3,1,1,n) = ima*zkp*zkp*zdhp(n)*zjp(n)

        dsigs(1,1,2,n) = ima*zkm*zjm(n)*zhm(n)
        dsigs(3,1,2,n) = ima*zkm*zkm*zdjm(n)*zhm(n)

c
c  compute zjsigs
c
        r1 = -1.0d0/n/(n+1.0d0)
        zjsigs(1,1,1,n) = r1*ima*epm(1)*a2p/b3p
        zjsigs(2,2,1,n) = -r1*ima*epm(2)*a2m/b3p

        zjsigs(1,2,2,n) = r1*ima*epm(2)*a2m/b3m
        zjsigs(2,1,2,n) = r1*ima*epm(1)*a2p/b3m

c
c  compute zmsigs
c
        zmsigs(1,3,1,n) = r1*ima*zmpm(1)*b2p/a3p
        zmsigs(2,4,1,n) = -r1*ima*zmpm(2)*b2m/a3p
       
        zmsigs(1,4,2,n) = r1*ima*zmpm(2)*b2m/a3m
        zmsigs(2,3,2,n) = r1*ima*zmpm(1)*b2p/a3m

c
c  compute esigs = i \om b3 \mu asig*jsig + 
c     a2 D e_{j}^{t} + a3 bsig*msig
c
c  and hsigs = -i \om a3p \ep asig*msig + b2 D e_{j}^{t} +
c     b3 bsig*jsig
c
        do j=1,4
          do l=1,3
            do k=1,2
              esigs(l,j,1,n) = esigs(l,j,1,n) + 
     1          ima*om*zmpm(1)*b3p*asigs(l,k,1,n)*zjsigs(k,j,1,n)
              esigs(l,j,1,n) = esigs(l,j,1,n) + 
     1          a3p*bsigs(l,k,1,n)*zmsigs(k,j,1,n)
              esigs(l,j,2,n) = esigs(l,j,2,n) + 
     1           ima*om*zmpm(2)*b3m*asigs(l,k,2,n)*zjsigs(k,j,2,n)
              esigs(l,j,2,n) = esigs(l,j,2,n) + 
     1          a3m*bsigs(l,k,2,n)*zmsigs(k,j,2,n)


              hsigs(l,j,1,n) = hsigs(l,j,1,n) - 
     1          ima*om*epm(1)*a3p*asigs(l,k,1,n)*zmsigs(k,j,1,n)
              hsigs(l,j,1,n) = hsigs(l,j,1,n) + 
     1          b3p*bsigs(l,k,1,n)*zjsigs(k,j,1,n)
              hsigs(l,j,2,n) = hsigs(l,j,2,n) - 
     1           ima*om*epm(2)*a3m*asigs(l,k,2,n)*zmsigs(k,j,2,n)
              hsigs(l,j,2,n) = hsigs(l,j,2,n) + 
     1          b3m*bsigs(l,k,2,n)*zjsigs(k,j,2,n)
            enddo
          enddo
        enddo
        esigs(1:3,1,1,n) = esigs(1:3,1,1,n) + a2p*dsigs(1:3,1,1,n)
        esigs(1:3,2,2,n) = esigs(1:3,2,2,n) + a2m*dsigs(1:3,1,2,n)
        hsigs(1:3,3,1,n) = hsigs(1:3,3,1,n) + b2p*dsigs(1:3,1,1,n)
        hsigs(1:3,4,2,n) = hsigs(1:3,4,2,n) + b2m*dsigs(1:3,1,2,n)
      enddo
      


      return
      end
c
c      
c     
c   
c
c

      subroutine jhfuns(njh, zk, fjvals, fjder, fhvals, fhder, frjvals, 
     1   frjder, frhvals, frhder)
c
c  Compute the spherical bessel and hankel functions, along with
c  their ricatti analogs at a given wave number
c 
c  Input arguments:
c    - njh: integer
c        maximum order of bessel/hankel functions 
c    - zk: complex *16
c        argument of the bessel/hankel functions
c 
c  Output arguments:
c    - fjvals: complex *16 (0:njh)
c        spherical bessel functions
c    - fjder: complex *16(0:njh)
c        derivatives of spherical bessel functions
c    - fhvals: complex *16(0:njh)
c        spherical hankel functions
c    - fhder: complex *16(0:njh)
c        derivatives of spherical hankel functions
c    - frjvals: complex *16 (0:njh)
c        Ricatti bessel functions
c    - frjder: complex *16(0:njh)
c        derivatives of ricatti bessel functions
c    - frhvals: complex *16(0:njh)
c        ricatti hankel functions
c    - frhder: complex *16(0:njh)
c        derivatives of ricatti hankel functions
c
c
c
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
      
      print *, "njh=", njh

      call besseljs3d(njh, zk, rscale, fjvals, ifder, fjder)
      call prin2('fjvals=*',fjvals,2*(njh+1))
      call prin2('fjder=*',fjder,2*(njh+1))
      call h3dall(njh, zk, rscale, fhvals, ifder, fhder)
      call prin2('fjvals=*',fhvals,2*(njh+1))
      call prin2('fjder=*',fhder,2*(njh+1))
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



