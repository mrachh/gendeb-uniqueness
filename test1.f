      implicit real *8 (a-h,o-z)
      real *8, allocatable :: argkp(:,:),argkm(:,:)
      integer, allocatable :: isin(:,:),isuccess(:,:)
      
      done = 1
      pi = atan(done)*4
      call prini(6,13)
      call prin2('pi=*',pi,1)
c
c  nlat should be odd
c
      nlat = 101
      allocate(argkp(nlat,nlat),argkm(nlat,nlat),isin(nlat,nlat))
      do i=1,nlat
        do j=1,nlat
          argkp(j,i) = (i-1.0d0)/(nlat-1)*pi
          argkm(j,i) = (j-1.0d0)/(nlat-1)*pi
          call test_isin(argkp(j,i),argkm(j,i),isin(j,i))
          if(isin(j,i).eq.1) then
            alpha = argkm(j,i) - argkp(j,i)
            if(alpha.gt.0.and.alpha.le.pi/2.and.argkp(j,i).gt.0.and.
     1        argkp(j,i).lt.pi/2) thet = argkp(j,i)
            if(alpha.gt.0.and.alpha.le.pi/2.and.argkp(j,i).ge.pi/2.and.
     1        argkp(j,i).lt.pi) thet = argkp(j,i)-pi + 2*alpha
            if(alpha.ge.-pi/2.and.alpha.lt.0.and.argkp(j,i).gt.pi/2.and.
     1        argkp(j,i).lt.pi) thet = argkp(j,i)-pi
            if(alpha.ge.-pi/2.and.alpha.lt.0.and.argkp(j,i).gt.-alpha
     1        .and.argkp(j,i).le.pi/2) thet = argkp(j,i)+ 2*alpha
            if(alpha.eq.0) thet = 0

            if(abs(argkp(j,i)).le.1.0d-16.and.argkm(j,i).lt.pi/2)
     1         thet=0
            if(abs(argkp(j,i)).le.1.0d-16.and.argkm(j,i).gt.pi/2)
     1         thet=pi
            if(abs(argkm(j,i)).le.1.0d-16) thet = -argkp(j,i)
            if(abs(argkm(j,i)-pi).le.1.0d-16) thet = pi-argkp(j,i)
            if(abs(argkp(j,i)-pi).le.1.0d-16.and.argkm(j,i).lt.pi/2)
     1         thet=pi
            if(abs(argkp(j,i)-pi).le.1.0d-16.and.argkm(j,i).gt.pi/2)
     1         thet=0
             
            thet2 = thet
            call comp_thet(argkp(j,i),argkm(j,i),thet)
            if(abs(thet-thet2).gt.1.0d-16) then
cc              call prinf('something bad happeneed=*',i,0)
cc              call prin2('argkp=*',argkp(j,i),1)
cc              call prin2('argkm=*',argkm(j,i),1)
cc              call prin2('thet=*',thet,1)
            endif



c
c  now that theta is defined, test all inequalities
c
c
            isuc = 1
            r1 = -thet + argkp(j,i) + 2*alpha
            if(r1+pi.le.1.0d-15) r1 = r1 + 2*pi
            if(r1-2*pi.ge.-1.0d-15) r1 = r1-2*pi
            if(r1.lt.0.or.r1.gt.pi) then
              call prinf('cond 1 failed for*',i,0)
              call prin2('argkp=*',argkp(j,i),1)
              call prin2('argkm=*',argkm(j,i),1)
              call prin2('alpha=*',alpha,1)
              call prin2('thet=*',thet,1)
              call prin2('r1=*',r1,1)
              stop
            endif

            r1 = thet + argkp(j,i)
            if(r1+pi.le.1.0d-15) r1 = r1 + 2*pi
            if(r1-2*pi.ge.-1.0d-15) r1 = r1-2*pi
            if(r1.lt.0.or.r1.gt.pi) then
              call prinf('cond 2 failed for*',i,0)
              call prin2('argkp=*',argkp(j,i),1)
              call prin2('argkm=*',argkm(j,i),1)
              call prin2('alpha=*',alpha,1)
              call prin2('thet=*',thet,1)
              call prin2('r1=*',r1,1)
              stop
            endif

            r1 = -thet + argkp(j,i)
            if(r1+pi.le.1.0d-15) r1 = r1 + 2*pi
            if(r1-2*pi.ge.-1.0d-15) r1 = r1-2*pi
            if(r1.lt.0.or.r1.gt.pi) then
              call prinf('cond 3 failed for*',i,0)
              call prin2('argkp=*',argkp(j,i),1)
              call prin2('argkm=*',argkm(j,i),1)
              call prin2('alpha=*',alpha,1)
              call prin2('thet=*',thet,1)
              call prin2('r1=*',r1,1)
              stop
            endif
          endif
        enddo
      enddo

      stop
      end


      
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




