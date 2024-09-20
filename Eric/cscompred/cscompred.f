      PROGRAM CSCOMPRED
      IMPLICIT none
 
      real*8 ein(20000),e(20000),ep(20000),th(20000),nu(20000)
      real*8 q2(20000),w2(20000),cs(20000),css(20000),cserr(20000)
      real*8 flux,x(20000),eps,norm(20),rat,erat,sys,eb,theta
      real*8 kappa,sin2,tan2,t1(20000),t2(20000),t3(20000)
      real*8 csmod(20000),A(20000),Z(20000),foc(20000)
      real*8 ev,epv,q2v,w2v,epsv,kappav,fluxv,nuccs,nuccstot
      real*8 f1,f2,fl,f1qe,flqe,r,rqe,sigt,sigl,sigm,f2qe,f2mec  
      real*8 alpha,pi,pi2,mp,mp2,res,veff,chi2,q2min,q2max,w2max
      real*8 xvalc(45) / 
c     & 0.94187E-01,0.10367E+02,0.13523E+00,0.68783E+01,0.76964E+00,
c     & 0.74171E+00,0.20269E+01,0.20269E+01,0.69603E+00,-.41246E+01,
c     & 0.98677E+00,0.97242E+00,0.10338E+01,0.98935E+00,0.10000E+01,
c     & 0.10041E+01,0.97018E+00,0.10073E+01,0.98889E+00,0.99529E+00,
c     & 0.10000E+01,0.99757E+00,0.10040E+01,0.10191E+01,0.10083E+01,
c     & 0.81949E+00,0.00000E+00,-.95810E+00,0.90514E+00,0.18136E+01,
c     & 0.21799E+01,0.21799E+01,0.27720E+01,0.47927E+00,0.22800E+00,
c     & 0.10133E-01,0.26823E+00,0.38365E-01,0.71421E-01,0.88152E-01,
c     & 0.74219E-01,0.82864E-01,0.29490E+00,0.66002E+00,0.81059E+01 /
     & 0.88963E-01,0.10663E+02,0.12860E+00,0.67797E+01,0.79358E+00,
     & 0.76437E-01,0.87115E+01,0.18087E+01,0.74344E+00,-.37314E+01,
     & 0.99011E+00,0.97805E+00,0.10320E+01,0.99223E+00,0.10000E+01,
     & 0.10115E+01,0.97463E+00,0.10072E+01,0.98889E+00,0.99346E+00,
     & 0.10138E+01,0.99819E+00,0.10030E+01,0.10126E+01,0.10046E+01,
     & 0.80164E+00,0.11295E-05,-.10046E+01,0.96819E+00,0.20656E+01,
     & 0.22689E+01,0.22699E+01,0.31164E+01,0.53181E+00,0.21558E+00,
     & 0.96730E-02,0.26133E+00,0.33069E-01,0.56432E-01,0.68785E-02,
     & 0.89007E-02,0.70430E-01,0.25760E+00,-.11541E-01,0.98134E-02 /
      

      integer i,j,k,ntot,set(20000),type/1/ 
      real*4 x4,a4,emcfac,fitemc
      logical thend/.false./
      logical coulomb/.true./
      logical new,wr/.true./
      LOGICAL GOODFIT/.true./  


      character*40 filename


      read(5,*) eb,theta
      

c     filename = '12C-dec21.dat'
c      filename = '12C-feb22.dat'
       filename = '12C-april23.dat'
      
      

      w2max = 4.0

      mp = .938272
      mp2 = mp*mp

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi

      ntot = 0
      i=0

      do i=1,15
        norm(i) = xvalc(10+i)
      enddo

      do i=16,20
         norm(i) = 1.0
      enddo

*****    Read in data   *****

      open(unit=15,file=filename,status='old')
      do i=1,19
        read(15,*)    !!!  Read header  !!!
      enddo

      i = 1
      dowhile(.not.thend)

       read(15,*,end=2000) Z(i),A(i),e(i),th(i),nu(i),cs(i),cserr(i),
     &   set(i)      

       cs(i) = cs(i)/12.0
       cserr(i) = cserr(i)/12.0
       
       ep(i) = e(i)-nu(i)

c       e(i) = e(i)-0.0015
       
       sys = 0.03
       if(set(i).EQ.8.OR.set(i).EQ.11) then
         sys = 0.025
       elseif(set(i).EQ.12) then
         sys = 0.005
       endif
       cserr(i) = sqrt(cserr(i)*cserr(i)+sys*sys*cs(i)*cs(i))

      
       nu(i) = e(i)-ep(i)
       q2(i) = 4.*e(i)*ep(i)*sin(th(i)*pi/180./2.)*sin(th(i)*pi/180./2.)
       w2(i) = mp2+2.*mp*nu(i)-q2(i)
       x(i) = q2(i)/2./mp/nu(i)

       ntot = i

c       write(6,*) i,w2(i),q2(i),cs(i),cserr(i),set(i),ntot
       i=i+1 

      enddo

 2000 thend = .true.


CCC///    Now do sorting in whatever variable    ///CCC
       k = 0
       chi2 = 0.0
       do j=1,ntot
         wr = .false.            !!! if true then include data  !!!
         new = .true.
         wr = .true.
         k = k+1

         if(coulomb) then
          call vcoul(A(j),Z(j),veff)
          foc(j) = 1.0D0 + veff/e(j)
          ev = e(j) + veff     
          epv = ep(j) + veff
         endif
         sin2 = sin(th(j)*pi/180./2.)*sin(th(j)*pi/180./2.)
         tan2 = sin2/(1.-sin2)        
        
         kappa = abs(w2(j)-mp2)/2./mp
        
         eps = 1./(1. + 2.*(nu(j)*nu(j)+q2(j))/q2(j)*tan2)
         flux = alpha*kappa/(2.*pi*pi*q2(j))*ep(j)/e(j)/(1.-eps) 
         
         cs(j) = norm(set(j))*cs(j)
         cserr(j) = norm(set(j))*cserr(j)

         cs(j) = cs(j)/flux/1000.0*abs(w2(j)-mp2)
         cserr(j) = cserr(j)/flux/1000.0*abs(w2(j)-mp2)


         
c         cs(j) = cs(j)/flux/1000.0
c         cserr(j) = cserr(j)/flux/1000.0

c         cs(j) = cs(j)*abs(w2(j)-mp2)/8.0d0/pi2/alpha/0.3894e3
c         cserr(j) = cserr(j)*abs(w2(j)-mp2)/8.0d0/pi2/alpha/0.3894e3


         q2v = 4.0*ev*epv*sin2
         epsv = 1./(1. + 2.*(nu(j)*nu(j)+q2v)/q2v*tan2)
         w2v = mp2+2*mp*nu(j)-q2v
         kappav = abs(w2v-mp2)/2./mp
         fluxv = alpha*kappav/(2.*pi*pi*q2v)*epv/ev/(1.-epsv) 
         
c         cserr(j) = sqrt(cserr(j)**2+0.04*0.04*cs(j)*cs(j))

c         call f1f2in09(Z(j),A(j),q2(j),w2(j),xvalc,f1,f2,r)  !!! includes MEC !!!
c         call f1f2qe09(Z(j),A(j),q2(j),w2(j),xvalc,f1qe,f2qe)

         cs(j) = cs(j)/1000.0
         cserr(j) = cserr(j)/1000.0

c         cs(j) = cs(j)*abs(w2v-mp2)/8.0d0/pi2/alpha/0.3894e3
c         cserr(j) = cserr(j)*abs(w2v-mp2)/8.0d0/pi2/alpha/0.3894e3

         
         call csfitcomp(w2v,q2v,A(j),Z(j),xvalc,type,sigt,sigl)

         sigt = 0.3894e3*8.0d0*pi2*alpha/abs(w2v-mp2)*sigt
         sigL =  0.3894e3*8.0d0*pi2*alpha/abs(w2v-mp2)*sigL
         
         sigm = fluxv*(sigt+epsv*sigl)
c         sigm = (sigt+epsv*sigl)
         sigm = sigm/12.0

         sigm = sigm/fluxv/1000.0*abs(w2v-mp2)
         
         nuccstot = 0.0
         
         do k=2,21
          call nuccs12cs(Z(j),A(j),ev,epv,th(j),k,nuccs)
          nuccstot = nuccstot+nuccs
         enddo
         nuccstot = nuccstot/1000.0

         nuccstot = nuccstot/fluxv*abs(w2v-mp2)/1000.0

c         nuccstot = 0.0  !!! check
         
         sigm = sigm+nuccstot
         
         sigm = sigm*foc(j)*foc(j)

         res = cs(j)-sigm
         chi2 = res*res/cserr(j)/cserr(j)

         rat = cs(j)/sigm
         erat = cserr(j)/cs(j)*rat

c         if(abs(e(j)-eb).LT.0.003.AND.abs(theta-th(j)).LT.02) 
c     &      write(6,*) e(j),ep(j),th(j),w2(j),cs(j),cserr(j),sigm
         
c         write(6,*) e(j),th(j),w2(j),w2v,q2(j),q2v,cs(j),sigm,rat
           
c         endif
c         if(set(j).EQ.6) wr = .false.

         if(wr.AND.abs(eb-e(j)).LT.0.001.AND.abs(th(j)-theta).LT.0.03)
     &      write(6,3000) e(j),ep(j),th(j),w2(j),q2(j),
     &            eps,cs(j),cserr(j),sigm,rat,erat,set(j)
           
       enddo

 3000 format(6f9.4,5e13.4,1i4)
      end






