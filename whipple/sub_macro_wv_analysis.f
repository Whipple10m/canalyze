c May 2000
c this subroutine contains the subroutines it calls
c it uses CERN library for the Gaussian distributions in the
c noise simulations (from EARTH link with -L/home2/cern/99/lib -lmathlib) 


c OUTPUT: data_adc(inner_pix) with ADC values for the selected 
c                            pixels and 0 for the rest (old adc_sel)


      subroutine sub_macro_wv_analysis(inner_pix,!number of data
     &                                 data_adc,!data array (inner_pix)
     &                                 pixel_size,!in arcmin
     &                                 ped,!pedestal array (inner_pix)
     &                                 ped_disp,!pedest disp arr (inner_pix)
     &                                 x,y,!coord arrays (inner_pix)--degs
     &                                 iflag,!intg-error flag
c     &                                 adc_sel,!selected adc(inner_pix)--intg
     &                                 ind_noise_sim) !ind noise simul done


      implicit none
      integer*4 inner_pix
c      integer*2 adc(inner_pix)
      real*4 data_adc(inner_pix)
      real*4 pixel_size
      real*4 ped(inner_pix)
      real*4 ped_disp(inner_pix)
      real*4 x(inner_pix),y(inner_pix)
      integer*4 iflag
c      integer*2 adc_sel(inner_pix)
      integer*4 i
c      real*4 mhw_arcmin_mh1
      real*4 mhw_arcmin_mh2
      real*4 mhw_arcmin_mh3
      real*4 mhw_arcmin_mh4
      real*4 mhw_arcmin_mh5
c      real*4 wv_adc_mh1(inner_pix)
      real*4 wv_adc_mh2(inner_pix)
      real*4 wv_adc_mh3(inner_pix)
      real*4 wv_adc_mh4(inner_pix)
      real*4 wv_adc_mh5(inner_pix)
      integer*4 ind_noise_sim


c real numbers requiered for the MexHat convolution
c      do i=1,inner_pix
c         data_adc(i)=real(adc(i))
c      enddo


c getting MexHat wv coeffs for 5 different widths for the data       

c      mhw_arcmin_mh1=pixel_size
c      call sub_conv_mexhat_arr(inner_pix,data_adc,x,y,
c     &                         pixel_size,mhw_arcmin_mh1,wv_adc_mh1)

      mhw_arcmin_mh2=2.*pixel_size
      call sub_conv_mexhat_arr(inner_pix,data_adc,x,y,
     &                         pixel_size,mhw_arcmin_mh2,wv_adc_mh2)
      mhw_arcmin_mh3=3.*pixel_size
      call sub_conv_mexhat_arr(inner_pix,data_adc,x,y,
     &                         pixel_size,mhw_arcmin_mh3,wv_adc_mh3)
      mhw_arcmin_mh4=4.*pixel_size
      call sub_conv_mexhat_arr(inner_pix,data_adc,x,y,
     &                         pixel_size,mhw_arcmin_mh4,wv_adc_mh4)
      mhw_arcmin_mh5=5.*pixel_size
      call sub_conv_mexhat_arr(inner_pix,data_adc,x,y,
     &                         pixel_size,mhw_arcmin_mh5,wv_adc_mh5)


      
c calling the subroutine to do the wavelet analysis


      call sub_wv_analysis_prob_wvcoeffs(inner_pix,
     &                       pixel_size,
     &                       data_adc,
c     &                       wv_adc_mh1,
     &                       wv_adc_mh2,
     &                       wv_adc_mh3,
     &                       wv_adc_mh4,
     &                       wv_adc_mh5,
     &                       ped,
     &                       ped_disp,
     &                       x,y,
c     &                       adc_sel,
     &                       ind_noise_sim)


c     subtracting pedestal for input image
      do i=1,inner_pix
         data_adc(i)=data_adc(i)-ped(i)
         if (data_adc(i).lt.0.) then
            data_adc(i) = 0.0
         endif
      enddo

      return

      end


c********************SUBROUTINES of sub_macro_wv_analysis*************c

      subroutine sub_conv_mexhat_arr(inner_pix,!number of data
     &                               data, !data array (inner_pix)---real
     &                               x,y, !corrd arrays (inner_pix) in degs
     &                               pix_size, !pixel_size in arcmin
     &                               mhw_arcmin, !MexHat width in arcmin
     &                               wv) !wv_coeffs array (inner_pix)

c convolution of an hexagonal camara array with the Mexican Hat 
c wavelet characterized by the width 



      implicit none
      integer*4 inner_pix
      integer*4 i,j
      real*4 x_start
      real*4 data(inner_pix)
      real*4 x(inner_pix)
      real*4 y(inner_pix)
      integer*4 pix(inner_pix)
      real*4 mhw_arcmin
      real*4 mhw
      real*4 pix_size
      real*4 pi
      real*4 wv(inner_pix)
      real*4 const
      real*4 dist2
      real*4 sum
      integer*4 int_nwidth




      do i=1,inner_pix
            x(i)=x(i)*60./pix_size
            y(i)=y(i)*60./pix_size
            pix(i)=i
      enddo



c convolving with the Mexican Hat 

      mhw=mhw_arcmin/pix_size


      pi=3.141592654e0
      do i=1,inner_pix
         wv(i)=0.
         const=1./(mhw*sqrt(2.*pi))
         do j=1,inner_pix
            dist2=(x(i)-x(j))**2.+(y(i)-y(j))**2.
            sum=data(j)*(2.-(dist2/mhw**2.))/exp(dist2/(2.*mhw**2.))
            wv(i)=wv(i)+sum
         enddo
         wv(i)=const*wv(i)
      enddo

      do i=1,inner_pix
            x(i)=x(i)*pix_size/60.
            y(i)=y(i)*pix_size/60.
      enddo


 100   continue
      return

      end





      subroutine sub_wv_analysis_prob_wvcoeffs(
     &                       inner_pix,!number of data
     &                       pix_size,!arcmin
     &                       data_adc,!data arr(inner_pix)
c     &                       wv_adc_mh1,!wv data arr(inner_pix)
     &                       wv_adc_mh2,!wv data arr(inner_pix)
     &                       wv_adc_mh3,!wv data arr(inner_pix)
     &                       wv_adc_mh4,!wv data arr(inner_pix)
     &                       wv_adc_mh5,!wv data arr(inner_pix)
     &                       ped,!pedestal arr (inner_pix)
     &                       ped_disp,!pedest disp arr (inner_pix)
     &                       x,y,!coord arrays (inner_pix)--degs
c     &                       adc_sel,!selected ADC (inner_pix)--intg
     &                       ind_noise_sim)!noise sims and distrib only once
c wv analysis of a signal plus noise map based on the distribution 
c function of the wv coeffs at each scale at each pixel obtained from
c several noise simulations (all produced from the same pedestal
c disperision map that is in the signal plus noise map)

      implicit none
      integer*4 inner_pix
      real*4 pixel_size
c      integer*2 adc(inner_pix)
      real*4 data_adc(inner_pix)
c      real*4 wv_adc_mh1(inner_pix)
      real*4 wv_adc_mh2(inner_pix)
      real*4 wv_adc_mh3(inner_pix)
      real*4 wv_adc_mh4(inner_pix)
      real*4 wv_adc_mh5(inner_pix)
      real*4 ped(inner_pix),ped_disp(inner_pix)
      real*4 x(inner_pix),y(inner_pix)
c      integer*2 adc_sel(inner_pix)
      real*4 noise(inner_pix)
      real*4 wv(inner_pix)
      integer*4 nsim
      parameter (nsim=300) !number of noise simulations to 
!                           generate wv distributions
c      real*4 wv_mh1(inner_pix,nsim),save wv_mh1
      real*4 wv_mh2(inner_pix,nsim),save wv_mh2
      real*4 wv_mh3(inner_pix,nsim),save wv_mh3
      real*4 wv_mh4(inner_pix,nsim),save wv_mh4
      real*4 wv_mh5(inner_pix,nsim),save wv_mh5
      integer*4 i,j
      real*4 pix_size
      integer*4 ind_change
      real*4 width_arcmin
      integer*4 int_width
      integer*4 ind_noise_sim

      if (ind_noise_sim.eq.1) then
      
      do i=1,nsim
         
c         print*,'!!!!!!!!!!!!!!SIMULATION NUMBER ----->',i
         
      call sub_get_noise_arr(inner_pix,ped,ped_disp,noise)

! R1=pix_size
c      width_arcmin=1.0*pix_size
c      print*,'convolution with MexHat width',width_arcmin
c      call sub_conv_mexhat_arr(inner_pix,noise,x,y,
c     &                         pix_size,width_arcmin,wv)
c      do j=1,inner_pix
c         wv_mh1(j,i)=wv(j)
c      enddo
! R2=2*pix_size
      width_arcmin=2.0*pix_size
c      print*,'convolution with MexHat width',width_arcmin
      call sub_conv_mexhat_arr(inner_pix,noise,x,y,
     &                         pix_size,width_arcmin,wv)
      do j=1,inner_pix
         wv_mh2(j,i)=wv(j)
      enddo
! R3=3.*pix_size
      width_arcmin=3.0*pix_size
c      print*,'convolution with MexHat width',width_arcmin
      call sub_conv_mexhat_arr(inner_pix,noise,x,y,
     &                         pix_size,width_arcmin,wv)
      do j=1,inner_pix
         wv_mh3(j,i)=wv(j)
      enddo
! R4=4.*pix_size
      width_arcmin=4.0*pix_size
c      print*,'convolution with MexHat width',width_arcmin
      call sub_conv_mexhat_arr(inner_pix,noise,x,y,
     &                         pix_size,width_arcmin,wv)
      do j=1,inner_pix
         wv_mh4(j,i)=wv(j)
      enddo
! R5=5.*pix_size
      width_arcmin=5.0*pix_size
c      print*,'convolution with MexHat width',width_arcmin
      call sub_conv_mexhat_arr(inner_pix,noise,x,y,
     &                         pix_size,width_arcmin,wv)
      do j=1,inner_pix
         wv_mh5(j,i)=wv(j)
      enddo


      enddo


      endif


      call sub_get_prob_pixel_arr(inner_pix,nsim,
c     &                            wv_mh1,
     &                            wv_mh2,wv_mh3,wv_mh4,wv_mh5,
     &                            data_adc,
c     &                            wv_adc_mh1,
     &                            wv_adc_mh2,wv_adc_mh3,
     &                            wv_adc_mh4,wv_adc_mh5,
c     &                            adc_sel,
     &                            ind_noise_sim)



 100  continue
      end





c****SUBROUTINES (inside sub_wv_analysis_prob_wvcoeffs*************c

      subroutine sub_get_noise_arr(inner_pix,!number of data
     &                         ped,!pedestal array (inner_pix)
     &                         ped_disp,!pedestal deviation arr (inner_pix)
     &                         noise)! noise arr (inner_pix)
      implicit none
      integer*4 inner_pix
      integer*4 i 
      real*4 ped(inner_pix),ped_disp(inner_pix)
      real*4 noise(inner_pix)
      integer*4 int_noise
c      real*4 xdummy,gauss,x_gauss
      real*4 gauss,x_gauss
      integer*4 xdummy,save
c      external gauss
      real*8 gasdev3
      external gasdev3
      data xdummy /-13/

c we will generate noise only in the inner pixels      
c noise will be generated at ecah pixel from a Gaussian
c distribution of mean ped(i) and dispersion ped_disp(i)


c      call ranstart

      do i=1,inner_pix
c            x_gauss=gauss(xdummy) !glenn's function
            x_gauss=real(gasdev3(xdummy))
            noise(i)=ped(i)+x_gauss*ped_disp(i)
            int_noise=int(noise(i))
            noise(i)=real(int_noise)
      enddo

c      call ranend

      return

      end
      
c***********c
c functions c
c***********c
ccccccccccc gaussian generator from numerical recipes

c-----*-------------------------------------------------------------

      FUNCTION GASDEV3(IDUM)
      implicit double precision (a-h,o-z)
      DATA ISET/0/
c      print*,'beg gasdev3'
      IF (ISET.EQ.0) THEN
1       V1=2.*RAN3(IDUM)-1.
        V2=2.*RAN3(IDUM)-1.
        R=V1**2+V2**2
        IF(R.GE.1.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV3=V2*FAC
        ISET=1
      ELSE
        GASDEV3=GSET
        ISET=0
      ENDIF
c      print*,'end of gasdev3'
      RETURN
      END

c-----*-------------------------------------------------------------

      function ran3(idum)
      implicit double precision (a-h,o-z)
      integer idum
c     REAL MBIG,MSEED,MZ
      integer mbig, mseed, mz
      real ran3, fac
c     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      parameter (mbig = 1000000000, mseed = 161803398, mz = 0, fac = 1.
     & / mbig)
      integer i, iff, ii, inext, inextp, k
c     REAL mj,mk,ma(55)
      integer mj, mk, ma(55)
      save ma, inextp, inext, iff
      data iff / 0 /
c 13 "ran3.for"
      if ((idum .lt. 0) .or. (iff .eq. 0)) then
      iff = 1
      mj = mseed - iabs(idum)
      mj = mod(mj,mbig)
      ma(55) = mj
      mk = 1
      do 11 i = 1, 54
      ii = mod(21 * i,55)
      ma(ii) = mk
      mk = mj - mk
      if (mk .lt. mz) mk = mk + mbig
      mj = ma(ii)
   11 continue
      do 13 k = 1, 4
      do 12 i = 1, 55
      ma(i) = ma(i) - ma(1 + mod(i + 30,55))
      if (ma(i) .lt. mz) ma(i) = ma(i) + mbig
   12 continue
   13 continue
      inext = 0
      inextp = 31
      idum = 1
      end if
      inext = inext + 1
      if (inext .eq. 56) inext = 1
      inextp = inextp + 1
      if (inextp .eq. 56) inextp = 1
      mj = ma(inext) - ma(inextp)
      if (mj .lt. mz) mj = mj + mbig
      ma(inext) = mj
      ran3 = mj * fac
c      print*,'end of ran3'
      return 
c  (C) Copr. 1986-92 Numerical Recipes Software 2114.
c 45 "ran3.for"
      end

      subroutine sub_get_prob_pixel_arr(inner_pix,!number of data
     &       nsim,!number of noise simulations
c     &       wv_mh1,
     &       wv_mh2,wv_mh3,wv_mh4,wv_mh5,!wv noise (inner_pix,nsim)
     &       data_adc,!data array (inner_pix)
c     &       wv_adc_mh1,
     &       wv_adc_mh2,wv_adc_mh3,wv_adc_mh4,wv_adc_mh5,!wv data
c     &       adc_sel,!output selected pixels (inner_pix)--intg
     &       ind_noise_sim)!calculating only the first time 

c get the probability for each pixel from the wvcoeffs 
c distributions calculated at each scale by the subroutine
c distrib_wvcoeffs_pixel !TO LINK WITH !!!!

      implicit none
      integer*4 inner_pix
      integer*4 nsim
c      real*4 wv_mh1(inner_pix,nsim)
      real*4 wv_mh2(inner_pix,nsim)
      real*4 wv_mh3(inner_pix,nsim)
      real*4 wv_mh4(inner_pix,nsim)
      real*4 wv_mh5(inner_pix,nsim)
c      integer*2 adc(inner_pix)
      real*4 data_adc(inner_pix)
c      real*4 wv_adc_mh1(inner_pix)
      real*4 wv_adc_mh2(inner_pix)
      real*4 wv_adc_mh3(inner_pix)
      real*4 wv_adc_mh4(inner_pix)
      real*4 wv_adc_mh5(inner_pix)
c      integer*2 adc_sel(inner_pix)
c outputs for the subroutine
c      real*4 xdist_mh1(inner_pix,nsim),save xdist_mh1
      real*4 xdist_mh2(inner_pix,nsim),save xdist_mh2
      real*4 xdist_mh3(inner_pix,nsim),save xdist_mh3
      real*4 xdist_mh4(inner_pix,nsim),save xdist_mh4
      real*4 xdist_mh5(inner_pix,nsim),save xdist_mh5
c      real*4 ndist_mh1(inner_pix,nsim),save ndist_mh1
      real*4 ndist_mh2(inner_pix,nsim),save ndist_mh2
      real*4 ndist_mh3(inner_pix,nsim),save ndist_mh3
      real*4 ndist_mh4(inner_pix,nsim),save ndist_mh4
      real*4 ndist_mh5(inner_pix,nsim),save ndist_mh5
c      real*4 val_mh1
      real*4 val_mh2
      real*4 val_mh3
      real*4 val_mh4
      real*4 val_mh5
      real*4 val_mh5_array(inner_pix)
      real*4 val,val_signal
      real*4 val_signal_array(inner_pix)
c      real*4 val_out(inner_pix)
      integer*4 i,k,j
c      real*4 prob_mh1(inner_pix)
      real*4 prob_mh2(inner_pix)
      real*4 prob_mh3(inner_pix)
      real*4 prob_mh4(inner_pix)
      real*4 prob_mh5(inner_pix)
      integer*4 noise_count,signal_count,real_signal_count
      real*4 mean_mh5,mean2_mh5,rms_mh5
      integer*4 count_out,count_out_signal
      integer*4 ind_prob_mh2
      integer*4 ind_prob_mh3
      integer*4 ind_prob_mh4
      integer*4 ind_prob_mh5
      integer*4 ind_noise_sim


c      open(15,file='prob_wv_pixel_scale.dat',status='unknown',
c     &        access='append')




      if (ind_noise_sim.eq.1) then

      call distrib_wvcoeffs_pixel_arr(inner_pix,nsim,
c     &                                wv_mh1,
     &                                wv_mh2,wv_mh3,
     &                                wv_mh4,wv_mh5,
c     &                                ndist_mh1,xdist_mh1,
     &                                ndist_mh2,xdist_mh2,
     &                                ndist_mh3,xdist_mh3,
     &                                ndist_mh4,xdist_mh4,
     &                                ndist_mh5,xdist_mh5)
      
      endif

      noise_count=0
      signal_count=0
      real_signal_count=0
      mean_mh5=0.
      mean2_mh5=0.
      do i=1,inner_pix
c         val=data_adc(i)
c         val_mh1=wv_adc_mh1(i)
         val_mh2=wv_adc_mh2(i)
         val_mh3=wv_adc_mh3(i)
         val_mh4=wv_adc_mh4(i)
         val_mh5=wv_adc_mh5(i)
         val_mh5_array(i)=val_mh5
         mean_mh5=mean_mh5+val_mh5
         mean2_mh5=mean2_mh5+val_mh5**2.
c get the probabilities of val at each scale
c         prob_mh1(i)=0.
         prob_mh2(i)=0.
         prob_mh3(i)=0.
         prob_mh4(i)=0.
         prob_mh5(i)=0.
c         do k=1,nsim
c            if (val_mh1.le.xdist_mh1(i,k)) then
c               prob_mh1(i)=ndist_mh1(i,k)
c               goto 101
c            else
c               continue
c            endif
c         enddo
c 101     continue
         do k=1,nsim
            if (val_mh2.le.xdist_mh2(i,k)) then
               prob_mh2(i)=ndist_mh2(i,k)
               goto 102
            else
               continue
            endif
         enddo
 102     continue
         do k=1,nsim
            if (val_mh3.le.xdist_mh3(i,k)) then
               prob_mh3(i)=ndist_mh3(i,k)
               goto 103
            else
               continue
            endif
         enddo
 103     continue
         do k=1,nsim
            if (val_mh4.le.xdist_mh4(i,k)) then
               prob_mh4(i)=ndist_mh4(i,k)
               goto 104
            else
               continue
            endif
         enddo
 104     continue
         do k=1,nsim
            if (val_mh5.le.xdist_mh5(i,k)) then
               prob_mh5(i)=ndist_mh5(i,k)
               goto 105
            else
               continue
            endif
         enddo
 105     continue

         if ((prob_mh2(i).lt.0.05.or.
     &           prob_mh2(i).gt.0.95).and.
     &          (prob_mh3(i).eq.0.0).and.
     &          (prob_mh4(i).eq.0.0).and.
     &          (prob_mh5(i).eq.0.0)) then 

           signal_count=signal_count+1
c           val_out(i)=data_adc(i)

        else

c           val_out(i)=0.
           data_adc(i)=0.

        endif

c           adc_sel(i)=int(val_out(i))


c (May-16-2000) writting the information about the probability of
c               the image wavelet coefficients at each pixel at
c               ecah scale

c         if ((prob_mh2(i).lt.0.05.or.
c     &           prob_mh2(i).gt.0.95)) then
c            ind_prob_mh2=1
c         else 
c            ind_prob_mh2=0
c         endif
c         if (prob_mh3(i).eq.0.0) then
c            ind_prob_mh3=1
c         else 
c            ind_prob_mh3=0
c         endif
c         if (prob_mh4(i).eq.0.0) then
c            ind_prob_mh4=1
c         else 
c            ind_prob_mh4=0
c         endif
c         if (prob_mh5(i).eq.0.0) then 
c           ind_prob_mh5=1
c         else 
c            ind_prob_mh5=0
c         endif
c         
c         write(15,*)adc_sel(i),ind_prob_mh2,ind_prob_mh3,
c     &                         ind_prob_mh4,ind_prob_mh5




      enddo

c      close(15)
      

      return
      end








      subroutine distrib_wvcoeffs_pixel_arr(inner_pix,!number of data
     &            nsim,!number of noise simulations
c     &            wv_mh1,
     &            wv_mh2,wv_mh3,
     &            wv_mh4,wv_mh5,!wv noise (inner_pix,nsim)
c     &            ndist_mh1,xdist_mh1, !probability params/pixel (inner_pix)
     &            ndist_mh2,xdist_mh2,
     &            ndist_mh3,xdist_mh3,
     &            ndist_mh4,xdist_mh4,
     &            ndist_mh5,xdist_mh5)

      

c determines the wv coeffs distribution at each pixel
c at each scale. 


      implicit none
      integer*4 inner_pix
      integer*4 nsim
      integer*4 i,k,k1,k2
c      real*4 wv_mh1(inner_pix,nsim)
      real*4 wv_mh2(inner_pix,nsim)
      real*4 wv_mh3(inner_pix,nsim)
      real*4 wv_mh4(inner_pix,nsim)
      real*4 wv_mh5(inner_pix,nsim)
      real*4 a(nsim)
      real*4 delta
c      real*4 wv_mh1_min,wv_mh1_max
      real*4 wv_mh2_min,wv_mh2_max
      real*4 wv_mh3_min,wv_mh3_max
      real*4 wv_mh4_min,wv_mh4_max
      real*4 wv_mh5_min,wv_mh5_max
      real*4 x(nsim+1),xmed(nsim+1)
c      real*4 xdist_mh1(inner_pix,nsim)
      real*4 xdist_mh2(inner_pix,nsim)
      real*4 xdist_mh3(inner_pix,nsim)
      real*4 xdist_mh4(inner_pix,nsim)
      real*4 xdist_mh5(inner_pix,nsim)
      real*4 n(nsim+1),nsum
c      real*4 ndist_mh1(inner_pix,nsim)
      real*4 ndist_mh2(inner_pix,nsim)
      real*4 ndist_mh3(inner_pix,nsim)
      real*4 ndist_mh4(inner_pix,nsim)
      real*4 ndist_mh5(inner_pix,nsim)


      



      do i=1,inner_pix
c pixel fixed!
c            do k=1,nsim
c               a(k)=wv_mh1(i,k)
c            enddo
c            call sort(nsim,a)
c            do k=1,nsim
c               wv_mh1(i,k)=a(k) !sorted in ascending order
c            enddo
            do k=1,nsim
               a(k)=wv_mh2(i,k)
            enddo
            call sort(nsim,a)
            do k=1,nsim
               wv_mh2(i,k)=a(k) !sorted in ascending order
            enddo
            do k=1,nsim
               a(k)=wv_mh3(i,k)
            enddo
            call sort(nsim,a)
            do k=1,nsim
               wv_mh3(i,k)=a(k) !sorted in ascending order
            enddo
            do k=1,nsim
               a(k)=wv_mh4(i,k)
            enddo
            call sort(nsim,a)
            do k=1,nsim
               wv_mh4(i,k)=a(k) !sorted in ascending order
            enddo
            do k=1,nsim
               a(k)=wv_mh5(i,k)
            enddo
            call sort(nsim,a)
            do k=1,nsim
               wv_mh5(i,k)=a(k) !sorted in ascending order
            enddo
c generating the cumulative distributions
c mh1
c            wv_mh1_min=wv_mh1(i,1)
c            wv_mh1_max=wv_mh1(i,nsim)
c            delta=(wv_mh1_max-wv_mh1_min)/real(nsim)
c            do k=1,nsim+1
c               if (k.eq.1) then
c                  x(k)=wv_mh1_min
c               else 
c                  x(k)=x(k-1)+delta
c               endif
c            enddo
c            do k=1,nsim+1
c               n(k)=0
c            enddo
c            do k1=1,nsim
c               do k2=2,nsim+1
c                  if (wv_mh1(i,k1).ge.x(k2-1).and.
c     &                wv_mh1(i,k1).lt.x(k2)) then
c                     n(k2-1)=n(k2-1)+1
c                     goto 10
c                  else
c                     continue
c                   endif
c               enddo
c 10            continue
c            enddo
c            nsum=0.
c            do k=1,nsim+1
c               nsum=nsum+n(k)
c            enddo
c            n(1)=n(1)/nsum
c            do k=2,nsim+1
c               n(k)=(n(k)/nsum)+n(k-1)
c            enddo
c    the x of the distribution would be
c            do k=2,nsim+1
c               xmed(k-1)=x(k-1)+((x(k)-x(k-1))/2.)
c            enddo
cc the distribution at pixel (ij) would be
c            do k=1,nsim
c               xdist_mh1(i,k)=xmed(k)
c               ndist_mh1(i,k)=n(k)
c            enddo
c
c mh2
            wv_mh2_min=wv_mh2(i,1)
            wv_mh2_max=wv_mh2(i,nsim)
            delta=(wv_mh2_max-wv_mh2_min)/real(nsim)
            do k=1,nsim+1
               if (k.eq.1) then
                  x(k)=wv_mh2_min
               else 
                  x(k)=x(k-1)+delta
               endif
            enddo
            do k=1,nsim+1
               n(k)=0
            enddo
            do k1=1,nsim
               do k2=2,nsim+1
                  if (wv_mh2(i,k1).ge.x(k2-1).and.
     &                wv_mh2(i,k1).lt.x(k2)) then
                     n(k2-1)=n(k2-1)+1
                     goto 20
                  else
                     continue
                  endif
               enddo
 20            continue
            enddo
            nsum=0.
            do k=1,nsim+1
               nsum=nsum+n(k)
            enddo
            n(1)=n(1)/nsum
            do k=2,nsim+1
               n(k)=(n(k)/nsum)+n(k-1)
            enddo
c    the x of the distribution would be
            do k=2,nsim+1
               xmed(k-1)=x(k-1)+((x(k)-x(k-1))/2.)
            enddo
c the distribution at pixel (ij) would be
            do k=1,nsim
               xdist_mh2(i,k)=xmed(k)
               ndist_mh2(i,k)=n(k)
            enddo
c
c mh3
            wv_mh3_min=wv_mh3(i,1)
            wv_mh3_max=wv_mh3(i,nsim)
            delta=(wv_mh3_max-wv_mh3_min)/real(nsim)
            do k=1,nsim+1
               if (k.eq.1) then
                  x(k)=wv_mh3_min
               else 
                  x(k)=x(k-1)+delta
               endif
            enddo
            do k=1,nsim+1
               n(k)=0
            enddo
            do k1=1,nsim
               do k2=2,nsim+1
                  if (wv_mh3(i,k1).ge.x(k2-1).and.
     &                wv_mh3(i,k1).lt.x(k2)) then
                     n(k2-1)=n(k2-1)+1
                     goto 30
                  else
                     continue
                  endif
               enddo
 30            continue
            enddo
            nsum=0.
            do k=1,nsim+1
               nsum=nsum+n(k)
            enddo
            n(1)=n(1)/nsum
            do k=2,nsim+1
               n(k)=(n(k)/nsum)+n(k-1)
            enddo
c    the x of the distribution would be
            do k=2,nsim+1
               xmed(k-1)=x(k-1)+((x(k)-x(k-1))/2.)
            enddo
c the distribution at pixel (ij) would be
            do k=1,nsim
               xdist_mh3(i,k)=xmed(k)
               ndist_mh3(i,k)=n(k)
            enddo
c
c mh4
            wv_mh4_min=wv_mh4(i,1)
            wv_mh4_max=wv_mh4(i,nsim)
            delta=(wv_mh4_max-wv_mh4_min)/real(nsim)
            do k=1,nsim+1
               if (k.eq.1) then
                  x(k)=wv_mh4_min
               else 
                  x(k)=x(k-1)+delta
               endif
            enddo
            do k=1,nsim+1
               n(k)=0
            enddo
            do k1=1,nsim
               do k2=2,nsim+1
                  if (wv_mh4(i,k1).ge.x(k2-1).and.
     &                wv_mh4(i,k1).lt.x(k2)) then
                     n(k2-1)=n(k2-1)+1
                     goto 40
                  else
                     continue
                  endif
               enddo
 40            continue
            enddo
            nsum=0.
            do k=1,nsim+1
               nsum=nsum+n(k)
            enddo
            n(1)=n(1)/nsum
            do k=2,nsim+1
               n(k)=(n(k)/nsum)+n(k-1)
            enddo
c    the x of the distribution would be
            do k=2,nsim+1
               xmed(k-1)=x(k-1)+((x(k)-x(k-1))/2.)
            enddo
c the distribution at pixel (ij) would be
            do k=1,nsim
               xdist_mh4(i,k)=xmed(k)
               ndist_mh4(i,k)=n(k)
            enddo
c
c mh5
            wv_mh5_min=wv_mh5(i,1)
            wv_mh5_max=wv_mh5(i,nsim)
            delta=(wv_mh5_max-wv_mh5_min)/real(nsim)
            do k=1,nsim+1
               if (k.eq.1) then
                  x(k)=wv_mh5_min
               else 
                  x(k)=x(k-1)+delta
               endif
            enddo
            do k=1,nsim+1
               n(k)=0
            enddo
            do k1=1,nsim
               do k2=2,nsim+1
                  if (wv_mh5(i,k1).ge.x(k2-1).and.
     &                wv_mh5(i,k1).lt.x(k2)) then
                     n(k2-1)=n(k2-1)+1
                     goto 50
                  else
                     continue
                  endif
               enddo
 50            continue
            enddo
            nsum=0.
            do k=1,nsim+1
               nsum=nsum+n(k)
            enddo
            n(1)=n(1)/nsum
            do k=2,nsim+1
               n(k)=(n(k)/nsum)+n(k-1)
            enddo
c    the x of the distribution would be
            do k=2,nsim+1
               xmed(k-1)=x(k-1)+((x(k)-x(k-1))/2.)
            enddo
c the distribution at pixel (ij) would be
            do k=1,nsim
               xdist_mh5(i,k)=xmed(k)
               ndist_mh5(i,k)=n(k)
            enddo
c end loop pixel fixed!
      enddo






      return

      end






      subroutine sort(nsim,a)
      dimension a(nsim)
      n=nsim
      l=n/2+1
      ir=n
 10   continue
      if (l.gt.1) then
         l=l-1
         rra=a(l)
      else
         rra=a(ir)
         a(ir)=a(1)
         ir=ir-1
         if (ir.eq.1) then
            a(1)=rra
            return
         endif
      endif
      i=l
      j=l+l
 20   if (j.le.ir) then
         if (j.lt.ir) then
            if (a(j).lt.a(j+1))j=j+1
         endif
         if (rra.lt.a(j)) then
            a(i)=a(j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
      goto 20
      endif
      a(i)=rra
      goto 10
      end
















c*********END SUBROUTINES inside sub_wv_analysis_prob_wvcoeffs**********c




c**********END SUBROUTINES of sub_macro_wv_analysis*************c
