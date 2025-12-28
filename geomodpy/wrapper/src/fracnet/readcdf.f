      subroutine readcdf(datafl,ivr,ivrw,nd,MAXD,vr,cdf)
c-----------------------------------------------------------------------
c
c                   Read Data and Establish a CDF
c                   *****************************
c
c
c Original: C.V. Deutsch                              Date: October 1991
c-----------------------------------------------------------------------
      real      vr(MAXD),cdf(MAXD),var(20)
      character datafl*40
c
c Open the data file:
c
      lin  = 3
      nd   = 0
      tcdf = 0.0
      open(lin,file=datafl,status='OLD')
      read(lin,*,err=96)
      read(lin,*,err=96) nvari
      do 1 i=1,nvari
 1    read(lin,*,err=96)
c
c Read all the data until the end of the file:
c
      do 2 i=1,MAXD
            read(lin,*,end=3,err=96) (var(j),j=1,nvari)
            vrt = var(ivr)
            if(vrt.ge.-9999.0.and.vrt.le.9999.0) then
                  nd     = nd + 1
                  vr(nd) = vrt
                  if(ivrw.le.0) then
                        cdf(nd) = 1.0
                  else
                        cdf(nd) = var(ivrw)
                  endif
                  tcdf = tcdf + cdf(nd)
            endif
 2    continue
      write(*,*) ' WARNING: May have left some data out ',MAXD
 3    close(lin)
      if(tcdf.le.0.00001) go to 96
c
c Sort the data:
c
      call sortem(1,nd,vr,1,cdf,c,d,e,f,g,h)
c
c Initialize CDF, loop over all data and cutoffs computing CDF:
c
      tcdf   = 1.0 / tcdf
      cdf(1) = cdf(1) * tcdf
      do 4 i=2,nd
            cdf(i) = cdf(i-1) + cdf(i) * tcdf
 4    continue
      return
c
c Error in an Input File Somewhere:
c
 96   stop 'ERROR in data file!'
      end

