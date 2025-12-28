	program fracnet  
c-----------------------------------------------------------------------
c                  3D simulation of fracture networks
c                  **********************************
c Author: E. Gringarten
c-----------------------------------------------------------------------
	include 	'fracnet.inc'
c
c Input/Output units used:
c
	lin  = 1
	lout = 2
	ldbg = 3
c
c Read the parameters and data:
c
	call readparm
c
c Call fracnet for the simulation(s):
c
	call fracnetm
c
c Finished:
c
	write(*,9998) VERSION
 9998 	format(/' FRACNET Version: ',f5.3, ' Finished'/)

	stop
	end


	subroutine readparm
c-----------------------------------------------------------------------
c
c                  Initialization and Read Parameters
c                  **********************************
c-----------------------------------------------------------------------
	include 	'fracnet.inc'
	real*8 		p,acorni
	character*40 	prefl,shafl,orifl,datafl,str
	logical		testfl
c
c Note VERSION number:
c
      	write(*,9999) VERSION
 9999 	format(/' FRACNET Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,40
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ')str='fracnet.par                      '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=95) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      write(*,*) ' Starting to read input parameters'

      read(lin,*,err=95)          i_dat
      write(*,*) ' conditioning data: ',i_dat

      read(lin,'(a40)',err=95) datafl
      call chknam(datafl,40)
      write(*,*) ' data file: ',datafl

      read(lin,'(a40)',err=95) denfl
      call chknam(denfl,40)
      write(*,*) ' density file: ',denfl

      read(lin,*,err=95)       nxd, nyd, nzd
      write(*,*) ' density grid: ',nxd, nyd, nzd

      read(lin,*,err=95)       iden
      write(*,*) ' column      : ',iden

      read(lin,'(a40)',err=95) angfl
      call chknam(angfl,40)
      write(*,*) ' strikes file: ',angfl

      read(lin,*,err=95)       iang,iangw
      write(*,*) ' strike, wt:   ',iang,iangw

      read(lin,'(a40)',err=95) dipfl
      call chknam(dipfl,40)
      write(*,*) ' dip file: ',dipfl

      read(lin,*,err=95)       idip,idipw
      write(*,*) ' dip, wt:   ',idip,idipw

      read(lin,*,err=95)          i_spd
      write(*,*) ' condition spacing: ',i_spd

      read(lin,'(a40)',err=95) spdfl
      call chknam(spdfl,40)
      write(*,*) ' spacing file: ',spdfl

      read(lin,*,err=95)       ispd,ispdw
      write(*,*) ' spacing, wt:   ',ispd,ispdw

c      read(lin,*,err=95)       i_corr
c      write(*,*) ' correlation: ',i_corr

      read(lin,*,err=95)          i_ms
      write(*,*) ' min. spacing: ',i_ms

      read(lin,*,err=95)          i_lmax, lmax
      write(*,*) ' max. length:  ',i_lmax, lmax

      read(lin,*,err=95)          i_len
      write(*,*) ' condition length: ',i_len

      read(lin,'(a40)',err=95) lenfl
      call chknam(spdfl,40)
      write(*,*) ' length file: ',lenfl

      read(lin,*,err=95)       ilen,ilenw
      write(*,*) ' length, wt:   ',ilen,ilenw

      read(lin,*,err=95)          i_pre
      write(*,*) ' existing primary set: ',i_pre

      read(lin,*,err=95)          i_hie
      write(*,*) ' hierarchical model: ',i_hie

      read(lin,'(a40)',err=95) prefl
      call chknam(prefl,40)
      write(*,*) ' primary set file: ',prefl

      read(lin,*,err=95)          i_prob, prob
      write(*,*) ' prob.of abuting:  ',i_prob, prob

      read(lin,*,err=95)         i_maxgrow
      write(*,*) ' grow to max extend', i_maxgrow

      read(lin,*,err=95)          i_sha
      write(*,*) ' existing shales: ',i_sha

      read(lin,'(a40)',err=95) shafl
      call chknam(shafl,40)
      write(*,*) ' shale file: ',shafl

      read(lin,*,err=95)          i_sprob, sprob
      write(*,*) ' prob.of stopping:  ',i_sprob, sprob

      read(lin,*,err=95)          i_ori
      write(*,*) ' orientation field: ',i_ori

      read(lin,'(a40)',err=95) orifl
      call chknam(orifl,40)
      write(*,*) ' orientation file: ',orifl

      read(lin,*,err=95)       iori
      write(*,*) ' column      : ',iori

      read(lin,'(a40)',err=95) outfl
      call chknam(outfl,40)
      write(*,*) ' output file: ',outfl

c      read(lin,*,err=95)       i_format
c      write(*,*) ' output format: ',i_format

      read(lin,*,err=95)       ixv(1)
      write(*,*) ' random number seed: ',ixv(1)
      do i=1,1000
             p = acorni(idum)
      end do


      read(lin,*,err=95)       nsim
      write(*,*) ' nummber of simulations: ',nsim

      read(lin,*,err=95)       nx,xsiz
      write(*,*) ' nx, xsiz:    ',nx,xsiz
      inx = nx / nxd
      xmn = real(inx)/2.

      read(lin,*,err=95)       ny,ysiz
      write(*,*) ' ny, ysiz:    ',ny,ysiz
      iny = ny / nyd
      ymn = real(iny)/2.

      read(lin,*,err=95)       nz,zsiz
      write(*,*) ' nz, zsiz:    ',nz,zsiz
      inz = nz / nzd
      zmn = real(inz)/2.

      nxy = nx*ny
      nxyz = nx*ny*nz

      write(*,*)

      close(lin)

c
c Perform some quick error checking:
c
      testfl = .false.
      if(nx.gt.MAXX.or.ny.gt.MAXY.or.nz.gt.MAXZ) then
            write(*,*) 'ERROR: available grid size: ',MAXX,MAXY,MAXZ
            write(*,*) '       you have asked for : ',nx,ny,nz
            testfl = .true.
      end if
      if(testfl) STOP

c
c Establish cdf's
c
	if(i_ori.eq.0) then
           call readcdf(angfl,iang,iangw,nang,MAXANG,ang,angcdf)
        endif
        call readcdf(dipfl,idip,idipw,ndip,MAXDIP,dip,dipcdf)
        if(i_spd.eq.1) then
           call readcdf(spdfl,ispd,ispdw,nspd,MAXSPD,spd,spdcdf)
        endif
        if(i_len.eq.1) then
           call readcdf(lenfl,ilen,ilenw,nlen,MAXL,len,lencdf)
        endif
c
c Read density data
c
        open(lin,file=denfl,status='old')
        read(lin,*)
        read(lin,*) njunk
        do i= 1, njunk
                read(lin,*)
        end do
        do k = 1, nzd
                do j = 1, nyd
                        do i = 1, nxd
                        read(lin,*) (junk(ic),ic=1,iden-1),n_d(i,j,k)
                        end do
                end do
        end do
        close(lin)
c
c Read orientation data
c
	PI = 3.14159265359
	if(i_ori.eq.1) then
	open(lin,file=orifl,status='old')
 	read(lin,*)
	read(lin,*) njunk
	do i= 1, njunk
		read(lin,*)
	end do	
	do k = 1, nzd
		do j = 1, nyd
			do i = 1, nxd
			read(lin,*) (junk(ic),ic=1,iori-1),ori(i,j,k)
			ori(i,j,k) = ori(i,j,k) * PI / 180.
			end do
		end do
	end do
	close(lin)
	endif
c
c Read primary data 
c
        if(i_pre.eq.1) then
        open(lin,file=prefl,status='old')
        read(lin,*)
        read(lin,*) njunk
        do i= 1, njunk
                read(lin,*)
        end do
        do k = 1, nz 
                do j = 1, ny 
                        do i = 1, nx 
			ind = i + (j-1)*nx + (k-1)*nxy
                        read(lin,*) ipre(ind)
                        end do
                end do
        end do
        close(lin)
        endif

c
c Read shale data
c
        if(i_sha.eq.1) then
        open(lin,file=shafl,status='old')
        read(lin,*)
        read(lin,*) njunk
        do i= 1, njunk
                read(lin,*)
        end do
        do k = 1, nzd
		kk = k * inz
                do j = 1, nyd
                        do i = 1, nxd
                        read(lin,*) njunk     
			do jj = (j-1)*iny+1, j*iny
			do ii = (i-1)*inx+1, i*inx
                        ind = ii + (jj-1)*nx + (kk-1)*nxy
                        isha(ind) = njunk 
c			ioff(ind) = njunk
                        end do
                        end do
                        end do
                end do
        end do
        close(lin)
        endif
c
c Read conditioning data
c
	if(i_dat.eq.1) then
	do i=1,nxyz
 		ival(i) =  0
	end do
	open(lin,file=datafl,status='old')
	read(lin,*)
        read(lin,*) njunk
        do i= 1, njunk
                read(lin,*)
        end do 
	nd = 0 
	do idata = 1, maxdat
	read(lin,*,err=6) ix,iy,iz,ico,angcond,dipcond
		ind=ix+(iy-1)*nx+(iz-1)*nxy
c		ival(ind)=ico
c		print*, ix,iy,iz,ival(ind),angcond,dipcond
                call getindx(nxd,xmn,real(inx),real(ix),id,testind)
                call getindx(nyd,ymn,real(iny),real(iy),jd,testind)
                call getindx(nzd,zmn,real(inz),real(iz),kd,testind)
c		do kd = 1, nzd
c                do jd = 1, nyd
c                do id = 1, nxd
c		if(ival(ind).eq.1) call generate(1)
 		if(ico.eq.1) call generate(1)
		if(ico.eq.0) ioff(ind) = 1
c		end do
cc		end do
c		end do
		nd = nd+1
	end do
	
 6	close(lin)
	print*, 'Number of data ',nd
	endif	
	


      return
 95   stop 'ERROR in parameter file'

      end


        subroutine fracnetm
c-----------------------------------------------------------------------
c
c                  Fracture simulation routine
c                  ***************************
c------------------------------------------------------------------------
	include		'fracnet.inc'
	real*8 		p,acorni

c
c   Loop over all simulations
c

	do isim = 1, nsim
		write(*,*) 'Working on simulation: ',isim
c
c   Initialise all points to 0
c
c		do i = 1, nxyz
c			ival(i) = 0
c		end do
c
c   Loop over all density blocks
c
        	do kd = 1, nzd
		iz = (kd-1) * inz + (inz/2)
		if(nz.eq.1) iz = 1
		do jd = 1, nyd
		iy = (jd-1) * iny + (iny/2)
		if(ny.eq.1) iy = 1
		do id = 1, nxd
		write(*,*) 'in density block ',id,jd,kd,n_d(id,jd,kd)
c
c   determine number of points to add or remove on scanline
c

                 	ncount = 0

                 	do i = 1, inx
                    	ind=(id-1)*inx+i+(iy-1)*nx+(iz-1)*nxy
                    	if(ival(ind).eq.1) ncount = ncount + 1
 			end do

                 	i_until = abs(n_d(id,jd,kd) - ncount)
                	if(ncount.eq.n_d(id,jd,kd)) then
                    		itest = 0
                   		goto 10 
                 	endif
                 	if(ncount.lt.n_d(id,jd,kd)) then
                    		itest = 1
                 	endif
                 	if(ncount.gt.n_d(id,jd,kd)) then
                    		itest = 2
                 	endif

c
c   generate or erase fractures
c

			do loop = 1, i_until

30			p = acorni(idum)
			ix = (id-1)*inx + int(p*inx) + 1

			if(itest.eq.1) then
				ind=ix+(iy-1)*nx+(iz-1)*nxy
                                if(ival(ind).eq.1) goto 30
                                if(ioff(ind).eq.1) goto 30
c
c   if conditioning to spacing distribution
c
                       if(ncount.eq.0.and.loop.ne.1) ncount=1
c                      if(i_spd.eq.1) then
                       if(i_spd.eq.1.and.ncount.ne.0) then
                          nspaces=0
                          mark1=-1
                          mark2=-1
                          flip=acorni(idum)    
                          spdraw=acorni(idum)    
c                         if(flip.ge.0.5) then
                          irlast = (int(ix/inx)+1)*inx
                          do 5000 iscan=ix,irlast
                             ind = iscan + (iy-1)*nx +(iz-1)*nxy
                             if(ival(ind).eq.1) then
                                mark1 = iscan
                                goto 5010
                             endif
                             nspaces=nspaces+1
 5000                     continue
 5010                     illast = (int(ix/inx))*inx+1
                          do 5020 iscan=ix,illast,-1
                             ind = iscan + (iy-1)*nx +(iz-1)*nxy
                             if(ival(ind).eq.1) then
                                mark2 = iscan
                                goto 5030
                             endif
                             nspaces=nspaces+1
 5020                     continue
 5030                     if(nspaces.le.spd(1)) goto 30
                          do 5040 ij = 1,nspd
                             if(spdraw.le.spdcdf(ij)) then
                                if(spd(ij).eq.nspaces) goto 5060
                                if(spd(ij).gt.nspaces) then
c                                  ij = 0
                                   spadraw=acorni(idum)    
                                   goto 5040
                                endif
                                nsps=int(spd(ij))
                                goto 5050
                             endif
 5040                     continue
 5050                     continue
                          if(mark1.eq.-1.and.mark2.eq.-1) goto 5060
                          if(flip.lt.0.5) then
                             ix = mark2 + nsps
                             if(mark2.eq.-1) ix = mark1 - nsps
                          else
                             ix = mark1 - nsps
                             if(mark1.eq.-1) ix = mark2 + nsps
                          endif
                       endif
c                      endif
 5060                  continue
                                ind=ix+(iy-1)*nx+(iz-1)*nxy
                                if(ival(ind).eq.1) goto 30
                                if(ioff(ind).eq.1) goto 30


				call generate(0)
			endif
			if(itest.eq.2) then
				if(i_maxgrow.eq.1) goto 20
				ind=ix+(iy-1)*nx+(iz-1)*nxy
				if(ival(ind).eq.0) goto 30
				call erase
			endif

 20			end do
		
 10		end do
		end do
		end do

c
c   Write results to file
c
	call output

        end do

	return
	end

        subroutine generate(inflag)
c-----------------------------------------------------------------------
c
c                  Generate fractures 
c                  ******************
c------------------------------------------------------------------------
        include         'fracnet.inc'
        real*8          p,acorni

	ind = ix  + (iy-1)*nx + (iz-1)*nxy
c
c  Draw strike and dip angles
c
	PI = 3.14159265359
	p = acorni(idum)
	if(i_ori.eq.0) then
	do i = 1,nang
		if(p.le.angcdf(i)) then
			angle(ind) = ang(i)
			angle(ind) = angle(ind) * PI / 180.
			goto 10
		endif
	end do
10	continue
	endif

	p = acorni(idum)
	do i = 1,ndip
		if(p.le.dipcdf(i)) then
			dipang(ind) = dip(i)
			dipang(ind) = dipang(ind) * PI / 180.
			goto 20
		endif
	end do
20	continue

	if(inflag.eq.1) then
		angle(ind)=angcond
		angle(ind) = angle(ind) * PI / 180.
		dipang(ind)=dipcond
		dipang(ind) = dipang(ind) * PI / 180.
	endif

c
c  Draw length if required       
c

	if(i_len.eq.1) then
	do i = 1,nlen
		if(p.le.lencdf(i)) then
                        length = len(i)
                        goto 25
		endif
	end do
25	continue
	endif
c
c  Pick y-coord start and end of fractures
c
	p = acorni(idum)
	j0 = iy - int(p*iny)
	if(j0.le.0) j0 = 1
	p = acorni(idum)
	jn = ny - iny/2 + int(p*iny)
	if(jn.gt.ny) jn = ny
c
c  Consider maximum length
c
	if(i_lmax.eq.1) then
		p = acorni(idum)
		j0 = iy - int(p*lmax)
		if(j0.le.0) j0 = 1
		jn = j0 + lmax
		if(jn.gt.ny) jn = ny
	endif
c
c  If conditioning to length distribution
c
	if(i_len.eq.1) then
	  	p = acorni(idum)
                j0 = iy - int(p*length)
		if(j0.le.0) j0 = 1
                jn = j0 + length
                if(jn.gt.ny) jn = ny
	endif

	
c
c  Grow fractures to their full extent
c
	if(i_maxgrow.eq.1) then
		j0=1 
		jn=ny
	endif

 	do k = (kd-1)*inz+1, nz

	psha=acorni(idum)
	phie=acorni(idum)

	ix2=ix
	iy2=iy
	iz2=iz
	do j=iy-1,j0,-1
		call grow(1,iflag,iflagsha)
		if(iflag.eq.1) goto 30
		if(i_ori.eq.1) then
		if(mod(j,iny).eq.0) then
			ix2 = i
			iy2 = j
		endif
		endif
  	end do
30	ix2=ix 
	iy2=iy
	iz2=iz
	do j=iy,jn     
		call grow(1,iflag,iflagsha)
		if(iflag.eq.1) goto 40
		if(i_ori.eq.1) then
		if(mod(j,iny).eq.0) then
			ix2 = i
			iy2 = j
		endif
		endif
  	end do

	if(iflagsha.eq.1) then 
		return 
	endif

40	end do

	return
	end


        subroutine erase
c-----------------------------------------------------------------------
c
c                  Erase fractures
c                  ***************
c------------------------------------------------------------------------
        include         'fracnet.inc'
        real*8          p,acorni


c
c Pick start of erasing point
c
        p = acorni(idum)
        j0 = iy - int(p*iny)
        if(j0.le.0) j0 = 1
        jn = ny
 
        do k = (kd-1)*inz+1, nz 

        phie = acorni(idum)
        psha = acorni(idum)

        ind = ix + (iy-1)*nx + (iz-1)*nxy

        ix2=ix0(ind)
        iy2=iy0(ind)
        iz2=iz0(ind)
        do j=iy0(ind)-1,j0,-1
                call grow(0,iflag,iflagsha)
                if(iflag.eq.1) goto 30
                if(i_ori.eq.1) then
                if(mod(j,iny).eq.0) then
                        ix2 = i
                        iy2 = j
                endif
                endif
        end do

30      ix2=ix0(ind)
        iy2=iy0(ind)
        iz2=iz0(ind)
        do j=iy0(ind),jn
                call grow(0,iflag,iflagsha)
                if(iflag.eq.1) goto 40
                if(i_ori.eq.1) then
                if(mod(j,iny).eq.0) then
                        ix2 = i
                        iy2 = j
                endif
                endif
        end do

 	if(iflagsha.eq.1) return

40      end do

        return
        end


        subroutine grow(inflag,iflag,iflagsha)
c-----------------------------------------------------------------------
c
c                  Grow fractures
c                  **************
c------------------------------------------------------------------------
        include         'fracnet.inc'
        real*8          p,acorni

	iflag=0
	iflagsha=0
	ind = ix  + (iy-1)*nx + (iz-1)*nxy
	if(inflag.eq.0) indang = ix0(ind) + (iy0(ind)-1)*nx
     +                                    + (iz0(ind)-1)*nxy
	if(inflag.eq.1) indang = ind
	if(i_ori.eq.1) then
                call getindx(nxd,xmn,real(inx),real(ix2),ixd,testind)
                call getindx(nyd,ymn,real(iny),real(iy2),iyd,testind)
                call getindx(nzd,zmn,real(inz),real(iz2),izd,testind)
                angle(indang) = ori(ixd,iyd,izd)
        endif
	
        i = ix2 + int(tan(angle(indang))*(j-iy2)*ysiz)
     +          + int(tan(dipang(indang))*(k-iz2)*zsiz)
        if(i.le.0.or.i.gt.nx) then
		iflag=1
		return
	endif
        ind = i  + (j-1)*nx + (k-1)*nxy

        if(i_hie.eq.1.and.ipre(ind).eq.1) then
		if(i_prob.eq.1.and.phie.gt.prob) then
			iflag = 0
		else
			iflag = 1
			return
		endif
	endif
        if(ioff(ind).eq.1) then
		iflag = 1
		return
	endif
        if(i_sha.eq.1.and.isha(ind).eq.1) then
                if(i_sprob.eq.1.and.psha.gt.sprob) then
                        iflagsha = 0
                else
                        iflagsha = 1
                        return
                endif
	endif

c
c  Generate
c
        if(inflag.eq.1) then
		ival(ind) = ival(ind)+1
		do ii=1,i_ms
			if(i+ii.le.nx) ioff(ind+ii)=1
			if(i-ii.ge.1) ioff(ind-ii)=1
		end do
	endif
			
c
c  Erase
c
        if(inflag.eq.0.and.j.ge.j0) then
		ival(ind) = 0
		do ii=1,i_ms
			if(i+ii.le.nx) ioff(ind+ii)=0
			if(i-ii.ge.1) ioff(ind-ii)=0
		end do
	endif

	if(inflag.eq.1) then
		ix0(ind) = ix
		iy0(ind) = iy
		iz0(ind) = iz
	endif

        return
        end


        subroutine output
c-----------------------------------------------------------------------
c
c                  Output results to file
c                  **********************
c------------------------------------------------------------------------
        include         'fracnet.inc'

        open(lout,file=outfl,status='UNKNOWN')
        write(lout,9997)
9997    format('Output from FRACNET',/,'1',/,'code')

	
	write(*,*) 'Writing simulation results to: ', outfl

	do i = 1, nxyz
		if(ival(i).lt.0) ival(i)=0
		if(i_pre.eq.1) ival(i) = max(ival(i),ipre(i))
		if(i_sha.eq.1.and.isha(i).eq.1) ival(i) = 2
		write(lout,'(i2)') ival(i)
	end do

	return
	end
