c
c	Static dimensioning for FRACNET
c

c
c User Adjustable Parameters:
c
	parameter(MAXX = 400, MAXY = 400 , MAXZ = 300,
     +		  MAXANG =  40, MAXDIP = 10,  MAXL = 300,
     +            MAXSPD = 10, MAXCOL = 10,
     +            MAXXD = 20, MAXYD = 20, MAXZD = 10,
     +            MAXDAT = 1400)
c
c Fixed Parameters:
c
	parameter(MXYZ=MAXX*MAXY*MAXZ,
     +		  EPSLON=1.0e-20, VERSION=1.00)
c
c ACORN parameters:
c
      parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
c
c Variable Declaration:
c
	integer		ival(MXYZ),ipre(MXYZ),ioff(MXYZ),
     +			junk(MAXCOL),isha(MXYZ),
     +			ix0(MXYZ),iy0(MXYZ),iz0(MXYZ),
     +                  n_d(MAXXD,MAXYD,MAXZD) 
	character*40 	denfl,outfl,angfl,dipfl,spdfl,lenfl
	real 		ang(MAXANG),angcdf(MAXANG),angle(MXYZ),
     +	      		dip(MAXDIP),dipcdf(MAXDIP),dipang(MXYZ),
     +	      		len(MAXL),lencdf(MAXL),
     +	      		spd(MAXSPD),spdcdf(MAXSPD),
     +			ori(MAXXD,MAXYD,MAXZD)
			

c
c Common Blocks:
c

	common /iaco/ 	ixv(MAXOP1)

	common /generl/	nx,ny,nz,xsiz,ysiz,zsiz,xmn,ymn,zmn,nxy,
     +			nxyz,lin,lout,n_d,
     + 			outfl,junk,inx,iny,inz,nsim,ival,ix,iy,iz,
     +			id,jd,kd,i,j,k,ix2,iy2,iz2,ix0,iy0,iz0,
     +			j0,jn,ipre,ioff,isha,lmax,prob,sprob,
     + 			phie,psha,angcond,dipcond

	common /flag/	i_pre,i_hie,i_maxgrow,i_lmax,i_sha,i_ms,
     +  		i_prob,i_sprob,i_len,i_spd,i_dat

	common /dens/	nxd,nyd,nzd,iden,denfl 

	common /dist/	angfl,iang,iangw,nang,ang,angcdf,
     +			dipfl,idip,idipw,ndip,dip,dipcdf,
     +			spdfl,ispd,ispdw,nspd,spd,spdcdf,
     +			lenfl,ilen,ilenw,nlen,len,lencdf,
     +			angle,dipang,i_ori,ori
