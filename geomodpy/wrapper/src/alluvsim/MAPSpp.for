c
c Created from GSLIB code in January 2005 to perform post processing 
c with MAPS selection.
c



      module  geostat

      parameter(MAXCAT=32,EPSLON=1.0e-10,UNEST=-99.0,VERSION=1.000)

      integer test,nd,nx,ny,nz,nxy,nxyz,nsim,lin,lout,ncat,
     +        cat(MAXCAT)
     
      real    xsiz,ysiz,zsiz,xmn,ymn,zmn,pdf(MAXCAT),mpdf(MAXCAT)
     
      integer,allocatable :: sim(:,:,:),cdat(:,:,:),ixpath(:),iypath(:),
     +                       izpath(:),ixmaps(:),iymaps(:),izmaps(:),
     +                       ixw(:),iyw(:),izw(:),icatw(:),imisw(:)
      real,allocatable    :: cwt(:,:,:),wtpath(:),wtmaps(:),
     +                       wnode(:),anode(:),adata(:)
      
      end module



      program main
c-----------------------------------------------------------------------
c
c
c
c
c
c
c
c
c AUTHOR: C.V. Deutsch                                DATE: January 2005
c-----------------------------------------------------------------------
      use       geostat
      parameter(MV=512)
      real      var(MV),tvar(MV)
      real*8    p,acorni
      character wellfl*512,datafl*512,outfl*512,str*512
      logical   testfl,inflag
c
c ACORN parameters:
c
      parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      common /iaco/   ixv(MAXOP1)
c
c Fortran unit numbers needed:
c
      lin  = 1
      lout = 2
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' MAPSpp Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'MAPSpp.par          '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'MAPSpp.par          ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a512)',err=98) wellfl
      call chknam(wellfl,512)
      write(*,*) ' well data file = ',wellfl(1:40)

      read(lin,*,err=98) icx,icy,icz,icf
      write(*,*) ' columns = ',icx,icy,icz,icf

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' initial image data file = ',datafl(1:40)

      read(lin,*,err=98) icol
      write(*,*) ' column for variable = ',icol

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=98) nsim
      write(*,*) ' number of simulations = ',nsim

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' X grid specification = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' Y grid specification = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' Z grid specification = ',nz,zmn,zsiz
      nxy  = nx*ny
      nxyz = nx*ny*nz

      read(lin,*,err=98) ixv(1)
      write(*,*) ' random number seed = ',ixv(1)
      nloop = int(ixv(1))/10
      do i=1,nloop
             p = acorni(idum)
      end do

      read(lin,*,err=98) nxt,nyt,nzt
      write(*,*) ' Window size = ',nxt,nyt,nzt

      read(lin,*,err=98) wexp
      write(*,*) ' Weighting exponent = ',wexp

      read(lin,*,err=98) nxm,nym,nzm
      write(*,*) ' MAPS size = ',nxm,nym,nzm

      close(lin)
c
c Some memory allocation:
c
      allocate(sim(nx,ny,nz),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(cdat(nx,ny,nz),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(cwt(nx,ny,nz),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      sim  = -99
      cdat = -99
      cwt  = -99
c
c Make sure the initial image exists:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*)
            write(*,*) 'No initial image - nothing to do'
            stop
      end if
c
c Read the initial image:
c
      write(*,*)
      write(*,*) 'Reading initial image'
      open(lin,file=datafl,status='OLD')
      read(lin,*,err=99)
      read(lin,*,err=99) nvari
      do i=1,nvari
            read(lin,*,err=99)
      end do
      ncat = 0
      cat  = 0
      pdf  = 0.
      do k=1,nz
      do j=1,ny
      do i=1,nx
            read(lin,*,end=99,err=99) (var(l),l=1,nvari)
            icat = int(var(icol)+0.5)
            do ic=1,ncat
                  if(icat.eq.cat(ic)) then
                        ind = ic
                        go to 2
                  end if      
            end do
            ncat = ncat + 1
            if(ncat.ge.MAXCAT) then
                  write(*,*) ' Too many categories in training',ncat
                  stop
            end if      
            cat(ncat) = icat
            ind = ncat
 2          continue
            sim(i,j,k) = ind
            pdf(ind)   = pdf(ind) + 1
            if(sim(i,j,k).lt.1.or.sim(i,j,k).gt.MAXCAT) then
                  write(*,*) ' problem with facies at ',i,j,k
                  stop
            end if
      end do
      end do
      end do
      close(lin)
      do ic=1,ncat
            pdf(ic) = pdf(ic) / real(nx*ny*nz)
            write(*,*)    '  category ',cat(ic),' proportion ',pdf(ic)
      end do
c
c Make sure the well data exists:
c
      inquire(file=wellfl,exist=testfl)
      if(.not.testfl) then
            write(*,*)
            write(*,*) 'No well data - nothing to do'
            stop
      end if
c
c Read the well data:
c
      write(*,*)
      write(*,*) 'Reading well data'
c
c     Find out how many we need to allocate for:
c
      open(lin,file=wellfl,status='OLD')
      read(lin,*,err=99)
      read(lin,*,err=99) nvari
      do i=1,nvari
            read(lin,*,err=99)
      end do
      maxwd = 0
 3    read(lin,*,end=4,err=99) (var(j),j=1,nvari)
      if(var(icf).lt.tmin.or.var(icf).ge.tmax) go to 3
      maxwd = maxwd + 1
      go to 3
 4    continue      
      rewind lin
      read(lin,*,err=99)
      read(lin,*,err=99) nvari
      do i=1,nvari
            read(lin,*,err=99)
      end do
c
c     Allocate as required:
c
      allocate(ixw(maxwd),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(iyw(maxwd),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(izw(maxwd),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(icatw(maxwd),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(imisw(maxwd),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
c
c     Now, read them for real:
c
      nwd   = 0
      nmiss = 0
 5    read(lin,*,end=6,err=99) (var(j),j=1,nvari)
      if(var(icf).lt.tmin.or.var(icf).ge.tmax) go to 5
      ifacies = var(icf)
      do ic=1,ncat
            if(ifacies.eq.cat(ic)) then
                  ind = ic
                  go to 7
            end if      
      end do
      go to 5
 7    continue      
      call getindx(nx,xmn,xsiz,var(icx),ii,inflag)
      if(.not.inflag) go to 5
      ix = max(min(ii,nx),1)
      call getindx(ny,ymn,ysiz,var(icy),ii,inflag)
      if(.not.inflag) go to 5
      iy = max(min(ii,ny),1)
      call getindx(nz,zmn,zsiz,var(icz),ii,inflag)
      if(.not.inflag) go to 5
      iz = max(min(ii,nz),1)
      nwd = nwd + 1
      ixw(nwd)   = ix
      iyw(nwd)   = iy
      izw(nwd)   = iz
      icatw(nwd) = ind
      imisw(nwd) = 0
      if(icatw(nwd).ne.sim(ix,iy,iz)) then
            imisw(nwd) = 1
            nmiss      = nmiss + 1
      end if
      go to 5
 6    continue
      close(lin)
      pmiss = real(nmiss)/real(nwd)*100
      write(*,101) nwd,nmiss,pmiss
 101  format(/,' Number of acceptable well data  = ',i8,
     +       /,' Number of them that mismatch    = ',i8,' =',f8.3,'%')
      if(nmiss.lt.1) then
            write(*,*) ' no well data to match!'
            stop
      end if      
c
c Calculate a spiral path:
c
      write(*,*)
      write(*,*) 'Calculating spiral path from data'
      npath = 0
      maxsize = 2*(nxt+1)*2*(nyt+1)*2*(nzt+1)
      allocate(ixpath(maxsize),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(iypath(maxsize),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(izpath(maxsize),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(wtpath(maxsize),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      do ixt=-nxt,nxt
      do iyt=-nyt,nyt
      do izt=-nzt,nzt
            xdis = real(ixt)/real(nxt)
            ydis = real(iyt)/real(nyt)
            zdis = real(izt)/real(nzt)
            dist = sqrt(xdis*xdis+ydis*ydis+zdis*zdis)
            if(dist.lt.1.0) then
                  npath = npath + 1
                  ixpath(npath) = ixt
                  iypath(npath) = iyt
                  izpath(npath) = izt
                  wtpath(npath) = (1.0-dist)**wexp
            end if
      end do
      end do
      end do
      write(*,102) npath
 102  format(/,' Number of locations on path  = ',i8)
c
c Setup MAPS weighting:
c
      write(*,*)
      write(*,*) 'Calculating MAPS weights'
      sumwts = 0
      nmaps  = 0
      maxsize = 2*(nxt+1)+2*(nyt+1)+2*(nzt+1)
      allocate(ixmaps(maxsize),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(iymaps(maxsize),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(izmaps(maxsize),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(wtmaps(maxsize),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      do ixm=-nxm,nxm
      do iym=-nym,nym
      do izm=-nzm,nzm
            xdis = real(ixm)/real(nxm)
            ydis = real(iym)/real(nym)
            zdis = real(izm)/real(nzm)
            dist = sqrt(xdis*xdis+ydis*ydis+zdis*zdis)
            if(dist.lt.1) then
                  nmaps = nmaps + 1
                  ixmaps(nmaps) = ixm
                  iymaps(nmaps) = iym
                  izmaps(nmaps) = izm
                  wtmaps(nmaps) = (1.2-dist)**1.5
                  sumwts = sumwts + wtmaps(nmaps)
            end if
      end do
      end do
      end do
      do i=1,nmaps
            wtmaps(i) = 0.8 * wtmaps(i) / sumwts
      end do
      write(*,103) nmaps
 103  format(/,' Number of locations in MAPS template = ',i8)
c
c Figure out a path:
c
      write(*,*)
      write(*,*) 'Figure out a path through the grid nodes'
      do iwd=1,nwd
            i1 = ixw(iwd)
            j1 = iyw(iwd)
            k1 = izw(iwd)
            do ip=1,npath
                  i2 = i1 + ixpath(ip)
                  j2 = j1 + iypath(ip)
                  k2 = k1 + izpath(ip)
                  if(i2.ge.1.and.i2.le.nx.and.
     +               j2.ge.1.and.j2.le.ny.and.
     +               k2.ge.1.and.k2.le.nz) then
                        if(cdat(i2,j2,k2).gt.0) then
                              if(wtpath(ip).gt.cwt(i2,j2,k2)) then
                                    cdat(i2,j2,k2) = iwd
                                    cwt(i2,j2,k2)  = wtpath(ip)
                              end if
                              go to 8
                        end if
                        cdat(i2,j2,k2) = iwd
                        cwt(i2,j2,k2)  = wtpath(ip)
 8                      continue
                  end if
            end do
      end do
c
c Unflag the cells that are close to data that match the well data:
c
      ntot = 0
      do k=1,nz
      do j=1,ny
      do i=1,nx
            if(cdat(i,j,k).ge.1) ntot = ntot + 1
      end do
      end do
      end do
      ncon = 0
      do k=1,nz
      do j=1,ny
      do i=1,nx
            if(cdat(i,j,k).ge.1) then
                  index = cdat(i,j,k)
                  if(imisw(index).eq.0) then
                        cdat(i,j,k) = -99
                        cwt(i,j,k)  = -99
                  end if      
            end if
            if(cdat(i,j,k).ge.1) ncon = ncon + 1
      end do
      end do
      end do
      write(*,104) ntot,ncon
 104  format(/,' Number of cells close to wells  = ',i8,
     +       /,' Number close to mismatch cells  = ',i8)
      if(ncon.lt.1) then
            write(*,*) ' no mismatches to fix!'
            go to 90
      end if
c
c Visit the locations close to the wells first:
c
      allocate(wnode(ncon),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(anode(ncon),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      allocate(adata(ncon),stat=test)
      if(test.ne.0) stop ' ERROR: memory allocation failure!'
      ind = 0
      do k=1,nz
      do j=1,ny
      do i=1,nx
            if(cdat(i,j,k).ge.1) then
                  index = i + (j-1)*nx+(k-1)*nx*ny
                  ind = ind + 1
                  anode(ind) = real(index)
                  adata(ind) = real(cdat(i,j,k))
                  wnode(ind) = -cwt(i,j,k)
            end if      
      end do
      end do
      end do
      if(ind.ne.ncon) write(*,*) ' ERROR: ',ind,ncon
      call sortem(1,ncon,wnode,2,anode,adata,d,e,f,g,h)
      nxy = ny * ny



c
c MAIN PART OF CODE: The MAPS Post Processor
c
      write(*,*)
      write(*,*) 'Running the MAPS Post Processing'
      do icon=1,ncon
c
c Get a location to consider:
c
      thewt = -wnode(icon)*((real(acorni(idum))*0.2)+0.9)
      index =  int(anode(icon)+0.5)
      iwdat =  int(adata(icon)+0.5)
      iz = int((index-1)/nxy) + 1
      iy = int((index-(iz-1)*nxy-1)/nx) + 1
      ix = index - (iz-1)*nxy - (iy-1)*nx
      if(thewt.ge.0.99) then
            sim(ix,iy,iz) = icatw(iwdat)
            go to 20
      end if
c
c Get MAPS weights:
c
      mpdf = 0.
      do im=1,nmaps
            i2 = ix + ixmaps(im)
            j2 = iy + iymaps(im)
            k2 = iz + izmaps(im)
            if(i2.ge.1.and.i2.le.nx.and.
     +         j2.ge.1.and.j2.le.ny.and.
     +         k2.ge.1.and.k2.le.nz) then
                  icat = sim(i2,j2,k2)
                  mpdf(icat) = mpdf(icat)+wtmaps(im)
            end if
      end do      
c
c Add in the well weight:
c
      icat = icatw(iwdat)
      mpdf(icat) = mpdf(icat) + thewt
c
c Draw a value to assign to this grid node:
c
      iold    = sim(ix,iy,iz)
      isimcat = 1
      pdfmax  = mpdf(1)
      do i=2,ncat
            if(mpdf(i).gt.pdfmax) then
                  isimcat = i
                  pdfmax  = mpdf(i)
            end if
      end do
      sim(ix,iy,iz) = isimcat
c
c Keep looping:
c
 20   continue
      end do



c
c Write final image:
c
 90   write(*,*)
      write(*,*) 'Writing output file'
      open(lout,file=outfl)
      write(lout,200)
 200  format('Output from MAPSpp',/,'1',/,'Category')
      do k=1,nz
      do j=1,ny
      do i=1,nx
            itest = sim(i,j,k)
            icat  = cat(itest)
            if(abs(icat).lt.9) then
                  write(lout,'(i2)') icat
            else
                  write(lout,'(i0)') icat
            end if
      end do
      end do
      end do
      close(lout)
c
c Finished:
c
      close(lout)
      write(*,9998) VERSION
 9998 format(/' MAPSpp Version: ',f5.3, ' Finished'/)
      stop
 97   stop 'ERROR in training image file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end



      double precision function acorni(idum)
c-----------------------------------------------------------------------
c
c Fortran implementation of ACORN random number generator of order less
c than or equal to 12 (higher orders can be obtained by increasing the
c parameter value MAXORD).
c
c
c NOTES: 1. The variable idum is a dummy variable. The common block
c           IACO is used to transfer data into the function.
c
c        2. Before the first call to ACORN the common block IACO must
c           be initialised by the user, as follows. The values of
c           variables in the common block must not subsequently be
c           changed by the user.
c
c             KORDEI - order of generator required ( must be =< MAXORD)
c
c             MAXINT - modulus for generator, must be chosen small
c                      enough that 2*MAXINT does not overflow
c
c             ixv(1) - seed for random number generator
c                      require 0 < ixv(1) < MAXINT
c
c             (ixv(I+1),I=1,KORDEI)
c                    - KORDEI initial values for generator
c                      require 0 =< ixv(I+1) < MAXINT
c
c        3. After initialisation, each call to ACORN generates a single
c           random number between 0 and 1.
c
c        4. An example of suitable values for parameters is
c
c             KORDEI   = 10
c             MAXINT   = 2**30
c             ixv(1)   = an odd integer in the (approximate) range 
c                        (0.001 * MAXINT) to (0.999 * MAXINT)
c             ixv(I+1) = 0, I=1,KORDEI
c
c
c
c Author: R.S.Wikramaratna,                           Date: October 1990
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      common/iaco/ ixv(MAXOP1)
      do i=1,KORDEI
            ixv(i+1)=(ixv(i+1)+ixv(i))
            if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
      end do
      acorni=dble(ixv(KORDEI+1))/MAXINT
      return
      end



      subroutine chknam(str,len)
c-----------------------------------------------------------------------
c
c                   Check for a Valid File Name
c                   ***************************
c
c This subroutine takes the character string "str" of length "len" and
c removes all leading blanks and blanks out all characters after the
c first blank found in the string (leading blanks are removed first).
c
c
c
c-----------------------------------------------------------------------
      parameter (MAXLEN=512)
      character str(MAXLEN)*1
c
c Find first two blanks and blank out remaining characters:
c
      do i=1,len-1
            if(str(i)  .eq.' '.and.
     +         str(i+1).eq.' ') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 2
            end if
      end do
 2    continue
c
c Look for "-fi" for file
c
      do i=1,len-2
            if(str(i)  .eq.'-'.and.
     +         str(i+1).eq.'f'.and.
     +         str(i+2).eq.'i') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 3
            end if
      end do
 3    continue
c
c Look for "\fi" for file
c
      do i=1,len-2
            if(str(i)  .eq.'\'.and.
     +         str(i+1).eq.'f'.and.
     +         str(i+2).eq.'i') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 4
            end if
      end do
 4    continue
c
c Return with modified file name:
c
      return
      end



      subroutine sortem(ib,ie,a,iperm,b,c,d,e,f,g,h)
c-----------------------------------------------------------------------
c
c                      Quickersort Subroutine
c                      **********************
c
c This is a subroutine for sorting a real array in ascending order. This
c is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
c in collected algorithms of the ACM.
c
c The method used is that of continually splitting the array into parts
c such that all elements of one part are less than all elements of the
c other, with a third part in the middle consisting of one element.  An
c element with value t is chosen arbitrarily (here we choose the middle
c element). i and j give the lower and upper limits of the segment being
c split.  After the split a value q will have been found such that 
c a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
c performs operations on the two segments (i,q-1) and (q+1,j) as follows
c The smaller segment is split and the position of the larger segment is
c stored in the lt and ut arrays.  If the segment to be split contains
c two or fewer elements, it is sorted and another segment is obtained
c from the lt and ut arrays.  When no more segments remain, the array
c is completely sorted.
c
c
c INPUT PARAMETERS:
c
c   ib,ie        start and end index of the array to be sorteda
c   a            array, a portion of which has to be sorted.
c   iperm        0 no other array is permuted.
c                1 array b is permuted according to array a
c                2 arrays b,c are permuted.
c                3 arrays b,c,d are permuted.
c                4 arrays b,c,d,e are permuted.
c                5 arrays b,c,d,e,f are permuted.
c                6 arrays b,c,d,e,f,g are permuted.
c                7 arrays b,c,d,e,f,g,h are permuted.
c               >7 no other array is permuted.
c
c   b,c,d,e,f,g,h  arrays to be permuted according to array a.
c
c OUTPUT PARAMETERS:
c
c    a      = the array, a portion of which has been sorted.
c
c    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
c
c NO EXTERNAL ROUTINES REQUIRED:
c
c-----------------------------------------------------------------------
      dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
c
c The dimensions for lt and ut have to be at least log (base 2) n
c
      integer   lt(64),ut(64),i,j,k,m,p,q
c
c Initialize:
c
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
c
c If this segment has more than two elements  we split it
c
 10   if (j-i-1) 100,90,15
c
c p is the position of an arbitrary element in the segment we choose the
c middle element. Under certain circumstances it may be advantageous
c to choose p at random.
c
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
c
c Start at the beginning of the segment, search for k such that a(k)>t
c
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
c
c Such an element has now been found now search for a q such that a(q)<t
c starting at the end of the segment.
c
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
c
c a(q) has now been found. we interchange a(q) and a(k)
c
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
c
c Update q and search for another pair to interchange:
c
      q = q-1
      go to 20
 50   q = k-1
 60   continue
c
c The upwards search has now met the downwards search:
c
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
c
c The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
c store the position of the largest segment in lt and ut
c
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
c
c Update m and split the new smaller segment
c
 80   m = m+1
      go to 10
c
c We arrive here if the segment has  two elements we test to see if
c the segment is properly ordered if not, we perform an interchange
c
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
c
c If lt and ut contain more segments to be sorted repeat process:
c
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      return
      end



      subroutine getindx(n,min,siz,loc,index,inflag)
c-----------------------------------------------------------------------
c
c     Gets the coordinate index location of a point within a grid
c     ***********************************************************
c
c
c n       number of "nodes" or "cells" in this coordinate direction
c min     origin at the center of the first cell
c siz     size of the cells
c loc     location of the point being considered
c index   output index within [1,n]
c inflag  true if the location is actually in the grid (false otherwise
c         e.g., if the location is outside then index will be set to
c         nearest boundary
c
c
c
c-----------------------------------------------------------------------
      integer   n,index
      real      min,siz,loc
      logical   inflag
c
c Compute the index of "loc":
c
      index = int( (loc-min)/siz + 1.5 )
c
c Check to see if in or out:
c
      if(index.lt.1) then
            index  = 1
            inflag = .false.
      else if(index.gt.n) then
            index  = n
            inflag = .false.
      else
            inflag = .true.
      end if
c
c Return to calling program:
c
      return
      end



      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='MAPSpp.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for MAPSpp',/,
     +       '                  *********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('Fluvsim/fluvsim01.dat          ',
     +       '-file with well conditioning data')
      write(lun,12)
 12   format('2  3  4  5                     ',
     +       '-  columns for X, Y, Z, facies')
      write(lun,13)
 13   format('-1.0       1.0e21              ',
     +       '-  trimming limits')
      write(lun,14)
 14   format('Fluvsim/fluvsim01t.out         ',
     +       '-file with inital image           ')          
      write(lun,15)
 15   format('1                              ',
     +       '-  column for categorical variable')                
      write(lun,16)
 16   format('1                              ',
     +       '-debugging level: 0,1,2,3,4')
      write(lun,17)
 17   format('MAPSpp.dbg                 ',
     +       '-file for debugging output')
      write(lun,18)
 18   format('MAPSpp.out                 ',
     +       '-file for simulation output')
      write(lun,19)
 19   format('1                              ',
     +       '-number of realizations')
      write(lun,20)
 20   format('100 25.0  50.0                 ',
     +       '-nx,xmn,xsiz')
      write(lun,21)
 21   format('100 25.0  50.0                 ',
     +       '-ny,ymn,ysiz')
      write(lun,22)
 22   format('1   0.005  0.01                ',
     +       '-nz,zmn,zsiz')
      write(lun,23)
 23   format('69069                          ',
     +       '-random number seed')
      write(lun,24)
 24   format('10                             ',
     +       '-number of loops')
      write(lun,25)
 25   format('10  10  4                      ',
     +       '-template size: nx, ny, nz')
      write(lun,26)
 26   format('1.0  0.25                      ',
     +       '-factor: multiple point and two point')



      close(lun)
      return
      end
