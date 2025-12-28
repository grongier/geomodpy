C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Developed in the Ph.D. Dissertation:                                 %
C                                                                      %
C The Integration of Geological Information into Geostatistical Models %
C                                                                      %
C Dr. Michael J. Pyrcz                                                 %
C Centre for Computational Geostatistics                               %
C University of Alberta                                                %
C                                                                      %
C Supervised by Dr. Clayton V. Deutsch                                 %
C Convocation Fall, 2004                                               %
C                                                                      %
C The programs in this collection are research code.  It is hoped that %
C they will be useful and that they may be adapted to a variety of     %
C applications.  See accompanying documention for a discussion of code %
C quality, limitations and ideas for future work.                      %
c                                                                      %
c Pyrcz, M.J. and Deutsch, C.V., in press, Alluvsim: a Conditional     %
c Event-based Fluvial Model, Computers & Geosciences.                  %
c                                                                      %
c Pyrcz, M.J., The Integration of Geologic Information into            %
c Geostatistical Models, Ph.D. Thesis, University of Alberta, 2004,    %
c 250 pages.                                                           %
c                                                                      %
c Appreciation to Fueng Zabel for her work to improve conditioning     %
c with a fitness match applied to new channels for early rejection,    %
c well conditioning summary and bug fixes in 2005.                     %
c                                                                      %
C No author or distributor accepts responsibility to anyone for the    %
C consequences of using them or for whether they serve any particular  %
C purpose or work at all, unless he says so in writing.  Everyone is   %
C granted permission to copy, modify and redistribute the programs in  %
C collection, but only under the condition that this original work is  %
C cited.                                                               %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Additional tips on using this program.
c
c 1.      Migration edge effects - the migration translation is set to zero 
c     at the source location and linearly increases to the calculated 
c     migration at 10% of the centerline length.  This was done to 
c     prevent the source location from moving or forming discontinuities 
c     (edge effects).  One way to deal with this would be to start the 
c     centerline outside the model.  These effects are more apparent if 
c     the model is set to migration 100% for many time steps.  
c     Infrequent major avulsions help mitigate these edge effects.
c 2.      Number of time steps and number of levels - the model is meant to 
c     run until the NTG is reached at each assigned level.  The maximum 
c     number of time steps is set to allow the user to stop the model 
c     early.  If no channels are placed vertically near the level of 
c     conditioning (within the user assigned z tolerance) then the 
c     conditioning will fail.  If too few levels are assigned to fill 
c     the model sufficiently to reach NTG the model will be stuck and 
c     run until maximum number of time steps.
c 3.      Well conditioning ease is improved by (1) having frequent 
c     avulsions as meander sets are moved together and this is difficult 
c     to condition, (2) making sure the trend models do not contradict 
c     the conditioning well data, (3) smaller elements (relative to well 
c     spacing), and (4) many events relative to the number of well 
c     intercepts to honor.  Conditioning only approximately honors net 
c     intervals - post processing with the included MAPSpp can help 
c     with exact reproduction of elements identified in wells.  Note 
c     the well report provide a mismatch percentage.
c 4.      Trend reproduction is soft.  Horizontal trends are "honored" by 
c     building a database of candidate centerlines and weighted drawing.  
c     The weights are the average trend value encountered along the 
c     path.  Therefore, a very specific trend model will not be 
c     reproduced (low local values can be averaged out by high local 
c     values elsewhere on the path).  Note: trend values closer to 0.0 
c     will make channel events at the location while larger values 
c     will encourage channels more likely.  Vertical trends are 
c     controlled by scaling the target NTG for each level.  Reasonable 
c     number of levels is required to reproduce the vertical trend.
c 5.      The NTG is honored within a single element volume.  When the target 
c     NTG is reached the model aggrades or terminates (if at the top 
c     level).  In this version, no attempt is made to scale the element 
c     sizes for precise NTG match.
c

      module streamsimMod

c-----------------------------------------------------------------------
c
c                 Alluvsim
c                 Event-based Fluvial Simulation
c                 ************************************************
c 
c 
c The event-based fluvial simulation is based on the streamline 
c building blocks. Alluvsim is tailored to fluvial models.  The building 
c blocks include elements:
c   
c   CH     - channel fills
c   LA     - lateral accretion
c   LV     - levees
c   CS     - crevasse channel and splay
c   FF(CH) - abandoned channel
c   FF     - overbank fines
c
c
c-----------------------------------------------------------------------
c
c The operators include:
c
c agradation - with simplified incision
c avulsion (within and proximal) - isolated and braiding streamlines  
c initialization - realistichannel streamline
c migration - realistimeander evolution
c cutoff - check for chute and neck cutoffs
c 
c The migration is calculated from the bank retreat fluvial model.  
c Based on procedures from:
c
c Howard, A.D., 1992, Lowland Floodplain Rivers: Geomorphological 
c Perspectives, chapter Modeling channel migration and floodplain 
c sedimentation in meandering streams, pages 1-37. John Whiley and Sons.
c 
c Sun, T., Meakin, P. and Jossang, T., 1996, A Simulation Model for 
c Meandering Rivers, Water Resources Research, Vol. 32, No. 9., pages 
c 2937-2954. 
c 
c The streamline initialization is based on a simplification of the 
c disturbed periodic model from:
c
c Ferguson, R.I., 1976, Disturbed PeriodiModel for River Meanders, 
c Earth Surface Processes, Vol. 1, pp. 337-347.
c
c AUTHOR: Michael J. Pyrcz                             DATE: 2002-2004
c-----------------------------------------------------------------------
c 
c External Subroutines:
c
c allocater       - allocate dynamiarrays
c avulsioninside  - avulse a streamline inside model
c azimuth              - calculate an azimuth for a line segment
c buildCHtable    - build streamline lookup table
c calc_levee           - calculate a levee geometry along a streamline
c calc_lobe            - calculate a lobe (SH, FS, PD) geometry along a 
c streamline           - may be applied for deepwater settings
c calcusb              - calculate the near bank velocity along a streamline
c copystream           - copy streamline
c curvature2           - calculate the streamline properties splines
c genabandonedchannel - claculate an abandoned channel geometry
c genchannel          - generate a channel geometry
c lookupstream        - draw a streamline from the streamline table
c migrate              - apply meander migration to a streamline
c morphCHendpts   - correct streamline to honor well intercept
c movwinsmooth    - apply moving window smoothing to a property array
c movwinsmoothcol - allpy moving window smoothing to a single array col.
c neckcutoff           - check for neck cutoffs and apply
c offset               - shift a node by an offset (distance and azimuth)
c onedrf              - calculate a 1D RF 
c probdraw             - draw from a weighted distribution
c repulseCHendpts - correct streamline to remove unwarranted intercept
c segmentwell     - format conditioning data
c
c-----------------------------------------------------------------------
c
c External GSLIB Subroutines:
c
c acorni          - random number generator
c chknam          - removes spaces from file names
c gauinv          - compute inverse standard normal distribution
c getinx          - calculate location index
c locate          - find location in an ascending array
c resc                 - linear interpolation
c 
c Deutsch, C.V. and Journel, A.G., GSLIB: Geostatistical Software 
c Library and Users Guide, Oxford University Press, New York, 1992,  
c 335 pages.  
c 
c----------------------------------------------------------------------- 
c
c External Numerical Recipes Subroutines (user must add these routines).
c 
c spline          - setup a cubic spline function
c splint          - spline interpolate
c 
c William H. Press, et al, Numerical Recipes in Fortran 90 : The Art of 
c Parallel ScientifiComputing (Fortran Numerical Recipes,  Vol 2), 
c Cambridge University Press, 1992, 571 pages.
c
c-----------------------------------------------------------------------

c
c Hydraulic Parameters
c
      
      real                     h0,us0,pb,Cf,scour_factor,gradient,g,
     +                          step,sinu,Q

c
c Geometric Parameters
c

      
      real      mCHdepth,stdevCHdepth
      real    mCHazi,stdevCHazi
      real    mCHwdratio,stdevCHwdratio
      real    mCHsource,stdevCHsource
      integer nCHdraw,ndiscr,nCHcor
      real    mLVdepth,stdevLVdepth,stdevCHdepth2
        real    mCHsinu,stdevCHsinu
      real    mLVwidth,stdevLVwidth
      real    mLVheight,stdevLVheight
      real    mLVasym,stdevLVasym
      real    mLVthin,stdevLVthin
      real    mCSlength,stdevCSlength
      real    mCSnum,stdevCSnum
      real    mCSnumlobe,stdevCSnumlobe
      real    mCSsource,stdevCSsource,CSsource
      real    mFFCHprop,stdevFFCHprop
      real    mdistMigrate,stdevdistMigrate
      real    mCSLOLL,stdevCSLOLL
      real    mCSLOWW,stdevCSLOWW
      real    mCSLOl,stdevCSLOl
      real    mCSLOw,stdevCSLOw
      real    mCSLO_hwratio,stdevCSLO_hwratio
      real    mCSLO_dwratio,stdevCSLO_dwratio

c
c Table parameters
c

      integer    ndis0,ntime,ndata_spline_table,ndata_chamat,
     +           ndata_chaprop,ndata_chadraw,ndata_chadrawprop,
     +           max_assoc,max_withinassoc ,ndata_wellheader,CHCSTEMP
      real,allocatable :: chamat(:,:,:),chaprop(:,:)

c
c Schedule Parameters
c

      real    probAvulOutside,probAvulInside
      integer CHcurrent,CHdraw,nlevel,nassoc
      real,allocatable :: level(:,:),schedule(:)

c
c NTG Parameters
c

      integer NTGcount,FFCHcount,NTGcurrent
      real    NTGtarget,actualNTG

c
c Facies Codes
c

      integer kFF,kFFCH,kCH,kLA,kLV,kCS
      real kFF_local,kFFCH_local,kCH_local,kLA_local,kLV_local,kCS_local
      real color_incr
      real,allocatable   :: facies(:,:,:)

c
c Grid and Files
c
  
        character                outfl*512,facfl*512,parfl*512,wellfl*512,
     +                           horifl*512,vertfl*512,ihsfl*512,fitfl*512,
     +                         upoutfl*512,upfacfl*512
      integer                  nx,ny,nz
      real                     xmn,ymn,zmn
      real                     xmin,ymin,zmin,xsiz,ysiz,zsiz,xmax,ymax,
     +                         zmax,ymid
      integer                  test,lout,lfac,lpar,lihs,lupout,lupfac,
     +                         lfit               
       
c
c IHS surfaces
c

      logical                  write_surfaces
      integer                  lsurf
      real,allocatable    ::   ihs(:,:,:,:)                
                       
c
c Temporary Arrays
c

      real,allocatable    ::   temp1dmat(:),testmat(:,:,:)                                      

c
c Horizontal and Vertical Trends 
c

      real,allocatable   ::    chadraw(:,:,:),chadrawprop(:,:)
      real,allocatable   ::    chadrawwt(:),horimat(:,:),vertmat(:)

c
c Spline arrays
c
      real,allocatable   ::    spline_table(:,:,:),spline_table2(:,:,:)
      real,allocatable   ::    spline_s(:)
      real,allocatable   ::    spline_d(:),spline_2d(:)
      real,allocatable   ::    spline_a(:),spline_2a(:)
      real,allocatable   ::    spline_x(:),spline_2x(:)
      real,allocatable   ::    spline_y(:),spline_2y(:)
      real,allocatable   ::    spline_w(:),spline_2w(:)
      real,allocatable   ::    spline_t(:),spline_2t(:)
      real,allocatable   ::    spline_c(:),spline_2c(:)
      real,allocatable   ::    spline_i(:),spline_2i(:)
      real,allocatable   ::    spline_z(:),spline_2z(:)
      real,allocatable   ::    spline_h(:),spline_2h(:)

c
c Data Conditioning 
c

      logical                  honor_well
      integer                  ndata,nwell,nCHdata,iwellcond    
      integer                  wcol,xcol,ycol,ztcol,zbcol,fcol
      integer                  node_buffer,totalViolate
      real                     xanis,yanis,zanis,ztol,buffer
      integer,allocatable::    wellpointer(:),welltemp(:),well_log(:) 
      real,allocatable::       percentViolate(:)  
      real,allocatable   ::    dataorder(:),assoc_tab(:,:)
      real,allocatable   ::    assoc_cond(:,:),well(:,:)
      real,allocatable   ::    welldata(:,:),wellheader(:,:)
      real,allocatable   ::    wellheader2(:,:),wellinterval(:,:)


c
c Program Constants
c

      parameter(PI=3.14159,VERSION=1.0,g=9.81,big_number=10.0e10)

c
c Random Number Seed
c

      common /iaco/   ixv(13)

c
c Module Finished
c

      end module

c-----------------------------------------------------------------------

      program main

c-----------------------------------------------------------------------

      use    streamsimMod
      implicit none

      integer itype
      logical rhonor_well
c
c Input/Output units used:
c
      
      lpar = 1
      lout = 2
      lfac = 3
      lfit = 4
      lihs = 7
      lupout = 8
      lupfac = 9
      
      kFF=-1
      kFFCH=0
      kCS=1
      kLV=2
      kLA=3
      kCH=4

c
c The total number of table parameters
c

      ndata_spline_table = 12
      ndata_chamat = 12
      ndata_chaprop = 23      
      ndata_chadraw = 3
      ndata_chadrawprop = 7
      ndata_wellheader = 57

c
c Read the parameters and trends
c
      
      call readparm

c
c Build a suite of potential channel splines (x,y,w)
c Generate candidate streamlines 
c
      call buildCHtable

c
c Check whether to perform well conditioning or not
c
      if(iwellcond.eq.1) then 

c To process conditioning data and prepare active arrays

          call segmentwell
      endif 
c
c Honor only soft data for unconditional run but also honor well data for conditioning run
c
      call alluvsim 

c
c Perform postprocessing to fine tune the match to all intervals
c
      call postprocessing

c
c Construct the model
c
      call constructmodelAll_noColorAdjust  

c
c Check the well intercepts
c
      call checkwell(rhonor_well)

c
c Construct the model again to adjust color
c
      call constructmodelAll


c Write out the model
c      
      itype = 0
      call writeout(itype)

c
c Finished:
c
      write(*,9998) VERSION
 9998 format(/' Alluvsim Version: ',f5.3, ' Finished'/)
      stop
      end

      subroutine readparm
c-----------------------------------------------------------------------
c
c                  Initialization and Read Parameters
c                  **********************************
c
c The input parameters and data are read in from their files. Some quick
c error checking is performed.
c
c AUTHOR: Michael J. Pyrcz                               DATE: 2002-2004
c-----------------------------------------------------------------------
      use       streamsimMod
      implicit none

      integer   i,ix,iy,iz,nvari,htcol,vtcol,lvert,lhori,lwell
      integer   ilevel,inflag,idum,idata
      real      var(50)
      real*8    p,acorni

      character str*512
      logical   testfl,testind

      lvert = 11
      lhori = 12
      lwell = 13

c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' Alluvsim Version: ',f5.3/)
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
      if(str(1:1).eq.' ') str(1:20) = 'alluvsim.par           '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'alluvsim.par           ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif

167   continue

      open(lpar,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lpar,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c

      read(lpar,'(a512)',err=98) wellfl 
      call chknam(wellfl,512)
      write(*,*) ' well file = ',wellfl(1:40) 

      read(lpar,*,err=98) wcol, xcol, ycol, ztcol, zbcol, fcol
      write(*,*) 'wcol, xcol, ycol, ztcol, zbcol, fcol = ',
     +      wcol, xcol, ycol, ztcol, zbcol, fcol

      read(lpar,*,err=98) xanis,yanis,zanis
      write(*,*) 'xanis,yanis,zanis = ',xanis,yanis,zanis

      read(lpar,*,err=98) buffer,ztol
      write(*,*) 'buffer,ztol = ',buffer, ztol

      read(lpar,'(a512)',err=98) horifl 
      call chknam(horifl,512)
      write(*,*) ' horizontal trend file = ',horifl(1:40) 

      read(lpar,*,err=98) htcol
      write(*,*) 'htcol = ',htcol

      read(lpar,'(a512)',err=98) vertfl 
      call chknam(vertfl,512)
      write(*,*) ' vertical trend file = ',vertfl(1:40) 

      read(lpar,*,err=98) vtcol
      write(*,*) 'vtcol = ',vtcol

      read(lpar,*,err=98) ntime      
      write(*,*) 'ntime = ',ntime

      max_assoc = ntime ! set to maximum possible to ensure 
      max_withinassoc = ntime ! used to size array to tract related groups of events

      read(lpar,*,err=98) nlevel                  
   
      allocate(level(nlevel,2),stat = test)
      if(test.ne.0)then
       write(*,*)'Level Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if

      backspace(lpar)

      read(lpar,*,err=98) ,nlevel,(level(i,1),i=1,nlevel) 
      write(*,*) ' nlevel,(level(i,1),i=1,nlevel) = ',
     +               nlevel,(level(i,1),i=1,nlevel) 

      read(lpar,*,err=98) NTGtarget,Cf,scour_factor,gradient,Q
      write(*,*) ' NTGtarget,Cf,scour_factor,gradient,Q = ',
     +               NTGtarget,Cf,scour_factor,gradient,Q  

      read(lpar,*,err=98) nCHdraw,ndiscr,nCHcor
      write(*,*) ' nCHdraw,ndiscr,nCHcor = ',nCHdraw,ndiscr,nCHcor  

      read(lpar,*,err=98) probAvulOutside,probAvulInside
      write(*,*) ' probAvulOutside,probAvulInside = ',
     +               probAvulOutside,probAvulInside
      if(probAvulOutside+probAvulInside.gt.1.00001) then
       stop 'Sum of Avulsion Probabilities Greater than 1.0'
      end if
     
      read(lpar,*,err=98) mCHazi,stdevCHazi
      write(*,*) ' mCHazi,stdevCHazi = ',mCHazi,stdevCHazi

      read(lpar,*,err=98) mCHsource,stdevCHsource
      write(*,*) ' mCHsource,stdevCHsource = ',mCHsource,stdevCHsource

      if(stdevCHsource.lt.0.0) then 
       write(*,*) ' Appling Uniform Source Distribution' 
      end if

      read(lpar,*,err=98) mCHdepth,stdevCHdepth,stdevCHdepth2
      write(*,*) ' mCHdepth,stdevCHdepth,stdevCHdepth2 = ',
     +               mCHdepth,stdevCHdepth,stdevCHdepth2

      read(lpar,*,err=98) mCHwdratio,stdevCHwdratio
      write(*,*) ' mCHwdratio,stdevCHwdratio = ',
     +               mCHwdratio,stdevCHwdratio 

      read(lpar,*,err=98) mCHsinu,stdevCHsinu
      write(*,*) ' mCHsinu,stdevCHsinu = ',mCHsinu,stdevCHsinu 

      read(lpar,*,err=98) mLVdepth,stdevLVdepth
      write(*,*) ' mLVdepth,stdevLVdepth = ',mLVdepth,stdevLVdepth 

      read(lpar,*,err=98) mLVwidth,stdevLVwidth
      write(*,*) ' mLVwidth,stdevLVwidth = ',mLVwidth,stdevLVwidth 

      read(lpar,*,err=98) mLVheight,stdevLVheight
      write(*,*) ' mLVheight,stdevLVheight = ',mLVheight,stdevLVheight 

      read(lpar,*,err=98) mLVasym,stdevLVasym 
      write(*,*) ' mLVasym,stdevLVasym = ',mLVasym,stdevLVasym 

      read(lpar,*,err=98) mLVthin,stdevLVthin 
      write(*,*) ' mLVthin,stdevLVthin = ',mLVthin,stdevLVthin 

      read(lpar,*,err=98) mCSnum,stdevCSnum
      write(*,*) ' mCSnum,stdevCSnum = ',mCSnum,stdevCSnum

      read(lpar,*,err=98) mCSnumlobe,stdevCSnumlobe
      write(*,*) ' mCSnumlobe,stdevCSnumlobe = ',
     +             mCSnumlobe,stdevCSnumlobe

      read(lpar,*,err=98) mCSsource,stdevCSsource
      write(*,*) ' mCSsource,stdevCSsource = ',mCSsource,stdevCSsource

      read(lpar,*,err=98) mCSLOLL,stdevCSLOLL
      write(*,*) ' mCSLOLL,stdevCSLOLL = ',mCSLOLL,stdevCSLOLL

      read(lpar,*,err=98) mCSLOWW,stdevCSLOWW
      write(*,*) ' mCSLOWW,stdevCSLOWW = ',mCSLOWW,stdevCSLOWW

      read(lpar,*,err=98) mCSLOl,stdevCSLOl
      write(*,*) ' mCSLOl,stdevCSLOl = ',mCSLOl,stdevCSLOl

      read(lpar,*,err=98) mCSLOw,stdevCSLOw
      write(*,*) ' mCSLOw,stdevCSLOw = ',mCSLOw,stdevCSLOw

      read(lpar,*,err=98) mCSLO_hwratio,stdevCSLO_hwratio
      write(*,*) ' mCSLO_hwratio,stdevCSLO_hwratio = ',
     +               mCSLO_hwratio,stdevCSLO_hwratio
      
      read(lpar,*,err=98) mCSLO_dwratio,stdevCSLO_dwratio
      write(*,*) ' mCSLO_dwratio,stdevCSLO_dwratio = ',
     +               mCSLO_dwratio,stdevCSLO_dwratio

      read(lpar,*,err=98) mFFCHprop,stdevFFCHprop
      write(*,*) ' mFFCHprop,stdevFFCHprop = ',mFFCHprop,stdevFFCHprop

      read(lpar,*,err=98) mdistMigrate,stdevdistMigrate
      write(*,*) ' mdistMigrate,stdevdistMigrate = ',
     +          mdistMigrate,stdevdistMigrate

      read(lpar,*,err=98) nx,xmn,xsiz
      write(*,*) ' X grid specification = ',nx,xmn,xsiz

      read(lpar,*,err=98) ny,ymn,ysiz
      write(*,*) ' Y grid specification = ',ny,ymn,ysiz

      read(lpar,*,err=98) nz,zmn,zsiz
      write(*,*) ' Z grid specification = ',nz,zmn,zsiz

      read(lpar,*,err=98) ixv(1),color_incr
      write(*,*) ' Random number seed, color_incr = ',ixv(1),color_incr

      do i=1,1000
             p = acorni(idum)
      end do

      read(lpar,'(a512)',err=98) facfl 
      call chknam(facfl,512)
      write(*,*) ' facies output file = ',facfl(1:40) 
      
      read(lpar,'(a512)',err=98) outfl 
      call chknam(outfl,512)
      write(*,*) ' streamline output file = ',outfl(1:40)    
      
      read(lpar,'(a512)',err=98) fitfl 
      call chknam(fitfl,512)
      write(*,*) ' measure of fitness output file = ',fitfl(1:40)          
      close(lpar)

c
c IHS accretionary surfaces
c

      write_surfaces = .false.

c
c Read in the well data
c
      iwellcond = 1            
      inquire(file=wellfl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR file ',wellfl(1:40),' does not exist!'
            write(*,*) '  no well file, assuming no well conditioning'   
              iwellcond = 0 
            nwell = 0
            go to 50
      end if
      open(lwell,file=wellfl,status='OLD') 
      read(lwell,*,err=99)
      read(lwell,*,err=99) nvari
      do i=1,nvari
            read(lwell,*,err=99)
      end do
      ndata=0
301   read(lwell,*,end=302,err=99)
      ndata=ndata+1
      goto 301
302   continue
c
      allocate(well(ndata,6),stat = test)
      if(test.ne.0)then
       write(*,*)'Well Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      well = 0
c
      allocate(wellpointer(ndata),stat = test)
      if(test.ne.0)then
       write(*,*)'WellPointer Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(welltemp(ndata),stat = test)
      if(test.ne.0)then
       write(*,*)'WellTemp Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      rewind(lwell)
      read(lwell,*,err=99)
      read(lwell,*,err=99) nvari
      do i=1,nvari
            read(lwell,*,err=99)
      end do


      do idata = 1, ndata
         read(lwell,*,err=99) (var(i),i=1,nvari)
            well(idata,1)=var(wcol)
       well(idata,2)=var(xcol)
       well(idata,3)=var(ycol)
            well(idata,4)=var(ztcol)
            well(idata,5)=var(zbcol)
       well(idata,6)=var(fcol)
      end do
      close(lwell)
  50  continue

c
c Read in the horizontal trend
c

      allocate(horimat(nx,ny),stat = test)
      if(test.ne.0)then
       write(*,*)'Horizontal Trend Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      horimat = 1.0
c
      inquire(file=horifl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR file ',horifl(1:40),' does not exist!'
            write(*,*) '  assuming no horizontal trend'
            do iy = 1, ny
             do ix = 1, nx
              horimat(ix,iy) = 1.0
             end do
            end do 
            goto 51
      end if
      open(lhori,file=horifl,status='OLD') 
      read(lhori,*,err=95)
      read(lhori,*,err=95) nvari
      do i=1,nvari
            read(lhori,*,err=95)
      end do
      do iy = 1, ny
       do ix = 1, nx
        read(lhori,*,err=95) (var(i),i=1,nvari)
        horimat(ix,iy) = var(htcol)
       end do
      end do
      close(lhori)
51    continue
      

c
c Read in the vertical trend
c

      allocate(vertmat(nz),stat = test)
      if(test.ne.0)then
       write(*,*)'Vertical Trend Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      vertmat = 0.0
c
      inquire(file=vertfl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR file ',vertfl(1:40),' does not exist!'
            write(*,*) '  assuming no vertical trend'
            do iz = 1, nz
             vertmat(iz)=1.0
            end do
            goto 52
      end if
      open(lvert,file=vertfl,status='OLD') 
      read(lvert,*,err=96)
      read(lvert,*,err=96) nvari
      do i=1,nvari
            read(lvert,*,err=96)
      end do
      do iz = 1, nz
       read(lvert,*,err=96) (var(i),i=1,nvari)
       vertmat(iz)=var(vtcol)
      end do
      close(lvert)
52    continue

c
c Convert vertical trend to cumulative
c

      do iz = 2, nz
       vertmat(iz)=vertmat(iz)+vertmat(iz-1)
      end do
      
c
c Pick the cumulatives for each level
c

      do ilevel = 1, nlevel
       call getindx(nz,zmn,zsiz,level(ilevel,1),iz,inflag)
       level(ilevel,2) = vertmat(iz)
      end do

c
c Calculate the NTG for each level (missing top is divided equally)
c

      NTGcount = int(NTGtarget*real(nx*ny*nz))     
      do ilevel = 1, nlevel
       level(ilevel,2) = real(NTGcount)*level(ilevel,2)/level(nlevel,2)
      end do      

c
c Find the needed parameters:
c
       
      step = min(xsiz,ysiz)

      us0 = ((g*Q*gradient)/(mCHdepth*mCHwdratio*Cf))**(1.0/3.0) ! average
      h0 = Q/(mCHdepth*mCHwdratio*us0)

      xmin = xmn-0.5*xsiz
      ymin = ymn-0.5*ysiz
      zmin = zmn-0.5*zsiz

      xmax = xmin+xsiz*real(nx)
      ymax = ymin+ysiz*real(ny)
      zmax = zmin+zsiz*real(nz)

      ymid = (ymax+ymin)/2.0

      step = (xsiz+ysiz)/2.0
      ndis0 = (((xmax-xmin)+(ymax-ymin))/2.0)/step*2

      node_buffer = int(real(buffer)/real(step))

c
c Checks
c
 
      if(mCHsource.gt.ymax.or.mCHsource.lt.ymin) then
       write(*,*)'mCHsource should not be set outside the model 
     +  (>zmax or <zmin).'
       stop
      end if

c
c Allocate arrays
c

      call allocater ! all arguments are in streamsimMOD module 

c
c Return to the main program:
c
      return

c
c Error in an Input File Somewhere:
c

 95   stop 'ERROR in horizontal trend file!'
 96   stop 'ERROR in vertical trend file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in well data file!'
      end



      subroutine alluvsim
c-----------------------------------------------------------------------
c
c           Generate Events
c           ***************
c
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c Modifications: Fueng Zabel                             DATE: 2004-2005 
c-----------------------------------------------------------------------
      use       streamsimMod
      implicit none

c
c General Variables
c

      integer   idum,ierr,ilevel,idis,ndis,nstart,avulsion_node,regrid,
     +           extrapolate,wellNum,prevWellNum,iSum,iCount
      real      color_current,val
      real      CHazi,CHsinu,CHdepth,distMigrate,CHmaxhalfwidth
      real      CSnum,CSnumlobe,CSsinu
      integer   iCSnum,CSnode,iCSnumlobe
      integer   CHlevel,iblklevel,inflagilevel,iCHdata,iblkCHinterval,
     +          nCHtoMatch,iCHtoMatch(10),iMatch,iwell(10),
     +          maxCHdrawn,rCase
      real      ztop,zbot,t1,t2
      real      curvature,CSazi,CSLOLL,CSLOWW,CSLOw,CSLOl
      real      CSLO_dwratio,CSLO_hwratio,xp,cutoff_tol
      real*8    p,acorni
      logical   cutoff,first_on_level,iFound,iPostprocess

c
c Streamline Parameters
c

      real   CHdwratio,CHwdratio,CHlength
      real   LVwidth,LVheight,LVdepth,LVasym,LVthin
      real   FFCHprop

c
c Assume the temporary CS streamline storage location
c

      CHCStemp = ntime

c
c Loop over the levels
c --------------------
 
      kFFCH_local = kFFCH
      kFF_local = kFF
      kCS_local = kCS
      kLV_local = kLV
      kLA_local = kLA
      kCH_local = kCH
      color_current = 0.0

      NTGcurrent = 0 ! count of net cells
      FFCHcount = 0  ! count of FF(CH) cells in last channel
      CHcurrent = 1  
      nassoc = 0
      actualNTG = 0 

      do ilevel = 1, nlevel
       first_on_level = .true.
       write(*,*) 'Working on Level #', ilevel

c Prepare extra information for well conditioning case
c -----------------------------------------------------

       if(iwellcond.eq.1) then ! for well conditioning
           nCHtoMatch=0
           iCount=1
           iSum = int(wellheader(iCount,6))
c
c Check if there is well interval to match for this level 
c 
             call getindx(nz,zmn,zsiz,level(ilevel,1),
     +       iblklevel,inflagilevel) ! block number for this ilevel
           
           do iCHdata=1,nCHdata
c
c Count net intervals to match for this level
                  
             ztop = welldata(iCHdata,4)
               zbot = welldata(iCHdata,5)
               iblkCHinterval = welldata(iCHdata,8) 

c Count number of wells.
             if(iCHdata.gt.iSum)then
                iCount = iCount + 1
                iSum = iSum + int(wellheader(iCount,6))
             end if

               if(abs(iblklevel-iblkCHinterval).le.1)then
               nCHtoMatch = nCHtoMatch + 1
               iCHtoMatch(nCHtoMatch) = iCHdata
               iwell(nCHtoMatch) = iCount
               end if           
           end do
        end if 
c
c Draw an initial channel
c
      
       call probdraw(chadrawwt,nCHdraw,CHdraw,val)              
       call lookupstream(ilevel,CHdraw,CHcurrent)
       chaprop(21,CHcurrent) = nassoc                 
       regrid = 1
       extrapolate = 1
       call curvature2(CHcurrent,extrapolate,regrid)

c
c Generate streamlines until reaching layer NTG 
c ---------------------------------------------
c

       do while(real(NTGcurrent-FFCHcount).lt.level(ilevel,2)) 
         write(*,*) '    Event #', CHcurrent


        if(CHcurrent.eq.78) then
          continue
        end if

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        FFCHprop = min(1.0,max(0.0,xp*stdevFFCHprop+mFFCHprop))
        chaprop(13,CHcurrent) = FFCHprop

c
c Decide to avulse or migrate  
c

        p = real(acorni(idum))


c
c Check for neck cutoff
c
        if(.not.first_on_level) then
         CHmaxhalfwidth = chaprop(17,CHcurrent)
         cutoff_tol= CHmaxhalfwidth*1.0
         cutoff=.false.
         call neckcutoff(cutoff_tol,cutoff,CHcurrent-1,CHcurrent)
         if(cutoff) then
          assoc_tab(nassoc,1) = assoc_tab(nassoc,1) + 1 ! update association tab
          assoc_tab(nassoc,assoc_tab(nassoc,1)+1) = CHcurrent 
          chaprop(21,CHcurrent) = nassoc  
          regrid = 1
          extrapolate = 1
          call curvature2(CHcurrent,extrapolate,regrid)
         end if
        end if

ccccchamat
c
c CASE 0:Neck cut off in the last step
c
cccc

        if(cutoff.and.(CHcurrent.gt.1)) then
         chaprop(22,CHcurrent-1) = kFFCH_local
         cutoff =.false.

c
c Include the cut off and oxbow lake as separate streamlines       
c


ccccchamat
c
c EVENT CASE 1:Avulse Proximal of Model
c
cccc

        else if(p.lt.probAvulOutside.or.first_on_level) then
         write(*,*) 'Operation: Proximal Avulsion'
         write(*,*) '  Element: FF(CH)'
         FFCHprop = chaprop(13,CHcurrent)
         if(.not.first_on_level.and.CHcurrent.ne.1) then 
          call genabandonedchannel(FFCHprop,CHcurrent-1,kFFCH_local,
     +      kCH_local)
          chaprop(22,CHcurrent-1) = kFFCH_local    
         end if
         

c
c Draw a new streamline
c

         if(iwellcond.eq.1) then
c            
c For well conditioning, select streamline based on Accepting/Rejecting rules. 
c Otherwise just simply draw streamline if there is no well data.
c         
              call selectstreamline(ilevel,nCHtoMatch,iwell,
     +                              iCHtoMatch,rCase)
                  
         else ! Alluvsim with no well conditioning
             call probdraw(chadrawwt,nCHdraw,CHdraw,val)       
             call lookupstream(ilevel,CHdraw,CHcurrent)             
         endif ! end well conditioning check

             nassoc = nassoc + 1
             assoc_tab(nassoc,1) = 1
             assoc_tab(nassoc,2) = CHcurrent
             chaprop(21,CHcurrent) = nassoc
             regrid = 1
             extrapolate = 1
             call curvature2(CHcurrent,extrapolate,regrid)
        
cccc
c
c EVENT CASE 2: Avulse Inside the Model
c
cccc

        else if(p.lt.probAvulOutside+probAvulInside) then     
         write(*,*) 'Operation: Avulsion'
         write(*,*) '  Element: FF(CH)'
         FFCHprop = chaprop(13,CHcurrent)
         call genabandonedchannel(FFCHprop,CHcurrent-1,kFFCH_local,
     +    kCH_local)
         chaprop(22,CHcurrent-1) = kFFCH_local
         CHsinu = chaprop(4,CHcurrent)
         CHazi = chaprop(1,CHcurrent)
         avulsion_node = -1 ! initialize as null - subroutine draws a node
         CHlength = -1.0
         call avulsioninside(avulsion_node,CHsinu,CHazi,CHlength,
     +    CHcurrent-1,CHcurrent)
         assoc_tab(nassoc,1) = assoc_tab(nassoc,1) + 1 ! update association tab
         assoc_tab(nassoc,assoc_tab(nassoc,1)+1) = CHcurrent  
         chaprop(21,CHcurrent) = nassoc    
         regrid = 1
         extrapolate = 1
         call curvature2(CHcurrent,extrapolate,regrid)  

cccc
c
c EVENT CASE 3: Migrate the Channel
c
cccc

        else 
         write(*,*) 'Operation: Migration'
         p = real(acorni(idum))
         call gauinv(p,xp,ierr)
         distMigrate = max(0.0,xp*stdevdistMigrate+mdistMigrate)
         chaprop(5,CHcurrent) = LVdepth
         write(*,*) '  Element: LA'
         call genchannel(CHcurrent-1,kLA_local) ! place LA before migration
         chaprop(22,CHcurrent-1) = kLA_local
         call calcusb(CHcurrent-1)
         call migrate(distMigrate,CHcurrent-1,CHcurrent)
         assoc_tab(nassoc,1) = assoc_tab(nassoc,1) + 1 ! update association tab
         assoc_tab(nassoc,assoc_tab(nassoc,1)+1) = CHcurrent
         chaprop(21,CHcurrent) = nassoc  
         regrid = 1
         extrapolate = 1
         call curvature2(CHcurrent,extrapolate,regrid)
        end if ! end event check
chamat
c Draw the architectural element geometric parameters
c ---------------------------------------------------

c
c NOTE: CHazi and CHsinu are applied in the CH table construction
c

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        CHdepth = max(0.0,xp*stdevCHdepth+mCHdepth)
        chaprop(2,CHcurrent) = CHdepth

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        CHwdratio = max(0.0,xp*stdevCHwdratio+mCHwdratio)
        CHdwratio = 1.0/CHwdratio
        chaprop(3,CHcurrent) = CHdwratio ! store CH depth:width ratio

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        LVdepth = max(0.0,xp*stdevLVdepth+mLVdepth)
        chaprop(5,CHcurrent) = LVdepth

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        LVwidth = max(0.0,xp*stdevLVwidth+mLVwidth)
        chaprop(6,CHcurrent) = LVwidth

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        LVheight = max(0.0,xp*stdevLVheight+mLVheight)
        chaprop(7,CHcurrent) = LVheight

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        LVasym = max(0.0,xp*stdevLVasym+mLVasym)
        chaprop(8,CHcurrent) = LVasym

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        LVthin = max(0.0,xp*stdevLVthin+mLVthin)
        chaprop(9,CHcurrent) = LVthin

c CSlength is drawn for each random walker

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        CSnum = max(0.0,xp*stdevCSnum+mCSnum)
        chaprop(11,CHcurrent) = CSnum

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        CSnumlobe = max(0.0,xp*stdevCSnumlobe+
     +   mCSnumlobe)
        chaprop(14,CHcurrent) = CSnumlobe

        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        CSsource = max(0.001,xp*stdevCSsource+CSsource)
        chaprop(12,CHcurrent) = CSsource
c
c Generate the architectural elements
c -----------------------------------
c
          call constructmodel(CHcurrent)  
c
c Copy the current streamline to the next
c ---------------------------------------
c

        CHcurrent = CHcurrent+1
        if(CHcurrent.eq.ntime-1) then
         write(*,*) 'WARNING: Reached ntime before NTG reached'
         actualNTG = actualNTG + (NTGcurrent-FFCHcount)
         actualNTG = real(actualNTG/(nx*ny*nz))
         goto 911
        end if       
        call copystream(CHcurrent-1,CHcurrent)
        regrid = 1
        extrapolate = 1
        call curvature2(CHcurrent,extrapolate,regrid)
        first_on_level=.false.
        continue
       end do ! end while loop NTG check for a level

c
c Sum actual NTG for all levels
c      
       actualNTG = actualNTG + (NTGcurrent-FFCHcount) 

      end do  ! end loop over all levels

      actualNTG = real(actualNTG/(nx*ny*nz))

911   continue
      return
      end

      subroutine writeout(rtype)
c-----------------------------------------------------------------------
c
c                  Write out the streamlines and Facies Model
c                  ******************************************
c 
c  Write out the streamlines and Facies Model
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      use       streamsimMod
      implicit none

c
c Local Variables
c

      integer   ix,iy,iz,idis,iCH,itime,lfil
      real      x,y,z,color,i,c,d

c
c Looked Up Variables
c

      integer   ndis

c
c Passed Variables
c

      integer   rtype

c
c Write out the facies model
c

      if(rtype.eq.0) then
       lfil=lfac
       write(*,*) ' Writing out facies model '
       open(lfil,file=facfl,status='UNKNOWN')
       write(lfil,*) 'ALLUVSIM_Facies_Output'
      else
       lfil=lupfac
       write(*,*) ' Writing out updated facies model '
       open(lfil,file=upfacfl,status='UNKNOWN')
       write(lfil,*) 'ALLUVSIM_Facies_Output'
      end if
      write(lfil,*) 1
      write(lfil,*) 'Facies'
      do iz = 1, nz
       do iy = 1, ny
        do ix = 1, nx
         write(lfil,711) facies(ix,iy,iz)
        end do
       end do
      end do
      close(lfil)
      
711   format(f5.2)

c
c Write out the streamlines
c

      if(rtype.eq.0) then
       lfil=lout
       write(*,*) ' Writing out event centerlines '
       open(lfil,file=outfl,status='UNKNOWN')
       write(lfil,*) 'ALLUVSIM_Event_Centerlines_Output'
      else
       lfil=lupout
       write(*,*) ' Writing out update event centerlines '
       open(lfil,file=upoutfl,status='UNKNOWN')
       write(lfil,*) 'ALLUVSIM_Updated_Event_Centerlines_Output'     
      end if
      write(lfil,*) 8
      write(lfil,*) 'streamline_number'
      write(lfil,*) 'x'
      write(lfil,*) 'y'
      write(lfil,*) 'z'
      write(lfil,*) 'facies'
      write(lfil,*) 'azi'
      write(lfil,*) 'curvature'
      write(lfil,*) 'dcurvature'

      do iCH = 1, CHcurrent-1
       ndis = chaprop(15,iCH)
       color = chaprop(22,iCH)
       do idis = 1, ndis
        x = chamat(idis,1,iCH)
        y = chamat(idis,2,iCH)
        z = chamat(idis,11,iCH)
        i = chamat(idis,9,iCH)
        c = chamat(idis,6,iCH)
        d = chamat(idis,8,iCH)
        write(lfil,712) iCH,x,y,z,color,i,c,d
       end do
      end do
      close(lout)

712   format(i5,1x,7(f12.5,1x))

c
c Write out the measure of fitness
c
      lfil=lfit
      write(*,*) ' Writing out measure of fitness '
      open(lfil,file=fitfl,status='UNKNOWN')
      write(lfil,*) 'ALLUVSIM_Fitness_Output'      
      write(lfil,*) 3
      write(lfil,*) 'random number seed'
       write(lfil,*) 'actual NTG'
      write(lfil,*) 'percent well violation'

      write(lfil,713) ixv(1),actualNTG,totalViolate
      close(lfil)

713   format(i9,' ',f5.2,' ',i5)

c
c Return to the main program:
c

      continue
      return
      end

      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='alluvsim.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for ALLUVSIM',/,
     +       '                  ***********************',/,/,
     +       'START OF PARAMETERS:')
      write(lun,31)
 31   format('welldata.dat                        ',
     +  '-file with well data')
      write(lun,32)
 32   format('1 2 3 4 7 9                       ',
     +  '- wcol,xcol,ycol,ztcol,zbcol,fcol')
      write(lun,33)
 33   format('1.0 1.0 1.0                         ',
     +  '-xanis,yanis,zanis')
      write(lun,34)
 34   format('10 3.0                              ',
     +  '- buffer, ztol')
      write(lun,11)
 11   format('horitrend.dat                       ',
     +  '-file with the horizontal trend')
      write(lun,12)
 12   format('1                                   ',
     +  '- htcol')
      write(lun,13)
 13   format('verttrend.dat                       ',
     +  '-file with the vertical trend')
      write(lun,14)
 14   format('1                                   ',
     +  '- vtcol')
      write(lun,85)
 85   format('100                                 ',
     +  '- ntime')      
      write(lun,15)
 15   format('3 7.0 13.0 17.0                     ',
     +  '- nlevel, level elevations')
      write(lun,16)
 16   format('0.05 0.0036 10.0 0.001 0.50         ',
     +  '- NTGtarget,Cf,Scour_Factor,gradient,Q')
      write(lun,17)
 17   format('100 10 10                           ',
     +  '- CHndraw,ndiscr,nCHcor')
      write(lun,18)
 18   format('0.3 0.3                             ',
     +  '- probAvulOutside,probAvulInside')
      write(lun,19)
 19   format('90.0 1.0                            ',
     +  '- CH element: mCHazi,stdevCHazi')
      write(lun,61)
 61   format('500.0 50.0                          ',
     +  '-  mCHsource,stdevCHsource')
      write(lun,20)
 20   format('2.0 0.5 0.2                         ',
     +  '-  mCHdepth,stdevCHdepth,stdevCHdepth2')
      write(lun,21)
 21   format('10.0 2.0                            ',
     +  '-  mCHwdratio,stdevCHwdratio')
      write(lun,22)
 22   format('1.3 0.2                             ',
     +  '-  mCHsinu,stdevCHsinu')
      write(lun,23)
 23   format('1.0 0.1                             ',
     +  '- LV Element: mLVdepth,stdevLVdepth')
      write(lun,24)
 24   format('40.0 5.0                            ',
     +  '-  mLVwidth,stdevLVwidth')
      write(lun,25)
 25   format('0.5 0.2                             ',
     +  '-  mLVheight,stdevLVheight')
      write(lun,26)
 26   format('0.3 0.1                             ',
     +  '-  mLVasym,stdevLVasym')
      write(lun,27)
 27   format('0.3 0.1                             ',
     +  '-  mLVthin,stdevLVthin')
      write(lun,29)
 29   format('2 2                                 ',
     +  '- CS Element: mCSnum,stdevCSnum')
      write(lun,30)
 30   format('3 2                                 ',
     +  '-  mCSnumlobe,stdevCSnumlobe')
      write(lun,81)
 81   format('50.0 20.0                           ',
     +  '-  mCSsource,stdevCSsource')
      write(lun,82)
 82   format('200.0 50.0                          ',
     +  '-  mCSLOLL,stdevCSLOLL')
      write(lun,83)
 83   format('30.0 10.0                           ',
     +  '-  mCSLOWW,stdevCSLOWW')
      write(lun,84)
 84   format('100.0 20.0                          ',
     +  '-  mCSLOl,stdevCSLOl')
      write(lun,35)
 35   format('20.0 10.0                           ',
     +  '-  mCSLOw,stdevCSLOw')
      write(lun,36)
 36   format('0.03 0.05                           ',
     +  '-  mCSLO_hwratio,stdevCSLO_hwratio')
      write(lun,37)
 37   format('0.02 0.05                           ',
     +  '-  mCSLO_dwratio,stdevCSLO_dwratio')
      write(lun,38)
 38   format('0.3  0.1                            ',
     +  '- FFCH Element: mFFCHprop,stdevFFCHprop')
      write(lun,39)
 39   format('50.0 20.0                           ',
     +  '- mdistMigrate,stdevdistMigrate')
      write(lun,40)
 40   format('100 5.0 10.0                        ',
     +  '-nx,xmn,xsiz')
      write(lun,41)
 41   format('100 5.0 10.0                        ',
     +  '-ny,ymn,ysiz')
      write(lun,42)
 42   format('20 0.5 1.0                          ',
     +  '-nz,zmn,zsiz')
      write(lun,43)
 43   format('69069 0.05                          ',
     +  '-random number seed, color_incr')
      write(lun,44)
 44   format('alluvsim.out                       ',
     +  '-file for output facies file')
      write(lun,45)
 45   format('alluvline.out                      ',
     +  '-file for output updated streamlines')
      write(lun,46)
 46   format('alluvsim_condreport.out                   ',
     +  '-file for conditioning report')
      write(lun,47)
 47   format('alluvlineupdate.out                  ',
     +  '-file for output updated streamlines')
 
      close(lun)
      return
      end


c-----------------------------------------------------------------------
c
c     Program Specific Subroutines 
c-----------------------------------------------------------------------


      subroutine allocater
    
c-----------------------------------------------------------------------
c
c     Generate the channel facies model
c     *********************************
c
c Generate a realistic channel cross section based on Deutsch and Wang.
c Reuires a channel spline.
c
c INPUT VARIABLES:
c
c   ndis             the current number of streamline discretizations
c   chamat           the streamline location, width and thalweg location
c
c OUTPUT VARIABLES:
c
c   facies           a realistic channel geometry fit to the streamline.
c
c
c
c EXTERNAL REFERENCES: azimuth    computes the azimuth
c                      spline     calculates spline
c                      splint     spline interpolation 
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none
  
c
c Conditioning Arrays
c

      allocate(assoc_tab(max_assoc,max_withinassoc),stat = test) 
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      assoc_tab=0.0
c
      allocate(assoc_cond(max_assoc,ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      assoc_cond=0.0
c

c
c Allocate arrays   
c 

      allocate(facies(nx,ny,nz),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      facies = kFF
c
      allocate(chamat(ndis0+1,ndata_chamat,0:ntime),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      chamat = 0.0
c
      allocate(spline_table(ndis0+1,ndata_spline_table,0:ntime),
     + stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      spline_table = 0.0
c
      allocate(spline_table2(ndis0+1,ndata_spline_table,0:ntime),
     + stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      spline_table2 = 0.0
c
      allocate(temp1dmat(ndis0),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      temp1dmat = 0.0
c
      if(write_surfaces) then
       allocate(ihs(nx,ny,0:ntime,2),stat = test)
       if(test.ne.0)then
        write(*,*)'ERROR 5: Allocation failed',
     +  ' due to insufficient memory.'
        stop
       end if
       ihs = -1
      end if
c
      allocate(chaprop(ndata_chaprop,ntime),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      chaprop = -9.9
c
      allocate(schedule(0:ntime),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      schedule = 1.0
c

c
c Horizontal Trend Draw arrays
c

      allocate(chadraw(ndis0,ndata_chadraw,nCHdraw),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      chadraw = -9.9
c
      allocate(chadrawprop(ndata_chadrawprop,nCHdraw),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      chadrawprop = -9.9
c
      allocate(chadrawwt(nCHdraw),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      chadrawwt = 0.0

c
c Spline arrays
c
      allocate(spline_s(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_d(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2d(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_a(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2a(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_x(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_w(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_y(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_t(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_c(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_i(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_z(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_h(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2x(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2y(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2w(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2t(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2c(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2i(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2z(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(spline_2h(ndis0+1),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c

c
c Test matrix
c 
 
      allocate(testmat(nx,ny,nz),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      testmat = 0.0
c

 
      continue                   

c
c Finished:
c
      return
      end




      subroutine avulsioninside(ianode,asinu,local_azi,tlength,
     + ctime1,ctime2)
    
c-----------------------------------------------------------------------
c
c     Generate a realistic channel streamline avulsion inside the model
c     ***************************************
c
c The subroutine generates a realistic channel streamline based on the 
c disturbed periodic model.
c
c INPUT VARIABLES:
c
c   ianode        the avulsion node (<1 then draw)   
c   asinu         the sinuosity    
c   local_azi     the principal channel direction (0 - y positive)
c
c GLOBAL VARIABLE:
C
c   sinu          approximate sinuosity of the initial channel
c
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none

c
c Local Variables
c
 
      integer   idis,idum,ierr
      real      delta_width,xp,ang,ang1,ang2
      real      b1,b2,k,h,m,s,val,phi
      real      mCHhalfwidth,CHwdratio,length
      real*8    x,y,cx,cy,p
      real,allocatable :: temp1Dmati(:)

c
c Looked Up Variables
c 
     
      integer   ndis
      real      CHdepth,CHdwratio,Chelev

c
c Passed Variables
c

      integer   ianode,ctime1,ctime2
      real      local_azi,asinu,tlength
 
c
c Function Calls
c

      real      resc
      real*8    acorni

c
c Test
c

      if(ctime1.eq.67) then
       continue
      end if

c
c Look up ndis from lookup table and pass on to the new channel
c
     
      write(*,*) ' Started: AvulsionInside streamline ',ctime1,
     + ' to ',ctime2
      ndis=chaprop(15,ctime1)
      CHdepth = chaprop(2,ctime1)
      CHdwratio = chaprop(3,ctime1)
      CHwdratio = 1.0/CHdwratio
      CHelev = chamat(1,11,ctime1) ! assume constant elevation

c
c Allocate arrays   
c

c
      allocate(temp1DMati(ndis0),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      temp1DMati=0.0

c
c Draw the avulsion location 
c

512   if(ianode.le.1) then ! if failed draw avulsion location again
       do idis=2,ndis-1
         temp1DMati(idis)=chamat(idis,6,ctime1)
       end do
       call probdraw(temp1DMati,ndis,ianode,val)
      end if

c
c Pass the angle at the avulsion node if not specified
c

      if(local_azi.lt.-999.0) then
       local_azi = chamat(ianode,9,ctime1)
      end if

c
c Disturbed periodic model - note observation showed a linear 
c relationship between s and sinuosity
c

      k = 0.3
      h = 0.8
      m = (450-local_azi)*(16.33/360.0)
      s = resc(1.0,2.0,1.0,13.0,asinu)

      phi = asind(h)
      b1 = 2.0*exp(-k*h)*cosd(k*cosd(phi))
      b2 = -1.0*exp(-2.0*k*h)

c
c Construct random N{0,1} vector and the weighting function
c
      
      do idis = 1,ndis0
       p = real(acorni(idum))
       call gauinv(p,xp,ierr)
       temp1dMati(idis) = xp*s+m
      end do
 
      ang1 = 450-local_azi
      ang2 = 450-local_azi

c
c Copy the nodes up to the avulsion node
c
      do idis = 1, ianode
       chamat(idis,1,ctime2) = chamat(idis,1,ctime1)
       chamat(idis,2,ctime2) = chamat(idis,2,ctime1)
       chamat(idis,5,ctime2) = chamat(idis,5,ctime1)
       chamat(idis,11,ctime2) = CHelev
      end do

      cx = chamat(ianode,1,ctime2)
      cy = chamat(ianode,2,ctime2)

c
c Avulsion: Stop criteria - length
c

      if(tlength.gt.0.0) then
       length = 0.0
       tlength = tlength
       do idis = ianode+1,ndis0
        ang= b1*ang1+b2*ang2+temp1dMati(idis)
        x = cx+step*cosd(ang)
        y = cy+step*sind(ang)        
        chamat(idis,1,ctime2) = x
        chamat(idis,2,ctime2) = y
        chamat(idis,11,ctime2) = CHelev
        length=length+sqrt((x-cx)**2.0+(y-cy)**2.0)
        if(length.gt.tlength) then
         ndis = idis
         goto 511
        end if
        ang2 = ang1
        ang1 = ang
        cx = x
        cy = y
       end do
      else

c
c Avulsion: Stop criteria - escape model
c

       do idis = ianode+1,ndis0
        ang= b1*ang1+b2*ang2+temp1dMati(idis)
        x = cx+step*cosd(ang)
        y = cy+step*sind(ang)
        if(x.lt.xmin.or.x.gt.xmax.or.y.lt.ymin.or.y.gt.ymax) then
         if(idis.lt.ndis0/10) then
           ianode = -1
           goto 512
         end if 
         ndis = idis-1
         goto 511
        end if        
        chamat(idis,1,ctime2) = x
        chamat(idis,2,ctime2) = y
        chamat(idis,11,ctime2) = CHelev
        ang2 = ang1
        ang1 = ang
        cx = x
        cy = y
       end do
      end if


511   continue
      chaprop(15,ctime2)=ndis

c
c Calculate a 1D G{} RF for the channel width
c

      call onedrf(nCHcor,ndis0,temp1dMati,CHdepth,
     +  stdevCHdepth2) 

c
c Convert to half width
c   
      
      do idis = 1, ndis
       temp1dMati(idis) = temp1dMati(idis)*CHwdratio*0.5 ! CHwdratio constant
      end do

      delta_width = chamat(ianode,5,ctime2)-temp1dMati(ianode)
      do idis = ianode+1, ndis
       chamat(idis,5,ctime2) = temp1dMati(idis)+delta_width
      end do

      deallocate(temp1dMati)

c
c Finished:
c
      
      write(*,*) ' Finished: AvulsionInside streamline ',ctime1,
     + ' to ',ctime2
      return

      end


      subroutine azimuth(x1,x2,y1,y2,azi)
    
c-----------------------------------------------------------------------
c
c     Calculate the azimuth
c     *********************
c
c Given the end points of a line extending from (x1,y1) to (x2,y2) this
c subroutine returns the azimuth with 0.0 being in the y positive 
c direction and increasing angle clockwise. 
c
c
c INPUT VARIABLES:
c
c   x1            difference in x   
c   y1            difference in y 
c   x1            difference in x   
c   y1            difference in y
c
c Author: Michael Pyrcz                                  DATE: 2001-2002 
c-----------------------------------------------------------------------
      
      implicit none

c
c Local Variables
c

      real      PI

c
c Passed Variables
c

      real      x1,x2,y1,y2,di,dj,azi
      parameter(PI=3.14159)
      
c
c Calculate the delta x and delta y
c   

      di = x2 - x1
      dj = y2 - y1

c
c Check for 0 or 180
c
      
      if(di.eq.0.0) then
       if(dj.gt.0.0) azi = 0.0
       if(dj.lt.0.0) azi = 180.0
       goto 566
      end if

c
c First quandrant
c

      if(di.gt.0.0.and.dj.ge.0.0) then
       azi = 90.0-180*atan(dj/di)/pi
       goto 566
      end if 

c
c Second quandrant
c
      if(di.lt.0.0.and.dj.ge.0.0) then
       azi = 270.0-180.0*atan(dj/di)/pi
       goto 566
      end if

c
c Third quandrant 
c
      if(di.lt.0.0.and.dj.lt.0.0) then
       azi = 270.0-180.0*atan(dj/di)/pi
       goto 566
      end if

c
c Fourth quandrant
c

      if(di.gt.0.0.and.dj.lt.0.0) then
       azi = 90.0-180.0*atan(dj/di)/pi
       goto 566
      end if

566   continue

c
c Finished:
c
      return
      end



      subroutine buildCHtable
    
c-----------------------------------------------------------------------
c
c     Generate a realistic channel streamline 
c     ***************************************
c
c The subroutine generates a realistic channel streamline based on the 
c disturbed periodic model.
c
c MODIFIED TO WRITE TO TEMP CALIBRATION TABLE
c
c INPUT VARIABLES:
c
c
c GLOBAL VARIABLE:
C
c   sinu          approximate sinuosity of the initial channel
c
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none

c
c Local Variables
c

      integer   i,ix,iy,iCHdraw,cdraw,idum,ierr,idis,inflag
      real      x,y,x0,y0,cx,cy,azi,ang,b1,b2,k,h,m,s,phi,xp
      real      ang1,ang2
      real      CHazi,CHsinu,CHdepth,CHdwratio,CHwdratio,mCHhalfwidth
      real*8    p
      real,allocatable :: temp1Dmati(:)

c
c Looked Up Variables
c

      integer   ndis

c
c Called Functions
c

      real      resc
      real*8    acorni
      real      cosd

c
c Allocate arrays   
c

c
      allocate(temp1DMati(ndis0),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c

c
c Generate a 1D Gaussian Random function
c
      
      k = 0.3
      h = 0.8
      do iCHdraw = 1, nCHdraw
       p = real(acorni(idum))
       call gauinv(p,xp,ierr)
       CHazi = xp*stdevCHazi+mCHazi
       chadrawprop(1,iCHdraw) = CHazi
       p = real(acorni(idum))
       call gauinv(p,xp,ierr)
       CHsinu = min(1.9,max(1.1,xp*stdevCHsinu+mCHsinu))
       chadrawprop(2,iCHdraw) = CHsinu
       x0 = 0.0
       p = real(acorni(idum))
       if(stdevCHsource.gt.0.0) then
        call gauinv(p,xp,ierr)
        y0 = min(ymax,max(ymin,xp*stdevCHsource+mCHsource))
       else
        y0 = ymin+p*(ymax-ymin) ! if source stdev<0 apply uniform
       end if

       m = (450-CHazi)*(16.33/360.0)
       s = resc(1.0,2.0,1.0,13.0,CHsinu)

       phi = asind(h)
       b1 = 2.0*exp(-k*h)*cosd(k*cosd(phi))
       b2 = -1.0*exp(-2.0*k*h)

c
c Construct random N{0,1} vector and the weighting function
c
      
512    do i = 1,ndis0
        p = real(acorni(idum))
        call gauinv(p,xp,ierr)
        temp1dMati(i) = xp*s+m
       end do
 
        ndis = ndis0
       ang1 = 450-CHazi
       ang2 = 450-CHazi
       cx = x0
       cy = y0
       chadraw(1,1,iCHdraw) = cx
       chadraw(1,2,iCHdraw) = cy
       do i = 2,ndis0
        ang= b1*ang1+b2*ang2+temp1dMati(i)
        x = cx+step*cosd(ang)
        y = cy+step*sind(ang)
        if(x.gt.xmax.or.x.lt.xmin.or.y.gt.ymax.or.y.lt.ymin) then
         if(i.lt.ndis0/10) goto 512 ! regenerate if short streamline
         ndis = i-1
         goto 511
        end if        
        chadraw(i,1,iCHdraw) = x
        chadraw(i,2,iCHdraw) = y
        ang2 = ang1
        ang1 = ang
        cx = x
        cy = y
       end do
511    continue
       chadrawprop(3,iCHdraw) = ndis

c
c Calculate a 1D G{} RF for the channel width
c

       p = real(acorni(idum))
       call gauinv(p,xp,ierr)
       CHdepth = max(0.0,xp*stdevCHdepth+mCHdepth)
       chadrawprop(6,iCHdraw) = CHdepth
       temp1dMati = 0.0
       call onedrf(nCHcor,ndis,temp1dMati,CHdepth,
     +  stdevCHdepth2)

       p = real(acorni(idum))
       call gauinv(p,xp,ierr) 
       CHwdratio = max(0.0,xp*stdevCHwdratio+mCHwdratio)
       CHdwratio = 1.0/CHwdratio
       chadrawprop(7,iCHdraw) = CHdwratio ! store CH depth:width ratio

c
c Convert to half width
c   
       do idis = 1, ndis
        temp1dMati(idis) = temp1dMati(idis)*CHwdratio*0.5 ! CHwdratio constant
       end do

       do idis = 1, ndis
        chadraw(idis,3,iCHdraw) = max(0.01,temp1dMati(idis))
       end do
      end do 
 
      deallocate(temp1dMati)

c
c Calculate weight of each path
c

      do iCHdraw = 1, nCHdraw
       chadrawprop(4,iCHdraw) = 0.0
      end do

      do iCHdraw = 1, nCHdraw
       ndis = chadrawprop(3,iCHdraw)
       do idis = 1,ndis 
        call getindx(nx,xmn,xsiz,chadraw(idis,1,iCHdraw),ix,inflag)
        call getindx(ny,ymn,ysiz,chadraw(idis,2,iCHdraw),iy,inflag)
        chadrawprop(4,iCHdraw)=chadrawprop(4,iCHdraw)+horimat(ix,iy)
       end do
       if(ndis.gt.0) then
         chadrawprop(4,iCHdraw)=chadrawprop(4,iCHdraw)/ndis ! calculate the avg.
       else
         chadrawprop(4,iCHdraw) = 0
       end if
      end do 

c
c Copy the weights to a vector so it can be handled by probdraw
c

      do iCHdraw = 1, nCHdraw
       chadrawwt(iCHdraw) = chadrawprop(4,iCHdraw)
      end do

c
c Finished:
c
      return
      end



      subroutine calc_levee(nstart,nend,ctime,k)
c------------------------------------------------------------------------
c
c    Generate LV element
c    *******************
c
c Given model dimensions this member function constructs a channel levee 
c element assosiate with the channel spline.
c
c INPUT VARIABLES
c
c node_start        the proximal extend of the levee (idis) - most likely 
c                   set to idis=0
c node_end          the distal extent of the levee (idis)
c ctime             the CHcurrent
c k                 the levee category code (real)
c
c GLOBAL INPUT VARIABLES
c
c levee_width       the levee width (currrently constant modified by 
c                   curvature and distal thinning 
c levee_height      maximum levee hieght above channel datum (top of 
c                   channel) 
c levee_depth       maximum depth of incision of the levee below the 
c                   channel datum (top of channel)
c levee_asymmetry   fraction 0.0 to 1.0 for the strength in levee 
c                   asymmetry on bends 0.0 = symmetrical, and 1.0 
c                   doubled on cutbank and none on point bar. 
c levee_thin        distal thinning of levees 0.0 to 1.0.  1.0 = thin to 
c                   0 at node_end, 0.0 = no thinning  
c
c GLOBAL OUTPUT
c
c facies(nx,ny,nz)  place LV element within the facies model
c
c Author: Michael Pyrcz                                  DATE: 2003-2004  
c-----------------------------------------------------------------------

      use     nr, only : splint
      use     streamsimMod
      implicit none

c
c Local Variables
c
      integer ix,iy,iz,idis,jdis,idiscr,cnode
      real    max_search,val
      real    x_loc,y_loc,z_loc,distance,close_distance,
     +         close_distance_datum,channel_elevation
      real    half_width,start_dist,end_dist,delta_dist,s_loc,s_close
      real    curvature,x_inter,y_inter,x_value_0,y_value_0
      real    x_value,y_value,azi1,azi2,dazi
      real    CHelev,CHhalfwidth,levee_width_scale

c
c Looked Up Variables
c

      integer ndis,spline_count
      real    levee_depth,levee_width,levee_height,levee_asym,
     +         levee_thin
      real    levee_top,levee_bottom,factor,factor_asym,factor_thin
      real    maxCHcurvature,maxCHhalfwidth

c
c Passed Variables
c

      integer nstart,nend,ctime
      real    k

c
c Look up associated parameters
c

      write(*,*) ' Started: LV streamline ',ctime
      ndis = chaprop(15,ctime)
      spline_count = chaprop(16,ctime)
      maxCHcurvature = chaprop(18,ctime)
      maxCHhalfwidth = chaprop(17,ctime)
 
      levee_depth = chaprop(5,ctime)
      levee_width = chaprop(6,ctime) 
      levee_width_scale = levee_width/6.0   ! this is how the LV equation works            
      levee_height = chaprop(7,ctime)
      levee_asym = chaprop(8,ctime)
      levee_thin = chaprop(9,ctime)
      max_search = maxCHhalfwidth+levee_width*1.5 ! add extra for interpolation

c
c Check the start and end nodes
c
      nstart=nstart+1
      nstart = max(1,nstart)
      nend = min(ndis,nend)

c
c Copy the spline vectors to the spline table
c

      do idis = 1, spline_count
       spline_x(idis)=spline_table(idis,1,ctime)
       spline_y(idis)=spline_table(idis,2,ctime)
       spline_s(idis)=spline_table(idis,3,ctime)
       spline_t(idis)=spline_table(idis,4,ctime)
       spline_w(idis)=spline_table(idis,5,ctime)
       spline_c(idis)=spline_table(idis,6,ctime)
       spline_d(idis)=spline_table(idis,8,ctime)
       spline_i(idis)=spline_table(idis,9,ctime)
       spline_z(idis)=spline_table(idis,11,ctime)            
      end do
      do idis = 1, spline_count
       spline_2x(idis)=spline_table2(idis,1,ctime)
       spline_2y(idis)=spline_table2(idis,2,ctime)
       spline_2t(idis)=spline_table2(idis,4,ctime)
       spline_2w(idis)=spline_table2(idis,5,ctime)
       spline_2c(idis)=spline_table2(idis,6,ctime)
       spline_2d(idis)=spline_table2(idis,8,ctime)
       spline_2i(idis)=spline_table2(idis,9,ctime)
       spline_2z(idis)=spline_table2(idis,11,ctime)            
      end do

c
c Loop Over All Model Areal Locations 
c
      
      do iy=1, ny
       y_loc=real(iy)*ysiz+ymn  
       do ix = 1, nx
        x_loc=real(ix)*xsiz+xmn   
        close_distance = big_number
   
c
c Find the Nearest Control Node
c

        do idis = nstart, nend
        distance=(chamat(idis,1,ctime)-x_loc)**2.0+
     +    (chamat(idis,2,ctime)-y_loc)**2
         if(distance<close_distance) then
          cnode = idis
          close_distance = distance
         end if
        end do
        close_distance = sqrt(max(0.001,close_distance)) ! calculate the dist 
      
c
c Refine the Search
c

        if(close_distance.lt.max_search) then      
         jdis=max(nstart,cnode-1)
         start_dist = chamat(jdis,3,ctime)
         jdis=min(nend-1,cnode+1)
         end_dist = chamat(jdis,3,ctime)
         delta_dist = end_dist-start_dist
         distance = big_number
         close_distance = big_number

         do idiscr=0,ndiscr
          s_loc=start_dist+(real(idiscr)/real(ndiscr))*delta_dist
          x_inter = splint(spline_s,spline_x,spline_2x,s_loc,spline_count)
          y_inter = splint(spline_s,spline_y,spline_2y,s_loc,spline_count)
          distance = (x_inter-x_loc)**2.0+(y_inter-y_loc)**2.0
          if(distance.lt.close_distance) then
           s_close=s_loc
           close_distance=distance
          end if
         end do
         close_distance=sqrt(max(close_distance,0.001)) ! prevent a math error 
         x_value = splint(spline_s,spline_x,spline_2x,s_close,spline_count)
         y_value = splint(spline_s,spline_y,spline_2y,s_close,spline_count)
         CHelev = splint(spline_s,spline_z,spline_2z,s_close,spline_count)
         CHhalfwidth = splint(spline_s,spline_w,spline_2w,s_close,spline_count)
         curvature = splint(spline_s,spline_c,spline_2c,s_close,spline_count)

c
c Levee thinning - factor_thin is centered on 1.0 with range 1+levee_thin to 1-levee_thin
c

         factor_thin=(1.0+levee_thin)-2.0*(s_close-spline_s(1))/
     +      (spline_s(spline_count)-spline_s(1))*levee_thin ! prop along spline

c
c Improved Left / Right Determination
c
         x_value_0 = splint(spline_s,spline_x,spline_2x,s_close-0.001,spline_count)
         y_value_0 = splint(spline_s,spline_y,spline_2y,s_close-0.001,spline_count)
         call azimuth(x_value_0,x_value,y_value_0,y_value,azi1)        
         call azimuth(x_value,x_loc,y_value,y_loc,azi2)
         dazi=azi2-azi1
         if(abs(dazi)>180.0) then
          dazi=(azi1-180.0)-azi2
         end if

c
c Determine if Cutbank or Pointbar
c
     
         if((dazi.le.0.0.and.curvature.le.0.0).or.
     +     (dazi.gt.0.0.and.curvature.gt.0.0)) then ! cutbank side
          factor_asym=1.0+levee_asym*
     +     abs(curvature/maxCHcurvature)
         else
          factor_asym=1.0-levee_asym*
     +       abs(curvature/maxCHcurvature)
         end if
 
c
c Levee Cross Section
c
         levee_top=-2.0 ! prevent the rare escapes 
         levee_bottom=-1.0 ! prevents the rare escapes 
         factor=factor_asym*factor_thin ! combine asymmetry and thinning  
         close_distance_datum = close_distance - CHhalfwidth ! set origin at CH 
         if(close_distance.lt.0.0) then ! check for case inside the CH  
          levee_top=CHelev
          levee_bottom=CHelev-levee_depth       
         else
          levee_top=levee_height*(close_distance_datum/
     +     (levee_width_scale*factor))*exp(-1.0*(close_distance_datum)/
     +     (levee_width_scale*factor))+CHelev
          levee_bottom=CHelev-levee_depth*((levee_width*factor-
     +     close_distance)/(levee_width))
          end if
         if(levee_top.gt.levee_bottom) then
c          write(*,*) levee_top,levee_bottom
          do iz = 1, nz
           z_loc=zmn+real(iz)*zsiz
           if(z_loc.ge.levee_bottom.and.z_loc.le.levee_top) then 
            if(facies(ix,iy,iz).lt.1.0) then
             NTGcurrent=NTGcurrent+1 ! increment  NTG counter
             end if  
            facies(ix,iy,iz)=k 
           end if ! if withn levee geometry
          end do ! loop over all iz
         end if ! levee has positive thickness
        end if ! if location within 2 times the levee length
       end do ! loop over all ix
      end do ! loop over all iy

c
c Finished:
c
      write(*,*) ' Finished: LV streamline ',ctime
      return
      end




      subroutine calc_lobe(node_start,lobe_LL,lobe_WW,lobe_l,lobe_w,
     + lobe_hwratio,lobe_dwratio,ctime,k)

c------------------------------------------------------------------------
c
c    Generate CS,FS,SH or DA element
c    ****************************
c
c Given model dimensions this member function constructs CS,FS,SH or DA 
c elements assosiate with the channel spline.
c
c INPUT VARIABLES
c
c node_start        the proximal extend of the levee (idis) - most likely 
c                   set to idis=0
c lobe_LL           the length of the lobe
c lobe_WW           the maximum width of the lobe
c lobe_l            the length of the lobed at maximum width        
c lobe_w            the width of the lobe at the proximal edge
c lobe_hwratio      the height to width ratio of the lobe
c lobe_dwratio      the depth to width ratio of the lobe
c ctime             the CHcurrent
c kLV               the levee category code (real)
c
c GLOBAL INPUT VARIABLES
c
c levee_width       the levee width (currrently constant modified by 
c                   curvature and distal thinning 
c levee_height      maximum levee hieght above channel datum (top of 
c                   channel) 
c levee_depth       maximum depth of incision of the levee below the 
c                   channel datum (top of channel)
c levee_asymmetry   fraction 0.0 to 1.0 for the strength in levee 
c                   asymmetry on bends 0.0 = symmetrical, and 1.0 
c                   doubled on cutbank and none on point bar. 
c levee_thin        distal thinning of levees 0.0 to 1.0.  1.0 = thin 
c                   to 0 at node_end, 0.0 = no thinning  
c
c GLOBAL OUTPUT
c
c facies(nx,ny,nz)  place CS,FS,SH or DA element within the facies model
c
c Author: Michael Pyrcz                                  DATE: 2003-2004  
c-----------------------------------------------------------------------

      use nr, only : splint
      use streamsimMod
      implicit none

c
c Local Variables
c

      integer ix,iy,iz
      integer idis,idiscr,close_node,jdis,FS_local_count,net_local_count
      real    x_loc,y_loc,z_loc,s_loc,distance,close_distance,lobe_s,
     +           x_inter,y_inter,s_close,x_value,y_value
      real    lobe_top,lobe_bottom,start_dist,end_dist,delta_dist,
     +           lobe_datum,half_width,y_function
      real    CHhalfwidth

c
c Passed Variables
c

      integer node_start,ctime
      real    lobe_LL,lobe_WW,lobe_l,lobe_w,lobe_hwratio,lobe_dwratio,k

c
c Looked Up Variables
c

      integer ndis,spline_count
      real    maxCHcurvature,maxCHhalfwidth 

c
c Check the start and end nodes
c

      write(*,*) ' Started: Lobe streamline ',ctime
      node_start = Int(max(0,node_start))
      ndis = chaprop(15,ctime)
      spline_count = chaprop(16,ctime)
      maxCHcurvature = chaprop(18,ctime)
      maxCHhalfwidth = chaprop(17,ctime)

      if(k.eq.1.1) then
       continue
      end if

c
c Copy the spline vectors to the spline table
c

      do idis = 1, spline_count
       spline_x(idis)=spline_table(idis,1,ctime)
       spline_y(idis)=spline_table(idis,2,ctime)
       spline_s(idis)=spline_table(idis,3,ctime)
       spline_t(idis)=spline_table(idis,4,ctime)
       spline_w(idis)=spline_table(idis,5,ctime)
       spline_c(idis)=spline_table(idis,6,ctime)
       spline_d(idis)=spline_table(idis,8,ctime)
       spline_i(idis)=spline_table(idis,9,ctime)
       spline_z(idis)=spline_table(idis,11,ctime)            
      end do
      do idis = 1, spline_count
       spline_2x(idis)=spline_table2(idis,1,ctime)
       spline_2y(idis)=spline_table2(idis,2,ctime)
       spline_2t(idis)=spline_table2(idis,4,ctime)
       spline_2w(idis)=spline_table2(idis,5,ctime)
       spline_2c(idis)=spline_table2(idis,6,ctime)
       spline_2d(idis)=spline_table2(idis,8,ctime)
       spline_2i(idis)=spline_table2(idis,9,ctime)
       spline_2z(idis)=spline_table2(idis,11,ctime)            
      end do

c
c  Lobe Parameters
c

      lobe_s=chamat(node_start,3,ctime)
      lobe_l=lobe_l+lobe_s ! translate to start location of lobe / distance  
      lobe_LL=Lobe_LL+lobe_s
      if(lobe_w.lt.0.0) then
       lobe_w = chamat(node_start,5,ctime) ! set the start width to CHhalfwidth 
      else
       lobe_w=lobe_w/2.0 ! divide the w width by 2 - half width are assumed 
      end if
      lobe_WW=lobe_WW/2.0 ! divide the W Width by 2 - half widths are assumed 


c
c Loop Over All Model Areal Locations 
c

      do iy = 1, ny
       y_loc=real(iy)*ysiz+ymn  
       do ix = 1, nx  
        x_loc=real(ix)*xsiz+xmn
        close_distance = big_number
      
c
c Find the Nearest Control Node
c

        do idis = node_start, ndis
         distance=(chamat(idis,1,ctime)-x_loc)**2.0+
     +   (chamat(idis,2,ctime)-y_loc)**2.0
         if(distance.lt.close_distance) then
          close_node = idis
          close_distance = distance
         end if
        end do
        close_distance = sqrt(max(0.001,close_distance))
      
c
c Refine the Search
c

        if(close_distance.le.lobe_WW) then      
         jdis=max(node_start,close_node-1)
         start_dist = chamat(jdis,3,ctime)
         jdis=min(ndis,close_node+1)
         end_dist = chamat(jdis,3,ctime)
         delta_dist = end_dist-start_dist
        
         distance = big_number
         close_distance = big_number

         do idiscr=0,ndiscr
          s_loc=start_dist+(float(idiscr)/float(ndiscr))*delta_dist
          x_inter = splint(spline_s,spline_x,spline_2x,s_loc,spline_count)
          y_inter = splint(spline_s,spline_y,spline_2y,s_loc,spline_count)          
          distance = (x_inter-x_loc)**2.0+(y_inter-y_loc)**2.0
          if(distance.lt.close_distance) then
           s_close=s_loc
           close_distance=distance
          end if
         end do
         close_distance=sqrt(max(close_distance,0.001)) ! prevent a math error 
         x_value = splint(spline_s,spline_x,spline_2x,s_close,spline_count)
         y_value = splint(spline_s,spline_y,spline_2y,s_close,spline_count)
         lobe_datum = splint(spline_s,spline_z,spline_2z,s_close,spline_count)      
         CHhalfwidth = splint(spline_s,spline_w,spline_2w,s_close,spline_count)
 
c
c Lobe Cross Section
c
         y_function=0.0
         lobe_top=-2.0
         lobe_bottom=-1.0
         if(s_close.gt.lobe_s.and.s_close.le.lobe_l) then 
          y_function=lobe_WW-(lobe_WW-lobe_w)*(1.0-(s_close-lobe_s)
     +       /(lobe_l-lobe_s))**2.0 ! new parabolic model
         end if
         if(s_close.gt.lobe_l.and.s_close.le.lobe_LL) then 
          y_function=lobe_WW*sqrt(1.0-((s_close-lobe_l)/
     +       (lobe_LL-lobe_l))**2.0)
         end if
         if(s_close.le.lobe_LL) then
          lobe_top=lobe_datum+y_function*lobe_hwratio*
     +       (1.0-((close_distance**2.0)/(y_function**2.0)))
          lobe_bottom=lobe_datum-y_function*lobe_dwratio*
     +     (1.0-((close_distance**2.0)/(y_function**2.0)))       
         end if
         if(lobe_top.gt.lobe_bottom) then
          do iz = 1, nz
           z_loc=zmn+real(iz)*zsiz
           if(z_loc.ge.lobe_bottom.and.z_loc.le.lobe_top) then 
            if(facies(ix,iy,iz).lt.1.0) then
             NTGcurrent=NTGcurrent+1 ! increment counter if FF facies replaced
            end if
            facies(ix,iy,iz)=k 
            testmat(ix,iy,iz)=k
           end if ! if within the lobe element
          end do ! loop over all iz
         end if ! if lobe has positive thickness
        end if ! if within W_lobe distance
       end do ! loop over all ix
      end do ! loop over all iy

c
c Finished:
c

      write(*,*) ' Finished: Lobe streamline ',ctime
      return
      end



      subroutine calcusb(ctime)
c-----------------------------------------------------------------------
c
c      Calculated the near bank velocity
c      *********************************
c
c Apply one time step to a river channel.  The equations and concepts 
c are based on A.D. Howard, 1992, and Sun et. al., 1996.  
c
c
c GLOBAL VARIABLES:
c
c   ndis             the number of streamline discretizations
c   chamat           the 2D stream line with x in col 1, and y in col 2
c   
c
c OUTPUT VARIABLES:
c
c   ->chamat(ndis,10,ctime)   the near bank velocity
c
c
c EXTERNAL REFERENCES: 
c       
c Author: Michael Pyrcz                                  DATE: 2003-2004                
c-----------------------------------------------------------------------
 
      use       streamsimMod
      implicit none

c
c Local Variables
c
 
      integer   j,idis,start
      real      Csi,ds,dCsds,usip,part1,part2,part3,part4
      real      inte,mCHhalfwidth

c
c Looked Up Variables
c 
     
      integer   ndis
      real      CHdwratio,factor

c
c Passed Variables
c

      integer   ctime

c
c Retreive the number of control nodes
c

      write(*,*) ' Started: USB streamline ',ctime
      ndis=chaprop(15,ctime)
      CHdwratio=chaprop(3,ctime)
      mCHhalfwidth=mCHdepth/CHdwratio

c
c Solve the recusive equation (Sun et al., 1996, EQN 15, p. 2940) 
c

      chamat(1,10,ctime) = 0.0 
      do idis = 2, ndis
       Csi = chamat(idis,6,ctime)
       ds = chamat(idis,7,ctime)
       dCsds = chamat(idis,8,ctime)
       usip = chamat(idis-1,10,ctime)
       part1 = -1.0*us0*Csi
       part2 = mCHhalfwidth*Cf/us0 ! should this be doubled (full width)?
       part3 = (us0**4.0)/(g*h0**2.0)
       part4 = (scour_factor+2.0)*((us0**2.0)/h0)
       ds = 0.0
       inte = 0.0
       start = max(1,idis-30)
       do j = idis, start, -1
        ds = ds + chamat(j,7,ctime) 
        inte = inte+exp(-2.0*Cf*ds/h0)*chamat(j,6,ctime)
       end do  
       chamat(idis,10,ctime) = part1+part2*(part3+part4)*inte
       if(idis.lt.ndis/10) then
         factor = real(idis-2)/real(ndis/10-(idis-2))
         chamat(idis,10,ctime) = chamat(idis,10,ctime) * factor
       end if 
      end do
      continue

c
c Finished:
c

      write(*,*) ' Finished: USB streamline ',ctime
      return
      end




      subroutine checkwell(rcheck)
    
c-----------------------------------------------------------------------
c
c     Check the Well Reproduction 
c     ***************************
c
c Given a well data set and a final model - check to see if the net 
c intervals are honored without unwarranted intercepts.
c
c
c INPUT VARIABLES:
c
c
c OUTPUT VARIABLES:
c rcheck = false if there is unwarranted intercepts
c percentViolate = percent violation for each well
c totalViolate = total percent violation for all wells
c
c Author: Michael Pyrcz                                  DATE: 2003-2004
c-----------------------------------------------------------------------
      
      use streamsimMOD
      implicit none

      logical  rcheck
      integer  iwell,iCHdata,rviolate,iNonNetBlk
      integer  wellnum,wellix,welliy,iz,inflag
      real     z,topz,botz,sum

      rcheck = .true.
      rviolate = 0
      sum = 0


      allocate(percentViolate(nwell),stat = test)
      if(test.ne.0)then
       write(*,*)'Percent Violate Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if

c
c Loop over vertical wells
c
   
      do iwell = 1, nwell
       wellnum = wellheader(iwell,1)
       wellix = wellheader(iwell,2)
       welliy = wellheader(iwell,3)

c
c Calculate well log
c

       well_log = 0
       iNonNetBlk = nz
       do iCHdata = 1, nCHdata
        if(int(welldata(iCHdata,1)).eq.wellnum) then
         topz = welldata(iCHdata,4)
         botz = welldata(iCHdata,5)
         do iz = 1, nz
          z = real(iz-1)*zsiz+zmn
          if(z.gt.botz.and.z.lt.topz) then 
           well_log(iz) = 1
           iNonNetBlk = iNonNetBlk-1 
          end if
         end do
        end if
       end do
   
c
c Check well log with model
c

       do iz = 1, nz
        write(*,*) 'wellix welliy iz ',wellix,welliy,iz
        write(*,*) 'facies = ',facies(wellix,welliy,iz)
        if(facies(wellix,welliy,iz).eq.kCH.or.
     +   facies(wellix,welliy,iz).eq.kLA.or.
     +   facies(wellix,welliy,iz).eq.kLV.or.
     +   facies(wellix,welliy,iz).eq.kCS.or.
     +   facies(wellix,welliy,iz).eq.kFFCH) then 
         if(well_log(iz).ne.1) then
          rcheck = .false.
          rviolate = rviolate + 1
         end if
        else if(facies(wellix,welliy,iz).ne.kCH.and.
     +   facies(wellix,welliy,iz).ne.kLA.and.
     +   facies(wellix,welliy,iz).ne.kLV.and.
     +   facies(wellix,welliy,iz).ne.kCS.and.
     +   facies(wellix,welliy,iz).ne.kFFCH) then
         if(well_log(iz).ne.0) then
          rcheck = .false.
          rviolate = rviolate + 1
         end if
        end if         
       end do ! loop over all z blocks
c
c Calculate percent of blocks that violates the well intercept
c
      percentViolate(iwell) = (real(rviolate)/iNonNetBlk)*100
      sum = sum + percentViolate(iwell)
      end do ! loop over vertical wells
      totalViolate = 0.
      if(nwell > 0) then
        totalViolate = sum/nwell
      end if
c
c Finished:
c
      return
      end


      subroutine constructmodel(itime)
    
c-----------------------------------------------------------------------
c
c     To construct current streamline and architectual elements
c     *********************************************************
c
c The subroutine is called to construct current streamline with attached 
c  architectural elements.
c
c
c INPUT VARIABLES:
c      itime = current channel constructed
c
c Author: Fueng Zabel                                    DATE: 2004-2005
c-----------------------------------------------------------------------
      use       streamsimMod
      implicit none
c
c General Variables
c
      integer   idum,ierr,ilevel,idis,ndis,nstart,regrid,extrapolate
      real      color_current,val
      real      CHdepth
      real      CSnum,CSnumlobe,CSsinu
      integer   iCSnum,CSnode,iCSnumlobe
      real      curvature,CSazi,CSLOLL,CSLOWW,CSLOw,CSLOl
      real      CSLO_dwratio,CSLO_hwratio,xp
      real*8    p,acorni
c
c Conditioning 
c
      integer   itime
c
c Streamline Parameters
c
      real   CHdwratio
      real   LVwidth,LVheight,LVdepth,LVasym,LVthin

      write(*,*) 'Constructing Streamline #', itime
      CHdepth = chaprop(2,itime)
      CHdwratio = chaprop(3,itime)
      LVdepth = chaprop(5,itime)
      LVwidth = chaprop(6,itime)
      LVheight = chaprop(7,itime)
      LVasym = chaprop(8,itime)
      LVthin = chaprop(9,itime)
      CSnum = chaprop(11,itime)
      CSnumlobe = chaprop(14,itime)
      CSsource = chaprop(12,itime)

c
c      Generate the architectural elements
c
c Generate CS again / should have minimal impact on NTG 
c 
      write(*,*) '  Element: CS'
      CSsinu = 1.2 ! currently hard coded as a low sinuosity
      if(CSnum.gt.0.0) then
        do iCSnum = 1, int(CSnum)
         ndis = chaprop(15,itime)
         do idis = 1, ndis
          temp1Dmat(idis) = abs(chamat(idis,6,itime))
         end do
         CSnode = 1  ! initialize variable
         call probdraw(temp1Dmat,ndis,CSnode,val) 
         write(*,*) 'Element: CS test'      
         CSazi = chamat(CSnode,9,itime)
         curvature = chamat(CSnode,6,itime)
         if(curvature.gt.0.0) then
          CSazi = CSazi-90.0
         else
          CSazi = CSazi+90.0
         end if
         CHdepth = chamat(CSnode,5,itime)*2.0*CHdwratio
         do iCSnumlobe = 1, int(CSnumlobe)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOLL = max(0.001,xp*stdevCSLOLL+mCSLOLL)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOWW = max(0.001,xp*stdevCSLOWW+mCSLOWW)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOw = max(0.001,xp*stdevCSLOw+mCSLOw)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOl = max(0.001,xp*stdevCSLOl+mCSLOl)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLO_hwratio = max(0.001,xp*stdevCSLO_hwratio+mCSLO_hwratio)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLO_dwratio = max(0.001,xp*stdevCSLO_dwratio+mCSLO_dwratio)
          call avulsioninside(CSnode,CSsinu,CSazi,CSLOLL*1.5,itime,
     +     CHCStemp) ! set the branch length for stop criteria
          regrid = 0
          extrapolate = 0
          call curvature2(CHCStemp,extrapolate,regrid) ! CS branches use temp 
          if(kCS_local.eq.1.1) then
           continue
          end if
          call calc_lobe(CSnode,CSLOLL,CSLOWW,CSLOl,CSLOw,
     +     CSLO_hwratio,CSLO_dwratio,CHCStemp,kCS_local)
         end do
        end do
        continue
      end if

c
c Generate LV
c     

      if(LVwidth.gt.0.0.and.LVdepth+LVheight.gt.0.0) then       
        write(*,*) '  Element: LV'
        nstart=1
        ndis = chaprop(15,itime)
        call calc_levee(nstart,ndis,itime,kLV_local)
      end if       

c
c
c Generate CH
c
c
      write(*,*) '  Element: CH'
      call genchannel(itime,kCH_local)
      chaprop(22,itime) = kCH_local

c
c Return to the main program:
c

      continue
      return
      end


      subroutine constructmodelAll
    
c-----------------------------------------------------------------------
c
c     To construct streamline and facies models
c     ******************************************
c
c The subroutine is called to contruct overall streamline and facies 
c models with color scale adjustment.
c
c Author: Fueng Zabel                                    DATE: 2004-2005
c-----------------------------------------------------------------------
      use       streamsimMod
      implicit none
c
c General Variables
c
      integer   idum,ierr,ilevel,idis,ndis,nstart,regrid,extrapolate
      real      color_current,val
      real      CHdepth
      real      CSnum,CSnumlobe,CSsinu
      integer   iCSnum,CSnode,iCSnumlobe
      real      curvature,CSazi,CSLOLL,CSLOWW,CSLOw,CSLOl
      real      CSLO_dwratio,CSLO_hwratio,xp
      real*8    p,acorni
      real code
c
c Conditioning 
c
      integer   itime
c
c Streamline Parameters
c
      real   CHdwratio
      real   LVwidth,LVheight,LVdepth,LVasym,LVthin
      real   FFCHprop

c
c Reset the architectural element model
c

      facies = kFF 

       kFFCH_local = kFFCH
      kFF_local = kFF
      kCS_local = kCS
      kLV_local = kLV
      kLA_local = kLA
      kCH_local = kCH
      color_current = 0.0

      do itime = 1, CHcurrent-1
       write(*,*) 'Constructing Streamline #', itime
       CHdepth = chaprop(2,itime)
       CHdwratio = chaprop(3,itime)
       LVdepth = chaprop(5,itime)
       LVwidth = chaprop(6,itime)
       LVheight = chaprop(7,itime)
       LVasym = chaprop(8,itime)
       LVthin = chaprop(9,itime)
       CSnum = chaprop(11,itime)
       CSnumlobe = chaprop(14,itime)
       CSsource = chaprop(12,itime)

c
c      Generate the architectural elements
c
c Generate CS again / should have minimal impact on NTG 
c 
       write(*,*) '  Element: CS'
       CSsinu = 1.2 ! currently hard coded as a low sinuosity
       if(CSnum.gt.0.0) then
        do iCSnum = 1, int(CSnum)
         ndis = chaprop(15,itime)
         do idis = 1, ndis
          temp1Dmat(idis) = abs(chamat(idis,6,itime))
         end do
         CSnode = 1  ! initialize variable
         call probdraw(temp1Dmat,ndis,CSnode,val) 
         write(*,*) 'Element: CS test'      
         CSazi = chamat(CSnode,9,itime)
         curvature = chamat(CSnode,6,itime)
         if(curvature.gt.0.0) then
          CSazi = CSazi-90.0
         else
          CSazi = CSazi+90.0
         end if
         CHdepth = chamat(CSnode,5,itime)*2.0*CHdwratio
         do iCSnumlobe = 1, int(CSnumlobe)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOLL = max(0.001,xp*stdevCSLOLL+mCSLOLL)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOWW = max(0.001,xp*stdevCSLOWW+mCSLOWW)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOw = max(0.001,xp*stdevCSLOw+mCSLOw)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOl = max(0.001,xp*stdevCSLOl+mCSLOl)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLO_hwratio = max(0.001,xp*stdevCSLO_hwratio+mCSLO_hwratio)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLO_dwratio = max(0.001,xp*stdevCSLO_dwratio+mCSLO_dwratio)
          call avulsioninside(CSnode,CSsinu,CSazi,CSLOLL*1.5,itime,
     +     CHCStemp) ! set the branch length for stop criteria
          regrid = 0
          extrapolate = 0
          call curvature2(CHCStemp,extrapolate,regrid) ! CS branches use temp 
          if(kCS_local.eq.1.1) then
           continue
          end if
          call calc_lobe(CSnode,CSLOLL,CSLOWW,CSLOl,CSLOw,
     +     CSLO_hwratio,CSLO_dwratio,CHCStemp,kCS_local)
         end do
        end do
        continue
       end if

c
c Generate LV
c     

       if(LVwidth.gt.0.0.and.LVdepth+LVheight.gt.0.0) then       
        write(*,*) '  Element: LV'
        nstart=1
        ndis = chaprop(15,itime)
        call calc_levee(nstart,ndis,itime,kLV_local)
       end if       

c
c Generate CH
c

       write(*,*) '  Element: CH'
       code = chaprop(22,itime)
       if(code.eq.KFFCH) then
            FFCHprop = chaprop(13,itime)
            call genabandonedchannel(FFCHprop,itime,code+color_current,KCH+color_current)
       else
            call genchannel(itime,code+color_current)
       end if
c 
c Adjust the color scales
c
       color_current = color_current+color_incr
       if(color_current.gt.0.5) then
         color_current = 0.0
       end if
       KFFCH_local = KFFCH+color_current
       KCS_local = KCS+color_current
       KLV_local = KLV+color_current
       KLA_local = KLA+color_current
       KCH_local = KCH+color_current
      end do !end construct each channel  

c
c Return to the main program:
c

      continue
      return
      end


      subroutine constructmodelAll_noColorAdjust
    
c-----------------------------------------------------------------------
c
c     To construct streamline and facies models
c     ******************************************
c
c The subroutine is called to contruct overall streamline and facies 
c models without color scale adjustment.
c
c
c Author: Fueng Zabel                                    DATE: 2004-2005
c-----------------------------------------------------------------------
      use       streamsimMod
      implicit none

c
c General Variables
c
      integer   idum,ierr,ilevel,idis,ndis,nstart,regrid,extrapolate
      real      val
      real      CHdepth
      real      CSnum,CSnumlobe,CSsinu
      integer   iCSnum,CSnode,iCSnumlobe
      real      curvature,CSazi,CSLOLL,CSLOWW,CSLOw,CSLOl
      real      CSLO_dwratio,CSLO_hwratio,xp
      real*8    p,acorni
      real code
c
c Conditioning 
c
      integer   itime
c
c Streamline Parameters
c
      real   CHdwratio
      real   LVwidth,LVheight,LVdepth,LVasym,LVthin
      real   FFCHprop

c
c Reset the architectural element model
c

      facies = kFF



       kFFCH_local = kFFCH
      kFF_local = kFF
      kCS_local = kCS
      kLV_local = kLV
      kLA_local = kLA
      kCH_local = kCH

      do itime = 1, CHcurrent-1
       write(*,*) 'Constructing Stremline #', itime
       CHdepth = chaprop(2,itime)
       CHdwratio = chaprop(3,itime)
       LVdepth = chaprop(5,itime)
       LVwidth = chaprop(6,itime)
       LVheight = chaprop(7,itime)
       LVasym = chaprop(8,itime)
       LVthin = chaprop(9,itime)
       CSnum = chaprop(11,itime)
       CSnumlobe = chaprop(14,itime)
       CSsource = chaprop(12,itime)

c Generate the architectural elements
c
c
c Generate CS again / should have minimal impact on NTG 
c 
       write(*,*) '  Element: CS'
       CSsinu = 1.2 ! currently hard coded as a low sinuosity
       if(CSnum.gt.0.0) then
        do iCSnum = 1, int(CSnum)
         ndis = chaprop(15,itime)
         do idis = 1, ndis
          temp1Dmat(idis) = abs(chamat(idis,6,itime))
         end do
         CSnode = 1  ! initialize variable
         call probdraw(temp1Dmat,ndis,CSnode,val) 
         write(*,*) 'Element: CS test'      
         CSazi = chamat(CSnode,9,itime)
         curvature = chamat(CSnode,6,itime)
         if(curvature.gt.0.0) then
          CSazi = CSazi-90.0
         else
          CSazi = CSazi+90.0
         end if
         CHdepth = chamat(CSnode,5,itime)*2.0*CHdwratio
         do iCSnumlobe = 1, int(CSnumlobe)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOLL = max(0.001,xp*stdevCSLOLL+mCSLOLL)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOWW = max(0.001,xp*stdevCSLOWW+mCSLOWW)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOw = max(0.001,xp*stdevCSLOw+mCSLOw)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLOl = max(0.001,xp*stdevCSLOl+mCSLOl)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLO_hwratio = max(0.001,xp*stdevCSLO_hwratio+mCSLO_hwratio)
          p = real(acorni(idum))
          call gauinv(p,xp,ierr)
          CSLO_dwratio = max(0.001,xp*stdevCSLO_dwratio+mCSLO_dwratio)
          call avulsioninside(CSnode,CSsinu,CSazi,CSLOLL*1.5,itime,
     +     CHCStemp) ! set the branch length for stop criteria
          regrid = 0
          extrapolate = 0
          call curvature2(CHCStemp,extrapolate,regrid) ! CS branches use temp 
          if(kCS_local.eq.1.1) then
           continue
          end if
          call calc_lobe(CSnode,CSLOLL,CSLOWW,CSLOl,CSLOw,
     +     CSLO_hwratio,CSLO_dwratio,CHCStemp,kCS_local)
         end do
        end do
        continue
       end if

c
c Generate LV
c    

       if(LVwidth.gt.0.0.and.LVdepth+LVheight.gt.0.0) then       
        write(*,*) '  Element: LV'
        nstart=1
        ndis = chaprop(15,itime)
        call calc_levee(nstart,ndis,itime,kLV_local)
       end if       

c
c Generate CH
c
c

       write(*,*) '  Element: CH'
       code = chaprop(22,itime)
       if(code.eq.KFFCH) then
            FFCHprop = chaprop(13,itime)
            call genabandonedchannel(FFCHprop,itime,code,KCH)
       else
            call genchannel(itime,code)
       end if
      end do !end construct each channel  

c
c Return to the main program:
c

      continue
      return
      end


      subroutine copystream(time1,time2)
    
c-----------------------------------------------------------------------
c
c     Make a copy of Streamline
c     *************************
c
c Copies a streamline
c
c
c INPUT VARIABLES:
c
c   ndis             the current number of streamline discretizations
c   chamat           the streamline location, width and thalweg location
c
c OUTPUT VARIABLES:
c
c   facies           a realistic channel geometry fit to the streamline.
c
c
c
c EXTERNAL REFERENCES: azimuth    computes the azimuth
c                      spline     calculates spline
c                      splint     spline interpolation 
c
c Created by: Michael J. Pyrcz                        March, 2004
c
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none

c
c Local Variables
c
 
      integer idis,idata

c
c Looked Up Variables
c 
     
      integer   ndis,spline_count


c
c Passed Variables
c

      integer   time1,time2  

c
c Look up ndis and spline_count from lookup table
c

      write(*,*) ' Started: Copy streamline ',time1,' to ',time2
      ndis=chaprop(15,time1)
      spline_count=chaprop(16,time1)

c
c Copy the channel property matrix 
c

      do idata = 1, ndata_chaprop
       chaprop(idata,time2)=chaprop(idata,time1)
      end do

c
c Copy the spline look up table entry
c
   
      do idata = 1, ndata_spline_table
       do idis = 1, spline_count
        spline_table(idis,idata,time2)=spline_table(idis,idata,time1)
        spline_table2(idis,idata,time2)=spline_table2(idis,idata,time1)
       end do
      end do

c
c Copy the control node matrix
c

      do idata = 1, ndata_chamat
       do idis = 1, ndis
        chamat(idis,idata,time2) = chamat(idis,idata,time1)
       end do
      end do

c
c Finished:
c

      write(*,*) ' Finished: Copy streamline ',time1,' to ',time2
      return
      end


      subroutine curvature2(ctime,cextrapolate,cregrid)
c-----------------------------------------------------------------------
c
c Second attempt to get stable results (regrid and then get properties)      
c
c Caulculate the curvature, dCs/ds, thalweg, and construct splines and 
c regrid to equal spacing with the original number of discretizations.
c
c
c INPUT VARIABLES:
c
c   ndis0            the original number of streamline discretizations
c   chamat           the current streamline location
c
c OUTPUT VARIABLES:
c
c   chamat           the regridded channel streamline with curvature etc.
c                    and regridded
c
c
c
c EXTERNAL REFERENCES: azimuth    computes the azimuth
c                      spline     calculates spline
c                      splint     spline interpolation 
c
c Author: Michael Pyrcz                                  DATE: 2003-2004            
c-----------------------------------------------------------------------
 
       use       nr, only : spline,splint
      use       streamsimMod
      implicit none

c
c Local Variables
c

      integer   idis,ispline,spline_count,nwin
      real      x,y
      real      x1,x2,y1,y2,x3,y3,cx,cy
      real      dist,dist0,distance,ds
      real      yp1x,ypnx,yp1y,ypny
      real      yxA,yyA,yxB,yyB
      real      azi,dazi,dazi1,dazi2,azi1,azi2,maxcurve
      real      curv1,curv2
      real      length,length1A,length2A,length1B,length2B
      real      CHarea,b,mCHhalfwidth
      real      thalweg,width,tmin
      real      avgCHelev
      integer   ix,iy
 
c
c Looked Up Variables
c 
     
      integer   ndis
      real      CHdepth,CHdwratio,maxCHhalfwidth

c
c Passed Variables
c

      integer   ctime,cextrapolate,cregrid

c
c Begin
c

      write(*,*) ' Started: Curvature streamline ',ctime
      tmin = -900.0
c     nwin = 10
      nwin = int((real(nx+ny)/2.0)/10.0) ! smoothing window for porperties

      if(ctime.eq.10) then
       continue
      end if


c
c Look up the current ndis in the chamat matrix
c

      ndis = chaprop(15,ctime)
      CHdepth = chaprop(2,ctime) 
      CHdwratio = chaprop(3,ctime)

c
c Find the maximum channel half width
c
      maxCHhalfwidth = 0.0
      do idis = 1, ndis
       maxCHhalfwidth=max(maxCHhalfwidth,chamat(idis,5,ctime))
      end do
      chaprop(17,ctime) = maxCHhalfwidth

c
c Calculate segment lengths (points have migrated since last gridding)
c       
             
      spline_count = 1      
      length = 0.0
      spline_s(spline_count) = 0.00
      spline_x(spline_count) = xmin
      spline_y(spline_count) = chamat(1,2,ctime)
      spline_w(spline_count) = chamat(1,5,ctime)
      spline_z(spline_count) = chamat(1,11,ctime)

      do idis = 2, ndis
       spline_count = spline_count+1 
       x1 = chamat(idis-1,1,ctime)
       y1 = chamat(idis-1,2,ctime)
       x2 = chamat(idis,1,ctime)
       y2 = chamat(idis,2,ctime)     
       length = length + sqrt((x2-x1)**2.0+(y2-y1)**2.0)
       spline_s(spline_count) = length
       spline_x(spline_count) = chamat(idis,1,ctime)
       spline_y(spline_count) = chamat(idis,2,ctime)
       spline_w(spline_count) = chamat(idis,5,ctime)
       spline_z(spline_count) = chamat(idis,11,ctime)
       if(x2.gt.xmax) goto 111 
      end do 

c
c Always extrapolate to the right hand side - more stable / less flexible
c

      if(cextrapolate.eq.1) then
       x3 = xmax+0.1*(xmax-xmin)
       y3 = y2
       spline_count = spline_count + 1
       spline_x(spline_count) = x3
       spline_y(spline_count) = y3
       spline_w(spline_count) = spline_w(spline_count-1)
          spline_z(spline_count) = spline_z(spline_count-1)  
       spline_s(spline_count) = length + sqrt((x3-x2)**2.0+(y3-y2)**2.0)
      end if

c
c Set the boundary conditions
c

      if(cregrid.ne.1) then
       continue
      end if

 
111   yp1x = 0.0
      ypnx = 0.0      
      yp1y = 0.0
      ypny = 0.0

c
c Calculate the azimuth only
c

      do ispline = 2, spline_count
            yxA = spline_x(ispline-1)
       yxB = spline_x(ispline)
       yyA = spline_y(ispline-1)
       yyB = spline_y(ispline)
       call azimuth(yxA,yxB,yyA,yyB,azi)
       spline_i(ispline) = azi
c      write(*,*) azi
      end do
      spline_i(1) = spline_i(2)

      call movwinsmooth_azimuth(spline_count,nwin,tmin,spline_i)
      continue

c
c Calculate the curvature
c
 
      do ispline = 2, spline_count
       dist = spline_s(ispline)-spline_s(ispline-1)
       azi1 = spline_i(ispline-1)
       azi2 = spline_i(ispline)
        if(abs(azi2-azi1).lt.abs(azi2-(azi1+360.0))) then
        dazi = azi2-azi1
       else 
        dazi = azi2-(azi1+360.0)
       end if
       spline_c(ispline) = dazi/dist 
      end do
      spline_c(1) = spline_c(2) 
      nwin = 10
      call movwinsmooth(spline_count,nwin,tmin,spline_c)
      continue

c
c Calculate the dCsi/ds 
c

      do ispline = 2, spline_count
       dist = spline_s(ispline)-spline_s(ispline-1)
       curv1 = spline_c(ispline-1)
       curv2 = spline_c(ispline)
       spline_d(ispline) = (curv2-curv1)/dist
      end do
      spline_d(1)=spline_d(2)
      nwin = 10
      call movwinsmooth(spline_count,nwin,tmin,spline_d)
      continue

c
c Calculate the relative thalweg along the channel
c

      maxcurve = 0.0
      do ispline = 2, spline_count
       maxcurve = max(abs(spline_c(ispline)),maxcurve)
      end do
      chaprop(18,ctime)=maxcurve
      do ispline = 1, spline_count
       if(spline_c(ispline).lt.0.0) then
        spline_t(ispline) = 0.5-0.25*abs(spline_c(ispline))/
     +  maxcurve
       else if(spline_c(ispline).ge.0.0) then
        spline_t(ispline) = 0.5+0.25*abs(spline_c(ispline))/
     +  maxcurve
       else
        spline_t(ispline) = 0.5
       end if
      end do

c
c Calculate channel depth accounting for conservation of hydraulic area (flow)
c
  
      mCHhalfwidth = (CHdepth/CHdwratio)/2.0
      b = 1 ! thalweg in ceter of channel (a=0.5)
      CHarea = 4.0*CHdepth*mCHhalfwidth*2.0*
     + (1.0/(b+1.0)-1.0/(2.0*b+1.0))
      chaprop(19,ctime) = CHarea
      do ispline = 1, spline_count
       width = 2.0*spline_w(ispline)
       thalweg = spline_t(ispline)
       if(thalweg.gt.0.5) then
        thalweg = 1.0 - thalweg ! area calculation is a symmetric problem
       end if
       b = -log(2.0)/log(thalweg)
       spline_h(ispline)=Charea/(4.0*width* ! eqn from integration of channel eqn
     +  (1.0/(b+1.0)-1.0/(2.0*b+1.0))) 
      end do 

c
c Set all splines
c

      call spline(spline_s,spline_x,yp1x,ypnx,spline_2x,spline_count)
      call spline(spline_s,spline_y,yp1y,ypny,spline_2y,spline_count)
      call spline(spline_s,spline_w,yp1x,ypnx,spline_2w,spline_count)
      call spline(spline_s,spline_t,yp1x,ypnx,spline_2t,spline_count)
      call spline(spline_s,spline_c,yp1x,ypnx,spline_2c,spline_count)
      call spline(spline_s,spline_d,yp1x,ypnx,spline_2d,spline_count)
      call spline(spline_s,spline_i,yp1x,ypnx,spline_2i,spline_count)
      call spline(spline_s,spline_z,yp1x,ypnx,spline_2z,spline_count)
      call spline(spline_s,spline_h,yp1x,ypnx,spline_2h,spline_count)

c
c Update the table entry with the number of discritizations
c      
      
      chaprop(16,ctime)=spline_count

c
c Copy the spline vectors to the spline table
c

      do ispline = 1, spline_count
       spline_table(ispline,1,ctime)=spline_x(ispline)
       spline_table(ispline,2,ctime)=spline_y(ispline)
       spline_table(ispline,3,ctime)=spline_s(ispline)
       spline_table(ispline,4,ctime)=spline_t(ispline)
       spline_table(ispline,5,ctime)=spline_w(ispline)
       spline_table(ispline,6,ctime)=spline_c(ispline)
       spline_table(ispline,8,ctime)=spline_d(ispline)
       spline_table(ispline,9,ctime)=spline_i(ispline)
       spline_table(ispline,11,ctime)=spline_z(ispline)
       spline_table(ispline,12,ctime)=spline_h(ispline)                   
      end do

      do ispline = 1, spline_count
       spline_table2(ispline,1,ctime)=spline_2x(ispline)
       spline_table2(ispline,2,ctime)=spline_2y(ispline)
       spline_table2(ispline,4,ctime)=spline_2t(ispline)
       spline_table2(ispline,5,ctime)=spline_2w(ispline)
       spline_table2(ispline,6,ctime)=spline_2c(ispline)
       spline_table2(ispline,8,ctime)=spline_2d(ispline)
       spline_table2(ispline,9,ctime)=spline_2i(ispline)
       spline_table2(ispline,11,ctime)=spline_2z(ispline)            
       spline_table2(ispline,12,ctime)=spline_2h(ispline)            
      end do


c
c Regrid the discretized chamat arrays 
c
  
      distance = 0.0
      ds = 0.0
      if(cregrid.eq.1) then
       
c
c Reset ndis in the chamat lookup table
c

       ndis = ndis0
       chaprop(15,ctime) = ndis
       dist0 = 0.0
       do idis = 1, ndis
        distance = real(idis-1)*(spline_s(spline_count)/real(ndis-1))
        chamat(idis,1,ctime) = splint(spline_s,spline_x,spline_2x,distance,spline_count)       
        chamat(idis,2,ctime) = splint(spline_s,spline_y,spline_2y,distance,spline_count)
        chamat(idis,3,ctime) = distance
        chamat(idis,4,ctime) = splint(spline_s,spline_t,spline_2t,distance,spline_count)
        chamat(idis,5,ctime) = splint(spline_s,spline_w,spline_2w,distance,spline_count)
        chamat(idis,6,ctime) = splint(spline_s,spline_c,spline_2c,distance,spline_count)
        chamat(idis,7,ctime) = distance-dist0
        chamat(idis,8,ctime) = splint(spline_s,spline_d,spline_2d,distance,spline_count)
        chamat(idis,9,ctime) = splint(spline_s,spline_i,spline_2i,distance,spline_count)
        chamat(idis,11,ctime) = splint(spline_s,spline_z,spline_2z,distance,spline_count)
        chamat(idis,12,ctime) = splint(spline_s,spline_h,spline_2h,distance,spline_count)
        dist0 = distance
       end do
      else  
       ndis = spline_count 
       chaprop(15,ctime) = ndis           
       cx = spline_x(1)
       cy = spline_y(1)
       do idis = 1, ndis
        x = spline_x(idis)
        y = spline_y(idis)
        ds = sqrt((x-cx)**2.0+(y-cy)**2.0)
        distance = distance + ds
        chamat(idis,1,ctime) = spline_x(idis)       
        chamat(idis,2,ctime) = spline_y(idis)
        chamat(idis,3,ctime) = distance
        chamat(idis,4,ctime) = spline_t(idis)
        chamat(idis,5,ctime) = spline_w(idis)
        chamat(idis,6,ctime) = spline_c(idis)
        chamat(idis,7,ctime) = ds
        chamat(idis,8,ctime) = spline_d(idis)
        chamat(idis,9,ctime) = spline_i(idis)
        chamat(idis,11,ctime) = spline_z(idis)
        chamat(idis,12,ctime) = spline_h(idis)
        cx = x
        cy = y 
       end do
      end if

      if(write_surfaces) then

c
c Calculate the average CH elevation
c
 
       avgCHelev = 0.0
       do idis = 1, ndis
        avgChelev = avgCHelev + chamat(idis,11,ctime)
       end do
       avgCHelev = avgCHelev/real(ndis)
       do iy = 1, ny
        do ix = 1, nx        
         ihs(ix,iy,ctime,1) = avgCHelev
         ihs(ix,iy,ctime,2) = 1
        end do
       end do
      end if

c
c Finished:
c

      write(*,*) ' Finished: Curvature event centerline ',ctime    
      return
      end



      subroutine genabandonedchannel(mud_prop,ctime,lkFFCH,lkCH)
    
c-----------------------------------------------------------------------
c
c     Generate the channel facies model
c     *********************************
c
c Generate a realistic channel cross section based on Deutsch and Wang.
c Reuires a channel spline.
c
c INPUT VARIABLES:
c
c   ndis             the current number of streamline discretizations
c   chamat           the streamline location, width and thalweg location
c
c GLOBAL INPUT VARIABLES
c
c
c
c OUTPUT VARIABLES:
c
c   facies           a realistic channel geometry fit to the streamline.
c
c
c
c EXTERNAL REFERENCES: azimuth    computes the azimuth
c                      spline     calculates spline
c                      splint     spline interpolation 
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      
      use       nr, only : splint
      use       streamsimMod
      implicit none

c
c Local Variables
c
 
      integer   i,ix,iy,iz,idis,jdis 
      real      x,y,z
      real      azi1,azi2,dazi
      real      sds,eds,ds0,sx,sx0,sy,sy0,ds
      real      dist,dist0,cdis,distance
      real      CHelev,CHhalfwidth,CHdepth,CHbot,FFCHbot
      real      w,WW,by,cy,a

c
c Looked Up Variables
c 

      integer   ndis,spline_count
      real      maxCHhalfwidth,CHdwratio
     
c
c Passed Variables
c

      integer   ctime
      real      lkCH,lkFFCH,mud_prop
   
c
c Begin
c

      write(*,*) ' Started: FF(CH) streamline ',ctime
      FFCHcount = 0
      schedule(ctime)=1 ! indicate channel abandonment in schedule

c
c Copy the spline arrays into the required vectors
c

c
c Get the number of spline discretizations
c      
      
      ndis = chaprop(15,ctime)
      spline_count=chaprop(16,ctime)
      maxCHhalfwidth=chaprop(17,ctime)
      CHdwratio = chaprop(3,ctime)

c
c Copy the spline vectors to the spline table
c

      do idis = 1, spline_count
       spline_x(idis)=spline_table(idis,1,ctime)
       spline_y(idis)=spline_table(idis,2,ctime)
       spline_s(idis)=spline_table(idis,3,ctime)
       spline_t(idis)=spline_table(idis,4,ctime)
       spline_w(idis)=spline_table(idis,5,ctime)
       spline_c(idis)=spline_table(idis,6,ctime)
       spline_d(idis)=spline_table(idis,8,ctime)
       spline_i(idis)=spline_table(idis,9,ctime)
       spline_z(idis)=spline_table(idis,11,ctime)            
      end do
      do idis = 1, spline_count
       spline_2x(idis)=spline_table2(idis,1,ctime)
       spline_2y(idis)=spline_table2(idis,2,ctime)
       spline_2t(idis)=spline_table2(idis,4,ctime)
       spline_2w(idis)=spline_table2(idis,5,ctime)
       spline_2c(idis)=spline_table2(idis,6,ctime)
       spline_2d(idis)=spline_table2(idis,8,ctime)
       spline_2i(idis)=spline_table2(idis,9,ctime)
       spline_2z(idis)=spline_table2(idis,11,ctime)            
      end do

c
c Generate the channel object
c

      if(ctime.eq.0) then
       continue
      end if

      do iy = 1, ny
       y = ymin + real(iy)*ysiz
       do ix = 1, nx
        x = xmin + real(ix)*xsiz        
        distance = big_number

c
c Find the nearest control node
c

        do idis = 1, ndis
         dist = (chamat(idis,1,ctime)-x)**2+
     +          (chamat(idis,2,ctime)-y)**2
         if(dist.lt.distance) then
          cdis = idis
          distance = dist
         end if
        end do
        distance = sqrt(distance)

c
c Refine search: Find the nearest location on the spline
c

        if(distance.lt.maxCHhalfwidth*2.0) then
         jdis = max(1,int(cdis-1))
         sds = chamat(jdis,3,ctime)
         jdis = min(ndis,int(cdis+1))
         eds = chamat(jdis,3,ctime)
         dist = big_number
         do i = 1, ndiscr
          ds0 = sds + (real(i)/real(ndiscr))*(eds-sds)
          sx = splint(spline_s,spline_x,spline_2x,ds0,spline_count)
          sy = splint(spline_s,spline_y,spline_2y,ds0,spline_count)
          dist0 = (sx-x)**2+(sy-y)**2
          if(dist0.lt.dist) then
           ds = ds0
           dist = dist0
          end if
         end do
711         dist = sqrt(dist)
         sx = splint(spline_s,spline_x,spline_2x,ds,spline_count)
         sy = splint(spline_s,spline_y,spline_2y,ds,spline_count)
         CHelev = splint(spline_s,spline_z,spline_2z,ds,spline_count)
         CHhalfwidth = splint(spline_s,spline_w,spline_2w,ds,spline_count)
         a = splint(spline_s,spline_t,spline_2t,ds,spline_count)
         sx0 = splint(spline_s,spline_x,spline_2x,ds-0.001,spline_count)
         sy0 = splint(spline_s,spline_y,spline_2y,ds-0.001,spline_count)
         call azimuth(sx0,sx,sy0,sy,azi1)
         call azimuth(sx,x,sy,y,azi2)

c
c Improved left / right determination
c
   
         dazi = azi2-azi1
         if(abs(dazi).gt.180.0) then
          dazi = (azi1-180.0)-azi2
         end if

         CHdepth = CHdwratio*CHhalfwidth*2.0
         WW = CHhalfwidth*2.0 
         if(dazi.lt.0.0) then ! left side
          w = -1.0*dist+CHhalfwidth
         else ! right side
          w = dist+CHhalfwidth
         end if
         a = min(0.999,max(0.001,a))
         if(w.lt.0.0.or.w.gt.WW) goto 712
         if(a.lt.0.5) then
          by = -1.0*log(2.0)/log(a)
          CHbot = CHelev - 4.0*CHdepth*(w/WW)**by*(1.0-(w/WW)**by)
          continue
         else
          if(a.gt.1.0) then
           continue
          end if
          cy = -1.0*log(2.0)/log(1.0-a)
          CHbot = CHelev - 4.0*CHdepth*(1.0-w/WW)**cy*
     +       (1.0-(1.0-w/WW)**cy)
         end if
         FFCHbot = CHelev-(CHelev-CHbot)*mud_prop

c
c IHS accretionary Surfaces
c

         if(write_surfaces) then
          ihs(ix,iy,ctime,1) = chbot
          ihs(ix,iy,ctime,2) = 0
         end if
         do iz = 1, nz
          z = zmin + real(iz)*zsiz
          if(z.ge.CHbot.and.z.le.CHelev) then
           if(z.gt.FFCHbot) then ! abandoned channel
            FFCHcount = FFCHcount + 1
            if(facies(ix,iy,iz).ge.1.0) then
             NTGcurrent = NTGcurrent-1
            end if
            facies(ix,iy,iz) = lkFFCH
           else                  ! channel fill              
            if(facies(ix,iy,iz).lt.1.0) then
             NTGcurrent = NTGcurrent+1
            end if 
            facies(ix,iy,iz) = lkCH 
           end if ! decide between CH and FFCH
          end if  ! if within channel
         end do   ! loop over all iz       
712     end if    ! check if distance within maxCHhalfwidth 
       end do     ! loop over all ix
      end do      ! loop over all iy
   
      continue                

c
c Finished:
c

      write(*,*) ' Finished: FF(CH) streamline ',ctime
      return
      end



      subroutine genchannel(ctime,k)
    
c-----------------------------------------------------------------------
c
c     Generate the channel facies model
c     *********************************
c
c Generate a realistic channel cross section based on Deutsch and Wang.
c Reuires a channel spline.
c
c INPUT VARIABLES:
c
c   ndis             the current number of streamline discretizations
c   chamat           the streamline location, width and thalweg location
c
c OUTPUT VARIABLES:
c
c   facies           a realistic channel geometry fit to the streamline.
c
c
c
c EXTERNAL REFERENCES: azimuth    computes the azimuth
c                      spline     calculates spline
c                      splint     spline interpolation 
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      
      use     nr, only : splint
      use     streamsimMod
      implicit none

c
c Local Variables
c
 
      integer ix,iy,iz,idis,jdis,ifac,i,cdis
      real    x,y,z,sx,sy,sx0,sy0,dist,dist0,distance,ds,sds,eds,ds0
      real    azi1,azi2,dazi
      real    w,WW,a,by,cy
      real    CHelev,CHhalfwidth,CHdepth,CHbot,FFCHbot

c
c Looked Up Variables
c 
            
      integer spline_count,ndis
      real    maxCHhalfwidth,CHdwratio,FFCHprop

c
c Passed Variables
c

      integer ctime
      real    k

c
c Begin 
c

      write(*,*) ' Started: CH streamline ',ctime
      FFCHcount = 0

cccc
c
c Copy the spline arrays into the required vectors
c
cccc

c
c Get the number of spline discretizations
c      
      
      CHdepth=chaprop(2,ctime)
      CHdwratio=chaprop(3,ctime)
      ndis=chaprop(15,ctime)
      spline_count=chaprop(16,ctime)
      maxCHhalfwidth=chaprop(17,ctime)
      FFCHprop = chaprop(13,ctime)

c
c Copy the spline vectors to the spline table
c

      do idis = 1, spline_count
       spline_x(idis)=spline_table(idis,1,ctime)
       spline_y(idis)=spline_table(idis,2,ctime)
       spline_s(idis)=spline_table(idis,3,ctime)
       spline_t(idis)=spline_table(idis,4,ctime)
       spline_w(idis)=spline_table(idis,5,ctime)
       spline_c(idis)=spline_table(idis,6,ctime)
       spline_d(idis)=spline_table(idis,8,ctime)
       spline_i(idis)=spline_table(idis,9,ctime)
       spline_z(idis)=spline_table(idis,11,ctime)            
      end do
      do idis = 1, spline_count
       spline_2x(idis)=spline_table2(idis,1,ctime)
       spline_2y(idis)=spline_table2(idis,2,ctime)
       spline_2t(idis)=spline_table2(idis,4,ctime)
       spline_2w(idis)=spline_table2(idis,5,ctime)
       spline_2c(idis)=spline_table2(idis,6,ctime)
       spline_2d(idis)=spline_table2(idis,8,ctime)
       spline_2i(idis)=spline_table2(idis,9,ctime)
       spline_2z(idis)=spline_table2(idis,11,ctime)            
      end do

c
c Generate the channel object
c

      if(ctime.eq.0) then
       continue
      end if

      do iy = 1, ny
       y = ymin + real(iy)*ysiz
       do ix = 1, nx
        x = xmin + real(ix)*xsiz        
        distance = big_number

c
c Find the nearest control node
c

        do idis = 1, ndis
         dist = (chamat(idis,1,ctime)-x)**2+
     +          (chamat(idis,2,ctime)-y)**2
         if(dist.lt.distance) then
          cdis = idis
          distance = dist
         end if
        end do
      
        distance = sqrt(distance)

c
c Refine search: Find the nearest location on the spline
c

        if(distance.lt.maxCHhalfwidth*2.0) then
         jdis = max(1,int(cdis-1))
         sds = chamat(jdis,3,ctime)
         jdis = min(ndis,int(cdis+1))
         eds = chamat(jdis,3,ctime)
         dist = big_number
         do i = 1, ndiscr
          ds0 = sds + (real(i)/real(ndiscr))*(eds-sds)
          sx = splint(spline_s,spline_x,spline_2x,ds0,spline_count)
          sy = splint(spline_s,spline_y,spline_2y,ds0,spline_count)
          dist0 = (sx-x)**2+(sy-y)**2
          if(dist0.lt.dist) then
           ds = ds0
           dist = dist0
c           goto 711
          end if
         end do
711         dist = sqrt(dist)
         sx = splint(spline_s,spline_x,spline_2x,ds,spline_count)
         sy = splint(spline_s,spline_y,spline_2y,ds,spline_count)
         chelev = splint(spline_s,spline_z,spline_2z,ds,spline_count)
         CHhalfwidth = splint(spline_s,spline_w,spline_2w,ds,spline_count)
         a = splint(spline_s,spline_t,spline_2t,ds,spline_count)
         sx0 = splint(spline_s,spline_x,spline_2x,ds-0.001,spline_count)
         sy0 = splint(spline_s,spline_y,spline_2y,ds-0.001,spline_count)
         call azimuth(sx0,sx,sy0,sy,azi1)
         call azimuth(sx,x,sy,y,azi2)

c
c Improved left / right determination
c
   
         dazi = azi2-azi1
         if(abs(dazi).gt.180.0) then
          dazi = (azi1-180.0)-azi2
         end if

         CHdepth = CHdwratio*CHhalfwidth*2.0
         WW = CHhalfwidth*2.0 
         if(dazi.lt.0.0) then ! left side
          w = -1.0*dist+CHhalfwidth
         else ! right side
          w = dist+CHhalfwidth
         end if
         a = min(0.999,max(0.001,a))
         if(w.lt.0.0.or.w.gt.WW) goto 712
         if(a.lt.0.5) then
          by = -1.0*log(2.0)/log(a)
          chbot = chelev - 4.0*CHdepth*(w/WW)**by*(1.0-(w/WW)**by)
          continue
         else
          if(a.gt.1.0) then
           continue
          end if
          cy = -1.0*log(2.0)/log(1.0-a)
          chbot = chelev - 4.0*CHdepth*(1.0-w/WW)**cy*
     +     (1.0-(1.0-w/WW)**cy)
         end if
         FFCHbot = CHelev-(CHelev-CHbot)*FFCHprop
         if(write_surfaces) then
          ihs(ix,iy,ctime,1) = chbot
          ihs(ix,iy,ctime,2) = 0
         end if
         do iz = 1, nz
          z = zmin + real(iz)*zsiz
          
c
c Erode "material" directly above the channel
c
            
               if(z.gt.chelev) then ! CH erodes residual LV in migration
           if(facies(ix,iy,iz).ge.1.0) then
            NTGcurrent=NTGcurrent-1
            facies(ix,iy,iz)=kFF
           end if
          end if
          if(z.gt.chbot.and.z.lt.chelev) then
           if(facies(ix,iy,iz).lt.1.0) then
            NTGcurrent=NTGcurrent+1
           end if
           facies(ix,iy,iz) = k
c
c Count potential FFCH if the channel avulsed 
c

           if(z.gt.FFCHbot) then
            FFCHcount=FFCHcount+1
           end if
            end if
         end do ! end do iz        
712     end if
       end do ! end do iy
      end do !end do ix  
   
      continue                   

c
c Finished:
c
 
      write(*,*) ' Finished: CH streamline ',ctime
      return
      end


      subroutine genchannellocal(rix,riy,rtime,rclosenode,riztop,
     +      rizbot)
    
c-----------------------------------------------------------------------
c
c     Generate the channel facies model
c     *********************************
c
c Generate a realistic channel cross section based on Deutsch and Wang.
c Requires a channel spline....
c
c INPUT VARIABLES:
c
c  rix = x location of well 
c  riy = y location of well 
c  rtime= closest streamline to particular CH interval in a well 
c
c OUTPUT VARIABLES:
c
c  rclosenode = closest node
c  riztop = channel top
c  rizbot = channel bottom
c
c GLOBAL VARIABLES:
c
c chaprop = channel property array 
c
c EXTERNAL REFERENCES: azimuth    computes the azimuth
c                      spline     calculates spline
c                      splint     spline interpolation 
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      
      use       nr, only : splint
      use       streamsimMOD
      implicit none

c
c Local Variables
c

      integer   idis,jdis,ispline,nwin,ni,cdis,inflagztop,inflagzbot
      real      x,y,z,w,c,d,t,i,h,WW
      real      sds,eds,ds0,sx,sy,CHhalfwidth,a,sx0,sy0
      real      by,cy
      real      dist,dist0,distance,ds
      real      yp1x,ypnx,yp1y,ypny
      real      azi,dazi,azi1,azi2
      real      CHelev,Chbot

c
c Looked Up Variables
c 
     
      integer   ndis,spline_count
      real      CHdepth,CHdwratio,maxCHhalfwidth

c
c Passed Variables
c

      integer   rix,riy,rtime,rclosenode,riztop,rizbot

c
c Begin
c  

c
c Check if channel is adandoned 
c

      CHelev = 0.0
      chbot = 0.0
      x = real(rix-1)*xsiz+xmn ! well stored now as ix, iy
      y = real(riy-1)*ysiz+ymn
      ni = 5

cccc
c
c Copy the spline arrays into the required vectors
c
cccc

c
c Get the number of spline discretizations
c      
      
      CHdepth=chaprop(2,rtime)
      CHdwratio=chaprop(3,rtime)
      ndis=chaprop(15,rtime)
      spline_count=chaprop(16,rtime)
      maxCHhalfwidth=chaprop(17,rtime)

c
c Copy the spline vectors to the spline table
c

      do idis = 1, spline_count
       spline_x(idis)=spline_table(idis,1,rtime)
       spline_y(idis)=spline_table(idis,2,rtime)
       spline_s(idis)=spline_table(idis,3,rtime)
       spline_t(idis)=spline_table(idis,4,rtime)
       spline_w(idis)=spline_table(idis,5,rtime)
       spline_c(idis)=spline_table(idis,6,rtime)
       spline_d(idis)=spline_table(idis,8,rtime)
       spline_i(idis)=spline_table(idis,9,rtime)
       spline_z(idis)=spline_table(idis,11,rtime)            
      end do
      do idis = 1, spline_count
       spline_2x(idis)=spline_table2(idis,1,rtime)
       spline_2y(idis)=spline_table2(idis,2,rtime)
       spline_2t(idis)=spline_table2(idis,4,rtime)
       spline_2w(idis)=spline_table2(idis,5,rtime)
       spline_2c(idis)=spline_table2(idis,6,rtime)
       spline_2d(idis)=spline_table2(idis,8,rtime)
       spline_2i(idis)=spline_table2(idis,9,rtime)
       spline_2z(idis)=spline_table2(idis,11,rtime)            
      end do

c
c Generate the channel object
c

      distance = big_number

c
c Find the nearest control node
c

      do idis = 1, ndis
       dist = (chamat(idis,1,rtime)-x)**2+
     +        (chamat(idis,2,rtime)-y)**2
       if(dist.lt.distance) then
        cdis = idis
        distance = dist
       end if
      end do
      distance = sqrt(distance)

c
c Refine search: Find the nearest location on the spline
c

      if(distance.lt.maxCHhalfwidth*2.0) then
       jdis = max(1,int(cdis-1))
       sds = chamat(jdis,3,rtime)
       jdis = min(ndis,int(cdis+1))
       eds = chamat(jdis,3,rtime)
       dist = big_number
       do i = 1, ndiscr
        ds0 = sds + (real(i)/real(ndiscr))*(eds-sds)
        sx = splint(spline_s,spline_x,spline_2x,ds0,spline_count)
        sy = splint(spline_s,spline_y,spline_2y,ds0,spline_count)
        dist0 = (sx-x)**2+(sy-y)**2
        if(dist0.lt.dist) then
         ds = ds0
         dist = dist0
c         goto 711
        end if
       end do
711       dist = sqrt(dist)
       sx = splint(spline_s,spline_x,spline_2x,ds,spline_count)
       sy = splint(spline_s,spline_y,spline_2y,ds,spline_count)
       CHelev = splint(spline_s,spline_z,spline_2z,ds,spline_count)
       CHhalfwidth = splint(spline_s,spline_w,spline_2w,ds,spline_count)
       CHdepth = splint(spline_s,spline_h,spline_2h,ds,spline_count)
       a = splint(spline_s,spline_t,spline_2t,ds,spline_count)
       sx0 = splint(spline_s,spline_x,spline_2x,ds-0.001,spline_count)
       sy0 = splint(spline_s,spline_y,spline_2y,ds-0.001,spline_count)
       call azimuth(sx0,sx,sy0,sy,azi1)
       call azimuth(sx,x,sy,y,azi2)

c
c Improved left / right determination
c
       dazi = azi2-azi1
       if(abs(dazi).gt.180.0) then
        dazi = (azi1-180.0)-azi2
       end if

       WW = CHhalfwidth*2.0 
       if(dazi.lt.0.0) then ! left side
        w = -1.0*dist+CHhalfwidth
       else ! right side
        w = dist+CHhalfwidth
       end if
       a = min(real(0.999),max(real(0.001),a))
       if(w.lt.0.0.or.w.gt.WW) then
        CHbot = CHelev
        goto 712
       end if
       if(a.lt.0.5) then
        by = -1.0*log(2.0)/log(a)
        chbot = CHelev - 4.0*CHdepth*(w/WW)**by*(1.0-(w/WW)**by)
        continue
       else
        if(a.gt.1.0) then
         continue
        end if
        cy = -1.0*log(2.0)/log(1.0-a)
        chbot = CHelev - 4.0*CHdepth*(1.0-w/WW)**cy*
     +  (1.0-(1.0-w/WW)**cy)
       end if
      end if

712      continue
      if(CHelev.eq.0.0.or.CHelev-chbot.le.0.0) then
       call getindx(nz,zmn,zsiz,chelev,riztop,inflagztop) 
       call getindx(nz,zmn,zsiz,chbot,rizbot,inflagzbot) 
       rizbot=rizbot+1
      else
       call getindx(nz,zmn,zsiz,CHelev,riztop,inflagztop) 
       call getindx(nz,zmn,zsiz,chbot,rizbot,inflagzbot) 
      end if 
      rclosenode = cdis
      continue                   

c
c Finished:
c
      return
      end



      subroutine locatenode(rx,ry,rz,rtime,rassoc,rdis,rfound)
c-----------------------------------------------------------------------
c
c      Find the nearest node to a data location
c      ****************************************
c
c
c
c GLOBAL VARIABLES:
c
c   ndis             the number of streamline discretizations
c   chamat           the 2D stream line with x in col 1, and y in col 2
c   
c
c OUTPUT VARIABLES:
c
c  
c EXTERNAL REFERENCES:
c
c 
c Author: Michael Pyrcz                                  DATE: 2003-2004  
c-----------------------------------------------------------------------
 
      use       streamsimMod
      implicit none

c
c Local Variables
c

      integer   idis,itime,stime,etime,iassoc
      real      mindist,dist
      real      dz,deltax,deltay,deltaz

c
c Looked Up Variables
c 
     
      integer   ndis

c
c Passed Variables
c

      logical   rfound
      integer   rdis,rtime,rassoc
      real      rx,ry,rz,rztol

c
c Begin
c
      rtime = -1        
      rassoc = -1        
      rdis = -1

c
c Check if the channel is identified
c

      if(rtime.gt.0) then ! case: specific streamline speicified
       stime = rtime
       etime = rtime
      else if(rassoc.gt.0) then ! case: streamline association specified
       stime = assoc_tab(rassoc,2)
       etime = assoc_tab(rassoc,assoc_tab(rassoc,1)+1)
      else ! global search
       stime = 1
       etime =  CHcurrent-1
      end if

c
c Find the nearest node to conditioning data
c

      rfound = .false.
      mindist=big_number
      do itime = stime, etime
       ndis = int(chaprop(15,itime))
       iassoc = int(chaprop(21,itime))
       do idis = 1, ndis
        dz = (rz-chamat(idis,11,itime))
        if(abs(dz).lt.ztol) then
         deltax = (rx-chamat(idis,1,itime))/xanis 
         deltay = (ry-chamat(idis,2,itime))/yanis 
         deltaz = dz/zanis      
         dist = deltax**2+deltay**2+deltaz**2
         if(dist.lt.mindist.and.assoc_cond(iassoc,idis).eq.0) then 
          mindist=dist
          rtime = itime
          rassoc = iassoc
          rdis = idis
          rfound = .true.
         end if
        end if
       end do
      end do

c
c Finished:
c

      return
      end



      subroutine lookupstream(clevel,dtime,ctime)
    
c-----------------------------------------------------------------------
c
c     Copy a Streamline from Look Up Table
c     ************************************
c
c Copies a streamline from the look up table into active arrays.
c
c
c INPUT VARIABLES:
c
c   clevel           the current level for elevation look up   
c   dtime            the streamline look up table entry
c   ctime            the storage location for streamline
c
c OUTPUT VARIABLES:
c
c Author: Michael Pyrcz                                  DATE: 2003-2004
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none

c
c Local Variables
c
 
      integer idis

c
c Looked Up Variables
c 
     
      integer   ndis

c
c Passed Variables
c

      integer   clevel,dtime,ctime  

c
c Look up ndis and spline_count from lookup table
c
      
      write(*,*) ' Started: Look-up streamline ',dtime,' to ',ctime
      ndis = chadrawprop(3,dtime)

c
c Pass drawn streamline properties
c
         
      chaprop(1,ctime)=chadrawprop(1,dtime)  ! azi
      chaprop(4,ctime)=chadrawprop(2,dtime)  ! sinuosity
      chaprop(15,ctime)=chadrawprop(3,dtime) ! ndis
      chaprop(19,ctime)=chadrawprop(5,dtime) ! CHarea
      chaprop(2,ctime)=chadrawprop(6,dtime)  ! CHdepth
      chaprop(3,ctime)=chadrawprop(7,dtime)  ! CHdwratio

c
c Copy drawn streamline into chamat
c

      do idis = 1, ndis
       chamat(idis,1,ctime) = chadraw(idis,1,dtime)
       chamat(idis,2,ctime) = chadraw(idis,2,dtime)
       chamat(idis,5,ctime) = chadraw(idis,3,dtime)
       chamat(idis,11,ctime) = level(clevel,1) ! constant elev 
      end do

c
c Finished:
c

      write(*,*) ' Finished: Look-up streamline ',dtime,' to ',ctime
      return
      end


      subroutine migrate(tmigrate,time1,time2)
    
c-----------------------------------------------------------------------
c
c     Migrate a Streamline with the Bank Retreat Model
c     ************************************************
c
c Generate a realistic channel meander migration step based on the
c procedure from:
c
c Howard, A.D., 1992, Lowland Floodplain Rivers: Geomorphological 
c Perspectives, chapter Modeling channel migration and floodplain 
c sedimentation in meandering streams, pages 1-37. John Whiley and Sons.
c
c Sun, T., Meakin, P. and Jossang, T., 1996, A Simulation Model for 
c Meandering Rivers, Water Resources Research, Vol. 32, No. 9., pages 
c 2937-2954. 
c
c INPUT VARIABLES:
c
c   ndis             the current number of streamline discretizations
c   chamat           the streamline location, width and thalweg location
c
c OUTPUT VARIABLES:
c
c   facies           a realistic channel geometry fit to the streamline.
c
c
c
c EXTERNAL REFERENCES: azimuth    computes the azimuth
c                      spline     calculates spline
c                      splint     spline interpolation 
c
c Author: Michael Pyrcz                                  DATE: 2003-2004
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none

c
c Local Variables
c
 
      integer   idis
      real      x1,y1,x2,y2,ang,dist
      real      maxmigrate,scale_migrate

c
c Looked Up Variables
c 
     
      integer   ndis

c
c Passed Variables
c

      integer   time1,time2
      real      tmigrate       
c
c Look up ndis from lookup table and pass on to the new channel
c

      write(*,*) ' Started: Migrate event ',time1,' to ',time2
      ndis=chaprop(15,time1)
      chaprop(15,time2)=ndis

c
c Restandardize the migration to keep regular per time step
c
      maxmigrate = 0.0
      do idis = 2, ndis
       maxmigrate = max(maxmigrate,abs(chamat(idis,10,time1)))
      end do
      scale_migrate = tmigrate/maxmigrate
      do idis = 2, ndis
       chamat(idis,10,time1) = chamat(idis,10,time1)*scale_migrate
      end do

      write(*,*) 'Migration Coefficient Scalar', scale_migrate

c
c Calculated the by-point migration and apply (used uniform E now)
c

      chamat(1,1,time2) = chamat(1,1,time1)
      chamat(1,2,time2) = chamat(1,2,time1)
      maxmigrate = 0.0
      do idis = 2,ndis
       ang = chamat(idis,9,time1) - 90.0
       x1 = chamat(idis,1,time1)
       y1 = chamat(idis,2,time1)
       dist = chamat(idis,10,time1)
       call offset(x1,y1,ang,dist,x2,y2)
       chamat(idis,1,time2) = x2
       chamat(idis,2,time2) = y2
      end do
      continue                   

c
c Inherrent the width and vertical profile from the last time step
c

      do idis = 1, ndis
       chamat(idis,5,time2)=chamat(idis,5,time1)
       chamat(idis,11,time2)=chamat(idis,11,time1)
      end do

c
c Finished:
c

      write(*,*) ' Finished: Migrate event ',time1,' to ',time2
      return
      end


      subroutine morphCHendpts(rtime,rdis,rix,riy,rztop,rzbot)
    
c-----------------------------------------------------------------------
c
c     Morph a channel belt 
c     ********************
c
c The subroutine morphs the channel streamline and preserves the 
c endpoints to prevent extrapolation problems in curvature2 subroutine.
c
c 
c
c INPUT VARIABLES:
c
c rtime    the current channel belt 
c rdis     the closest node in the CH element of the channel belt 
c rix      the ix location of the conditioning data
c riy      the iy location of the conditioning data
c rztop    the top of the well contact 
c rzbot    the bottom of the well contact
c
c GLOBAL VARIABLE:
C
c   sinu          approximate sinuosity of the initial channel
c
c Author: Michael Pyrcz                                  DATE: 2003-2004
c-----------------------------------------------------------------------
      
      use       streamsimMOD
      implicit none

c
c Local Variables
c

      logical   converge,no_prev_cond
      integer   i,idis,this_cdis,lowcond,highcond,itime,closenode
      integer   inflagz,iztop_target,izbot_target
      integer   ithick,ithick_target,iztop,izbot,ithick_error
      integer   itol_lower,itol_upper,irepo
      integer   extrapolate,regrid
      real      x,y,xchi,ychi,azi,cstep,cstepx,cstepy,x0,y0,x1,y1
      real      init_xCHi,init_yCHi
      real      iterstartx,iterstarty,stepx,stepy,iterstep,iter
      real      deltax,deltay,deltaz,distance,small_number
      real,allocatable :: shiftmat(:)
      real      test3(6,1000)

c
c Looked Up Variables
c 
     
      integer   ndis,cassoc

c
c Passed Variables
c

      integer   rtime,rdis,rix,riy,riztop,rizbot 
      real      rztop,rzbot

c
c Begin
c  

c
c Prevent any regridding / required for conditioning tables
c

      regrid = 0
      extrapolate = 0

c
c Initialize variables
c
 
      small_number = (xsiz+ysiz)/2.0*0.01
      converge = .false.
      x0 = 0.0
      y0 = 0.0
      x1 = 0.0
      y1 = 0.0
      azi = 0.0

c
c Look up the current and target locations
c

      ndis = chaprop(15,rtime)
      if (ndis.lt.ndis0) then
          ndis = ndis0
      end if
      cassoc = chaprop(21,rtime)

c
      allocate(shiftmat(ndis),stat = test)

      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
     
c
c Well x and y
c

      x = real(rix-1)*xsiz+xmn
      y = real(riy-1)*ysiz+ymn

c
c Column of CH cells
c

      call getindx(nz,zmn,zsiz,rztop,iztop_target,inflagz) ! z block of ztop of CH interval
      call getindx(nz,zmn,zsiz,rzbot,izbot_target,inflagz) ! z block of zbot of CH interval
      ithick_target = iztop_target - izbot_target + 1 !thickness equivalent to no. of z blocks

c
c Find the conditioned locations on this channel system and move in end locations
c

      no_prev_cond=.true.
      lowcond = 1
      highcond = ndis
      do idis = rdis-1, 1, -1      ! node from closest node to well intercept toward the first
       if(assoc_cond(cassoc,idis).eq.1) then
        lowcond = idis ! move in the low end location 
        no_prev_cond=.false.
       end if
      end do
      do idis = rdis+1, ndis ! node from the closest node to well intercept toward the end
       if(assoc_cond(cassoc,idis).eq.1) then
        highcond = idis ! move in the high end location
        no_prev_cond=.false.
       end if
      end do
      shiftmat = 0.0 ! initialize the shift matrix

c
c Set up the parabolic shift function - goes to 0 at the end points
c
      i=1
      do idis = 1, ndis

       if(idis.lt.rdis) then
        shiftmat(idis)=1.0-(real((idis-1)-(rdis-1))/real(rdis-1))**2.0
       else 
        shiftmat(idis)=1.0-(real(idis-rdis)/real(ndis-rdis))**2.0
       end if
      end do

      itol_lower = 0
      itol_upper = 0

      irepo = 50

      xCHi = chamat(rdis,1,rtime) ! the nearest location
      yCHi = chamat(rdis,2,rtime)
      init_xCHi = xCHi
      init_yCHi = yCHi
      call azimuth(xCHi,x,yCHi,y,azi) ! angle from current location to well

c
c Iterate to get the correct horizontal location, the channel is shifted towards the 
c the well until the correct thickness is found.  The convergence criteria
c is 1/4 cell size - or below the resolution of the model.
c
        write(*,*) "cstep  this_error  azi  xchi  ychi xwell ywell"   
c      step = (xsiz+ysiz)/2.0 ! first step is the average cell size
       cstep = (xsiz*real(nx)/20.0+ysiz*real(ny)/20.0)
      do while(.not.converge)
       continue
       call genchannellocal(rix,riy,rtime,closenode,iztop,izbot)
       ithick = iztop - izbot + 1 ! zero case is accounted for
       ithick_error=ithick-ithick_target
       if(ithick_error.le.itol_upper.and.
     +    ithick_error.ge.itol_lower) then
        write(*,*) 'Assoc,Centerline #',cassoc,rtime,' converge!'
        converge = .true.
       else
        xCHi = chamat(rdis,1,rtime) ! update the nearest location
        yCHi = chamat(rdis,2,rtime)
        call azimuth(xCHi,x,yCHi,y,azi) ! angle from current location to well
        if(ithick_error.gt.0.0) then
         cstep=-abs(cstep)*0.99 ! move away from well
        else
         cstep=abs(cstep)*0.99 ! move towards well
        end if
        call offset(x0,y0,azi,cstep,x1,y1) 
        cstepx = x1-x0
        cstepy = y1-y0
        xCHi=xCHi+cstepx
        yCHi=yCHi+cstepy
        do idis = 1, ndis
         chamat(idis,1,rtime) = chamat(idis,1,rtime)+cstepx* ! shift all nodes in streamline in x direction
     +      shiftmat(idis)
         chamat(idis,2,rtime) = chamat(idis,2,rtime)+cstepy* ! shift all nodes in streamline in Y direction
     +      shiftmat(idis)         
        end do
        call curvature2(rtime,extrapolate,regrid) 
        write(*,670) cstep,ithick_error,azi,xchi,ychi,x,y,iztop,izbot
      if(i.gt.100) then
      continue
      end if
      test3(1,i) = cstep
      test3(2,i) = xchi
      test3(3,i) = ychi
      test3(4,i) = x
      test3(5,i) = y
      test3(6,i) = azi
      i=i+1
        if(abs(cstep).lt.small_number) then
         converge = .true.
         write(*,*) 'Channel #',cassoc,' failed to converge!'
        end if
       end if ! end checking ithick_error
      end do ! end do while checking convergence
670   format(f10.3,1x,i9,5(f10.3,1x),2(i5,1x))    
      call curvature2(rtime,extrapolate,regrid)
c
c Now adjust the rest of the streamline association
c
 
      deltax = xCHi-init_xCHi ! how far the nearest node is shifted 
      deltay = yCHi-init_yCHi

      do i = 2, assoc_tab(cassoc,1)+1 ! now adjust the entire streamline assoc
       itime = assoc_tab(cassoc,i)
       if(itime.ne.rtime) then  
        ndis=chaprop(15,itime)
        do idis = 1, ndis
         chamat(idis,1,itime) = chamat(idis,1,itime)+deltax*
     +       shiftmat(idis)
         chamat(idis,2,itime) = chamat(idis,2,itime)+deltay*
     +       shiftmat(idis)         
        end do
        call curvature2(itime,extrapolate,regrid)
       end if 
      end do 

c
c Correct the vertical position for the streamline association
c
      
      deltaz = rztop-chamat(rdis,11,rtime) 

c
c Shift the vertical position - now shifted after horizontal set
c

      if(rtime.eq.3) then
       continue
      end if

      if(.not.no_prev_cond) then ! with previous conditioning
       do i = 2, assoc_tab(cassoc,1)+1 
        itime = assoc_tab(cassoc,i)
        ndis=chaprop(15,itime)
        do idis=1,ndis
         chamat(idis,11,itime) = chamat(idis,11,itime)+
     +    deltaz*shiftmat(idis)
        end do
        call curvature2(itime,extrapolate,regrid)
       end do
      else ! if no other conditioning just shift the channel 
       do i = 2, assoc_tab(cassoc,1)+1 
        itime = assoc_tab(cassoc,i)
        ndis=chaprop(15,itime)       
        do idis=1,ndis
        chamat(idis,11,itime) = chamat(idis,11,itime)+
     +   deltaz
        end do
        call curvature2(itime,extrapolate,regrid)       
       end do
      end if

      deallocate(shiftmat)

c
c Finished:
c
      return
      end



      subroutine movwinsmooth(ndim,nsiz,tmin,valmat)
c-----------------------------------------------------------------------
c
c                    Moving Window Smoothing Subroutine
c                    *******************************************
c
c Apply equal weighted moving window smoothing to a column in a matrix.
c
c INPUT VARIABLES:
c
c   ndim             the number of rows in the matrix
c   nsiz             the size of the moving window
c   valmat           the matrix
c   icol             the column to be smoothed 
c   ncol             the number of columns in the matrix
c
c OUTPUT VARIABLES:
c
c   valmat           original matrix with specificied column smoothed
c
c
c Author: Michael Pyrcz                                  DATE: 2001-2002 
c-----------------------------------------------------------------------

      use      streamsimMod

      integer   ix
      real      valmat(*),sumwt,wt,sum,tmin
      real,allocatable :: tempmat6(:)

c
      allocate(tempmat6(ndim),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if

c
c Apply the moving window smoothing
c
   
      do ix = 1, ndim
         sumwt = 0.0
         sum = 0.0
          do iix = -nsiz, nsiz
           if(ix+iix.ge.1.and.ix+iix.le.ndim) then
            if(valmat(ix+iix).gt.tmin) then
             wt = real((nsiz-abs(iix)+1)/real(nsiz+1))
             sum = sum + valmat(ix+iix)*wt
             sumwt = sumwt+wt
            end if
           end if
          end do
       if(sumwt.gt.0.0) then
        tempmat6(ix) = sum/sumwt
       else
        tempmat6(ix) = -999.9
       end if
      end do

c
c Copy the results over the original array
c

       do ix = 1, ndim
        valmat(ix) = tempmat6(ix)
       end do

      deallocate(tempmat6)
      

c
c Finished:
c
      return
      end



      subroutine movwinsmooth_azimuth(ndim,nsiz,tmin,valmat)
c-----------------------------------------------------------------------
c
c                    Moving Window Smoothing Subroutine for Azimuths
c                    ***********************************************
c
c Apply equal weighted moving window smoothing to a column in a matrix.
c
c INPUT VARIABLES:
c
c   ndim             the number of rows in the matrix
c   nsiz             the size of the moving window
c   valmat           the matrix
c   icol             the column to be smoothed 
c   ncol             the number of columns in the matrix
c
c OUTPUT VARIABLES:
c
c   valmat           original matrix with specificied column smoothed
c
c
c Author: Michael Pyrcz                                  DATE: 2001-2002 
c-----------------------------------------------------------------------

      use      streamsimMod

      integer   ix
      real      valmat(*),sumwt,wt,sum,tmin
      real      x_loc,y_loc,x1,y1,azi,azi0,azi1,azil,azih
      real,allocatable :: tempmat6(:)

c
      allocate(tempmat6(ndim),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if

c
c Apply the moving window smoothing
c
   

      do ix = 1, ndim
          x_loc = 0.0
          y_loc = 0.0
          do iix = -nsiz, nsiz
           if(ix+iix.ge.1.and.ix+iix.le.ndim) then
            if(valmat(ix+iix).gt.tmin) then
             wt = real((nsiz-abs(iix)+1)/real(nsiz+1))
             azi = valmat(ix+iix)
             call offset(x_loc,y_loc,azi,wt,x1,y1)
             x_loc = x1
             y_loc = y1
            end if
           end if
          end do
       if(x_loc.ne.0.0.or.y_loc.ne.0.0) then
         call azimuth(0.0,x_loc,0.0,y_loc,azi)
         tempmat6(ix) = azi 
       else
        tempmat6(ix) = -999.9
       end if
      end do

c
c Ensure continuity for spline operation
c

      do ix = 2, ndim
        azi0 = tempmat6(ix-1)
        azi1 = tempmat6(ix)
        azil = azi1-360.0
        azih = azi1+360.0 
        if(abs(azi0-azil).lt.abs(azi0-azi1)) then
          tempmat6(ix) = azil
        end if
        if(abs(azi0-azih).lt.abs(azi0-azi1)) then
          tempmat6(ix) = azih
        end if
      end do

c
c Copy the results over the original array
c

       do ix = 1, ndim
        valmat(ix) = tempmat6(ix)
       end do

      deallocate(tempmat6)
      

c
c Finished:
c
      return
      end




      subroutine neckcutoff(ctol,bcutoff,ctime1,ctime2)
    
c-----------------------------------------------------------------------
c
c     Check for Neck Cut Offs
c     ***********************
c
c Checks for a neck cutoff and removes oxbow lake segment
c
c INPUT VARIABLES:
c
c   ctol              the distance between control nodes for cut off
c   ctime1            the channel checked
c   ctime2            the channel with orbow lake removed
c
c OUTPUT VARIABLES:
c
c   bcutoff           true - cutoff occured, false - none
c   chamant(ndis,11,ctime2)  the modified channel
c
c
c
c EXTERNAL REFERENCES: azimuth    computes the azimuth
c                      spline     calculates spline
c                      splint     spline interpolation 
c
c Author: Michael Pyrcz                                  DATE: 2003-2004
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none

c
c Local Variables
c
 
      integer   j,dis_thresh,idis,jdis,idata,count
      real      thresh,xi,xj,yi,yj,cdist

c
c Looked Up Variables
c 
     
      integer   ndis
      real      maxCHhalfwidth

c
c Passed Variables
c

      integer   ctime1,ctime2
      real      ctol
      logical   bcutoff
             
c
c Begin
c

      write(*,*) ' Started: Cutoff event ',ctime1,' to ',ctime2
      if(ctime1.eq.20) then
        continue
      end if

      ndis = chaprop(15,ctime1)
      maxCHhalfwidth = chaprop(17,ctime1)

c
c Copy the streamline across and work on the copy
c

      do idis = 1, ndis0
       do idata = 1, ndata_chamat
         chamat(idis,idata,ctime2) = chamat(idis,idata,ctime1)
       end do
      end do

      do idis = 1, ndis0
       do idata = 1, ndata_spline_table
        spline_table(idis,idata,ctime2) = 
     +   spline_table(idis,idata,ctime1)
        spline_table2(idis,idata,ctime2) = 
     +   spline_table2(idis,idata,ctime1)
       end do
      end do

      do idata = 1, ndata_chaprop
       chaprop(idata,ctime2) = chaprop(idata,ctime1)
      end do

c
c Check for neck cut off and remove
c

      thresh = (ctol)**2 ! to save time work in square root
      dis_thresh = int((ctol)/
     + (chamat(ndis,3,ctime2)/ndis))+2 ! set as #nodes for tolerance + 2 node 
      bcutoff = .false.
435   do idis = 1, ndis
       do jdis = idis+dis_thresh, ndis
        xi = chamat(idis,1,ctime2)
        yi = chamat(idis,2,ctime2)
        xj = chamat(jdis,1,ctime2)
        yj = chamat(jdis,2,ctime2)
        cdist = (xi-xj)**2+(yi-yj)**2
        if(cdist.lt.thresh) then
         bcutoff = .true.
         count = 1
         continue
         do j = jdis+1,ndis
          chamat(idis+count,1,ctime2) = chamat(j,1,ctime2)
          chamat(idis+count,2,ctime2) = chamat(j,2,ctime2)
          count = count+1
         end do
         ndis = ndis - (jdis-idis)
         goto 435
        end if 
       end do
      end do

c
c Store reduced ndis 
c

      chaprop(15,ctime2)=ndis

c
c Finished:
c

      write(*,*) ' Finished: Cutoff event ',ctime1,' to ',ctime2
      return
      end


      subroutine offset(x1,y1,ang,dist,x2,y2)
    
c-----------------------------------------------------------------------
c
c     Calculate an Offset
c     *********************
c
c Given an initial point, angle (0.0-360.0 with 0 on y+), and a distance
c claculate the new point.
c
c
c INPUT VARIABLES:
c
c   x1            inital x   
c   y1            inital y 
c   ang           angle of offset   
c   dist          distance offset
c
c
c Author: Michael Pyrcz                                  DATE: 2002-2003
c-----------------------------------------------------------------------
      
      implicit none

c
c Local Variables
c
 
      real     PI

c
c Passed Variables
c

      real     x1,x2,y1,y2,dist,ang,rang,tang

c
c Begin
c

      PI=3.14159

c
c Test that ang is >0.0 and <360.0
c
      tang = ang

40    if(tang.lt.0.0) then
       tang = tang+360.0
       goto 40
      end if 

41    if(tang.gt.360.0) then
       tang = tang-360.0
       goto 41
      end if



c
c First quandrant
c

      if(tang.ge.0.0.and.tang.lt.90.0) then
       rang = tang*pi/180.0
       x2 = x1+dist*sin(rang)
       y2 = y1+dist*cos(rang)
       goto 566
      end if 

c
c Second quandrant
c
      if(tang.ge.90.0.and.tang.lt.180.0) then
       tang = tang-90.0
       rang = tang*pi/180.0
       x2 = x1+dist*cos(rang)
       y2 = y1-dist*sin(rang)
       goto 566
      end if

c
c Third quandrant 
c
      if(tang.ge.180.0.and.tang.lt.270.0) then
       tang = tang-180.0
       rang = tang*pi/180.0
       x2 = x1-dist*sin(rang)
       y2 = y1-dist*cos(rang)
       goto 566
      end if

c
c Fourth quandrant
c

      if(tang.ge.270.0.and.tang.lt.360.0) then
       tang = tang-270.0
       rang = tang*pi/180.0
       x2 = x1-dist*cos(rang)
       y2 = y1+dist*sin(rang)
       goto 566
      end if

566   continue

c
c Finished:
c
      return
      end



      subroutine ONEDRF(l,dim,valmat,tmean,tstdev)
c-----------------------------------------------------------------------
c
c      General a 1D Gaussian Random Field with Triangular Covariance
c      *************************************************************
c
c The 1D RF is generated by moving average based on a trangular weighting
c function.  Results are restandardize to N{tmean,tstdev}.
c
c
c INPUT VARIABLES:
c
c   l                the range of the trangular variogram in cells
c   dim              the number of cells in the generated array
c   valmat           the 1D RF
c
c OUTPUT VARIABLES:
c
c   valmat           the 1D RF
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
 
      implicit none

c
c Local Variables
c
 
      integer   i,j,ierr,test,idum
      real      xp,last,mean,var,stdev
      real,allocatable :: r1Dmat(:),w1Dmat(:)
      real*8    p,acorni

c
c Looked Up Variables
c 
     
c
c Passed Variables
c

      integer   l,dim
      real      tmean,tstdev
      real valmat(*)



      allocate(r1DMat(-l:dim+l),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if

      allocate(w1DMat(-l:l),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if

c
c Construct random N{0,1} vector and the weighting function
c
      
      do i = -l,dim+l
       p = real(acorni(idum))
       call gauinv(p,xp,ierr)
       r1dMat(i) = xp
      end do

      w1DMat(0) = 1.0/real(l)
      do i = 1,l
       w1Dmat(i) = w1DMat(0) - real(i)*(w1DMat(0)/(real(l)))
       w1DMat(-i) = w1DMat(i)
      end do

c
c Apply the weighting function as moving window
c             

      do i = 1,dim
       valmat(i) = 0.0
       do j = -l,l
        valmat(i) = valmat(i)+r1DMat(i+j)*w1DMat(j)
       end do
      end do

c
c Standardize to restore N{0,1}
c
      var = 0.0
      mean = 0.0

      do i = 1,dim
       mean = mean+valmat(i)
       var = var+valmat(i)**2
      end do
      
      mean = mean/real(dim)
      var = var/real(dim)-mean**2
      
      stdev = sqrt(max(var,0.001))  
      
      do i = 1,dim         
       valmat(i) = (valmat(i)-mean)*(tstdev/stdev)+tmean
      end do 
 
      deallocate(r1DMat)
      deallocate(w1DMat)      

c
c Finished:
c
      return
      end



      subroutine postprocessing 
    
c-----------------------------------------------------------------------
c
c     Postprocessing 
c     **************
c
c This subroutine performs postprocessing to all net intervals by moving 
c each intersected channel or attracting each closest channels 
c horizontally to match the thickness of each interval.  It is then 
c moved vertically to match the top of that interval.
c
c Author: Michael Pyrcz                                  DATE: 2003-2004 
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none

c
c Local Variables
c

      integer iwell
      integer i,iCHdata,iwelldata
      integer ixWell,iyWell,wellnum 
      integer closetime,closeassoc,closenode 
      real      ztop,zbot,x,y     
      logical iNodeFound 
      
      
      do iCHdata = 1, nCHdata       
          wellnum = welldata(iCHdata,1)
          ixWell = welldata(iCHdata,2) 
          iyWell = welldata(iCHdata,3)

          x = real(ixWell-1)*xsiz+xmn
          y = real(iyWell-1)*ysiz+ymn
          ztop = welldata(iCHdata,4)
          zbot = welldata(iCHdata,5)

c
c Find the nearest spline node
c
          call locatenode(x,y,ztop,closetime,    
     +               closeassoc,closenode,iNodeFound)          
c
          if(.not.iNodeFound) then
              write(*,*) 'Failed to find suitable centerline for'
              write(*,*) 'Data = ',iCHdata
              write(*,*) 'Likely caused by model set to stop early'
              go to 13
          end if
          
c
c Update welldata and wellheader 
c
          welldata(iCHdata,6) = closeassoc
          welldata(iCHdata,7) = closetime
          do iwell = 1, nwell
              if(wellnum.eq.wellheader(iwell,1)) then ! find the correct wellnum
                    wellheader(iwell,7) = wellheader(iwell,7)+1
                    wellheader(iwell,(wellheader(iwell,7)+7)) = closeassoc
                    wellheader2(iwell,(wellheader(iwell,7))) = closetime
                  go to 10
              endif
          enddo
c
c Paint the buffer and close node in assoc_cond
c
  10          do i=-node_buffer,node_buffer 
                if(closenode+i.ge.1.and.closenode+i.le.ndis0) then
                      assoc_cond(closeassoc,closenode+i) = 2
                endif
          enddo
c
c Move streamline horizontally and vertically
c                              
          call morphCHendpts(closetime,closenode,
     +                       ixWell,iyWell,ztop,zbot)

c
c Update the conditioning location
c
          assoc_cond(closeassoc,closenode) = 1
          assoc_cond(closeassoc,ndis0+1) 
     +                    = assoc_cond(closeassoc,ndis0+1)+1 

      enddo ! all intervals
 13   continue

c Finish
        write(*,*)'Finish postprocessing'
      write(*,*)'*********************'
      return
      end      


      subroutine probdraw(valmat,nrow,irow,val)
c-----------------------------------------------------------------------
c
c                    Select from a weighted distibution
c                    **********************************
c
c
c INPUT VARIABLES:
c
c   valmat           the 1D array with the weights
c   nrow             the number of values in the array
c
c
c OUTPUT VARIABLES:
c
c   irow             the randomly selected row
c   val              the value in the associated row
c
c
c Author: Michael Pyrcz                                  DATE: 2002-2003
c-----------------------------------------------------------------------
 
      implicit none

c
c Local Variables
c
 
      integer idum
      real    last
      real*8  p,acorni

c
c Passed Variables
c

      integer irow,nrow
      real    val,valmat(nrow),tempmat5(nrow)

c
c Convert array to cummulative
c
    
      last = 0.0
      do irow = 1, nrow
       tempmat5(irow) = valmat(irow) + last
       last = tempmat5(irow)
      end do

      p = real(acorni(idum)) ! check that this works at run time with local idum
      p = p*tempmat5(nrow)
      
      do irow = 1, nrow
       if(p.lt.tempmat5(irow)) then
        goto 7031
       end if
      end do
      irow = nrow

7031      continue
       
      val = valmat(irow)

c
c Finished:
c
      return
      end




      subroutine segmentwell
    
c-----------------------------------------------------------------------
c
c     Segment Well Intercepts into Multiple Arrays 
c     ********************************************
c
c Given an initial point, angle (0.0-360.0 with 0 on y+), and a distance
c claculate the new point.
c
c
c INPUT VARIABLES:
c
c   x1            inital x   
c   y1            inital y 
c   ang           angle of offset   
c   dist          distance offset
c
c Author: Michael Pyrcz                                  DATE: 2003-2004     
c-----------------------------------------------------------------------
      
      use streamsimMOD
      implicit none

      logical  inflagx,inflagy,inflagzt,inflagzb
      integer  i,idata,jdata,ix,iy,iztop,izbot
      integer  iCHdata,izt,izb,iwell,wellnum,idum,nWellCHdata(10)
      real     x1,x2,y1,y2,dist,ang,rang,tang
      real     d1,d2,d3,d4,d5,d6 ! dummy variables
      real*8   acorni
      integer,allocatable :: temp(:)

c
c Segment channels
c      

c     nWellCH = 0
c      nWellLA = 0
c      nWellLV = 0
c      nWellCS = 0
      
      nCHdata = 0
      nWell = 0

c      do iwell = 1, nwell
c       if(well(iwell,6).eq.kCH) then
c        nWellCH = nWellCH + 1
c       else if(well(iwell,6).eq.kLA) then
c        nWellLA = nWellLA + 1
c       else if(well(iwell,6).eq.kLV) then
c        nWellLV = nWellLV + 1
c       else if(well(iwell,6).eq.kCS) then
c        nWellCS = nWellCS + 1
c       end if 
c      end do

      write(*,*) '_______________________'
      write(*,*) 'In Subroutine segmentwell'
c
c Count the number of CH data (nCHdata) and number of wells (nWell)
c
     
      iCHdata = 0
      do idata = 1, ndata
       if(well(idata,4).le.zmax.and.well(idata,5).ge.zmin) then
        if(well(idata,6).eq.kCH.or.well(idata,6).eq.kLA) then
         call getindx(nx,xmn,xsiz,well(idata,2),ix,inflagx)
         call getindx(ny,ymn,ysiz,well(idata,3),iy,inflagy)
         if(inflagx.and.inflagy) then         
          nCHdata = nCHdata + 1
           do jdata = idata-1,1,-1
           if(well(idata,1).eq.well(jdata,1)) then 
            nWellCHdata(nWell) = nWellCHdata(nWell) + 1 !add no. of CH' interval for each well
              goto 43
           end if
          end do
          nWell = nWell + 1
            nWellCHdata(nWell) = 1
43        continue
            write(*,*) 'nWell= ',nWell,'nCHdata= ',nCHdata
            write(*,*) 'nWellCHdata(nWell)= ',nWellCHdata(nWell)
         end if
        end if         
       end if
      end do      
      
c
c Check if there is any conditioning to apply
c       

      if(nCHdata.le.0) then
       Stop 'No Conditioning In Model Space - No Updating Performed'
      end if

c
c Conditioning Arrays
c

c
      allocate(well_log(nz),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c      
      allocate(dataorder(nCHdata),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if      
c
      allocate(temp(nCHdata),stat = test)
      if(test.ne.0)then
       write(*,*)'ERROR 5: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
c
      allocate(wellheader(nWell,ndata_wellheader),stat = test)
      if(test.ne.0)then
       write(*,*)'Well Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      wellheader = 0.0
c
      allocate(wellheader2(nWell,ndata_wellheader),stat = test)
      if(test.ne.0)then
       write(*,*)'Well Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      wellheader2 = 0.0
c
      allocate(welldata(nCHdata,9),stat = test)
      if(test.ne.0)then
       write(*,*)'Well Array: Allocation failed',
     + ' due to insufficient memory.'
       stop
      end if
      welldata = 0.0
c

      iCHdata = 0
      iWell = 0
      do idata = 1, ndata
       if(well(idata,4).lt.zmax.and.well(idata,5).gt.zmin) then
        if(well(idata,6).eq.kCH.or.well(idata,6).eq.kLA) then
         call getindx(nx,xmn,xsiz,well(idata,2),ix,inflagx)
         call getindx(ny,ymn,ysiz,well(idata,3),iy,inflagy)
         call getindx(nz,zmn,zsiz,well(idata,4),izt,inflagzt)
         call getindx(nz,zmn,zsiz,well(idata,5),izb,inflagzb)
         if(inflagx.and.inflagy) then         
          iCHdata = iCHdata + 1
          welldata(iCHdata,1) = well(idata,1)
          welldata(iCHdata,2) = ix
          welldata(iCHdata,3) = iy
          welldata(iCHdata,4) = well(idata,4)
            welldata(iCHdata,5) = well(idata,5)
            welldata(iCHdata,8) = izt      !FZ 03/24/05
            welldata(iCHdata,9) = izb   !FZ 03/24/05
           do jdata = idata-1,1,-1
           if(well(idata,1).eq.well(jdata,1)) then
            goto 44
           end if
          end do
          iWell = iWell + 1
          do i = 1, 5
           wellheader(iWell,i) = welldata(iCHdata,i)
          end do      
44        continue
         end if
        end if         
       end if
      end do       

c
c Set the data order
c
      
      do idata = 1, nCHdata
       dataorder(idata) = idata
       temp(idata) = real(acorni(idum))
      end do
      call sortem(1,nCHdata,temp,1,dataorder,d1,d2,d3,d4,d5,d6) 
      deallocate(temp)

c
c Fill out the well header, iztop, izbot
c

      do iwell = 1, nWell
        wellnum = wellheader(iwell,1)
        do idata = 1, nData
          if(int(well(idata,1)).eq.int(wellnum)) then
            wellheader(iwell,4) = max(wellheader(iwell,4),
     +       well(idata,4))
            wellheader(iwell,5) = min(wellheader(iwell,5),
     +       well(idata,5))
                   wellheader(iwell,6) = nWellCHdata(iwell) ! FZ Assign no. of CH' intervals
              write(*,*) 'iwell= ',iwell  
              write(*,*) 'wellheader(iwell,6)= ',wellheader(iwell,6) 
          end if 
        end do
      end do

c
c Convert ztop and zbot to iztop and izbot in well header
c
 
      do iwell = 1, nWell
       call getindx(nz,zmn,zsiz,wellheader(iWell,4),iztop,inflagzt)
       wellheader(iWell,4)=iztop
       call getindx(nz,zmn,zsiz,wellheader(iWell,5),izbot,inflagzb)
       wellheader(iWell,5)=izbot
      end do      

c
c Finished:
c
      return
      end



      subroutine selectstreamline(ilevel,nCHtoMatch,iwell,iCHtoMatch,
     +                            rCase)
    
c-----------------------------------------------------------------------
c
c     Accept or Reject Generated Channel Objects 
c     ******************************************
c
c The subroutine checks if thickness of channel objects intersect with 
c the well within the thickness tolerance
c 
c
c INPUT VARIABLES:
c     ilevel = aggradation level
c     nCHtoMatch = number of intervals to match on ilevel
c     iwell = well number (1 to nwell)
c     iCHtoMatch = interval number to match on ilevel (1 to nCHtoMatch 
c      on ilevel)
c
c OUTPUT VARIABLES:
c     rCase = case number of the selected streamline
c
c GLOBAL VARIABLE:
c
c   chaprop      channel property array
c   welldata 
c   level
c   CHcurrent current channel number 
c   nwell     total number of wells        
c
c
c Author: Fueng Zabel                                    DATE: 2004-2005
c-----------------------------------------------------------------------
      
      use       streamsimMod
      implicit none

c
c Local Variables
c

      integer maxCHdrawn,CHlevel,iwell(10),ilevel,iCHtoMatch(10),iMatch
      integer i,prevCase,rCase,nCHtoMatch,keepCHdraw
      real    zbot,t1,t2,val
      logical iIntercept,iInterceptFound,iThickness,iNoUnwarrant
      
      maxCHdrawn = 1000 ! set up maximum no. of times allowing to draw
      CHlevel=1
      prevCase=0
      rCase=0
      keepCHdraw=0

     
   1      if (CHlevel.le.maxCHdrawn) then

             call probdraw(chadrawwt,nCHdraw,CHdraw,val)
c
c Lookup and copy streamline (with properties) from streamline table 
c into active arrays (chaprop and chamat)
c                  
          call lookupstream(ilevel,CHdraw,CHcurrent)
c
c Applying ACCEPTING/REJECTING RULES
c
c
c 1. Check for unwarranted intercept
c    -------------------------------
c
          iNounwarrant = .false.
          do i=1,nwell
              call unwarrantedintercept(i,iNounwarrant)
              if(.not.iNoUnwarrant) then
                  rCase = 4 ! with unwarranted intercept
                  CHlevel = CHlevel+1
                  write(*,*) 'rCase= ',rCase  
                    write(*,*) 'Case4. Found unwarranted intercept...',
     +                            'draw a new one'
                  go to 1
              endif  
          enddo ! all wells
c
c 2. Check for well intercept (if the drawn channel intersects with any 
c    net interval on ilevel)
c    ------------------------------------------------------------------
c
          if(nCHtoMatch.gt.0) then
              iIntercept = .false.
              iInterceptFound = .false.
              iThickness = .false.
              do iMatch = 1,nCHtoMatch 
                  call wellintercept(iwell(iMatch),iIntercept)
                  if(iIntercept) then
                      iInterceptFound = .true.
c
c Check if thickness of net matching > thickness of unwarranted net
c
                           zbot = welldata(iCHtoMatch(iMatch),5)
                      t1 = abs(level(ilevel,1) - zbot)
                        t2 = abs(zbot - chaprop(2,CHcurrent))

                      if(t1.gt.t2) then
                          iThickness = .true. 
                          go to 2
                      endif
                  endif ! end well intercept                           
              end do ! end all intervals to match        
          else ! no interval to match on ilevel
              go to 30   
          endif ! end interval to match 

   2      if(iInterceptFound.and.iThickness) then
              rCase = 1 ! pass thickness criteria with no unwarranted intercept
              write(*,*) 'rCase= ',rCase  
              write(*,*) 'Case1. Pass thickness criteria with ',
     +                                 'no unwarranted intercept'
              go to 30
          else if (iInterceptFound) then
              rCase = 2 ! intercept with no unwarranted intercept
              write(*,*) 'rCase= ',rCase  
              write(*,*) 'Case2. Intercept with no unwarranted ',
     +                                     'intercept...draw a new one'
          else
              rCase = 3 ! no intercept with no unwarranted intercept
              write(*,*) 'rCase= ',rCase  
              write(*,*) 'Case 3. No intercept with no unwarranted ',
     +                              'intercept...draw a new one'
          endif ! end assigning different case
          
c
c Keep the good one before drawing a new streamline if it does not reach 
c maximum times allowed to draw
c
          if (prevCase.eq.0.or.
     +        (prevCase.gt.0.and.rCase.le.prevCase)) then
              keepCHdraw = CHdraw
          endif
          prevCase = rCase
          CHlevel = CHlevel+1

      endif ! end maximum times allowed to draw checking

c
c 3. Take the best case if rCase=1 is not found
c    ----------------------------------------------
c
       CHdraw = keepCHdraw
       keepCHdraw = 0
       call lookupstream(ilevel,CHdraw,CHcurrent)
       
c Finish
  30      write(*,*)'Finish selectstreamline'
      write(*,*)'****************************************************'
      return
      end      


      subroutine unwarrantedintercept(iwell,iNotFound)
    
c-----------------------------------------------------------------------
c
c     Check unwarranted intercept
c     ***************************
c
c Purpose: To check for unwarranted intercept (to see if the channel 
c object intersects with non-net interval)
c         
c
c INPUT VARIABLES:
c     iwell = well number
c
c OUTPUT VARIABLES:
c      iNotFound = true if there is no unwarranted intercept
c                 false if there is unwarranted intercept
c
c Author: Fueng Zabel                                    DATE: 2004-2005
c-----------------------------------------------------------------------
      
      use streamsimMOD
      implicit none

      logical  iNotFound
      integer  iwell,iCHdata,ilevel,idis,ndis,nCHinterval
      integer  wellnum,wellix,welliy,welliz,iz,ixNode,iyNode,izNode
      real     z,topz,botz,xNode,yNode,zNode 

c
c Loop over vertical wells
c
   
      iNotFound = .true.
      wellnum = wellheader(iwell,1)
      wellix = wellheader(iwell,2) ! block no. in X direction
      welliy = wellheader(iwell,3)
      nCHinterval = wellheader(iwell,6) ! no. of well CH intervals

c Calculate well log
c

      well_log = 0
      do iCHdata = 1, nCHinterval
        if(int(welldata(iCHdata,1)).eq.wellnum) then
         topz = welldata(iCHdata,4)
         botz = welldata(iCHdata,5)
         do iz = 1, nz
          z = real(iz-1)*zsiz+zmn ! Z location for the bottom of iz block
          if(z.gt.botz.and.z.lt.topz) then
             welliz = iz 
           well_log(iz) = 1
          end if
         end do ! loop over all z blocks
        end if
      end do !loop over all well CH' interval

c
c Loop through all the nodes of a streamline and check if any node violate
c non-net intercept
  
       ndis = chaprop(15,CHcurrent)
       do idis=1,ndis
            xNode = chamat(idis,1,CHcurrent)
            yNode = chamat(idis,2,CHcurrent)
            zNode = chamat(idis,11,CHcurrent)
            ixNode = int((xNode-xmn)/xsiz) + 1
            iyNode = int((yNode-ymn)/ysiz) + 1
            izNode = int((zNode-zmn)/zsiz) + 1
            if((ixNode.eq.wellix).and.(iyNode.eq.welliy).and.
     +       (well_log(izNode).eq.0)) then
                  iNotFound = .false. 
                  go to 10
            endif
       end do ! loop over all nodes in streamline
c
c Finished:
c
  10  write(*,*) 'Finish unwarrantedintercept'
      write(*,*) '***********************************************'
      return
      end


      subroutine wellintercept(iwell,iFound)
    
c-----------------------------------------------------------------------
c
c     Check well intersection
c     ************************
c
c Purpose: To check if the channel object intersects with the well net 
c interval for the same aggradation level 
c
c
c INPUT VARIABLES:
c     iwell = well number
c
c OUTPUT VARIABLES:
c      iFound = true if CH object is found to be in the same block as 
c      well CH' interval false if not found
c
c Author: Michael Pyrcz                                  DATE: 2003-2004    
c-----------------------------------------------------------------------
      
      use streamsimMOD
      implicit none

      logical  iFound
      integer  iwell,iCHdata,ilevel,idis,ndis,nCHinterval
      integer  wellnum,wellix,welliy,welliz,iz,ixNode,iyNode,izNode
      real     z,topz,botz,xNode,yNode,zNode 

c
c Loop over vertical wells
c
   
      iFound = .false.
      wellnum = wellheader(iwell,1)
      wellix = wellheader(iwell,2) ! block no. in X direction
      welliy = wellheader(iwell,3)
      nCHinterval = wellheader(iwell,6) 
c
c Calculate well log
c

      well_log = 0
      do iCHdata = 1, nCHinterval
        if(int(welldata(iCHdata,1)).eq.wellnum) then
         topz = welldata(iCHdata,4)
         botz = welldata(iCHdata,5)
         do iz = 1, nz
          z = real(iz-1)*zsiz+zmn ! Z location for the bottom of iz block
          if(z.gt.botz.and.z.lt.topz) then
             welliz = iz 
           well_log(iz) = 1
          end if
         end do ! loop over all z blocks
        end if
      end do !loop over all well CH' interval

c
c Loop through all the nodes of a streamline and check if any node is in the same 
c block as well net interval
  
       ndis = chaprop(15,CHcurrent)
       do idis=1,ndis
            xNode = chamat(idis,1,CHcurrent)
            yNode = chamat(idis,2,CHcurrent)
            zNode = chamat(idis,11,CHcurrent)
            ixNode = int((xNode-xmn)/xsiz) + 1
            iyNode = int((yNode-ymn)/ysiz) + 1
            izNode = int((zNode-zmn)/zsiz) + 1
            if((ixNode.eq.wellix).and.(iyNode.eq.welliy).and.
     +       (well_log(izNode).eq.1)) then
                  iFound = .true. 
                  go to 10
            endif
       end do ! loop over all nodes in streamline
c
c Finished:
c
  10  write(*,*) 'Finish wellintercept'
      write(*,*) '********************'
      return
      end


c-----------------------------------------------------------------------
c
c     General Utility Subroutines 
c-----------------------------------------------------------------------



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
c From Deutsch and Journel, 1999, GSLIB: Geostatistical Software Library
c   and User's Guide. Oxford University Press, New Yourk, 2nd Edition. 
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

      subroutine gauinv(p,xp,ierr)
c-----------------------------------------------------------------------
c
c Computes the inverse of the standard normal cumulative distribution
c function with a numerical approximation from : Statistical Computing,
c by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95.
c
c
c
c INPUT/OUTPUT:
c
c   p    = double precision cumulative probability value: dble(psingle)
c   xp   = G^-1 (p) in single precision
c   ierr = 1 - then error situation (p out of range), 0 - OK
c
c
c-----------------------------------------------------------------------
      real*8 p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim,p
      save   p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim
c
c Coefficients of approximation:
c
      data lim/1.0e-10/
      data p0/-0.322232431088/,p1/-1.0/,p2/-0.342242088547/,
     +     p3/-0.0204231210245/,p4/-0.0000453642210148/
      data q0/0.0993484626060/,q1/0.588581570495/,q2/0.531103462366/,
     +     q3/0.103537752850/,q4/0.0038560700634/
c
c Check for an error situation:
c
      ierr = 1
      if(p.lt.lim) then
            xp = -1.0e10
            return
      end if
      if(p.gt.(1.0-lim)) then
            xp =  1.0e10
            return
      end if
      ierr = 0      
c
c Get k for an error situation:
c
      pp   = p
      if(p.gt.0.5) pp = 1 - pp
      xp   = 0.0
      if(p.eq.0.5) return
c
c Approximate the function:
c
      y  = dsqrt(dlog(1.0/(pp*pp)))
      xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) /
     +               ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
      if(real(p).eq.real(pp)) xp = -xp
c
c Return with G^-1(p):
c
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
c From Deutsch and Journel, 1999, GSLIB: Geostatistical Software Library
c   and User's Guide. Oxford University Press, New Yourk, 2nd Edition. 
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



      real function resc(xmin1,xmax1,xmin2,xmax2,x111)
      real*8 rsc
c
c Simple linear rescaling (get a value in coordinate system "2" given
c a value in "1"):
c
      rsc  = dble((xmax2-xmin2)/(xmax1-xmin1))
      resc = xmin2 + real( dble(x111 - xmin1) * rsc )
      return
      end







c-----------------------------------------------------------------------
c
c     To Be Added by User
c-----------------------------------------------------------------------
c
c   spline.for - Press,W.H., Flannery, B. P., Teukolsky, S.A.,& 
c     Vetterling,W. T. 1992,Numerical Recipes in Fortran 90: The Art of 
c     Scientific Computing (Cambridge: Cambridge Univ. Press)
c
c   Used for spline interpolation between cetenterline control nodes.
c
c   splint.for - Press,W.H., Flannery, B. P., Teukolsky, S.A.,& 
c     Vetterling,W. T. 1992,Numerical Recipes in Fortran 90: The Art of 
c     Scientific Computing (Cambridge: Cambridge Univ. Press)
c
c   Used to set up spline interpolation between cetenterline control 
c     nodes.
