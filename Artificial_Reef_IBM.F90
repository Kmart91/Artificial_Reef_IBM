!     Last change: KM 2-10-22

MODULE allvar

  IMPLICIT NONE
  INTEGER, PARAMETER:: imodf=1000      !the number of individuals per age-class
  INTEGER,PARAMETER:: inyears = 10
  INTEGER, PARAMETER ::  nspecies = 6    !*NEW - inserted variable nspecies so that itotf is dimensioned for the exact number of species/individuals MDC - 1apr09
  INTEGER, PARAMETER :: nages = 10       !*NEW - Number of age groups (allocates appropriate space for vul and ksat)
  INTEGER, PARAMETER:: itotf=imodf*inyears*nspecies    !total model individuals  !*FIX - inserted nspecies
  INTEGER idum                        ! random number seed
  INTEGER, PARAMETER:: nrow=100       ! number of rows
  INTEGER, PARAMETER:: ncol=100       ! number of columns
  INTEGER, PARAMETER:: nhab=4         ! number of different habitat types
  INTEGER, PARAMETER:: npyr=3         ! number of pyramid cells
  INTEGER, PARAMETER:: numnatcells = 128  ! number of natural bank cells to represent Baker bank
  REAL*8  xweight(itotf),xgrowth(itotf),xworth(itotf)  ! weight in g ww Not-used worth
  INTEGER xhab(itotf), xpyr(itotf), xnat(itotf)        !  0 if not and 1 if on rig during daytime movement
  INTEGER xxcell(itotf),xycell(itotf) !  x cell number and y cell number - lower left corner is 1,1
  INTEGER xage(itotf)                 ! age in years 1=1 to 12 months
  INTEGER xalive(itotf)               ! 0 alive 1 means worth too small
  INTEGER xspecies(itotf)             !  1=red snapper 2=grunts 3=grouper 4=bluefish 5=jacks 6=ga
  REAL*8 xdailymove(itotf)           ! total daily movement (m)
  REAL*8 xcmort(itotf)             ! daily cumulative mortality experienced
  REAL*8 grow                         ! g ww pred/g wwpred/hour
  INTEGER ifish, inight, ihour        ! id 1=night and 2=day  hours from 1 to 24
  INTEGER iday, jday, iyear           ! model day  julian day  and year counter
  REAL*8 light_hours                  ! number of sunlight hours in a day, eqn9001 units are initially in minutes
  INTEGER nrigs, rigmeth              !  number of cells that have rigs
  INTEGER pyrdens(npyr)               ! density of pyramids. will adjust this value depending on the number of pyramids in each grid in the simulation
  REAL*8 pyrweight(npyr), rigweight, natweight  ! weightings for each habitat type
  REAL*8  cellsize                    ! size of cells in meters
  INTEGER zhab(ncol,nrow)             ! 0 = rig cell.  1 = benthic cell. 2 = pyramid cell. 3 - natural bank - used to multiply by prey densities
  INTEGER habtype                     ! 1-benthic,0=rig, 2=pyr, 3=nat
  INTEGER habloc(nhab,ncol,nrow)      ! creates a variable to assign zhab specific to the habitat type at a certain location
  INTEGER rigcol(ncol*nrow),rigrow(ncol*nrow)  ! column and row number of the ith rig
  integer pyrcol(ncol*nrow),pyrrow(ncol*nrow)  ! column number of ith pyramid and row number of ith pyramid
  integer natcol(ncol*nrow),natrow(ncol*nrow)  ! column number of ith nat cell and row number of ith nat cell
  integer colnearrig(ncol,nrow), rownearrig(ncol,nrow) ! column number of nearest rig cell and row number
  integer colnearpyr(ncol,nrow), rownearpyr(ncol,nrow) ! column and row number of nearest pyramid cell
  integer colnearnat(ncol,nrow), rownearnat(ncol,nrow) ! column and row number of nearest pyramid cell
  REAL*8 distnearrig(ncol,nrow), distnearpyr(ncol,nrow), distnearnat(ncol,nrow)   ! distance in meters to center of nearest rig cell
  REAL*8 cumhours, cumdays
  INTEGER, PARAMETER:: nprey=5       ! number of prey types
  REAL*8 kprey(ncol,nrow,nprey) ! carrying capacity (g ww/m2) in logistic for each prey type by cell
  REAL*8 rprey(ncol,nrow,nprey) ! population growth (1/hour) in logistic for each prey type by cell
  REAL*8 zprey(ncol,nrow,nprey)  ! actual prey density (g ww/m2)
  REAL*8 kkprey(nprey),rrprey(nprey) ! temporary storage of kprey and rprey
  REAL*8 xeat(itotf,nprey),geat(nprey) ! g prey/g predator/day per hour??
  REAL*8 cumeat(ncol,nrow,nprey)  ! sum of prey type nprey eaten in each cell, summed over all fish inthe that cell g prey/day -SCALED
  REAL*8 watertemp      ! water tmperature in degrees C
  REAL*8 watertempbycol(ncol,nrow)
  REAL*8 wtemp(itotf)
  REAL*8 depth(ncol)
  INTEGER xstage(itotf) ! aggregated ages 1=ages 1 and 2, 2=ages 3 and 4
  REAL*8 totcon(itotf)  ! g of all prey eaten/g snapper/hour, sum of the con(jj) by prey type over an hour (reset every hour to 0)
  REAL*8 dailyp(itotf)  ! ratio of actual daily consumption (sum of con(jj)) / gcmax
  REAL*8 vul(nspecies,nprey,10), ksat(nspecies,nprey,10) ! parameter of functional response - species, prey type, life stage
  REAL*8 con(nprey)  ! g prey type eaten/g snapper/day in a hour?
  INTEGER idaval(itotf)  ! array of available ids for new fish
  REAL*8 totpop(nspecies)  ! initial of age-1 individuals per rig by species
  REAL*8 nmort, nmorth,pmort,pmorth
  REAL*8 fmort, fmorth
  REAL*8 initnum(nspecies,inyears)
  REAL*8 daydist(inyears),nightdist(inyears)
  REAL*8 initwt(nspecies,inyears)  ! initial weight by age in g wet weight
  REAL*8 dailyc(itotf)   ! g all of prey eaten/g snapper/day - sum of totcon over hours of geeding in a day
  REAL*8 xxloc(itotf),xyloc(itotf)
  REAL*8 xdistnearrig(itotf)  ! distance in meters from individuals's location to the cetner of the nearest rig cell
  REAL*8 xdistnearpyr(itotf)
  REAL*8 xdistnearnat(itotf)
  INTEGER xcolnearrig(itotf),xrownearrig(itotf)  ! column and row of nearest rig cell
  INTEGER xcolnearpyr(itotf), xrownearpyr(itotf)
  INTEGER xcolnearnat(itotf), xrownearnat(itotf)
  REAL*8 yield(nspecies),yieldnum(nspecies),yieldwt(nspecies) ! summed hourly yield for each day - biomass, numbers, average weight per ind, by species
  REAL*8 naturalp(nspecies), naturalpnum(nspecies) ! summed hourly natural mortality losses by day - biomass and numbers
  REAL*8 growthp(nspecies)  ! summed hourly growth production for each day by species - g ww predator per day
  REAL*8 ener(nprey)  ! calories per g ww by prey type
  REAL*8 rec(nspecies)  ! initial number of new recruits (age-1) each year by species
  REAL*8 en(nspecies)   ! energy density of the three fish species
  INTEGER ipinfish,icroaker,ijack,ibluefish,sumfish
  REAL*8 gcmax,avgden,egest,exc,resp,egg
  REAL*8 predavg   ! instantenous mortality rate (1/hour)
  INTEGER xriskhr(itotf)  ! number of hours each individual spent off the rig during daylight for each day
  INTEGER xuniqrig(itotf), xuniqpyr(itotf), xuniqnat(itotf) ! number of times rig cell changes from previous daytime rig cell
  INTEGER xxlastrig(itotf), xxlastpyr(itotf), xxlastnat(itotf) ! column number of last daytime rig cell visited
  INTEGER xylastrig(itotf), xylastpyr(itotf), xylastnat(itotf) ! row number of last daytime rig cell visited
  REAL*8  xdistcentrig(itotf), xdistcentpyr(itotf), xdistcentnat(itotf)     ! distance of fish from the previous daytime rig cell
  INTEGER i,j,k
  INTEGER meth         ! Method used to place rigs on grid (read from input file)
  REAL*8 :: ga(nspecies),gb(nspecies),gtmax(nspecies),gtopt(nspecies),gtheta1(nspecies)  ! consumption bioenergetics terms
  REAL*8 :: ra(nspecies),rb(nspecies),rtmax(nspecies),rtopt(nspecies),rtheta2(nspecies),act(nspecies)  ! respiration bioenergetics terms
  REAL*8 :: fa(nspecies),fb(nspecies),fg(nspecies)  !egestion bioenergetics terms
  REAL*8 :: ua(nspecies),ub(nspecies),ug(nspecies)  !excretion bioenergetics terms
  REAL*8 :: movebio(ncol,nrow),bluebio,jackbio, gabio, bio(nspecies,ncol,nrow)
  REAL*8 :: onto_mov
  INTEGER, PARAMETER :: iout=1    !when iout = 1 reduces the amount of written files when running production size model
  INTEGER, PARAMETER :: newrigs=0  !when newrigs > 0 choose to add (newrigs = x) number of platforms at some set interval or year determined at the first of the year.
  INTEGER, PARAMETER :: remrigs=0  !when remrigs > 0 choose to remove (remrigs = x) number of platforms at some set interval or year determined at the first of the year.
  !REAL*8 eqn9001
END MODULE allvar
!-------------------------------------------------------------------------------------------------------
!***** new numerical recipes stuff below here
! you can ignore these three modules as they are only needed for numerical recipes routines
         MODULE nrtype
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
	INTEGER, PARAMETER :: SP = KIND(1.0)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: LGT = KIND(.true.)
	REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
	REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
	REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
	REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
	REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
	REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
	REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
	REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
	TYPE sprs2_sp
		INTEGER(I4B) :: n,len
		REAL(SP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_sp
	TYPE sprs2_dp
		INTEGER(I4B) :: n,len
		REAL(DP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_dp
END MODULE nrtype
!********************************************************************************
MODULE nrutil
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
	INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
	INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
	INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
	INTEGER(I4B), PARAMETER :: NPAR_POLY=8
	INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
	INTERFACE
		FUNCTION rank(indx)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		INTEGER(I4B), DIMENSION(size(indx)) :: rank
		END FUNCTION rank
	END INTERFACE
	INTERFACE swap
		MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
			swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
			masked_swap_rs,masked_swap_rv,masked_swap_rm
	END INTERFACE
	INTERFACE assert_eq
		MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
	END INTERFACE
	INTERFACE arth
		MODULE PROCEDURE arth_r, arth_d, arth_i
	END INTERFACE
CONTAINS
	SUBROUTINE swap_i(a,b)
	INTEGER(I4B), INTENT(INOUT) :: a,b
	INTEGER(I4B) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_i
!BL
	SUBROUTINE swap_r(a,b)
	REAL(SP), INTENT(INOUT) :: a,b
	REAL(SP) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_r
!BL
	SUBROUTINE swap_rv(a,b)
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
	REAL(SP), DIMENSION(SIZE(a)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_rv
!BL
	SUBROUTINE swap_c(a,b)
	COMPLEX(SPC), INTENT(INOUT) :: a,b
	COMPLEX(SPC) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_c
!BL
	SUBROUTINE swap_cv(a,b)
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
	COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_cv
!BL
	SUBROUTINE swap_cm(a,b)
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
	COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_cm
!BL
	SUBROUTINE swap_z(a,b)
	COMPLEX(DPC), INTENT(INOUT) :: a,b
	COMPLEX(DPC) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_z
!BL
	SUBROUTINE swap_zv(a,b)
	COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
	COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_zv
!BL
	SUBROUTINE swap_zm(a,b)
	COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
	COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_zm
!BL
	SUBROUTINE masked_swap_rs(a,b,mask)
	REAL(SP), INTENT(INOUT) :: a,b
	LOGICAL(LGT), INTENT(IN) :: mask
	REAL(SP) :: swp
	if (mask) then
		swp=a
		a=b
		b=swp
	end if
	END SUBROUTINE masked_swap_rs
!BL
	SUBROUTINE masked_swap_rv(a,b,mask)
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(a)) :: swp
	where (mask)
		swp=a
		a=b
		b=swp
	end where
	END SUBROUTINE masked_swap_rv
!BL
	SUBROUTINE masked_swap_rm(a,b,mask)
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
	LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
	where (mask)
		swp=a
		a=b
		b=swp
	end where
	END SUBROUTINE masked_swap_rm

	SUBROUTINE assert1(n1,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, INTENT(IN) :: n1
	if (.not. n1) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert1'
	end if
	END SUBROUTINE assert1
!BL
	SUBROUTINE assert2(n1,n2,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, INTENT(IN) :: n1,n2
	if (.not. (n1 .and. n2)) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert2'
	end if
	END SUBROUTINE assert2
!BL
	SUBROUTINE assert3(n1,n2,n3,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, INTENT(IN) :: n1,n2,n3
	if (.not. (n1 .and. n2 .and. n3)) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert3'
	end if
	END SUBROUTINE assert3
!BL
	SUBROUTINE assert4(n1,n2,n3,n4,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, INTENT(IN) :: n1,n2,n3,n4
	if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert4'
	end if
	END SUBROUTINE assert4
!BL
	SUBROUTINE assert_v(n,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	LOGICAL, DIMENSION(:), INTENT(IN) :: n
	if (.not. all(n)) then
		write (*,*) 'nrerror: an assertion failed with this tag:', &
			string
		STOP 'program terminated by assert_v'
	end if
	END SUBROUTINE assert_v
!BL
	FUNCTION assert_eq2(n1,n2,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2
	INTEGER :: assert_eq2
	if (n1 == n2) then
		assert_eq2=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq2'
	end if
	END FUNCTION assert_eq2
!BL
	FUNCTION assert_eq3(n1,n2,n3,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2,n3
	INTEGER :: assert_eq3
	if (n1 == n2 .and. n2 == n3) then
		assert_eq3=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq3'
	end if
	END FUNCTION assert_eq3
!BL
	FUNCTION assert_eq4(n1,n2,n3,n4,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2,n3,n4
	INTEGER :: assert_eq4
	if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
		assert_eq4=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq4'
	end if
	END FUNCTION assert_eq4
!BL
	FUNCTION assert_eqn(nn,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, DIMENSION(:), INTENT(IN) :: nn
	INTEGER :: assert_eqn
	if (all(nn(2:) == nn(1))) then
		assert_eqn=nn(1)
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eqn'
	end if
	END FUNCTION assert_eqn
!BL
	SUBROUTINE nrerror(string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	write (*,*) 'nrerror: ',string
	STOP 'program terminated by nrerror'
	END SUBROUTINE nrerror
!BL
	FUNCTION arth_r(first,increment,n)
	REAL(SP), INTENT(IN) :: first,increment
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(n) :: arth_r
	INTEGER(I4B) :: k,k2
	REAL(SP) :: temp
	if (n > 0) arth_r(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_r(k)=arth_r(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_r(k)=arth_r(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
	END FUNCTION arth_r
!BL
	FUNCTION arth_d(first,increment,n)
	REAL(DP), INTENT(IN) :: first,increment
	INTEGER(I4B), INTENT(IN) :: n
	REAL(DP), DIMENSION(n) :: arth_d
	INTEGER(I4B) :: k,k2
	REAL(DP) :: temp
	if (n > 0) arth_d(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_d(k)=arth_d(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_d(k)=arth_d(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
	END FUNCTION arth_d
!BL
	FUNCTION arth_i(first,increment,n)
	INTEGER(I4B), INTENT(IN) :: first,increment,n
	INTEGER(I4B), DIMENSION(n) :: arth_i
	INTEGER(I4B) :: k,k2,temp
	if (n > 0) arth_i(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_i(k)=arth_i(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_i(k)=arth_i(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
	END FUNCTION arth_i
!BL
END MODULE nrutil

module nr1

	INTERFACE indexx
		SUBROUTINE indexx_sp(arr,index)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
		END SUBROUTINE indexx_sp
		SUBROUTINE indexx_i4b(iarr,index)
		USE nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
		END SUBROUTINE indexx_i4b
	END INTERFACE

END module nr1
!---- of modules for new numerical recipes
!------------------------------------------------------------------------------------------------
! ---- main program
program redsnapper
USE allvar
REAL*8 ran1
REAL*8 eqn9001
REAL*8 eqn8003, predh
                                                        ! Format number
OPEN(UNIT=7,FILE='rigloc1.out',STATUS='unknown')         ! 1001
open(UNIT=8,FILE='fishind1.out',STATUS='unknown')        ! Not actually being used
open(UNIT=9,FILE='neargrid1.out',STATUS='unknown')       ! 1002
open(UNIT=10,FILE='onefish1.out',STATUS='unknown')       ! 1000
open(UNIT=11,FILE='zprey1.out',STATUS='unknown')         ! 1006
open(UNIT=12,FILE='wtemp_daylight1.out',STATUS='unknown')! 1101
open(UNIT=13,FILE='age1.out',STATUS='unknown')           ! 1009
open(UNIT=14,FILE='meanwt1.out',STATUS='unknown')        ! 1010
open(UNIT=15,FILE='meanprey1.out',STATUS='unknown')      ! 1011
open(UNIT=16,FILE='biomass1.out',STATUS='unknown')       ! 1009
open(UNIT=17,FILE='bgsnapper1.out',STATUS='unknown')     ! 1012
open(UNIT=18,FILE='bioenergetics1.out',STATUS='unknown') ! 1050
open(UNIT=19,FILE='bghsnapper1.out',STATUS='unknown')    ! 1012
open(UNIT=20,FILE='yielddaily1.out',STATUS='unknown')    ! 1012
open(UNIT=21,FILE='yield1.out',STATUS='unknown')         ! 1012
open(UNIT=22,FILE='bgpin1.out',STATUS='unknown')      ! 1012
open(UNIT=23,FILE='bgcroak1.out',STATUS='unknown')     ! 1012
open(UNIT=24,FILE='bghpin1.out',STATUS='unknown')     ! 1012
open(UNIT=25,FILE='bghcroak1.out',STATUS='unknown')    ! 1012
open(UNIT=26,FILE='riskhr1.out',STATUS='unknown')        ! 1010
open(UNIT=27,FILE='rsinput.dat',STATUS='old')           ! 1051
open(UNIT=28,FILE='onefish51.out',STATUS='unknown')      ! 1052
open(UNIT=29,FILE='bioen1.dat',STATUS='old')            ! 1053 consumption
OPEN(UNIT=30,FILE='bioen2.dat',STATUS='old')            ! 1054 respiration
open(UNIT=31,FILE='bioen3.dat',STATUS='old')            ! 1055 egestion and excretion
open(UNIT=32,FILE='dailymove1.out',STATUS='unknown')     ! 1056
open(UNIT=33,FILE='weeklyrigs1.out',STATUS='unknown')    ! 1057
open(UNIT=53,FILE='weeklypyrs1.out',STATUS='unknown')    ! 1057
open(UNIT=63,FILE='weeklynats1.out',STATUS='unknown')    ! 1057
open(UNIT=34,FILE='foragedistrig1.out',STATUS='unknown')    ! 1058
open(UNIT=54,FILE='foragedistpyr1.out',STATUS='unknown')    ! 1058
open(UNIT=64,FILE='foragedistnat1.out',STATUS='unknown')    ! 1058
open(UNIT=35,FILE='bgblue1.out',STATUS='unknown')      ! 1012
open(UNIT=36,FILE='bgjack1.out',STATUS='unknown')     ! 1012
open(UNIT=37,FILE='bghblue1.out',STATUS='unknown')     ! 1012
open(UNIT=38,FILE='bghjack1.out',STATUS='unknown')    ! 1012
OPEN(UNIT=39,FILE='rigloc21.out',STATUS='unknown')         ! 1001
OPEN(UNIT=40,FILE='randcheck1.out',STATUS='unknown')     !1059
OPEN(UNIT=41,FILE='randcheck21.out',STATUS='unknown')     !1060
OPEN(UNIT=42,FILE='predcheck1.out',STATUS='unknown')      !1061
OPEN(UNIT=43,FILE='cummort1.out',STATUS='unknown')      !1062
open(UNIT=44,FILE='bgga1.out',STATUS='unknown')    ! 1012
open(UNIT=45,FILE='bghga1.out',STATUS='unknown')    ! 1012
open(UNIT=46,FILE='natreef.dat',STATUS='old')   ! 2001
open(UNIT=47,FILE='pyrreef.dat',STATUS='old')   ! 2002
open(UNIT=148,FILE='rigreefact.dat',STATUS='old')   ! 2012
open(UNIT=149,FILE='rigreefhigh.dat',STATUS='old')   ! 2012
open(UNIT=150,FILE='rigreefmed.dat',STATUS='old')   ! 2012
open(UNIT=151,FILE='rigreeflow.dat',STATUS='old')   ! 2012
open(UNIT=48,FILE='gawt_age1.out',STATUS='unknown')  ! 2003
open(UNIT=49,FILE='gctemp_check.out',STATUS='unknown') ! 2004
open(UNIT=101,FILE='costben.out',STATUS='unknown') ! 2020
!NEW - READ statement to read in model options specied from an input file and used in a batch run  MDC 16apr09
!READ (UNIT = 27, FMT = 1051) totpop(1), ipinfish, totpop(2), icroaker, totpop(3),ibluefish, totpop(4),ijack, meth, nrigs, &
!                           & intvl, bffr

!1051 FORMAT(1X, F4.0, 1X, I1, 1X, F4.0, 1X, I1, 1X, F4.0, 1X, I1, F3.0, 1X, I1, 1X, F3.0, 1X, I1, 1X, I3, 1X, I2, 1X, I1)

!PRINT *, nspecies, totpop(1), ipinfish, totpop(2), icroaker, totpop(3), meth, nrigs, intvl, bffr
!PRINT *

READ (UNIT = 29, FMT = 1053) ga(1),ga(2),ga(3),ga(4),ga(6),gb(1),gb(2),gb(3),gb(4),gb(6),gtmax(1),gtmax(2),gtmax(3),gtmax(4),gtmax(6),gtopt(1), &
                           & gtopt(2),gtopt(3),gtopt(4),gtopt(6),gtheta1(1),gtheta1(2),gtheta1(3),gtheta1(4),gtheta1(6)
                !ga        gb        gtmax & opt  gtheta1   ra         rb        rtmax & opt rtheta2
1053 FORMAT(1X,5(F6.4,1X),5(F6.3,1X),10(F4.1,1X),5(F4.2,1X))


READ (UNIT = 30,FMT =1054) ra(1),ra(2),ra(3),ra(4),ra(6),rb(1),rb(2),rb(3),rb(4),rb(6),rtmax(1),rtmax(2),rtmax(3),rtmax(4),rtmax(6),rtopt(1), &
                           & rtopt(2),rtopt(3),rtopt(4),rtopt(6),rtheta2(1),rtheta2(2),rtheta2(3),rtheta2(4),rtheta2(6),act(1),act(2),act(3),act(4),act(6)

1054 FORMAT(1X,5(F6.4,1X),5(F6.3,1X),10(F4.1,1X),5(F4.2,1X),5(F4.2,1X))

READ (UNIT = 31, FMT = 1055) fa(1),fa(2),fa(3),fa(4),fa(6),fb(1),fb(2),fb(3),fb(4),fb(6),fg(1),fg(2),fg(3),fg(4),fg(6),ua(1),ua(2),ua(3),ua(4),ua(6), &
                           & ub(1),ub(2),ub(3),ub(4),ub(6),ug(1),ug(2),ug(3),ug(4),ug(6)
                !fa        fb         fg         ua         ub         ug
1055 FORMAT(1X,5(F5.3,1X),5(F6.3,1X),5(F5.3,1X),5(F5.3,1X),5(F4.2,1X),5(F6.3,1X))

READ(UNIT = 46, FMT=2001) natcol(1),natcol(2),natcol(3),natcol(4),natcol(5),natcol(6),natcol(7),natcol(8),natcol(9),natcol(10),natcol(11),natcol(12),natcol(13),natcol(14),natcol(15),natcol(16),natcol(17),natcol(18), &
                           & natcol(19),natcol(20),natcol(21),natcol(22),natcol(23),natcol(24),natcol(25),natcol(26),natcol(27),natcol(28),natcol(29),natcol(30), &
                           & natcol(31),natcol(32),natcol(33),natcol(34),natcol(35),natcol(36),natcol(37),natcol(38),natcol(39),natcol(40),natcol(41),natcol(42),natcol(43), natcol(44),natcol(45),natcol(46),natcol(47),natcol(48), &
                           & natcol(49),natcol(50),natcol(51),natcol(52),natcol(53),natcol(54),natcol(55),natcol(56),natcol(57),natcol(58),natcol(59),natcol(60), &
                           & natcol(61),natcol(62),natcol(63),natcol(64),natcol(65),natcol(66),natcol(67),natcol(68),natcol(69),natcol(70),natcol(71),natcol(72),natcol(73),natcol(74),natcol(75),natcol(76),natcol(77),natcol(78), &
                           & natcol(79),natcol(80),natcol(81),natcol(82),natcol(83),natcol(84),natcol(85),natcol(86),natcol(87),natcol(88),natcol(89),natcol(90), &
                           & natcol(91),natcol(92),natcol(93),natcol(94),natcol(95),natcol(96),natcol(97),natcol(98),natcol(99),natcol(100),natcol(101),natcol(102),natcol(103),natcol(104),natcol(105),natcol(106),natcol(107),natcol(108), &
                           & natcol(109),natcol(110),natcol(111),natcol(112),natcol(113),natcol(114),natcol(115),natcol(116),natcol(117),natcol(118),natcol(119),natcol(120), &
                           & natcol(121),natcol(122),natcol(123),natcol(124),natcol(125),natcol(126),natcol(127),natcol(128), &
                           & natrow(1),natrow(2),natrow(3),natrow(4),natrow(5),natrow(6),natrow(7),natrow(8),natrow(9),natrow(10),natrow(11),natrow(12),natrow(13),natrow(14),natrow(15),natrow(16),natrow(17),natrow(18), &
                           & natrow(19),natrow(20),natrow(21),natrow(22),natrow(23),natrow(24),natrow(25),natrow(26),natrow(27),natrow(28),natrow(29),natrow(30), &
                           & natrow(31),natrow(32),natrow(33),natrow(34),natrow(35),natrow(36),natrow(37),natrow(38),natrow(39),natrow(40),natrow(41),natrow(42),natrow(43), natrow(44),natrow(45),natrow(46),natrow(47),natrow(48), &
                           & natrow(49),natrow(50),natrow(51),natrow(52),natrow(53),natrow(54),natrow(55),natrow(56),natrow(57),natrow(58),natrow(59),natrow(60), &
                           & natrow(61),natrow(62),natrow(63),natrow(64),natrow(65),natrow(66),natrow(67),natrow(68),natrow(69),natrow(70),natrow(71),natrow(72),natrow(73),natrow(74),natrow(75),natrow(76),natrow(77),natrow(78), &
                           & natrow(79),natrow(80),natrow(81),natrow(82),natrow(83),natrow(84),natrow(85),natrow(86),natrow(87),natrow(88),natrow(89),natrow(90), &
                           & natrow(91),natrow(92),natrow(93),natrow(94),natrow(95),natrow(96),natrow(97),natrow(98),natrow(99),natrow(100),natrow(101),natrow(102),natrow(103),natrow(104),natrow(105),natrow(106),natrow(107),natrow(108), &
                           & natrow(109),natrow(110),natrow(111),natrow(112),natrow(113),natrow(114),natrow(115),natrow(116),natrow(117),natrow(118),natrow(119),natrow(120), &
                           & natrow(121),natrow(122),natrow(123),natrow(124),natrow(125),natrow(126),natrow(127),natrow(128) 
2001 FORMAT(1x,256(i2,1x))
                   
READ(UNIT = 47, FMT=2002) pyrcol(1),pyrcol(2),pyrcol(3),pyrrow(1),pyrrow(2),pyrrow(3)
2002 FORMAT(1x,6(i2,1x))

PRINT *, 'Which rig method?'
PRINT *, '1 - Actual, 2 - High, 3 - Medium, 4 - Low'
READ *, rigmeth
IF(rigmeth.eq.1)THEN
nrigs = 12
READ(UNIT = 148, FMT=2012) rigcol(1), rigcol(2), rigcol(3), rigcol(4), rigcol(5), rigcol(6), rigcol(7), rigcol(8), rigcol(9), rigcol(10), &
                            &   rigcol(11), rigcol(12), &
                            &   rigrow(1), rigrow(2), rigrow(3), rigrow(4), rigrow(5), rigrow(6), rigrow(7), rigrow(8), rigrow(9), rigrow(10), &
                            &   rigrow(11), rigrow(12)
2012 FORMAT(1x,24(i2,1x))
                            
ELSE IF(rigmeth.eq.2)THEN
nrigs = 20
     READ(UNIT = 149, FMT=2022) rigcol(1), rigcol(2), rigcol(3), rigcol(4), rigcol(5), rigcol(6), rigcol(7), rigcol(8), rigcol(9), rigcol(10), &
                            &   rigcol(11), rigcol(12), rigcol(13), rigcol(14), rigcol(15), rigcol(16), rigcol(17), rigcol(18), rigcol(19), rigcol(20), &
                            &   rigrow(1), rigrow(2), rigrow(3), rigrow(4), rigrow(5), rigrow(6), rigrow(7), rigrow(8), rigrow(9), rigrow(10), &
                            &   rigrow(11), rigrow(12), rigrow(13), rigrow(14), rigrow(15), rigrow(16), rigrow(17), rigrow(18), rigrow(19), rigrow(20)
2022 FORMAT(1x,40(i2,1x))
                            
ELSE IF(rigmeth.eq.3)THEN
nrigs = 16
READ(UNIT = 150, FMT=2032) rigcol(1), rigcol(2), rigcol(3), rigcol(4), rigcol(5), rigcol(6), rigcol(7), rigcol(8), rigcol(9), rigcol(10), &
                            &   rigcol(11), rigcol(12), rigcol(13), rigcol(14), rigcol(15), rigcol(16), &
                            &   rigrow(1), rigrow(2), rigrow(3), rigrow(4), rigrow(5), rigrow(6), rigrow(7), rigrow(8), rigrow(9), rigrow(10), &
                            &   rigrow(11), rigrow(12), rigrow(13), rigrow(14), rigrow(15), rigrow(16)
2032 FORMAT(1x,32(i2,1x))
                            
ELSE IF(rigmeth.eq.4)THEN
nrigs = 9
READ(UNIT = 151, FMT=2042) rigcol(1), rigcol(2), rigcol(3), rigcol(4), rigcol(5), rigcol(6), rigcol(7), rigcol(8), rigcol(9), rigrow(1), rigrow(2), rigrow(3), &
                            &   rigrow(4), rigrow(5), rigrow(6), rigrow(7), rigrow(8), rigrow(9)
2042 FORMAT(1x,18(i2,1x))
END IF


!Check the correctness of the bioenergetics data file read in
do i=1,nspecies
PRINT *,i,ga(i),gb(i),gtmax(i),gtopt(i),gtheta1(i),ra(i),rb(i),rtmax(i),rtopt(i),rtheta2(i),act(i),fa(i),fb(i),fg(i),ua(i),ub(i),ug(i)
END do 
pause

! Starting off with all of the cells as benthic cells
do i=1,ncol
   do j=1,nrow
      zhab(i,j)=1
   end do
end do

! manually entering then densities and weightings for each of the pyramid, rig, and nat cells KM 5-11
PRINT *, 'Enter the Pyramid Density for Pyr Site 1'
READ *, pyrdens(1)

PRINT *, 'Enter the Pyramid Density for Pyr Site 2'
READ *, pyrdens(2)

PRINT *, 'Enter the Pyramid Density for Pyr Site 3'
READ *, pyrdens(3)

PRINT *, 'Enter rig distance weightings'
READ *, rigweight

PRINT *, 'Enter natural bank distance weightings'
READ *, natweight

PRINT *, 'Enter pyramid distance weightings for sites 1, 2, and 3'
READ *, pyrweight(1)
READ *, pyrweight(2)
READ *, pyrweight(3)


! check the correctness of the reef files  
do i=1,numnatcells
    PRINT*, i,natcol(i),natrow(i)
    zhab(natcol(i),natrow(i))=3
end do


do i=1,npyr
   PRINT*, i,pyrcol(i),pyrrow(i),pyrdens(i)
   zhab(pyrcol(i),pyrrow(i))=2
end do


do i=1,nrigs
   PRINT*, i,rigcol(i),rigrow(i)
   zhab(rigcol(i),rigrow(i))=0
end do

! Use either read from file (above) or read from screen here
 !PRINT *,'enter total pop of red snapper on each rig in numbers'  ! number of age-1 individual snapper in population
 !READ *,totpop(1)
 totpop(1)=1000
 
 !PRINT *,'enter 0 for no pinfish and 1 for pinfish'  ! 1 means the second species is grunts
 !READ *,ipinfish
 ipinfish=1

 IF(ipinfish.eq.1)then
  !PRINT *,'enter total pop of pinfish on each rig in numbers'
  !READ *,totpop(2)
  totpop(2)=1000
 endif
!
 !PRINT *,'enter 0 for no Atl croaker and 1 for Atl croaker'  ! 1 means the third species is grouper
 !READ *,icroaker
 icroaker=1
!
 IF(icroaker.eq.1)then
  !PRINT *,'enter total pop of Atl croaker on each rig in numbers'
  !READ *,totpop(3)
  totpop(3)=1000
 endif
!
!!*NEW - entered setup information to include a 4th IBM species ... bluefish
  !PRINT *,'enter 0 for no bluefish and 1 for bluefish'  ! 1 means the fourth species is bluefish
  !READ *,ibluefish
  ibluefish=1
!
 IF(ibluefish.eq.1)then
  !PRINT *,'enter total pop of bluefish on each rig in numbers'
  !READ *,totpop(4)
  totpop(4)=1000
 endif
!
  !PRINT *,'enter 0 for no jack and 1 for jack'  ! 1 means the fifth species is jack
  !READ *,ijack
  ijack=1
!
 IF(ijack.eq.1)then
  !PRINT *,'enter total pop of jack on each rig in numbers'
  !READ *,totpop(5)
  totpop(5)=1000
 endif
!
!PRINT *, 'Enter the total pop of Greater Amberjack on each rig in numbers'
!READ *, totpop(6)
totpop(6)=1000

call parameters   ! assign values to many model parameters
call setupgrid        ! set up the grid and environmental variables

sumfish=0
call initfish(1)        ! set up initial population for snapper
IF(ipinfish.eq.1)call initfish(2)   ! if grunts, then set up initial pop for grunts
IF(icroaker.eq.1)call initfish(3)  ! if grouper, set up initial pop for ggouper
IF(ibluefish.eq.1)call initfish(4) ! if bluefish, set up initial pop for bluefish
IF(ijack.eq.1)call initfish(5)
call initfish(6) ! set up initial fish for greater amberjack

!do i=1,itotf
!if (MOD(i,10).eq.0) then
!  PRINT *, i,xspecies(i),xweight(i),xworth(i)
!  pause
!END if
!END do

cumeat=0.0       ! for first call to preyupdate
call preyupdate  ! update prey to get started


do 100 iyear=1,10   ! loop over years

!if (MOD(iyear,10).eq.0.and.newrigs.gt.0)call addrig
!if (iyear.eq.5.and.remrigs.gt.0)call subrig   ! subtract remrigs at year 5
!if (iyear.eq.5.and.newrigs.gt.0)call addrig   ! add newrigs at year 5

  do 200 iday=1,365    ! loop over days in a year
     cumdays = (iyear-1)*365 + iday   ! days from 1 to end of run for plotting 1,2,3 ......,678,....,1254, etc.
     jday=iday+181        ! iday=1 is jday=182 which is July 1
     IF(jday.gt.365)jday=jday-365  ! wrap around Dec 31 for jday

     PRINT *,'iyear iday jday',iyear,iday,jday

      dailyc=0.0
      yield=0.0  ! summed hourly biomass for each day due to fishing mortality
      yieldnum=0.0; yieldwt=0.0
      naturalp=0.0; growthp=0.0; naturalpnum=0.0
      xriskhr=0
      xdailymove=0.0
      xcmort=0.0

      

      watertemp=eqn8003(dfloat(jday))   ! cal function to get water temperature for today

      ! based off of Figure 6 from Turner et al. 2017 - KM 3-4.
      ! This equation will predict the water temperature for each cell and assumes that the far left of the grid is 10 m and that each cell is a 5 m depth contour
      ! This equation assumed that there is an increase in 0.13 m with every 100 m depth (i.e. every cell)
      tempcalib = -3.801*(LOG(0.13*ncol/2+20))+37.098  ! this is the temperature function created from Turner et al. 2017. This is to calibrate the temperature so that it represents realistic temperatures that revolve around the depth at the center of the grid. This is also assuming that the first column represents 20 m depth
      b=37.098+(watertemp-tempcalib)    
      do i=1,ncol
         do j=1,nrow
      depth(i)=0.13*i+20 ! 0.13 represents a realistic increase in depth every 100 m (or every cell)
      watertempbycol(i,j)=-3.801*LOG(depth(i))+b
         end do
      end do

!      IF(iout.ge.1)then
!      WRITE(12,1007)cumdays,iyear,iday,jday,watertemp,light_hours
! 1007 FORMAT(1x,f12.4,1x,3(i7,1x),f12.4,1x,f12.4,1x)
!      endif
                

      do 300 ihour=1,24   ! loop over hours
          cumhours = (iyear-1)*365*24 + (iday-1)*24 + ihour   ! cumulative hours for plotting
          cumeat=0.0   ! zero matrix for every hour since prey is updated after each hour
          ! Determining the number of light hours in a day. Equation is based off of https://gml.noaa.gov/grad/solcalc/calcdetails.html   KM 1-28-22
          light_hours = eqn9001(dfloat(jday)) ! determining the number of light hours in a day throughout the year 
          IF(ihour.le.light_hours)THEN
             inight=2
             ELSE
             inight=1
          END IF
          IF(iyear.eq.1)THEN
             WRITE(12,1101)cumdays,iyear,iday,jday,ihour,watertemp,light_hours,inight
1101         FORMAT(1x,f12.4,1x,4(i7,1x),f12.4,1x,f12.4,1x,i7,1x)
          END IF
         do 400 ifish=1,itotf    ! loop over all fish
           !IF(xspecies.eq.1) nhood=eqn9002(dfloat(xweight)) ! equation for the change in foraging area with weight
           IF(xalive(ifish).eq.0)then   ! if alive then continue
             wtemp(ifish)=watertempbycol(xxcell(ifish),xycell(ifish))
                 
             IF(inight.eq.2)then   ! if day then forage
              call moveday      ! move hourly during day
             if (xspecies(ifish).ne.5) then  ! no bioenergetics for jacks strictly movement
                call growth(1,xxcell(ifish),xycell(ifish))   ! grow the individual in location from last hour

                  IF(xspecies(ifish).eq.4)grow=0.0  ! bluefish don't grow but they do consume/compete
                 dailyc(ifish)=dailyc(ifish)+totcon(ifish)    ! g prey/g pred summed over hours in a day
                 dailyp(ifish)=dailyc(ifish)/gcmax
                 !                                                 (g prey/g fish)*(g fish/ind)*(num of pop ind)
                 growthp(xspecies(ifish)) = growthp(xspecies(ifish)) + grow*xweight(ifish)*xworth(ifish) ! growth production

                 xweight(ifish)=xweight(ifish)+grow*xweight(ifish) ! update body weight (g ww per ind) of model individual     
                
                 do k=1,nprey   ! loop over prey types
                   xeat(ifish,k)=geat(k) ! store g prey/g pred eaten per hour of prey type k for individual ifish
                   cumeat(xxcell(ifish),xycell(ifish),k)=cumeat(xxcell(ifish),xycell(ifish),k)+&
                     & geat(k)*xweight(ifish)*xworth(ifish)   ! g prey/g pred/h to g prey by the population worth of all species
                 end do   
                 
                 IF(xweight(ifish).lt.1.0)then   ! if less than one gram wet weight, then starved
                  xweight(ifish)=1.0
                  xalive(ifish)=2    ! starved flag
                 endif !end if for xweight.lt.1.0              
              endif  ! for xsp ne 5
            endif  ! for inight=2

             if(inight.eq.1)then  ! if night then move towards a reef cell 
              if(xspecies(ifish).le.3.or.xspecies(ifish).eq.6) then
               if(xhab(ifish).eq.1.and.xspecies(ifish).le.3) call movenight    ! move to nearest rig, pyr, or nat cell - no growth but mort does occurs (no end if statement needed in this form)
                    call respiration   ! respiration for this hour
                avgden=(ener(1)+ener(2)+ener(3)+ener(4)+ener(5))/5.0   ! default if nothing eaten, otherwise would get zero
!               g ww/ind    = g ww/ind  -   g equivalents of prey respired in an hour * (cal/g prey)/(cal/g pred) * g ww/ind
                xweight(ifish)=xweight(ifish)-(resp/24.0)*(avgden/en(xspecies(ifish)))*xweight(ifish) ! decrement weight for hourly respiration
!               g ww          =    g ww - loss due to respiration but times pop worth of each individual
                growthp(xspecies(ifish))=growthp(xspecies(ifish))-(resp/24.0)*avgden/en(xspecies(ifish))*xweight(ifish)*xworth(ifish)

                IF(xhab(ifish).eq.1.and.xspecies(ifish).le.3)xriskhr(ifish)=xriskhr(ifish)+1   ! number of hours off rig during nighttime
              end if ! ifish.le.3
              !bluefish & jacks are going to move the same regardless of time of day, move according to biomass of species 1-3
              !bluefish & jacks do not respire
              IF(xspecies(ifish).ge.4) call moveday

            endif ! for inight=1

             IF(iyear.eq.9.AND.(iday.gt.200.and.iday.lt.220))then
                 WRITE(18,1050)iyear,iday,ihour,inight,ifish,xspecies(ifish),xage(ifish),xhab(ifish),xxcell(ifish),xycell(ifish),&
                          &  xxloc(ifish),xyloc(ifish),xweight(ifish),gcmax, &
                          &  totcon(ifish),dailyc(ifish),dailyp(ifish),resp/12.0,egest,exc,egg/12.0,grow
             1050    FORMAT(1x,10(i5,1x),12(f12.4,1x))
              endif !end if for write statement unit 18

             if(xspecies(ifish).le.3.or.xspecies(ifish).eq.6)call mortality  

          ENDIF    ! end if for xalive=0

           IF(iout.eq.1)then
           !changed year to 5 from 10 .... 23 june 2009.  Needs to be changed back with full model.
            IF(iyear.ge.1.and.(iday.eq.1.or.iday.eq.182).and.MOD(ifish,100).eq.0)then
      !     IF(MOD(iyear,5).eq.0.AND.MOD(ifish,100).eq.0.and.(iday.eq.10.or.iday.eq.11.or.iday.eq.12.or.iday.eq.13))then
           !IF(iyear.eq.2.AND.MOD(ifish,100).eq.0.and.(iday.eq.10.or.iday.eq.11.or.iday.eq.12.or.iday.eq.13))then
                
            WRITE(10,1000)cumhours,iyear,iday,jday,ihour,inight,xage(ifish),xspecies(ifish),ifish,xalive(ifish),xxcell(ifish),&
                 &   xycell(ifish),xxloc(ifish),xyloc(ifish),xhab(ifish),xweight(ifish),(geat(kk),kk=1,nprey)  ! geat is g prey/g pred/hour
            
      !      WRITE(10,1000)cumhours,iyear,iday,jday,ihour,inight,xspecies(ifish),ifish,xalive(ifish),xxcell(ifish),&
      !           &   xycell(ifish),xxloc(ifish),xyloc(ifish),xhab(ifish),xweight(ifish),(geat(kk),kk=1,nprey)  ! geat is g prey/g pred/hour
            
           endif !end if write statement unit 10
           1000      FORMAT(1x,e12.3,1x,7(i5,1x),i10,1x,3(i4,1x),1x,2(f12.4,1x),i3,1x,10(f12.4,1x))

!           IF(MOD(iyear,20).eq.0.AND.MOD(ifish,100).eq.0.and.MOD(iday,30).eq.0) then
!           WRITE(28,1052)cumhours,iyear,iday,jday,ihour,inight,xspecies(ifish),ifish,xalive(ifish),xxcell(ifish),&
!           &   xycell(ifish),xrig(ifish),xweight(ifish),(geat(kk),kk=1,nprey),gcmax,totcon(ifish)  ! geat is g prey/g pred/hour
!           endif !end if for write unit 28
!           1052      FORMAT(1x,e12.3,1x,6(i5,1x),i10,1x,3(i4,1x),1x,i3,1x,10(f12.4,1x))

           IF(iyear.eq.10.and.iday.eq.1.and.(xspecies(ifish).eq.6.or.xspecies(ifish).le.3)) THEN
              WRITE(48,2003) iyear,iday,xspecies(ifish),xage(ifish),xweight(ifish),xalive(ifish)
           2003 FORMAT(1x,4(i5,1x),f12.4,1x,i5,1x)
           END IF !for GA weights
           
           ! problem with gctemp failing, checking to see what the numbers are like the day before and the day of the failure
           !IF(xspecies(ifish).ne.5.and.(iday.eq.29.or.iday.eq.30).and.MOD(ifish,100).eq.0) then
              ! WRITE(49,2004) iyear, iday, ihour, xspecies(ifish), ifish, gctemp
!2004           FORMAT(1x,5(i5,1x),f12.4,1x)     
              ! end if !for gctemp
           endif !iout.eq.1

400      continue   ! end of loop over fish

!  WRITE(21,1012)cumhours,iyear,iday,jday,ihour,inight,yield(1),yieldnum(1),yieldwt(1), &   ! temporary to see hourly yield
!        & yield(2),yieldnum(2),yieldwt(2),yield(3),yieldnum(3),yieldwt(3)
!  1012  FORMAT(1x,e12.4,1x,5(i5,1x),9(e12.5,1x))

         call preyupdate
         call bioout

300   continue    ! end of loop over hours

      call dailyout

200 continue   ! end of loop over days

    call age   ! age individuals for species 1,2,and 3, remove old age individuals, and introduce imodf
               ! new model individuals corresponding to recruits (new age1)

100 continue   ! end of loop over years

end program
!****************************************************
subroutine parameters
USE allvar
cellsize=100   ! meters

! 5 prey types  1-zoop 2-crabs 3-shrimp 4-pelagic-fish 5-benthic-fish

kkprey(1)=6.7           ! zoop   g ww/m2
kkprey(2)=4.0           ! crabs  tons/km2   1000 kg/mt * 1000 g/kg km2 to m2 1000.0*1000.0/(1000.0*1000.0)
kkprey(3)=4.0           ! shrimp
kkprey(4)=15.0          ! pelagic fish
kkprey(5)=15.0          ! benthic fish

rrprey(1)=17.3/365.0/24.0  ! P/B in per year to per hour - used 24 because preyupdate called 24 times per day 12? KAR
rrprey(2)=4.0/365.0/24.0
rrprey(3)=4.0/365.0/24.0
rrprey(4)=4.0/365.0/24.0
rrprey(5)=4.0/365.0/24.0

! Ask Matt about the units of these values - should be in J/g and might be but the units below just could say the wrong thing
ener(1) = 3511.0       ! cal/gdw by prey type were originally given
ener(2) = 3138.0       ! average ed = 4087 cal/gdw
ener(3) = 3894.0
ener(4) = 4947.0
ener(5) = 4947.0

! species 1 - red snapper, 2 - pinfish, 3 - Atlantic croaker, 4 - Bluefish, and 5 - Jack spp
! prey 1 - zooplankton, 2 - crabs, 3 - shrimp, 4 - pelagic fish, and 5 - benthic fish
! ages 1 - 10

vul(1,1,1)=0.3;   vul(1,1,2)=0.3;   vul(1,1,3)=0.3;   vul(1,1,4)=0.3;   vul(1,1,5)=0.3   ! zoop being eaten by life stage 1 to 5
vul(1,2,1)=0.1;   vul(1,2,2)=0.1;   vul(1,2,3)=0.1;   vul(1,2,4)=0.1;   vul(1,2,5)=0.1   ! crabs being eaten by stage 2
vul(1,3,1)=0.0;   vul(1,3,2)=0.0;   vul(1,3,3)=0.0;   vul(1,3,4)=0.0;   vul(1,3,5)=0.0   ! shrimp being eaten by stage 3
vul(1,4,1)=0.15;  vul(1,4,2)=0.15;  vul(1,4,3)=0.15;  vul(1,4,4)=0.15;  vul(1,4,5)=0.15   ! pelagic fish being eaten by stage 4
vul(1,5,1)=0.03;  vul(1,5,2)=0.03;  vul(1,5,3)=0.03;  vul(1,5,4)=0.03;  vul(1,5,5)=0.03   ! benthic fish being eaten by stage 5

vul(1,1,6)=0.3 ;  vul(1,1,7)=0.3;   vul(1,1,8)=0.3;   vul(1,1,9)=0.3;   vul(1,1,10)=0.3   ! zoop being eaten by life stage 1 to 5
vul(1,2,6)=0.1 ;  vul(1,2,7)=0.1;   vul(1,2,8)=0.1;   vul(1,2,9)=0.1;   vul(1,2,10)=0.1   ! crabs being eaten by stage 2
vul(1,3,6)=0.0 ;  vul(1,3,7)=0.0;   vul(1,3,8)=0.0;   vul(1,3,9)=0.0;   vul(1,3,10)=0.0   ! shrimp being eaten by stage 3
vul(1,4,6)=0.15;  vul(1,4,7)=0.15;  vul(1,4,8)=0.15;  vul(1,4,9)=0.15;  vul(1,4,10)=0.15   ! pelagic fish being eaten by stage 4
vul(1,5,6)=0.03;  vul(1,5,7)=0.03;  vul(1,5,8)=0.03;  vul(1,5,9)=0.03;  vul(1,5,10)=0.03   ! benthic fish being eaten by stage 5

! species 2
vul(2,1,1)=1.0;   vul(2,1,2)=1.0;   vul(2,1,3)=1.0;   vul(2,1,4)=1.0;   vul(2,1,5)=1.0   ! zoop
vul(2,2,1)=1.0;   vul(2,2,2)=1.0;   vul(2,2,3)=1.0;   vul(2,2,4)=1.0;   vul(2,2,5)=1.0   ! crabs
vul(2,3,1)=1.0;   vul(2,3,2)=1.0;   vul(2,3,3)=1.0;   vul(2,3,4)=1.0;   vul(2,3,5)=1.0   ! shrimp
vul(2,4,1)=0.0;   vul(2,4,2)=0.0;   vul(2,4,3)=0.0;   vul(2,4,4)=0.0;   vul(2,4,5)=0.0 ! pelagic fish
vul(2,5,1)=0.0;   vul(2,5,2)=0.0;   vul(2,5,3)=0.0;   vul(2,5,4)=0.0;   vul(2,5,5)=0.0   ! benthic fish

vul(2,1,6)=1.0;   vul(2,1,7)=1.0;   vul(2,1,8)=1.0;   vul(2,1,9)=1.0;   vul(2,1,10)=1.0   ! zoop
vul(2,2,6)=1.0;   vul(2,2,7)=1.0;   vul(2,2,8)=1.0;   vul(2,2,9)=1.0;   vul(2,2,10)=1.0   ! crabs
vul(2,3,6)=1.0;   vul(2,3,7)=1.0;   vul(2,3,8)=1.0;   vul(2,3,9)=1.0;   vul(2,3,10)=1.0   ! shrimp
vul(2,4,6)=0.0;   vul(2,4,7)=0.0;   vul(2,4,8)=0.0;   vul(2,4,9)=0.0;   vul(2,4,10)=0.0   ! pelagic fish
vul(2,5,6)=0.0;   vul(2,5,7)=0.0;   vul(2,5,8)=0.0;   vul(2,5,9)=0.0;   vul(2,5,10)=0.0   ! benthic fish

! species 3
vul(3,1,1)=0.0;   vul(3,1,2)=0.0;   vul(3,1,3)=0.0;   vul(3,1,4)=0.0;   vul(3,1,5)=0.0   ! zoop
vul(3,2,1)=0.6;   vul(3,2,2)=0.6;   vul(3,2,3)=0.6;   vul(3,2,4)=0.6;   vul(3,2,5)=0.6   ! crabs
vul(3,3,1)=0.4;   vul(3,3,2)=0.4;   vul(3,3,3)=0.4;   vul(3,3,4)=0.4;   vul(3,3,5)=0.4   ! shrimp
vul(3,4,1)=0.05;  vul(3,4,2)=0.05;  vul(3,4,3)=0.05;  vul(3,4,4)=0.05;  vul(3,4,5)=0.05   ! pelagic fish
vul(3,5,1)=0.15;  vul(3,5,2)=0.15;  vul(3,5,3)=0.15;  vul(3,5,4)=0.15;  vul(3,5,5)=0.15   ! benthic fish

vul(3,1,6)=0.0;   vul(3,1,7)=0.0;   vul(3,1,8)=0.0;   vul(3,1,9)=0.0;   vul(3,1,10)=0.0   ! zoop
vul(3,2,6)=0.6;   vul(3,2,7)=0.6;   vul(3,2,8)=0.6;   vul(3,2,9)=0.6;   vul(3,2,10)=0.6   ! crabs
vul(3,3,6)=0.4;   vul(3,3,7)=0.4;   vul(3,3,8)=0.4;   vul(3,3,9)=0.4;   vul(3,3,10)=0.4   ! shrimp
vul(3,4,6)=0.05;  vul(3,4,7)=0.05;  vul(3,4,8)=0.05;  vul(3,4,9)=0.05;  vul(3,4,10)=0.05   ! pelagic fish
vul(3,5,6)=0.15;  vul(3,5,7)=0.15;  vul(3,5,8)=0.15;  vul(3,5,9)=0.15;  vul(3,5,10)=0.15   ! benthic fish

! species 4 ... currently not called during growth - movement only
vul(4,1,1)=0.0;   vul(4,1,2)=0.0;   vul(4,1,3)=0.0;   vul(4,1,4)=0.0;   vul(4,1,5)=0.0   ! zoop
vul(4,2,1)=0.0;   vul(4,2,2)=0.0;   vul(4,2,3)=0.0;   vul(4,2,4)=0.0;   vul(4,2,5)=0.0   ! crabs
vul(4,3,1)=0.0;   vul(4,3,2)=0.0;   vul(4,3,3)=0.0;   vul(4,3,4)=0.0;   vul(4,3,5)=0.0   ! shrimp
vul(4,4,1)=0.2;   vul(4,4,2)=0.2;   vul(4,4,3)=0.2;   vul(4,4,4)=0.2;   vul(4,4,5)=0.2   ! pelagic fish
vul(4,5,1)=0.2;   vul(4,5,2)=0.2;   vul(4,5,3)=0.2;   vul(4,5,4)=0.2;   vul(4,5,5)=0.2   ! benthic fish

vul(4,1,6)=0.0;   vul(4,1,6)=0.0;   vul(4,1,7)=0.0;   vul(4,1,8)=0.0;   vul(4,1,10)=0.0   ! zoop
vul(4,2,6)=0.0;   vul(4,2,6)=0.0;   vul(4,2,7)=0.0;   vul(4,2,8)=0.0;   vul(4,2,10)=0.0   ! crabs
vul(4,3,6)=0.0;   vul(4,3,6)=0.0;   vul(4,3,7)=0.0;   vul(4,3,8)=0.0;   vul(4,3,10)=0.0   ! shrimp
vul(4,4,6)=0.2;   vul(4,4,6)=0.2;   vul(4,4,7)=0.2;   vul(4,4,8)=0.2;   vul(4,4,10)=0.2   ! pelagic fish
vul(4,5,6)=0.2;   vul(4,5,6)=0.2;   vul(4,5,7)=0.2;   vul(4,5,8)=0.2;   vul(4,5,10)=0.2   ! benthic fish

! Species 6 - Greater Amberjack
vul(6,1,1)=0.0;    vul(6,1,2)=0.0;    vul(6,1,3)=0.0;    vul(6,1,4)=0.0;    vul(6,1,5)=0.0     ! zoop being eaten by life stage 1 to 5
vul(6,2,1)=0.0;    vul(6,2,2)=0.0;    vul(6,2,3)=0.0;    vul(6,2,4)=0.0;    vul(6,2,5)=0.0     ! crabs being eaten by stage 2
vul(6,3,1)=0.1;    vul(6,3,2)=0.1;    vul(6,3,3)=0.1;    vul(6,3,4)=0.1;    vul(6,3,5)=0.1     ! shrimp being eaten by stage 3
vul(6,4,1)=0.2;    vul(6,4,2)=0.2;    vul(6,4,3)=0.2;    vul(6,4,4)=0.2;    vul(6,4,5)=0.2     ! pelagic fish being eaten by stage 4
vul(6,5,1)=0.25;   vul(6,5,2)=0.25;   vul(6,5,3)=0.25;   vul(6,5,4)=0.25;   vul(6,5,5)=0.25    ! benthic fish being eaten by stage 5

vul(6,1,6)=0.0;    vul(6,1,7)=0.0;    vul(6,1,8)=0.0;    vul(6,1,9)=0.0;    vul(6,1,10)=0.0     ! zoop being eaten by life stage 1 to 5
vul(6,2,6)=0.0;    vul(6,2,7)=0.0;    vul(6,2,8)=0.0;    vul(6,2,9)=0.0;    vul(6,2,10)=0.0     ! crabs being eaten by stage 2
vul(6,3,6)=0.1;    vul(6,3,7)=0.1;    vul(6,3,8)=0.1;    vul(6,3,9)=0.1;    vul(6,3,10)=0.1     ! shrimp being eaten by stage 3
vul(6,4,6)=0.2;    vul(6,4,7)=0.2;    vul(6,4,8)=0.2;    vul(6,4,9)=0.2;    vul(6,4,10)=0.2     ! pelagic fish being eaten by stage 4
vul(6,5,6)=0.25;   vul(6,5,7)=0.25;   vul(6,5,8)=0.25;   vul(6,5,9)=0.25;   vul(6,5,10)=0.25    ! benthic fish being eaten by stage 5


! half-saturation parameters of  functional response - all in units g ww/m2
! ksat's have to be a number greater than zero because they are used as a denominator in the functional response
! species 1 - red snapper
ksat(1,1,1)=155.0;   ksat(1,1,2)=160.0;   ksat(1,1,3)=165.0;   ksat(1,1,4)=170.0;  ksat(1,1,5)=175.0   ! zoop being eaten by five life stages
ksat(1,2,1)=180.0;   ksat(1,2,2)=155.0;   ksat(1,2,3)=160.0;   ksat(1,2,4)=165.0;  ksat(1,2,5)=170.0   !
ksat(1,3,1)=5.0;     ksat(1,3,2)=5.0;     ksat(1,3,3)=5.0;     ksat(1,3,4)=5.0;    ksat(1,3,5)=5.0   !
ksat(1,4,1)=120.0;   ksat(1,4,2)=120.0;   ksat(1,4,3)=120.0;   ksat(1,4,4)=120.0;  ksat(1,4,5)=120.0     !
ksat(1,5,1)=130.0;   ksat(1,5,2)=130.0;   ksat(1,5,3)=130.0;   ksat(1,5,4)=130.0;  ksat(1,5,5)=130.0   !

ksat(1,1,6)=180.0;   ksat(1,1,7)=185.0;   ksat(1,1,8)=190.0;   ksat(1,1,9)=195.0;  ksat(1,1,10)=200.0   ! zoop being eaten by five life stages
ksat(1,2,6)=175.0;   ksat(1,2,7)=180.0;   ksat(1,2,8)=185.0;   ksat(1,2,9)=190.0;  ksat(1,2,10)=195.0   !
ksat(1,3,6)=5.0;     ksat(1,3,7)=5.0;     ksat(1,3,8)=5.0;     ksat(1,3,9)=5.0;    ksat(1,3,10)=5.0   !
ksat(1,4,6)=120.0;   ksat(1,4,7)=120.0;   ksat(1,4,8)=120.0;   ksat(1,4,9)=120.0;  ksat(1,4,10)=120.0     !
ksat(1,5,6)=130.0;   ksat(1,5,7)=130.0;   ksat(1,5,8)=130.0;   ksat(1,5,9)=130.0;  ksat(1,5,10)=130.0   !

! species 2
ksat(2,1,1)=196.0;   ksat(2,1,2)=196.5;  ksat(2,1,3)=197.0;   ksat(2,1,4)=197.5;  ksat(2,1,5)=198.0   !
ksat(2,2,1)=196.0;   ksat(2,2,2)=196.5;  ksat(2,2,3)=197.0;   ksat(2,2,4)=197.5;  ksat(2,2,5)=198.0   !
ksat(2,3,1)=196.0;   ksat(2,3,2)=196.5;  ksat(2,3,3)=197.0;   ksat(2,3,4)=197.5;  ksat(2,3,5)=198.0   !
ksat(2,4,1)=20.0;    ksat(2,4,2)=20.0;   ksat(2,4,3)=20.0;    ksat(2,4,4)=20.0;   ksat(2,4,5)=20.0   !
ksat(2,5,1)=20.0;    ksat(2,5,2)=20.0;   ksat(2,5,3)=20.0;    ksat(2,5,4)=20.0;   ksat(2,5,5)=20.0   !
!
ksat(2,1,6)=198.5;  ksat(2,1,7)=198.5;  ksat(2,1,8)=198.5;   ksat(2,1,9)=198.5;  ksat(2,1,10)=198.5   !
ksat(2,2,6)=198.5;  ksat(2,2,7)=198.5;  ksat(2,2,8)=198.5;   ksat(2,2,9)=198.5;  ksat(2,2,10)=198.5   !
ksat(2,3,6)=198.5;  ksat(2,3,7)=198.5;  ksat(2,3,8)=198.5;   ksat(2,3,9)=198.5;  ksat(2,3,10)=198.5   !
ksat(2,4,6)=20.0;   ksat(2,4,7)=20.0;   ksat(2,4,8)=20.0;    ksat(2,4,9)=20.0;   ksat(2,4,10)=20.0   !
ksat(2,5,6)=20.0;   ksat(2,5,7)=20.0;   ksat(2,5,8)=20.0;    ksat(2,5,9)=20.0;   ksat(2,5,10)=20.0   !

! species 3
ksat(3,1,1)=5.0;     ksat(3,1,2)=5.0;     ksat(3,1,3)=5.0;    ksat(3,1,4)=5.0;    ksat(3,1,5)=5.0   !
ksat(3,2,1)=145.0;   ksat(3,2,2)=147.0;   ksat(3,2,3)=149.0;  ksat(3,2,4)=151.0;  ksat(3,2,5)=153.0   !
ksat(3,3,1)=145.0;   ksat(3,3,2)=147.0;   ksat(3,3,3)=149.0;  ksat(3,3,4)=151.0;  ksat(3,3,5)=153.0   !
ksat(3,4,1)=145.0;   ksat(3,4,2)=147.0;   ksat(3,4,3)=149.0;  ksat(3,4,4)=151.0;  ksat(3,4,5)=153.0   !
ksat(3,5,1)=145.0;   ksat(3,5,2)=147.0;   ksat(3,5,3)=149.0;  ksat(3,5,4)=151.0;  ksat(3,5,5)=153.0   !

ksat(3,1,6)=5.0;     ksat(3,1,7)=5.0;     ksat(3,1,8)=5.0;    ksat(3,1,9)=5.0;    ksat(3,1,10)=5.0   !
ksat(3,2,6)=154.0;   ksat(3,2,7)=155.0;   ksat(3,2,8)=156.0;  ksat(3,2,9)=157.0;  ksat(3,2,10)=158.0   !
ksat(3,3,6)=154.0;   ksat(3,3,7)=155.0;   ksat(3,3,8)=156.0;  ksat(3,3,9)=157.0;  ksat(3,3,10)=158.0   !
ksat(3,4,6)=154.0;   ksat(3,4,7)=155.0;   ksat(3,4,8)=156.0;  ksat(3,4,9)=157.0;  ksat(3,4,10)=158.0   !
ksat(3,5,6)=154.0;   ksat(3,5,7)=155.0;   ksat(3,5,8)=156.0;  ksat(3,5,9)=157.0;  ksat(3,5,10)=158.0   !

! species 4
ksat(4,1,1)=5.0;     ksat(4,1,2)=5.0;    ksat(4,1,3)=5.0;    ksat(4,1,4)=5.0;    ksat(4,1,5)=5.0   !
ksat(4,2,1)=5.0;     ksat(4,2,2)=5.0;    ksat(4,2,3)=5.0;    ksat(4,2,4)=5.0;    ksat(4,2,5)=5.0   !
ksat(4,3,1)=5.0;     ksat(4,3,2)=5.0;    ksat(4,3,3)=5.0;    ksat(4,3,4)=5.0;    ksat(4,3,5)=5.0   !
ksat(4,4,1)=100.0;   ksat(4,4,2)=100.0;  ksat(4,4,3)=100.0;  ksat(4,4,4)=100.0;  ksat(4,4,5)=100.0   !
ksat(4,5,1)=100.0;   ksat(4,5,2)=100.0;  ksat(4,5,3)=100.0;  ksat(4,5,4)=100.0;  ksat(4,5,5)=100.0   !

ksat(4,1,6)=5.0;    ksat(4,1,7)=5.0;    ksat(4,1,8)=5.0;    ksat(4,1,9)=5.0;    ksat(4,1,10)=5.0   !
ksat(4,2,6)=5.0;    ksat(4,2,7)=5.0;    ksat(4,2,8)=5.0;    ksat(4,2,9)=5.0;    ksat(4,2,10)=5.0   !
ksat(4,3,6)=5.0;    ksat(4,3,7)=5.0;    ksat(4,3,8)=5.0;    ksat(4,3,9)=5.0;    ksat(4,3,10)=5.0   !
ksat(4,4,6)=100.0;  ksat(4,4,7)=100.0;  ksat(4,4,8)=100.0;  ksat(4,4,9)=100.0;  ksat(4,4,10)=100.0   !
ksat(4,5,6)=100.0;  ksat(4,5,7)=100.0;  ksat(4,5,8)=100.0;  ksat(4,5,9)=100.0;  ksat(4,5,10)=100.0   !

! species 6
ksat(6,1,1)=5.0;    ksat(6,1,2)=5.0;    ksat(6,1,3)=5.0;    ksat(6,1,4)=5.0;    ksat(6,1,5)=5.0   !
ksat(6,2,1)=5.0;    ksat(6,2,2)=5.0;    ksat(6,2,3)=5.0;    ksat(6,2,4)=5.0;    ksat(6,2,5)=5.0   !
ksat(6,3,1)=20.0;   ksat(6,3,2)=21.0;   ksat(6,3,3)=22.0;   ksat(6,3,4)=23.0;   ksat(6,3,5)=24.0   !
ksat(6,4,1)=101.0;  ksat(6,4,2)=102.0;  ksat(6,4,3)=103.0;  ksat(6,4,4)=104.0;  ksat(6,4,5)=105.0   !
ksat(6,5,1)=110.0;  ksat(6,5,2)=120.0;  ksat(6,5,3)=130.0;  ksat(6,5,4)=140.0;  ksat(6,5,5)=150.0   !

ksat(6,1,6)=5.0;    ksat(6,1,7)=5.0;    ksat(6,1,8)=5.0;    ksat(6,1,9)=5.0;    ksat(6,1,10)=5.0   !
ksat(6,2,6)=5.0;    ksat(6,2,7)=5.0;    ksat(6,2,8)=5.0;    ksat(6,2,9)=5.0;    ksat(6,2,10)=5.0   !
ksat(6,3,6)=25.0;   ksat(6,3,7)=26.0;   ksat(6,3,8)=27.0;   ksat(6,3,9)=28.0;   ksat(6,3,10)=29.0   !
ksat(6,4,6)=106.0;  ksat(6,4,7)=107.0;  ksat(6,4,8)=108.0;  ksat(6,4,9)=109.0;  ksat(6,4,10)=110.0   !
ksat(6,5,6)=160.0;  ksat(6,5,7)=170.0;  ksat(6,5,8)=180.0;  ksat(6,5,9)=190.0;  ksat(6,5,10)=200.0   !


nmort=0.09               ! minimum instantaneous annual natural mortality rate experienced (set at background noise). Natural mortality based on SEDAR 52 for RS: KM 2-1-22
nmorth=nmort/365.0/24.0  ! convert to per hour instantenous rate
fmort=0.1                ! annual fishing rate per year instantenous          - could be made age dependent
fmorth=fmort/365.0/12.0  ! fishing is only imposed during 12 hours of daylight - convert to hourly
pmort=0.5                ! maximum possible additional predation mortality (jacks and blues) that can occur
pmorth=pmort/365.0/12.0  ! made hourly

daydist=200.0          ! snapper move 200 m per hour when heading back to rig starting at sunset or when inight switched to 1
nightdist=200.0        ! snapper move 200  meters per hour when foraging during daytime hours (when inight switched to 2)
                       ! typical weight at age of red snapper in g ww at the end of the age - age-1 means 12 months

!NEW - Makes age class intiwt (therefore recruitment weight), species specific
!Red Snapper              Pinfish                   Atlantic Croaker          Bluefish
initwt(1,1)=111.5;        initwt(2,1)=36.87;        initwt(3,1)=88.71;        initwt(4,1)=100.7
initwt(1,2)=416.6;        initwt(2,2)=58.43;        initwt(3,2)=181.86;       initwt(4,2)=450.6
Initwt(1,3)=961.2;        initwt(2,3)=78.85;        initwt(3,3)=286.30;       initwt(4,3)=800.0
initwt(1,4)=1738.8;       initwt(2,4)=96.61;        initwt(3,4)=388.82;       initwt(4,4)=1150.2
initwt(1,5)=2717.5;       initwt(2,5)=111.25;       initwt(3,5)=481.81;       initwt(4,5)=1500.0
initwt(1,6)=3854.2;       initwt(2,6)=122.91;       initwt(3,6)=561.98;       initwt(4,6)=1850.7
initwt(1,7)=5102.7;       initwt(2,7)=131.97;       initwt(3,7)=628.76;       initwt(4,7)=2200.7
initwt(1,8)=6420.1;       initwt(2,8)=138.91;       initwt(3,8)=683.06;       initwt(4,8)=2550.0
initwt(1,9)=7768.7;       initwt(2,9)=142.91;       initwt(3,9)=726.46;       initwt(4,9)=2900.0
initwt(1,10)=12121.2;     initwt(2,10)=146.91;      initwt(3,10)=760.73;      initwt(4,10)=3000.0

!Jack species             Greater Amberjack - Manooch & Potts 1997
initwt(5,1)=100.7;        initwt(6,1)=1448.1
initwt(5,2)=450.6;        initwt(6,2)=3595.8
initwt(5,3)=800.0;        initwt(6,3)=6229.4
initwt(5,4)=1150.2;       initwt(6,4)=8993.4
initwt(5,5)=1500.0;       initwt(6,5)=11653.9
initwt(5,6)=1850.7;       initwt(6,6)=14080.3
initwt(5,7)=2200.7;       initwt(6,7)=16214.8
initwt(5,8)=2550.0;       initwt(6,8)=18045.9
initwt(5,9)=2900.0;       initwt(6,9)=19588.3
initwt(5,10)=3000.0;      initwt(6,10)=20870.3

END subroutine

!----------------------------------------------------*
subroutine setupgrid
USE allvar

CALL setrig

! determine column and row and distance of nearst rig cell to each cell on the grid
colnearrig=0; rownearrig=0   ! start with cell numbers of nearest cell at zeroes
do i=1,ncol
   do j=1,nrow
       olddistrig=nrow*ncol*cellsize        ! use maximum distance possible for initial value prior to search
       do k=1,nrigs
           xside=ABS(rigcol(k)-i)*cellsize
           yside=ABS(rigrow(k)-j)*cellsize
           distrig=SQRT(xside**2 + yside**2)
           
       ! added this in to check to see if the model was actually using the other rigs but it's not - KM 5-1
       !    WRITE(40,1059)iyear,iday,ihour,inight,k,rigcol(k), rigrow(k),colnearrig(i,j),rownearrig(i,j),distrig,&
       !                & distnearrig(i,j),olddistrig
       !  1059   FORMAT(1x,9(i6,1x),3(e12.2,1x)) !replicate once, 11 values, i6 - create a column where it'll take integers of 6 places, then create 3 that are out to 12 decimals. so basically put you are writing out the file and formatting it with the 1060 format from above
   
           IF(distrig.lt.olddistrig)then
               distnearrig(i,j)=distrig           ! distance in meters from nearest rig cell to cell i,j
               colnearrig(i,j)=rigcol(k)      ! column number of nearest rig cell to cell i,j
               rownearrig(i,j)=rigrow(k)      ! row number of nearest rig cell to cell i,j
               olddistrig=distrig
           endif
       END do
  end do
end do

! determine columnc and row and distance of nearest pyramid cell to each cell on the grid
colnearpyr=0; rownearpyr=0   ! start with cell numbers of nearest cell at zeroes
do i=1,ncol
   do j=1,nrow
       olddistpyr=nrow*ncol*cellsize        ! use maximum distance possible for initial value prior to search
       do k=1,npyr
           xside=ABS(pyrcol(k)-i)*cellsize
           yside=ABS(pyrrow(k)-j)*cellsize
           distpyr=SQRT(xside**2 + yside**2)
           IF(distpyr.lt.olddistpyr)then
               distnearpyr(i,j)=distpyr           ! distance in meters from nearest rig cell to cell i,j
               colnearpyr(i,j)=pyrcol(k)      ! column number of nearest rig cell to cell i,j
               rownearpyr(i,j)=pyrrow(k)      ! row number of nearest rig cell to cell i,j
               olddistpyr=distpyr
           endif
       END do
  end do
end do

! determine columnc and row and distance of nearest nat cell to each cell on the grid
colnearnat=0; rownearnat=0   ! start with cell numbers of nearest cell at zeroes
do i=1,ncol
   do j=1,nrow
       olddistnat=nrow*ncol*cellsize        ! use maximum distance possible for initial value prior to search
       do k=1,numnatcells
           xside=ABS(natcol(k)-i)*cellsize
           yside=ABS(natrow(k)-j)*cellsize
           distnat=SQRT(xside**2 + yside**2)
           IF(distnat.lt.olddistnat)then
               distnearnat(i,j)=distnat           ! distance in meters from nearest rig cell to cell i,j
               colnearnat(i,j)=natcol(k)      ! column number of nearest rig cell to cell i,j
               rownearnat(i,j)=natrow(k)      ! row number of nearest rig cell to cell i,j
               olddistnat=distnat
           endif
       END do
  end do
end do

! here is the code for making prey r and K values cell specific based on distance to nearest rig cell
do i=1,ncol
   do j=1,nrow
      do k=1,nprey
         kprey(i,j,k)=kkprey(k)
         rprey(i,j,k)=rrprey(k)
         IF(distnearrig(i,j)/cellsize.lt.4.0.or.distnearpyr(i,j)/cellsize.lt.4.0.or.distnearnat(i,j)/cellsize.lt.4.0)then   ! option to set prey K and r to cell-specific values - eg near rig
           kprey(i,j,k)=kkprey(k)*1.0
           rprey(i,j,k)=rrprey(k)*1.0
         endif
      end do
   end do
end do

zprey=kprey   ! initial densities of all 5 prey types set to their carrying capacity (g ww/m2)

do i=1,ncol   ! write out the column and row and distance to nearest rig cell for every cell on the grid
  do j=1,nrow
     WRITE(9,1002)i,j,distnearrig(i,j),distnearpyr(i,j),distnearnat(i,j)
1002  FORMAT(1x,2(i6,1x),3(f12.5,1x))
  end do
end do


end subroutine
!------------------------------------------------------*
REAL*8 FUNCTION distmorth(ii,jj)
USE allvar
INTEGER ii,jj
REAL*8 tt1,m1,dmortc
! function that adjust mortality rate for distance of cell column ii and row jj to its nearest rig cell

xsp=xspecies(ifish)

dmortc=0.09/365/24  !(constrains mortality to some value for the distance function) hourly rate!!!

IF (distnearrig(ii,jj).le.distnearpyr(ii,jj).and.distnearrig(ii,jj).lt.distnearnat(ii,jj)) then ! If the fish is closer to a rig, pyramid, or nat, then apply that distance - KM 3-11
tt1=1.114 - 1.2003e-7*distnearrig(ii,jj)**2
else if (distnearpyr(ii,jj).lt.distnearrig(ii,jj).and.distnearpyr(ii,jj).le.distnearnat(ii,jj))then
tt1=1.114 - 1.2003e-7*distnearpyr(ii,jj)**2
else if (distnearnat(ii,jj).le.distnearrig(ii,jj).and.distnearnat(ii,jj).lt.distnearpyr(ii,jj))then
tt1=1.114 - 1.2003e-7*distnearnat(ii,jj)**2
end if

IF(distnearrig(ii,jj).gt.1000.0.or.distnearpyr(ii,jj).gt.1000.0.or.distnearnat(ii,jj).gt.1000.0) tt1=1.0
IF(tt1.lt.1.0)tt1=1.0

m1=dmortc  !red snapper
!scale natural mortality up for atl croaker and pinfish and greater amberjack
if (xsp.eq.2) then  !Nelson 2002 listed a yearly mortality rate of 0.78/yr
     m1=dmortc*6.0
end if
if (xsp.eq.3) then
     m1=dmortc*6.0
end if
if (xsp.eq.6) then
     m1=dmortc*2.8 ! SEDAR 70 estimated a natural mortality of 0.28/y average across ages for greater amberjack - KM 12/23
end if

distmorth=m1*tt1

!distmorth=0.0     ! override

return
end
!-----------------------------------------------------*

subroutine nearrig
USE allvar
! determine nearest rig cell for each fish for each hour

  colnearrig=0  ! start at zero
  rownearrig=0

  olddistrig=nrow*ncol*cellsize  ! largest distance that is possible to start with
       do k=1,nrigs  ! search over the list of rig cells
           xside=(rigcol(k))*cellsize + cellsize/2.0  ! assume rigs at center of cell - meters from left edge of grid
           yside=(rigrow(k))*cellsize + cellsize/2.0  ! meters from bottom of grid
           opp=(yside-xyloc(ifish))  ! meters of vertical line from current location to cetner of k-th rig cell
           adj=(xside-xxloc(ifish))  ! meters of horizontal line from current location to center of k-th rig cell
           distrig=SQRT(opp**2 + adj**2)*rigweight   ! WEIGHTED distance in m to the nearest rig cell
       !if (ifish.eq.23.and.iyear.eq.1.and.ihour.eq.14) then
       !  WRITE(41,1060)iyear,iday,ihour,inight,ifish,xspecies(ifish),k,xycell(ifish),xxcell(ifish),xcolnearrig(ifish),xrownearrig(ifish),distrig,&
       !                & xdistnearrig(ifish),olddistrig
       !  1060   FORMAT(1x,11(i6,1x),3(e12.2,1x)) !replicate once, 11 values, i6 - create a column where it'll take integers of 6 places, then create 3 that are out to 12 decimals. so basically put you are writing out the file and formatting it with the 1060 format from above
       !end if
       IF(distrig.lt.olddistrig)then ! see if shortest so far in search through the list of rig cells
               xdistnearrig(ifish)=distrig  ! actual distance of shortest distance
               xcolnearrig(ifish)=rigcol(k)  ! column number of the rig cell
               xrownearrig(ifish)=rigrow(k)   ! row number of the rig cell
               olddistrig=distrig  ! replace shortest so far with the new shortest vaale
           endif




       END do

end subroutine
!-----------------------------------------------------*

subroutine nearpyr
USE allvar
! determine nearest pyramid cell for each fish for each hour

  colnearpyr=0  ! start at zero
  rownearpyr=0

  olddistpyr=nrow*ncol*cellsize  ! largest distance that is possible to start with
       do k=1,npyr  ! search over the list of rig cells
           xside=(pyrcol(k))*cellsize + cellsize/2.0  ! assume rigs at center of cell - meters from left edge of grid
           yside=(pyrrow(k))*cellsize + cellsize/2.0  ! meters from bottom of grid
           opp=(yside-xyloc(ifish))  ! meters of vertical line from current location to cetner of k-th rig cell
           adj=(xside-xxloc(ifish))  ! meters of horizontal line from current location to center of k-th rig cell
           distpyr=SQRT(opp**2 + adj**2)*pyrweight(k)  ! WEIGHTED distance in m to the nearest kth pyr cell
            IF(distpyr.lt.olddistpyr)then ! see if shortest so far in search through the list of rig cells
               xdistnearpyr(ifish)=distpyr  ! actual distance of shortest distance
               xcolnearpyr(ifish)=pyrcol(k)  ! column number of the rig cell
               xrownearpyr(ifish)=pyrrow(k)   ! row number of the rig cell
               olddistpyr=distpyr  ! replace shortest so far with the new shortest vaale
           endif

!       if (ifish.eq.23.and.iyear.eq.1.and.iday.ge.300) then
!         WRITE(41,1060)iyear,iday,ihour,inight,ifish,xsp,k,xycell(ifish),xxcell(ifish),xcolnear(ifish),xrownear(ifish),dist,&
!                       & xdistnear(ifish),oldist
!         1060   FORMAT(1x,11(i6,1x),3(e12.2,1x)) !replicate once, 11 values, i6 - create a column where it'll take integers of 6 places, then create 3 that are out to 12 decimals. so basically put you are writing out the file and formatting it with the 1060 format from above
!       end if


       END do

end subroutine
!-----------------------------------------------------*

subroutine nearnat
USE allvar
! determine nearest rig cell for each fish for each hour

  colnearnat=0  ! start at zero
  rownearnat=0

  olddistnat=nrow*ncol*cellsize  ! largest distance that is possible to start with
       do k=1,numnatcells  ! search over the list of rig cells
           xside=(natcol(k))*cellsize + cellsize/2.0  ! assume rigs at center of cell - meters from left edge of grid
           yside=(natrow(k))*cellsize + cellsize/2.0  ! meters from bottom of grid
           opp=(yside-xyloc(ifish))  ! meters of vertical line from current location to cetner of k-th rig cell
           adj=(xside-xxloc(ifish))  ! meters of horizontal line from current location to center of k-th rig cell
           distnat=SQRT(opp**2 + adj**2)*natweight  ! WEIGHTED distance in m to the nearest nat cell
            IF(distnat.lt.olddistnat)then ! see if shortest so far in search through the list of rig cells
               xdistnearnat(ifish)=distnat  ! actual distance of shortest distance
               xcolnearnat(ifish)=natcol(k)  ! column number of the rig cell
               xrownearnat(ifish)=natrow(k)   ! row number of the rig cell
               olddistnat=distnat  ! replace shortest so far with the new shortest vaale
           endif

!       if (ifish.eq.23.and.iyear.eq.1.and.iday.ge.300) then
!         WRITE(41,1060)iyear,iday,ihour,inight,ifish,xsp,k,xycell(ifish),xxcell(ifish),xcolnear(ifish),xrownear(ifish),dist,&
!                       & xdistnear(ifish),oldist
!         1060   FORMAT(1x,11(i6,1x),3(e12.2,1x)) !replicate once, 11 values, i6 - create a column where it'll take integers of 6 places, then create 3 that are out to 12 decimals. so basically put you are writing out the file and formatting it with the 1060 format from above
!       end if

       END do

end subroutine

!----------------------------------------------------

!*FIX - inserted code to make call to onrig flexible.  Initialization and recruitment was problematic  MDC 2apr2009

subroutine onrig(kk)
USE allvar
! determine if individual is on a rig cell this hour
     xhab(kk)=1  ! assume off rig so can call from movenight and moveday
       do k=1,nrigs   ! search list of rigs
          IF(xxcell(kk).eq.rigcol(k).and.xycell(kk).eq.rigrow(k))then   ! column and row correspond to a rig cell
             xhab(kk)=0   ! switch value of flag to denote now on a rig cell
! *NEW - Inserted IF THEN statement below.  Otherwise no original code was modified or deleted.
! In the original code when the fish is on a rig in the daytime, the fish is placed at center
! During the evening the fish does not get relocated to the center of the rig cell.  MDC 17March09

             IF(inight.eq.1)then  !*NEW - see note above. ! changes to inight..eq.1 based on the KM changes in night/day movement
             xxloc(kk)=(rigcol(k))*cellsize + cellsize/2.0   ! put fish who reach a rig during nighttime at center of rig cell
             xyloc(kk)=(rigrow(k))*cellsize + cellsize/2.0
             END if
          endif                                                  
       END do

if (inight.eq.1.and.xhab(kk).eq.0) then 
    if (xxcell(kk).eq.xxlastrig(kk).and.xycell(kk).eq.xylastrig(kk)) then
      xuniqrig(kk)=xuniqrig(kk)
  !    print *, 'no new rig', xuniqrig(kk)
  !    pause
    else
      xuniqrig(kk)=xuniqrig(kk)+1
   !   print *, 'new rig added', xuniqrig(kk)
   !   pause
    end if

    xxlastrig(kk)=xxcell(kk)
    xylastrig(kk)=xycell(kk)
end if

end subroutine onrig

!----------------------------------------------------
! didn't include last nat and unique nat because we don't care about that. just want to make sure that they know where the nat is and that they are generally staying there
subroutine onnat(kk)
USE allvar
! determine if individual is on a rig cell this hour
       do k=1,numnatcells   ! search list of rigs
          IF(xxcell(kk).eq.natcol(k).and.xycell(kk).eq.natrow(k))then   ! column and row correspond to a nat cell
             xhab(kk)=3   ! switch value of flag to denote now on a nat cell
             IF(inight.eq.1)then  !move fish to the center of the cell at night - KM 2/22
             xxloc(kk)=(natcol(k))*cellsize + cellsize/2.0   ! put fish who reach a nat during daytime at center of nat cell
             xyloc(kk)=(natrow(k))*cellsize + cellsize/2.0
             END IF
          END IF                                                  
       END do
if (inight.eq.1.and.xhab(kk).eq.3) then ! changed to 1 KM 2/22
    if (xxcell(kk).eq.xxlastnat(kk).and.xycell(kk).eq.xylastnat(kk)) then
      xuniqnat(kk)=xuniqnat(kk)
  !    print *, 'no new rig', xuniqrig(kk)
  !    pause
    else
      xuniqnat(kk)=xuniqnat(kk)+1
   !   print *, 'new rig added', xuniqrig(kk)
   !   pause
    end if

    xxlastnat(kk)=xxcell(kk)
    xylastnat(kk)=xycell(kk)
end if
end subroutine onnat

!----------------------------------------------------

!DEL - commented out subroutine onrig and substituted the onrig(kk) version.  So that xrig can be updated in the initial conditions and within the loop
!CP2 --- removed and placed at the bottom of the program with full remarks - !MDC - 2apr2009

!*FIX - inserted code to make call to onrig flexible.  Initialization and recruitment was problematic  MDC 2apr2009

subroutine onpyr(kk)
USE allvar
! determine if individual is on a rig cell this hour
       do k=1,npyr   ! search list of rigs
          IF(xxcell(kk).eq.pyrcol(k).and.xycell(kk).eq.pyrrow(k))then   ! column and row correspond to a pyr cell
             xhab(kk)=2   ! switch value of flag to denote now on a pyr cell
             IF(inight.eq.1)then  !*NEW - see note above. changed to 1 KM 2/22
             xxloc(kk)=(pyrcol(k))*cellsize + cellsize/2.0   ! put fish who reach a pyr during daytime at center of pyr cell
             xyloc(kk)=(pyrrow(k))*cellsize + cellsize/2.0
             END if
          endif                                                  
       END do

if (inight.eq.1.and.xhab(kk).eq.2) then ! changed to 1 KM 2/22
    if (xxcell(kk).eq.xxlastpyr(kk).and.xycell(kk).eq.xylastpyr(kk)) then
      xuniqpyr(kk)=xuniqpyr(kk)
    else
      xuniqpyr(kk)=xuniqpyr(kk)+1
    end if

    xxlastpyr(kk)=xxcell(kk)
    xylastpyr(kk)=xycell(kk)
end if

end subroutine onpyr

!---------------------------------------------------------------------
subroutine initfish(xsp)
USE allvar
REAL*8 num2,totpop1
REAL*8 ran1
INTEGER xsp
! set up initial populations of fish
initnum=0.0
inight=2   !*NEW - initialized inight = 2, so that call onrig, which sets initial xrig status, works correctly
           !       onrig requires a day/night status.  inight resets within the hour DO Loop, initialization will only happen once.  MDC 20March09. Changed to inight=2 for initialization because onrig status is now in moveday instead of movenight - KM 3-16

!DEL - commented out next line.  Total population depends on rig number.  Making tests of rig number impossible to test MDC 29apr09
!num2=nrigs*totpop(xsp)   ! initial age-1 individuals per rig times number of rigs

!INS - next line (removal of nrigs*).  Makes totpop a constant so that we can test the effect of adding rigs against a constant pop size
num2=totpop(xsp)   ! initial age-1 individuals per rig times number of rigs

initnum(xsp,1)=num2
rec(xsp)=num2   ! recruitment (number of age-1) by species

do i=2,inyears
   initnum(xsp,i)=initnum(xsp,i-1)*EXP(-4.0*nmort-fmort) ! initial number by age in stable age distribution 4 x natural as a guess for movement mort
END do

! totpop1 is the total initial number of individuals of each species - initial total pop - not used
totpop1=0.0
do i=1,inyears
  totpop1=totpop1+initnum(xsp,i)
end do

do i=1,inyears
 do j=1,imodf
    jj=(i-1)*imodf+j+sumfish    ! this jj continues over the three species so individuals of all species are in a single list
    xxcell(jj)=INT(ran1(idum)*ncol)+1    ! random column
    xycell(jj)=INT(ran1(idum)*nrow)+1     ! random row
! convert cell locations to continuous values of meters from origin - randomly inside the cell
    xxloc(jj) = (xxcell(jj)-1)*cellsize + ran1(idum)*cellsize
    xyloc(jj) = (xycell(jj)-1)*cellsize + ran1(idum)*cellsize
!    xweight(jj)=initwt(i)  ! initial weight from weight-at-age schedule - made it species specific MDC
    xweight(jj)=initwt(xsp,i)  ! initial weight from weight-at-age schedule
    xgrowth(jj)=0.1        ! does not matter as growth is called immediately

    xuniqrig(jj)=0
    xxlastrig(jj)=-1
    xylastrig(jj)=-1

    xuniqpyr(jj)=0
    xxlastpyr(jj)=-1
    xylastpyr(jj)=-1

    xuniqnat(jj)=0
    xxlastnat(jj)=-1
    xylastnat(jj)=-1

    !*DEL - Commented out xrig(jj) = 0.  Why start all fish arbitrarily at 0.  Why not find their exact status?  MDC 20March09
    !xrig(jj)=0             ! start with xrig equal zero - overwritten if actually started on a rig cell.  See *FIX below
    !*FIX - put in onrig(jj) so that onrig can be used in initfish, the loop and the call AGE  MDC 2apr2009
    call onrig(jj)
    call onpyr(jj)
    call onnat(jj)

    xage(jj)=i             ! age is i - 1,2,3,...
     IF(xage(jj).le.2)xstage(jj)=1  ! assign stage value based on age
     IF(xage(jj).ge.3.and.xage(jj).le.4)xstage(jj)=2
     IF(xage(jj).ge.5.and.xage(jj).le.6)xstage(jj)=3
     IF(xage(jj).ge.7.and.xage(jj).le.8)xstage(jj)=4
     IF(xage(jj).ge.9.and.xage(jj).le.10)xstage(jj)=5
    xworth(jj)=initnum(xsp,i)/float(imodf)   ! initial worth is pop number of age-1 divided by the number of model individuals
    xalive(jj)=0        ! alive flag 0 = alive
    xspecies(jj)=xsp   ! species number (1=snapper, 2=grunts, 3=grouper, 4 bluefish)
    xriskhr(jj)=0   ! start risk (number of hours of daylight not on a rig cell) at zero
   end do
END do

sumfish = jj  ! store running total of new individuals so can add to it for the next species --sumfish is used to
              ! start jj for the next species stopped in the loop, jj=(i-1)*imodf + sumfish

!  PRINT *,'sumfish i jj',sumfish,i,jj
!  pause

end subroutine initfish

!*******************************************************

subroutine growth(iflag,icol,jrow)
USE allvar
REAL*8 ran1
REAL*8 cnum
REAL*8 ka,kb,g1,g2,l1,l2
REAL*8 ck1,ck4,cto,cq,ctl,cm
INTEGER icol,jrow,iflag,xsp
!  called every hour during daytime hours (1 to 12)
         
    xsp=xspecies(ifish)         ! short hand species id
    is=xage(ifish)  ! shorthand for stage value of this individual

     if (xsp.eq.4) then
       !wisconsin consumption equation 3 for cold water species (bluefish)
       !hardwired for bluefish bioenergetics estimates

       tt5=1.0/(23.0 - 10.2)
       t5=tt5 * log((0.98*(1 - 0.156))/(0.156*(1.0 - 0.98)))
       t4=EXP(t5*(watertempbycol(icol,jrow) - 10.2))
       ka=(0.156*t4)/(1+0.156*(t4-1.0))

       tt7=1.0/(32.0 - 28.0)
       t7=tt7 * log((0.98*(1 - 0.850))/(0.850*(1.0 - 0.98)))
       t6=EXP(t7*(32.0 - watertempbycol(icol,jrow)))
       kb=(0.850*t6)/(1+0.850*(t6-1.0))

       !written in case we want to include more than 1 cold water species
       !tt5=1.0/(te2(xsp)-te1(xsp))
       !t5=tt5 * log((xk2(xsp)*(1-xk1(xsp)))/(xk1(xsp)*(1.0-xk2(xsp))))
       !t4=EXP(t5*(watertemp-te1(xsp)))
       !ka=(xk1(xsp)*t4)/(1+xk1(xsp)*(t4-1.0))

       !tt7=1.0/(te4(xsp)-te3(xsp))
       !t7=tt7 * log((xk3(xsp)*(1-xk4(xsp)))/(xk4(xsp)*(1.0-xk3(xsp))))
       !t6=EXP(t7*(te4(xsp)-watertemp))
       !kb=(xk4(xsp)*t6)/(1+xk4(xsp)*(t6-1.0))

       gctemp=ka*kb

!       ga(4)=0.5197  !bluefish a and b consumption coefficients hardwired rather than written in
!       gb(4)=-0.288  !hardwired because the fields for cold water species have more coefficients to deal with

       ga(4)=0.1732  !bluefish a and b consumption coefficients hardwired rather than written in
       gb(4)=-0.288  !hardwired because the fields for cold water species have more coefficients to deal with


     else if (xsp.eq.1) then
       ! creating a stepwise equation for optimal temperature with age ONLY for red snapper
       rsgtopt = -0.1111 * (is) + 27.5111
       ! wisconsin consumption equation 2 for warm water species (redsnapper, pinfish, croaker, amberjack) 
       v=(gtmax(xsp)-watertempbycol(icol,jrow))/(gtmax(xsp) - rsgtopt)
       z=(log(gtheta1(xsp)))*(gtmax(xsp) - rsgtopt)
       y=(log(gtheta1(xsp)))*(gtmax(xsp) - rsgtopt+2)
       x=((z**2)*((1+(sqrt((1+40)/y)))**2))*0.0025

       gctemp=(v**x)*EXP(x*(1-v))
       
     else 
       !wisconsin consumption equation 2 for warm water species (redsnapper, pinfish, croaker, amberjack) 
       v=(gtmax(xsp)-watertempbycol(icol,jrow))/(gtmax(xsp) - gtopt(xsp))
       z=(log(gtheta1(xsp)))*(gtmax(xsp) - gtopt(xsp))
       y=(log(gtheta1(xsp)))*(gtmax(xsp) - gtopt(xsp)+2)
       x=((z**2)*((1+(sqrt((1+40)/y)))**2))*0.0025

       gctemp=(v**x)*EXP(x*(1-v))
     end if


!    if (xsp.eq.2) then
!    PRINT *, 'ga, gb', ga(xsp),gb(xsp)
!    PRINT *, 'gtmax,gtopt,gtheta1', gtmax(xsp), gtopt(xsp),gtheta1(xsp)
!    PRINT *, 'comsumption gctemp and gctemp2', gctemp1, gctemp
!    pause
!    END if

    !*DEL - MDC 9june09, had to make the gcmax equation general so that species specific bioenergetics are incorporated
    gcmax = ga(xsp)*xweight(ifish)**(gb(xsp)) * gctemp

     !IF(is.eq.0)then    !ECHK - MDC 1apr09
     !PRINT *,'iday ihour ifish is xsp',iday,ihour,ifish,is,xsp
     !pause
     !endif

     cnum=0.0
     avgden=(ener(1)+ener(2)+ener(3)+ener(4)+ener(5))/5.0   ! simple average over all prey types if nothing is eaten
                                                            ! otherwise would an energy density of zero
      do jj=1,nprey                                        
        cnum = cnum + zprey(icol,jrow,jj)*vul(xsp,jj,is)/ksat(xsp,jj,is)  ! denominator term for functional response
      end do

      totcon(ifish)=0.0; geat=0.0; con=0.0; egest=0.0; exc=0.0; sda=0.0

      IF(iflag.eq.1.and.dailyc(ifish).gt.1.0*gcmax)go to 101   ! iflag=1 means really growing (not a call from movement) dailyc for movement?
      avgden=0.0    !NEW - inserted avgden=0.0.  Otherwise it was using avgden as calculated above.  So the running total below was out of whack.  Fish were growing massive within 6 days.
        do 100 jj=1,nprey  ! loop over prey types  - should this be done in random order?
          !con(jj) = gcmax * zhab(icol,jrow)* zprey(icol,jrow,jj)*vul(xsp,jj,is)/ksat(xsp,jj,is)/(1.0+cnum)  !this is where you make rig cells have no consumption!!!!
          con(jj) = gcmax * zprey(icol,jrow,jj)*vul(xsp,jj,is)/ksat(xsp,jj,is)/(1.0+cnum)  !this is where you make rig cells have no consumption!!!!
          geat(jj) = con(jj)   ! g of prey type jj eaten/g snapper in this hour (the most would be cmax in a single hour)
          totcon(ifish) = totcon(ifish) + con(jj)    ! g of allprey types eaten/g snapper in this hour
          avgden=avgden + con(jj)*ener(jj)   ! average energy density this hour based on actual diet

          !if (xsp.eq.2.and.iflag.eq.1) then
          !PRINT *,'iday ihour ifish jj con totcon avgden ener',iday,ihour,ifish,jj,con(jj),totcon(ifish),avgden,ener(jj)  !*ECHK - MDC 1apr09
          !pause
          !END if

          IF(iflag.eq.1.AND.((dailyc(ifish)+totcon(ifish)).gt.1.0*gcmax))go to 101 ! if reached Cmax skip - the remaining prey types
 100   continue

 !  *NOTE - Original coding is still used to sum avgden consumed (avgden equation above).  Original coding was not in the correct units. Should be cal/gprey.  MDC 20March09
 !  *NEW - Inserted IF THEN code below to divide avgden by totcon so that the denominator then becomes gprey and the units are cal/gprey.  MDC 20March09

  !    PRINT *,'totcon avgden',totcon(ifish),avgden  !ECHK - MDC 1apr09
  !    pause

      IF(totcon(ifish).gt.0.0)avgden=avgden/totcon(ifish)  ! *NEW - see note above

  !   PRINT *,'totcon ifish avgden en',totcon(ifish),avgden,en    !ECHK - MDC 1apr09
  !     pause

      !* NEW MDC - redsnapper uses a temperature dependent function to calculate egestion and excretion
      IF(xsp.eq.1.or.xsp.eq.6)THEN 
       p=totcon(ifish)/gcmax   ! p-value
       pe=fa(xsp)*watertempbycol(icol,jrow)**(fb(xsp)) *EXP(fg(xsp)*p)
       pff=0.0
       pf=((pe-0.1)/0.9) * (1-pff)+pff
       egest = pf*totcon(ifish)   ! egestion g prey/g fish this hour
       exc=ua(xsp)*watertempbycol(icol,jrow)**ub(xsp)*EXP(ug(xsp)*p)*(totcon(ifish)-egest) ! excretion g prey/g fish this hour
       sda = 0.0
      END IF

      !pinfish, and atl croaker use a constant proportion of total consumption to calculate egest and exc
      IF(xsp.eq.2.or.xsp.eq.3) THEN
       egest = fa(xsp)*totcon(ifish)
       exc= ua(xsp)*(totcon(ifish)-egest)
      END IF

101  CONTINUE   ! if reached cmax we skip to here and then should continue for respiration for each of the 12 hours of nighttime

 !      if (iflag.eq.1.and.iyear.eq.1.and.iday.eq.1.and.ihour.eq.6.and.MOD(ifish,100).eq.0) then
 !      PRINT *,'bioenergetics',xsp, ga,gb,gtmax,gtopt,gtheta1,ra,rb,rtmax,rtopt,rtheta2,fa,fb,fg,ua,ub,ug
 !      PRINT *,'gcmax', gcmax,totcon
 !      pause
 !      END if

       IF(iflag.eq.1.and.jday.eq.201)then      ! spawning on jday 201
         t1=1.0
         IF(xage(ifish).eq.1)t1=0.0   ! immature age
         IF(xage(ifish).eq.2)t1=0.5   ! 50% maturity at age-2
         IF(xsp.eq.1)egg=0.08*t1             ! *NOTE - eggs are 8% of body mass.  MDC 20March09
         IF(xsp.eq.2)egg=0.04*t1
         IF(xsp.eq.3)egg=0.08*t1
         IF(xsp.eq.6)egg=0.0325*t1  ! Based off of Harris et al. 2007 - KM 2-17-22
       else
         egg=0.0
       endif

       !   egg = 0.0 !ECHK - MDC 1apr09

       ! cal/g dw
       en(1)=4999.0   ! Shipley USA dissertation
       en(2)=5200.0 !4229.0   ! Thompson USA dissertation
       en(3)=5050.0 !4245.0   ! Craig croaker bioenergetics
       en(4)=7498.0
       en(5)=5500.0 ! was getting a divide by 0 so just put this in here to check - KM 3-11
       en(6)=5500.0 ! SEDAR70-RD-02-BOFFF all for meat and liver and around 5500 cal/g

       call respiration   ! hourly respiration for nighttime hours

       IF(iflag.eq.1.or.iflag.eq.2)then   ! compute hourly growth rate (g prey/g fish) for real and for call from movement
              ! dw/wdt in g ww pred/gww pred for this hour

    !       avgden=en     !ECHK - MDC 1apr09

    !       tt1= (totcon(ifish)-egest-exc-sda)*avgden/en
    !       tt2= (resp/24.0)*avgden/en
    !       tt3 = egg/12.0

    !       IF(jday.eq.201.and.iflag.eq.1.and.xage(ifish).ge.2)then
    !      PRINT *,'jday ihour ifish avgden en tt1 tt2 tt3 egg xweight'
    !      PRINT *,jday,ihour,ifish,avgden,en,tt1,tt2,tt3,egg,xweight(ifish)

    !      pause
    !      endif

       grow=(totcon(ifish)-egest-exc-sda)*avgden/en(xsp) - (resp/24.0)*avgden/en(xsp) - egg/12.0 ! *FIX

        !* NEW check values generated.
        if (ifish.eq.501.and.iflag.eq.1) then
         !p=totcon(ifish)/gcmax   ! p-value
         !m=egest+exc+sda
         !r=resp/24
         !PRINT *,'iyear iday ihour ifish',iyear,iday,ihour,ifish
         !Print *,'p',p,'totcon',totcon(ifish),'egest',egest,'exc',exc,'sda',sda,'mtn',m,'resp',r,'ed', en(xsp)
         !PRINT *, 'grow', grow
         !pause
        endif !for print check above
       endif !iflag.eq.1.or.iflag.eq.2

end subroutine growth

!-----------------------------------------------------*

subroutine respiration
USE allvar
INTEGER icol1,jrow1,xsp
REAL*8 goxcon
! called hourly

 xsp= xspecies(ifish)         ! short hand species id

    v=(rtmax(xsp)-wtemp(ifish))/(rtmax(xsp) - rtopt(xsp))
    z=(log(rtheta2(xsp)))*(rtmax(xsp) - rtopt(xsp))
    y=(log(rtheta2(xsp)))*(rtmax(xsp) - rtopt(xsp)+2)
    x=((z**2)*((1+(sqrt(1+(40/y))))**2))*0.0025

    grtemp=(v**x)*EXP(x*(1.0-v))

    !optional print statement
  !  if (xsp.eq.2) then
  !  PRINT *, 'ra, rb', ra(xsp),rb(xsp)
  !  PRINT *, 'rtmax,rtopt,rtheta2', rtmax(xsp), rtopt(xsp),rtheta2(xsp)
  !  PRINT *, 'respiration grtemp and grtemp2', grtemp1, grtemp
  !  PRINT *, 'act', act(xsp)
  !  pause
  !  END if

  if(xsp.eq.2.and.xage(ifish).ge.1)act(xsp)=1.0

  if (xsp.eq.1)goxcon=3.2387
  if (xsp.eq.2)goxcon=3.8558
  if (xsp.eq.3)goxcon=2.9511
  if (xsp.eq.4)goxcon=2.1593
  if (xsp.eq.6)goxcon=3.3750 ! 14140 J/gO2 from Brodie et al. 2016, then used cross multiplication with the value from Campbell equation 11 and the value for species 1 above (14140/x)=(13569/3.2387)

  resp = ra(xsp)*xweight(ifish)**(rb(xsp)) * grtemp * act(xsp) * goxcon !3.9686   ! g prey/g fish, x 2 is activity multipler, 0.7938 is g oxygen/g prey

    IF(inight.eq.1.and.xhab(ifish).eq.1)then  ! nighttime hours and not on a rig, pyr, or nat cell
       !CHANGE - original conditions  MDC 6may09
       IF(xspecies(ifish).eq.1.or.xspecies(ifish).eq.3)resp=resp*1.3  ! 30% increase in respiration for snapper and grouper
       IF(xspecies(ifish).eq.2.or.xspecies(ifish).eq.6)resp=resp*1.3   ! 30% increase for grunts and greater amberjack

       !NEW - use higher respiration penalty for grunts  MDC 6may09
       !IF(xspecies(ifish).eq.1.or.xspecies(ifish).eq.3)resp=resp*1.1  ! 30% increase in respiration for snapper and grouper
       !IF(xspecies(ifish).eq.2)resp=resp*1.6   ! 60% increase for grunts
    endif

end subroutine respiration

!-----------------------------------------------------*

subroutine moveday()
USE allvar
USE nrtype; USE nrutil; USE nr1
INTEGER, PARAMETER:: maxcells=100   ! maximum number of cells in a neighborhood
REAL*8 oldgrow, ran1, oldbio
INTEGER xxout(maxcells),xxout1(maxcells),rowk(maxcells),colk(maxcells),xsp,ia
REAL*8 xcenter, ycenter, opp, adj, theta, tt1,tt2,tt3
INTEGER xxc,xyc
REAL xx(maxcells)    ! note not double precision so matches the dummary array in the numerical recipes
REAL*8 eucdist     ! euclidean distance formula to calculate distance actually moved this hour
  
  xsp=xspecies(ifish)     ! shorthand species id
  ia=xage(ifish)          ! shorthand age value
  icol=xxcell(ifish)      ! shorthand column number
  jrow=xycell(ifish)      ! shorthand row number
  nhood=2                 ! plus/minus number of cells in x and y directions in searching neighborhood
  oldgrow=-99.0           ! dummy very small number that gets overwritten with first cell evaluated
  oldbio=-99.0
  xx=99.0   ! this is a trick - values will be 0 to 1 and so elements not needed in neighbohood will get 99 and be ranked last
  xxout=0; xxout1=0
  newcol=xxcell(ifish)    ! the column of the best cell in neighborhood - starts at current location
  newrow=xycell(ifish)    ! the row of best cell inthe neighborhood - starts at the curent location
  iii=1                   ! counter for number of viable cells in neighborhood
    do ii=-nhood,nhood   ! search nhood cells the left of current column and to right of current column
       ic=icol + ii      ! convert nhood to the actual column number on the grid
       iF(ic.le.0)ic = 1  ! truncate if went too far left
       IF(ic.ge.ncol)ic = ncol   ! truncate if went too far right
       do jj=-nhood,nhood   ! for each look at a column , look up and down the column (i.e., rows)
          jr=jrow + jj     ! convert to acutual row number based on current row number
          IF(jr.le.0)jr=1   ! truncate if too far down
          IF(jr.ge.nrow)jr = nrow  ! truncate if too far up
          rowk(iii)=jr     ! store row number in an array so later can know the row of the iii-th cell in the list
          colk(iii)=ic     ! store column number so can later can know the column number of the iii-th element in the list
          iii=iii+1        ! number of cells in actual neighborhood nhood=1 should be 9 cells; nhood=2 should be 26 cells
       end do              ! is this OK?  or should I not count cells that were off the grid instead of truncating their row and col?
    end do                 ! could do this with go to statements that skip rowk and colk and incrementing iii KAR


    do i=1,iii-1           ! put random numbers in an array whose length is the number of cells in neighborhood (i.e., iii-1 cells)
         xx(i)=ran1(idum)
    end do

!   call indexxold(iii-1,xx,xxout)  ! old version of indexx - did not work correctly but I am not sure why
    call indexx_sp(xx,xxout1)   ! new version from numerical receipes for fortran90 - xxout1 is integers
    xxout=rank(xxout1)          ! now xxout is the rank values of the real numbers in the array xx

    do 200 ii = 1,iii-1     ! search all cells in the neighborhood as a list rather than plus/minus number of rows and columns
       icell = xxout(ii)    ! icell is ii-th element in the list of ranks
       ic = colk(icell)     ! ic is the column number of the ii-th cell in the list
       jr = rowk(icell)       ! jr is the row number of the ii-th cells in the list

       !only call growth to move species 1,2,3
       if (xsp.le.3) then
        call growth(2,ic,jr)   ! the 2 tells growth subroutine that this call is not for real but for evaluating cells for movement
         IF(ihour.le.2)THEN !first 2 hours off the rig choose a random direction to swim
           newcol = ic     ! column number of best cell so far in the list
           newrow = jr     ! row number of best cells in the list so far
         else   ! all other hours try and maximize growth based on knowledge of your neighborhood
            ! Creating the stepwise decreases on the amount of prey difference that is needed to get the fish to move using the new variable onto_mov - KM 2/22 
           IF(xstage(ifish).eq.1)onto_mov = 1.20
           IF(xstage(ifish).eq.2)onto_mov = 1.15
           IF(xstage(ifish).eq.3)onto_mov = 1.10
           IF(xstage(ifish).eq.4)onto_mov = 1.05
           IF(xstage(ifish).eq.5)onto_mov = 1.02
           IF(grow.gt.onto_mov*oldgrow)THEN  ! note 1.2 to avoid 4.999 versus 5.0000 from attracting fish - must be at least 20% greater !!! CREATE NEW FUNCTION TO HAVE THIS DECREASE BY AGE - NEED NEW VARIABLE
             newcol = ic     ! column number of best cell so far in the list
             newrow = jr     ! row number of best cells in the list so far
             oldgrow=grow    ! actual growth rate projected for the best cell so far
           endif
         endif
       endif

       !NEW MDC, find highest prey density to move bluefish and jack around grid, called for both day and night movement
       !Removing the if statement comparing movebio to oldbio
       if (xsp.ge.4) then
        IF(movebio(ic,jr).gt.oldbio)THEN  ! note 1.2 to avoid 4.999 versus 5.0000 from attracting fish - must be at least 20% greater
            newcol = ic     ! column number of best cell so far in the list
            newrow = jr     ! row number of best cells in the list so far
            oldbio=movebio(ic,jr)    ! actual biomass for the best cell so far (previous hours biomass all species)
        endif  ! movebio
       end if ! xsp.ge.4

!       if (mod(ifish,50).eq.0.and.iyear.eq.10.and.iday.ge.300) then
!         WRITE(40,1059)iyear,iday,ihour,inight,ifish,xsp,icol,jrow,ii,icell,colk(icell),rowk(icell),newcol,newrow,grow,oldgrow
!         1059   FORMAT(1x,14(i6,1x),2(e12.6,1x))
!       end if

200  continue

!if (iyear.eq.5.and.xsp.eq.4.and.ifish.eq.1900) then
!   PRINT *, 'oldbio, x, y old',oldbio, xxcell(ifish),xycell(ifish), 'new x y',newcol,newrow
!   pause
!end if

!  determine center location of the desired cell - assume in the middle
     xcenter = (newcol-1)*cellsize + cellsize/2.0   ! this is the column location in meters from the lower left corner (0,0)
     ycenter = (newrow-1)*cellsize + cellsize/2.0   ! this is the row location in meters from the lower left corner
     opp=(ycenter-xyloc(ifish))       ! this is the length of the vertical line from the current location in meters to the center
     adj=(xcenter-xxloc(ifish))       ! this is the horizonatal line from the current location to the center in meters
     theta = atan(ABS(opp+0.001)/ABS(adj+0.001)) ! theta is the angle in radians

    IF(ran1(idum).lt.0.5)THEN   !  randomness on theta
        tt1=-1.0
    else
        tt1=1.0
    endif
    theta=theta + tt1*ran1(idum)*0.5   ! +- 0.5 out of 6 radians
     IF(opp.gt.0.0.and.adj.lt.0.0)theta=3.14-theta  ! theta is 0 to 2pie
     IF(opp.lt.0.0.and.adj.lt.0.0)theta=3.14+theta
     IF(opp.lt.0.0.and.adj.gt.0.0)theta=2.0*3.14-theta

! multiplier that is species specific - used below in this routine to create species differences in movement
     t1=1.0               ! snapper get a one so that inputted values are used for snapper
     IF(xsp.eq.2)t1=1.0   ! pinfish same as snapper
     IF(xsp.eq.3)t1=1.2   ! croaker wander more
     IF(xsp.eq.4)t1=1.5   ! 1.5 bluefish move fast
     IF(xsp.eq.5)t1=1.5   ! 1.5 jacks move fast
     IF(xsp.eq.6)t1=1.5   ! GA move fast

     IF(ran1(idum).lt.0.5)then !  randomness on daydist
        tt2=-1.0   ! move left
     else
        tt2=1.0     ! move right
     endif
    tt3 = daydist(ia) + tt2*daydist(ia)*0.3*ran1(idum)   ! +- 30% of inputted distance moved in an hour during the night (meters)
    IF(tt3.le.0.0)tt3=0.0   ! this is distance, nothing is never really negative. this is solved based on the if else statement from the previous 4 lines. can't have neg speed to they just stand still
    xdist = t1*tt3*COS(theta)   ! distance in meters moved in the x direction this hour
    ydist = t1*tt3*SIN(theta)

    eucdist = sqrt(xdist**2 + ydist**2)   ! distance from previous position (m)
    xdailymove(ifish)= xdailymove(ifish) + eucdist  ! running total of distance moved over a day (reset each day)

    xxl = xxloc(ifish) + xdist    ! temporary update of x location in meters from the lower left corner
    xyl = xyloc(ifish) + ydist    ! temporary update of the y location in meters

   IF(xxl.gt.-cellsize.and.xxl.lt.0.0)then  ! so new x location is negative but within a cellsize
          xdist=-xdist     ! THIS IS BAD CODING  KAR - switch sign so go the other direction - will use below for actual updating
   endif
   IF(xyl.GT.-cellsize.and.xyl.lt.0.0)ydist=-ydist     ! this is when you are off the edge but by less than one cell
                                                  ! no need to update xxc and xyc for switched xdist because truncation
                                                  ! when within one cell of the edge still puts in the edge cell
     xxc = INT(xxl/cellsize) + 1       ! new column number correcpodning to the new x location in meters
     xyc = INT(xyl/cellsize) + 1       ! new row number corresponding to the new y location in meters

    IF(xxc.le.0.or.xxc.GE.(ncol+1))then  ! now if more than cell off the grid inthe x direction (left or right) then switch sign
         xdist=-xdist
    endif
    IF(xyc.le.0.or.xyc.GE.(nrow+1))ydist=-ydist ! if more than a cell off the grid in the y direction (up or down) then switch sign

    xxloc(ifish) = xxloc(ifish) + xdist        ! really update positions because we are certain the new location is on the grid
    xyloc(ifish) = xyloc(ifish) + ydist

    xxcell(ifish) = INT(xxloc(ifish)/cellsize) + 1 ! now must update cell locations because xdist or ydist could have
    xycell(ifish) = INT(xyloc(ifish)/cellsize) + 1 ! changed enough to put you into a different cell
    
    ! distance to closest rig
    xcenterrig = (xxlastrig(ifish)-1)*cellsize + cellsize/2.0   ! this is the column location in meters from the lower left corner (0,0)
    ycenterrig = (xylastrig(ifish)-1)*cellsize + cellsize/2.0   ! this is the row location in meters from the lower left corner
    adjrig=(xcenterrig-xxloc(ifish))       ! this is the horizontal line from the current location to the center in meters
    opprig=(ycenterrig-xyloc(ifish))       ! this is the length of the vertical line from the current location in meters to the center
    xdistcentrig(ifish)=sqrt(adjrig**2 + opprig**2)  ! distance from center of previous nighttime rig, pyr, or nat cell (m)

    ! distance to closest pyramid
    xcenterpyr = (xxlastpyr(ifish)-1)*cellsize + cellsize/2.0   ! this is the column location in meters from the lower left corner (0,0)
    ycenterpyr = (xylastpyr(ifish)-1)*cellsize + cellsize/2.0   ! this is the row location in meters from the lower left corner
    adjpyr=(xcenterpyr-xxloc(ifish))       ! this is the horizontal line from the current location to the center in meters
    opppyr=(ycenterpyr-xyloc(ifish))       ! this is the length of the vertical line from the current location in meters to the center
    xdistcentpyr(ifish)=sqrt(adjpyr**2 + opppyr**2)  ! distance from center of previous nighttime rig, pyr, or nat cell (m)

    ! distance to closest nat cell
    xcenternat = (xxlastnat(ifish)-1)*cellsize + cellsize/2.0   ! this is the column location in meters from the lower left corner (0,0)
    ycenternat = (xylastnat(ifish)-1)*cellsize + cellsize/2.0   ! this is the row location in meters from the lower left corner
    adjnat=(xcenternat-xxloc(ifish))       ! this is the horizontal line from the current location to the center in meters
    oppnat=(ycenternat-xyloc(ifish))       ! this is the length of the vertical line from the current location in meters to the center
    xdistcentnat(ifish)=sqrt(adjnat**2 + oppnat**2)  ! distance from center of previous nighttime rig, pyr, or nat cell (m)

    if(jday.eq.182.and.xage(ifish).eq.1) THEN
       xdistcentrig(ifish)=0.0
       xdistcentpyr(ifish)=0.0
       xdistcentnat(ifish)=0.0
    endif

    call onrig(ifish)  ! check and see if fish is currently on a rig
    call onpyr(ifish)
    call onnat(ifish)

end subroutine moveday

!*****************************************************

subroutine movenight()
USE allvar
INTEGER xsp
REAL*8 ran1
INTEGER xxc,xyc
REAL*8 eucdist

! hourly move towards to nearest rig during the daylight hours (ihour is 13 to 24)
! see comments in movenight as many statements are the same as here with better comments

  xx=0.0; yy=0.0
  con=0.0; geat=0.0       ! zero these because do not call growth during nighttime movement
  xsp=xspecies(ifish)     ! shorthand for species id
  ia=xage(ifish)
  icol=xxcell(ifish)
  jrow=xycell(ifish)

  call nearrig   ! determine the distance in meters of a line from the current location to the center of the closest rig cell
  call nearpyr   ! determine the distance to the center of the closest pyramid cell
  call nearnat   ! determine the distance to the center of the closest nat cell

  ! can include weightings here to incorporate weighted habitat distances so that they can detect differences in pyramid density
  dist1=xdistnearrig(ifish)   ! shorthand distance in meters to nearest rig cell - ADDED DISTANCE WEIGHTINGS
  dist2=xdistnearpyr(ifish)   ! shorthand distance to nearest pyr cell
  dist3=xdistnearnat(ifish)   ! shorthand distance to nearest nat cell - ADDED DISTANCES WEIGHTINGS
  distrand = ran1(idum)
  
! have Matt check this KM 1-25-22. Making sure that all of my bases are covered if any of the fish happen to be perfectly between any two habitats. I am making the assumption that a fish will not be perfectly between all 3 hab types
  IF(dist1.lt.dist2.and.dist1.lt.dist3) THEN
  xcenter = (xcolnearrig(ifish)-1)*cellsize + cellsize/2.0 ! x distance in meters from lower left corner of center of nearest rig cell
  ycenter = (xrownearrig(ifish)-1)*cellsize + cellsize/2.0 ! y distance
  xhab(ifish) = 0
  ELSE IF(dist2.lt.dist1.and.dist2.lt.dist3) THEN
  xcenter = (xcolnearpyr(ifish)-1)*cellsize + cellsize/2.0 ! x distance in meters from lower left corner of center of nearest rig cell
  ycenter = (xrownearpyr(ifish)-1)*cellsize + cellsize/2.0 ! y distance
  xhab(ifish) = 2
  ELSE IF(dist3.lt.dist2.and.dist3.lt.dist1) THEN
  xcenter = (xcolnearnat(ifish)-1)*cellsize + cellsize/2.0 ! x distance in meters from lower left corner of center of nearest rig cell
  ycenter = (xrownearnat(ifish)-1)*cellsize + cellsize/2.0 ! y distance
  xhab(ifish) = 3
  ELSE IF((dist1.eq.dist2.or.dist1.eq.dist3).and.distrand.lt.0.5) THEN
  xcenter = (xcolnearrig(ifish)-1)*cellsize + cellsize/2.0 ! x distance in meters from lower left corner of center of nearest rig cell
  ycenter = (xrownearrig(ifish)-1)*cellsize + cellsize/2.0 ! y distance
  xhab(ifish) = 0
  ELSE IF((dist2.eq.dist1.or.dist2.eq.dist3).and.distrand.ge.0.5) THEN
  xcenter = (xcolnearpyr(ifish)-1)*cellsize + cellsize/2.0 ! x distance in meters from lower left corner of center of nearest rig cell
  ycenter = (xrownearpyr(ifish)-1)*cellsize + cellsize/2.0 ! y distance
  xhab(ifish) = 2
  ELSE IF(dist3.eq.dist2.and.distrand.lt.0.5) THEN
  xcenter = (xcolnearnat(ifish)-1)*cellsize + cellsize/2.0 ! x distance in meters from lower left corner of center of nearest rig cell
  ycenter = (xrownearnat(ifish)-1)*cellsize + cellsize/2.0 ! y distance
  xhab(ifish) = 3
  ELSE IF(dist3.eq.dist1.and.distrand.ge.0.5) THEN
  xcenter = (xcolnearnat(ifish)-1)*cellsize + cellsize/2.0 ! x distance in meters from lower left corner of center of nearest rig cell
  ycenter = (xrownearnat(ifish)-1)*cellsize + cellsize/2.0 ! y distance
  xhab(ifish) = 3
  ELSE
  PRINT *, "Have a fish that is equidistant between multiple habitat types. Need to add more scenarios in the IF statement in subroutine movenight()"
  STOP
  END IF
  
  opp=(ycenter-xyloc(ifish))  ! vertical line in meters from current location to center of nearest rig cell
  adj=(xcenter-xxloc(ifish))  ! horizontal line
  if (adj.eq.0) adj = 0.0000001 ! had to inclue b/c getting a floating point error "divide by zero" on the next line - KM 3/11

  theta = atan(ABS(opp)/ABS(adj))  ! angel in radians

  IF(ran1(idum).lt.0.5)then   ! randomness in theta
        tt1=-1.0      ! 50/50 to be negative or postive
  else
        tt1=1.0
  endif
  theta=theta  + tt1*ran1(idum)*0.25  ! +- 0.5 out of 2 pie (about 6) radians
  IF(opp.gt.0.0.and.adj.lt.0.0)theta=3.14-theta
  IF(opp.lt.0.0.and.adj.lt.0.0)theta=3.14+theta
  IF(opp.lt.0.0.and.adj.gt.0.0)theta=2.0*3.14-theta

  t1=1.0               ! species specific multiplier
  IF(xsp.eq.2)t1=0.8   ! 0.8 home slow
  IF(xsp.eq.3)t1=1.0   ! 1.0 home same
  IF(xsp.eq.6)t1=1.5   ! GA move fast - this should be same as moveday()

  IF(ran1(idum).lt.0.5)then
      tt2=-1.0
  else
      tt2=1.0
  endif

    tt3 = nightdist(ia) + tt2*nightdist(ia)*0.2*ran1(idum)   ! +- 20% of nightdistance
    IF(tt3.le.0.0)tt3=0.0
    xdist = t1*tt3*COS(theta)

    ydist = t1*tt3*sin(theta)

    eucdist = sqrt(xdist**2 + ydist**2)

    xdailymove(ifish)= xdailymove(ifish) + eucdist

    xxl = xxloc(ifish) + xdist   ! temporary new x location in meters from left edge
    xyl = xyloc(ifish) + ydist   ! temporary new y location in meters from bottom of grid

   IF(xxl.GT.-cellsize.and.xxl.lt.0.0)xdist=-xdist     ! THIS IS BAD CODING  KAR
   IF(xyl.GT.-cellsize.and.xyl.lt.0.0)ydist=-ydist     ! this is when you are off the edge but by less than one cell
                                                  ! no need to update xxc and xyc for switched xdist because truncation
                                                  ! when within one cell of the edge still puts in the edge cell
     xxc = INT(xxl/cellsize) + 1        ! wrapping around the column
     xyc = INT(xyl/cellsize) + 1        ! wrapping around the row 

    IF(xxc.le.0.or.xxc.GE.(ncol+1))xdist=-xdist    ! if you go too far them then move back towards the middle column
    IF(xyc.le.0.or.xyc.GE.(nrow+1))ydist=-ydist    ! bounce back towards the middle row

    xxloc(ifish) = xxloc(ifish) + xdist            ! distance in meters from the left wall
    xyloc(ifish) = xyloc(ifish) + ydist            ! distance in meters from the bottom row

    xxcell(ifish) = INT(xxloc(ifish)/cellsize) + 1 ! column update
    xycell(ifish) = INT(xyloc(ifish)/cellsize) + 1 ! row update

    call onrig(ifish)  ! update onrig status
    call onpyr(ifish)  ! update onpyr status
    call onnat(ifish)  ! update onnat status
end subroutine movenight

!**********************************************************

subroutine mortality
USE allvar
INTEGER xsp
REAL*8 m1,distmorth,m2,blueratio, jackratio, garatio, eqn8090,m3a,m3b,m3c,m3,m4,z
! could be day verus night dependent

   xsp=xspecies(ifish)

   m1=nmorth
   ! atlantic croaker and pinfish experience increased natural mortality
   if (xsp.eq.2) then  !Nelson 2002 listed a yearly mortality rate of 0.78/yr
     m1=nmorth*6.0
   end if
   if (xsp.eq.3) then
     m1=nmorth*6.0
   end if
   if (xsp.eq.6) then ! SEDAR 70
     m1=nmorth*3.0
   end if
   !extra predation mortality scaled by their current distance from nearest rig cell, unless located on rig cell
!   IF(inight.eq.2.and.xrig(ifish).eq.0)then   ! it is daylight and not on rig yet
   IF(xhab(ifish).eq.1)then   ! off rig, pyr, or nat
     m2=distmorth(xxcell(ifish),xycell(ifish)) !calculate additional natural mortality based on distance from closest rig
   else   ! on rig
!    m2=distmorth(xxcell(ifish),xycell(ifish)) !no refuge effect - or turn mortality on even though on a rig
     m2=0.0    !refuge effect - turn mortality off when on a rig
   endif
   
   !calculates the ratio of current cell predator density / gridwide avg. predator density
   if(iday.eq.1.and.ihour.eq.1)then  !hour 1 of day 1 biomass in previous hour is unknown so we can't scale mort.  Set to 0
      blueratio = 0.0
      jackratio = 0.0
      garatio = 0.0
   else
      blueratio = bio(4,xxcell(ifish),xycell(ifish))/bluebio ! cell by cell biomass of species 4 / gridwide average biomass species 4
      jackratio = bio(5,xxcell(ifish),xycell(ifish))/jackbio ! cell by cell biomass of species 5 / gridwide average biomass species 5
      garatio = bio(6,xxcell(ifish),xycell(ifish))/gabio     ! cell by cell biomass of species 6 / gridwide average biomass species 6
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !additional predation mortality scaled by  predator ratio above
!   !deathtrap scenario
!   m3a = pmorth*eqn8090(blueratio )    ! predation ratio multiplier (0-1) eqn8090, multiplied by pmorth (constrained)
!   m3b = pmorth*eqn8090(jackratio)    ! predation ratio multiplier (0-1) eqn8090, multiplied by pmorth (constrained)
!   m3 = m3a + m3b
!
!   !deathtrap scenario
!   if(inight.eq.1.and.xrig(ifish).eq.0)then
!     m3 = m3*0.0
!   endif
!   if(inight.eq.1.and.xrig(ifish).eq.1)then
!     m3 = m3*0.5
!   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !additional predation mortality scaled by  predator ratio above
   !refuge scenario

   if(inight.eq.1.and.xhab(ifish).eq.1)then  !nighttime/on benthic, apply mortality - KM
     m3a = pmorth*eqn8090(blueratio )    ! predation ratio multiplier (0-1) eqn8090, multiplied by pmorth (constrained)
     m3b = pmorth*eqn8090(jackratio)    ! predation ratio multiplier (0-1) eqn8090, multiplied by pmorth (constrained)
     m3c = pmorth*eqn8090(garatio)    ! predation ratio multiplier (0-1) eqn8090, multiplied by pmorth (constrained)
     m3 = m3a + m3b + m3c
   ! no refuge effect for younger fish at natural reef (not as much high relief nat reefs on the TX coast - KM 3-11
   ELSE IF(inight.eq.1.and.xsp.le.3.and.xage(ifish).le.4.and.xhab(ifish).eq.3)then
    m3a = pmorth*eqn8090(blueratio )    ! predation ratio multiplier (0-1) eqn8090, multiplied by pmorth (constrained)
     m3b = pmorth*eqn8090(jackratio)    ! predation ratio multiplier (0-1) eqn8090, multiplied by pmorth (constrained)
     m3c = pmorth*eqn8090(garatio)    ! predation ratio multiplier (0-1) eqn8090, multiplied by pmorth (constrained)
     m3 = m3a + m3b + m3c
   else  !nighttime/on rig, daytime both, provides refuge
     m3 = 0.0
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   if(m3.lt.0)m3=0.0  ! constrain the function to 0 ... otherwise you get increased survival
   ! Increasing fishing mortality for RS once they enter the directed fishery at age 2
   if(xhab(ifish).ne.1.and.xage(ifish).ge.2)then    !on habitat status to toggle refuge effect, species and age specific. Can change this to reflect Gallaway 2009 - might be age structure to fmort - KM
      m4=fmorth  !off rig fmorth applied
!   else if(xspecies(ifish).eq.1.and.xhab(ifish).ne.1.and.xage(ifish).le.1)then
!      m4=(0.57/365/light_hours)          !Based on the average fishing mortality estimates for age 0 and 1 RS reported in Gallaway et al. 2009 - KM 2-1-22
   else
      m4=0.0   ! full refuge - fmorth not applied on rig
   endif

   if(inight.eq.1)m4=m4*0.5  !night time fmorth cut in half - less fishing pressure at night (commercial fleet only)
   if(xsp.gt.1)m4=0.0  !no fishing mortality for pinfish and croaker and ga

   z = m1+m2+m3+m4

!   if(iday.eq.2.and.mod(ifish,100).eq.0)then
!     m1=m1*365*24; m2=m2*365*24; m3=m3*365*24; m4=m4*365*24
!     print *, 'species, night/day, onrig', xsp,inight,xrig(ifish)
!     print *, 'nmort,dmort,distance',m1,m2,distnear(xxcell(ifish),xycell(ifish))
!     print *, 'blues, jacks, predation', blueratio, jackratio,m3
!     print *, 'fishing mort', m4
!     print *, 'Z', m1+m2+m3+m4
!     pause
!     m1=m1/365/24; m2=m2/365/24; m3=m3/365/24; m4=m4/365/24
!   end if

   xcmort(ifish) = xcmort(ifish) + m1 + m2 + m3 + m4  ! total cumulative mortality in a day

   if(mod(iyear,10).eq.0.and.iday.lt.20.and.ifish.eq.100)then
     m1=m1*365*24; m2=m2*365*24; m3=m3*365*24; m4=m4*365*24; z=z*365*24

      WRITE(42,1061)cumhours,iyear,iday,ihour,inight,ifish,xsp,xhab(ifish),m1,distnearrig(xxcell(ifish),xycell(ifish)),&
      &             m2, blueratio, jackratio, garatio, m3, m4, z, xcmort(ifish)
      1061  FORMAT(1x,e12.5,7(i6,1x),9(f14.9,1x))

     m1=m1/365/24; m2=m2/365/24; m3=m3/365/24; m4=m4/365/24; z=z/365/24
   endif

   if(mod(iyear,10).eq.0.and.iday.le.20.and.ihour.eq.24)then
    xcmort(ifish) = xcmort(ifish)*365

      WRITE(43,1062)cumhours,iyear,iday,ihour,inight,ifish,xsp,xhab(ifish),xcmort(ifish)
      1062  FORMAT(1x,e12.5,7(i6,1x),f14.9)

    xcmort(ifish) = xcmort(ifish)/365
   end if

! partition natural and fishing mortality into yield (gww of each species, summed over individuals, in this hour)
   yield(xsp) = yield(xsp) + (m4/(m4+m1+m2+m3))*(1.-EXP(-m1-m2-m3-m4))*xworth(ifish)*xweight(ifish)
   yieldnum(xsp) = yieldnum(xsp) + (m4/(m4+m1+m2+m3))*(1.-EXP(-m1-m2-m3-m4))*xworth(ifish)

! partition natural and fishing mortality into natural losses
   naturalp(xsp) = naturalp(xsp) + ((m1+m2+m3)/(m1+m2+m3+m4))*(1.-EXP(-m1-m2-m3-m4))*xworth(ifish)*xweight(ifish)
   naturalpnum(xsp) = naturalpnum(xsp) + ((m1+m2+m3)/(m1+m2+m3+m4))*(1.-EXP(-m1-m2-m3-m4))*xworth(ifish)

! decrement worth based on total (natural plus fishing) mortality
   xworth(ifish)=xworth(ifish)*EXP(-m1-m2-m3-m4)
                   
   IF(xworth(ifish).lt.0.01)then   ! if worth get too small then remove the model individuals as it is worth so little
      xworth(ifish)=0.0
      xalive(ifish)=1
   endif

end subroutine mortality

!***********************************************************

!Predator ratio scaling function
!scales from 0-1 ... then used as a multiplier of maximum predation mortality you want

REAL*8 FUNCTION eqn8090(x)
!-----------------------------------------------------------
! TableCurve C:\Users\Kenny\Desktop\mortscale.f90 Jan 6, 2010 9:37:39 AM 
! C:\Users\Kenny\Desktop\mortscale.txt 
! X=  
! Y=  
! Eqn# 8090  AsymSig(a,b,c,d,e) 
! r2=0.9965032671209708D0 
! r2adj=0.979019602725825D0 
! StdErr=0.04181347198588728D0 
! Fval=142.4906193288693D0 
! a= -0.1426330958218728D0 
! b= 1.135386805451159D0 
! c= 17.64505768173239D0 
! d= 2.971982469222264D0 
! e= 0.2296318423765162D0 
! Constraints: d<>0,e>0 
!-----------------------------------------------------------
  REAL*8 x,y
  y=(-0.1426330958218728D0)+1.135386805451159D0/(1.0+&
   &DEXP(-(x-2.971982469222264D0*DLOG(2.0**(1.0/&
   &0.2296318423765162D0)-1.0)-&
   &17.64505768173239D0)/2.971982469222264D0))**&
   &0.2296318423765162D0 
  eqn8090=y
  RETURN
    END

!****************************************************************************
subroutine costben
USE allvar
REAL*8 weighting(nspecies,ncol,nrow)
REAL*8 comm_value(nspecies,ncol,nrow)
REAL*8 weighted_value(nspecies,ncol,nrow)
REAL*8 ecowt(ncol,nrow), totwt(ncol,nrow)
INTEGER xsp


rs_value = 12176300.
rs_landed = 2754861.
ga_value = 45695.
ga_landed = 22201.

! loop over all cells to find the cost benefit ratios for each cell
do i = 1,ncol
    do j = 1,nrow
 ! computing the weighted value for both RS and GA
            weighting(1,i,j) = rs_value/(rs_value+ga_value)
            comm_value(1,i,j) = rs_value/rs_landed
            weighted_value(1,i,j) = weighting(1,i,j)*comm_value(1,i,j)

            weighting(6,i,j) = ga_value/(ga_value+ga_value)
            comm_value(6,i,j) = ga_value/ga_landed
            weighted_value(6,i,j) = weighting(6,i,j)*comm_value(6,i,j)  

        ecowt(i,j) = (bio(1,i,j)*weighted_value(1,i,j)) + (bio(6,i,j)*weighted_value(6,i,j))  ! getting the total weighted biomass for each cell for RS and GA
        totwt(i,j) = bio(1,i,j)+bio(2,i,j)+bio(3,i,j)+bio(4,i,j)+bio(5,i,j)+bio(6,i,j)    ! getting the total overall biomass per cell

        ! costben output
        WRITE(101,2020) iyear,iday,jday,i,j,zhab(i,j),ecowt(i,j),totwt(i,j)
2020    FORMAT(1x,6(i5,1x),2(f12.4,1x))

    end do ! nrow
end do ! ncol   

end subroutine costben      
            
!****************************************************************************
subroutine preyupdate
USE allvar
REAL*8 ran1
! called hourly

!cumeat = 0 ! density depence switch

! loop over all cells
do i=1,ncol
  do j=1,nrow
     do k=1,nprey   ! loop over all prey types
        ! logistic update of each prey type with mort from fish eatening them
        zprey(i,j,k)=zprey(i,j,k) + zprey(i,j,k)*1.0*rprey(i,j,k)&
        &           *(1.0-zprey(i,j,k)/kprey(i,j,k)) -cumeat(i,j,k)/(cellsize**2)  ! per m2
        IF(zprey(i,j,k).le.0.001)zprey(i,j,k)=0.001   ! do not let them go to zero
     end do
     IF(iyear.gt.5.and.iday.eq.1)then
!     IF((iyear.EQ.1.or.MOD(iyear,10).eq.0).and.MOD(iday,6).eq.0.AND.(iday.gt.30.and.iday.lt.60))then
       !IF((iyear.EQ.1.or.MOD(iyear,5).eq.0).and.MOD(iday,30).eq.0) then ! switched to print out 5 yr info MDC
     !  IF(MOD(iday,30).eq.0.and.ihour.eq.1)then
       WRITE(11,1006)iyear,iday,jday,ihour,inight,i,j,(zprey(i,j,k),k=1,nprey)
1006   FORMAT(1x,7(i5,1x),10(f12.4,1x))
     ! endif
      endif
  end do
end do

end subroutine preyupdate

!*********************************************************

subroutine age
USE allvar
REAL*8 ran1
INTEGER xsp
! called once per year
! output here to get the last day of the last age class included

  idcount=0   ! number of newly freed up spots in the arrays due to removal of individuals from old age
  do i=1,itotf
     IF(xspecies(i).le.3.or.xspecies(i).eq.6)then ! species 4 does not age or recruit new (itotf = 1501-2000)
       xage(i)=xage(i)+1   ! one year older   1 is 12-24 months old  2=24 to 36
       IF(xage(i).le.2)xstage(i)=1  ! update stage flag for this individual
       IF(xage(i).ge.3.and.xage(i).le.4)xstage(i)=2
       IF(xage(i).ge.5.and.xage(i).le.6)xstage(i)=3
       IF(xage(i).ge.7.and.xage(i).le.8)xstage(i)=4
       IF(xage(i).ge.9.and.xage(i).le.10)xstage(i)=5

       IF(xage(i).eq.11)then    ! if reach 11 then die of old age and remove, freeing up the element for new fish
          idcount=idcount+1   ! one more spot available
          idaval(idcount)=i   ! actual loction inthe arrays that is available
       endif
     endif ! end xspecies(i).le.3
  end do

!  add new year-class for each species
itt=2   ! always add new individuals for snapper and GA
IF(ipinfish.eq.1)itt=itt+1   ! means 3 species in the model
IF(icroaker.eq.1)itt=itt+1  ! means 3 or 4 species in the model
! itt=2 red snapper and greater amberjack only; itt=3 snapper, ga, and either grunts or grouper; itt=4 all four species

  do j=1,imodf*itt   ! number of new individuals of all four species
      jj=idaval(j)   ! array element number that was occupied by an old age individual removed above
      idcount=idcount-1 ! one less in the counter of available spots
      IF(idcount.LE.-1)then
         PRINT *,'problem with idaval'
         pause
      endif

    xxcell(jj)=INT(ran1(idum)*ncol)+1    ! random initial placement in column
    xycell(jj)=INT(ran1(idum)*nrow)+1     ! random initial placement in row
    xxloc(jj) = (xxcell(jj)-1)*cellsize + ran1(idum)*cellsize  ! random meters location within cell in x direction
    xyloc(jj) = (xycell(jj)-1)*cellsize + ran1(idum)*cellsize  ! random meters location within cell in y direction
    xgrowth(jj)=0.1                ! arbitrary initial value for growth rate

    call onrig(jj)   ! update on rig status
    call onpyr(jj)
    call onnat(jj)

    xage(jj)=1       ! must be age-1
    xstage(jj)=1    ! since must be age-1 then must be stage 1
    xuniqrig(jj)=0
    xxlastrig(jj)=-1
    xylastrig(jj)=-1
    xuniqpyr(jj)=0
    xxlastpyr(jj)=-1
    xylastpyr(jj)=-1
    xuniqnat(jj)=0
    xxlastnat(jj)=-1
    xylastnat(jj)=-1

! assign species id (1=snapper, 2=grunts, 3=grouper, 6=greater amberjack)
    IF(j.le.imodf)xspecies(jj)=1  ! must always have snapper so first set of imodf new individuals is snapper

! assigning greater amberjack
    IF(j.gt.imodf.and.j.le.(2*imodf))xspecies(jj)=6 ! must have greater amberjack as 2 in case there aren't grunts or grouper

! so if second set of imodf number of new individuals and grunt are present then grunt
    IF((j.gt.(2*imodf).and.j.LE.(3*imodf)).and.ipinfish.EQ.1)xspecies(jj)=2

! second set of new ind but grouper and no grunts so these are grouper
    IF((j.gt.(2*imodf).and.j.LE.(3*imodf)).AND.(ipinfish.EQ.0.and.icroaker.eq.1))xspecies(jj)=3

! third set of imodf individuals - must have gruntsand grouper so these are grouper
    IF((j.gt.3*imodf).AND.(ipinfish.EQ.1.and.icroaker.eq.1))xspecies(jj)=3   ! fixed the imodf multiplier to read 2*imodf

    xweight(jj)=initwt(xspecies(jj),1)              ! initial weight by species for age 1, in g ww

    xworth(jj)=rec(xspecies(jj))/float(imodf)   ! start with constant recruitment (rec) every year!!!!!!!!
                                              ! so 1000 age-1 represented by 20 model ind means worth of each model ind is 50
    xalive(jj)=0   ! individual starts off alive
    xriskhr(jj)=0

  end do

end subroutine age

!********************************************************

subroutine dailyout
USE allvar
REAL*8 numage(nspecies,inyears),mwt(nspecies,inyears),mmwt(nspecies,inyears),zmean(nprey),totbio(nspecies),totp(nspecies)
REAL*8 mmxriskhr(nspecies,inyears),mxriskhr(nspecies,inyears)
REAL*8 mxdailymove(nspecies,inyears), mmxdailymove(nspecies,inyears)
REAL*8 mxuniqrig(nspecies,inyears), mmxuniqrig(nspecies,inyears)
REAL*8 mxuniqpyr(nspecies,inyears), mmxuniqpyr(nspecies,inyears)
REAL*8 mxuniqnat(nspecies,inyears), mmxuniqnat(nspecies,inyears)
REAL*8 mxdistcentrig(nspecies,inyears), mmxdistcentrig(nspecies,inyears)!, maxxdistcent(nspecies,inyears)
REAL*8 mxdistcentpyr(nspecies,inyears), mmxdistcentpyr(nspecies,inyears)!, maxxdistcent(nspecies,inyears)
REAL*8 mxdistcentnat(nspecies,inyears), mmxdistcentnat(nspecies,inyears)!, maxxdistcent(nspecies,inyears)
INTEGER xsp

numage=0.0; mwt=0.0; mmwt=0.0; mxriskhr=0.0;  mmxriskhr=0.0; mxdailymove=0.0; mmxdailymove=0.0
mxuniqrig=0.0; mmxuniqrig=0.0; mxuniqpyr=0.0; mmxuniqpyr=0.0; mxuniqnat=0.0; mmxuniqnat=0.0
mxdistcentrig=0.0; mmxdistcentrig=0.0; maxxdistcentrig=0.0
mxdistcentpyr=0.0; mmxdistcentpyr=0.0; maxxdistcentpyr=0.0
mxdistcentnat=0.0; mmxdistcentnat=0.0; maxxdistcentnat=0.0

do i=1,itotf
    xsp=xspecies(i)
    IF(xalive(i).eq.0)then
       ia=xage(i)   ! shorthand for age
       numage(xsp,ia)=numage(xsp,ia) + xworth(i)  ! numbers of pop individuals by species and age class
       mwt(xsp,ia) = mwt(xsp,ia) + xworth(i)*xweight(i)  ! biomass by species and age class
       mxriskhr(xsp,ia) = mxriskhr(xsp,ia) + xriskhr(i)*xworth(i)  ! hours at risk (off rig cells in daylight) by species and age
       mxdailymove(xsp,ia) = mxdailymove(xsp,ia) + xworth(i)*xdailymove(i)
       mxuniqrig(xsp,ia) = mxuniqrig(xsp,ia) + xworth(i)*xuniqrig(i)
       mxuniqpyr(xsp,ia) = mxuniqpyr(xsp,ia) + xworth(i)*xuniqpyr(i)
       mxuniqnat(xsp,ia) = mxuniqnat(xsp,ia) + xworth(i)*xuniqnat(i)
       mxdistcentrig(xsp,ia) = mxdistcentrig(xsp,ia) + xworth(i)*xdistcentrig(i)
       mxdistcentpyr(xsp,ia) = mxdistcentpyr(xsp,ia) + xworth(i)*xdistcentpyr(i)
       mxdistcentnat(xsp,ia) = mxdistcentnat(xsp,ia) + xworth(i)*xdistcentnat(i)
    endif
end do

!print *, mxuniqrig(1,1),mxuniqrig(2,2), mxuniqrig(3,3)
!pause

totbio=0.0
do i=1,inyears   ! loop over age classes to get averaged values from running sums computed above
  do j=1,nspecies   ! loop over species KAR - what if two or one species inthe run? - I think you see zeroes on output files
    IF(numage(j,i).gt.0.0)then
       mmwt(j,i)=mwt(j,i)/numage(j,i)   ! average weight of an individual over all ages by species
       mmxriskhr(j,i)=mxriskhr(j,i)/numage(j,i) ! average risk of an individual over all ages by species
       mmxdailymove(j,i)=mxdailymove(j,i)/numage(j,i)  ! average daily movement of an individual over all ages by species
       mmxuniqrig(j,i)=mxuniqrig(j,i)/numage(j,i)   ! average number of rig changes over a week (written out and reset after day 7)
       mmxuniqpyr(j,i)=mxuniqpyr(j,i)/numage(j,i)   ! average number of rig changes over a week (written out and reset after day 7)
       mmxuniqnat(j,i)=mxuniqnat(j,i)/numage(j,i)   ! average number of rig changes over a week (written out and reset after day 7)
       mmxdistcentrig(j,i)=mxdistcentrig(j,i)/numage(j,i)  ! average night distance moved away from rig by age and species
       mmxdistcentpyr(j,i)=mxdistcentpyr(j,i)/numage(j,i)  ! average night distance moved away from rig by age and species
       mmxdistcentnat(j,i)=mxdistcentnat(j,i)/numage(j,i)  ! average night distance moved away from rig by age and species
    else
       mmwt(j,i)=0.0
       mmxriskhr(j,i)=0
       mmxdailymove(j,i)=0.0
       mmxuniqrig(j,i)=0.0
       mmxuniqpyr(j,i)=0.0
       mmxuniqnat(j,i)=0.0
       mmxdistcentrig(j,i)=0.0
       mmxdistcentpyr(j,i)=0.0
       mmxdistcentnat(j,i)=0.0
    endif
    totbio(j)=totbio(j) + mwt(j,i)  ! not an average here because it summed to get total biomass by species over all cells
  end do
end do


! compute production and yield by species
totp=0.0
do i=1,nspecies
  IF(yieldnum(i).gt.0.0)then
      yieldwt(i)=yield(i)/yieldnum(i)   ! average weight per harvested ind = total biomass of yield / total number harvested
  else
      yieldwt(i)=0.0
  endif

! total production is yield in biomass plus loss due to natural mort and gains or losses due to groiwth in weight (g ww)
  totp(i)=yield(i)+naturalp(i)+growthp(i)    ! ignored naturalpnum - but I think it is correct to ignore changes in numbers ! should be growth plus mortality plsu remaining biomass  KAR
enddo                                        ! because the naturalp terms includes worth, which relfects mortality



!numage is number of pop individuals by age and species
WRITE(13,1009)cumhours,iyear,iday,jday,ihour,inight,(numage(1,j),j=1,inyears),(numage(2,j),j=1,inyears),&
     &  (numage(3,j),j=1,inyears), (numage(6,j),j=1,inyears)
!  mmwt is the average weight per individual fish (g ww) by age class
WRITE(14,1010)cumhours,iyear,iday,jday,ihour,inight,(mmwt(1,j),j=1,inyears),(mmwt(2,j),j=1,inyears),&
     & (mmwt(3,j),j=1,inyears), (mmwt(6,j),j=1,inyears)
! totbio is total biomass over all ages by species; mwt is biomass by age class
WRITE(16,1009)cumhours,iyear,iday,jday,ihour,inight,totbio(1),totbio(2),totbio(3),totbio(6),(mwt(1,j),j=1,inyears),&
     &  (mwt(2,j),j=1,inyears), (mwt(3,j),j=1,inyears), (mwt(6,j),j=1,inyears)
! totp is total production and semms OK  KAR
WRITE(20,1012)cumhours,iyear,iday,jday,ihour,inight,(totp(j),j=1,3),(yield(j),yieldnum(j),yieldwt(j),j=1,3),&
     &  (naturalp(j),j=1,3), (growthp(j),j=1,3),(naturalpnum(j),j=1,3)
! average risk defined as hours off rig cells during nighttime hours
WRITE(26,1010)cumhours,iyear,iday,jday,ihour,inight,(mmxriskhr(1,j),j=1,inyears),(mmxriskhr(2,j),j=1,inyears),&
     &  (mmxriskhr(3,j),j=1,inyears)
WRITE(32,1056)cumhours,iyear,iday,jday,ihour,inight,(mmxdailymove(1,j),j=1,inyears),(mmxdailymove(2,j),j=1,inyears),&
     &  (mmxdailymove(3,j),j=1,inyears), (mmxdailymove(4,j),j=1,inyears), (mmxdailymove(5,j),j=1,inyears), (mmxdailymove(6,j),j=1,inyears)
if(mod(iday,7).eq.0)then
WRITE(33,1057)cumhours,iyear,iday,jday,ihour,inight,(mmxuniqrig(1,j),j=1,inyears),(mmxuniqrig(2,j),j=1,inyears),&
     &  (mmxuniqrig(3,j),j=1,inyears), (mmxuniqrig(4,j),j=1,inyears), (mmxuniqrig(5,j),j=1,inyears), (mmxuniqrig(6,j),j=1,inyears)
WRITE(53,1057)cumhours,iyear,iday,jday,ihour,inight,&
     &  (mmxuniqpyr(1,j),j=1,inyears),(mmxuniqpyr(2,j),j=1,inyears),(mmxuniqpyr(3,j),j=1,inyears), (mmxuniqpyr(4,j),j=1,inyears), (mmxuniqpyr(5,j),j=1,inyears), (mmxuniqpyr(6,j),j=1,inyears)
WRITE(63,1057)cumhours,iyear,iday,jday,ihour,inight,&
     &  (mmxuniqnat(1,j),j=1,inyears),(mmxuniqnat(2,j),j=1,inyears),(mmxuniqnat(3,j),j=1,inyears), (mmxuniqnat(4,j),j=1,inyears), (mmxuniqnat(5,j),j=1,inyears), (mmxuniqnat(6,j),j=1,inyears)
  xuniqrig=0  ! only reset xuniqrig to 0 after a week
  xuniqpyr=0
  xuniqnat=0
end if
WRITE(34,1058)cumhours,iyear,iday,jday,ihour,inight,(mmxdistcentrig(1,j),j=1,inyears),(mmxdistcentrig(2,j),j=1,inyears),&
     &  (mmxdistcentrig(3,j),j=1,inyears),(mmxdistcentrig(4,j),j=1,inyears),(mmxdistcentrig(5,j),j=1,inyears),(mmxdistcentrig(6,j),j=1,inyears)
WRITE(54,1058)cumhours,iyear,iday,jday,ihour,inight,(mmxdistcentpyr(1,j),j=1,inyears),(mmxdistcentpyr(2,j),j=1,inyears),&
     &  (mmxdistcentpyr(3,j),j=1,inyears),(mmxdistcentpyr(4,j),j=1,inyears),(mmxdistcentpyr(5,j),j=1,inyears),(mmxdistcentpyr(6,j),j=1,inyears)
WRITE(64,1058)cumhours,iyear,iday,jday,ihour,inight,(mmxdistcentnat(1,j),j=1,inyears),(mmxdistcentnat(2,j),j=1,inyears),&
     &  (mmxdistcentnat(3,j),j=1,inyears),(mmxdistcentnat(4,j),j=1,inyears),(mmxdistcentnat(5,j),j=1,inyears),(mmxdistcentnat(6,j),j=1,inyears)

1009  FORMAT(1x,e12.5,1x,5(i5,1x),44(e12.5,1x))
1010  FORMAT(1x,e12.5,1x,5(i5,1x),40(f12.2,1x))
1012  FORMAT(1x,e12.5,1x,5(i5,1x),21(e12.5,1x))
1056  FORMAT(1x,e12.5,1x,5(i5,1x),60(f12.5,1x))
1057  FORMAT(1x,e12.5,1x,5(i5,1x),60(f12.5,1x))
1058  FORMAT(1x,e12.5,1x,5(i5,1x),60(f12.5,1x))

! total prey by type on the entire grid by each day
zmean=0.0
do i=1,ncol
  do j=1,nrow
     do k=1,nprey
         zmean(k)=zmean(k)+zprey(i,j,k)  ! sum of each prey type (g ww/m2) over the entire grid
     end do
  END do
END do

!  zmean/(ncol*nrow) is the overal grid average prey density (g ww/m2) by prey type k
WRITE(15,1011)cumhours,iyear,iday,jday,ihour,inight,(zmean(k)/(float(ncol*nrow)),k=1,nprey)
1011  FORMAT(1x,e12.5,1x,5(i5,1x),10(f12.2,1x))

end subroutine
!************************************************************
subroutine bioout
USE allvar
INTEGER xsp
! called every hour - biomass information

movebio=0.0
bluebio=0.0
jackbio=0.0
gabio = 0.0
bio=0.0
do i=1,itotf
    xsp=xspecies(i)
    IF(xalive(i).eq.0)then  ! if i-th individual is alive
        ii=xxcell(i)
        jj=xycell(i)
        bio(xsp,ii,jj) = bio(xsp,ii,jj) + xworth(i)*xweight(i) ! biomass (g ww) by species and cell
        bluebio=bluebio+bio(4,ii,jj) ! running sum of cell by cell bluefish biomass
        jackbio=jackbio+bio(5,ii,jj)
        gabio=gabio+bio(6,ii,jj)
        IF(xsp.le.3)movebio(ii,jj) = movebio(ii,jj) + xworth(i)*xweight(i) ! used to calculate cell by cell total biomass for species 1
    endif
end do

bluebio=bluebio/(ncol*nrow)  !gridwide average of bluefish biomass
jackbio=jackbio/(ncol*nrow)  !gridwide average of jack biomass
gabio=gabio/(ncol*nrow)      !gridwide average of ga biomass 

! write out biomass by species and cell during year five, every 60 days
!IF(mod(iday,5).eq.0.and.iyear.le.1)then
!IF(iyear.le.10.and.iday.eq.10)then
IF(MOD(iyear,5).eq.0.and.MOD(iday,20).eq.0)then
do i=1,ncol
  do j=1,nrow
     IF(ihour.eq.12.or.ihour.eq.24)then
       WRITE(17,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(1,i,j)
       WRITE(22,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(2,i,j)
       WRITE(23,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(3,i,j)
       WRITE(35,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(4,i,j)
       WRITE(36,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(5,i,j)
       WRITE(44,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(6,i,j)
     endif
     IF(iday.eq.240)then
       WRITE(19,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(1,i,j)
       WRITE(24,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(2,i,j)
       WRITE(25,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(3,i,j)
       WRITE(37,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(4,i,j)
       WRITE(38,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(5,i,j)
       WRITE(45,1012)cumhours,iyear,iday,jday,ihour,inight,i,j,bio(6,i,j)
     endif

1012  FORMAT(1x,e12.5,1x,7(i5,1x),e12.5,1x)
  end do
end do
endif

end subroutine
!********************************************************

SUBROUTINE setrig
USE allvar
IMPLICIT NONE

!DEL - commented out print and read statements.  Data is read from file 'input.dat'
!PRINT *,'Method: 1 = MU-775: pyramid field and randomly places rigs, 2 = RGV: random clusters of pyramids and randomly places rigs, 3 = POC: Places 2 large pyramid fields and randomly places rigs, 4 = placing random rigs'
!READ *, meth
!PRINT *, 'Number of rigs'
!READ *, nrigs


!IF (meth.eq.1)CALL pyrfieldandrandrig

!IF (meth.eq.2)CALL pyrclustandrandrig

!IF (meth.eq.3)CALL POCpyrclustandrandrig

!IF (meth.eq.4)CALL randrigonly

DO i = 1,ncol
   DO j = 1,nrow
    IF(zhab(i,j).ne.1) THEN
    write(7,1000)i,j,zhab(i,j)
    1000 FORMAT(1x,3(i3,1x))
    END IF

   END DO
END DO

END SUBROUTINE setrig

! This subroutine will simply call random rigs where natural reefs and pyramids are not place
SUBROUTINE randrigonly
USE allvar
IMPLICIT NONE
INTEGER rigcnt
REAL*8 ran1

rigcnt = 1

!Randomly places rigs on the right side of the grid
    DO i =1,nrigs
5678    rigcol(i)=ran1(idum)*ncol   ! new rig col
        rigrow(i)=ran1(idum)*nrow   ! new rig row
        !IF (rigcol(i).lt.(ncol/3)) GO TO 5678
        IF (zhab(rigcol(i),rigrow(i)).ne.1) then
            GO TO 5678  ! check and see if the cell is already a habitat type. If it is, then start over
    else
        zhab(rigcol(i),rigrow(i))=0
        PRINT *, i, rigcol(i), rigrow(i), zhab(rigcol(i),rigrow(i)) ! rig location not working for some reason and needed to check
    end if
        rigcnt = rigcnt + 1
        
    END DO

END SUBROUTINE randrigonly


! This subroutine creates an evenly spaces pyramid field on the left side of the grid and then randomly places several rigs on the right side of the grid
SUBROUTINE pyrfieldandrandrig
USE allvar
IMPLICIT NONE
INTEGER :: upperpyrcol, lowerrigcol
INTEGER :: bffr, intvl, pyrcnt, rigcnt
REAL*8 ran1

pyrcnt = 1
rigcnt = 1
upperpyrcol = ncol/3
lowerrigcol = ncol/5
bffr = 5
intvl = 2

! Setting up the pyramid field
DO i = bffr,nrow+bffr
    DO j = bffr,upperpyrcol-bffr

    if (MOD(i-bffr,intvl).eq.0.and.MOD(j-bffr,intvl).eq.0) then

        pyrcol(pyrcnt)=j
        pyrrow(pyrcnt)=i
        IF(zhab(pyrcol(pyrcnt),pyrrow(pyrcnt)).ne.1) THEN
           PRINT *, "Trying to place pyramid field over natural bank"
           STOP
        END IF
        zhab(pyrcol(pyrcnt), pyrrow(pyrcnt))=2

        if (pyrcnt.eq.npyr)GO TO 678

        pyrcnt = pyrcnt+1

    end if
    END DO
END DO

678 continue

!Randomly places rigs on the right side of the grid
    DO i =1,nrigs

5678    rigcol=ran1(idum)*ncol   ! new rig col
        rigrow=ran1(idum)*nrow   ! new rig row
        IF (rigcol(i).lt.(ncol/3)) GO TO 5678
        if (zhab(rigcol(i),rigrow(i)).ne.1)then
          GOTO 5678  ! check and see if the cell is already a rig
        else
          zhab(rigcol,rigrow)=0
        endif
        rigcnt = rigcnt + 1
    END DO

END SUBROUTINE pyrfieldandrandrig

! This subroutine creates several random clusters of pyramids and randomly places rigs 
SUBROUTINE POCpyrclustandrandrig
USE allvar
IMPLICIT NONE
INTEGER :: l, m
INTEGER :: pyrcnt, rigcnt, upperpyrcol  !counter and the upper column limit for pyramid placement
INTEGER :: pyrclustcol, pyrclustrow  !column and row of clusters of pyramids
INTEGER :: clustbffr, clustintvl, pyrbffr, pyrintvl !buffer and interval for the clusters and pyramids
INTEGER :: npyrclust, npyrinclust, numclust1, numclust2
REAL*8 ran1

pyrcnt = 1
numclust1 = 6
numclust2 = 3
npyrclust = 1
clustbffr = 1
clustintvl = 1
pyrbffr = 2
pyrintvl = 1
upperpyrcol = ncol/3
npyrinclust = npyr/numclust1

!Placing the Keeping it Wild pyramid clusters
  DO i = 1,numclust1  
     DO j = clustbffr,nrow/2+clustbffr
        DO k = clustbffr,upperpyrcol+clustbffr
              pyrclustcol=k
              pyrclustrow=j
              DO l = pyrclustrow, pyrclustrow+pyrbffr
                 DO m = pyrclustcol, pyrclustcol+pyrbffr
                    IF (MOD(m-pyrbffr,pyrintvl).eq.0.and.MOD(l-pyrbffr,pyrintvl).eq.0) then
                       pyrcol(pyrcnt)=m
                       pyrrow(pyrcnt)=l
                       IF(zhab(pyrcol(pyrcnt),pyrrow(pyrcnt)).ne.1) THEN
                          PRINT *, "Trying to place pyramid field over natural bank"
                          STOP
                       END IF
                       zhab(pyrcol(pyrcnt), pyrrow(pyrcnt))=2
                       pyrcnt = pyrcnt + 1
                       IF (pyrcnt.gt.8)GO TO 678
                    END IF
                 END DO
              END DO
           END DO
        END DO
678     pyrcnt =  1
!        IF (npyrclust.eq.numclust1) GO TO 4001
        npyrclust = npyrclust + 1
   END DO
4001 continue

! Placing the Shell pyramid clusters
npyrclust = 1
clustbffr = 2
clustintvl = 1
npyrinclust = npyr/numclust2
!  DO i = 1,numclust2
!     DO j = clustbffr,nrow+clustbffr
!        DO k = clustbffr,upperpyrcol+clustbffr
!           if (MOD(j-clustbffr,clustintvl).eq.0.and.MOD(k-clustbffr,clustintvl).eq.0) then
!              pyrclustcol=k
!              pyrclustrow=j
!              DO l = pyrclustrow, pyrclustrow+pyrbffr
!                 DO m = pyrclustcol, pyrclustcol+pyrbffr
!                    IF (MOD(l-pyrbffr,pyrintvl).eq.0.and.MOD(m-pyrbffr,pyrintvl).eq.0) then
!                       pyrcol(pyrcnt)=m
!                       pyrrow(pyrcnt)=l
!                       zhab(pyrcol(pyrcnt), pyrrow(pyrcnt))=2
!                       pyrcnt = pyrcnt + 1
!                       IF (npyrclust.eq.1.or.npyrclust.eq.2.and.pyrcnt.gt.30)GO TO 679
!                       IF (npyrclust.eq.3.and.pyrcnt.gt.20)GO TO 679
!                    END IF
!                 END DO
!              END DO
!           END IF
!679     pyrcnt =  1
!        END DO
!     END DO
!        IF (npyrclust.eq.numclust2) GO TO 4002
!        npyrclust = npyrclust + 1
!  END DO
!4002 continue

!Randomly places rigs
    DO i =1,nrigs

5678    rigcol=ran1(idum)*ncol   ! new rig col
        rigrow=ran1(idum)*nrow   ! new rig row
        IF (rigcol(i).lt.(ncol/3)) GO TO 5678

        if (zhab(rigcol(i),rigrow(i)).ne.1)then
          GOTO 5678  ! check and see if the cell is already a rig

        else

          zhab(rigcol,rigrow)=0

        endif
        rigcnt = rigcnt + 1
    END DO
END SUBROUTINE POCpyrclustandrandrig


! This subroutine creates two large pyramid fields similar to POC with a large number of pyramids very close together
SUBROUTINE pyrclustandrandrig
USE allvar
IMPLICIT NONE
INTEGER :: rem  !remainder
INTEGER :: pyrcnt, rigcnt, upperpyrcol  !counter and the upper column limit for pyramid placement
REAL(8) :: pyrclustcol, pyrclustrow  !column and row of randomly placed clusters of pyramids
REAL(8) :: npyrclust, numclust, bffr, intvl, thispyrcol(nrow*ncol), thispyrrow(nrow*ncol)
REAL*8 ran1

pyrcnt = 1
numclust = npyr/4
npyrclust = 1
bffr = 0.5
intvl = 0.5
upperpyrcol = ncol/3
rem = npyr - numclust*4
IF (rem .ne. 0) THEN
PRINT *, "Number of pyramids is not a multiple of 4"
STOP
END IF 

  DO i = 1,numclust
    pyrclustcol=ran1(idum)*upperpyrcol   ! new pyramid cluster col
    pyrclustrow=ran1(idum)*nrow   ! new pyramid cluster row
     DO j = pyrclustrow, pyrclustrow+bffr
        DO k = pyrclustcol, pyrclustcol+bffr
!           IF (MOD(k-bffr,intvl).eq.0.and.MOD(j-bffr,intvl).eq.0) then
              thispyrcol(pyrcnt)=k
              thispyrrow(pyrcnt)=j
              IF(zhab(thispyrcol(pyrcnt),thispyrrow(pyrcnt)).ne.1) THEN
                   PRINT *, "Trying to place pyramid field over natural bank"
                   STOP
              END IF
              zhab(thispyrcol(pyrcnt), thispyrrow(pyrcnt))=2
              pyrcnt = pyrcnt + 1
              IF (pyrcnt.gt.4)GO TO 678
 !          END IF
        END DO
     END DO
678     pyrcnt =  1
        IF (npyrclust.eq.numclust) GO TO 4001
        npyrclust = npyrclust + 1
   END DO
4001 continue

!Randomly places rigs
    DO i =1,nrigs

5678    rigcol=ran1(idum)*ncol   ! new rig col
        rigrow=ran1(idum)*nrow   ! new rig row
        IF (rigcol(i).lt.(ncol/3)) GO TO 5678

        if (zhab(rigcol(i),rigrow(i)).ne.1)then
          GOTO 5678  ! check and see if the cell is already a rig

        else

          zhab(rigcol,rigrow)=0

        endif
        rigcnt = rigcnt + 1
    END DO
END SUBROUTINE pyrclustandrandrig

!********************************************************
!subroutine to physically insert desired rig columns and rows MDC - 14april09
SUBROUTINE pinsert
USE allvar
IMPLICIT NONE

    !rigcol = (/44,36,28,20,12,4,40,32,24,16,8,44,36,28,20,12,4,40,32,24,16,8,44,36,28,20,&
    !12,4,40,32,24,16,8,44,36,28,20,12,4,40,32,24,16,8,44,36,28,20,12,4/)
    !rigrow = (/5,5,5,5,5,5,10,10,10,10,10,15,15,15,15,15,15,20,20,20,20,20,25,25,25,25,25&
    !25,30,30,30,30,30,35,35,35,35,35,35,40,40,40,40,40,45,45,45,45,45,45/)

    !DO i = 1,ncol
    !    DO j = 1,nrow
    !        DO k = 1,nrigs
    !            IF(rigcol(k).eq.i.and.rigrow(k).eq.j)THEN
    !            zhab(i,j) = 0
    !            END IF
    !        END DO
    !    END DO
    !END DO

END SUBROUTINE pinsert

!*********************************************************************************************

!NEW - subroutine uses LOCAL clsize and #rigs to place 1 cluster at the center of the grid
SUBROUTINE centclust
USE allvar
IMPLICIT NONE
INTEGER :: clsize,clsize2 !cluster size, and remainder cluster size to fill
INTEGER :: rem  !remainder
INTEGER :: center = nrow/2  !center of the grid
INTEGER :: cnt  !counter
INTEGER :: ic, jr  !adjusts clsize to the center of the grid
REAL(8) :: chknrigs

cnt = 1
chknrigs = nrigs
!PRINT *, '1st check chknrigs', chknrigs

    IF (nrigs.lt.9)THEN
    clsize2 = 1
    rem = nrigs
    GO TO 55
    END IF

    clsize = (SQRT(chknrigs)-1)/2.0
!    PRINT *, '2nd check clsize', clsize

44  IF (MOD((SQRT(chknrigs)-1),2.0).eq.0) THEN !checks to see if nrigs is an odd perfect square
    clsize = (SQRT(chknrigs)-1)/2.0
    rem = nrigs - chknrigs
        IF (rem.gt.0) THEN
        clsize2 = clsize + 1
        END IF
!    PRINT *, 'clsize, rem', clsize, rem
    ELSE
    chknrigs = chknrigs - 1
    GO TO 44
    END IF

! If there are more rigs than cells stop the program.
    IF (nrigs.eq.ncol*nrow) THEN
    PRINT *, 'rigs equal to or greater than cells available - no where to put food'
    stop
    END IF

    DO i = -clsize,clsize
        ic = center + i
        DO j = -clsize,clsize
        jr = center + j

        rigcol(cnt)=ic
        rigrow(cnt)=jr
        zhab(rigcol(cnt),rigrow(cnt)) = 0

        cnt = cnt + 1
        END DO
    END DO

55  IF (rem.gt.0) THEN
        DO i = -clsize2,clsize2
        ic = center + i
            DO j = -clsize2,clsize2
            jr = center + j

                IF (zhab(ic,jr).gt.0) THEN
                
                rigcol(cnt)=ic
                rigrow(cnt)=jr
                zhab(rigcol(cnt),rigrow(cnt)) = 0

                cnt = cnt + 1
                rem = rem-1
                END IF

                IF (rem.eq.0)GO TO 4000
                
            END DO
        END DO
    END IF
4000 continue

END SUBROUTINE centclust

!*********************************************************************************************

! uniform, sequentially spaced

SUBROUTINE uniform
USE allvar
IMPLICIT NONE
INTEGER :: cnt, bffr, intvl       ! Counter
INTEGER :: rem       ! Remainder check
INTEGER :: rowcount  ! Trims number of column places available in a row to 1 - square root of nrigs
!INTEGER :: intvl, bffr

cnt = 1

!DEL - commented out print and read statements.  Data is read from file 'input.dat'  MDC 16apr09
PRINT *, 'Set interval, buffer, and rowcount (0=buffer, 1=no buffer, boolean flag) - check the interval against total rigs desired'
READ *, intvl, bffr, rowcount

!rowcount = ((sqrt(float(nrigs)))-1)*intvl

!DO i = 1,ncol
!    DO j = 1,nrow
!
!    if (MOD(i,intvl).eq.0.and.MOD(j,intvl).eq.0) then
!
!        rigcol(cnt)=j-bffr
!        rigrow(cnt)=i-bffr
!        zhab(rigcol(cnt), rigrow(cnt))=0
!
!        if (cnt.eq.nrigs)GO TO 678
!
!        cnt = cnt+1
!
!    end if
!    END DO
!END DO

DO i = bffr,rowcount+bffr
    DO j = bffr,rowcount+bffr

    if (MOD(i-bffr,intvl).eq.0.and.MOD(j-bffr,intvl).eq.0) then

        rigcol(cnt)=j
        rigrow(cnt)=i
        zhab(rigcol(cnt), rigrow(cnt))=0

        if (cnt.eq.nrigs)GO TO 678

        cnt = cnt+1

    end if
    END DO
END DO

678 continue
END SUBROUTINE uniform
!********************************************************************************************

!NEW - subroutine to physically insert desired rig columns and rows MDC - 14april09
SUBROUTINE randrig
USE allvar
IMPLICIT NONE
REAL*8 ran1


    DO i =1,nrigs
    
5678    rigcol=ran1(idum)*ncol   ! new rig col
        rigrow=ran1(idum)*nrow   ! new rig col
        if (zhab(rigcol(i),rigrow(i)).eq.0)then

          GOTO 5678  ! check and see if the cell is already a rig

        else

          zhab(rigcol,rigrow)=0

        endif

    END DO

END SUBROUTINE randrig

!********************************************************************************************
! Subroutine to randomly add rigs to the grid and update rig info
SUBROUTINE addrig
USE allvar
IMPLICIT NONE
REAL*8 ran1

    DO i =nrigs+1,nrigs+newrigs
    
5678    rigcol(i)=ran1(idum)*ncol   ! new rig col
        rigrow(i)=ran1(idum)*nrow   ! new rig col
        if (zhab(rigcol(i),rigrow(i)).eq.0) GOTO 5678  ! check and see if the cell is already a rig
        zhab(rigcol(i),rigrow(i))=0   ! If the cell is open, make it a rig cell

    END DO

    nrigs=nrigs+newrigs ! update nrigs to include newly added rigs

    DO i = 1,nrigs
      write(39,1000)i,rigcol(i),rigrow(i),zhab(rigcol(i),rigrow(i))
      1000 FORMAT(1x,4(i3,1x))
    END DO

!    DO i = 1,ncol
!      DO j = 1,nrow
!
!      IF(zhab(i,j).eq.0) THEN
!      write(39,1000)i,j,rigcol(i),rigrow(j),zhab(i,j)
!      1000 FORMAT(1x,3(i3,1x))
!      END IF
!
!      END DO
!    END DO

END SUBROUTINE addrig

!********************************************************************************************
! Subroutine to randomly remove rigs from the grid and update rig info
SUBROUTINE subrig
USE allvar
INTEGER, PARAMETER:: maxcells=100   ! maximum number of cells in a neighborhood
INTEGER xxout(maxcells),xxout1(maxcells)
REAL*4 xx(maxcells)    ! note not double precision so matches the dummy array in the numerical recipes

DO i=1,remrigs

    DO j = 1,nrigs
      xx(j)=ran1(idum)
    END DO

    call indexx_sp(xx,xxout1)   ! new version from numerical receipes for fortran90 - xxout1 is integers
    xxout=rank(xxout1)          ! now xxout is the rank values of the real numbers in the array xx
    icell=xxout(1) ! icell is the first random integer in the list, also the rig chosen randomly to be removed

    DO k = 1,nrigs
      if (k.eq.icell)then
        zhab(rigcol(icell),rigrow(icell))=1 ! make habitat benthic again
        rigcol(k)=rigcol(k+1) ! move everybody else up in the list of rigs
        rigrow(k)=rigrow(k+1)
      end if
    END DO

    nrigs = nrigs-1  ! update nrigs to match current status

END DO


DO i = 1,nrigs
  write(39,1000)i,rigcol(i),rigrow(i),zhab(rigcol(i),rigrow(i))
  1000 FORMAT(1x,4(i3,1x))
END DO

!    DO i = 1,ncol
!      DO j = 1,nrow
!        IF(zhab(i,j).eq.0) THEN
!        write(39,1000)i,j,zhab(i,j)
!        1000 FORMAT(1x,3(i3,1x))
!        END IF
!      END DO
!    END DO

END SUBROUTINE subrig

!********************************************************
REAL(8) FUNCTION RAN1(IDUM)  ! generates a uniform 0 to 1 random number
INTEGER idum
DIMENSION R(97)
INTEGER, SAVE:: IX1, IX2, IX3	!KPE added declaration and SAVE
INTEGER:: J
PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
PARAMETER (M3=243000,IA3=4561,IC3=51349)
DATA IFF /0/
SAVE R				!KPE added SAVE
!      PRINT *,'idum=',idum

IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
    IFF=1
    IX1=MOD(IC1-IDUM,M1)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX2=MOD(IX1,M2)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX3=MOD(IX1,M3)
    DO 11 J=1,97
       IX1=MOD(IA1*IX1+IC1,M1)
       IX2=MOD(IA2*IX2+IC2,M2)
       R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11  CONTINUE
    IDUM=21
ENDIF

IX1=MOD(IA1*IX1+IC1,M1)
IX2=MOD(IA2*IX2+IC2,M2)
IX3=MOD(IA3*IX3+IC3,M3)
J=1+(97*IX3)/M3
IF(J.GT.97.OR.J.LT.1)THEN
   print *,'problem in ran1, J=',J
   PAUSE
ENDIF

RAN1=R(J)
R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
RETURN
END

!-----------------------------------------------------------
REAL*8 FUNCTION wtemp(x)
!-----------------------------------------------------------
! TableCurve C:\smelt-data\wtemp.f90 Oct 4, 2007 7:28:57 PM 
!  
! X=  
! Y=  
! Eqn# 8055  Beta(a,b,c,d,e,f) 
! r2=0.967502946789219D0 
! r2adj=0.9350058935784381D0 
! StdErr=1.10092368110586D0 
! Fval=41.680829234558D0 
! a= 14.18865521423792D0 
! b= 11.4758914163255D0 
! c= 241.8974614712954D0 
! d= 314.3232001489511D0 
! e= 2.926681291108284D0 
! f= 2.239830453865759D0 
! Constraints: x>(c*t-d*e+d)/t or x<(c*t+d*f-d)/t t=e+f-2,d>0,e>1,f>1 
!-----------------------------------------------------------
  REAL*8 x,y
  REAL*8 n 
  REAL*8 m 
  n=(2.239830453865759D0 - 1.0)/(2.926681291108284D0+&
   &2.239830453865759D0 - 2.0)
  m=(2.926681291108284D0 - 1.0)/(2.926681291108284D0+&
   &2.239830453865759D0 - 2.0)
  y=14.18865521423792D0+11.47589141632550D0*((x-&
   &241.8974614712954D0+314.3232001489511D0*m)/&
   &314.3232001489511D0)**(2.926681291108284D0 - 1.0)*&
   &(1.0-(x-241.8974614712954D0+314.3232001489511D0*&
   &m)/314.3232001489511D0)**&
   &(2.239830453865759D0 - 1.0)/(m**(2.926681291108284D0 - 1.0)*n**&
   &(2.239830453865759D0 - 1.0))

  wtemp=y

  RETURN
END

!-----------------------------------------------------------
REAL*8 FUNCTION eqn8003(x)
!-----------------------------------------------------------
! TableCurve C:\redsnapper\8003.f90 Oct 4, 2007 8:07:25 PM 
!  
! X=  
! Y=  
! Eqn# 8003  Gaussian(a,b,c,d) 
! r2=0.9440600890272728D0 
! r2adj=0.9160901335409092D0 
! StdErr=1.273866373332105D0 
! Fval=50.62897344371202D0 
! a= 13.52108336438569D0 
! b= 12.52503105456064D0 
! c= 231.5464390570307D0 
! d= 75.24108463989667D0 
! Constraint: d>0 
!-----------------------------------------------------------
  REAL*8 x,y
  REAL*8 n 
  n=(x-231.5464390570307D0)/75.24108463989667D0 
  y=13.52108336438569D0+12.52503105456064D0*DEXP(-0.5*n*n) 
  eqn8003=y
  RETURN
END
!------------------------------------------------------*
! Function where y represents the number of hours of sunlight in a day as a function of day. Where x=0 is Jan 1 and x=365 is Dec 31
REAL*8 FUNCTION eqn9001(x)
  REAL*8 y,x
  y=(-0.00000000000274D0*(x**6))+(0.00000000286590D0*(x**5))-&
  & (0.00000092186564D0*(x**4))+(0.00006351930143D0*(x**3))+&
  & (0.00720366012820D0*(x**2))+(0.59679045015946D0*(x))+631.1331924037070D0
  eqn9001=y/60.0
  RETURN
END
 
!--------------------------------------------
! Function describing the changing home range (or foraging area) with size, where y=nhood and x=fish weight
! This linear function has been taken from Figure 9 in Paraino & Szedlmayer 2014 and modified to be used with weight and not size with equation 1 from the results in Wilson & Nieland 2001
REAL*8 FUNCTION eqn9002(x)
  REAL*8 x,y
  y=0.00365*((x/(1.17*10**(-8)))**(1/3.04))-0.74
  eqn9002=y
  RETURN
END

!---------------------------------------------
SUBROUTINE INDEXXold(N,ARRIN,INDX)
      real*8 ARRIN(N)
!---changed foloowing from 2 to 4
      integer*4 INDX(N)
      DO 11 J=1,N
      INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
      IF(L.GT.1)THEN
      L=L-1
      INDXT=INDX(L)
      Q=ARRIN(INDXT)
      ELSE
      INDXT=INDX(IR)
      Q=ARRIN(INDXT)
      INDX(IR)=INDX(1)
      IR=IR-1
      IF(IR.EQ.1)THEN
       INDX(1)=INDXT
       RETURN
      ENDIF
      ENDIF
      I=L
      J=L+L
20    IF(J.LE.IR)THEN
      IF(J.LT.IR)THEN
       IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
      ENDIF
      IF(Q.LT.ARRIN(INDX(J)))THEN
       INDX(I)=INDX(J)
       I=J
       J=J+J
      ELSE
       J=IR+1
      ENDIF
      GO TO 20
      ENDIF
      INDX(I)=INDXT
      GO TO 10
      END
!***************************************************************************
! -- new numerical recipes
 FUNCTION rank(index)
 USE nrtype; USE nrutil, ONLY : arth
IMPLICIT NONE
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: index
INTEGER(I4B), DIMENSION(size(index)) :: rank
!Given index as output from the routine indexx, this routine returns a same-size array
!rank, the corresponding table of ranks.
rank(index(:))=arth(1,1,size(index))
END FUNCTION rank

	SUBROUTINE indexx_sp(arr,index)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	n=assert_eq(size(index),size(arr),'indexx_sp')
	index=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=index(j)
				a=arr(indext)
				do i=j-1,l,-1
					if (arr(index(i)) <= a) exit
					index(i+1)=index(i)
				end do
				index(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(index(k),index(l+1))
			call icomp_xchg(index(l),index(r))
			call icomp_xchg(index(l+1),index(r))
			call icomp_xchg(index(l),index(l+1))
			i=l+1
			j=r
			indext=index(l+1)
			a=arr(indext)
			do
				do
					i=i+1
					if (arr(index(i)) >= a) exit
				end do
				do
					j=j-1
					if (arr(index(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(index(i),index(j))
			end do
			index(l+1)=index(j)
			index(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (arr(j) < arr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_sp

	SUBROUTINE indexx_i4b(iarr,index)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	INTEGER(I4B) :: a
	INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	n=assert_eq(size(index),size(iarr),'indexx_sp')
	index=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=index(j)
				a=iarr(indext)
				do i=j-1,l,-1
					if (iarr(index(i)) <= a) exit
					index(i+1)=index(i)
				end do
				index(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(index(k),index(l+1))
			call icomp_xchg(index(l),index(r))
			call icomp_xchg(index(l+1),index(r))
			call icomp_xchg(index(l),index(l+1))
			i=l+1
			j=r
			indext=index(l+1)
			a=iarr(indext)
			do
				do
					i=i+1
					if (iarr(index(i)) >= a) exit
				end do
				do
					j=j-1
					if (iarr(index(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(index(i),index(j))
			end do
			index(l+1)=index(j)
			index(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (iarr(j) < iarr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_i4b

!*****************************************************

!CP1 --- removed from subroutine setupgrid

!! column and row numbers of rig cells
!rigcol(1)=2; rigrow(1)=2   ;  zhab(2,2)=0          ! no prey on rig cells
!rigcol(2)=45; rigrow(2)=45 ;  zhab(10,10)=0       ! *DEL - commented out this line because cell positions and zhab location does not agree. MDC 17March09
!rigcol(2)=45; rigrow(2)=45 ;  zhab(45,45)=0        ! *FIX - Changed to match.  45,45 = 45,45 MDC 17March2009
!rigcol(3)=15; rigrow(3)=15 ;  zhab(15,15)=0
!rigcol(4)=20; rigrow(4)=20 ; zhab(20,20)=0
!rigcol(5)=20; rigrow(5)=10 ; zhab(20,10)=0
!rigcol(6)=10; rigrow(6)=20 ; zhab(10,20)=0
!rigcol(7)=30; rigrow(7)=30 ; zhab(30,30)=0
!rigcol(8)=40; rigrow(8)=30 ; zhab(40,30)=0
!rigcol(9)=30; rigrow(9)=40 ; zhab(30,40)=0
!rigcol(10)=10; rigrow(10)=10 ; zhab(45,45)=0      ! *DEL - commented out this line because cell position and zhab location do not agree.  MDC 20MARCH09
!rigcol(10)=10; rigrow(10)=10 ; zhab(10,10)=0       ! *FIX - changed to match. 10,10 = 10,10  MDC 17March2009
!rigcol(11)=30; rigrow(11)=33 ; zhab(30,33)=0
!rigcol(12)=30; rigrow(12)=34 ; zhab(30,34)=0

!*****************************************************

