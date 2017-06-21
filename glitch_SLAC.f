c  ******************************************************************
c
c    Program GLITCH  --  Simultaneous reflection mapping program
c
c    This program finds the geometric conditions under which diffrac-
c  tion from secondary sets of crystal planes can occur simultaneously
c  with diffraction from a given set of planes.  This is presented
c  as a display of curves in a two-dimensional space whose
c  axes are the energy of the diffracting radiation and the azimuthal
c  angle phi of the crystal (measuring rotation of the crystal
c  around the given scattering vector).  Each curve represents
c  the locus of conditions under which a particular set of secondary
c  planes can diffract simultaneously with the given planes.
c
c  
c  Input:  This program looks for input concerning the crystal type,
c            primary reflection, and range of energy-phi space to 
c            examine in a text file called GLITCH.DAT   From
c            this file, the program takes the following data:
c
c          Crys_name is an identifier in a file of lattice parameters
c          T is the crystal temperature, in degrees K.
c          Hz(3) gives the (h,k,l) for the primary reflection.
c          Tm(3) is a reciprocal lattice vector that lies in the 
c                 scattering plane when phi=0.  Defines phi=0.
c          Emin and Emax are the limits for the energy range to be
c                 considered.  Values less than 100 will be interpreted
c                 as wavelength values in angstroms, and will be conver-
c                 ted to energy values internally.  Values greater than
c                 100 will be interpreted as energy values, in eV.
c          Phii and Phif are the limits for the phi range, in degrees.
c          Nphi is the number of phi values, between the phi limits,
c                 to consider while searching for glitches. 
c
c      While examining reciprocal lattice vectors for possible simultan-
c        eous reflections, this program normally stops when the magnitudes
c        of the reciprocal vectors become very large, since the reflec-
c        tions they could cause would be so weak as to have negligible 
c        effect on the primary beam.  The point at which the program 
c        stops is controlled by the constant Max, equal to the sum of the
c        squares of the indices for the largest reciprocal lattice 
c        vector that will be examined.  The default value of Max is
c        100, but when the program is run it offers the opportunity to
c        change the value of Max.
c
c
c  Output:     In addition to the graphical output mentioned above, this 
c           program identifies the (h,k,l) indices of each curve on the 
c           energy vs phi graph, so that the lines on the graph can be
c           properly identified.  
c
c           When the incident beam simultaneously satisfies the Bragg condi-
c           tion for the given reciprocal lattice vector P and for a
c           secondary reciprocal lattice vector H, the primary diffracted
c           beam will also satisfy the Bragg condition for the reciprocal
c           lattice vector H-P.  Thus a glitch usually involves at least
c           two secondary reciprocal lattice vectors, and conversely
c           a given secondary reciprocal lattice vector will usually
c           contribute to two glitches at different energies.
c               It is customary to describe a simultaneous reflection 
c           situation by counting the number of radiation beams involved.
c           A "three-beam" case can be described by the origin O=(0,0,0),
c           the given reciprocal lattice vector P, and a third reciprocal
c           lattice vector H, whose Bragg condition is satisfied simul-
c           taneously with that of P.  In solving the diffraction problem,
c           it will be necessary to consider not only O, P, and H, but also
c           -P, -H, H-P, and P-H.  However, the problem is completely
c           specified by O, P, and H.
c               This program finds every contributing vector to each glitch,
c           including forbidden reflections (which become "active" via the
c           detour or Umweg path).  Thus, the secondary reflection H may be
c           forbidden, if the complement H-P is allowed.
c
c
c
c  This program and its subroutines call on some routines from the MTRX
c    library of matrix arithmetic subroutines.
c  This program generates graphical output through the GRAPH subroutine
c    attached at the end of the program.  Since plotting libraries vary
c    from one machine to another, most users will find it necessary to
c    rewrite this subroutine.
c
c
c    The technique for finding multiple-reflection conditions used in this
c  program follows the approach of Rek, et al. (Z.U. Rek, G.S. Brown, and
c  T. Troxel, in "EXAFS and Near Edge Structure III", edited by K.O. Hodgson,
c  B. Hedman, and J.E. Penner-Hahn, Springer, 1984.)  Central to this tech-
c  nique are two equations relating the Bragg angles for the secondary and
c  primary reflections during simultaneous reflection:
c
c    sin(thetaS) = (hh/hhz)*sin(thetaP) = alph1*cos(thetaP)*cos(phi) +
c                          alph2*cos(thetaP)*sin(phi) +- alph3*sin(thetaP)
c
c  ThetaS is the secondary reflection Bragg angle, thetaP is the given
c  reflection Bragg angle, and the other variables are as defined in this
c  program.  The +- sign is + when the incident beam satisfies the secon-
c  dary Bragg condition, and - when the primary diffracted beam satisfies
c  the secondary Bragg condition.
c
c    VERSION 4.0      J. Arthur, SLAC, June 2017
c
c  ******************************************************************
      real*8 alatt(3),angle(3,3),hx(3),hy(3),hz(3),tm(3),emin,
     z       emax,phii,phif,hc,s,hhx,hhy,hhz,hh,h(3),elow(3,3),
     z       e(3,3),ee(3,3,3),eup(3,3),e0,phi,alph1,alph2,alph3,
     z       small,t,temp,tabl(512,5),dh,hpri(3),hsec(3),
     z       a,test1,test2,sinth(2),tanth,egl(2),eglch,dot,alamda,
     z       glitchlist(100,6)
      real*4 max,sumsq,x(2,512),y(2,512),xmin,xmax,ymin,ymax,
     z       xx(512),yy(512)
      complex*16 psi,fact
      character*25 elabl,philbl,crys_name
      character*1 ianser
      integer ihz(3),npts(2),ioptn,nchrs,itest,lim,lamflg,ncurv
      common /cc1/e,ee
      data small,s,hc/1.0d-3,1.7453293d-2,1.239852d4/
c
c  Read the input parameters from the file XTAL1.DAT
c 
      open(unit=1,file='GLITCH.DAT',action='READ',status='OLD',
     z      form='FORMATTED',err=1)
      goto 2
 1    write(*,1000)
1000  format(/' Error opening GLITCH.DAT')
      call exit
 2    continue
      read(1,*,end=5,err=5)crys_name
      read(1,*,end=5,err=5)t
      read(1,*,end=5,err=5)hz(1),hz(2),hz(3)
      read(1,*,end=5,err=5)tm(1),tm(2),tm(3)
      read(1,*,end=5,err=5)emin,emax
      read(1,*,end=5,err=5)phii,phif,nphi
      goto 7
 5    write(*,1010)
1010  format(/' Error reading GLITCH.DAT')
      call exit
 7    continue
c
c  Close the input data file
c
      close(1)
c
c  Determine value for Max
c
      max = 100.0
      write(*,1020)max
1020  format(/' Maximum h**2 + k**2 + l**2 =',f6.0,2x,
     z       'Change? (CR=N): ',$)
      read(*,1030)ianser
1030  format(a1)
      if(ianser .ne. 'Y' .and. ianser .ne. 'y')goto 8
      write(*,1040)
1040  format(' New value for Max: ',$)
      read(*,*)max
 8    lim = sqrt(max) + 1   !Largest h, k, l that need be considered
c
c  Set up crystal lattice common data.  (This subroutine
c    finds the crystal lattice parameters in a data file and
c    loads the common block with data useful in performing
c    vector products in the crystal lattice.)
c
      call lattice(crys_name,t,elow,eup,e0)
c
c  LATTICE calls the subroutine SETCRYS, which
c  provides crystal lattice info.  
c
c  Now find basis vectors Hx and Hy at right angles to Hz. 
c    (Hx is in the scattering plane when phi=0.) 
c
      call cross(hz,tm,hy,elow,e0)
      call cross(hy,hz,hx,elow,e0)
      hhz = dsqrt(dot(hz,hz,eup))
      hhx = dsqrt(dot(hx,hx,eup))
      hhy = dsqrt(dot(hy,hy,eup))
c
c  Prepare data for output plot, with proper units and precision
c
      elabl = (' Energy (eV)')
      if(emin .lt. 100.d0) elabl=(' Wavelength')
      philbl = (' Phi (deg)')
      xmin = emin		!plot package prefers single precision
      xmax = emax
      ymin = phii
      ymax = phif
c
c  If wavelength units were given, change to energy units
c 
      lamflg = 1
      if(emin .gt. 100.d0)goto 10
      lamflg = 0
      temp = emin
      emin = hc/emax
      emax = hc/temp
 10   ncurv = 0            !Number of glitches found
c
c  Set up loops through h,k,l values for possible glitches
c
      do 50 k1=1,2                           !For pos, neg h values
        do 49 i1=1,lim
          if(i1.eq.1.and.k1.eq.2)goto 49     !Don't re-count (0,k,l)
          ih = (i1-1)*(2*k1 - 3)
          do 45 k2=1,2                       !For pos, neg k values
            do 44 i2=1,lim
              if(i2.eq.1.and.k2.eq.2)goto 44 !Don't re-count (h,0,l)
              ik = (i2-1)*(2*k2 - 3)
              do 40 k3=1,2                   !For pos, neg l values
                do 39 i3=1,lim
                  if(i3.eq.1.and.k3.eq.2)goto 39   !Ditto (h,k,0)
                  il = (i3-1)*(2*k3 - 3)
                  sumsq = ih**2 + ik**2 + il**2
                  if(sumsq .gt. max)goto 39
                  if(sumsq .lt. small)goto 39   !Don't count (0,0,0)
                  do 12 j=1,3
                    ihz(j) = hz(j) + 0.1  !Find integer h,k,l for hz
 12                 if(hz(j).lt.0.d0) ihz(j) = hz(j) - 0.1
                  if(ih.eq.ihz(1) .and. ik.eq.ihz(2) .and.
     z               il.eq.ihz(3))goto 39     !Avoid primary reflec
c
c  Now (ih,ik,il) has been selected to test for possible glitch
c  First, see if it gives an allowed reflection
c
                  h(1) = dfloat(ih)
                  h(2) = dfloat(ik)
                  h(3) = dfloat(il)
                  dh = dsqrt(dot(h,h,eup))
                  alamda = hc/emax
                  call scatamp(h,fact)
c
c  The subroutine SCATAMP provides enough information about the 
c  scattering amplitude of H to determine whether or not
c  H is an allowed reflection.  If H is forbidden, we skip it
c  at this point
c
                  if(cdabs(fact) .lt. small)goto 39
                  hh = dsqrt(dot(h,h,eup))
                  alph1 = dot(h,hx,eup)/hh/hhx
                  alph2 = dot(h,hy,eup)/hh/hhy
                  alph3 = dot(h,hz,eup)/hh/hhz
c
c  Now loop through the phi range.  If a glitch condition occurs
c    within the desired energy limits, find the energy vs phi line
c    that describes it, and prepare to draw it on the glitch map.
c
                  npts(1) = 0   !# points in the line(s) to be drawn
                  npts(2) = 0
                  do 30 jj=1,nphi
                    phi = phii + (phif-phii)/nphi*float(jj-1)
                    phi = phi*s
                    a = alph1*dcos(phi) + alph2*dsin(phi)
                    if(dabs(hh/hhz - dabs(alph3)) .gt. small)goto 13
c
c   When hh/hhz = abs(alph3), a glitch occurs with primary
c   Bragg angle equal to 90 deg
c
                    sinth(1) = 1.d0
                    sinth(2) = 1.d0
                    goto 18
 13                 itest = 0
                    test1 = a/(hh/hhz - alph3)
                    test2 = a/(hh/hhz + alph3)
                    if(test1 .gt. 0.d0)itest = itest + 1
                    if(test2 .gt. 0.d0)itest = itest + 1
                    if(dabs(a) .lt. small)itest = itest + 10
                    if(itest .gt. 0)goto 14
c
c   ITEST=0  The glitch condition is impossible at any energy.
c
                    goto 30
 14                 continue
                    if(itest .ne. 2)goto 15
c
c   ITEST=2  There are two glitch energies.  At one energy, the  
c     incident beam satisfies the secondary Bragg condition.  At the 
c     other energy, the primary diffracted beam satisfies the 
c     secondary Bragg condition.
c
                    tanth = test1
                    sinth(1) = tanth/dsqrt(tanth**2 + 1.d0)
                    tanth = test2
                    sinth(2) = tanth/dsqrt(tanth**2 + 1.d0)
                    goto 18
 15                 continue
                    if(itest .eq. 1)goto 16
c
c   ITEST>2  The glitch condition occurs only at infinite energy.
c
                    goto 30
 16                 continue
c
c   ITEST=1  A single glitch condition exists at finite energy.
c
                    if(test1 .lt. 0.d0)goto 17
c
c   TEST1>0  The glitch occurs when the incident beam satisfies the
c     secondary Bragg condition for reciprocal lattice vector H.
c
c   TEST2>0  The glitch occurs when the primary diffracted beam 
c     satisfies the secondary Bragg condition for vector H-P.
c
                    tanth = test1
                    sinth(1) = tanth/dsqrt(tanth**2 + 1.d0)
                    sinth(2) = 0.d0
                    goto 18
 17                 tanth = test2
                    sinth(1) = 0.d0
                    sinth(2) = tanth/dsqrt(tanth**2 + 1.d0)
 18                 continue
c
c  At this point, one or possibly two glitches have been found
c  involving a particular combination of h,k,l and phi values.
c  Test to see if the energies of these glitches are within
c  the desired range, and if so store the E and phi values
c  for plotting.
c  
                    do 20 j=1,2
                      if(sinth(j) .lt. small)goto 20
                      eglch = hc*hhz/2.d0/sinth(j)
                      if(eglch.lt.emin .or. eglch.gt.emax)goto 20
                      egl(j) = eglch
                      npts(j) = npts(j) + 1
                      x(j,npts(j)) = egl(j)
                      y(j,npts(j)) = phi/s
 20                   continue
 30                 continue
c
c  End phi loop
c 
                 do 35 j=1,2
                    if(npts(j) .eq. 0)goto 35
c 
c  Get h,k,l values for reciprocal lattice vectors H and H-P
c
					if(j .eq. 1) then
						hpri(1) = h(1)
						hpri(2) = h(2)
						hpri(3) = h(3)
						hsec(1) = h(1) - hz(1)
						hsec(2) = h(2) - hz(2)
						hsec(3) = h(3) - hz(3)
					elseif(j .eq. 2) then
						hpri(1) = h(1) + hz(1)
						hpri(2) = h(2) + hz(2)
						hpri(3) = h(3) + hz(3)
						hsec(1) = h(1)
						hsec(2) = h(2)
						hsec(3) = h(3)
					endif
					if(ncurv .lt. 1) goto 33
c 
c  Check to see if new glitch is a duplicate of one already found
c  Most glitches should be found twice, once when H lies on the primary
c  beam Ewald sphere, and once when H-P lies on the scattered beam
c  Ewald sphere.  But we only want to count them once
c
					do 32 i=1,ncurv
						if(hpri(1).eq.glitchlist(i,1) .and. 
     z                  	hpri(2).eq.glitchlist(i,2) .and.
     z                      hpri(3).eq.glitchlist(i,3)) goto 35
  32                continue
  33				ncurv = ncurv + 1
 					if(ncurv .gt. 100) then
 					  write(*,1045)
1045  format(/' 100 glitch limit has been exceeded')
					  goto 51
 					endif
 					glitchlist(ncurv,1) = hpri(1)
 					glitchlist(ncurv,2) = hpri(2)
 					glitchlist(ncurv,3) = hpri(3)
 					glitchlist(ncurv,4) = hsec(1)
 					glitchlist(ncurv,5) = hsec(2)
 					glitchlist(ncurv,6) = hsec(3)					
c
c  Plot the new glitch curve
c
                    do 34 k=1,npts(j)
                      if(lamflg .eq. 0) x(j,k) = hc/x(j,k)
                      xx(k) = x(j,k)
  34                  yy(k) = y(j,k)
 					call graph(ncurv,npts(j),xx,yy,xmin,xmax,elabl,
     z                   ymin,ymax,philbl,hpri)
 35                 continue
c
c  End all the h,k,l loops
c
 39               continue  ! end l loops
 40             continue
 44           continue      ! end k loops
 45         continue
 49       continue          ! end h loops
 50     continue
 51   continue
c
c  Finally, print out the table of glitch h,k,l values
c
      write(*,1052)trim(crys_name),hz(1),hz(2),hz(3)
1052  format(/ a,2x,3f3.0)
      write(*,1053)trim(elabl),emin,emax
1053  format(a,1x,'range',3x,f6.0,2x,f6.0)
      write(*,1054)trim(philbl),phii,phif
1054  format(a,1x,'range',3x,f6.2,2x,f6.2)
      write(*,1055)ncurv,max
1055  format(i3,2x,'glitches found with h^2+k^2+l^2 <',f4.0/)
      write(*,1056)
1056  format(' Glitch H',7x,'Glitch H-P')
      do 55 i=1,ncurv
 55   write(*,1060)glitchlist(i,1),glitchlist(i,2),glitchlist(i,3),
     z         glitchlist(i,4),glitchlist(i,5),glitchlist(i,6)
1060  format(3f4.0,3x,3f4.0)
c
c  Close the plot Xwindow and quit
c
      call pgclos
      stop 
      end
c
c     -----------------------  subroutines:  ------------------------
c
      double precision function dot(a,b,e)
      real*8 a(3),b(3),e(3,3),x
      x=0.0d0
      do 10 i=1,3
      do 10 j=1,3
 10   x=x+a(i)*b(j)*e(i,j)
      dot=x
      return
      end
c
c  **************************************
c
      subroutine cross(a,b,c,elow,e0)
      real*8 a(3),b(3),c(3),x,elow(3,3),e,eup,ee,e0
      common /cc1/e(3,3),ee(3,3,3)
      do 10 l=1,3
        x=0.0d0
        do 20 i=1,3
        do 20 j=1,3
        do 20 k=1,3
 20     x=x+e0*a(i)*b(j)*elow(l,k)*ee(i,j,k)
 10   c(l)=x
      return
      end
c
c  ************************************
c
      subroutine lattice(crys_name,t,elow,eup,e0)
c
c  This subroutine finds the crystal lattice parameters in a
c  data file and produces matrices useful for
c  performing vector products in the crystal lattice.
c  It also loads the common block with matrices E and EE.
c
c  elow(3,3) contains dot products of the unit cell vectors on
c    each other
c  e(3,3) is the identity matrix
c  ee(3,3,3) is a unitary antisymmetric matrix used in cross products
c  eup(3,3) is the inverse of elow
c  e0 equals the reciprocal of the unit cell volume
c
c  Makes calls to SETCRYS, LINEQ, and DTMNT
c     (LINEQ and DTMNT are in the MTRX library)
c
      real*8 angle(3,3),alatt(3),elow(3,3),eup(3,3),eupp(3,3),e0,
     z       elow1(3,3),e(3,3),ee(3,3,3),t,s
      integer ierr
      character*20 crys_name
      data s/1.7453293d-2/       !conversion from deg to rad
      common/cc1/e,ee
      call setcrys(crys_name,t,alatt,angle,ierr)
c
c  SETCRYS takes the lattice data associated with CRYS_NAME
c  and loads it into that ALATT and ANGLE matrices
c
      if(ierr .eq. 0)goto 5
      write(*,999)
  999 format(/' Error with SETCRYS')
      return
c
   5  angle(2,1)=angle(1,2)
      angle(3,2)=angle(2,3)
      angle(1,3)=angle(3,1)
      do 10 i=1,3
        angle(i,i)=0.0d0
        do 20 j=1,3
          elow(i,j)=alatt(i)*alatt(j)*dcos(s*angle(i,j))
          do 30 k=1,3
            ee(i,j,k)=0.0d0
 30         if(i.ne.j.and.j.ne.k.and.k.ne.i) ee(i,j,k)=-1.0d0
 20       e(i,j)=0.0d0
 10     e(i,i)=1.0d0
      ee(1,2,3)=1.0d0
      ee(2,3,1)=1.0d0
      ee(3,1,2)=1.0d0
      call lineq(elow,e,eup,3,3,3,elow1,ioops)
      if(ioops .eq. 0) goto 39
      write(*,1000)
1000  format(/' Error in LINEQ')
      call exit
 39   continue
      do 40 i=1,3
      do 40 j=1,3
 40     eupp(i,j)=eup(i,j)
      call dtmnt(eupp,e0)
      e0=dsqrt(e0)
      return
      end
c
c  *******************************************
c
      subroutine setcrys(crys_name,t,alatt,angle,ierr)
c
c  This subroutine loads the matrices ALATT and ANGLE
c  with values appropriate for the crystal identified 
c  by CRYS_NAME
c
c  Only diamond, silicon, and germanium lattices are supported
c  by this simple version
c
      real*8 alatt(3),angle(3,3),t,temp_coef,a0
      character*20 crys_name
      if(crys_name .eq. 'diamond')goto 10
      if(crys_name .eq. 'silicon')goto 20
      if(crys_name .eq. 'germanium')goto 30
      write(*,1000)
1000  format(/ 'Crystal name not recognized')
      ierr = 1
      return
  10  a0 = 3.567d0
      temp_coef = 1.05d-6
      goto 50
  20  a0 = 5.431d0
      temp_coef = 2.6d-6
      goto 50
  30  a0 = 5.658d0
      temp_coef = 5.8d-6
  50  continue
      do 55 i=1,3
        do 55 j=1,3
  55    angle(i,j) = 0.0d0
      angle(2,3) = 90.0d0
      angle(3,1) = 90.0d0
      angle(1,2) = 90.0d0
      do 60 i=1,3
  60    alatt(i) = a0*(1.d0 + temp_coef*(t-293.d0))
      ierr = 0
      return
      end
c
c  ***************************
c
      subroutine scatamp(h,fact)
c
c  A simple routine to identify allowed reflections
c  in a diamond-cubic lattice
c
c  input h(3) is the reciprocal lattice vector
c
c  output fact is 0 for forbidden reflections, 1 for allowed
c
      real*8 h(3),test1,test2
      complex*16 fact
      test1 = abs(dmod((h(1)+h(2)+h(3)),4.d0)) +
     z      abs(dmod(h(1),2.d0)) + abs(dmod(h(2),2.d0)) + 
     z      abs(dmod(h(3),2.d0))
      test2 = abs(dmod(h(1),2.d0)) + abs(dmod(h(2),2.d0)) + 
     z      abs(dmod(h(3),2.d0))
      if(test1 .eq. 0.d0) goto 50
      if(test2 .eq. 3.d0) goto 50
      fact = (0.0d0,0.0d0)
      return
   50 fact = (1.0d0,0.0d0)
      return
      end
c
c  ******************************
c
      subroutine graph(ioptn,n,x,y,xmin,xmax,xlbl,ymin,ymax,ylbl,hkl)
c
c  This subroutine calls routines from the PGPLOT library
c  to generate graphics commands for an Xwindow 
c
c  ioptn is a flag for choosing among plotting options
c  n is the number of (x,y) points
c  x is a REAL*4 array containing abscissa values
c  y is a REAL*4 array containing ordinate values
c  xmin and xmax are the x-limits
c  ymin and ymax are the y-limits
c  xlabl and ylabl are labels for the x and y axes.
c    
c  Just for GLITCH: hkl is a 3-element REAL*8 vector
c     containing h,k,l values for the primary glitch 
c     reflection being plotted
c
      real*4 x(512),y(512),xmin,xmax,ymin,ymax
      real*8 hkl(3)
      integer pgopen,istat,ioptn
      character*25 xlbl,ylbl
      character*8 hklstr
      character*1 curschar
c
c  Make the h,k,l values into a text string, for printing
c
      write(hklstr,1000)int(hkl(1)),int(hkl(2)),int(hkl(3))
1000  format(i2,',',i2,',',i2)
c
c  If this is the first glitch to be plotted, open the plotting 
c     window
c
      if(ioptn .gt. 1) goto 10
      istat = pgopen('/xwindow')
      if(istat .le. 0) stop
      call pgenv(xmin,xmax,ymin,ymax,0,1)
      call pglab(xlbl,ylbl,'')
  10  continue
c
c  Now draw the glitch curve
c
      call pgline(n,x,y)
c
c  Finally, label the line with h,k,l value
c
      call pgptxt(x(n),y(n),90.0,1.0,hklstr)
      return
      end  