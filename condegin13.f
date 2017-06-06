* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   CALCULATION OF ELECTRON TRANSPORT COEFFICIENTS IN MAGNETIC FIELDS  *
*           for the case  of very strong electron degeneracy           *
*      (thermal averaging is completely neglected in this version).    *
*   Whenever need arbitrary degeneracy, use version `CONDUCT' instead! *
* Last revision: 01.02.2013.                                           *
*   For theoretical background and references see:                     *
*           http://www.ioffe.rssi.ru/astro/conduct/                    *
* Difference from generic CONDEGEN - allowance for nuclear form factor *
*         in the steplike nuclear profile approximation                *
*      Remarks and suggestions are welcome. Please send them to        *
*       Alexander Potekhin <palex@astro.ioffe.ru>                      *
*   This code stems from condegin08.f (v.06.10.2008)                   *
* Differences: - ion thermal conduction (CONDIN) is switched off,      *
*     because its implemented approximation might cause problems;      *
* - exponential suppression of conductivities at low T is switched off *
*     following the criticism by Chugunov (2012, Ast. Lett. 38, 25).   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **

*   ----------------------   MAIN block   ---------------------------  *
*       This is auxiliary MAIN program for input/output purposes.      *
*                You can change it or write your own.                  *
*                     (It may be commented-out.)                       *
*     Calculations are performed in the subroutine CONDEGIN (below).   *
*       Most of internal quantities are in the relativistic units      *
*         (\hbar = m_e = c = 1; default throughout the file)           *
*   -----------------------------------------------------------------  *
C%C      implicit double precision (A-H), double precision (O-Z)
C%C      save
C%C      character KEY
C%C      data AUM/1822.9/,AUD/15819.4/,DRIP/4.3d11/
C%C* NOTATIONS:
C%C*        AUM - atomic mass unit divided by the electron mass
C%C*        AUD - relativistic unit of density in g/cm^3
C%C*        DRIP - neutron drip density in g/cm^3
C%C      data UNISIG/7.763d20/,UNIKAP/2.778d15/
C%C* Rel.units of SIGMA and KAPPA expressed in CGS.
C%C      data BOHR/137.036/,PI/3.14159265/
C%C      write(*,'('' Charge and atomic mass of ions (Z,A): ''$)')
C%C      read*,Zion,CMI
C%C   20 continue
C%C      write(*,'('' Impurity parameter (effective Z): ''$)')
C%C      read*,Zimp
C%C   25 continue
C%C      write(*,'('' Magnetic field B12 (in 10^12 Gauss): ''$)')
C%C      read*,B12
C%C      B=B12/44.14 ! B is the magnetic field in relativistic units
C%C   50 continue
C%C      write(*,'('' Temperature T6 (in 10^6 K): ''$)')
C%C      read*,T6
C%C      TEMP=T6/5930. ! Temperature in mc^2
C%C  100 continue
C%C      write(*,'('' Density of ions rho (in g/cm^3): ''$)')
C%C      read*,RHO
C%C      RHO6=RHO/1.D6
C%C      if (RHO.gt.DRIP) then
C%C   10    write(*,'('' Total number of nucleons per nucleus: ''$)')
C%C         read*,CMI1
C%C        if (CMI1.lt.CMI) then
C%C           print*,'Should be >',CMI,' !  Try again.'
C%C           goto 10
C%C        endif
C%C      else
C%C         CMI1=CMI
C%C      endif
C%C      DENSI=RHO/(AUD*AUM*CMI1) ! number density of ions
C%C*  Call for the central subroutine which calculates the transport       *
C%C*  coefficients SIGMA,KAPPA,Q (in Relativistic units):
C%C      call CONDEGIN(TEMP,DENSI,B,Zion,CMI,CMI1,Zimp,
C%C     & RSIGMA,RTSIGMA,RHSIGMA,RKAPPA,RTKAPPA,RHKAPPA)
C%C* If you want to allow for ion thermal conduction in the approximation
C%C*  of Chugunov & Haensel (2007), then uncomment the next line:
C%CC      call CONDIN(TEMP,DENSI,Zion,CMI,CMI1,RKAPi)
C%C      RKAPPA=RKAPPA+RKAPi
C%C      RTKAPPA=RTKAPPA+RKAPi
C%C*   -------  CONVERSION TO ORDINARY PHYSICAL (CGSE) UNITS:   -------- *
C%C      SIGMA=RSIGMA*UNISIG ! SIGMA in s^{-1}
C%C      CKAPPA=RKAPPA*UNIKAP ! KAPPA in erg/(K cm s)
C%C      SIGMAT=RTSIGMA*UNISIG
C%C      CKAPPAT=RTKAPPA*UNIKAP
C%C      SIGMAH=RHSIGMA*UNISIG
C%C      CKAPPAH=RHKAPPA*UNIKAP
C%C*   ----------------------   OUTPUT:   ------------------------------ *
C%C      T=T6*1.D6
C%C      write(*,113)
C%C      write(*,111) RHO,T,SIGMA,CKAPPA,
C%C     *   SIGMAT,CKAPPAT,
C%C     *   SIGMAH,CKAPPAH
C%C      write(*,'('' New density? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 100
C%C      write(*,'('' New temperature? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 50
C%C      write(*,'('' New B? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 25
C%C      write(*,'('' New Zimp? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 20
C%C      stop
C%C  111 format(1PE12.4,E10.3,3(2X,2E10.2))
C%C  112 format('  Output of ''condegen'':',
C%C     *' rho in g/cc, sigma in 1/s, kappa in erg/K/cm/s'/
C%C     /F10.4,2F10.2,'  (values of B12, Z, A)')
C%C  113 format(29X,'longitudinal',
C%C     *'          transverse            off-diagonal'/
C%C     *'      rho       T   ',3('       sigma     kappa'))
C%C      end

* * * * * * * * * * * * *   Block CONDEGIN  * * * * * * * * * * * * * **
*  This subroutine calculates the electron electric conductivity tensor*
*                   in strongly degenerate matter                      *
*   --------------------------------------------------  Version 12.11.07
      subroutine CONDEGIN(TEMP,DENSI,B,Zion,CMI,CMI1,Zimp,
     & RSIGMA,RTSIGMA,RHSIGMA,RKAPPA,RTKAPPA,RHKAPPA)
* Input: TEMP - temperature, DENSI - number density of ions, B - magn.f.
*        Zion and CMI - ion charge and mass numbers,
*        CMI1 - number of nucleons per nucleus, Zimp -impurity parameter
*          (eff.Z) - that is, Z_{imp}^2 = < n_j (Z-Z_j)^2 > / n,
*         where Z_j, n_j are partial charges and densities of impurities
* Output: RSIGMA, RTSIGMA, RHSIGMA - 
*         longitudinal, transverse and off-diagonal conductivities
*         (all in the relativistic units: \hbar = m_e = c = 1)
      implicit double precision (A-H), double precision (O-Z)
      save
      data PI/3.14159265d0/
      data BOHR/137.036/
*   -------------------   RESTRICTIONS:   --------------------------   *
      if (TEMP.le.0..or.DENSI.le.0..or.B.lt.0..or.Zion.le.0.
     & .or.CMI.le.0.) stop'CONDEGIN: Non-positive input parameter'
      if (CMI1.lt.CMI) stop'CONDEGIN: Incorrect CMI1'
      if (Zion.lt..5) stop'CONDEGIN: Too small ion charge'
      if (CMI.lt.1.) stop'CONDEGIN: Too small ion mass'
      if (DENSI.gt.1.d6) stop'CONDEGIN: Too high density'
*   -----------------   PLASMA PARAMETERS   ------------------------   *
      DENS=DENSI*Zion ! DENS - number density of electrons
      SPHERION=(.75/PI/DENSI)**.3333333 ! Ion sphere radius
      GAMMA=Zion**2/BOHR/TEMP/SPHERION ! Ion coupling parameter
      XSR=(3.*PI**2*DENS)**.3333333 ! special relativity parameter 
*   XSR equals the non-magnetic Fermi momentum
      EF0=dsqrt(1.+XSR**2) ! non-mag.Fermi energy
      Q2e0=4./PI/BOHR*XSR*EF0 ! non-magn.e-screening at strong degen.
      CST=XSR**3/.75/PI*Zion/BOHR**2 ! =4\pi n_i(Ze^2)^2
      if (CMI.eq.CMI1) then ! outer envelope
         xnuc=.00155*(CMI/Zion)**.33333*XSR ! i.e. r_nuc=1.15 A^{1/3} fm
      else ! inner envelope
         xnuc=.00247*XSR ! i.e. r_nuc=1.83 Z^{1/3} fm (Itoh&Kohyama'83)
      endif ! xnuc=r_nuc/a_i - nucleus size parameter
*   -------------------  Chemical potential   ----------------------   *
      if (XSR**2.gt.4.d2*B) then ! non-quantizing case
         PCL=XSR
         Q2e=Q2e0
      else ! quantizing magn.field
         PM0=XSR**3/1.5/B ! p_F in strongly quantizing field
        if (PM0**2.le.2.*B) then ! strongly quantizing case
           PCL=PM0
           Q2e=B/PI/BOHR/PCL ! e-screening
        else ! weakly quantizing case - find p_F by iteration
           Pmax=PM0
           Pmin=XSR/2.
   22      PCL=(Pmax+Pmin)/2.
           NL=PCL**2/2./B
           SN=PCL
           SM=1./PCL
          do N=1,NL
             PN=dsqrt(PCL**2-2.*B*N) ! =p_n
             SN=SN+2.*PN
             SM=SM+2./PN
          enddo
           D=SN*B/2./PI**2 ! estimate n_e
          if (D.lt.DENS) then ! increase PCL
             Pmin=PCL
          else ! decrease PCL
             Pmax=PCL
          endif
          if (dabs(D-DENS).gt.1.d-4*DENS) goto 22 ! next iteration
           Q2e=B/PI/BOHR*SM
        endif
      endif
*   -------------------   Relaxation times    ----------------------   *
      call COULIN(PCL,XSR,GAMMA,B,Zion,CMI,Q2e,xnuc,
     *   CLeff,CLlong,CLtran,SN,THtoEL)
      E=dsqrt(1.+PCL**2) ! magn.Fermi energy
      TAU=PCL**3/E/4./PI/DENSI/(Zion/BOHR)**2/CLlong/SN
      GYROM=B/E ! gyrofrequency
      TAUt0=PCL**3/E/CST/CLtran*SN
      if (Zimp.gt.0.) call COUL99I(PCL,XSR,GAMMA,B,Q2e, ! incl.impurity
     *   CLeffI,CLlongI,CLtranI,SN)
      TAUlong=TAU*CLlong*Zion**2/(CLlong*Zion**2+CLlongI*Zimp**2)
      TAUt=TAUt0*CLtran*Zion**2/(CLtran*Zion**2+CLtranI*Zimp**2)
      TAUtran=TAUt/(1.+(TAUt*GYROM)**2)
      TAUhall=TAUt*GYROM*TAUtran
* Modification of thermal conductivity: inclusion of THtoEL (21.11.99)
      TAUlongT=TAU*CLlong*Zion**2/
     /  (CLlong*Zion**2*THtoEL+CLlongI*Zimp**2)
      TAUtT=TAUt0*CLtran*Zion**2/(CLtran*Zion**2*THtoEL+CLtranI*Zimp**2)
      TAUtranT=TAUtT/(1.+(TAUt*GYROM)**2)
      TAUhallT=TAUtT*GYROM*TAUtranT
*   ----------------------------------------------------------------   *
*   Longitudinal transport coefficients:
      C=SN*PCL**3/3./PI**2/E/BOHR ! common factor
      RSIGMA=C*TAUlong
      RTSIGMA=C*TAUtran
      RHSIGMA=C*TAUhall
* Find thermal conductivity from the Wiedemann-Franz law:
      CTH=C*PI**2*TEMP/3.*BOHR
      if (B.eq.0.) then ! corrected 12.11.07
         call TAUEESY(XSR,TEMP,TAUEE) ! eff.e-e relax.time
         EECOR=TAUEE/(TAUlongT+TAUEE)
      else
         EECOR=1.
      endif
      RKAPPA=CTH*TAUlongT*EECOR
      RTKAPPA=CTH*TAUtranT*EECOR
      RHKAPPA=CTH*TAUhallT*EECOR
      return
      end

*  ================   EFFECTIVE COULOMB LOGARITHM  ==================  *
      subroutine COULIN(PCL,XSR,GAMMA,B,Zion,CMI,Q2e,xnuc,
     *   CLeff,CLlong,CLtran,SN,THtoEL)
*                                                       Version 24.02.00
*                                                     Corrected 23.05.07
*   Difference from COUL99 - inclusion of xnuc
*   Input: PCL - non-magnetic electron momentum \equiv \sqrt(E^2-1),
*          Q2e - squared electron screening wavenumber(IN RELATIV.UNITS)
*          XSR = p_F/mc - relativity (density) parameter,
*          GAMMA - Coulomb coupling parameter of ions,
*          B - magnetic field,
*          Zion - mean charge of the ion,
*          CMI - mean atomic weight,
*          xnuc - ion radius divided by Wigner-Seitz cell radius
*   Output: CLlong, CLtran - eff.Coulomb log.,
*           SN = N_e(E)/N_0(E) = (3/2)(eB\hbar/c)\sum_{ns} p_n/p_0^3
      implicit double precision (A-H), double precision (O-Z)
      save
      data Uminus1/2.78/,Uminus2/12.973/,AUM/1822.9/,BOHR/137.036/
* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
*        Uminus1,Uminus2 - dimensionless frequency moments of phonons
*        AUM - atomic mass unit divided by the electron mass
*        BOHR - radius of the first Bohr orbit in the rel.units
      data PI/3.14159265/
*   ----------------------   Preliminaries   -----------------------   *
      DENS=XSR**3/3./PI**2 ! number density of electrons
      DENSI=DENS/Zion ! number density of ions (rel.)
      SPHERION=(.75/PI/DENSI)**.3333333 ! Ion sphere radius
      Q2icl=3.*GAMMA/SPHERION**2 ! squared Debye screening momentum
      ECL=sqrt(1.+PCL**2) ! Energy
      VCL=PCL/ECL ! Velocity
      PM2=(2.*PCL)**2 ! squared max.momentum transfer
      TRP=Zion/GAMMA*sqrt(CMI*AUM*SPHERION/3./BOHR) ! =T/T_p
      BORNCOR=VCL*Zion*PI/BOHR ! first non-Born correction
*   ---------------------   Non-magnetic fit   ---------------------   *
      C=(1.+.06*GAMMA)*dexp(-dsqrt(GAMMA))
      Q2s=(Q2icl*C+Q2e)*dexp(-BORNCOR) ! eff.scr.wavenumber in rel.un.
      XS=Q2s/PM2 ! eff.screening param.
      R2W=Uminus2/Q2icl*(1.+.3333*BORNCOR)
      XW=R2W*PM2 ! eff. Debye-Waller param.
** Modification WITH FINITE SIZES OF NUCLEI; xnuc=r_{nuc}/a_i
      XW1=14.7327*xnuc**2 ! =4(9\pi/4)^{2/3} x_{nucl}^2 =coeff.at q^2
      XW1=XW1*(1.+.3333*BORNCOR)*(1.+Zion/13.*dsqrt(xnuc))
      CL=COULAN2(XS,XW,VCL,XW1)
      A0=1.683*sqrt(PCL/CMI/Zion) ! zero-vibr.param.(Baiko&Yakovlev95)
      VIBRCOR=exp(-A0/4.*Uminus1*exp(-9.1*TRP)) ! corr.for zero-vibr.
      T0=.19/Zion**.16667 ! critical T/T_p parameter
      G0=TRP/sqrt(TRP**2+T0**2)*(1.+(Zion/125.)**2) ! mod.10.01.99
      GW=G0*VIBRCOR
      CLeff=CL*GW ! 1st FIT (for non-magnetic electrical conductivity)
      G2=TRP/sqrt(.0081+TRP**2)**3
      THtoEL=1.+G2/G0*(1.+BORNCOR*VCL**3)*.0105*(1.-1./Zion)*
     *  (1.d0+xnuc**2*dsqrt(2.d0*Zion))
Commented-out is the obsolete "normal-processes" correction:
C      TRU=TRP*3.*VCL*BOHR/Zion**.3333333 ! T/Tu
C      if (TRU.lt.20.) then ! correction for dying-out umklapp processes
C         CLhigh=CLeff
C         EU=dexp(-1.d0/TRU)
C* CLlowK, CLlowS corrected on 23.05.07 following A.Chugunov's remark
C         CLlowK=15.918*sqrt(XSR/CMI/Zion)*TRP**3 ! low-T lim.for kappa
C         CLlowS=CLlowK/VCL/BOHR*PI**2/.75*TRP**2 ! low-T lim.for sigma
C         CLeff=CLhigh*EU+CLlowS*(1.-EU)
C         THtoEL=(CLhigh*THtoEL*EU+CLlowK*(1.-EU))/CLeff
C      endif
      if (PCL**2.gt.4.d2*B) then ! Non-magnetic case
         CLlong=CLeff
         CLtran=CLeff
         SN=1.d0
         goto 50
      endif
*   -----------------------   Magnetic fit   -----------------------   *
      ENU=PCL**2/2.d0/B
      NL=ENU
      SN=0.
      do N=0,NL
         PB=dsqrt(ENU-N) ! =p_n/sqrt(2b)
         SN=SN+PB
        if (N.ne.0) SN=SN+PB
      enddo
      SN=SN*1.5d0*B*dsqrt(2.d0*B)/PCL**3
      if (ENU.le.1.d0) then ! Exact calculation     
         Xis=Q2s/2./B ! Screening parameter, scaled magnetically
         ZETA=R2W*2.*B ! magn.scaled exponential coefficient
         Xi=2.*PCL**2/B
         Xsum=Xi+Xis
         Q2M=(EXPINT(Xsum,1)-
     -     dexp(-ZETA*Xi)*EXPINT((1.+ZETA)*Xsum,1))/Xsum
         CLlong=(PCL*VCL/B)**2*Q2M/1.5*GW
         QtranM=(1.+Xsum)*EXPINT(Xsum,0)-1.-dexp(-ZETA*Xi)*
     *      ((1.+(1.+ZETA)*Xsum)*EXPINT((1.+ZETA)*Xsum,0)-1.)
         QtranP=(1.+Xis)*EXPINT(Xis,0)-1.-
     -   ((1.+(1.+ZETA)*Xis)*EXPINT((1.+ZETA)*Xis,0)-1.)
         Q=(ECL**2*QtranP+QtranM)*B/PCL**2 ! Q(E,b)
         CLtran=.375*Q/ECL**2*GW
      else
*   Preliminaries:
         DNU=ENU-NL
         XS1=(dsqrt(XS)+1./(2.+XW/2.))**2
         PN=dsqrt(2.*B*DNU)
         SQB=dsqrt(B)
         X=dmax1(PN/SQB,1.d-10)
*   Longitudinal:
        if (XW.lt..01) then
           EXW=1.
        elseif (XW.gt.50.) then
           EXW=1./XW
        else
           EXW=(1.d0-dexp(-XW))/XW
        endif
         A1=(30.-15.*EXW-(15.-6.*EXW)*VCL**2)/
     /     (30.-10.*EXW-(20.-5.*EXW)*VCL**2)
         Q1=.25*VCL**2/(1.-.667*VCL**2)
         DLT=SQB/PCL*(A1/X-sqrt(X)*(1.5-.5*EXW+Q1)+
     +     (1.-EXW+.75*VCL**2)/(1.+VCL**2)*(X-sqrt(X))/NL)
         Y1=1./(1.+DLT)
         CL0=dlog(1.d0+1.d0/XS1)
         P2=CL0*(.07+.2*EXW)
         Y2=1.5*CL0*(X**3-X/3.)/(NL+.75/(1.+2.*B)**2*X**2)+P2*X
         PY=1.+.06*CL0**2/NL**2
         DT=dsqrt(PY*Y1**2+Y2**2) ! ratio of relax.times
         CLlong=CLeff/DT
*   Transverse:
         DB=1./(1.+.5/B)
         CL1=XS1*CL0
         P1=.8*(1.+CL1)+.2*CL0
         P2=1.42-.1*DB+sqrt(CL1)/3.
         P3=(.68-.13*DB)*CL1**.165
         P4=(.52-.1*DB)*sqrt(sqrt(CL1))
         DLT=SQB/PCL*(P1/X**2*SQB/PCL+P3*alog(NL+0.)/X-
     -     (P2+P4*alog(NL+0.))*sqrt(X))
         CLtran=CLeff*(1.+DLT)
      endif
   50 return
      end

      subroutine COUL99I(PCL,XSR,GAMMA,B,Q2e,
     *   CLeff,CLlong,CLtran,SN) ! IMPURITY
*                                                       Version 29.10.99
*  This is a simplified version of COUL99 for el.-impurity scattering
*   Input: XSR = p_F/mc - relativity (density) parameter,
*          PCL - non-magnetic electron momentum \equiv \sqrt(E^2-1),
*          GAMMA - Coulomb coupling parameter of ions,
*          B - magnetic field,
*          Q2e - squared electron screening wavenumber
*    (ALL IN THE RELATIVISTIC UNITS) 
*   Output: CLlong, CLtran - eff.Coulomb log.,
*           SN = N_e(E)/N_0(E) = (3/2)(eB\hbar/c)\sum_{ns} p_n/p_0^3
      implicit double precision (A-H), double precision (O-Z)
      save
* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
*        BOHR - radius of the first Bohr orbit in the rel.units
      data PI/3.14159265/,XW/1.d99/ ! XW=infinity
*   ----------------------   Preliminaries   -----------------------   *
      DENS=XSR**3/3./PI**2 ! number density of electrons
      ECL=sqrt(1.+PCL**2) ! Energy
      VCL=PCL/ECL ! Velocity
      PM2=(2.*PCL)**2 ! squared max.momentum transfer
*   ---------------------   Non-magnetic fit   ---------------------   *
      C=(1.+.06*GAMMA)*dexp(-dsqrt(GAMMA))
      Q2s=Q2e ! eff.scr.wavenumber in rel.un.
      XS=Q2s/PM2 ! eff.screening param.
      CLeff=COULAN2(XS,XW,VCL,0.d0) ! 1st FIT (non-magn.el.conductivity)
      if (PCL**2.gt.4.d2*B) then ! Non-magnetic case
         CLlong=CLeff
         CLtran=CLeff
         SN=1.d0
         goto 50
      endif
*   -----------------------   Magnetic fit   -----------------------   *
      ENU=PCL**2/2.d0/B
      NL=ENU
      SN=0.
      do N=0,NL
         PB=dsqrt(ENU-N) ! =p_n/sqrt(2b)
         SN=SN+PB
        if (N.ne.0) SN=SN+PB
      enddo
      SN=SN*1.5d0*B*dsqrt(2.d0*B)/PCL**3
      if (ENU.le.1.d0) then ! Exact calculation     
         Xis=Q2s/2./B ! Screening parameter, scaled magnetically
         Xi=2.*PCL**2/B
         Xsum=Xi+Xis
         Q2M=EXPINT(Xsum,1)
         CLlong=(PCL*VCL/B)**2*Q2M/1.5
         QtranM=(1.+Xsum)*EXPINT(Xsum,0)-1.
         QtranP=(1.+Xis)*EXPINT(Xis,0)-1.
         Q=(ECL**2*QtranP+QtranM)*B/PCL**2 ! Q(E,b)
         CLtran=.375*Q/ECL**2
      else
*   Preliminaries:
         DNU=ENU-NL
         XS1=(dsqrt(XS)+1./(2.+XW/2.))**2
         PN=dsqrt(2.*B*DNU)
         SQB=dsqrt(B)
         X=dmax1(PN/SQB,1.d-10)
*   Longitudinal:
        if (XW.lt..01) then
           EXW=1.
        elseif (XW.gt.50.) then
           EXW=1./XW
        else
           EXW=(1.d0-dexp(-XW))/XW
        endif
         A1=(30.-15.*EXW-(15.-6.*EXW)*VCL**2)/
     /     (30.-10.*EXW-(20.-5.*EXW)*VCL**2)
         Q1=.25*VCL**2/(1.-.667*VCL**2)
         DLT=SQB/PCL*(A1/X-sqrt(X)*(1.5-.5*EXW+Q1)+
     +     (1.-EXW+.75*VCL**2)/(1.+VCL**2)*(X-sqrt(X))/NL)
         Y1=1./(1.+DLT)
         CL0=dlog(1.d0+1.d0/XS1)
         P2=CL0*(.07+.2*EXW)
         Y2=1.5*CL0*(X**3-X/3.)/(NL+.75/(1.+2.*B)**2*X**2)+P2*X
         PY=1.+.06*CL0**2/NL**2
         DT=dsqrt(PY*Y1**2+Y2**2) ! ratio of relax.times
         CLlong=CLeff/DT
*   Transverse:
         DB=1./(1.+.5/B)
         CL1=XS1*CL0
         P1=.8*(1.+CL1)+.2*CL0
         P2=1.42-.1*DB+sqrt(CL1)/3.
         P3=(.68-.13*DB)*CL1**.165
         P4=(.52-.1*DB)*sqrt(sqrt(CL1))
         DLT=SQB/PCL*(P1/X**2*SQB/PCL+P3*alog(NL+0.)/X-
     -     (P2+P4*alog(NL+0.))*sqrt(X))
         CLtran=CLeff*(1.+DLT)
      endif
   50 return
      end

      function COULAN2(XS,XW0,V,XW1)
*  ------   Analytic expression for Coulomb logarithm - Version 23.05.00
*   XS=(q_s/2p)^2, where p - momentum, q_s - eff.scr.momentum
*   XW0=u_{-2} (2p/\hbar q_D)^2, where u_{-2}=13, q_D^2=3\Gamma/a_i^2
*   V=p/(mc)
*   XW1=s1*(2p/\hbar)^2, s1 \approx r_{nuc}^2
      implicit double precision (A-H), double precision (O-Z)
      save
      data EPS/1.d-2/,EPS1/1.D-3/,EULER/0.5772156649 d0/
      if(XS.lt.0..or.XW0.lt.0..or.V.lt.0..or.XW1.lt.0.)stop'COULAN2'
      do I=0,1
        if (I.eq.0) then
           XW=XW0+XW1
           B=XS*XW
        else ! to do the 2nd term
          XW=XW1
          B=XS*XW
        endif
        if (I.eq.0.or.KEY.eq.2) then ! 23.05.00: for KEY=2 re-check
* Check applicability of asymptotes:
          if (XW.lt.EPS) then
             KEY=1
             goto 50
          endif
          if (XW.gt.1./EPS.and.B.gt.1./EPS) then
             KEY=2
          elseif (XS.lt.EPS1.and.B.lt.EPS1/(1.+XW)) then
             KEY=3
          else
             KEY=4
          endif
        endif
   50   continue
         EA=dexp(-XW)
         E1=1.-EA
        if (KEY.ne.1) E2=(XW-E1)/XW
        if (KEY.eq.1) then
           CL0=dlog((XS+1.)/XS)
           CL1=.5*XW*(2.-1./(XS+1.)-2.*XS*CL0)
           CL2=.5*XW*(1.5-3.*XS-1./(XS+1.)+3.*XS**2*CL0)
        elseif (KEY.eq.2) then
           CL0=dlog(1.d0+1./XS)
           CL1=(CL0-1.d0/(1.+XS))/2.
           CL2=(2.*XS+1.)/(2.*XS+2.)-XS*CL0
        elseif (KEY.eq.3) then
           CL1=.5*(EA*EXPINT(XW,0)+dlog(XW)+EULER)
           CL2=.5*E2
        elseif (KEY.eq.4) then
           CL0=dlog((XS+1.)/XS)
           EL=EXPINT(B,0)-EXPINT(B+XW,0)*EA
           CL1=.5*(CL0+XS/(XS+1.)*E1-(1.+B)*EL)
           CL2=.5*(E2-XS*XS/(1.+XS)*E1-2.*XS*CL0+XS*(2.+B)*EL)
        else
           stop'COULAN2:invalid KEY'
        endif
        if (I.eq.0) then ! 1st term calculated
           COULAN2=CL1-V**2*CL2
          if (XW1.lt.EPS1) return ! don't calculate the 2nd term
        else ! 2nd term calculated
           COULAN2=COULAN2-(CL1-V**2*CL2)
        endif
      enddo
      return
      end
*   ==================  AUXILIARY SUBROUTINE  ================   *
      function EXPINT(XI,L) ! = e^XI E_{L+1}(XI)
      implicit double precision (A-H), double precision (O-Z)
      save
      data GAMMA/.5772156649D0/,Nrep/21/
      if (XI.ge.1.) then ! continued fraction
         CL=L
         CI=Nrep
         C=0.
         do 11 I=Nrep,1,-1
         C=CI/(XI+C)
         C=(CL+CI)/(1.+C)
   11    CI=CI-1.
         Q0=1./(XI+C)
      else ! power series
         PSI=-GAMMA
         do 21 K=1,L
   21    PSI=PSI+1./K ! Psi(L+1)
         Q0=0.
         CMX=1. ! (-XI)^M/M!
         CL=L
         CM=-1.
         do 22 M=0,Nrep
         CM=CM+1.
        if (M.ne.0) CMX=-CMX*XI/CM ! (-XI)^M/M!
        if (M.ne.L) then ! Actually DQ=-deltaQ0
           DQ=CMX/(CM-CL)
        else
           DQ=CMX*(dlog(XI+1.d-20)-PSI)
        endif
         Q0=Q0-DQ
   22    continue
         Q0=exp(XI)*Q0
      endif
   50 continue
      EXPINT=Q0
      return
      end

*   ----------------------------------------------------

      subroutine TAUEESY(X,TEMP,TAUEE)
      implicit double precision (A-H), double precision (O-Z)
* Relaxation time of electron-electron collisions
*   according to Shternin & Yakovlev (2006),
*   corrected in the nondegenerate regime so as to match Lampe (1968)
* X - rel.parameter, DENS - electron density, TEMP-temperature [rel.un.]
*                                                      Version 17.07.06
      E=dsqrt(1.+X**2)
      V=X/E
      Y=.0963913/TEMP*X*dsqrt(V) ! .096391=2\sqrt{\alpha/\pi}
      C1=.123636+.016234*V**2
      C2=.0762+.05714*V**4
      A=12.2+25.2*V**3
      C=A*dexp(C1/C2)
      YV=Y*V
      CIL=dlog(1.d0+128.56/(37.1*Y+10.83*Y**2+Y**3))*
     * (.1587-.02538/(1.+.0435*Y))/V
      CIT=V**3*dlog(1.d0+C/(A*YV+YV**2))*
     *  (2.404/C+(C2-2.404/C)/(1.+.1*YV))
      CILT=V*dlog(1.d0+C/(A*Y+10.83*YV**2+YV**(8.d0/3.d0)))*
     *  (18.52*V**2/C+(C2-18.52*V**2/C)/(1.+.1558*Y**(1.d0-.75d0*V)))
      FI=CIL+CIT+CILT
      FREQ=.00021381*X*Y*dsqrt(V)*FI ! .00021381=6\alpha^{3/2}/\pi^{5/2}
* Correction for partial degeneracy:
      if (X.gt..001) then
         THETA=TEMP/(E-1.d0) ! degeneracy
      else
         THETA=2.*TEMP/X**2
      endif
      T=25.*THETA
      TAUEE=(1.+T+.4342*dsqrt(THETA)*T**2)/(1.+T**2)/FREQ
      return
      end

      subroutine CONDIN(TEMP,DENSI,Zion,CMI,CMI1,RKAPi)
* ion thermal conductivities in the inner crust, Chugunov & Haensel'07
* Input: TEMP - temperature, DENSI - number density of ions
*        Zion and CMI - ion charge and mass numbers
*        CMI1 - number of nucleons per nucleus
* Output: RKAPi - ion thermal conduction
* All quantities are in the relativistic units: \hbar = m_e = c = 1.
*                                                       Version 06.10.08
      implicit double precision (A-H), double precision (O-Z)
      save
      data Uminus1/2.8/,Uminus2/13./,AUM/1822.9/,BOHR/137.036/
* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
*        Uminus1,Uminus2 - dimensionless frequency moments of phonons
*        AUM - atomic mass unit divided by the electron mass
*        BOHR - radius of the first Bohr orbit in the rel.units
      parameter(PI=3.14159265d0,EPS=1.d-8,TINY=1.d-99)
*   ----------------------   Preliminaries   -----------------------   *
      DENS=DENSI*Zion ! number density of electrons
      XSR=(3.d0*PI**2*DENS)**.3333333 ! x_r - density parameter
      PCL=XSR ! classical Fermi momentum; equality due to B=0.
      VCL=PCL/dsqrt(1.d0+PCL**2)
      SPHERION=(.75/PI/DENSI)**.3333333 ! Ion sphere radius
      QBZ=PCL*(2./Zion)**.3333333 ! q_{BZ}
      GAMI=Zion**2/(BOHR*SPHERION*TEMP)
      TRP=Zion/GAMI*dsqrt(CMI*AUM*SPHERION/3.d0/BOHR) ! =T/T_pi
      OMPI=TEMP/TRP ! ion plasma frequency
*      if (GAMI.lt.TINY**.2) pause'CONDIN: GAMI is too low: kappa>HUGE'
*   ---------------------- reduced ion heat capacity:
      call HLfit8(1.d0/TRP,F,U,CV,S,1)
*   ---------------------- ion-electron:
      A0=1.683*sqrt(PCL/CMI/Zion) ! zero-vibr.param.(Baiko&Yakovlev95)
      WDW=A0*(.5*Uminus1*exp(-9.1*TRP)+Uminus2*TRP) ! DW factor (B&Y'95)
      if (CMI.eq.CMI1) then ! outer envelope
         xnuc=.00155*(CMI/Zion)**.33333*XSR ! i.e. r_nuc=1.15 A^{1/3} fm
      else ! inner envelope
         xnuc=.00247*XSR ! i.e. r_nuc=1.83 Z^{1/3} fm (Itoh&Kohyama'83)
      endif ! xnuc=r_nuc/a_i - nucleus size parameter
      W=WDW+43.*xnuc**2
      F=.014+.03/(1.d0+dexp(.2/TRP)) ! CH'07, Eq.(45)
      Y=QBZ/(2.*PCL)
      WY2=W*Y**2
      if (W.gt.EPS) then
         EXPW=dexp(-W)
         EXPWY2=dexp(-WY2)
         CL=(EXPINT(WY2,0)*EXPWY2-EXPINT(W,0)*EXPW-
     -     VCL*(EXPWY2-EXPW)/W)/2.
      else ! small-W asymptote
         CL=dlog(1.d0/Y)-VCL*(1.d0-Y**2)/2.
      endif
      FLEN=8.32d5*SPHERION/(1.d0+PCL**2)/CL*F/Zion*
     *  dsqrt((XSR/1.00884)**3*CMI/Zion) ! Eq.(43)
      CS=OMPI/(3.*QBZ) ! sound speed
      RKAPie=CV*DENSI*CS*FLEN/3. ! Eq.(42) of CH'07
*   ---------------------- ion-ion:
      if (TRP.lt.1./dlog(1./TINY)) then
*         pause'CONDI: T is too low:kappa>HUGE'
         RKAPi=RKAPie
      else
         RKAP0=OMPI*DENSI*SPHERION**2
         RKAPii=RKAP0*dsqrt(1./GAMI**5/
     /     dlog(2.d0+.57735/dsqrt(GAMI)**3)**2+
     +     0.16+(GAMI/77.)**2*dexp(.666667/TRP)) ! Eq.(26) of CH'07
*   ---------------------- total ion:
         RKAPi=1./(1./RKAPii+1./RKAPie)
      endif
      return
      end
      
      subroutine HLfit8(eta,F,U,CV,S,LATTICE)
*                                                       Version 15.02.08
* Baiko, Potekhin, & Yakovlev (2001). Stems from HLfit v.20.03.01
* Fit to thermal part of the thermodynamic functions.
* Zero-point lattice quantum energy 1.5u_1\eta INCLUDED (unlike HLfit)
* Input: eta=Tp/T
* Output: F and U (normalized to NkT),
*   CV and S (normalized to Nk) in the HL model for bcc Coulomb lattice
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(EPS=1.d-5)
      if (LATTICE.eq.1) then ! bcc lattice
         CLM=-2.49389 ! 3*ln<\omega/\omega_p>
         U1=.5113875
         ALPHA=.265764
         BETA=.334547
         GAMMA=.932446
         A1=.1839
         A2=.593586
         A3=.0054814
         A4=5.01813d-4
         A6=3.9247d-7
         A8=5.8356d-11
         B0=261.66
         B2=7.07997
         B4=.0409484
         B5=.000397355
         B6=5.11148d-5
         B7=2.19749d-6
         C9=.004757014
         C11=.0047770935
      elseif (LATTICE.eq.2) then ! fcc lattice
         CLM=-2.45373
         U1=.513194
         ALPHA=.257591
         BETA=.365284
         GAMMA=.9167070
         A1=.0
         A2=.532535
         A3=.0
         A4=3.76545d-4
         A6=2.63013d-7
         A8=6.6318d-11
         B0=303.20
         B2=7.7255
         B4=.0439597
         B5=.000114295
         B6=5.63434d-5
         B7=1.36488d-6
         C9=.00492387
         C11=.00437506
      else
         stop'HLfit: unknown lattice type'
      endif
      if (eta.gt.1./EPS) then ! asymptote of Eq.(13) of BPY'01
         U=3./(C11*eta**3)
         F=-U/3.
         CV=4.*U
        goto 50
      elseif (eta.lt.EPS) then ! Eq.(17) of BPY'01
         F=3.*dlog(eta)+CLM-1.5*U1*eta+eta**2/24. 
         U=3.-1.5*U1*eta+eta**2/12.
         CV=3.-eta**2/12.
         goto 50
      endif
      B9=A6*C9
      B11=A8*C11
      UP=1.+A1*eta+A2*eta**2+A3*eta**3+A4*eta**4+A6*eta**6+A8*eta**8
      DN=B0+B2*eta**2+B4*eta**4+B5*eta**5+B6*eta**6+
     +  B7*eta**7+B9*eta**9+B11*eta**11
      EA=dexp(-ALPHA*eta)
      EB=dexp(-BETA*eta)
      EG=dexp(-GAMMA*eta)
      F=dlog(1.d0-EA)+dlog(1.d0-EB)+dlog(1.-EG)-UP/DN ! thermal free energy/NT
      UP1=A1+
     + 2.*A2*eta+3.*A3*eta**2+4.*A4*eta**3+6.*A6*eta**5+8.*A8*eta**7
      UP2=2.*A2+6.*A3*eta+12.*A4*eta**2+30.*A6*eta**4+56.*A8*eta**6
      DN1=2.*B2*eta+4.*B4*eta**3+5.*B5*eta**4+6.*B6*eta**5+
     +  7.*B7*eta**6+9.*B9*eta**8+11.*B11*eta**10.
      DN2=2.*B2+12.*B4*eta**2+20.*B5*eta**3+30.*B6*eta**4+
     +  42.*B7*eta**5+72.*B9*eta**7+110.*B11*eta**9
      U=ALPHA*EA/(1.d0-EA)+BETA*EB/(1.d0-EB)+GAMMA*EG/(1.d0-EG)-
     -  (UP1*DN-DN1*UP)/DN**2 ! int.en./NT/eta
      CV=ALPHA**2*EA/(1.d0-EA)**2+BETA**2*EB/(1.d0-EB)**2+
     +  GAMMA**2*EG/(1.d0-EG)**2+
     +  ((UP2*DN-DN2*UP)*DN-2.*(UP1*DN-DN1*UP)*DN1)/DN**3 ! cV/eta^2
      U=U*eta
      CV=CV*eta**2
   50 continue
      S=U-F
* Add zero-point lattice energy:
      E0=1.5*U1*eta
      U=U+E0
      F=F+E0
      return
      end