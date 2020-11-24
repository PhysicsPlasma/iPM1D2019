module ModuleSEE
    use ModuleParticleBoundary
    use ModuleParticleOne
    Use Constants
    Use ModuleParticleBundle
    implicit none
    TYPE GasTemp
        real(8):: MGas,TGas
    end Type
contains
    Subroutine  Selectron(PB,PBDO)
        Implicit none
        class(ParticleBundle),intent(inout) :: PB
        Type(ParticleBoundaryOne),intent(inout) :: PBDO
        Type(ParticleOne) :: ParticleTemp        
        Integer :: i,j,n3,n4
        Type(GasTemp)::SEEGas=(GasTemp(ElectronMass,3.0*eVtoK))
        Real(8) :: XRandom=0.1,VFactor,Residue,GammaR
        VFactor = 1.d0/PB%VFactor

        !==========================================================================================================
        do i = 1,PBDO%CountMinOne           !×ó¼«°å
            if(PB%Mass==ElectronMass) THEN  !Secondary electron emmision of electron
                PBDO%Gamma = Gamma_Caculate(PBDO%PBLower%PO(i)%Energy(PBDO%PBLower%Mass,PBDO%PBLower%VFactor))
            ELSE                            !Secondary electron emmision of ion
                PBDO%Gamma = 0.1
            end if
            n3 = int(PBDO%Gamma)
            GammaR = PBDO%Gamma-n3
            Call Random_NUMBER(R)
            if (R<GammaR) then
                n3 = n3+1
            end if
            do j = 1, n3
                call ParticleTemp%AccInpInit()
                call ParticleTemp%VelMaxInit(SEEGas%MGas,SEEGas%TGas)
                call ParticleTemp%VelRes(VFactor)
                ParticleTemp%Vx = abs(ParticleTemp%Vx)
                Call Random_NUMBER(R)
                ParticleTemp%X = PBDO%XMin+R*ParticleTemp%Vx
                call PB%AddOne(ParticleTemp)
                Call PBDO%PBLower%AddOne(ParticleTemp)
            end do
            PBDO%SEECountMinOne = PBDO%SEECountMinOne+n3
        end do
        !==========================================================================================================
        do i = 1,PBDO%CountMaxOne           !ÓÒ¼«°å
            if(PB%Mass==ElectronMass) THEN  !Secondary electron emmision of electron
                PBDO%Gamma = Gamma_Caculate(PBDO%PBUpper%PO(i)%Energy(PBDO%PBUpper%Mass,PBDO%PBUpper%VFactor))
            ELSE                            !Secondary electron emmision of ion
                PBDO%Gamma = 0.1
            end if
            n4 = int(PBDO%Gamma)
            GammaR = PBDO%Gamma-n4
            Call Random_NUMBER(R)
            if (R<GammaR) then
                n4 = n4+1
            end if
            do j = 1, n4
                call ParticleTemp%AccInpInit()
                call ParticleTemp%VelMaxInit(SEEGas%MGas,SEEGas%TGas)
                call ParticleTemp%VelRes(VFactor)
                ParticleTemp%Vx = -abs(ParticleTemp%Vx)
                Call Random_NUMBER(R)
                ParticleTemp%X = PBDO%XMax+R*ParticleTemp%Vx
                call PB%AddOne(ParticleTemp)
                Call PBDO%PBUpper%AddOne(ParticleTemp)
            end do
            PBDO%SEECountMaxOne = PBDO%SEECountMaxOne+n4
        end do
        !==========================================================================================================

        return
      end subroutine Selectron

    ! function Gamma_Caculate(Energy)!Reference :: The effect of electron processes on metal walls in magnetized microdischarges
    !     implicit none
    !     Real(8),intent(in) :: Energy
    !     Real(8) :: Gamma_caculate
    !     Real(8) :: Ee,Em,E0,Gamma_m,s
    !     Em = 262.0  !eV
    !     E0 = 150.0  !eV
    !     Ee = Energy/JtoeV
    !     Gamma_m = 1.06
    !     s = 1.35
    !     Gamma_caculate = Gamma_m*s*(Ee/Em)/(s-1+(Ee/Em)**s)
    !     return
    ! end function Gamma_Caculate
  
    function Gamma_Caculate(Energy)!Reference :: Simulations of multipactor-assisted breakdown in radio frequency plasmas
	   implicit none
		Real(8),intent(in) :: Energy
		Real(8) :: Gamma_caculate
		Real(8) :: Ep,Em,Gamma_m
		Em = 400.d0  !eV
		Ep = Energy/JtoeV
		Gamma_m = 2.4d0
		Gamma_caculate = Gamma_m*exp*exp*(Ep/Em)*Dexp(-2.0*sqrt(Ep/Em))
		return
   end function Gamma_Caculate
end module ModuleSEE