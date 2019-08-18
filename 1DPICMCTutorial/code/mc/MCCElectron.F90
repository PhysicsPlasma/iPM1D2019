Module ModuleMCCElectron
    Use ModuleTypeMCC
    Use MCCEnergyKai 
    Implicit none
    contains 
   !  ReactionType is definition of the reactions:
   !   -1-The injected particle are removed because the generated particle is no longer included in simulations . 
   !   0-no collision
   !   1~99 for electrons  (1~10 Isotropic : 1-Elastic; 2-Exitation; 3- Ionization; 4-Attachment; 5-Dissociation)
   !                                    (11-20 Anisotropic Ar: 11-Elastic; 12-Exitation; 13- Ionization;)
    
    
    Subroutine SelectCollisionElectron(MCPO,SO,GO,RO)   
         Implicit none
         Type(MCCParticleOne),intent(inout) :: MCPO
         Type(SpecyOne),intent(in) :: SO
         Type(GasOne),intent(in) :: GO
         Type(ReactionOne),intent(in) :: RO

         Select case(RO%ReactionType)
                  Case(0_4) ! No collisions occur.
                  Case(11_4)  ! Isotropic Elastic collisions.
                          Call PostCollisionVelocity(MCPO,CosKaiVahedi)
                  Case(12_4)  ! Isotropic Excitation collisions occur.
                          Call  ElectronExcitation(MCPO,RO,CosKaiVahedi)
                  Case(13_4)  ! Isotropic Ionization collisions.
                          Call ElectronIonization(MCPO,SO,GO,RO,CosKaiVahedi,CreatEnergyVehadi)
                  Case(14_4)   ! Attachment collisions.
                          Call ElectronAttachment(MCPO,GO,RO)
                  Case(105_4)
                          !Call DissociationElectron(MCPO,Gas,TempParticle,RO%Threshold,IsotropicCosKai,IsotropicEnergy)

              End select
        return
        
      contains
                Subroutine ElectronExcitation(MCPO,RO,CosTheta)
                    Implicit none
                    Type(MCCParticleOne),intent(inout) :: MCPO
                     Type(ReactionOne),intent(in) :: RO
                     Real(8),external :: CosTheta
                    Real(8) :: VFactor
                    VFactor=Dsqrt(1.d0-RO%Threshold/MCPO%Energy)
                    MCPO%Energy=MCPO%Energy-RO%Threshold
                    Call MCPO%VelocityUpdater(VFactor)
                    Call PostCollisionVelocity(MCPO,CosTheta)
                    return
                End subroutine ElectronExcitation
      !         
                Subroutine ElectronIonization(MCPO,SO,GO,RO,CosTheta,CreateEnergy)
               ! This is the  subroutine for mainly the ionization collision of electrons. After the  collision, a new electron  and a new ion
                ! will be created. 
                    Implicit none
                    Type(MCCParticleOne),intent(inout) :: MCPO
                    Type(SpecyOne),intent(in) :: SO
                    Type(GasOne),intent(in) :: GO
                   Type(ReactionOne),intent(in) :: RO
                    Real(8),external :: CosTheta,CreateEnergy
                    Type(MCCParticleOne) :: TempMCPO
                    Type(ParticleOne),Target :: TempPO
                    Real(8) :: VFactor,EnergyTemp,EnergyOld,EnergyNew
                    ! Particle data are temperally stored.
                    TempMCPO=MCPO
                    ! Injected electron's velocity is determinated below.
                    
                    MCPO%ParticleCreation=.True.
                    EnergyTemp=MCPO%Energy-RO%Threshold
                    EnergyOld=CreateEnergy(EnergyTemp)
                    VFactor=Dsqrt(EnergyOld/MCPO%Energy)
                    Call MCPO%VelocityUpdater(VFactor)
                    Call PostCollisionVelocity(MCPO,CosTheta)
                    
                    ! Created electron's velocity is determinated below. Energy is conserved.
                     EnergyNew=EnergyTemp-EnergyOld
                     TempPO=TempMCPO%POI
                     TempMCPO%POI=>TempPO
                     VFactor=Dsqrt(EnergyNew/TempMCPO%Energy)
                     Call TempMCPO%VelocityUpdater(VFactor)
                     Call PostCollisionVelocity(TempMCPO,CosTheta)
                     MCPO%NPONew=2 
                     !Write(*,*)"SS,Ionization"
                     Allocate(MCPO%PON(2))
                     MCPO%PON(1)%ParticleOne=TempMCPO%POI
                     Call MCPO%PON(1)%IndexInit(SO%SpecyIndex)
                     
                     !Call MCPO%PON(2)%ParticleOne%Copy(TempMCPO%POI)
                     !Call MCPO%PON(2)%ParticleOne%VelMaxInit(GO%Mass,GO%Temperature)
                     MCPO%PON(2)%ParticleOne=MCPO%POT
                     Call MCPO%PON(2)%IndexInit(RO%Resultant)
                    return
                End subroutine ElectronIonization
                
                Subroutine ElectronAttachment(MCPO,GO,RO)
                    Type(MCCParticleOne),intent(inout) :: MCPO
                    Type(GasOne),intent(in) :: GO
                    Type(ReactionOne),intent(in) :: RO    
                    MCPO%ParticleAnnihilation=.True.
                    MCPO%ParticleCreation=.True.
                    MCPO%NPONew=1
                    Allocate(MCPO%PON(1))
                    MCPO%PON(1)%ParticleOne=MCPO%POT
                    Call MCPO%PON(1)%IndexInit(RO%Resultant)
                    return
                End subroutine ElectronAttachment
      ENd  subroutine   SelectCollisionElectron       
    End   Module ModuleMCCElectron                
                
      !
      !          Subroutine DissociationElectron(MCPO,Gas,OutputParticle,Threshold,CosTheta,CreateEnergy)
      !              Implicit none
      !              Type(MCCParticleOne),intent(inout) :: MCPO
      !              Type(GasType),intent(in) :: Gas
      !              Type(MCCParticleOne) :: ParticleTemp
      !              Type(ParticleOne) ::  OutputParticle(2)
      !              Real(8),intent(in) :: Threshold
      !              Real(8),external :: CosTheta,CreateEnergy 
      !              Real(8) :: VFactor,EnergyTemp
      !              EnergyTemp=MCPO%Energy-RO%Threshold
      !              VFactor=Dsqrt(EnergyTemp/MCPO%Energy)
      !              Call MCPO%VelocityUpdater(VFactor)
      !              Call PostCollisionVelocity(MCPO,CosKaiVahedi)
      !              OutputParticle(1)=MCPO%PhaseSpace
      !              OutputParticle(2)=MCPO%PhaseSpace  
      !              !Call Maxwellian(Gas,OutputParticle(1)) 
      !              !Call Maxwellian(Gas,OutputParticle(2)) 
      !              return
      !          End subroutine DissociationElectron


    
    

        !Abstract    Interface
        !    Function CosChi(Energy) 
        !        Real(8) :: CosChi,Energy
        !    end Function CosChi
        !    
        !    Function PartitionEnergy(Energy) 
        !        Real(8) :: PartitionEnergy,Energy
        !    end Function PartitionEnergy
        !End Interface          
        !Procedure(CosChi), pointer:: CosChiPointer => NULL()
        !Procedure(PartitionEnergy), pointer:: PartitionEnergyPointer => NULL()
        !
        !Select case(RO%ReactionType)
        !    Case(111)
        !      CosChiPointer => EArCosKai
        !    Case(112)
        !      CosChiPointer => EArCosKai
        !    Case(113)
        !      CosChiPointer => EArCosKai
        !      PartitionEnergyPointer =>EArEnergy
        !    Case DEFAULT
        !      CosChiPointer => IsotropicCosKai
        !      PartitionEnergyPointer => IsotropicEnergy
        !End Select    