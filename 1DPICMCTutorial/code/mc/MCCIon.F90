Module ModuleMCCIon
    Use ModuleTypeMCC
    Use MCCEnergyKai
     Implicit none
     contains
     subroutine SelectCollisionIon(MCPO,SO,GO,RO)      
         Implicit none
         Type(MCCParticleOne),intent(inout) :: MCPO
         Type(SpecyOne),intent(in) :: SO
         Type(GasOne),intent(in) :: GO
         Type(ReactionOne),intent(in) :: RO

         Select case (RO%ReactionType)
                  !Case(0_4)
                  Case(101_4)
                          Call PostCollisionVelocity(MCPO,IsotropicCosKai)
                  Case(102_4)
                          MCPO%POI=MCPO%POT
                          !Call PostVelocityIonElastic(InputParticle,InputGas,ReactiveCosKai)
                          !Call AddParticle(InputParticle%PhaseSpace,OuputParticle(RO%Reactant))
                  Case(103_4)  ! ChargeExchange
                           Call MCPO%POI%Copy(MCPO%POT)
                  Case(110_4)
                           MCPO%ParticleAnnihilation=.True.
                  Case(111_4)
                          Call  PostVelocityIonReactive(MCPO,RO)
                  Case(112_4)
                          Call  PostVelocityIonReactive(MCPO,RO)
                          MCPO%ParticleAnnihilation=.True.
                          MCPO%ParticleCreation=.True.
                          MCPO%NPONew=1
                           Allocate(MCPO%PON(1))
                           MCPO%PON(1)%ParticleOne=MCPO%POT
                           Call MCPO%PON(1)%IndexInit(RO%Resultant)
                          !Call MCCParticleIndex(Ns,RO%Resultant,1_4,ParticleIndex(1))
                  Case(-1_4)
          end select
        return
        contains
           Subroutine PostVelocityIonReactive(MCPO,RO)
               implicit none
               Type(MCCParticleOne),intent(inout) :: MCPO
               Type(ReactionOne),intent(in) :: RO 
               Real(8) :: VFactor,MassRatioA,MassRatioB
               
              VFactor=Dsqrt(1.d0-RO%Threshold/MCPO%Energy)
              MCPO%Energy=MCPO%Energy-RO%Threshold
              Call MCPO%VelocityUpdater(VFactor)
               
               MassRatioA=MCPO%MassRatio
               MassRatioB=1.d0-MassRatioA
               
                Associate (Gx=>MCPO%Gx,Gy=>MCPO%Gy,Gz=>MCPO%Gz,Gper=>MCPO%Gper,G=>MCPO%G,&
                    Vx=>MCPO%POI%Vx,Vy=>MCPO%POI%Vy,Vz=>MCPO%POI%Vz)
                    Vx=MassRatioB*Vx+MassRatioA*MCPO%POT%Vx+MassRatioB*Gx
                    Vy=MassRatioB*Vz+MassRatioA*MCPO%POT%Vy+MassRatioB*Gy
                    Vz=MassRatioB*Vz+MassRatioA*MCPO%POT%Vz+MassRatioB*Gz
                      return
                 ENd Associate 
            end subroutine PostVelocityIonReactive
  end  subroutine SelectCollisionIon  
end Module ModuleMCCIon

!
