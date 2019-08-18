Module ModuleTypeMCC
    Use ModuleParticleOne
    Use ModuleSpecyOne
    !Use GasModule
    Implicit none
                   !  ReactionType is definition of the reactions:
                   !   -1-The injected particle are removed because the generated particle is no longer included in simulations . 
                   !   0-no collision
                   !   100~199 for electron-neutral collision  (100~110 Base : 101-Elastic; 102-Exitation; 103- Ionization; 104-Attachment; 105-Dissociation)
                   !                                    (111-120 Anisotropic Ar: 111-Elastic; 112-Exitation; 113- Ionization;)
                   !   200~299 for  Ion-neutral collision (201~210 non-reactive : 201-isotropic; 202-ChargeExchange; 203-Anisotropic)
                   !                                    (211~220 reactive : 211-original particle; 212-new particle)    
                Type  ReactionOne    ! Define a collosion process.
                     Integer(4) :: Reactant,ReactionType,Resultant
                     Real(8) ::  Threshold
                End Type  ReactionOne
       
                  !  Mcc bundle for all collision process of one particle withe one gas. 
!                       ! Model is the index for the different MCC models:
!                       ! 0-no collision; 1-elelctrons,the gas moculars are stationary; 2-nonreactive ions; 3-reactive ions;
                Type MCCBundle
                    Integer(4) ::  Model,NReaction=0,NSigma=0
                    Real(8) :: EnergyMin,EnergyInterval,EnergyMax
                    Real(8) :: CollisionRatio,SigmaMax
                    Type(SpecyOne),Pointer :: SO=>Null()
                    Type(GasOne),Pointer :: GO=>Null()
                    Type(ReactionOne),Allocatable :: Reaction(:)
                    Real(8),Allocatable ::  Probility(:,:)
                    contains
                    procedure :: Dump=>DumpMCCBundle
                End Type MCCBundle

                    Type MCCParticleOne
                          Integer(4) :: ReactionIndex=0
                          Logical :: ParticleAnnihilation=.False.,ParticleCreation=.False.
                          Type(ParticleOne),Pointer :: POI=>Null()
                          Real(8) :: Miu,MassRatio,Energy,Beta=0.d0
                          Real(8) :: Gx,Gy,Gz,Gper,G
                          Type(ParticleOne) :: POT
                          Integer(4) :: NPONew=0
                          Type(ParticleOneIndex),Allocatable :: PON(:)
                    contains
                          procedure :: Select=>SelectMCCParticleOne
                          procedure :: Updater=>UpdateMCCParticleOne
                          procedure :: VelocityUpdater=>UpdateVelocityMCCParticleOne
                    End Type  MCCParticleOne
           
    contains         
                Subroutine SelectMCCParticleOne(MCPO,PO)
                     Implicit none
                     Class(MCCParticleOne),intent(inout) :: MCPO
                     Type(ParticleOne),Target,intent(in) :: PO
                     MCPO%POI=>PO
                    return  
              End subroutine SelectMCCParticleOne
    
    
              Subroutine UpdateMCCParticleOne(MCPO,SO,GO)
                     Implicit none
                     Class(MCCParticleOne),intent(inout) :: MCPO
                     Type(SpecyOne),intent(in) :: SO
                     Type(GasOne),intent(in) :: GO
                     Real(8) :: Gtemp,G2
                      Call  MCPO%POT%Copy(MCPO%POI)
                      MCPO%Miu=SO%Mass*GO%Mass/(SO%Mass+GO%Mass)
                      MCPO%MassRatio=GO%Mass/(SO%Mass+GO%Mass)
                      Call MCPO%POT%VelMaxInit(GO%Mass,GO%Temperature)
                     Associate (Gx=>MCPO%Gx,Gy=>MCPO%Gy,Gz=>MCPO%Gz,Gper=>MCPO%Gper,G=>MCPO%G,Energy=>MCPO%Energy)
                          Gx=MCPO%POI%Vx-MCPO%POT%Vx
                          Gy=MCPO%POI%Vy-MCPO%POT%Vy
                          Gz=MCPO%POI%Vz-MCPO%POT%Vz
                          Gtemp=Gy*Gy+Gz*Gz
                          Gper=Dsqrt(Gtemp)
                          G2=Gtemp+Gx*Gx
                          G=Dsqrt(G2)
                          Energy=0.5d0*MCPO%Miu*G2/JtoeV
                          MCPO%POT%Ax=MCPO%POT%Ax*SO%Mass/GO%Mass
                          MCPO%POT%Ay=MCPO%POT%Ay*SO%Mass/GO%Mass
                          MCPO%POT%Az=MCPO%POT%Az*SO%Mass/GO%Mass
                      End Associate
                    return  
              End subroutine UpdateMCCParticleOne
              
                 Subroutine UpdateVelocityMCCParticleOne(MCPO,VFactor)
                    Implicit none
                    Class(MCCParticleOne),intent(inout) :: MCPO
                     Real(8),intent(in) :: VFactor
                     MCPO%Gx=VFactor*MCPO%Gx
                     MCPO%Gy=VFactor*MCPO%Gy
                     MCPO%Gz=VFactor*MCPO%Gz
                     MCPO%Gper=VFactor*MCPO%Gper
                     MCPO%G=VFactor*MCPO%G
                     !MCPO%Energy=VFactor*VFactor*MCPO%Energy
                     MCPO%POI%Vx=VFactor*MCPO%POI%Vx
                     MCPO%POI%Vy=VFactor*MCPO%POI%Vy
                     MCPO%POI%Vz=VFactor*MCPO%POI%Vz
                     !MCPO%Energy=VFactor*VFactor*MCPO%Energy
                     return
                 End  subroutine UpdateVelocityMCCParticleOne
                 
             subroutine PostCollisionVelocity(MCPO,CosKai)
               Implicit none
               Type(MCCParticleOne),intent(inout) :: MCPO
               Real(8),external ::  CosKai
               real(8) :: hx,hy,hz
               real(8) :: APhi,CosPhi,SinPhi
               real(8) :: CosTheta,FCosTheta,SinTheta
               real(8) :: Energy1,Energy2
               
               
               Associate (Gx=>MCPO%Gx,Gy=>MCPO%Gy,Gz=>MCPO%Gz,Gper=>MCPO%Gper,G=>MCPO%G,Energy=>MCPO%Energy,&
                   Vx=>MCPO%POI%Vx,Vy=>MCPO%POI%Vy,Vz=>MCPO%POI%Vz,MassRatio=>MCPO%MassRatio)
               !Energy1=Vx*Vx+Vy*Vy+Vz*Vz
                  If(Abs(Gper)<MinReal)  then
                           Call MCPO%POI%VelRanInit(G)
               else
                       CosTheta=CosKai(MCPO%Energy)
                       FcosTheta=1.d0-cosTheta
                       SinTheta=Dsqrt(1.d0-cosTheta*cosTheta)

                       CALL RANDOM_NUMBER(R)
                       APhi=2.d0*Pi*R
                       CosPhi=DCos(APhi)
                       SinPhi=DSin(APhi)
                       hx=gper*CosPhi
                       hy=-(gx*gy*CosPhi+g*gz*SinPhi)/gper
                       hz=-(gx*gz*CosPhi-g*gy*SinPhi)/gper
                       Vx=Vx-MassRatio*(gx*FcosTheta+hx*SinTheta)
                       Vy=Vy-MassRatio*(gy*FcosTheta+hy*SinTheta)
                       Vz=Vz-MassRatio*(gz*FcosTheta+hz*SinTheta)
               End if 
               !Energy2=Vx*Vx+Vy*Vy+Vz*Vz
              End Associate
              return
             end subroutine PostCollisionVelocity
             
             Subroutine DumpMCCBundle(MCB)
               Implicit none
              Class(MCCBundle),intent(inout) :: MCB
              !Integer(4),intent(in) :: Mode 
               Character(len=99) :: Filename
               Real(8) :: NEnergy,Energy
               Integer(4),save :: i,j,Timer=0
               Timer=Timer+1
               Associate(Emin=>MCB%EnergyMin,Eint=>MCB%EnergyInterval,Nr=>MCB%NReaction,Ns=>MCB%NSigma)
                   Write(filename,*) "MCCB",Timer,".dat"
                   Open (10,file=filename)
                        do i=1,Ns
                            Energy=10.d0**(Emin+dble(i-1)*Eint)
                            Write(10,FMt="(*(es21.14,1x))") Energy, (MCB%Probility(j,i),j=1,Nr)
                        ENd do
               ENd Associate
               Close(10)
               return
             End Subroutine DumpMCCBundle
             
         subroutine SelectProbility(MCPO,MCB)
           Implicit None 
           Type(MCCParticleOne),intent(inout) :: MCPO
           Type(MCCBundle),intent(in) :: MCB
           Integer(4) :: i,N
           Real(8) :: Energy,S1,S2,Residue
           Real(8) :: TempProbility(MCB%NReaction)
           !Real(8),Allocatable :: 
           Associate(Energy=>MCPO%Energy,Emin=>MCB%EnergyMin,Eint=>MCB%EnergyInterval,Nr=>MCB%NReaction,Ns=>MCB%NSigma,Probility=>MCB%Probility,EnergyMax=>MCB%EnergyMax)
                TempProbility=0.d0
                N=Ceiling((Log10(Energy)-Emin)/Eint)
                !IF (Energy>16.d0) Then
                !    Write(*,*)"SS"
                ! ENd If   
                If (N<1) Then
                    TempProbility=Probility(1:Nr,1)
                ELse If (N<Ns) Then
                     S1=Dble(N)-(Log10(Energy)-Emin)/Eint
                     S2=1.d0-S1
                     !If(S1<0.d0.or.S1>1.d0.or.S2<0.d0.or.S2>1.d0) Write(*,*) "ERRRRRR S1S2SelectProbility"
                     TempProbility(1:Nr)=Probility(1:Nr,N)*S1+Probility(1:Nr,N+1)*S2
                Else
                     TempProbility(1:Nr)=Probility(1:Nr,Ns)*Energy/EnergyMax
                End If
                
                DO i=1,Nr
                       If (Energy<MCB%Reaction(i)%Threshold) Then
                            TempProbility(i)=0.d0
                       End If
                  End Do
                
                 Call Random_Number(R)
                 MCPO%ReactionIndex=0
                 Do i=1,Nr
                       If (R<TempProbility(i)) Then
                           MCPO%ReactionIndex=i
                           !if (i==3) Then
                           !    Write(*,*)"SS"
                           ! ENd If
                           exit
                       End If
                  End Do
            End Associate  
         End  subroutine SelectProbility
         
      subroutine SelectProbilityReactive(MCPO,MCB)
          Implicit none
          Type(MCCParticleOne),intent(inout) :: MCPO
          Type(MCCBundle),intent(in) :: MCB
                   Call Random_Number(R)
                  MCPO%Beta=MCB%GO%BetaMax*Dsqrt(R)
                  If (MCPO%Beta<=1.d0) then   !  Reactive collisions
                           If(MCPO%Energy<=MCB%Reaction(2)%Threshold) then
                                !  Reactive elastic collision with isotropic scattering
                                MCPO%ReactionIndex=1_4
                           Else
                                Call SelectProbility(MCPO,MCB)
                           End If
                 Else        ! Nonreactive Elastic collision with anisotropic scattering
                          MCPO%ReactionIndex=1_4
                 End If
                 
           return
        end subroutine  SelectProbilityReactive 
   
    End Module ModuleTypeMCC
