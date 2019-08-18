Module ModuleMCCInitialization
     Use ModuleMCCPublic 
     Use  ModuleControlFlow
     
     Implicit none

     !Integer(4),parameter :: NsMCCTemp=4_4
     !Type(MCCParticleOne),Allocatable,private :: TempMCCParticle(:)
     !Type(MCCBundle),Allocatable :: MCCBundleGlobal(:,:)
     
     Contains
            subroutine  MCCBundleInit(CF,SO,GO)
               Implicit none
               Class(ControlFlow),intent(in) :: CF
               Type(SpecyOne),intent(in) :: SO(0:CF%Ns)
               Type(GasOne),intent(in) :: GO(CF%Ng)

               Associate(Ns=>CF%Ns,Ng=>CF%Ng)
                   Allocate(MCCBundleGlobal(0:Ns,Ng))
                   Call MCCBundleElelctronInit(CF,SO(0),GO(1:Ng),MCCBundleGlobal(0,1:Ng))
                   Call MCCBundleIonInit(CF,SO(1:Ns),GO(1:Ng),MCCBundleGlobal(1:Ns,1:Ng))
               End Associate
               Call MCCBundleGlobal(0,1)%Dump
               Return
            End subroutine  MCCBundleInit
    
    
   subroutine  MCCBundleElelctronInit(CF,SO,GO,MCCB)
             Implicit none
             Class(ControlFlow), intent(in) :: CF
             Type(SpecyOne),intent(in),Target  :: SO
             Type(GasOne),intent(in),Target :: GO(CF%Ng)
             Type(MCCBundle),intent(inout)  :: MCCB(1:CF%Ng)
             
             Type(SigmaNormalized) :: SN 
             Integer(4) ::  i
             Logical :: Inited
             Associate(Ns=>CF%Ns,Ng=>CF%Ng)
                 Do i=1,Ng
                        !MCCBundleGlobal(0,i)
                        MCCB(i)%SO=>SO
                        MCCB(i)%GO=>GO(i)
                        
                        Call SigmaNormalization(SN,SO,GO(i),Inited)
                        !Call UpdateReactionIndex(SN,GO(i))
                        Call ProbilityNonReactive(MCCB(i),SN,SO,GO(i))
                        MCCB(i)%CollisionRatio=1.d0-DExp(-MCCB(i)%SigmaMax*CF%dt)
                 End Do
             End Associate
             return
   End subroutine MCCBundleElelctronInit

   
      subroutine  MCCBundleIonInit(CF,SO,GO,MCCB) 
             Implicit none
             Class(ControlFlow), intent(in) :: CF
             Type(SpecyOne),intent(in),Target  :: SO(1:CF%Ns)
             Type(GasOne),intent(in),Target  :: GO(1:CF%Ng)
             Type(MCCBundle),intent(inout)  :: MCCB(1:CF%Ns,1:CF%Ng)
             Type(SigmaNormalized) :: SN 
             Integer(4) ::  i,j
             Logical :: Inited
             Associate(Ns=>CF%Ns,Ng=>CF%Ng)
                     Do i=1,Ng
                            Do j=1,Ns
                                  MCCB(j,i)%SO=>SO(j)
                                  MCCB(j,i)%GO=>GO(i)
                                  Call UpdateReactionIndex(SN,GO(i),GO(SO(j)%GasIndex))
                                  If (SO(j)%GasIndex==i) Then
                                      Call SigmaNormalization(SN,SO(j),GO(i),Inited)
                                      If (Inited) Then
                                           If(SN%Model==3) Then
                                               Call ProbilityReactive(MCCB(j,i),SN)
                                           Else   
                                           Call ProbilityNonReactive(MCCB(j,i),SN,SO(j),GO(i))
                                           ENd If
                                       Else    
                                           Call  ProbilityHardShpere(MCCB(j,i),SN,SO(j),GO(i))
                                      ENd If
                                  Else
                                      !MCCB(j,i)%Model=0
                                      Call  ProbilityHardShpere(MCCB(j,i),SN,SO(j),GO(i))
                                  ENd IF
                            MCCB(j,i)%CollisionRatio=1.d0-DExp(-MCCB(j,i)%SigmaMax*CF%dt)
                            End Do
                     End DO
             End Associate
             return
      End subroutine MCCBundleIonInit
      
      
      

      
      subroutine  ProbilityNonReactive(MCB,SN,SO,GO)
                Implicit none
                Type(MCCBundle),intent(inout) :: MCB
                Type(SigmaNormalized),intent(in) :: SN 
                Type(SpecyOne),intent(in) :: SO
                Type(GasOne),intent(in) :: GO
                Real(8) :: NEnergy,Energy,V,Max,Miu
                Integer(4) :: i,j
                
                MCB%Model=SN%Model
                MCB%NReaction=SN%NReaction
                MCB%NSigma=SN%NSigma
                MCB%EnergyMin=SN%EnergyMin
                MCB%EnergyInterval=SN%EnergyInterval
                MCB%EnergyMax=10.d0**SN%EnergyMax
                
                !If(Allocated(SN%Reaction)) Deallocate(SN%Reaction)
                Allocate(MCB%Reaction(SN%NReaction))
                MCB%Reaction=SN%Reaction
                !If(Allocated(SN%Sigma)) Deallocate(SN%Sigma)
                Allocate(MCB%Probility(SN%NReaction,SN%NSigma))

                Associate(Sigma=>SN%Sigma,Emin=>MCB%EnergyMin,Eint=>MCB%EnergyInterval,Nr=>MCB%NReaction,Ns=>MCB%NSigma,Probility=>MCB%Probility)
                        !Probility=SN%Sigma
                        
                        Miu=SO%Mass*GO%Mass/(SO%Mass+GO%Mass)

                        do i=1,Ns
                             NEnergy=dble(i-1)
                             Energy=10.d0**(Emin+dble(i-1)*Eint)
                              V=DSqrt(2.d0*Energy*JtoeV/Miu)
                              do j=1,Nr
                                    Probility(j,i)=Sigma(j,i)*V*GO%InitDensity
                              end do
                        end do
                        
                      Max=0.d0  
                      do i=1,Ns
                           do j=2,Nr
                                Probility(j,i)=Probility(j,i)+Probility(j-1,i)
                            end do
                      end do
                      
                      Max=MAXVAL(Probility)
                      Probility=Probility/Max
                      MCB%SigmaMax=Max
                      Call MCB%Dump
                End Associate
             return 
      End subroutine  ProbilityNonReactive
      
      Subroutine  ProbilityHardShpere(MCB,SN,SO,GO)
                Implicit none
                Type(MCCBundle),intent(inout) :: MCB
                Type(SigmaNormalized),intent(inout) :: SN 
                Type(SpecyOne),intent(in) :: SO
                Type(GasOne) :: GO
                Real(8) :: V,SigmaMax,Max,Miu,NEnergy,Energy
                Integer(4) :: i,j
                Type(ReactionOne) :: RO(2)=(/ReactionOne(0,101,0,0.d0), ReactionOne(0,103,0,0.d0)/)
                If (GO%Ns==1) Then
                        MCB%Model=2
                        MCB%NReaction=2
                        MCB%NSigma=NSigmaIon
                        MCB%EnergyMin=SN%EnergyMin
                        MCB%EnergyInterval=SN%EnergyInterval
                        MCB%EnergyMax=10**SN%EnergyMax 

                        V=DSqrt(2.d0*MCB%EnergyMax*JtoeV/GO%Mass)
                        SigmaMax=2.d0*PI*0.5d0*(GO%Radius+SO%Radius)*0.5d0*(GO%Radius+SO%Radius)

                        do i=1,1
                              RO(i)%Reactant=SO%SpecyIndex
                              RO(i)%Resultant=SO%SpecyIndex
                        End Do
                
                        do i=2,2
                              RO(i)%Reactant=SO%SpecyIndex
                              RO(i)%Resultant=GO%IndexStart+1
                        End Do
                
                        !If(Allocated(SN%Reaction)) Deallocate(SN%Reaction)
                        Allocate(MCB%Reaction(MCB%NReaction))
                        !If(Allocated(SN%Sigma)) Deallocate(SN%Sigma)
                        Allocate(MCB%Probility(MCB%NReaction,MCB%NSigma))
                        MCB%Reaction=RO
                
                
                
                        Associate(Sigma=>SN%Sigma,Emin=>MCB%EnergyMin,Eint=>MCB%EnergyInterval,Nr=>MCB%NReaction,Ns=>MCB%NSigma,Probility=>MCB%Probility)
                        Probility(1,:)=0.5d0
                        Probility(2,:)=1.d0
                        Miu=SO%Mass*GO%Mass/(SO%Mass+GO%Mass)
                                do i=1,Ns
                                     NEnergy=dble(i-1)
                                     Energy=10.d0**(Emin+dble(i-1)*Eint)
                                      V=DSqrt(2.d0*Energy*JtoeV/Miu)
                                      do j=1,Nr
                                            Probility(j,i)=Probility(j,i)*SigmaMax*V*GO%InitDensity
                                      end do
                                end do
                        
                              Max=MAXVAL(Probility)
                              Probility=Probility/Max
                              MCB%SigmaMax=Max
                              Call MCB%Dump
                        End Associate
                Else
                        MCB%Model=SN%Model
                        MCB%NReaction=SN%NReaction
                        MCB%NSigma=SN%NSigma
                        MCB%EnergyMin=SN%EnergyMin
                        MCB%EnergyInterval=SN%EnergyInterval
                        MCB%EnergyMax=10.d0**SN%EnergyMax

                        V=DSqrt(2.d0*MCB%EnergyMax*JtoeV/GO%Mass)
                        SigmaMax=2.d0*PI*0.5d0*(GO%Radius+SO%Radius)*0.5d0*(GO%Radius+SO%Radius)

                        do i=1,1
                              RO(i)%Reactant=SO%SpecyIndex
                              RO(i)%Resultant=SO%SpecyIndex
                        End Do
                
                        !If(Allocated(SN%Reaction)) Deallocate(SN%Reaction)
                        Allocate(MCB%Reaction(MCB%NReaction))
                        !If(Allocated(SN%Sigma)) Deallocate(SN%Sigma)
                        Allocate(MCB%Probility(MCB%NReaction,MCB%NSigma))
                        MCB%Reaction=RO

                        Associate(Sigma=>SN%Sigma,Emin=>MCB%EnergyMin,Eint=>MCB%EnergyInterval,Nr=>MCB%NReaction,Ns=>MCB%NSigma,Probility=>MCB%Probility)
                        Probility(1,:)=1.d0
                        Miu=SO%Mass*GO%Mass/(SO%Mass+GO%Mass)
                                do i=1,Ns
                                     NEnergy=dble(i-1)
                                     Energy=10.d0**(Emin+dble(i-1)*Eint)
                                      V=DSqrt(2.d0*Energy*JtoeV/Miu)
                                      do j=1,Nr
                                            Probility(j,i)=Probility(j,i)*SigmaMax*V*GO%InitDensity
                                      end do
                                end do
                        
                              Max=MAXVAL(Probility)
                              Probility=Probility/Max
                              MCB%SigmaMax=Max
                              Call MCB%Dump
                        End Associate
                ENd If
             return 
      End subroutine  ProbilityHardShpere
      
      subroutine  ProbilityReactive(MCB,SN)
                implicit none
                Type(MCCBundle),intent(inout) :: MCB
                Type(SigmaNormalized),intent(inout) :: SN 
                Integer(4) :: i,j
                
                
                MCB%Model=SN%Model
                MCB%NReaction=SN%NReaction
                MCB%NSigma=SN%NSigma
                MCB%EnergyMin=SN%EnergyMin
                MCB%EnergyInterval=SN%EnergyInterval
                MCB%EnergyMax=10.d0**SN%EnergyMax
                
                Allocate(MCB%Reaction(SN%NReaction))
                MCB%Reaction=SN%Reaction
                !If(Allocated(SN%Sigma)) Deallocate(SN%Sigma)
                Allocate(MCB%Probility(SN%NReaction,SN%NSigma))
                
                MCB%Probility=SN%Sigma
 
               do i=1,MCB%NSigma
                        do j=2,MCB%NReaction
                            MCB%Probility(j,i)=MCB%Probility(j,i)+MCB%Probility(j-1,i)
                         end do
               end do
               Call MCB%Dump
              return
     End subroutine  ProbilityReactive
      
      !subroutine  CollisionRatioReactive(InputSpecy,InputGas,NGas,dt,CollisionRatio)
      !       implicit none
      !       Type(GasOne),intent(in) :: GO
      !       Type(SpecyOne),intent(in) :: SO
      !       Real(8),intent(in) :: NGas,dt
      !       Real(8),intent(inout) :: CollisionRatio
      !       Real(8) :: Miu,BetaMax,Radius
      !       Miu=InputSpecy%Mass*InputGas%MGas/(InputGas%MGas+InputSpecy%Mass)
      !       BetaMax=InputGas%BetaMax
      !       Radius=InputGas%Radius
      !       CollisionRatio=DSQRT(PI*Radius*a0*a0*a0*ElectronCharge*ElectronCharge/(Epsilon*Miu))*BetaMax*BetaMax*Ngas*dt
      !       return 
      !  End subroutine  CollisionRatioReactive
End Module ModuleMCCInitialization