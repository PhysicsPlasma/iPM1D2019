Module DiagnosticsMomentum
    Use ModuleGrid
    Use ModuleField
    Use ModuleParticleBundle
    Use ModuleParticleBoundary
    Implicit none
    Type ParticleMomentumOne!(Nx)
                        !Integer(4),Len :: Nx=100
                        !Integer(4) :: IOIndex=1
                        !Integer(4) :: XStart=0,XEnd=NxMax-1
                        Integer(4) :: Nx=NxMax,Timer=0
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        !Real(8) :: QdM
                        Real(8) ::  RhoOne(1:NxMax),ChiOne(1:NxMax)
                        Real(8) ::  JxOne(1:NxMax),JyOne(1:NxMax),JzOne(1:NxMax)
                        Real(8) ::  TOne(1:NxMax)
                   EndType ParticleMomentumOne
    contains

     Subroutine  DiagParticleFieldPeriod(GD,NSpecy,PB,FG,Mode)
         Implicit none
         Class(*),intent(inout)  :: GD
         Integer(4),intent(in) ::  NSpecy
         Type(ParticleBundle),intent(in) :: PB(0:NSpecy)
         Type(Field),intent(in) :: FG
         Integer(4),intent(in) ::  Mode
         !Integer(4),parameter :: NSpecyMax=2_4 
         Type(ParticleMomentumOne) ::  TempPMO
         Integer(4) :: i,j,Shift
         Real(8) :: HeatingRate(FG%Nx)
         Select Type (GD)
             Type is (Grid1D(*,*))
                 Select Case (Mode)

                    case(0)
                        Shift=1
                        Do i=0, NSpecy
                                      Call WeightingParticleMomentum(PB(i),TempPMO)
                                      Call GD%Update(GD%Nx,TempPMO%RhoOne,Shift)
                                      Call GD%Update(GD%Nx,TempPMO%JxOne,Shift)
                                      Call GD%Update(GD%Nx,TempPMO%TOne,Shift)
                                      HeatingRate=TempPMO%JxOne*FG%Ex
                                      Call GD%Update(GD%Nx,HeatingRate,Shift)
                        End do
                        Call GD%Update(GD%Nx,FG%Ex,Shift)
                        Call GD%Update(GD%Nx,FG%Phi,Shift)
                        Call GD%Update(GD%Nx,FG%Rho,Shift)
                        GD%Timer=GD%Timer+1
                     Case(1) 
                             Call GD%Rescale
                             Call GD%Dump(1)
                             GD%Value=0.d0
                             GD%Timer=0
                     Case(2)
                         Call GD%Dump(0)
                     case default
                         
                     End Select
                Type is (Grid2D(*,*,*))
                    Select Case (Mode)

                        case(0)
                            Shift=1
                             Do i=0, NSpecy
                                          Call WeightingParticleMomentum(PB(i),TempPMO)
                                          Call GD%Update(GD%Nx,TempPMO%RhoOne,Shift)
                                          Call GD%Update(GD%Nx,TempPMO%JxOne,Shift)
                                          Call GD%Update(GD%Nx,TempPMO%TOne,Shift)
                                          HeatingRate=TempPMO%JxOne*FG%Ex
                                          Call GD%Update(GD%Nx,HeatingRate,Shift)
                            End do
                            Call GD%Update(GD%Nx,FG%Ex,Shift)
                            Call GD%Update(GD%Nx,FG%Phi,Shift)
                            Call GD%Update(GD%Nx,FG%Rho,Shift)
                            GD%Timer=GD%Timer+1
                        Case(1) 
                                 Call GD%Rescale
                                 Call GD%Dump(1)
                                 GD%Value=0.d0
                                 GD%Timer=0
                         Case(2)
                             Call GD%Dump(0) 
                         case default
                     End Select     
                End select
        Return    
     End Subroutine  DiagParticleFieldPeriod

  subroutine WeightingParticleMomentum(PB,PMO)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleMomentumOne),intent(inout) :: PMO
                Real(8) :: RhoFactor,ChiFactor,JFactor,TFactor
                Real(8) :: S1,S2,Energy
                Integer(4) :: i,N
                PMO%RhoOne=0.d0
                PMO%JxOne=0.d0
                !PMO%JyOne=0.d0
                !PMO%JzOne=0.d0
                PMO%TOne=0.d0
                do i=1,PB%Npar
                   N=Ceiling(PB%PO(i)%X)

                   S1=Dble(N)-PB%PO(i)%X
                   S2=1.d0-S1
                   PMO%RhoOne(N)=PMO%RhoOne(N)+S1
                   PMO%RhoOne(N+1)=PMO%RhoOne(N+1)+S2
                   
                   PMO%JxOne(N)=PMO%JxOne(N)+S1*PB%PO(i)%Vx
                   PMO%JxOne(N+1)=PMO%JxOne(N+1)+S2*PB%PO(i)%Vx

                   Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
                
                   PMO%TOne(N)=PMO%TOne(N)+S1*Energy
                   PMO%TOne(N+1)=PMO%TOne(N+1)+S2*Energy
                end do
                
                RhoFactor=PB%Weight
                PMO%RhoOne=PMO%RhoOne*RhoFactor
                ChiFactor=0.5d0*PB%Charge/PB%Mass*PB%dt*PB%dt/Epsilon*PB%Charge
                PMO%ChiOne=PMO%RhoOne*ChiFactor
                JFactor=PB%Charge*PB%Weight*PB%VFactor
                PMO%JxOne=PMO%JxOne*JFactor
                TFactor=PB%Weight
                PMO%TOne=PMO%TOne*TFactor
                do i=1,PMO%Nx
                    if (PMO%RhoOne(i)>0.d0) then
                        PMO%TOne(i)=PMO%TOne(i)/PMO%RhoOne(i)
                    end if
                end do
                PMO%RhoOne(1)=2.d0*PMO%RhoOne(1)
                PMO%RhoOne(PMO%Nx)=2.d0*PMO%RhoOne(PMO%Nx)
                PMO%JxOne(1)=2.d0*PMO%JxOne(1)
                PMO%JxOne(PMO%Nx)=2.d0*PMO%JxOne(PMO%Nx)
                !PMO%TOne(1)=2.d0*PMO%TOne(1)
                !PMO%TOne(PMO%Nx)=2.d0*PMO%TOne(PMO%Nx)
                PMO%ChiOne(1)=2.d0*PMO%ChiOne(1)
                PMO%ChiOne(PMO%Nx)=2.d0*PMO%ChiOne(PMO%Nx)
                return
  end subroutine WeightingParticleMomentum
    End Module DiagnosticsMomentum
!!   
!!    
!!     
!!    subroutine WeightingParticleEDF(PB,PEDF)
!!                implicit none
!!                Type(ParticleBundle),intent(in) :: PB
!!                Type(ParticleEDF),intent(inout) :: PEDF
!!                Real(8) :: Energy,Frac
!!                Integer(4) :: i,N
!!                
!!                PEDF%EDF=0.d0
!!                PEDF%EDFNormalized=0.d0
!!                
!!                Frac=1.d0/DBLE(PB%Npar)
!!                do i=1,PB%Npar
!!                   Call CalEnergy (PB%Mass,PB%PO(i),Energy)
!!                   Energy=Energy*PB%VFactor*PB%VFactor/JtoeV
!!                   N=Ceiling( Energy/PEDF%EnergyInterval)
!!                   If (N>=1.and.N<PEDF%Ne) Then
!!                           PEDF%EDF(N)=PEDF%EDF(N)+PB%Weight/PB%XMax
!!                           PEDF%EDFNormalized(N)=PEDF%EDFNormalized(N)+Frac
!!                   End IF
!!                end do
!!                return
!!    end subroutine WeightingParticleEDF
!!    
!!        subroutine WeightingParticleCollisionRate(PB,PCR,MCCB)
!!                implicit none
!!                Type(ParticleBundle),intent(in) :: PB
!!                Type(ParticleCollisionRate),intent(inout) :: PCR
!!                Type(MCCBundle),intent(in) ::  MCCB
!!                Real(8) :: EnergyInterval,Energy,S1,S2
!!                Real(8) :: TempProbility(MCCB%NReaction)
!!                Integer(4) :: i,j,Nc,Np,Index,Upper,Center,Lower,NReaction
!!
!!                PCR%CollisionRate=0.d0
!!                NReaction=MCCB%NReaction
!!                PCR%NReaction=NReaction
!!                do j=1,PB%Npar
!!                   Call CalEnergy (PB%Mass,PB%PO(j),Energy)
!!                   Energy=Energy*PB%VFactor*PB%VFactor/JtoeV
!!                   TempProbility=0.d0 
!!                   EnergyInterval=MCCB%EnergyInterval
!!                   
!!                   Nc=Int(Energy/EnergyInterval)
!!                    If(Nc<2) then
!!                         TempProbility=MCCB%Probility(1:NReaction)
!!                    Else if (Nc<MCCB%NSigma) then
!!                         Index=Nc*NReaction
!!                         do i=1,NReaction
!!                               Center=Index+i
!!                               Upper=Center+NReaction
!!                               Lower=Center-NReaction
!!                               If(Energy<=MCCB%Reaction(i)%Threshold) Then
!!                                   TempProbility(i:NReaction)=0.d0
!!                                   exit       
!!                              else
!!                                   S1=(Energy-(dble(Nc)*EnergyInterval))/EnergyInterval
!!                                   S2=1.d0-S1 
!!                                   TempProbility(i)=S1*MCCB%Probility(Lower)+S2*MCCB%Probility(Upper)
!!                              end If
!!                          end do
!!                   else
!!                          Index=(MCCB%NSigma-1)*NReaction 
!!                          TempProbility=MCCB%Probility(Index+1:Index+NReaction)
!!                   End if
!!                   TempProbility=TempProbility*MCCB%SigmaMax*MCCB%Ng
!!                   Do i=NReaction,2
!!                         If ((TempProbility(i)-TempProbility(i-1))>=MinReal) Then
!!                              TempProbility(i)=TempProbility(i)-TempProbility(i-1)
!!                         End If
!!                   End do
!!                   
!!                   Np=Ceiling(PB%PO(j)%X)
!!                   Do  i=1, NReaction
!!                          PCR%CollisionRate(Np,i)=TempProbility(i)
!!                   End DO
!!                end do
!!                return
!!        end subroutine WeightingParticleCollisionRate
!!        !
!!        subroutine WeightingParticleFlux(PE)
!!                implicit none
!!                Type(ParticleElectrode),intent(inout)  :: PE
!!                Real(8) :: Energy
!!                Integer(4) :: i,Timer,NAva
!!
!!                PE%Flux=0.d0
!!                PE%EnergyFlux=0.d0
!!                NAva=PE%Period/PE%NFlux
!!                do i=1,PE%PBLower%Npar
!!                   Call CalEnergy (PE%PBLower%Mass,PE%PBLower%PO(i),Energy)
!!                   Energy=Energy*PE%PBLower%VFactor*PE%PBLower%VFactor/JtoeV
!!                   !Write(*,*) Energy,"Lo"
!!                   Timer=Int(PE%PBLower%PO(i)%Az)/NAva+1
!!                   If(Timer>=1.and.Timer<=PE%NFLux) Then
!!                       PE%Flux(Timer,1)=PE%Flux(Timer,1)+PE%PBLower%Weight
!!                       PE%EnergyFlux(Timer,1)=PE%EnergyFlux(Timer,1)+PE%PBLower%Weight*Energy
!!                   End IF
!!                End do
!!                
!!                do i=1,PE%PBUpper%Npar
!!                   Call CalEnergy (PE%PBUpper%Mass,PE%PBUpper%PO(i),Energy)
!!                   Energy=Energy*PE%PBUpper%VFactor*PE%PBUpper%VFactor/JtoeV
!!                   !Write(*,*) Energy,"Up"
!!                   Timer=Int(PE%PBUpper%PO(i)%Az)/NAva+1
!!                   If(Timer>=1.and.Timer<=PE%NFLux) Then
!!                       PE%Flux(Timer,2)=PE%Flux(Timer,2)+PE%PBUpper%Weight
!!                       PE%EnergyFlux(Timer,2)=PE%EnergyFlux(Timer,2)+PE%PBUpper%Weight*Energy
!!                   End IF
!!                End do
!!                PE%Flux=PE%Flux/(dble(NAva)*PE%PBLower%dt)
!!                PE%EnergyFlux=PE%EnergyFlux/(dble(NAva)*PE%PBLower%dt)
!!                return
!!        end subroutine WeightingParticleFlux
!!  
!!    
!!                        Integer(4),Parameter,Private :: NeMax=1000
!!                    Type ParticleEDF
!!                        Integer(4) :: Ne=NeMax
!!                        Real(8) :: EnergyInterval,Mass
!!                        Real(8) :: EDF(NeMax),EDFNormalized(NeMax)
!!                    EndType ParticleEDF
!!                    
!!                    Type ParticleCollisionRate
!!                        Integer(4) :: Nx=NxMax
!!                        Integer(4) :: NReaction
!!                        Real(8) :: CollisionRate(NxMax,NReactionMax)
!!                    EndType ParticleCollisionRate
!!                    
!!                    Integer(4),Parameter,Private :: NfMax=1000
!!                    Type ParticleElectrode
!!                        Integer(4) :: Timer=0,Period,NFlux=100              
!!                        Real(8) :: Dt
!!                        Type(ParticleBundle) :: PBLower,PBUpper
!!                        Real(8) :: Flux(1:NfMax,1:2)=0.d0,EnergyFlux(1:NfMax,1:2)=0.d0
!!    EndType ParticleElectrode
!!    
!!    
!!             subroutine DiagParticleTestParticle(PB,FG,Mode)
!!             Implicit none
!!            Type(ParticleBundle),intent(in) :: PB
!!            Type(Field),intent(in) :: FG
!!            Integer(4),intent(in) ::  Mode
!!            Type(ParticleBundle),save :: TestPB,TempPB
!!            Type(ParticleOne),save :: TempPO
!!            Integer(4) :: i
!!            Real(8) :: Energy
!!            Select Case (Mode)
!!                       Case(-1)
!!                       TempPB=PB
!!                       TempPB%Npar=1
!!                       TempPB%Timer=0
!!                       TempPB%PO(1)%X=1.d0
!!                       TempPB%PO(1)%Vx=0.d0
!!                       TempPB%PO(1)%Vy=0.d0
!!                       TempPB%PO(1)%Vz=0.d0
!!                       TempPB%PO(1)%Ax=0.d0
!!                       TempPB%PO(1)%Ay=0.d0
!!                       TempPB%PO(1)%Az=0.d0
!!                       TestPB=TempPB
!!                       TestPB%Timer=1
!!                       
!!                       
!!                       Case(0)
!!                           If (TempPB%Timer<=NParMax) Then
!!                               If ((TempPB%PO(1)%X>TempPB%XMin).and.(TempPB%PO(1)%X<TempPB%XMax)) Then
!!                               Call ParticleMove(TempPB,FG)
!!                               TempPO=TempPB%PO(1)
!!                               Call CalEnergy (TestPB%Mass,TempPO,Energy)
!!                               Energy=Energy*PB%VFactor*PB%VFactor/JtoeV
!!                               TempPO%Az=Energy
!!                                Call AddParticle(TempPO,TestPB)
!!                                End If
!!                               TempPB%Timer=TempPB%Timer+1
!!                               TestPB%Timer=TestPB%Timer+1
!!                           End IF
!!                        Case(1)
!!                            Call DumpParticle(TestPB,1)
!!                        Case(2)
!!                         case default
!!                         End Select
!!           return
!!        end subroutine DiagParticleTestParticle  
!!    
!!    
!!    
!!        subroutine DiagParticleElectrode(PE,PB,FG,Mode)
!!             Implicit none
!!            Type(ParticleElectrode),intent(inout)  :: PE
!!            Type(ParticleBundle),intent(in) :: PB
!!            Type(Field),intent(in) :: FG
!!            !Type(ParticleBoundaryOne),intent(in) :: PBDO
!!            Integer(4),intent(in) ::  Mode
!!            Type(ParticleBundle),save :: TempPB
!!            Type(Grid1D(Nx=1000,Ns=2)),save :: GD_EDF
!!            Type(Grid1D(Nx=737,Ns=4)),save :: GD_Flux
!!            Integer(4) :: i,ParticleTimer
!!            
!!            !Allocate(Type(Grid1D(Nx=1000,Ns=2))::GD_Flux)
!!            
!!            Select Case (Mode)
!!                 Case(-1)
!!                        PE%Timer=0
!!                        PE%Period=FG%Period
!!                        PE%NFlux=FG%Period
!!                        PE%Dt=FG%Dt
!!                        PE%PBLower=PB
!!                        PE%PBLower%NPar=0
!!                        PE%PBLower%Timer=0
!!                        PE%PBUpper=PB
!!                        PE%PBUpper%NPar=0
!!                        PE%PBUpper%Timer=0
!!                        PE%Flux=0.d0
!!                        PE%EnergyFLux=0.d0
!!                        Case(0)
!!                            TempPB=PB
!!                            Call ParticleMove(TempPB,FG)
!!                            ParticleTimer=Mod(PE%Timer,PE%Period)
!!                            do i=TempPB%NPar,1,-1
!!                                 If (TempPB%PO(i)%X<=TempPB%XMin) then
!!                                     TempPB%PO(i)%Ay=dble(ParticleTimer)
!!                                     Call AddParticle(TempPB%PO(i),PE%PBLower)
!!                                 else If(TempPB%PO(i)%X>=TempPB%XMax) then
!!                                     TempPB%PO(i)%Ay=dble(ParticleTimer)
!!                                     Call AddParticle(TempPB%PO(i),PE%PBUpper)
!!                                 end if
!!                            end do
!!                            PE%Timer=PE%Timer+1
!!                       Case(1)
!!                                 !Allocate(Type(Grid1D(100,2))::GD_Flux)
!!                                 Call DumpParticle(PE%PBLower,1)
!!                                 Call DumpParticle(PE%PBUpper,1)
!!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBLower,-1)
!!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBLower,0)
!!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBLower,1)
!!                                 
!!                                 
!!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBUpper,-1)
!!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBUpper,0)
!!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBUpper,1)
!!                                 
!!                                 Call DiagParticleFlux(GD_Flux,PE,-1)
!!                                 Call DiagParticleFlux(GD_Flux,PE,0)
!!                                 Call DiagParticleFlux(GD_Flux,PE,1)
!!                                 !
!!                                 !PE%Timer=0
!!                                 !PE%PBUpper%NPar=0
!!                                 !PE%PBUpper%Timer=0 
!!                                 !PE%PBUpper%NPar=0
!!                                 !PE%PBUpper%Timer=0 
!!                         Case(2)
!!                         case default
!!                         End Select
!!           return
!!        end subroutine DiagParticleElectrode
!!        
!!        subroutine DiagParticleFlux(GD,PE,Mode)
!!            Implicit none
!!            Type(Grid1D(*,*))  :: GD
!!            Type(ParticleElectrode),intent(inout)  :: PE
!!            Integer(4),intent(in) ::  Mode
!!            Integer(4) :: Shift,i
!!                     Select Case (Mode)
!!                        Case(-1)
!!                            Call GridInitialization(GD,PE%PBLower%Period,PE%PBLower%Dx,PE%PBLower%Dt)
!!                        Case(0)
!!                            Shift=1
!!                            PE%NFlux=100 !GD%Nx
!!                            Call WeightingParticleFlux(PE)
!!                            Do i=1,2
!!                                Call GD%Update(PE%Nflux,PE%Flux(1:PE%NFlux,i),Shift)
!!                                Call GD%Update(PE%Nflux,PE%EnergyFlux(1:PE%NFlux,i),Shift)
!!                            END DO
!!                            GD%Timer=GD%Timer+1
!!                         Case(1) 
!!                                 Call GD%Rescale
!!                                 Call GD%Dump(1)
!!                                 GD%Value=0.d0
!!                                 GD%Timer=0
!!                         case default
!!                         End Select
!!            return
!!        end subroutine DiagParticleFlux
!!
!!       subroutine DiagParticleCollisionRateOne(GD,PB,MCCB,Mode)
!!            Implicit none
!!            Class(*),intent(inout)  :: GD
!!            Type(ParticleBundle),intent(in) :: PB
!!            Type(MCCBundle),intent(in) ::  MCCB
!!            Integer(4),intent(in) ::  Mode
!!            Type(ParticleCollisionRate) :: PCR
!!            Integer(4) :: Shift,i
!!            Select Type (GD)
!!                  Type is (Grid1D(*,*))
!!                     Select Case (Mode)
!!                        Case(-1)
!!                            Call GridInitialization(GD,PB%Period,PB%Dx,PB%Dt)
!!                        case(0)
!!                            Shift=1
!!                            Call WeightingParticleCollisionRate(PB,PCR,MCCB)
!!                            do i=1,MCCB%NReaction
!!                                 Call GD%Update(GD%Nx,PCR%CollisionRate(:,i),Shift)
!!                            End Do
!!                            GD%Timer=GD%Timer+1
!!                         Case(1) 
!!                                 Call GD%Rescale
!!                                 Call GD%Dump(1)
!!                                 GD%Value=0.d0
!!                                 GD%Timer=0
!!                         Case(2)
!!                             Call GD%Dump(0)
!!                         case default
!!                         End Select
!!                    Type is (Grid2D(*,*,*))
!!                        Select Case (Mode)
!!                            Case(-1)
!!                                Call GridInitialization(GD,PB%Period,PB%Dx,PB%Dt)
!!                            case(0)
!!                                Shift=1
!!                                Call WeightingParticleCollisionRate(PB,PCR,MCCB)
!!                                do i=1,MCCB%NReaction
!!                                      Call GD%Update(GD%Nx,PCR%CollisionRate(:,i),Shift)
!!                                End Do
!!                                GD%Timer=GD%Timer+1
!!                            Case(1) 
!!                                     Call GD%Rescale
!!                                     Call GD%Dump(1)
!!                                     GD%Value=0.d0
!!                                     GD%Timer=0
!!                             Case(2)
!!                                 Call GD%Dump(0)
!!                             case default
!!                         End Select     
!!                    End select
!!            return
!!        end subroutine DiagParticleCollisionRateOne
!!
!!        subroutine DiagParticleEDFOne(GD,PB,Mode)
!!            Implicit none
!!            Class(*),intent(inout)  :: GD
!!            Type(ParticleBundle),intent(in) :: PB
!!            Integer(4),intent(in) ::  Mode
!!            Type(ParticleEDF) :: PEDF
!!            Integer(4) :: Shift
!!
!!            If (PB%SpecyIndex==0) Then
!!                PEDF%EnergyInterval=0.1d0
!!            Else
!!                PEDF%EnergyInterval=1.d0
!!            End If
!!            PEDF%Mass=PB%Mass
!!
!!            Select Type (GD)
!!                  Type is (Grid1D(*,*))
!!                     Select Case (Mode)
!!                        Case(-1)
!!                            Call GridInitialization(GD,PB%Period,PEDF%EnergyInterval,PB%Dt)
!!                        case(0)
!!                            
!!                            PEDF%Ne=GD%Nx
!!                            Shift=1
!!                            Call WeightingParticleEDF(PB,PEDF)
!!                            Call GD%Update(GD%Nx,PEDF%EDF,Shift)
!!                            Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
!!                            GD%Timer=GD%Timer+1
!!                         Case(1) 
!!                                 Call GD%Rescale
!!                                 Call GD%Dump(1)
!!                                 GD%Value=0.d0
!!                                 GD%Timer=0
!!                         Case(2)
!!                             Call GD%Dump(0)
!!                         case default
!!                         End Select
!!                    Type is (Grid2D(*,*,*))
!!                        Select Case (Mode)
!!                            Case(-1)
!!                                Call GridInitialization(GD,PB%Period,PEDF%EnergyInterval,PB%Dt)
!!                            case(0)
!!                                Shift=1
!!                                Call WeightingParticleEDF(PB,PEDF)
!!                                Call GD%Update(GD%Nx,PEDF%EDF,Shift)
!!                                Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
!!                                GD%Timer=GD%Timer+1
!!                            Case(1) 
!!                                     Call GD%Rescale
!!                                     Call GD%Dump(1)
!!                                     GD%Value=0.d0
!!                                     GD%Timer=0
!!                             Case(2)
!!                                 Call GD%Dump(0)
!!                             case default
!!                         End Select     
!!                    End select
!!            return
!!        end subroutine DiagParticleEDFOne
!!
!!        !    ! subroutine DiagTrajectoryINit(InputParticle,InputNx,InputPeriod)
!!!    !    Implicit none
!!!    !     Integer(4),intent(in) ::  InputNx
!!!    !     Type(ParticleBundle),intent(in) ::  InputParticle
!!!    !     Type(PICParticlePara),intent(in) :: InputPICParticlePara
!!!    !    
!!!    !     Type(ParticleBundle)::ParticleTempBundle
!!!    !     Type(ParticleOne) :: ParticleTemp
!!!    !     Real(8),save :: XRandom=2.0d0,VFactor,Xmin,Xmax
!!!    !     Integer(4),parameter:: NRandom=8
!!!    !     Type(Gas),parameter :: TraEGas=(Gas(1,0,9.1095d-31,0.026d0*11605.d0,0.d0,0.d0)) 
!!!    !      
!!!    !      
!!!    !    Integer(4),intent(in) ::  InputPeriod
!!!    !    Integer(4),save :: i,j,nx,Period,Timer,Npar1
!!!    !   ! Type(Grid),save :: TempGrid
!!!    !    Type(ParticleBundle),save:: trajectory
!!!    !       nx=InputNx
!!!    !    
!!!    !       trajectory%XFactor=InputParticle%XFactor
!!!    !       trajectory%VFactor=InputParticle%VFactor
!!!    !       trajectory%Charge=InputParticle%Charge
!!!    !       trajectory%Weight=InputParticle%Weight
!!!    !       ParticleTempBundle=InputParticle
!!!    !   
!!!    !        Xmin=dble(1-1)
!!!    !        Xmax=dble(Nx-1)
!!!    !    
!!!    !        Timer=0
!!!    !        NPar1=1
!!!    !        Nx=InputNx
!!!    !        Period=InputPeriod
!!!    !        !TempGrid%Nx=1
!!!    !        !TempGrid%Ny=Period
!!!    !        !TempGrid%Value=0.d0
!!!    !        ParticleTempBundle%NPar=1
!!!    !        trajectory%NPar=2*period
!!!    !        Write(trajectory%Name,*) 'eTrajectory'
!!!    !        
!!!    !        
!!!    !      ! CALL RANDOM_NUMBER(R)
!!!    !            ParticleTemp%X=63.d0!XMin+XRandom
!!!    !            ParticleTemp%Ax=0.0
!!!    !            ParticleTemp%Ay=0.0
!!!    !            ParticleTemp%Az=0.0
!!!    !        ParticleTemp%Vx=0.d0!InputParticle%PO(NRandom)%Vx
!!!    !        ParticleTemp%Vy=0.d0!InputParticle%PO(NRandom)%Vy
!!!    !        ParticleTemp%Vz=0.d0!InputParticle%PO(NRandom)%Vz
!!!    !        !Call Maxwellian(TraEGas,ParticleTemp)    
!!!    !        trajectory%PO(Npar1)=ParticleTemp
!!!    !        
!!!    !        ParticleTempBundle%PO(1)=ParticleTemp
!!!    !        Return 
!!!    !Entry DiagTrajectory(InputPICParticlePara,InputNx)
!!!    !   If (trajectory%PO(Npar1)%X>XMin.and.trajectory%PO(Npar1)%X<XMax) then
!!!    !           Call PushES1DEt(ParticleTempBundle,Nx,InputPICParticlePara%Ex,InputPICParticlePara%Ey,InputPICParticlePara%Bx,InputPICParticlePara%By)
!!!    !   endif
!!!    !       Timer=Timer+1
!!!    !        NPar1=NPar1+1
!!!    !        trajectory%PO(Npar1)=ParticleTempBundle%PO(1)
!!!    !        Return 
!!!    !   Entry  DiagTrajectoryFinal()
!!!    !     Call DumpParticle(trajectory)
!!!    !    return
!!!    !   end subroutine DiagTrajectoryINit
     !
     !    Subroutine  DiagUpdaterParticleFieldPeriod(GD,Ns,PB,FG)
     !    Implicit none
     !    Class(*),intent(inout)  :: GD
     !    Integer(4),intent(in) ::  Ns
     !    Type(ParticleBundle),intent(in) :: PB(0:Ns)
     !    Type(Field),intent(in) :: FG
     !    !Integer(4),parameter :: NsMax=2_4 
     !    Type(ParticleMomentumOne) ::  TempPMO
     !    Integer(4) :: i,j,Shift
     !    Real(8) :: HeatingRate(FG%Nx)
     !    Select Type (GD)
     !        Type is (Grid1D(*,*))
     !                   Shift=1
     !                   Do i=0, Ns
     !                                 Call WeightingParticleMomentum(PB(i),TempPMO)
     !                                 Call GD%Update(GD%Nx,TempPMO%RhoOne,Shift)
     !                                 Call GD%Update(GD%Nx,TempPMO%JxOne,Shift)
     !                                 Call GD%Update(GD%Nx,TempPMO%TOne,Shift)
     !                                 HeatingRate=TempPMO%JxOne*FG%Ex
     !                                 Call GD%Update(GD%Nx,HeatingRate,Shift)
     !                   End do
     !                   Call GD%Update(GD%Nx,FG%Ex,Shift)
     !                   Call GD%Update(GD%Nx,FG%Phi,Shift)
     !                   Call GD%Update(GD%Nx,FG%Rho,Shift)
     !                   GD%Timer=GD%Timer+1
     !           Type is (Grid2D(*,*,*))
     !                       Shift=1
     !                        Do i=0, Ns
     !                                     Call WeightingParticleMomentum(PB(i),TempPMO)
     !                                     Call GD%Update(GD%Nx,TempPMO%RhoOne,Shift)
     !                                     Call GD%Update(GD%Nx,TempPMO%JxOne,Shift)
     !                                     Call GD%Update(GD%Nx,TempPMO%TOne,Shift)
     !                                     HeatingRate=TempPMO%JxOne*FG%Ex
     !                                     Call GD%Update(GD%Nx,HeatingRate,Shift)
     !                       End do
     !                       Call GD%Update(GD%Nx,FG%Ex,Shift)
     !                       Call GD%Update(GD%Nx,FG%Phi,Shift)
     !                       Call GD%Update(GD%Nx,FG%Rho,Shift)
     !                       GD%Timer=GD%Timer+1
     !           End select
     !   Return    
     !End Subroutine  DiagUpdaterParticleFieldPeriod