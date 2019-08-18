!Module DiagnosticsElectrodeParticles
!    Use ModuleGrid
!    Use ModuleField
!    Use ModuleParticleBundle
!    Use DiagnosticsEEPF
!    Implicit none
!                   Integer(4),Parameter,Private :: NfMax=1000
!                    Type ParticleElectrode
!                        Integer(4) :: Timer=0,Period,NFlux=100              
!                        Real(8) :: Dt
!                        Type(ParticleBundle) :: PBLower,PBUpper
!                        Real(8) :: Flux(1:NfMax,1:2)=0.d0,EnergyFlux(1:NfMax,1:2)=0.d0
!                    EndType ParticleElectrode
!
!    contains
!             subroutine DiagParticleElectrode(PE,PB,FG,Mode)
!             Implicit none
!            Type(ParticleElectrode),intent(inout)  :: PE
!            Type(ParticleBundle),intent(in) :: PB
!            Type(Field),intent(in) :: FG
!            !Type(ParticleBoundaryOne),intent(in) :: PBDO
!            Integer(4),intent(in) ::  Mode
!            Type(ParticleBundle),save :: TempPB
!            Type(Grid1D(Nx=1000,Ns=2)),save :: GD_EDF
!            Type(Grid1D(Nx=737,Ns=4)),save :: GD_Flux
!            Integer(4) :: i,ParticleTimer
!            
!            !Allocate(Type(Grid1D(Nx=1000,Ns=2))::GD_Flux)
!            
!            Select Case (Mode)
!                 Case(-1)
!                        !PE%Timer=0
!                        !!PE%Period=FG%Period
!                        !!PE%NFlux=FG%Period
!                        !PE%Dt=FG%Dt
!                        !PE%PBLower=PB
!                        !PE%PBLower%NPar=0
!                        !PE%PBLower%Timer=0
!                        !PE%PBUpper=PB
!                        !PE%PBUpper%NPar=0
!                        !PE%PBUpper%Timer=0
!                        !PE%Flux=0.d0
!                        !PE%EnergyFLux=0.d0
!                        Case(0)
!                            TempPB=PB
!                            Call ParticleMove(TempPB,FG)
!                            ParticleTimer=Mod(PE%Timer,PE%Period)
!                            !do i=TempPB%NPar,1,-1
!                            !     If (TempPB%PO(i)%X<=TempPB%XMin) then
!                            !         TempPB%PO(i)%Ay=dble(ParticleTimer)
!                            !         Call AddParticle(TempPB%PO(i),PE%PBLower)
!                            !     else If(TempPB%PO(i)%X>=TempPB%XMax) then
!                            !         TempPB%PO(i)%Ay=dble(ParticleTimer)
!                            !         Call AddParticle(TempPB%PO(i),PE%PBUpper)
!                            !     end if
!                            !end do
!                            PE%Timer=PE%Timer+1
!                       Case(1)
!                                 !Allocate(Type(Grid1D(100,2))::GD_Flux)
!                                 Call PE%PBLower%Dump(1)
!                                 Call PE%PBUpper%Dump(1)
!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBLower,-1)
!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBLower,0)
!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBLower,1)
!                                 
!                                 
!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBUpper,-1)
!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBUpper,0)
!                                 Call DiagParticleEDFOne(GD_EDF,PE%PBUpper,1)
!                                 
!                                 Call DiagParticleFlux(GD_Flux,PE,-1)
!                                 Call DiagParticleFlux(GD_Flux,PE,0)
!                                 Call DiagParticleFlux(GD_Flux,PE,1)
!                                 !
!                                 !PE%Timer=0
!                                 !PE%PBUpper%NPar=0
!                                 !PE%PBUpper%Timer=0 
!                                 !PE%PBUpper%NPar=0
!                                 !PE%PBUpper%Timer=0 
!                         Case(2)
!                         case default
!                         End Select
!           return
!        end subroutine DiagParticleElectrode
!        
!        subroutine DiagParticleFlux(GD,PE,Mode)
!            Implicit none
!            Type(Grid1D(*,*))  :: GD
!            Type(ParticleElectrode),intent(inout)  :: PE
!            Integer(4),intent(in) ::  Mode
!            Integer(4) :: Shift,i
!                     Select Case (Mode)
!                        Case(-1)
!                            !Call GridInitialization(GD,PE%PBLower%Period,PE%PBLower%Dx,PE%PBLower%Dt)
!                        Case(0)
!                            Shift=1
!                            PE%NFlux=100 !GD%Nx
!                            Call WeightingParticleFlux(PE)
!                            Do i=1,2
!                                Call GD%Update(PE%Nflux,PE%Flux(1:PE%NFlux,i),Shift)
!                                Call GD%Update(PE%Nflux,PE%EnergyFlux(1:PE%NFlux,i),Shift)
!                            END DO
!                            GD%Timer=GD%Timer+1
!                         Case(1) 
!                                 Call GD%Rescale
!                                 Call GD%Dump(1)
!                                 GD%Value=0.d0
!                                 GD%Timer=0
!                         case default
!                         End Select
!            return
!        end subroutine DiagParticleFlux     
!        
!                subroutine WeightingParticleFlux(PE)
!                implicit none
!                Type(ParticleElectrode),intent(inout)  :: PE
!                Real(8) :: Energy
!                Integer(4) :: i,Timer,NAva
!
!                PE%Flux=0.d0
!                PE%EnergyFlux=0.d0
!                NAva=PE%Period/PE%NFlux
!                do i=1,PE%PBLower%Npar
!                   Call CalEnergy (PE%PBLower%Mass,PE%PBLower%PO(i),Energy)
!                   Energy=Energy*PE%PBLower%VFactor*PE%PBLower%VFactor/JtoeV
!                   !Write(*,*) Energy,"Lo"
!                   Timer=Int(PE%PBLower%PO(i)%Az)/NAva+1
!                   If(Timer>=1.and.Timer<=PE%NFLux) Then
!                       PE%Flux(Timer,1)=PE%Flux(Timer,1)+PE%PBLower%Weight
!                       PE%EnergyFlux(Timer,1)=PE%EnergyFlux(Timer,1)+PE%PBLower%Weight*Energy
!                   End IF
!                End do
!                
!                do i=1,PE%PBUpper%Npar
!                   Call CalEnergy (PE%PBUpper%Mass,PE%PBUpper%PO(i),Energy)
!                   Energy=Energy*PE%PBUpper%VFactor*PE%PBUpper%VFactor/JtoeV
!                   !Write(*,*) Energy,"Up"
!                   Timer=Int(PE%PBUpper%PO(i)%Az)/NAva+1
!                   If(Timer>=1.and.Timer<=PE%NFLux) Then
!                       PE%Flux(Timer,2)=PE%Flux(Timer,2)+PE%PBUpper%Weight
!                       PE%EnergyFlux(Timer,2)=PE%EnergyFlux(Timer,2)+PE%PBUpper%Weight*Energy
!                   End IF
!                End do
!                PE%Flux=PE%Flux/(dble(NAva)*PE%PBLower%dt)
!                PE%EnergyFlux=PE%EnergyFlux/(dble(NAva)*PE%PBLower%dt)
!                return
!        end subroutine WeightingParticleFlux
!End Module DiagnosticsElectrodeParticles