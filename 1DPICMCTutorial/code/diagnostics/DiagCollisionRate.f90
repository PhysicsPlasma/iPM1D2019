!!
    Module DiagnosticsCollisionRate
    Use ModuleGrid
    Use ModuleParticleBundle
    Use ModuleTypeMCC
    Implicit none
    Integer(4),Parameter  :: NReactionMax=100
                    Type ParticleCollisionRate
                        Integer(4) :: Nx=NxMax
                        Integer(4) :: NReaction
                        Real(8) :: CollisionRate(NxMax,NReactionMax)
                    EndType ParticleCollisionRate
    Type(ParticleCollisionRate), save, allocatable ::  TempPCRGlobal [:]
    Type(ParticleCollisionRate), save, allocatable ::  TempPCRLocal [:]
    contains
    subroutine DiagParticleCollisionRateInitilalization()
        allocate (TempPCRGlobal [*])
        allocate (TempPCRLocal [*])
    end subroutine DiagParticleCollisionRateInitilalization
    subroutine DiagParticleCollisionRateOne(GD,PB,MCCB,Mode)
            Implicit none
            Class(*),intent(inout)  :: GD
            Type(ParticleBundle),intent(in) :: PB
            Type(MCCBundle),intent(in) ::  MCCB
            Integer(4),intent(in) ::  Mode
            ! Type(ParticleCollisionRate) :: PCR
            Integer(4) :: Shift,i,k
            Select Type (GD)
                  Type is (Grid1D(*,*))
                     Select Case (Mode)
                        Case(-1)
                            !Call GD%Init(CF)
                        case(0)
                            Shift=1
                            TempPCRGlobal%CollisionRate=0.d0
                            Call WeightingParticleCollisionRate(PB,TempPCRLocal,MCCB)
                            sync all
                            do k=1,imageSize
                                TempPCRGlobal%CollisionRate=TempPCRGlobal%CollisionRate+TempPCRLocal[k]%CollisionRate
                            end do
                            sync all
                            do i=1,MCCB%NReaction
                                 Call GD%Update(GD%Nx,TempPCRGlobal%CollisionRate(:,i),Shift)
                            End Do
                            GD%Timer=GD%Timer+1
                         Case(1) 
                                 Call GD%Rescale()
                                 Call GD%Dump(1)
                                 Call GD%Reset()
                         Case(2)
                             Call GD%Dump(0)
                         case default
                         End Select
                    Type is (Grid2D(*,*,*))
                        Select Case (Mode)
                            Case(-1)
                                !Call GridInitialization(GD,PB%Period,PB%Dx,PB%Dt)
                            case(0)
                                Shift=1
                                TempPCRGlobal%CollisionRate=0.d0
                                Call WeightingParticleCollisionRate(PB,TempPCRLocal,MCCB)
                                sync all
                                do k=1,imageSize
                                    TempPCRGlobal%CollisionRate=TempPCRGlobal%CollisionRate+TempPCRLocal[k]%CollisionRate
                                end do
                                sync all
                                do i=1,MCCB%NReaction
                                    Call GD%Update(GD%Nx,TempPCRGlobal%CollisionRate(:,i),Shift)
                                End Do
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
            return !sort()
    end subroutine DiagParticleCollisionRateOne
    
       subroutine WeightingParticleCollisionRate(PB,PCR,MCCB)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleCollisionRate),intent(inout) :: PCR
                Type(MCCBundle),intent(in) ::  MCCB
                Real(8) :: EnergyInterval,Energy,S1,S2
                Real(8) :: TempProbility(MCCB%NReaction)
                Integer(4) :: i,j,Nc,Np,Index,Upper,Center,Lower,NReaction

                PCR%CollisionRate=0.d0
                NReaction=MCCB%NReaction
                PCR%NReaction=NReaction
                do j=1,PB%Npar
                    
                   Energy=PB%PO(j)%Energy(PB%Mass,PB%VFactor)/JtoeV
             Associate(Emin=>MCCB%EnergyMin,Eint=>MCCB%EnergyInterval,Nr=>MCCB%NReaction,Ns=>MCCB%NSigma,Probility=>MCCB%Probility,EnergyMax=>MCCB%EnergyMax)
                TempProbility=0.d0
                Nc=Int((Log10(Energy)-Emin)/Eint)
                If (Nc<1) Then
                    TempProbility=Probility(1:Nr,1)
                ELse If (Nc<Ns) Then
                     S2=(Log10(Energy)-Emin)/Eint-Dble(Nc)
                     S1=1.d0-S2
                     TempProbility(1:Nr)=Probility(1:Nr,Nc)*S1+Probility(1:Nr,Nc+1)*S2
                Else
                     TempProbility(1:Nr)=Probility(1:Nr,Ns)*Energy/EnergyMax
                End If
                
                   TempProbility=TempProbility*MCCB%SigmaMax*MCCB%GO%Density
                   Do i=NReaction,2
                         If ((TempProbility(i)-TempProbility(i-1))>=MinReal) Then
                              TempProbility(i)=TempProbility(i)-TempProbility(i-1)
                         End If
                   End do
            End Associate  

                   Np=Ceiling(PB%PO(j)%X)
                   Do  i=1, NReaction
                          PCR%CollisionRate(Np,i)=PCR%CollisionRate(Np,i)+TempProbility(i)
                   End DO
                end do
                return
        end subroutine WeightingParticleCollisionRate

 End Module DiagnosticsCollisionRate 