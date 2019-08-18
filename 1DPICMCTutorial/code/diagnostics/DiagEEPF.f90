Module DiagnosticsEEPF
    Use ModuleGrid
    Use ModuleParticleBundle
    Implicit none

                   Integer(4),Parameter,Private :: NeMax=1000
                    Type ParticleEDF
                        Integer(4) :: Ne=NeMax
                        Real(8) :: EnergyInterval,Mass
                        Real(8) :: EDF(NeMax),EDFNormalized(NeMax)
                    EndType ParticleEDF
    contains
             subroutine DiagParticleEDFOne(GD,PB,Mode)
            Implicit none
            Class(*),intent(inout)  :: GD
            Type(ParticleBundle),intent(in) :: PB
            Integer(4),intent(in) ::  Mode
            Type(ParticleEDF) :: PEDF
            Integer(4) :: Shift

            If (PB%SO%SpecyIndex==0) Then
                PEDF%EnergyInterval=0.1d0
            Else
                PEDF%EnergyInterval=1.d0
            End If
            PEDF%Mass=PB%Mass

            Select Type (GD)
                  Type is (Grid1D(*,*))
                     Select Case (Mode)
                        Case(-1)
                            !Call GridInitialization(GD,PB%Period,PEDF%EnergyInterval,PB%Dt)
                        case(0)
                            
                            PEDF%Ne=GD%Nx
                            Shift=1
                            Call WeightingParticleEDF(PB,PEDF)
                            Call GD%Update(GD%Nx,PEDF%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
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
                            Case(-1)
                                !Call GridInitialization(GD,PB%Period,PEDF%EnergyInterval,PB%Dt)
                            case(0)
                                Shift=1
                                Call WeightingParticleEDF(PB,PEDF)
                                Call GD%Update(GD%Nx,PEDF%EDF,Shift)
                                Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
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
            return
             end subroutine DiagParticleEDFOne
             
        subroutine WeightingParticleEDF(PB,PEDF)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleEDF),intent(inout) :: PEDF
                Real(8) :: Energy,Frac
                Integer(4) :: i,N
                
                PEDF%EDF=0.d0
                PEDF%EDFNormalized=0.d0
                
                Frac=1.d0/DBLE(PB%Npar)
                do i=1,PB%Npar
                   Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
                   N=Ceiling( Energy/PEDF%EnergyInterval)
                   If (N>=1.and.N<PEDF%Ne) Then
                           !PEDF%EDF(N)=PEDF%EDF(N)+PB%Weight/PB%XMax
                           PEDF%EDFNormalized(N)=PEDF%EDFNormalized(N)+Frac
                   End IF
                end do
                return
    end subroutine WeightingParticleEDF     
End Module DiagnosticsEEPF