Module DiagnosticsTestParticle
    Use ModuleGrid
    Use ModuleField
    Use ModuleParticleBundle
    Use ModuleParticleBoundary
    Implicit none
       
    contains
         subroutine DiagParticleTestParticle(PB,PBDO,FG,Mode)
             Implicit none
            Type(ParticleBundle),intent(in) :: PB
            Type(ParticleBoundaryOne),intent(inout) :: PBDO
            Type(Field),intent(in) :: FG
            Integer(4),intent(in) ::  Mode
            Type(ParticleBundle),save :: TestPB,TempPB
            Type(ParticleOne),save :: TempPO
            Integer(4) :: i
            Real(8) :: Energy
            Select Case (Mode)
                       Case(-1)
                       TempPB=PB
                       TempPB%Npar=1
                       !TempPB%Timer=0
                       TempPB%PO(1)%X=1.d0
                       !!!这里修改粒子的初始位置，注意一次只能跟踪一个粒子
                       !!!粒子位置和速度都是归一化后的，记住，需要自己算一下归一化
                       TempPB%PO(1)%Vx=0.d0
                       TempPB%PO(1)%Vy=0.d0
                       TempPB%PO(1)%Vz=0.d0
                       TempPB%PO(1)%Ax=0.d0
                       TempPB%PO(1)%Ay=0.d0
                       TempPB%PO(1)%Az=0.d0
                       TestPB=TempPB
                       !TestPB%Timer=1
 
                       Case(0)
                           !If (TempPB%Timer<=NParMax) Then
                          TempPO=TempPB%PO(1)
                          If (TempPO%X<=PBDO%XMin.or.TempPO%X>=PBDO%XMax) Then
                               Call TestPB%Dump(3)
                               Stop
                          Else!Call ParticleMove(TempPB,FG)
                               Call TempPB%MoveEM(FG)
                               TempPO=TempPB%PO(1)
                               Energy=TempPO%Energy(TempPB%Mass,TempPB%Vfactor)
                               TempPO%Az=Energy/JtoeV
                              Call TestPB%Addone(TempPO)
                              Write(*,*) TestPB%Npar
                          End If

                               !Call CalEnergy (TestPB%Mass,TempPO,Energy)
                               !Energy=Energy*PB%VFactor*PB%VFactor/JtoeV
                              
                                !End If
                               !TempPB%Timer=TempPB%Timer+1
                               !TestPB%Timer=TestPB%Timer+1
                           !End IF
                        Case(1)
                            Call TestPB%Dump(4)
                        Case(2)
                         case default
                         End Select
           return
        end subroutine DiagParticleTestParticle  
  End Module DiagnosticsTestParticle