Module ModuleOneStep
      Use ModuleParticleBundle
      Use ModuleSpecyOne
      Use ModuleMCCPublic
      Use ModuleOneStepField
      Use ModuleMCCInitialization
      
      !Use MoveModule
!      Use Diagnostics
      Implicit none
      ! This section defines the particles.
                  Type(ControlFlow) :: ControlFlowGlobal
      
                  Type(ParticleBundle), save, Allocatable       ::  ParticleLocal(:) [:]
                  Type(ParticleBoundaryOne), save, Allocatable  ::  ParticleBDOneLocal(:) [:]
                  !Type(ParticleBoundary),save  :: ParticleBDGlobal
                
    contains
    Subroutine AllInitilalization()
             Implicit none
             Integer(4) :: i
             Call InitializationControlFlow(ControlFlowGlobal)
             Call GasInit(ControlFlowGlobal)
             Call InitializationField(ControlFlowGlobal)
             Allocate(ParticleLocal(0:ControlFlowGlobal%Ns)[*])
		Allocate(ParticleBDOneLocal(0:ControlFlowGlobal%Ns)[*])
             DO i=0,ControlFlowGlobal%Ns
                    Call ParticleLocal(i)%AllInit(SpecyGlobal(i),ControlFlowGlobal)
                    Call ParticleBDOneLocal(i)%AllInit(ParticleLocal(i),ControlFlowGlobal)
             End do
             Call MCCBundleInit(ControlFlowGlobal,SpecyGlobal,GasGlobal)
            Return  
    End Subroutine AllInitilalization

    Subroutine OneStep()
             !  Use ,only : R
               Implicit none
               !Integer(4) :: Ntemp                   !Ntemp=ParticleGlobal(0)%NPar
            !   Real(8)::GamaE=0.2d0,GamaI=0.2d0
                                    !Write(*,*) ParticleBDOneGlobal%CountMinOne,ParticleBDOneGlobal%CountMin,ParticleGlobal%NPar
               Integer(4) :: i,j
               do i=0,ControlFlowGlobal%Ns
                  !    Call ParticleLocal(i)%MoveEM(FieldGlobal)
                     Call ParticleLocal(i)%MoveES(FieldGlobal)
                     !Call ParticleMove(ParticleGlobal(i),FieldGlobal)
                     Call ParticleAborption(ParticleLocal(i),ParticleBDOneLocal(i))
                     !Call Selectron(ParticleGlobal(0),ParticleBDOneGlobal(i))
                     !Call WeightingOne(ParticleGlobal(i),FieldOneGlobal(i))
                     Call ParticleLocal(i)%WeightP2C(FieldOneLocal(i))
               end do
               sync all

               ! 合并场
			do j=0,ControlFlowGlobal%Ns
				FieldOneGlobal(j)%RhoOne = 0.d0
				FieldOneGlobal(j)%ChiOne = 0.d0
				do i = 1, imageSize
					FieldOneGlobal(j)%RhoOne = FieldOneGlobal(j)%RhoOne + FieldOneLocal(j)[i]%RhoOne
					FieldOneGlobal(j)%ChiOne = FieldOneGlobal(j)%ChiOne + FieldOneLocal(j)[i]%ChiOne
				end do
			end do

			! ! 合并上下边界粒子数
			! ParticleBDOneLocal%CountMinGlobal = 0
			! ParticleBDOneLocal%CountMaxGlobal = 0
			
			! do i = 1, imageSize
			! 	ParticleBDOneLocal%CountMinGlobal = ParticleBDOneLocal%CountMinGlobal + ParticleBDOneLocal[i]%CountMin
			! 	ParticleBDOneLocal%CountMaxGlobal = ParticleBDOneLocal%CountMaxGlobal + ParticleBDOneLocal[i]%Countmax
			! end do

			sync all
                  
               Call FieldOneStep(ControlFlowGlobal%Ns,FieldOneGlobal,FieldGlobal,FieldBoundaryGlobal,FieldSolverGlobal)
               !Write(*,*) "Period Before",ParticleGlobal%NPar
               Call MCC(ControlFlowGlobal%Ns,ControlFlowGlobal%Ng,ParticleLocal,SpecyGlobal,GasGlobal,MCCBundleGlobal) 

               sync all

               ! 合并总粒子数
                  ParticleLocal%NParAll = 0

                  do i = 1, imageSize
                        ParticleLocal%NParAll = ParticleLocal%NParAll + ParticleLocal[i]%NPar
                  end do

                  sync all

               !Write(*,*) "Period After",ParticleGlobal%NPar
              return
    End  subroutine OneStep
    
        Subroutine OneStepRestart()
               !Use FileIO 
               Implicit none
               Integer(4) :: i
               !Write(*,*) Period,ParticleGlobal%NPar,"Period"!,FieldBoundaryGlobal%Qmin,FieldBoundaryGlobal%Qmax
               Call DumpFieldSolver(FieldSolverGlobal)
               Call DumpField(FieldGlobal)
               Call DumpFieldOne(ControlFlowGlobal%Ns,FieldOneGlobal)
               Call DumpFieldOne(ControlFlowGlobal%Ns,FieldOneLocal,imageRank)
               do i=0,ControlFlowGlobal%Ns
                     Call ParticleLocal(i)%Dump()
               End do
               Call FieldBoundayFinalization(FieldBoundaryGlobal)
              return
        End  subroutine OneStepRestart

        Subroutine MergeAndSplit()
            Implicit none
          integer(4) :: m, k
          
                ! 粒子合并与分裂
                DO m=0,ControlFlowGlobal%Ns
                      If (ParticleLocal(m)%NParAll > int(1.1 * ParticleLocal(m)%ParticleNumberAll)) then

                            sync all

                            if (1 == imageRank) then
                                  Write(*,*) "before:",(ParticleLocal(k)%NParAll,k=0,ControlFlowGlobal%Ns)
                                  open (10,file='1 Weightingchange.dat',position='APPEND') 
                                  Write(10,FMt="(*(es21.14,1x))") ControlFlowGlobal%Timer*ControlFlowGlobal%Dt,(ParticleLocal(k)%Weight*ControlFlowGlobal%Dx,k=0,ControlFlowGlobal%Ns)
                                  close(10)
                            end if
                            
                  ! 聚变所集群上，main中好像不能使用下面的循环。。。
                  ! Call ParticleBundleNormalizationMergeAll(ControlFlowGlobal, ParticleLocal)
                            do k=0,ControlFlowGlobal%Ns
                                  if (ParticleLocal(k)%NParAll > 0) then
                                        Call ParticleBundleNormalization(ParticleLocal(k), 0.75d0)
                                  end if
                            end do

                            if (1 == imageRank) then
                                  open (10,file='1 Weightingchange.dat',position='APPEND') 
                                  Write(10,FMt="(*(es21.14,1x))") ControlFlowGlobal%Timer*ControlFlowGlobal%Dt,(ParticleLocal(k)%Weight*ControlFlowGlobal%Dx,k=0,ControlFlowGlobal%Ns)
                                  close(10)
                                  Write(*,*) "after:", (ParticleLocal(k)%NParAll,k=0,ControlFlowGlobal%Ns)
                            end if

                            sync all
                      End If
                end do

          return
      End Subroutine MergeAndSplit

    End Module ModuleOneStep
    
    
       

    !

    !
    !

        
        
        
        

