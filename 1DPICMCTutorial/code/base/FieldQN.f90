Module ModuleFieldQN
     Use ModuleControlFlow
     Use ModuleVector1DX
     Use ModuleParticleBundle
     Use ModuleField
     Implicit none
                  !Type FieldQN !(Nx)\
                  !      Integer(4) :: Nsheath=0,NL=0,NR=NxMax
                  !      Integer(4) :: Nx=NxMax
                  !      Real(8) :: Dx=Inputdx,Dt=Inputdt
                  !      Real(8) ::  Ex(1:NxMax)=0.d0,Ey(1:NxMax)=0.d0,Ez(1:NxMax)=0.d0
                  !      Real(8) ::  Bx(1:NxMax)=0.d0,By(1:NxMax)=0.d0,Bz(1:NxMax)=0.d0
                  !      Real(8) ::  Rho(1:NxMax),Phi(1:NxMax)
                  !      Real(8) ::  Chi(1:NxMax)
                  !      contains
                  !           procedure :: Dump=>DumpField
                  !           procedure :: Load=>LoadField
                  ! EndType Field
                   
                    Type FieldSolverQNSheath
                        Integer(4) :: Nx=NxMax,Nxs=0,Nsheath=0
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        Real(8),Allocatable :: Source(:)
                        Real(8),Allocatable :: Solve(:)
                        Real(8),Allocatable :: CoeA(:),CoeB(:),CoeC(:)
                        contains
                             procedure :: Init=>InitFieldSolverQNSheath
                             procedure :: UpdateNsheath=>UpdateNsheathFieldSolverQNSheath
                             !procedure :: UpdateSource=>UpdateSourceFieldSolverQNSheath
                        EndType FieldSolverQNSheath
                        
                    Type FieldSolverQNBulk
                        Integer(4) :: Nx=NxMax,NL=1,NR=NxMax,Nbulk=0
                        !,NL=0,NR=NxMax
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        Real(8) :: Ni(1:NxMax)=0.d0,Vi(1:NxMax)=0.d0,Te(1:NxMax)=0.d0
                        Real(8) :: Vbohm(1:NxMax)=0.d0
                        Real(8) :: Eb(1:NxMax)=0.d0
                        !Real(8),Allocatable :: Eb(:)
                    contains
                          !procedure :: Init=>InitFieldSolverQNBulk
                          !procedure :: Update=>UpdateOneStepFieldSolverQNBulk
                          procedure :: UpdateE=>UpdateEFieldFieldSolverQNBulk
                          procedure :: UpdatenT=>UpdatenTFieldSolverQNBulk
                        EndType FieldSolverQNBulk
                        
                    Type ParticleMomentumOne!(Nx)
                        !Integer(4),Len :: Nx=100
                        !Integer(4) :: IOIndex=1
                        !Integer(4) :: XStart=0,XEnd=NxMax-1
                        Integer(4) :: Nx=NxMax
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        Real(8) :: Vbohm(1:NxMax)=0.d0
                        !Real(8) :: QdM
                        Real(8) ::  N(1:NxMax),V(1:NxMax),T(1:NxMax)
                   EndType ParticleMomentumOne

 
    contains
    
    
    
    subroutine InitFieldSolverQNSheath(FSQNS)
                    Implicit none
                    Class(FieldSolverQNSheath),intent(inout) :: FSQNS
                    If(allocated(FSQNS%Source)) Then
                       Deallocate(FSQNS%Source)
                       Deallocate(FSQNS%Solve)
                       Deallocate(FSQNS%CoeA)
                       Deallocate(FSQNS%CoeB)
                       Deallocate(FSQNS%CoeC)
                    End If
                    !
                    !If (FSQNS%Nsheath/=0) Then
                    !FSQNS%Nxs=2*FSQNS%Nsheath+1
                    Allocate(FSQNS%Source(FSQNS%Nxs))
                    Allocate(FSQNS%Solve(FSQNS%Nxs))
                    Allocate(FSQNS%CoeA(FSQNS%Nxs)) 
                    Allocate(FSQNS%CoeB(FSQNS%Nxs))
                    Allocate(FSQNS%CoeC(FSQNS%Nxs))
                    !End If!
                    
                    return
    end subroutine InitFieldSolverQNSheath
    
        !subroutine InitFieldSolverQNBulk(FSQNB)
        !            Implicit none
        !            Class(FieldSolverQNBulk),intent(inout) :: FSQNB
        !            If(allocated(FSQNB%Eb)) Then
        !               Deallocate(FSQNB%EB)
        !            End If
        !            Allocate(FSQNB%Eb(FSQNB%NL:FSQNB%NR))
        !            return
        !end subroutine InitFieldSolverQNBulk
    
            subroutine UpdateNsheathFieldSolverQNSheath(FSQNS,FSQNB)
                    Implicit none
                    Class(FieldSolverQNSheath),intent(inout) :: FSQNS
                    Type(FieldSolverQNBulk),intent(inout) :: FSQNB
                    Integer(4) :: i,NSheathMax=100
                    
                    FSQNS%Nsheath=NSheathMax
                    FSQNS%Nxs=2*(FSQNS%Nsheath-1)+1
                    FSQNB%NL=FSQNS%Nsheath+1 
                    FSQNB%NBulk=FSQNB%Nx-2-2*FSQNS%Nsheath
                    FSQNB%NR=FSQNB%NL+FSQNB%NBulk
                    
                    Call FSQNS%Init
                    return
            end subroutine UpdateNsheathFieldSolverQNSheath
            
        !subroutine UpdateNsheathFieldSolverQNSheath(FSQNS,FS)
        !            Implicit none
        !            Class(FieldSolverQNSheath),intent(inout) :: FSQNS
        !            Type(FieldSolverQNBulk),intent(in) :: FSQNB
        !            Integer(4) :: i,NSheathNew=0,NSheathMax=0
        !            !Logical@:: 
        !            
        !            NSheathMax=(FSQNS%Nx+1)/2                    
        !            NSheathNew=0
        !            
        !            
        !
        !            Do i=NSheathMax,1
        !               If(ABS(FSQNB%Vi(i))>ABS(FSQNB%Vbohm(i))) Then
        !                 NSheathNew=NSheathMax-i
        !                 Exit
        !               ENd if
        !            ENd do
        !            
        !            If (NSheathNew>FSQNS%Nsheath) Then
        !               FSQNS%Nsheath=NSheathNew
        !            End If
        !            
        !            Call FSQNS%Init
        !            return
        !            
        !End subroutine UpdateNsheathFieldSolverQNSheath
            

        
        subroutine UpdateEFieldFieldSolverQNBulk(FSQNB)
                    Implicit none
                    Class(FieldSolverQNBulk),intent(inout) :: FSQNB
                    Real(8) :: i,nTn(1:FSQNB%Nx),nTc(1:FSQNB%Nx+1),GradY(1:FSQNB%Nx),GradZ(1:FSQNB%Nx)
                    
                    nTn=FSQNB%Ni*FSQNB%Te

                    Call Interpolation1DX(nTc, nTn)
                    
                    Call Grad1DX(nTn, GradY, GradZ, nTc, FSQNB%dx)
                    
                    FSQNB%Eb=-1.d0*nTn/FSQNB%Ni/ElectronCharge
                    
                    !Do i=FSQNB%NL,FSQNB%NR
                    !   
                    !End DO
                    return
        end subroutine UpdateEFieldFieldSolverQNBulk

        subroutine UpdatenTFieldSolverQNBulk(FSQNB,Ns,PB)
                    Implicit none
                    Class(FieldSolverQNBulk),intent(inout) :: FSQNB
                    Integer(4),intent(in) :: Ns
                    Type(ParticleBundle),intent(in) :: PB(0:Ns)
                    Type(ParticleMomentumOne):: PMO(0:Ns)
                    Integer(4) :: i
                    
                    FSQNB%Te=0.d0
                    FSQNB%Ni=0.d0
                    FSQNB%Vi=0.d0
                    
                    Do i=0,Ns
                       Call WeightingParticleMomentum(PB(i),PMO(i))
                    End Do
                    
                    FSQNB%Te=PMO(0)%T
                    
                    Do i=1,Ns
                       PMO(i)%Vbohm=DSQRT(PMO(0)%T/PB(i)%Mass)
                       FSQNB%Ni=FSQNB%Ni+PMO(i)%N
                       FSQNB%Vi=FSQNB%Vi+PMO(i)%V
                       FSQNB%Vbohm= FSQNB%Vbohm+PMO(i)%Vbohm
                    End Do
                   
                    FSQNB%Vi=FSQNB%Vi/dble(Ns)
                    FSQNB%Vbohm=FSQNB%Vbohm/dble(Ns)
 
                    return
        end subroutine UpdatenTFieldSolverQNBulk
        
        

        
        subroutine WeightingParticleMomentum(PB,PMO)
            Implicit none
            Type(ParticleBundle),intent(in) :: PB
            Type(ParticleMomentumOne),intent(inout) :: PMO
                Real(8) :: NFactor,VFactor,TFactor
                Real(8) :: S1,S2,Energy
                Integer(4) :: i,N
                PMO%N=0.d0
                PMO%V=0.d0
                PMO%T=0.d0
                do i=1,PB%Npar
                   N=Ceiling(PB%PO(i)%X)

                   S1=Dble(N)-PB%PO(i)%X
                   S2=1.d0-S1
                   PMO%N(N)=PMO%N(N)+S1
                   PMO%N(N+1)=PMO%N(N+1)+S2
                   
                   PMO%V(N)=PMO%V(N)+S1*PB%PO(i)%Vx
                   PMO%V(N+1)=PMO%V(N+1)+S2*PB%PO(i)%Vx

                   !Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)
                   !Energy=0.5d0*PB%Mass*(PB%PO(i)%Vx*PB%PO(i)%Vx)*PB%VFactor*PB%VFactor
                   Energy=PB%Mass*(PB%PO(i)%Vx*PB%PO(i)%Vx)*PB%VFactor*PB%VFactor
                   PMO%T(N)=PMO%T(N)+S1*Energy
                   PMO%T(N+1)=PMO%T(N+1)+S2*Energy
                end do
                
                NFactor=PB%Weight
                PMO%N=PMO%N*NFactor
                VFactor=PB%Weight*PB%VFactor
                PMO%V=PMO%V*VFactor
                TFactor=PB%Weight
                PMO%T=PMO%T*TFactor
                do i=1,PMO%Nx
                    if (PMO%N(i)>0.d0) then
                        PMO%V(i)=PMO%V(i)/PMO%N(i)
                        PMO%T(i)=PMO%T(i)/PMO%N(i)
                    end if
                end do
                PMO%N(1)=2.d0*PMO%N(1)
                PMO%N(PMO%Nx)=2.d0*PMO%N(PMO%Nx)
                
                !PMO%V(1)=2.d0*PMO%V(1)
                !PMO%V(PMO%Nx)=2.d0*PMO%V(PMO%Nx)
                !PMO%T(1)=2.d0*PMO%T(1)
                !PMO%T(PMO%Nx)=2.d0*PMO%T(PMO%Nx)
                !PMO%ChiOne(1)=2.d0*PMO%ChiOne(1)
                !PMO%ChiOne(PMO%Nx)=2.d0*PMO%ChiOne(PMO%Nx)
                return
  end subroutine WeightingParticleMomentum
        
        
        
        
        
        
    
        !subroutine InitFieldSolverQNBulk(FSQNB,FSQNS)
        !            Implicit none
        !            Class(FieldSolverQNBulk) ::
        !            Type(FieldSolverQNSheath),intent(inout) :: FS
        !            Integer(4) :: Mode
        !            Character(len=99) :: Filename
        !            Integer(4):: i,j,NameIndex
        !            Integer(4),save :: NameTimer=1
        !            If (Mode==0) Then
        !                       NameIndex=DefaultNameIndex
        !            else
        !                       NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
        !                       NameTimer=NameTimer+1
        !            End If 
        !     
        !            Write(filename,*) NameIndex,"FieldSolver",".dat"
        !            Write(*,*) "Saving ",filename," Please wait..."
        !            open (10,file=filename)
        !            !Write(10,*)  FS%Ns,FS%Dx,FS%Dt
        !            do i=1,FS%Ns
        !                   Write(10,FMt="(*(es21.14,1x))") dble(i-1)*FS%Dx,FS%Solve(i),FS%Source(i),FS%CoeA(i),FS%CoeB(i),FS%CoeC(i)!,FG%Bx(i),FG%By(i)
        !            end do
        !            close(10)
        !            Write(*,*) "Save ",filename,"Complete!"  
        !            return
        !     end subroutine DumpFieldSolver

             !subroutine DumpFieldSolver(FS,Mode)
             !       Implicit none
             !       Class(FieldSolver),intent(inout) :: FS
             !       Integer(4) :: Mode
             !       Character(len=99) :: Filename
             !       Integer(4):: i,j,NameIndex
             !       Integer(4),save :: NameTimer=1
             !       If (Mode==0) Then
             !                  NameIndex=DefaultNameIndex
             !       else
             !                  NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
             !                  NameTimer=NameTimer+1
             !       End If 
             !
             !       Write(filename,*) NameIndex,"FieldSolver",".dat"
             !       Write(*,*) "Saving ",filename," Please wait..."
             !       open (10,file=filename)
             !       !Write(10,*)  FS%Ns,FS%Dx,FS%Dt
             !       do i=1,FS%Ns
             !              Write(10,FMt="(*(es21.14,1x))") dble(i-1)*FS%Dx,FS%Solve(i),FS%Source(i),FS%CoeA(i),FS%CoeB(i),FS%CoeC(i)!,FG%Bx(i),FG%By(i)
             !       end do
             !       close(10)
             !       Write(*,*) "Save ",filename,"Complete!"  
             !       return
             !end subroutine DumpFieldSolver
             !
             ! subroutine DumpFieldOne(Ns,FO,Mode)
             !       Implicit none
             !       Integer(4),intent(in) :: Ns
             !       Type(FieldOne),intent(inout) :: FO(0:Ns)
             !       Integer(4) :: Mode
             !       Character(len=99) :: Filename
             !       Integer(4):: i,j,NameIndex
             !       Integer(4),save :: NameTimer=1
             !       If (Mode==0) Then
             !                  NameIndex=DefaultNameIndex
             !       else
             !                  NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
             !                  NameTimer=NameTimer+1
             !       End If 
             !
             !       Write(filename,*) NameIndex,"FieldOne",".dat"
             !       Write(*,*) "Saving ",filename," Please wait..."
             !       open (10,file=filename)
             !       !Write(10,*)  FO%Nx,FO%Dx,FO%Dt
             !       do i=1,FO(0)%Nx
             !              Write(10,FMt="(*(es21.14,1x))") FO(0)%Dx*dble(i-1),(ABS(FO(j)%N(i)/ElectronCharge),j=0,Ns)
             !       end do
             !       close(10)
             !       Write(*,*) "Save ",filename,"Complete!"  
             !       return
             ! end subroutine DumpFieldOne    
           
 End Module ModuleFieldQN