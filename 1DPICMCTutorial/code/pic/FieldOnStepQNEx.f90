Module ModuleOneStepFieldQNEx
  Use ModuleField
  Use ModuleFieldBoundary
  Use Numrical
  Use ModuleFieldQN
  Use ModuleOneStepFieldEx

     
     Type(FieldSolverQNSheath) :: FieldSolverQNSheathGlobal
     Type(FieldSolverQNBulk) :: FieldSolverQNBulkGlobal
                  
    contains
         
    
         subroutine FieldOneStepQNEx(Ns,FO,FG,FB,FS,FSQNS,FSQNB,PB)
            Implicit none
            Integer(4),intent(in) :: Ns
            Type(FieldOne),intent(in) :: FO(0:Ns)
            Type(Field),intent(inout) :: FG
            Type(FieldSolver), intent(inout):: FS
            Type(FieldBoundary), intent(inout):: FB
            Type(FieldSolverQNSheath), intent(inout):: FSQNS
            Type(FieldSolverQNBulk),intent(inout) :: FSQNB
            !Integer(4),intent(in) :: Np
            Type(ParticleBundle),intent(in) :: PB(0:Ns)
            
            Call AccumulationField(FG,Ns,FO)
            Call FB%Updater
            Call UpdaterCoeFieldSolver(FS,FG,FB)
            Call UpdaterFieldSolverQNS(FSQNS,FSQNB,FS)
            !Call UpdaterFieldSolver(FS)
            Call UpdaterField(FG,FS,FB)
            Call UpdateFieldSolverQNBulk(FSQNB,Ns,PB,FG)
            
            !Call UpdaterFieldSolver(FS)
            
            return
         End subroutine FieldOneStepQNEx
         
        
        !Call UpdaterField(FG,FS,FB) 
         
        
         

        !subroutine AccumulationField(FG,Ns,FO)
        !    Implicit none
        !    Type(Field),intent(inout) :: FG
        !    Integer(4),intent(in) :: Ns
        !    Type(FieldOne),intent(in) :: FO(0:Ns)
        !    Integer(4) :: i!,Nx
        !    FG%Rho=0.d0
        !    FG%Chi=0.d0
        !    !Nx=FG%Nx
        !            do i=0,Ns
        !                    FG%Rho=FG%Rho+FO(i)%RhoOne
        !                    FG%Chi=FG%Chi+FO(i)%ChiOne
        !            End do
        !    return
        !End subroutine AccumulationField
         
        subroutine UpdateFieldSolverQNBulk(FSQNB,Ns,PB,FG)
                    Implicit none
                    Class(FieldSolverQNBulk),intent(inout) :: FSQNB
                    Integer(4),intent(in) :: Ns
                    Type(ParticleBundle),intent(in) :: PB(0:Ns)
                    Type(Field),intent(inout) :: FG
                    Integer(4) :: i
                    Call FSQNB%UpdatenT(Ns,PB)
                    Call FSQNB%UpdateE
                    DO i=FSQNB%NL-1,FSQNB%NR+1
                    !DO i=FSQNB%NL,FSQNB%NR
                    !DO i=101,156
                       FG%Ex(i)=FSQNB%Eb(i)
                    ENd DO
                    
                    return
        end subroutine UpdateFieldSolverQNBulk
        
        subroutine UpdaterFieldSolverQNS(FSQNS,FSQNB,FS)
            Implicit none
           Class(FieldSolverQNSheath), intent(inout):: FSQNS
           Type(FieldSolverQNBulk),intent(inout) :: FSQNB
           Type(FieldSolver), intent(inout):: FS
           Integer(4) :: i,Nsheath!,NL,NR
           
           Call FSQNS%UpdateNsheath(FSQNB)
           Nsheath=FSQNS%Nsheath
           !NL=FSQNS%NL
           !NR=FSQNS%NR
           !Fs%Nxs=2*Nsheath+1
           
           !FSQNS%CoeA(Nsheath)=FS%CoeA(FS%Ns/2)
           !FSQNS%CoeB(Nsheath)=FS%CoeB(FS%Ns/2)
           !FSQNS%CoeC(Nsheath)=FS%CoeC(FS%Ns/2)
           FSQNS%CoeA(Nsheath)=1.d0
           FSQNS%CoeB(Nsheath)=-2.d0
           FSQNS%CoeC(Nsheath)=1.d0
           FSQNS%Source(Nsheath)=0.d0
           
           do i=1,Nsheath-1
              FSQNS%CoeA(i)=FS%CoeA(i)
              FSQNS%CoeB(i)=FS%CoeB(i)
              FSQNS%CoeC(i)=FS%CoeC(i)
              FSQNS%Source(i)=FS%Source(i)
              FSQNS%CoeA(FSQNS%Nxs-i+1)=FS%CoeA(FS%Ns-i+1)
              FSQNS%CoeB(FSQNS%Nxs-i+1)=FS%CoeB(FS%Ns-i+1)
              FSQNS%CoeC(FSQNS%Nxs-i+1)=FS%CoeC(FS%Ns-i+1)
              FSQNS%Source(FSQNS%Nxs-i+1)=FS%Source(FS%Ns-i+1)
            ENd Do

           Call Tridag(FSQNS%CoeA,FSQNS%CoeB,FSQNS%CoeC,FSQNS%Source,FSQNS%Solve,FSQNS%Nxs)
           
           FS%Solve(FSQNB%NL:FSQNB%NR)=FSQNS%Solve(Nsheath)
           
           
           do i=1,Nsheath-1
              FS%Solve(i)=FSQNS%Solve(i)
              FS%Solve(FS%Ns-i+1)=FSQNS%Solve(FSQNS%Nxs-i+1)
            ENd Do
           return
        End subroutine UpdaterFieldSolverQNS
    
        !  subroutine UpdaterFieldSolver(FS)
        !    Implicit none
        !   Class(FieldSolver), intent(inout):: FS
        !   Call Tridag(FS%CoeA,FS%CoeB,FS%CoeC,FS%Source,FS%Solve,Fs%Ns)
        !   return
        !End subroutine UpdaterFieldSolver

        
        !subroutine UpdaterField(FG,FS,FB)
        !    Implicit none
        !    Type(Field),intent(inout) :: FG
        !    Type(FieldSolver), intent(in):: FS
        !    Type(FieldBoundary), intent(in):: FB
        !    Integer(4) :: i
        !    Associate (Nx=>FG%Nx)
        !        FG%Phi(2:Nx-1)=FS%Solve(1:Nx-2)
        !        FG%Phi(1)=FB%V1
        !        FG%Phi(Nx)=FB%V2
        !        FG%Ex(1)=(3.d0*FG%Phi(1)-4.d0*FG%Phi(2)+FG%Phi(3))/(2.d0*FG%dx)
        !        do i=2,Nx-1
        !            FG%Ex(i)=(FG%Phi(i-1)-FG%Phi(i+1))/(2.d0*FG%dx)
        !        end do
        !        FG%Ex(Nx)=-1.d0*(FG%Phi(Nx-2)-4.d0*FG%Phi(Nx-1)+3.d0*FG%Phi(Nx))/(2.d0*FG%dx)
        !    End Associate
        !    return
        !End subroutine UpdaterField
        
        !subroutine UpdaterCoeFieldSolver(FS,FG,FB)
        !    Implicit none
        !    Type(FieldSolver), intent(inout):: FS
        !    Type(Field), intent(in) :: FG
        !    Type(FieldBoundary), intent(in) :: FB
        !    Integer(4) :: i
        !    Real(8) :: GeoFactor
        !     Real(8) :: ATemp, BTemp,CTemp
        !    Associate (Ns=>FS%Ns,Nx=>FG%Nx,dx=>FG%dx)
        !            do i=1,FS%Ns
        !                     ATemp=1.d0+0.5d0*(FG%Chi(i)+FG%Chi(i+1))
        !                      !ATemp=1.d0+Max(Chi(i),Chi(i+1))
        !                     FS% CoeA(i)=ATemp
        !                      !CTemp=1.d0+Max(Chi(i),Chi(i+1))
        !                      CTemp=1.d0+0.5d0*(FG%Chi(i+1)+FG%Chi(i+2))
        !                      FS%CoeC(i)=CTemp
        !          
        !                      BTemp=-(ATemp+CTemp)
        !                      FS%CoeB(i)=BTemp
        !                End do
        !
        !            GeoFactor=-dx*dx/Epsilon
        !            FS%Source(1:Nx-2)=FG%Rho(2:Nx-1)*GeoFactor
        !            FS%Source(1)=FS%Source(1)-FB%V1*FS%CoeA(1)
        !            FS%Source(Nx-2)=FS%Source(Nx-2)-FB%V2*FS%CoeC(Nx-2)
        !    End Associate
        !    return
        !End subroutine UpdaterCoeFieldSolver
        
End Module ModuleOneStepFieldQNEx