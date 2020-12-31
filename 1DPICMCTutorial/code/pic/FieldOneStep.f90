Module ModuleOneStepField
  Use ModuleField
  Use ModuleFieldBoundary
  Use Numrical
     
     Type(Field) :: FieldGlobal
     Type(FieldSolver) :: FieldSolverGlobal
     Type(FieldBoundary) :: FieldBoundaryGlobal
     Type(FieldOne) :: FieldOneGlobalCO[2,*]
     Type(FieldOne),Allocatable :: FieldOneGlobal(:)
                  
    contains
         subroutine InitializationField(CF)
            Implicit none
            Type(ControlFlow),intent(inout) :: CF
            Logical ::  Status
            Call FieldBoundaryGlobal%Init(CF)
            Call FieldGlobal%Load(Status)
            Allocate(FieldOneGlobal(0:CF%Ns))
            !Allocate(FieldOneGlobalCO(0:CF%Ns)[*])
         End subroutine InitializationField
         
    
         subroutine FieldOneStep(Ns,FO,FG,FB,FS)
            Implicit none
            Integer(4),intent(in) :: Ns
            Type(FieldOne),intent(in) :: FO(0:Ns)
            Type(Field),intent(inout) :: FG
            Type(FieldSolver), intent(inout):: FS
            Type(FieldBoundary), intent(inout):: FB
            Call AccumulationField(FG,Ns,FO)
            Call FB%Updater
            Call UpdaterCoeFieldSolver(FS,FG,FB)
            Call UpdaterFieldSolver(FS)
            Call UpdaterField(FG,FS,FB)
            return
         End subroutine FieldOneStep
         
         subroutine FieldOneStepCO(Ns,FO,FG,FB,FS)
            Implicit none
            Integer(4),intent(in) :: Ns
            Type(FieldOne),intent(inout) :: FO(0:Ns)
            !Type(FieldOne),intent(in) :: FOCO(0:Ns)
            Type(Field),intent(inout) :: FG
            Type(FieldSolver), intent(inout):: FS
            Type(FieldBoundary), intent(inout):: FB
            Integer(4) :: i,j
            
            do i=0,Ns
            FO(i)%RhoOne=0.d0
            FO(i)%ChiOne=0.d0
            do j=1,num_images()
               FO(i)%RhoOne=FO(i)%RhoOne+FieldOneGlobalCO[i,j]%RhoOne
               FO(i)%ChiOne=FO(i)%ChiOne+FieldOneGlobalCO[i,j]%ChiOne
            ENd do
            ENd DO
            !Call AccumulationFieldCO(FG,Ns,FO,FOCO)
            Call AccumulationField(FG,Ns,FO)
            Call FB%Updater
            Call UpdaterCoeFieldSolver(FS,FG,FB)
            Call UpdaterFieldSolver(FS)
            Call UpdaterField(FG,FS,FB)
            return
         End subroutine FieldOneStepCO
         


        subroutine AccumulationField(FG,Ns,FO)
            Implicit none
            Type(Field),intent(inout) :: FG
            Integer(4),intent(in) :: Ns
            Type(FieldOne),intent(in) :: FO(0:Ns)
            Integer(4) :: i!,Nx
            FG%Rho=0.d0
            FG%Chi=0.d0
            !Nx=FG%Nx
                    do i=0,Ns
                            FG%Rho=FG%Rho+FO(i)%RhoOne
                            FG%Chi=FG%Chi+FO(i)%ChiOne
                    End do
            return
    End subroutine AccumulationField
    
          subroutine UpdaterFieldSolver(FS)
            Implicit none
           Class(FieldSolver), intent(inout):: FS
           Call Tridag(FS%CoeA,FS%CoeB,FS%CoeC,FS%Source,FS%Solve,Fs%Ns)
           return
        End subroutine UpdaterFieldSolver

        
        subroutine UpdaterField(FG,FS,FB)
            Implicit none
            Type(Field),intent(inout) :: FG
            Type(FieldSolver), intent(in):: FS
            Type(FieldBoundary), intent(in):: FB
            Integer(4) :: i
            Associate (Nx=>FG%Nx)
                FG%Phi(2:Nx-1)=FS%Solve(1:Nx-2)
                FG%Phi(1)=FB%V1
                FG%Phi(Nx)=FB%V2
                FG%Ex(1)=(3.d0*FG%Phi(1)-4.d0*FG%Phi(2)+FG%Phi(3))/(2.d0*FG%dx)
                do i=2,Nx-1
                    FG%Ex(i)=(FG%Phi(i-1)-FG%Phi(i+1))/(2.d0*FG%dx)
                end do
                FG%Ex(Nx)=-1.d0*(FG%Phi(Nx-2)-4.d0*FG%Phi(Nx-1)+3.d0*FG%Phi(Nx))/(2.d0*FG%dx)
            End Associate
            return
        End subroutine UpdaterField
        
        subroutine UpdaterCoeFieldSolver(FS,FG,FB)
            Implicit none
            Type(FieldSolver), intent(inout):: FS
            Type(Field), intent(in) :: FG
            Type(FieldBoundary), intent(in) :: FB
            Integer(4) :: i
            Real(8) :: GeoFactor
             Real(8) :: ATemp, BTemp,CTemp
            Associate (Ns=>FS%Ns,Nx=>FG%Nx,dx=>FG%dx)
                    do i=1,FS%Ns
                             ATemp=1.d0+0.5d0*(FG%Chi(i)+FG%Chi(i+1))
                              !ATemp=1.d0+Max(Chi(i),Chi(i+1))
                             FS% CoeA(i)=ATemp
                              !CTemp=1.d0+Max(Chi(i),Chi(i+1))
                              CTemp=1.d0+0.5d0*(FG%Chi(i+1)+FG%Chi(i+2))
                              FS%CoeC(i)=CTemp
                  
                              BTemp=-(ATemp+CTemp)
                              FS%CoeB(i)=BTemp
                        End do

                    GeoFactor=-dx*dx/Epsilon
                    FS%Source(1:Nx-2)=FG%Rho(2:Nx-1)*GeoFactor
                    FS%Source(1)=FS%Source(1)-FB%V1*FS%CoeA(1)
                    FS%Source(Nx-2)=FS%Source(Nx-2)-FB%V2*FS%CoeC(Nx-2)
            End Associate
            return
        End subroutine UpdaterCoeFieldSolver
        
End Module ModuleOneStepField