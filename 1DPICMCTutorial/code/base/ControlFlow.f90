Module ModuleControlFlow
     Use Constants
     Implicit none
              Integer(4),Parameter :: NxMax=257
              Real(8),parameter,private :: ZLength=0.02d0
              Real(8),parameter :: Inputdx=ZLength/dble(NxMax-1)
              Real(8),parameter:: Inputdt=0.25d-10  !4.d-12
              
            Integer(4),Parameter  ::  DefaultNameIndex=10000
            Integer(4),Parameter  ::  DefaultNameIndexInit=20000
            Integer(4),Parameter  ::  ModeMultiplier=1000
 
     !  Ns--NSpecy,Ng---NGas
     Type ControlFlow
           Real(8)  :: Dx=Inputdx,Dt=Inputdt
           Integer(4) :: ParticlePerGrid=100
           Real(8)  :: InitDensity=1.d16

           Integer(4) :: Ns=0,Ng=0
           Integer(4) :: Nx=NxMax,NxL=0,NxU=NxMax-1
           Integer(4) :: Timer=0,Period=0
           Integer(4) :: NRun=0,NDiagShort=1,NDiagLong=0
           Logical :: ReStartParticles=.TRUE.
           !contains
              !procedure :: Init=>InitializationControlFlow
           End Type ControlFlow
       contains     
        Subroutine InitializationControlFlow(CF)
               Type(ControlFlow) :: CF
                Logical :: Alive
                Character(len=99) :: Filename
                NAMELIST /ControlFlow/ CF
                Filename="./input/controlflow.txt"
                  Inquire(file=Filename,exist=alive)
                   If(alive) then
                       OPEN(20,FILE=Filename)
                       Read(20,NML=ControlFlow)
                       Close (20)
                   Else
                      Write(*,*)  "ControlFlow Load", Trim(Filename),"ERROR! The program will abort!"
                      !Stop
                   ENd If
        End subroutine InitializationControlFlow
End Module ModuleControlFlow