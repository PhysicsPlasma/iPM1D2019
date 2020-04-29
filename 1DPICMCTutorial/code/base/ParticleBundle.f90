Module ModuleParticleBundle
    Use ModuleParticleOne
    Use ModuleSpecyOne
    use h5fortran, only: hdf5_file
   Implicit none
   type(hdf5_file) :: h5f
   Type :: ParticleBundle
         Integer(4) :: NPar=0,NParNormal=1000
         Real(8) :: XFactor,VFactor
         Logical :: LXScaled=.False.,LVScaled=.False.
         Real(8) :: Charge,Mass,Weight
         Real(8) :: Dx,Dt
         Type(SpecyOne),Pointer :: SO
         Type(ParticleOne),Allocatable :: PO(:)
         Contains
         Procedure :: AllInit=>AllInitializationParticleBundle
         Procedure :: AddOne=>AddParticleOneParticleBundle
         Procedure :: DelOne=>DelParticleOneParticleBundle
         !Procedure :: NumRes2=>ParticleBundleNumberRescale2
         Procedure :: PosRes=>PositionRescaleParticleBundle
         Procedure :: VelRes=>VelocityRescaleParticleBundle
         Procedure :: MoveES=>MoveElectrostaticParticleBundle
         Procedure :: MoveEM=>MoveElectromagneticParticleBundle
         Procedure :: WeightP2C=>WeightP2CParticleBundle
         Procedure :: Dump=>DumpParticleBundle
         Procedure :: Load=>LoadParticleBundle
         Procedure :: Norm=>ParticleBundleNormalization
         !Procedure :: RescaleFieldOne=>RescaleFieldOneParticleBundle
         End Type ParticleBundle
         
         Type :: ParticleBundleIndex
         Integer(4) :: NPar=0,NParNormal=10000
         Type(ParticleOneIndex),Allocatable :: POI(:)
         Contains
         Procedure :: Init=>InitializationParticleBundleIndex
         Procedure :: AddOne=>AddParticleOneParticleBundleIndex
         Procedure :: DelOne=>DelParticleOneParticleBundleIndex
         !Procedure :: NumRes2=>ParticleBundleNumberRescale2
         !Procedure :: Dump=>DumpParticleBundleIndex
         End Type ParticleBundleIndex

  contains
!DIR$ ATTRIBUTES FORCEINLINE ::  MoveElectrostaticParticleBundle,  MoveElectromagneticParticleBundle, WeightP2CParticleBundle      
           Subroutine AllInitializationParticleBundle(PB,SO,CF)
             Class(ParticleBundle), intent(inout) :: PB
             Type(SpecyOne), intent(in),Target :: SO 
             Class(ControlFlow), intent(in) :: CF
             Integer(4) :: i,NPArMax
             Real(8) :: XFactor,VFactor
             Logical ::  Status
              PB%dx=CF%dx
              PB%dt=CF%dt
              PB%XFactor=PB%dx
              PB%VFactor=PB%dx/PB%dt
              
              PB%SO=>SO
              PB%Charge=PB%SO%Charge
              PB%Mass=PB%SO%Mass
              PB%Weight=CF%InitDensity/Dble(CF%ParticlePerGrid)

              PB%NParNormal=CF%ParticlePerGrid*(CF%Nx-1)*PB%SO%InitDensity
              PB%NPar=PB%NParNormal
              
              NPArMax=Ceiling(2.5*PB%NParNormal)
              If(Allocated(PB%PO)) Deallocate(PB%PO)
              Allocate(PB%PO(NPArMax))
                      If (CF%ReStartParticles) Then
                              Do i=1,PB%NPAR
                              PB%LXScaled=.True.
                              Call PB%PO(i)%PosInit(Dble(CF%NxL),Dble(CF%NxU))
                              PB%LVScaled=.True.
                              Call PB%PO(i)%VelMaxInit(PB%SO%Mass,PB%SO%InitTemperature)
                              VFactor=1.d0/PB%VFactor
                              Call PB%PO(i)%VelRes(VFactor)
                              Call PB%PO(i)%AccInpInit()
                              End Do
                      Else
                          Call LoadParticleBundle(PB,Status)
                          If (Status) Then
                              Do i=1,PB%NPAR
                                  PB%LXScaled=.True.
                                  Call PB%PO(i)%PosInit(Dble(CF%NxL),Dble(CF%NxU))
                                  PB%LVScaled=.True.
                                  Call PB%PO(i)%VelMaxInit(PB%SO%Mass,PB%SO%InitTemperature)
                                  VFactor=1.d0/PB%VFactor
                                  Call PB%PO(i)%VelRes(VFactor)
                                  Call PB%PO(i)%AccInpInit()
                              End Do
                          Else
                              XFactor=1.d0/PB%XFactor
                              Call PB%PosRes()
                              VFactor=1.d0/PB%VFactor
                              Call PB%VelRes()
                          End If
                      End If
                      
          End Subroutine AllInitializationParticleBundle

          Subroutine  AddParticleOneParticleBundle(PB,PO)
           Implicit none
               Class(ParticleBundle),intent(inout) ::  PB
               Class(ParticleOne),intent(in) ::  PO
                PB%NPar=PB%NPar+1
                PB%PO(PB%NPar)=PO
           return   
           End  Subroutine  AddParticleOneParticleBundle
   
           Subroutine  DelParticleOneParticleBundle(PB,NDel)
           Implicit none
               Class(ParticleBundle),intent(inout) ::  PB
               Integer(4),intent(in) ::  NDel
               PB%PO(NDel)=PB%PO(PB%NPar)
               PB%NPar=PB%NPar-1
           return   
           End  Subroutine  DelParticleOneParticleBundle
           
           subroutine PositionRescaleParticleBundle(PB)
             Class(ParticleBundle), intent(inout) :: PB
             Integer(4) :: i
             Real(8) :: XFactor
             If (PB%LXScaled) Then
                 XFactor=PB%XFactor
                 PB%LXScaled=.False.
             Else
                 XFactor=1.d0/PB%XFactor
                PB%LXScaled=.True.
              End If
              Do i=1,PB%NPar
                    Call PB%PO(i)%PosRes(XFactor)
              End Do
           end subroutine PositionRescaleParticleBundle
           
            subroutine VelocityRescaleParticleBundle(PB)
             Class(ParticleBundle), intent(inout) :: PB
             Integer(4) :: i
             Real(8) :: VFactor
             If (PB%LVScaled) Then
                 VFactor=PB%VFactor
                 PB%LVScaled=.False.
             Else
                 VFactor=1.d0/PB%VFactor
                 PB%LVScaled=.True.
              End If
              Do i=1,PB%NPar
                    Call PB%PO(i)%VelRes(VFactor)
              ENd DO
            end subroutine VelocityRescaleParticleBundle
            
          subroutine MoveElectrostaticParticleBundle(PB,FG)
             Class(ParticleBundle), intent(inout) :: PB
             Class(Field), intent(inout) :: FG
             Integer(4) :: i
             Real(8) :: Ex,Ey,Ez,EFactor
             EFactor=PB%Charge/PB%Mass*PB%dt/(PB%dx/PB%dt)
              Do i=1,PB%NPar
                    Call PB%PO(i)%WeightC2PES(FG,Ex)
                    Ex=Ex*EFactor
                    Call PB%PO(i)%MoveES(Ex)
              End Do
           End subroutine MoveElectrostaticParticleBundle  

          Subroutine MoveElectromagneticParticleBundle(PB,FG)
             Class(ParticleBundle), intent(inout) :: PB
             Type(Field), intent(in) :: FG
             Integer(4) :: i
             Real(8) :: Ex,Ey,Ez,Bx,By,Bz,EFactor,BFactor
             EFactor=PB%Charge/PB%Mass*PB%dt/(PB%dx/PB%dt)
             BFactor=PB%Charge/PB%Mass*PB%dt
              Do i=1,PB%NPar
                    Call PB%PO(i)%WeightC2PEM(FG,Ex,Ey,Ez,Bx,By,Bz)
                    Ex=Ex*EFactor
                    Ey=Ey*EFactor
                    Ez=Ez*EFactor
                    Bx=Bx*BFactor
                    By=By*BFactor
                    Bz=Bz*BFactor
                    Call PB%PO(i)%MoveEM(Ex,Ey,Ez,Bx,By,Bz)
              End Do
          End subroutine MoveElectromagneticParticleBundle
          
          subroutine WeightP2CParticleBundle(PB,FO)
             Class(ParticleBundle), intent(inout) :: PB
             Type(FieldOne), intent(inout) :: FO
             Integer(4) :: i
             Real(8) :: RhoFactor,ChiFactor
             FO%RhoOne=0.d0
             FO%ChiOne=0.d0
              Do i=1,PB%NPar
                    Call PB%PO(i)%WeightP2C(FO)
              End Do
              RhoFactor=PB%Charge*PB%Weight
              FO%RhoOne=FO%RhoOne*RhoFactor
              ChiFactor=0.5d0*PB%Charge/PB%Mass*PB%dt*PB%dt/Epsilon
              FO%ChiOne=FO%RhoOne*ChiFactor
              Return
          End subroutine WeightP2CParticleBundle
          
          subroutine DumpParticleBundle(PB,Mode)
              implicit none
              Class(ParticleBundle),intent(inout) :: PB
              Type(ParticleBundle) :: TempPB
              Integer(4),intent(in)  :: Mode
              Character(len=99) :: Filename
              Integer(4) :: i,NameIndex
              Real(8) :: XFactor,YFactor
              Integer(4) :: a1(2)
              real(8) :: a2(3)
              real(8),allocatable :: a3(:,:)
              Integer(4),save :: NameTimer=1
              If (Mode==0) Then
                         NameIndex=DefaultNameIndex
              else
                         NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
                         NameTimer=NameTimer+1
              End If 
              TempPB=PB
              Call TempPB%PosRes
              Call TempPB%VelRes
              Write(Filename,*) NameIndex,trim(TempPB%SO%Name),".h5"
              !Write(Filename,*)trim(adjustl(TempPB%SO%Name)),".h5"
              Filename=adjustl(Filename)
              Write(*,*) 'Saving ', trim(TempPB%SO%Name), TempPB%NPar,' Please Wait...' 
              allocate(a3(TempPB%NPar,7))
              a1(1)=TempPB%NPar
              a1(2)=TempPB%NParNormal
              a2(1)=TempPB%Charge
              a2(2)=TempPB%Mass
              a2(3)=TempPB%Weight
              do i=1,a1(1)
                a3(i,1)=TempPB%PO(i)%X
                a3(i,2)=TempPB%PO(i)%Vx
                a3(i,3)=TempPB%PO(i)%Vy
                a3(i,4)=TempPB%PO(i)%Vz
                a3(i,5)=TempPB%PO(i)%Ax
                a3(i,6)=TempPB%PO(i)%Ay
                a3(i,7)=TempPB%PO(i)%Az
              End do
              call h5f%initialize(Filename, status='new',action='w')
              call h5f%write('/x',a1)
              call h5f%write('/y',a2)
              call h5f%write('/z',a3)
              call h5f%finalize()
              Write(*,*) 'Saving ', trim(TempPB%SO%Name),' Complete!'  
              return
          end subroutine  DumpParticleBundle
 
          subroutine LoadParticleBundle(PB,Status)
              implicit none
              Class(ParticleBundle),intent(inout) :: PB
              Logical,intent(inout) ::  Status
              Logical :: alive
              Character(len=99) :: Filename
              Integer(4) :: i,NameIndex,NPArMax
              Integer(4) :: b1(2)
              real(8) :: b2(3)
              real(8),allocatable :: b3(:,:)
              NameIndex=DefaultNameIndex
              Write(Filename,*) NameIndex,trim(PB%SO%Name),".h5"
              !Write(Filename,*)trim(adjustl(PB%SO%Name)),".h5"
              Filename=adjustl(Filename)
              inquire(file=Filename,exist=alive)
              If(alive) then
                     call h5f%initialize(Filename, status='old',action='r')
                     call h5f%read('/x',b1)
                     call h5f%read('/y',b2)
                     PB%NPar=b1(1)
                     PB%NParNormal=b1(2)
                     PB%Charge=b2(1)
                     PB%Mass=b2(2)
                     PB%Weight=b2(3)
                     allocate(b3(PB%NPar,7))
                     call h5f%read('/z',b3)
                     call h5f%finalize()
                     Write(*,*) 'Loading ', trim(PB%SO%Name), PB%NPar,' Please Wait...' 
                      NPArMax=Ceiling(2.5*PB%NParNormal)
                      !If(Allocated(PB%PO)) Deallocate(PB%PO)
                     ! Allocate(PB%PO(NPArMax))
                      do i=1,PB%NPar
                        PB%PO(i)%X=b3(i,1)
                        PB%PO(i)%Vx=b3(i,2)
                        PB%PO(i)%Vy=b3(i,3)
                        PB%PO(i)%Vz=b3(i,4)
                        PB%PO(i)%Ax=b3(i,5)
                        PB%PO(i)%Ay=b3(i,6)
                        PB%PO(i)%Az=b3(i,7)
                     End do
                     Status=.False.
                     Write(*,*) 'Loading ', trim(PB%SO%Name),' Complete!'   
               else
                     Write(*,*) 'Can not find the file for ', trim(PB%SO%Name),' the particle will be randomly initilalized.'    
                     Status=.True. 
               End if        
              return
          end subroutine  LoadParticleBundle
          
          Subroutine ParticleBundleNormalization(PB,NParNew) 
            Implicit none
            Class(ParticleBundle),intent(inout) :: PB
            Integer(4),intent(in) :: NParNew
            Type(ParticleOne) :: ParticleTemp
            Integer(4) :: i,NTemp,Ndiff,Index
            NTemp=PB%NPar
            PB%Weight=PB%Weight*dble(NTemp)/dble(NParNew)
            If(NParNew>NTemp) then
               Ndiff=NParNew-NTemp
               do  i=1,Ndiff
                     Index=Ceiling(R*NTemp)
                     ParticleTemp=PB%PO(i)
                     Call PB%AddOne(ParticleTemp)
               end do
           else
               Ndiff=NTemp-NParNew
               do i=1,NDiff
                     CALL RANDOM_NUMBER(R)
                     Index=Ceiling(R*PB%NPar)
                     Call PB%DelOne(Index)
               end do
           end if  
           return 
         end  Subroutine  ParticleBundleNormalization
          
           Subroutine InitializationParticleBundleIndex(PBI)
             Class(ParticleBundleIndex), intent(inout) :: PBI
              If(Allocated(PBI%POI)) Deallocate(PBI%POI)
              Allocate(PBI%POI(PBI%NParNormal))
              PBI%NPar=0
              return
          End Subroutine InitializationParticleBundleIndex

          Subroutine  AddParticleOneParticleBundleIndex(PBI,POI)
           Implicit none
               Class(ParticleBundleIndex),intent(inout) ::  PBI
               Type(ParticleOneIndex),intent(in) ::  POI
                PBI%NPar=PBI%NPar+1
                Call PBI%POI(PBI%NPar)%Copy(POI)
           return   
           End  Subroutine  AddParticleOneParticleBundleIndex
   
           Subroutine  DelParticleOneParticleBundleIndex(PBI,NDel)
           Implicit none
               Class(ParticleBundleIndex),intent(inout) ::  PBI
               Integer(4),intent(in) ::  NDel
               Call PBI%POI(NDel)%Copy(PBI%POI(PBI%NPar))
               PBI%NPar=PBI%NPar-1
           return   
           End  Subroutine  DelParticleOneParticleBundleIndex
          
          
       !   subroutine WeightParticleBundle(PB,FO)
             !Class(ParticleBundle), intent(inout) :: PB
       !      Class(FieldOne), intent(inout) :: FO
       !      Integer(4) :: i
       !       Do i=1,PB%NPar
       !             Call PB%PO(i)%WeightP2C(FO)
       !       End Do
       !   End subroutine WeightParticleBundle            
!RescaleFieldOneParticleBundle

  !  subroutine ParticleBundleMove2(PB)
             !Class(ParticleBundle), intent(inout) :: PB
       !      Integer(4) :: i
       !      Real(8) :: Ex=3.d0
       !       Associate(PO=>PB%PO)
       !           Select Type (PO)
       !               Type is (ParticleOne1DDI)
       !       Do i=1,PB%NPar
       !       PO(i)%Vx=PO(i)%Vx+0.5d0*Ex
       !       PO(i)%X=PO(i)%X+0.5d0*Ex
       !       PO(i)%Ax=0.5d0*(PO(i)%Ax+Ex)
       !       PO(i)%Vx=PO(i)%Vx+0.5d0*PO(i)%Ax
       !       PO(i)%X=PO(i)%X+PO(i)%Vx
       !       End Do
       !           ENd Select
       !           
       !       End  Associate
       !
       !   end subroutine ParticleBundleMove2
            
       !    Subroutine ParticleBundleNumberRescale2(PB,)
             !Class(ParticleBundle), intent(inout) :: PB 
       !      Integer(4),intent(in) :: NParNew
       !      Integer(4) :: i
       !      Real(8) :: VFactor
       !      If (PB%LVScaled) Then
       !          VFactor=PB%VFactor
       !          PB%LVScaled=.False.
       !      Else
       !          VFactor=1.d0/PB%VFactor
       !          PB%LVScaled=.True.
       !       End If
       !       Do i=1,PB%NPar
       !             Call PB%PO(i)%POVR(VFactor)
       !       ENd DO
       !    end subroutine ParticleBundleNumberRescale2
           
           
            
            
  
END Module ModuleParticleBundle