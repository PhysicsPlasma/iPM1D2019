Module ModuleMCCSigma
     Use Numrical
     Use ModuleTypeMCC
     Use ModuleSpecyOne
     Implicit none
     ! 1-linear 2- Log10
     Integer(4),parameter,private :: MCEnergyScale=1
     
     Real(8),parameter,private :: ElectronEnergyMin=-3.d0,ElectronEnergyMax=3.d0
     Real(8),parameter :: IonEnergyMin=-1.d0,IonEnergyMax=3.d0
     Integer(4),parameter :: NSigmaElectron=10001,NSigmaIon=1001
     !  The points of the Sigma data is simply devided into two kinds, 0-larger for Electron, 1-smaller for ion.
     Integer(4),parameter,private :: NSigmaRawMax=500_4,NDataMax=9
     
     !Real(8)  :: Threshold,MaxEnergy
     !Integer(4) :: NSigmaRaw
     
    Type  SigmaRaw
        Type(ReactionOne) :: Reaction
        Real(8) :: EnergySigma(1:2,1:NSigmaRawMax)=0.d0
    end Type  SigmaRaw
                              
    Type SigmaNormalized
        !Character(99) :: Name="He"
        Integer(4) :: Model,NReaction=0,NSigma=0
        Real(8) :: EnergyMin,EnergyInterval=0.d0,EnergyMax
        Type(ReactionOne),Allocatable :: Reaction(:)
        Real(8),Allocatable ::  Sigma(:,:)
     End Type SigmaNormalized
     
    contains

    
      subroutine  SigmaNormalization(SN,SO,GO,Inited)
                  Implicit none
                  Type(SigmaNormalized),Intent(out) :: SN !, intent(inout)
                  Type(SpecyOne),Intent(in)  :: SO
                   Type(GasOne),Intent(in) :: GO
                   Logical,Intent(out) :: Inited
                  Type(SigmaRaw),Allocatable :: SR(:)
                  Integer(4) :: Model,NReaction
                  Character(len=99) :: Filename
                  !Logical :: Alive
                  NAMELIST /MRList/ Model,NReaction
                  NAMELIST /SRList/ SR
                  
                  !GO%Name="Ar"
                  Filename="./gas/"//Trim(GO%Name)//"/"//Trim(SO%Name)//Trim(GO%Name)//".txt"

                  Inquire(file=Filename,exist=Inited)
                  If(Inited) then
                         OPEN(17,FILE=Filename)
                          Read(17,NML=MRList)
                          Call SigmaNormalizedInit(SN,Model,NReaction)
                          
                          Allocate(SR(NReaction))
                          !OPEN(18,FILE='test6.txt')
                          !WRITE(18,NML=SRList)
                          Read(17,NML=SRList)
                          close (17)
                            Select case(Model)
                              Case (1)
                                 Call SigmaNormalizedUpdate(SN,SR)
                              Case(2)
                                Call SigmaNormalizedUpdate(SN,SR)
                              Case(3)
                                Call SigmaNormalizedUpdateReactive(SN,SR,SO,GO)
                            End Select
                            Call UpdateReactionIndex(SN,GO)
                            Call SigmaNormalizedFinalization(SN)
                            
                          DeAllocate(SR)
                  End If
                  return
      End  subroutine
      
      
   Subroutine  SigmaNormalizedInit(SN,Model,NReaction)
                  Implicit none
                  Type(SigmaNormalized), intent(inout) :: SN
                  Integer(4), intent(in) :: Model,NReaction
                  Real(8) ::  EnergyInterval,EnergyMin,EnergyMax
                  Integer(4) :: NSigma
                  SN%Model=Model
                  SN%NReaction=NReaction
                            Select case (SN%Model)
                                  Case(1_4)  ! For electrons
                                      EnergyMin=ElectronEnergyMin
                                      EnergyMax=ElectronEnergyMax
                                      EnergyInterval=(ElectronEnergyMax-ElectronEnergyMin)/dble(NSigmaElectron-1)
                                      SN%EnergyMin=EnergyMin
                                      SN%EnergyInterval=EnergyInterval
                                      SN%EnergyMax=EnergyMax 
                                      NSigma=NSigmaElectron
                                      SN%NSigma=NSigma
                                  Case(2_4)   ! For nonreactive ions
                                      EnergyMin=IonEnergyMin
                                      EnergyMax=IonEnergyMax 
                                      EnergyInterval=(IonEnergyMax-IonEnergyMin)/dble(NSigmaElectron-1)
                                      SN%EnergyMin=EnergyMin
                                      SN%EnergyInterval=EnergyInterval
                                      SN%EnergyMax=EnergyMax 
                                      NSigma=NSigmaElectron
                                      SN%NSigma=NSigma
                                  Case(3_4)   !  For reactive ions
                                      EnergyMin=IonEnergyMin
                                      EnergyMax=IonEnergyMax 
                                      EnergyInterval=(IonEnergyMax-IonEnergyMin)/dble(NSigmaIon-1)
                                      SN%EnergyMin=EnergyMin
                                      SN%EnergyInterval=EnergyInterval
                                      SN%EnergyMax=EnergyMax 
                                      NSigma=NSigmaIon
                                      SN%NSigma=NSigma
                                  End Select
                              If(allocated(SN%Reaction)) DeAllocate(SN%Reaction)
                              Allocate(SN%Reaction(NReaction))
                              If(allocated(SN%Sigma)) DeAllocate(SN%Sigma)
                              Allocate(SN%Sigma(NReaction,NSigma))
                   return
    End  subroutine  SigmaNormalizedInit
         
    subroutine  SigmaNormalizedUpdate(SN,SR)
                  Implicit none
                  Type(SigmaNormalized), intent(inout) :: SN
                  Type(SigmaRaw), intent(inout) :: SR(SN%NReaction)
                  Type(SigmaRaw) :: TempRaw
                  Integer(4) :: i,j,k
                  Integer(4) :: NSigmaRaw,NSigmaRawLoc(1)
                  Real(8) :: Threshold,MaxEnergy,Alfa
                  Real(8) :: NEnergy,Energy,LEnergy,UEnergy,S1,S2
                  
                  Associate(Emin=>SN%EnergyMin,Eint=>SN%EnergyInterval,Nr=>SN%NReaction,Ns=>SN%NSigma)
                   Do j=1,Nr
                                   TempRaw=SR(j)
                                   SN%Reaction(j)=TempRaw%Reaction
                                   
                                   Threshold=TempRaw%Reaction%Threshold
                                   NSigmaRawLoc= Maxloc(TempRaw%EnergySigma(1,:))
                                   NSigmaRaw=NSigmaRawLoc(1)
                                   !SS=Maxloc(TempRaw%EnergySigma(1,:))
                                   !TempRaw%EnergySigma(1,:)=0.d0
                                   MaxEnergy=TempRaw%EnergySigma(1,NSigmaRaw)
                                   Alfa= TempRaw%EnergySigma(2,NSigmaRaw)*MaxEnergy
                                   
                                   do i=1,NSigmaRaw
                                        If (TempRaw%EnergySigma(1,i)<10.d0**Emin) then
                                             TempRaw%EnergySigma(1,i)=Emin
                                        Else    
                                             TempRaw%EnergySigma(1,i)=(Log10(TempRaw%EnergySigma(1,i))-Emin)/Eint
                                         End if 
                                    End do
                                    do i=1,Ns
                                         NEnergy=dble(i-1)
                                         Energy=10.d0**(Emin+dble(i-1)*Eint)
                                          If (Energy<=Threshold) then
                                              SN%Sigma(j,i)=0.d0
                                          Else if (Energy>=MaxEnergy)  then
                                              SN%Sigma(j,i)=Alfa/Energy
                                          Else
                                             Call Locate(TempRaw%EnergySigma(1,:),NSigmaRaw,NEnergy,k)
                                             If (k<1) Then
                                                !Write(*,*) NEnergy,Energy,j,i,k
                                                k=1
                                             End If
                                              LEnergy=TempRaw%EnergySigma(1,k)
                                              UEnergy=TempRaw%EnergySigma(1,k+1)
                                              S1=(dble(i)-LEnergy)/(UEnergy-LEnergy)
                                              S2=1.d0-S1

                                              SN%Sigma(j,i)=S2*TempRaw%EnergySigma(2,k)+S1*TempRaw%EnergySigma(2,k+1)
                                          end if
                                    end do
                                    !SR(j)%EnergySigma(2,:)=SR(j)%EnergySigma(2,:)*1.d-20
                                    Call SigmaRawFinalization(SR(j))
                                    
                   End do
                   !SR(1)%EnergySigma(2,:)=SR(1)%EnergySigma(2,:)*1.d-20
                   !Call SigmaRawFinalization(SR(1))
                  END Associate
                Return
    End  subroutine  SigmaNormalizedUpdate
    
    Subroutine  SigmaNormalizedUpdateReactive(SN,SR,SO,GO)
                  Implicit none
                  Type(SigmaNormalized), intent(inout) :: SN
                  Type(SigmaRaw), intent(in) :: SR(SN%NReaction)
                  Type(SpecyOne),Intent(in)  :: SO
                  Type(GasOne),Intent(in) :: GO
                  Integer(4) :: NComplex,i,j
                  Real(8) :: S,NEnergy,Energy,PSum,PTemp(SN%NReaction)
                     
                    NComplex=SO%NAtom+GO%NAtom
                    S=(3.d0*dble(NComplex)-6.d0)/2.d0-1.d0
                    SN%Reaction=SR%Reaction
                     Associate(Emin=>SN%EnergyMin,Eint=>SN%EnergyInterval,Nr=>SN%NReaction,Ns=>SN%NSigma)
                        Do i=1,Ns
                              NEnergy=dble(i-1)
                              Energy=10.d0**(Emin+dble(i-1)*Eint)
                              Do j=1, Nr
                                    If(Energy>SN%Reaction(j)%Threshold)  then
                                      PTemp(j)=(Energy-SN%Reaction(j)%Threshold)**S
                                  else
                                      PTemp(j)=0.d0
                                  end if
                              End Do
                                PSum=SUM(PTemp)
                                SN%Sigma(:,i)=PTemp/PSum
                          ENd Do
                     END Associate       
        return
    end subroutine  SigmaNormalizedUpdateReactive
    
    Subroutine SigmaRawFinalization(SR)
               Implicit none
              ! Integer(4) 
              Type(SigmaRaw), intent(inout) :: SR
              !Integer(4),intent(in) :: Mode 
               Character(len=99) :: Filename
               Integer(4) :: NSigmaRaw,NSigmaRawLoc(1)
               Integer(4),save :: i,Timer=0
               Timer=Timer+1
               NSigmaRawLoc= Maxloc(SR%EnergySigma(1,:))
               NSigmaRaw=NSigmaRawLoc(1)
                   Write(filename,*) "SR",Timer,".dat"
                   Open (10,file=filename)
                        do i=1,NSigmaRaw
                            Write(10,FMt="(*(es21.14,1x))") SR%EnergySigma(:,i)
                        ENd do
                   Close(10)
               return
     End Subroutine SigmaRawFinalization
             
    Subroutine SigmaNormalizedFinalization(SN)
               Implicit none
              Type(SigmaNormalized), intent(inout) :: SN
              !Integer(4),intent(in) :: Mode 
               Character(len=99) :: Filename
               Real(8) :: NEnergy,Energy
               Integer(4),save :: i,j,Timer=0
               Timer=Timer+1
               Associate(Emin=>SN%EnergyMin,Eint=>SN%EnergyInterval,Nr=>SN%NReaction,Ns=>SN%NSigma)
                   Write(filename,*) "SN",Timer,".dat"
                   Open (10,file=filename)
                        do i=1,Ns
                            Energy=10.d0**(Emin+dble(i-1)*Eint)
                            Write(10,FMt="(*(es21.14,1x))") Energy, (SN%Sigma(j,i),j=1,Nr)
                        ENd do
               ENd Associate
               Close(10)
               return
     End Subroutine SigmaNormalizedFinalization
      
    subroutine SigmaNormalizationFunc(SN,ExFunc)
                   implicit none
                   Type(SigmaNormalized), intent(out) :: SN
                   Real(8),external :: ExFunc
                   Integer(4),save :: NSigma,Index,NReaction,i
                   Real(8) :: Energy
                       Index=Index+1 
                       NReaction=SN%NReaction
                       NSigma=SN%NReaction
                       do i=0,NSigma-1
                                          Energy=10.d0**(SN%EnergyMin+dble(i)*SN%EnergyInterval)
                                          SN%Sigma(Index,i)=ExFunc(Energy)
                       end do
                       return
    End  subroutine  SigmaNormalizationFunc
    
    Subroutine  UpdateReactionIndex(SN,GOT,GOI)
             Implicit none
             Type(SigmaNormalized),intent(inout) :: SN
             Type(GasOne),intent(in) :: GOT
             Type(GasOne),intent(in),Optional :: GOI
             Integer(4) :: i
             do i=1,SN%NReaction
                 If (Present(GOI)) Then
                     SN%Reaction(i)%Reactant=SN%Reaction(i)%Reactant+GOI%IndexStart
                 End If
                 SN%Reaction(i)%Resultant=SN%Reaction(i)%Resultant+GOT%IndexStart
             End do
             return
      End subroutine UpdateReactionIndex   

End module ModuleMCCSigma
