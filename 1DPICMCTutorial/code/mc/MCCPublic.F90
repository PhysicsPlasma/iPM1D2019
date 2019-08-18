Module ModuleMCCPublic
     Use ModuleTypeMCC
     Use ModuleParticleOne
     Use ModuleMCCElectron
     Use  ModuleMCCIon
     Use ModuleMCCSigma
     Use ModuleParticleBundle
     Implicit none

     Type(MCCBundle),Allocatable :: MCCBundleGlobal(:,:)
     Type(ParticleBundleIndex) :: PBIGlobal
     
     Integer(4),private :: NCollision
     Integer(4),Allocatable,private :: MCCParticleIndex(:)
     !Logical,Private :: ParticleAnnihilation=.False.
     Logical,Allocatable,private :: ParticleAnnihilationIndex(:)
     
     Type(MCCParticleOne),Allocatable,private :: MCCParticleBundle(:)
     

                !   Type MCCBundle
                !    Integer(4) ::  Model,NReaction=0,NSigma=0
                !    Real(8) :: EnergyMin,EnergyInterval,EnergyMax
                !    Real(8) :: CollisionRatio,SigmaMax
                !    Type(SpecyOne),Pointer :: SO
                !    Type(GasOne),Pointer :: GO
                !    Type(ReactionOne),Allocatable :: Reaction(:)
                !    Real(8),Allocatable ::  Probility(:,:)
                !    contains
                !    procedure :: Dump=>DumpMCCBundle
                !End Type MCCBundle


    contains
      Subroutine MCC(Ns,Ng,PB,SO,GO,MCCB) 
        Implicit none
        ! j is the Particle Index;
        Integer(4),intent(in) :: Ns,Ng
        Type(ParticleBundle),intent(inout) :: PB(0:Ns)
        Type(SpecyOne),intent(in) :: SO(0:Ns)
        Type(GasOne),intent(in) :: GO(Ng)
        Type(MCCBundle) :: MCCB(0:Ns,Ng)
        Integer(4) :: i,j,k,m,Index
        Integer(4) :: CounterParticleAnnihilationIndex,NParTemp
        Call PBIGlobal%Init

        Do i=1,Ng
              do j=0,Ns
                  If(MCCB(j,i)%Model/=0_4) then
                      CounterParticleAnnihilationIndex=0
                      NParTemp=0
                      Call SelectParticle(PB(j),MCCB(j,i)%CollisionRatio)
                      do k=1,NCollision
                           Call MCCParticleBundle(k)%POI%VelRes(PB(j)%VFactor)
                           Call MCCParticleBundle(k)%Updater(SO(j),GO(i))
                           !Call MCCParticleBundle(k)%VelocityUpdater(PB(j)%VFactor)
                           Select case (MCCB(j,i)%Model)
                                 Case(1_4)
                                   Call SelectProbility(MCCParticleBundle(k),MCCB(j,i))
                                     Index=MCCParticleBundle(k)%ReactionIndex
                                     If (Index>0) Then
                                         Call  SelectCollisionElectron(MCCParticleBundle(k),SO(j),GO(i),MCCB(j,i)%Reaction(Index))
                                     End If
                                 Case(2_4)
                                      Call SelectProbility(MCCParticleBundle(k),MCCB(j,i))
                                      Index=MCCParticleBundle(k)%ReactionIndex
                                      If (Index>0) Then
                                      Call  SelectCollisionIon(MCCParticleBundle(k),SO(j),GO(i),MCCB(j,i)%Reaction(Index))
                                      End If
                                 Case(3_4)
                                      Call  SelectProbilityReactive(MCCParticleBundle(k),MCCB(j,i))
                                      Index=MCCParticleBundle(k)%ReactionIndex
                                      If (Index>0) Then
                                      Call  SelectCollisionIon(MCCParticleBundle(k),SO(j),GO(i),MCCB(j,i)%Reaction(Index))
                                      End If    
                                 End select
                                 Call MCCParticleBundle(k)%POI%VelRes(1.d0/PB(j)%VFactor)
                                 MCCParticleBundle(k)%POI=>Null()
                                 If (MCCParticleBundle(k)%ParticleCreation) Then
                                     Do m=1,MCCParticleBundle(k)%NPONew
                                           Call PBIGlobal%AddOne(MCCParticleBundle(k)%PON(m))
                                     End Do
                                 ENd If
                      End do

                      
                      !If (ParticleAnnihilation) Then
                              Do k=1,NCollision
                                    If (MCCParticleBundle(k)%ParticleAnnihilation) Then
                                         ParticleAnnihilationIndex(MCCParticleIndex(k))=.True.    
                                    ENd If
                              End Do
                      
                              Do k=PB(j)%NPar,1,-1
                                    If (ParticleAnnihilationIndex(k)) Then
                                        Call PB(j)%DelOne(k)
                                    End If    
                              End Do
                      !End If
                      Deallocate(MCCParticleBundle)
                  Endif
              ENddo
           ENd Do
           Call MCCAddNewParticle(Ns,PB,PBIGlobal)
        return
     End  Subroutine MCC
          
            subroutine SelectParticle(PB,CollisionRatio)
                    Implicit none
                    Type(ParticleBundle),intent(inout) :: PB
                    Real(8),intent(in) :: CollisionRatio
                    Integer(4) :: i,N
                    Real(8) ::  Residue
                    
                    Associate(NPar=>PB%NPar)
                            If(Allocated(MCCParticleIndex)) Deallocate(MCCParticleIndex)
                            Allocate(MCCParticleIndex(NPar))
                            If(Allocated(ParticleAnnihilationIndex)) Deallocate(ParticleAnnihilationIndex)
                            Allocate(ParticleAnnihilationIndex(NPar))
                            ParticleAnnihilationIndex=.False.
                            
                            NCollision=Int(PB%NPar*CollisionRatio)
                            Residue=PB%NPar*CollisionRatio-Dble(NCollision)
                            Call Random_Number(R)
                            If (R<Residue) Then
                                NCollision=NCollision+1
                            End If
                            Call  SelectIndex(NPar,NCollision,MCCParticleIndex)
                            
                            If(Allocated(MCCParticleBundle)) Deallocate(MCCParticleBundle)
                            Allocate(MCCParticleBundle(NCollision))
                            
                            do i=1, NCollision
                                 Call MCCParticleBundle(i)%Select(PB%PO(MCCParticleIndex(i)))
                            end do
                       ENd Associate
                    return
            end subroutine   SelectParticle
            
            subroutine SelectIndex(Nindex,Nselect,POIndex)
                    Implicit none
                    Integer(4),intent(in) :: Nindex,Nselect
                    Integer(4),intent(inout) :: POIndex(Nindex)
                    Integer :: Temp

                    Integer(4) :: i,j,k
                    POIndex=(/(i,i=1,Nindex)/)
                    Do i=1,Nselect
                          Call Random_Number(R)
                          k=(i+1)+(Nindex-(i+1))*R
                          Temp=POIndex(i)
                          POIndex(i)=POIndex(k)
                          POIndex(k)=Temp
                          
                    End Do
                    !qsort()
                    return
            end subroutine   SelectIndex
            
             Subroutine MCCAddNewParticle(Ns,PB,PBI)
                    Implicit none
                    Integer(4) :: Ns
                    Type(ParticleBundle),intent(inout) :: PB(0:Ns)
                    Type(ParticleBundleIndex),intent(inout) :: PBI
                    Integer(4) :: i,Index
                    
                    !Associate(NPar=>PBI%NPar,POI=>PBI%POI,Index=>PBI%POI(i)%Index)
                            Do i=1,PBI%NPar
                                 Index=PBI%POI(i)%Index
                                 Call PBI%POI(i)%VelRes(1.d0/PB(Index)%VFactor)
                                 Call PB(Index)%AddOne(PBI%POI(i)%ParticleOne)
                             End Do
                   !ENd Associate
                    return
            End subroutine MCCAddNewParticle
            
End Module ModuleMCCPublic 

