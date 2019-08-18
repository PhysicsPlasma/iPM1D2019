!Module RecombinationModule
!     Use ParticleModule
!     Implicit none
!     contains
!       Subroutine  IonIonRecombinationRate(Nx,ParticlePositive,ParticleNegative,Weighting,Volume,dt)
!            Implicit none
!            Type(ParticleBundle),intent(inout) :: ParticlePositive,ParticleNegative
!            Real(8),intent(in) :: Weighting,Volume,dt
!            Real(8),parameter :: K0=2.7d-13 !G文献中 都是1d-13
!            Integer(4) :: i,Nx,N
!            Real(8) :: NumberPositive(Nx), NumberNegative(Nx),NumberRecombination(Nx),NumberRecombinationTemp(Nx)
!            
!            NumberPositive=0.d0
!            NumberNegative=0.d0
!            NumberRecombination=0.d0
!            NumberRecombinationTemp=0.d0
!            !Write (*,*)  ParticlePositive%Name,ParticlePositive%NPar,ParticleNegative%Name,ParticleNegative%NPar,'Before'
!            do i=1,ParticlePositive%NPar
!                       N=Ceiling(ParticlePositive%PhaseSpace(i)%X)
!                       NumberPositive(N)=NumberPositive(N)+1.d0   !求出每个网格内粒子数 N
!            End do
!           
!            do i=1,ParticleNegative%NPar
!                       N=Ceiling(ParticleNegative%PhaseSpace(i)%X) !同上
!                       NumberNegative(N)=NumberNegative(N)+1.d0 
!            End do
!           
!            NumberRecombination(1:Nx)=Weighting*K0/Volume*dt*NumberPositive(1:Nx)*NumberNegative(1:Nx) !Na*w*k*dt/V*Nb, Na*w*k*dt/V=n*v*sigame 注意(k=v*sigma)
!            
!            NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx)
!            Call Recombination(ParticlePositive,Nx,NumberRecombinationTemp) !处理复合的正离子-del
!            
!            NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx) !因为上面的调用函数已经改变了NumberRecombinationTemp值，需要重新赋值
!            Call Recombination(ParticleNegative,Nx,NumberRecombinationTemp)  !处理复合的负离子-del
!            !Write (*,*)  ParticlePositive%Name,ParticlePositive%NPar,ParticleNegative%Name,ParticleNegative%NPar,'After'       
!           return
!       End Subroutine  IonIonRecombinationRate
!          
!       Subroutine  Recombination(InputParticle,Nx,NumberRecombination)
!               Implicit none
!               Integer(4),intent(in) ::  Nx
!               Type(ParticleBundle),intent(inout) :: InputParticle
!               Real(8),intent(inout) :: NumberRecombination(Nx)
!               Integer(4) :: i,N,IndexBegin
!               Real(8) :: Ratio
!               
!                   CALL RANDOM_NUMBER(R)
!                   IndexBegin=Ceiling(R*dble(InputParticle%Npar))  !随机一个遍历起始点
!                     
!                   do i=IndexBegin,1,-1   !这里是分成两部分，将数组都遍历了一下，---改进：可以做一个总数，附加判断条件一旦总数达到 就不进行进一步的判断了
!                       N=Ceiling(InputParticle%PhaseSpace(i)%X) 
!                       If(NumberRecombination(N)>MinReal) Then
!                           If(NumberRecombination(N)>=1.d0)  Then
!                               NumberRecombination(N)=NumberRecombination(N)-1.d0
!                               Call DelParticle(i,InputParticle)
!                            Else
!                                 CALL RANDOM_NUMBER(R)
!                                 Ratio=NumberRecombination(N)
!                                 NumberRecombination(N)=NumberRecombination(N)-1.d0
!                                 If(R<Ratio)  Then
!                                    Call DelParticle(i,InputParticle)
!                                 End If 
!                             End if
!                        End if
!                 End do
!                 
!                i=IndexBegin+1
!                do while (i<=InputParticle%NPar)
!                       N=Ceiling(InputParticle%PhaseSpace(i)%X) 
!                       If(NumberRecombination(N)>MinReal) Then
!                           If(NumberRecombination(N)>=1.d0)  Then
!                               NumberRecombination(N)=NumberRecombination(N)-1.d0
!                               Call DelParticle(i,InputParticle)
!                            Else
!                                 CALL RANDOM_NUMBER(R)
!                                 Ratio=NumberRecombination(N)
!                                 NumberRecombination(N)=NumberRecombination(N)-1.d0
!                               If(R<Ratio)  Then
!                                 Call DelParticle(i,InputParticle)
!                                End If 
!                             End if
!                        End if
!                          i=i+1
!                 End do      
!               return
!          End Subroutine  Recombination
!
!  subroutine ElectronIonRecombinationRate(Nx,ParticlePositive,ParticleNegative,Weighting,Volume,dt) !这个函数没有进行处理，一般这个反应确实可以忽略，电子太少了
!        implicit none
!        Type(ParticleBundle),intent(inout) :: ParticlePositive,ParticleNegative
!            Real(8),intent(in) :: Weighting,Volume,dt
!            Real(8),parameter :: K0=1d-13
!            Integer(4) :: i,Nx,N
!            Real(8) :: NumberPositive(Nx), NumberNegative(Nx),NumberRecombination(Nx),NumberRecombinationTemp(Nx)
!     
!
!        return
!   end subroutine ElectronIonRecombinationRate
!End Module RecombinationModule