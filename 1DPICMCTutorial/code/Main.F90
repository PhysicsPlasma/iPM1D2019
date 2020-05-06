!program Hello_World
!  implicit none
!  integer :: i  ! Local variable
!  character(len=20) :: name[*] ! scalar coarray, one "name" for each image.
!  ! Note: "name" is the local variable while "name[<index>]" accesses the
!  ! variable in a specific image; "name[this_image()]" is the same as "name".
!
!  ! Interact with the user on Image 1; execution for all others pass by.
!  if (this_image() == 1) then   
!    write(*,'(a)',advance='no') 'Enter your name: '
!    read(*,'(a)') name
!
!    ! Distribute information to other images
!    do i = 2, num_images()
!      name[i] = name
!    end do
!  end if
!
!  sync all ! Barrier to make sure the data have arrived.
!
!  ! I/O from all images, executing in any order, but each record written is intact. 
!  write(*,'(3a,i0)') 'Hello ',trim(name),' from image ', this_image()
!  pause
!end program Hello_world
    
    !!!! This version is finished by 2016/12/09, the code can run well and it can be a produciable code, but still improvements are possible, please refer next version.
!! Note Due to the crosssections are different from the ones used in our code, there are some differences in the results.    
!    
    Program PIC_MCC_for_CCP
   Use ModuleOneStep
   Use ModuleDiagOneStep
   Implicit none
   Integer(4) :: i,j,k
   real(8) Cpu1,Cpu2
   Logical ::  Status
   integer :: img, nimgs
   
   
      Integer(4) :: NRun=1,NDiagShort=1,NDiagLong=0
   !Integer(4) :: NRun=10,NDiagShort=0,NDiagLong=0
   !Integer(4) :: NRun=10000,NDiagShort=200,NDiagLong=200
      Write(*,*) this_image(),"myid"  
      !pause
   Call AllInitilalizationCO()
!      img = this_image()
!   Write(*,*) "myid",img
!   
!
!   !DO j=1,NRun
!   !do i=1,ControlFlowGlobal%Period
!   !     Call OneStep()
!   !      If (ParticleGlobal(0)%Npar>ParticleGlobal(0)%NParNormal) then
!   !                  do k=0,1                
!   !                          Write(*,*) ParticleGlobal(k)%Npar,k,"before"
!   !                          Call ParticleBundleNormalization(ParticleGlobal(k),ParticleGlobal(k)%Npar/2)
!   !                          Write(*,*) ParticleGlobal(k)%Npar,k,"after" 
!   !                  end do
!   !              End If     
!   !ENd DO
!   
!   !open (10,position='append',file='ParticleNumber.dat')
!   !write (10,*) ParticleGlobal%NPar,ParticleGlobal%weight
!   !close (10)
!   
!   !do 
!   !Call OneStepRestart()
!   !Write(*,*) 'Period ',j,ParticleGlobal%NPar
!   
!
!   !
!   !
!   !!
      Call CPU_TIME(CPU1)
   Call ParticleGlobalCO(0)%Dump(0)
   Call CPU_TIME(CPU2)
   Write(*,*) 'Period ',CPU2-CPU1,'DumpTime1!!!',ParticleGlobalCO(0)%NPar
   
   
   Call CPU_TIME(CPU1)
   Call ParticleGlobalCO(0)%Load(Status)
   Call CPU_TIME(CPU2)
   Write(*,*) 'Period ',CPU2-CPU1,'LoadTime1!!!',ParticleGlobalCO(0)%NPar

   !pause

   !pause
   
      
   Call CPU_TIME(CPU1)
   Call ParticleGlobalCO(0)%Dump(1)
   Call CPU_TIME(CPU2)
   Write(*,*) 'Period ',CPU2-CPU1,'DumpTime2!!!',ParticleGlobalCO(0)%NPar
   
   Call CPU_TIME(CPU1)
   Call ParticleGlobalCO(0)%Load(Status)
   Call CPU_TIME(CPU2)
   Write(*,*) 'Period ',CPU2-CPU1,'LoadTime2!!!',ParticleGlobalCO(0)%NPar

   !
   !!ENd Do
   !Call CPU_TIME(CPU1)
   !Call ParticleGlobal(0)%Dump(0)
   !Call CPU_TIME(CPU2)
   !Write(*,*) 'Period ',CPU2-CPU1,'DumpTime1!!!',ParticleGlobal(0)%NPar
!   !
!   !
!   !
!   !Call CPU_TIME(CPU1)
!   !Call ParticleGlobal(0)%Load(Status)
!   !Call CPU_TIME(CPU2)
!   !Write(*,*) 'Period ',CPU2-CPU1,'LoadTime1!!!',ParticleGlobal(0)%NPar
!   !
!   !   Call CPU_TIME(CPU1)
!   !Call ParticleGlobal(0)%Dump(0)
!   !Call CPU_TIME(CPU2)
!   !Write(*,*) 'Period ',CPU2-CPU1,'DumpTime2!!!',ParticleGlobal(0)%NPar
!   !
!   !Call CPU_TIME(CPU1)
!   !Call ParticleGlobal(0)%Load(Status)
!   !Call CPU_TIME(CPU2)
!   !Write(*,*) 'Period ',CPU2-CPU1,'LoadTime2!!!',ParticleGlobal(0)%NPar
!   
!   
!   
!   !Write(*,*) 'Period ', j,' Complete!'  
!   !Call DiagInitilalization(ControlFlowGlobal)
!   !  do j=1,NDiagShort
!   !      do i=1,ControlFlowGlobal%Period
!   !          Call OneStep()
!   !          Call DiagOneStep()
!   !      End do
!   !  ENd Do
!   !  Call DiagOneStepFinal()
!   !!
!
pause
stop
end  Program
!
