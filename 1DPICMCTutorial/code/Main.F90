!!! This version is finished by 2016/12/09, the code can run well and it can be a produciable code, but still improvements are possible, please refer next version.
! Note Due to the crosssections are different from the ones used in our code, there are some differences in the results.    
    
    Program PIC_MCC_for_CCP
   Use ModuleOneStep
   Use ModuleDiagOneStep
   Implicit none
   Integer(4) :: i,j,k
   real(8) Cpu1,Cpu2
   Logical ::  Status
   
      Integer(4) :: NRun=1,NDiagShort=1,NDiagLong=0
   !Integer(4) :: NRun=10,NDiagShort=0,NDiagLong=0
   !Integer(4) :: NRun=10000,NDiagShort=200,NDiagLong=200
      
   Call AllInitilalization()
   

   !DO j=1,NRun
   !do i=1,ControlFlowGlobal%Period
   !     Call OneStep()
   !      If (ParticleGlobal(0)%Npar>ParticleGlobal(0)%NParNormal) then
   !                  do k=0,1                
   !                          Write(*,*) ParticleGlobal(k)%Npar,k,"before"
   !                          Call ParticleBundleNormalization(ParticleGlobal(k),ParticleGlobal(k)%Npar/2)
   !                          Write(*,*) ParticleGlobal(k)%Npar,k,"after" 
   !                  end do
   !              End If     
   !ENd DO
   
   !open (10,position='append',file='ParticleNumber.dat')
   !write (10,*) ParticleGlobal%NPar,ParticleGlobal%weight
   !close (10)
   
   !do 
   !Call OneStepRestart()
   !Write(*,*) 'Period ',j,ParticleGlobal%NPar
   

   !
   !
   !
   Call CPU_TIME(CPU1)
   Call ParticleGlobal(0)%Load(Status)
   Call CPU_TIME(CPU2)
   Write(*,*) 'Period ',CPU2-CPU1,'LoadTime1!!!',ParticleGlobal(0)%NPar
   
   If(Status) Write(*,*) "Load ERORRRR!!!!!!"
   
   Call CPU_TIME(CPU1)
   Call ParticleGlobal(0)%Dump(0)
   Call CPU_TIME(CPU2)
   Write(*,*) 'Period ',CPU2-CPU1,'DumpTime1!!!',ParticleGlobal(0)%NPar
   

   
   Call CPU_TIME(CPU1)
   Call ParticleGlobal(0)%Load(Status)
   Call CPU_TIME(CPU2)
   Write(*,*) 'Period ',CPU2-CPU1,'LoadTime2!!!',ParticleGlobal(0)%NPar
   
   Call CPU_TIME(CPU1)
   Call ParticleGlobal(0)%Dump(0)
   Call CPU_TIME(CPU2)
   Write(*,*) 'Period ',CPU2-CPU1,'DumpTime2!!!',ParticleGlobal(0)%NPar
   
   !ENd Do
   !Call CPU_TIME(CPU1)
   !Call ParticleGlobal(0)%Dump(0)
   !Call CPU_TIME(CPU2)
   !Write(*,*) 'Period ',CPU2-CPU1,'DumpTime1!!!',ParticleGlobal(0)%NPar
   !
   !
   !
   !Call CPU_TIME(CPU1)
   !Call ParticleGlobal(0)%Load(Status)
   !Call CPU_TIME(CPU2)
   !Write(*,*) 'Period ',CPU2-CPU1,'LoadTime1!!!',ParticleGlobal(0)%NPar
   !
   !   Call CPU_TIME(CPU1)
   !Call ParticleGlobal(0)%Dump(0)
   !Call CPU_TIME(CPU2)
   !Write(*,*) 'Period ',CPU2-CPU1,'DumpTime2!!!',ParticleGlobal(0)%NPar
   !
   !Call CPU_TIME(CPU1)
   !Call ParticleGlobal(0)%Load(Status)
   !Call CPU_TIME(CPU2)
   !Write(*,*) 'Period ',CPU2-CPU1,'LoadTime2!!!',ParticleGlobal(0)%NPar
   
   
   
   !Write(*,*) 'Period ', j,' Complete!'  
   !Call DiagInitilalization(ControlFlowGlobal)
   !  do j=1,NDiagShort
   !      do i=1,ControlFlowGlobal%Period
   !          Call OneStep()
   !          Call DiagOneStep()
   !      End do
   !  ENd Do
   !  Call DiagOneStepFinal()
   !!

pause
stop
end  Program

