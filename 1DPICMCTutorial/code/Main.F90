!!! This version is finished by 2016/12/09, the code can run well and it can be a produciable code, but still improvements are possible, please refer next version.
! Note Due to the crosssections are different from the ones used in our code, there are some differences in the results.    
    
    Program PIC_MCC_for_CCP
   Use ModuleOneStep
   Use ModuleDiagOneStep
   Implicit none
   Integer(4) :: i,j,k
   real(8) Cpu1,Cpu2,CPU3
   
   Integer(4) :: NRun=100000,NDiagShort=1000,NDiagLong=0,Per=100
      
   Call AllInitilalization()
   
   Call CPU_TIME(CPU1)
   DO j=1,NRun
        ! do i=1,ControlFlowGlobal%Period
        do i=1,Per
            Call OneStep()
            ! If (ParticleGlobal(0)%Npar>ParticleGlobal(0)%NParNormal) then
            If (ParticleGlobal(1)%Npar>2*ControlFlowGlobal%ParticlePerGrid*(NxMax-1) .or. ParticleGlobal(0)%Npar>2*ControlFlowGlobal%ParticlePerGrid*(NxMax-1)) then
                do k=0,1                
                        Write(*,*) ParticleGlobal(k)%Npar,k,"before"
                        Call ParticleBundleNormalization(ParticleGlobal(k),ParticleGlobal(k)%Npar/2)
                        Write(*,*) ParticleGlobal(k)%Npar,k,"after" 
                end do
            End If
        ENd DO
        open (10,position='append',file='ParticleNumber.dat')
        write (10,*) ParticleGlobal%NPar,ParticleGlobal%weight
        close (10)
        if ( mod(j,100)==0 ) then
            Call OneStepRestart()
            Call CPU_TIME(CPU2)
            write(*,'(*(A,I3))') 'time:',floor((CPU2-CPU1)/3600),'h',floor(((CPU2-CPU1)-floor((CPU2-CPU1)/3600)*3600)/60),'m',mod(int((CPU2-CPU1)),60),'s'
        end if
        !Call OneStepRestart()
        Write(*,*) 'Period ',j,ParticleGlobal%NPar
   ENd Do
   call CPU_TIME(Cpu2)
      !Write(*,*) 'Period ',CPU2-CPU1,ParticleGlobal%NPar
   !Write(*,*) 'Period ', j,' Complete!'  
   Call DiagInitilalization(ControlFlowGlobal)
   
    !  do j=1,NDiagShort
    !     do i=1,Per
    !     !do i=1,Per
    !          Call OneStep()
    !          Call DiagOneStep()
    !      End do
    !      Write(*,*) 'DiagPeriod ',j
    !      if ( mod(j,100)==0 ) then
    !         Call CPU_TIME(CPU2)
    !         write(*,'(*(A,I3))') 'time:',floor((CPU2-CPU1)/3600),'h',floor(((CPU2-CPU1)-floor((CPU2-CPU1)/3600)*3600)/60),'m',mod(int((CPU2-CPU1)),60),'s'
    !     end if
    !  ENd Do

    ! ControlFlowGlobal%Period=100000
   WRITE(*,*) 'diag total period:',ControlFlowGlobal%Period,',and diag start now:'
    do i=1,ControlFlowGlobal%Period
        Call OneStep()
        Call DiagOneStep()
        if ( mod(i,1000)==0 ) then
            Write(*,*) 'DiagPeriod ',i
            Call CPU_TIME(CPU3)
            write(*,'(*(A,I3))') 'time:',floor((CPU3-CPU1)/3600),'h',floor(((CPU3-CPU1)-floor((CPU3-CPU1)/3600)*3600)/60),'m',mod(int((CPU3-CPU1)),60),'s'
        end if
    End do        
        
     Call DiagOneStepFinal()

     Call CPU_TIME(CPU3)
     write(*,'(*(A,I3))') 'Total time:',floor((CPU3-CPU1)/3600),'h',floor(((CPU3-CPU1)-floor((CPU3-CPU1)/3600)*3600)/60),'m',mod(int((CPU3-CPU1)),60),'s'
     write(*,'(*(A,I3))') 'Run time:',floor((CPU2-CPU1)/3600),'h',floor(((CPU2-CPU1)-floor((CPU2-CPU1)/3600)*3600)/60),'m',mod(int((CPU2-CPU1)),60),'s'
     write(*,'(*(A,I3))') 'Diag time:',floor((CPU3-CPU2)/3600),'h',floor(((CPU3-CPU2)-floor((CPU3-CPU2)/3600)*3600)/60),'m',mod(int((CPU3-CPU2)),60),'s'

pause
stop
end  Program

