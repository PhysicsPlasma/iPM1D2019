!!! This version is finished by 2016/12/09, the code can run well and it can be a produciable code, but still improvements are possible, please refer next version.
! Note Due to the crosssections are different from the ones used in our code, there are some differences in the results.    
    
    Program PIC_MCC_for_CCP
   Use ModuleOneStep
   Use ModuleDiagOneStep
   use ModuleControlFlow
   Implicit none
   integer(4)  :: one[*]
   Integer(4) :: i,j,k
   integer(8) :: time_start, time_end
   
    !  Integer(4) :: NRun=1,NDiagShort=1,NDiagLong=0
   Integer(4) :: NRun=1000,NDiagShort=10,NDiagLong=0
   !Integer(4) :: NRun=10000,NDiagShort=200,NDiagLong=200
   
   imageRank = this_image()
   imageSize = num_images()
   if (1 == imageRank) then
        call system_clock(time_start)
    end if

    call random_seed()
   Call AllInitilalization()
   
   !Call CPU_TIME(CPU1)
   DO j=1,ControlFlowGlobal%NRun
        do i=1,ControlFlowGlobal%Period
                Call OneStep()
                call MergeAndSplit()   
        ENd DO
        sync all
        if (mod(j, 100) == 0) then
            Call OneStepRestart()
        end if
        if (1 == imageRank) then
            call system_clock(time_end)
            open (10,position='append',file='ParticleNumber.dat')
            write (10,*) ParticleLocal%NParAll,ParticleLocal%weight
            close (10)
            Write(*,*) 'Period ',j,ParticleLocal%NParAll
            Write(*, '(a30, f12.6)')       'CPU Time (s):         ', dble(time_end - time_start) / 1000000.d0
        end if
   ENd Do
   sync all
   !Call CPU_TIME(CPU2)
   !Write(*,*) 'Period ',CPU2-CPU1,ParticleGlobal%NPar
   !Write(*,*) 'Period ', j,' Complete!'  
   Call DiagInitilalization(ControlFlowGlobal)
     do j=1,ControlFlowGlobal%NDiagShort
         do i=1,ControlFlowGlobal%Period
             Call OneStep()
             Call DiagOneStep()
             call MergeAndSplit()
         End do
         if (1 == imageRank) then
            call system_clock(time_end)
            Write(*,*) 'Diag Period ',j,ParticleLocal%NParAll
            Write(*, '(a30, f12.6)')       'CPU Time (s):         ', dble(time_end - time_start) / 1000000.d0
        end if
     ENd Do
     sync all
     Call DiagOneStepFinal()
   !

stop
end  Program

