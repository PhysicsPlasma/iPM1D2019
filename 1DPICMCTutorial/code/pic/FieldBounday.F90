Module ModuleFieldBoundary
  Use ModuleControlFlow
  Use ModuleField
  Use ModuleParticleBoundary
  !Use ModuleParticleBoundary
  Implicit none
     Type FieldBoundary
                        !Integer(4) :: XStart=0,XEnd=1
                         Integer(4) :: FieldBoundaryModel=11
                         Integer(4) :: Timer=0,Period
                         Real(8) :: Dt=Inputdt
                         Real(8) :: Frequency(2)=(/13.56d6,2.d6/),Voltage(2)=(/200.d0,100.d0/)
                         Real(8) :: V1=0.d0,V2=0.d0
                         Real(8) :: Vdc=0.d0!-20.d0
                         Contains
                         procedure :: Init=>InitilalizationFieldBounday
                         procedure :: Updater=>UpdaterFieldBounday
     
     EndType FieldBoundary

    contains
Subroutine InitilalizationFieldBounday(FB,CF,name)
			Implicit none
			Class(FieldBoundary),intent(inout) :: FB
			Type(ControlFlow),intent(inout) :: CF
			Character(len=CharLenthMax), intent(in), optional :: name
			Real(8) :: Dt
			integer(4) :: i
			Character(len=CharLenthMax) :: Filename
			Logical :: alive

			Write(filename,*) "./restart/","FieldBoundary",".dat"
			if (present(name)) filename = name

			Inquire(file=filename,exist=alive)
			If(alive) then
				Open (10,file=filename)
					Read(10,*) FB%FieldBoundaryModel
					Read(10,*) FB%Timer,FB%Period
					Read(10,*) FB%V1,FB%V2,FB%Vdc
					Read(10,*) (FB%Frequency(i),i=1,3)
					Read(10,*) (FB%Voltage(i),i=1,3)
					Read(10,*) FB%V1
					Read(10,*) FB%V2
					Read(10,*) FB%Vdc
				Close(10)
			Else
				if (1 == imageRank) then
					Write(*,*) 'Can not find the file for ', filename,' the FieldBoundary will be set to Zero.'
				end if

				FB%Timer=0
				FB%Vdc=0.d0
				FB%V1=0.d0
				FB%V2=0.d0
			End If

			CF%Period = Int(1.d0/FB%Frequency(1)/CF%dt)
			CF%Timer = FB%Timer

			if (1 == imageRank) then
				write(*,*) "CF%Period=",CF%Period
			end if
			
			FB%Dt = CF%dt
			FB%Period = CF%Period

			return
		End Subroutine InitilalizationFieldBounday
     
     Subroutine UpdaterFieldBounday(FB)
            implicit none
            Class(FieldBoundary),intent(inout) :: FB

            Select case (FB%FieldBoundaryModel)
                      case (11)
                               FB%V1=FB%Voltage(1)*DSin(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))
                               FB%V2=0.d0
                               !Write (*,*) FB%V1,FB%V2
                      case (21)
                               FB%V1=FB%Voltage(1)*DCos(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))+FB%Voltage(2)*DCos(2*PI*FB%Frequency(2)*FB%dt*Dble(FB%Timer))
                               FB%V2=0.d0
                      case(22)
                               FB%V1=FB%Voltage(2)*DCos(2*PI*FB%Frequency(2)*FB%dt*Dble(FB%Timer))
                               FB%V2=FB%Voltage(1)*DCos(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))  
                      Case(31)
                              
              end  Select
  
            FB%Timer=FB%Timer+1
            return
    end subroutine  UpdaterFieldBounday

     
     Subroutine FieldBoundayFinalization(FB,name)
        Implicit none
        class(FieldBoundary),intent(in) :: FB
        Character(len=CharLenthMax), intent(in), optional :: name
        Character(len=CharLenthMax) :: Filename
        integer(4) :: i

        if (1 == imageRank) then
            Write(filename,*) "./restart/","FieldBoundary",".dat"
            if (present(name)) filename = name

            Open (10,file=filename)
                write(10,*) FB%FieldBoundaryModel
                write(10,*) FB%Timer,FB%Period
                write(10,*) FB%V1,FB%V2,FB%Vdc
                write(10,*) (FB%Frequency(i),i=1,3)
                write(10,*) (FB%Voltage(i),i=1,3)
                write(10,*) FB%V1
                write(10,*) FB%V2
                write(10,*) FB%Vdc
            Close(10)
        end if

        return
     End Subroutine FieldBoundayFinalization
     !
     Subroutine UpdateFieldBounday(FB,Ns,PBDO)
               Implicit none
               Type(FieldBoundary),intent(inout) :: FB
               Integer(4), intent(in) :: Ns
               Type(ParticleBoundaryOne),intent(inout) :: PBDO(0:Ns)
               Integer(4) :: i
               Real(8) :: NetCharge
               If(Mod(FB%Timer,FB%Period)==0) Then
                   NetCharge=0.d0
                   do i=0,Ns
                        NetCharge=NetCharge+PBDO(i)%CountMin*PBDO(i)%PBLower%SO%Charge*PBDO(i)%PBLower%Weight
                        PBDO(i)%CountMin=0
                   End Do
                   If (NetCharge<0.d0) then
                       FB%Vdc=FB%Vdc-1.d0
                   Else if (NetCharge>0.d0) Then
                       FB%Vdc=FB%Vdc+1.d0
                   ENd IF
               End IF
                  return
     End Subroutine UpdateFieldBounday
     
end  Module ModuleFieldBoundary

              
              