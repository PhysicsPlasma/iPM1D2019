Module ModuleField
     Use ModuleControlFlow
     Implicit none
                  Type Field !(Nx)
                        Integer(4) :: Nx=NxMax
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        Real(8) ::  Ex(1:NxMax)=0.d0,Ey(1:NxMax)=0.d0,Ez(1:NxMax)=0.d0
                        Real(8) ::  Bx(1:NxMax)=0.d0,By(1:NxMax)=0.d0,Bz(1:NxMax)=0.d0
                        Real(8) ::  Rho(1:NxMax),Phi(1:NxMax)
                        Real(8) ::  Chi(1:NxMax)
                        contains
                             procedure :: Dump=>DumpField
                             procedure :: Load=>LoadField
                   EndType Field
                   
                    Type FieldSolver
                        Integer(4) :: Ns=NxMax-2
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        Real(8) :: Source(1:NxMax-2)
                        Real(8) :: Solve(1:NxMax-2)
                        Real(8) :: CoeA(1:NxMax-2),CoeB(1:NxMax-2),CoeC(1:NxMax-2)
                        contains
                             procedure :: Dump=>DumpFieldSolver
                    EndType FieldSolver

                    Type FieldOne
                        Integer(4) :: Nx=NxMax
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        Real(8) ::  RhoOne(1:NxMax),ChiOne(1:NxMax)
                    EndType FieldOne
 
    contains
              subroutine DumpField(FG, name)
                     Implicit none

                     Class(Field), intent(inout) :: FG
                     Character(len=CharLenthMax), intent(in), optional :: name
                     Character(len=CharLenthMax) :: filename
                     Integer(4):: i
                     
                     if (1 == imageRank) then
                            Write(filename,*) "./restart/","FieldGlobal",".dat"
                            if (present(name)) filename = name
                            Write(*,*) "Saving ",filename," Please wait..."
                            open (10,file=filename)
                                   do i=1,FG%Nx
                                          Write(10,FMt="(*(es21.14,1x))") dble(i-1)*FG%Dx,FG%Ex(i),FG%Phi(i),FG%Rho(i),FG%Chi(i)
                                   end do
                            close(10)
                            Write(*,*) "Save ",filename,"Complete!"
                     end if

                     return
              end subroutine DumpField
            ! 
              subroutine LoadField(FG,Status,name)
			Implicit none
			Class(Field),intent(inout) :: FG
			Logical,intent(inout) ::  Status
			Character(len=CharLenthMax), intent(in), optional :: name
			Character(len=CharLenthMax) :: filename
			logical :: alive
			Integer(4) :: i,j
			Real(8) :: TempDx

			Write(filename,*) "./restart/","FieldGlobal",".dat"
			if (present(name)) filename = name

			Inquire(file=filename,exist=alive)
			If(alive) then
				if (1 == imageRank) then
					Write(*,*) "Loading ",filename," Please wait..."
				end if

				open (10,file=filename)
				do i=1,FG%Nx
					Read(10,*) TempDx,FG%Ex(i),FG%Phi(i),FG%Rho(i),FG%Chi(i)
				end do

				close(10)
				Status=.False.
				if (1 == imageRank) then
					Write(*,*) "Load ",filename,"Complete!"
				end if
			else
				if (1 == imageRank) then
					Write(*,*) 'Can not find the file for ', filename,' set for zero.'
				end if
				
				FG%Ex=0.d0
				FG%Phi=0.d0
				FG%Rho=0.d0
				FG%Chi=0.d0
				Status=.True. 
			End if    
			return
		end subroutine LoadField

             subroutine DumpFieldSolver(FS,name)
                     Implicit none
                     Class(FieldSolver),intent(inout) :: FS
                     Character(len=CharLenthMax), intent(in), optional :: name
                     Character(len=CharLenthMax) :: filename
                     Integer(4):: i
                     if (1 == imageRank) then
			       Write(filename,*) "./restart/","FieldSolver",".dat"
				if (present(name)) filename = name
                            Write(*,*) "Saving ",filename," Please wait..."
                            open (10,file=filename)
                    !Write(10,*)  FS%Ns,FS%Dx,FS%Dt
                                   do i=1,FS%Ns
                                          Write(10,FMt="(*(es21.14,1x))") dble(i-1)*FS%Dx,FS%Solve(i),FS%Source(i),FS%CoeA(i),FS%CoeB(i),FS%CoeC(i)!,FG%Bx(i),FG%By(i)
                                   end do
                            close(10)
                    Write(*,*) "Save ",filename,"Complete!"  
                     end if
                    return
             end subroutine DumpFieldSolver
             
             subroutine DumpFieldOne(Ns,FO,rank,name)
                     Implicit none
                     Integer(4),intent(in) :: Ns
                     Type(FieldOne),intent(inout) :: FO(0:Ns)
                     Integer(4), intent(in), optional :: rank
                     Character(len=CharLenthMax), intent(in), optional :: name
                     Character(len=CharLenthMax) :: filename
                     Integer(4):: i,j

                     if (present(rank) .and. rank >= 1) then
                            Write(filename, '(a, a, i3, a)') "./restart/","FieldOne ", rank, ".dat"
                            if (present(name)) filename = name
                            Write(*,*) "Saving ",filename," Please wait..."
                            open (10,file=filename)
                                   do i=1,FO(0)%Nx
                                          Write(10,FMt="(*(es21.14,1x))") FO(0)%Dx*dble(i-1),(ABS(FO(j)%RhoOne(i)/ElectronCharge),j=0,Ns)
                                   end do
                            close(10)
                            Write(*,*) "Save ",filename,"Complete!"
                     else
                            if (1 == imageRank) then
                                   Write(filename,*) "./restart/","FieldOne",".dat"
                                   if (present(name)) filename = name
                                   Write(*,*) "Saving ",filename," Please wait..."
                                   open (10,file=filename)
                                          do i=1,FO(0)%Nx
                                                 Write(10,FMt="(*(es21.14,1x))") FO(0)%Dx*dble(i-1),(ABS(FO(j)%RhoOne(i)/ElectronCharge),j=0,Ns)
                                          end do
                                   close(10)
                                   Write(*,*) "Save ",filename,"Complete!"
                            end if
                     end if
                     return
              end subroutine DumpFieldOne     
           
 End Module ModuleField 