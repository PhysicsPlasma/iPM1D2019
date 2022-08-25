Module ModuleGrid
    Use ModuleControlFlow
    Use Numrical
    Implicit none
    !  This Section is the definitions of the physicla parameters.
             Type :: Grid1D(Nx,Ns)
                        Integer(4),Len :: Nx=NxMax,Ns=1_4
                        Real(8) ::  Value(1:Nx,1:Ns)
                        Integer(4) :: Timer=0,Period=1
                        !Integer(4) :: NxAva=1
                        Real(8) :: Dx=0.d0
             contains
                        procedure :: Init=>InitializationGrid1D
                        procedure :: Dump=>DumpGrid1D
                        procedure :: Update=>UpdateGrid1D
                        procedure :: Rescale=>RescaleGrid1D
                        procedure :: Reset=>ResetGrid1D
                        
             End Type Grid1D
                   
                Type :: Grid2D(Nx,Ny,Ns)
                        Integer(4),Len :: Nx=NxMax,Ny=1000,Ns=1_4
                        Real(8) ::  Value(1:Nx,0:Ny,1:Ns)
                         Real(8) :: Dx=0.d0,Dy=0.d0
                         Integer(4) :: Timer=0,Period=1
                         Integer(4) :: NxAva=1,NyAva=1
                       contains
                        procedure :: Init=>InitializationGrid2D
                        procedure :: Dump=>DumpGrid2D
                        procedure :: Update=>UpdateGrid2D
                        procedure :: Rescale=>RescaleGrid2D
                        procedure :: Reset=>ResetGrid2D
                 End Type Grid2D

    Contains
    
                            
    
    subroutine InitializationGrid1D(GD,CF)
        Implicit none
        Class(Grid1D(*,*)),intent(inout)  :: GD
        Class(ControlFlow), intent(in) :: CF
                    GD%Timer=0
                    GD%Period=CF%Period
                    GD%Dx=CF%Dx
                    GD%Value=0.d0
        return
    end subroutine InitializationGrid1D
    
    
          
     subroutine DumpGrid1D(GD,Mode,rank)
        Implicit none
        Class(Grid1D(*,*)),intent(inout)  :: GD
        Integer(4),intent(in) :: Mode
        Integer(4),intent(in), optional :: rank
        Character(len=99) :: Filename
        Integer(4) :: i,j,k,NameIndex,Ns,dumpRank=1
        Integer(4),save :: NameTimer=1
        if (present(rank)) dumpRank = rank
        If (Mode==0) Then
                   NameIndex=DefaultNameIndex
        else
                   NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
                   NameTimer=NameTimer+1
        End If 
        if (imageRank == dumpRank) then
                Write(filename,*) "Grid1D",NameIndex,".dat"
                Write(*,*) "Saving ",Filename," Please wait..."
                open (10,file=filename)
                Ns=GD%Ns
                do i=1,GD%Nx
                       Write(10,FMt="(*(es21.14,1x))")  dble(i-1)*GD%dx,(GD%Value(i,k),k=1,Ns)
                end do
                close(10)
                Write(*,*) "Save ",Filename,"Complete!"
        end if
        return
     end subroutine DumpGrid1D
     
     Subroutine UpdateGrid1D(GD,Nx,A1D,Shift)
         Implicit none
         Class(Grid1D(*,*)),intent(inout)  :: GD
         Integer(4),intent(in) :: Nx
         Real(8),intent(in)  :: A1D(Nx)
         Real(8) :: TempValue(1:GD%Nx)
         Integer(4),intent(inout) :: Shift
         
            If (Shift>=1.and.Shift<=GD%Ns) Then
                Call AnyAvg(A1D(1:Nx), TempValue(1:GD%Nx))
                GD%Value(1:GD%Nx,Shift)=GD%Value(1:GD%Nx,Shift)+TempValue(1:GD%Nx)
            End If
             Shift=Shift+1
         Return    
     End Subroutine  UpdateGrid1D
     
     Subroutine RescaleGrid1D(GD)
         Implicit none
         Class(Grid1D(*,*)),intent(inout)  :: GD
                 If (GD%Timer>=1) Then
                       GD%Value=GD%Value/dble(GD%Timer)
                 End If
         Return    
     End Subroutine  RescaleGrid1D
     
    subroutine ResetGrid1D(GD)
        Implicit none
        Class(Grid1D(*,*)),intent(inout)  :: GD
                      GD%Value=0.d0
                      GD%Timer=0
        return
     end subroutine ResetGrid1D
     
      subroutine InitializationGrid2D(GD,CF)
        Implicit none
        Class(Grid2D(*,*,*)),intent(inout)  :: GD
        Class(ControlFlow), intent(in) :: CF
            GD%Timer=0
            GD%Period=CF%Period
            GD%Dx=CF%Dx
            GD%Dy=2.d0/Dble(GD%Ny)!FG%Dt
            GD%Value=0.d0
        return
     end subroutine InitializationGrid2D
     
    subroutine DumpGrid2D(GD,Mode,rank)
        Implicit none
        Class(Grid2D(*,*,*)),intent(inout)  :: GD
        Integer(4),intent(in) :: Mode
        Integer(4),intent(in), optional :: rank
        Character(len=99) :: Filename
        Integer(4) :: i,j,k,NameIndex,Ns,dumpRank=1
        Integer(4),save :: NameTimer=1
        if (present(rank)) dumpRank = rank
        If (Mode==0) Then
                   NameIndex=DefaultNameIndex
        else
                   NameIndex=DefaultNameIndexInit+ModeMultiplier*Mode+NameTimer
                   NameTimer=NameTimer+1
        End If 
        if (imageRank == dumpRank) then
                Write(filename,*) "Grid2D",NameIndex,".dat"
                Filename=Trim(filename)
                Write(*,*) "Saving ",Filename," Please wait..."
                open (10,file=filename)
                Ns=GD%Ns
                do i=0,GD%Ny
                     do j=1,GD%Nx
                       Write(10,FMt="(*(es21.14,1x))")  dble(j-1)*GD%dx,dble(i)*GD%dy,(GD%Value(j,i,k),k=1,Ns)
                     End do
                end do
                close(10)
                Write(*,*) "Save ",Filename,"Complete!"
        end if
         return
     end subroutine DumpGrid2D
     

        Subroutine UpdateGrid2D(GD,Nx,A1D,Shift)
         Implicit none
         Class(Grid2D(*,*,*)),intent(inout)  :: GD
         Integer(4),intent(in) :: Nx
         Real(8),intent(in)  :: A1D(Nx)
         Integer(4),intent(inout) :: Shift
         Integer(4) :: NZip,N0,N1,Index
                 N0=Mod(GD%Timer,GD%Period)
                 If (GD%Period>GD%Ny) Then
                         NZip=GD%Period/GD%Ny
                         N1=NZip*GD%Ny
                        If (Shift>=1.and.Shift<=GD%Ns) Then
                                 If (N0<N1) Then
                                    Index=N0/NZip
                                    Call AnyAvg(A1D(1:Nx), GD%Value(1:GD%Nx,Index,Shift))
                                Else
                                    Index=GD%Ny
                                    Call AnyAvg(A1D(1:Nx), GD%Value(1:GD%Nx,Index,Shift))
                                ENd If
                        End If
                 Else
                     Call AnyAvg(A1D(1:Nx), GD%Value(1:GD%Nx,Index,Shift))
                 End If
             Shift=Shift+1
         Return    
     End Subroutine  UpdateGrid2D
     
     Subroutine RescaleGrid2D(GD)
         Implicit none
         Class(Grid2D(*,*,*)),intent(inout)  :: GD
         Real(8) :: Spaces,Timers
                 If (GD%Timer>=1) Then
                     Timers=DBLE(GD%Timer/GD%Period+1)
                     If (GD%Period>GD%Ny) Then
                         Spaces=dble(GD%Period/GD%Ny)
                         GD%Value(:,0:GD%Ny-1,:)=GD%Value(:,0:GD%Ny-1,:)/(Spaces*Timers)
                         Spaces=dble(GD%Period-GD%Period/GD%Ny*GD%Ny)
                         GD%Value(:,GD%Ny,:)=GD%Value(:,GD%Ny,:)/(Spaces*Timers)
                     Else
                         Spaces=1.d0
                         
                         GD%Value(:,0:GD%Ny,:)=GD%Value(:,0:GD%Ny,:)/(Spaces*Timers)
                     End If
                 End IF
         Return    
     End Subroutine  RescaleGrid2D
     
      subroutine ResetGrid2D(GD)
        Implicit none
        Class(Grid2D(*,*,*)),intent(inout)  :: GD
                      GD%Value=0.d0
                      GD%Timer=0
        return
     end subroutine ResetGrid2D
ENd Module ModuleGrid          
