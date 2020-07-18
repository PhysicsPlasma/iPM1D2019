Module ModuleVector1DX
   implicit none
contains

   subroutine Interpolation1DX(Results, Raw)
      Implicit none
      Real(8), intent(in)  :: Raw(:)
      Real(8), intent(out)  :: Results(:)
      Integer(4) :: i, Ni, No

      Ni = Size(Raw)
      No = Size(Results)

      If (No == (Ni - 1)) Then   !!!!!!!!! This part is for Normal condition without extraplotation
         do i = 1, No
            Results(i) = (Raw(i + 1) + Raw(i))*0.5d0
         End do

      Else If (No == Ni) Then  !!!!!! This part is for Periodical condition
         do i = 1, No - 1
            Results(i) = (Raw(i + 1) + Raw(i))*0.5d0
         End do
         Results(No) = Results(1)

      Else If (No == Ni + 1) Then !!!!! This part is for Normal condition with extraplotation
         do i = 2, No - 1
            Results(i) = (Raw(i) + Raw(i - 1))*0.5d0
         End do
         Results(1) = 2.d0*Raw(2) - Raw(3)
         Results(No) = 2.d0*Raw(No - 1) - Raw(No - 2)
      Else
         WRite (*, *) "Error Interpolation!!"
         Stop
      End If
      Return
   ENd subroutine Interpolation1DX

   subroutine Grad1DX(GradX, GradY, GradZ, Scalar, dx)
      implicit none
      Real(8), intent(in)  :: Scalar(:), dx
      Real(8), intent(out)  :: GradX(:), GradY(:), GradZ(:)
      Integer(4) :: i, Ni, No
      Real(8) ::invdx
      Invdx = 1.d0/dx

      Ni = Size(Scalar)
      No = Size(GradX)
      GradY = 0.d0
      GradZ = 0.d0

      If (No == (Ni - 1)) Then   !!!!!!!!! This part is for Normal condition without extraplotation
         do i = 1, No
            GradX(i) = (Scalar(i + 1) - Scalar(i))*invdx
         End do

      Else If (No == Ni) Then  !!!!!! This part is for Periodical condition
         do i = 1, No - 1
            GradX(i) = (Scalar(i + 1) - Scalar(i))*invdx
         End do
         GradX(No) = GradX(1)

      Else If (No == Ni + 1) Then !!!!! This part is for Normal condition with extraplotation
         do i = 2, No - 1

            GradX(i) = (Scalar(i) - Scalar(i - 1))*invdx
         End do
         GradX(1) = 2.d0*GradX(2) - GradX(3)
         GradX(No) = 2.d0*GradX(No - 1) - GradX(No - 2)
      Else
         WRite (*, *) "Error grad!!"
         Stop
      End If
      Return
   ENd subroutine Grad1DX

   subroutine Laplacian1DX(LapS, Scalar, dx, Periodical)
      implicit none
      Real(8), intent(in)  :: Scalar(:), dx
      Real(8), intent(out)  :: LapS(:)
      Integer(4), Optional :: Periodical
      Integer(4) :: i, Ni, No
      Real(8), Allocatable :: GradX(:), GradY(:), GradZ(:)!,DivS(:)
      Ni = Size(Scalar)
      No = Size(LapS)
      If (Present(Periodical)) Then
         Periodical = 1
         Allocate (GradX(Ni))
         Allocate (GradY(Ni))
         Allocate (GradZ(Ni))
         Call Grad1DX(GradX, GradY, GradZ, Scalar, dx)
         Call Div1DX(LapS, GradX, GradY, GradZ, dx)
      Else
         If (No == (Ni - 2)) Then   !!!!!!!!! This part is for Normal condition without extraplotation
            Allocate (GradX(Ni - 1))
            Allocate (GradY(Ni - 1))
            Allocate (GradZ(Ni - 1))
            Call Grad1DX(GradX, GradY, GradZ, Scalar, dx)
            Call Div1DX(LapS, GradX, GradY, GradZ, dx)
         Else If (No == Ni) Then
            Allocate (GradX(Ni - 1))
            Allocate (GradY(Ni - 1))
            Allocate (GradZ(Ni - 1))
            Call Grad1DX(GradX, GradY, GradZ, Scalar, dx)
            Call Div1DX(LapS, GradX, GradY, GradZ, dx)
         Else
            WRite (*, *) "Error Laplacian1DX!!"
            Stop
         End If
      ENd If
   ENd subroutine Laplacian1DX

   !
   subroutine LaplacianPoisson1DX(LapS, Scalar, dx, Periodical)
      implicit none
      Real(8), intent(in)  :: Scalar(:), dx
      Real(8), intent(out)  :: LapS(:)
      Integer(4), Optional :: Periodical
      Integer(4) :: i, Ni, No
      Real(8) ::invdx
      invdx = 1.d0/dx
      Ni = Size(Scalar)
      No = Size(LapS)
      If (Present(Periodical)) Then
         Periodical = 1
         do i = 2, Ni - 1
            LapS(i) = (Scalar(i + 1) - 2.d0*Scalar(i) + Scalar(i - 1))*invdx*invdx
         End do
         LapS(1) = (Scalar(Ni - 1) - 2.d0*Scalar(1) + Scalar(2))*invdx*invdx
         LapS(Ni) = LapS(1)
      Else
         If (No == (Ni - 2)) Then   !!!!!!!!! This part is for Normal condition without extraplotation
            do i = 2, Ni - 1
               LapS(i - 1) = (Scalar(i + 1) - 2.d0*Scalar(i) + Scalar(i - 1))*invdx*invdx
            End do
         Else If (No == Ni) Then
            do i = 2, Ni - 1
               LapS(i) = (Scalar(i + 1) - 2.d0*Scalar(i) + Scalar(i - 1))*invdx*invdx
            End do
            LapS(1) = 2.d0*LapS(2) - LapS(3)
            LapS(Ni) = 2.d0*LapS(Ni - 1) - LapS(Ni - 2)
         Else
            WRite (*, *) "Error Laplacian1DX!!"
            Stop
         End If
      End IF
   End subroutine LaplacianPoisson1DX

   subroutine Div1DX(divV, VectorX, VectorY, VectorZ, dx)
      implicit none
      Real(8), intent(in)  :: VectorX(:), VectorY(:), VectorZ(:), dx
      Real(8), intent(out)  :: divV(:)
      Integer(4) :: i, Ni, No
      Real(8) ::invdx
      invdx = 1.d0/dx

      Ni = Size(VectorX)
      No = Size(divV)

      If (No == (Ni - 1)) Then   !!!!!!!!! This part is for Normal condition without extraplotation
         do i = 1, No
           divV(i) = (VectorX(i + 1) - VectorX(i))*invdx + (VectorY(i + 1) - VectorY(i))*invdx + (VectorZ(i + 1) - VectorZ(i))*invdx
         End do

      Else If (No == Ni) Then  !!!!!! This part is for Periodical condition
         do i = 1, No - 1
           divV(i) = (VectorX(i + 1) - VectorX(i))*invdx + (VectorY(i + 1) - VectorY(i))*invdx + (VectorZ(i + 1) - VectorZ(i))*invdx
         End do
         divV(No) = divV(1)

      Else If (No == Ni + 1) Then !!!!! This part is for Normal condition with extraplotation
         do i = 2, No - 1
           divV(i) = (VectorX(i) - VectorX(i - 1))*invdx + (VectorY(i) - VectorY(i - 1))*invdx + (VectorZ(i) - VectorZ(i - 1))*invdx
         End do
         divV(1) = 2.d0*divV(2) - divV(3)
         divV(No) = 2.d0*divV(No - 1) - divV(No - 2)
      Else
         WRite (*, *) "Error Div!!"
         Stop
      End If
      Return
   ENd subroutine Div1DX

   subroutine Curl1DX(CurlX, CurlY, CurlZ, VectorX, VectorY, VectorZ, dx)
      implicit none
      Real(8), intent(in)  :: VectorX(:), VectorY(:), VectorZ(:), dx
      Real(8), intent(out)  :: CurlX(:), CurlY(:), CurlZ(:)
      Integer(4) :: i, Ni, No
      Real(8) ::invdx
      invdx = 1.d0/dx

      Ni = Size(VectorX)
      No = Size(CurlX)

      If (No == (Ni - 1)) Then   !!!!!!!!! This part is for Normal condition without extraplotation
         do i = 1, No
            CurlX(i) = 0.d0
            CurlY(i) = -1.d0*(VectorZ(i + 1) - VectorZ(i))*invdx
            CurlZ(i) = (VectorY(i + 1) - VectorY(i))*invdx
         End do

      Else If (No == Ni) Then  !!!!!! This part is for Periodical condition
         do i = 1, No - 1
            CurlX(i) = 0.d0
            CurlY(i) = -1.d0*(VectorZ(i + 1) - VectorZ(i))*invdx
            CurlZ(i) = (VectorY(i + 1) - VectorY(i))*invdx
         End do
         CurlX(No) = CurlX(1)
         CurlY(No) = CurlY(1)
         CurlZ(No) = CurlZ(1)

      Else If (No == Ni + 1) Then !!!!! This part is for Normal condition with extraplotation
         do i = 2, No - 1
            CurlX(i) = 0.d0
            CurlY(i) = -1.d0*(VectorZ(i) - VectorZ(i - 1))*invdx
            CurlZ(i) = (VectorY(i) - VectorY(i - 1))*invdx
         End do
         CurlX(1) = 2.d0*CurlX(2) - CurlX(3)
         CurlY(1) = 2.d0*CurlY(2) - CurlY(3)
         CurlZ(1) = 2.d0*CurlZ(2) - CurlZ(3)
         CurlX(No) = 2.d0*CurlX(No - 1) - CurlX(No - 2)
         CurlY(No) = 2.d0*CurlY(No - 1) - CurlY(No - 2)
         CurlZ(No) = 2.d0*CurlZ(No - 1) - CurlZ(No - 2)
      Else
         WRite (*, *) "Error Curl!!"
         Stop
      End If
      Return
   ENd subroutine Curl1DX

   subroutine DivSymmTensorN2C(VectorX, VectorY, VectorZ, TensorXX, TensorYY, TensorZZ, TensorXY, TensorYZ, TensorXZ, dx)
      !subroutine DivSymmTensorN2C(Nx,TensorXX,TensorXY,TensorXZ,VectorX,dx)
      implicit none
      Real(8), intent(in)  :: TensorXX(:), TensorYY(:), TensorZZ(:), TensorXY(:), TensorYZ(:), TensorXZ(:), dx
      Real(8), intent(out)  :: VectorX(:), VectorY(:), VectorZ(:)
      Integer(4) :: i, Ni, No
      Real(8) ::invdx
      invdx = 1.d0/dx
      Ni = Size(TensorXX)
      No = Size(VectorX)

      If (No == (Ni - 1)) Then   !!!!!!!!! This part is for Normal condition without extraplotation
         do i = 1, No
  VectorX(i) = (TensorXX(i + 1) - TensorXX(i))*invdx + (TensorXY(i + 1) - TensorXY(i))*invdx + (TensorXZ(i + 1) - TensorXZ(i))*invdx
  VectorY(i) = (TensorXY(i + 1) - TensorXY(i))*invdx + (TensorYY(i + 1) - TensorYY(i))*invdx + (TensorYZ(i + 1) - TensorYZ(i))*invdx
  VectorZ(i) = (TensorXZ(i + 1) - TensorXZ(i))*invdx + (TensorYZ(i + 1) - TensorYZ(i))*invdx + (TensorZZ(i + 1) - TensorXZ(i))*invdx
         End do

      Else If (No == Ni) Then  !!!!!! This part is for Periodical condition
         do i = 1, No - 1
  VectorX(i) = (TensorXX(i + 1) - TensorXX(i))*invdx + (TensorXY(i + 1) - TensorXY(i))*invdx + (TensorXZ(i + 1) - TensorXZ(i))*invdx
  VectorY(i) = (TensorXY(i + 1) - TensorXY(i))*invdx + (TensorYY(i + 1) - TensorYY(i))*invdx + (TensorYZ(i + 1) - TensorYZ(i))*invdx
  VectorZ(i) = (TensorXZ(i + 1) - TensorXZ(i))*invdx + (TensorYZ(i + 1) - TensorYZ(i))*invdx + (TensorZZ(i + 1) - TensorXZ(i))*invdx
         End do
         VectorX(No) = VectorX(1)
         VectorY(No) = VectorY(1)
         VectorZ(No) = VectorZ(1)

      Else If (No == Ni + 1) Then !!!!! This part is for Normal condition with extraplotation
         do i = 2, No - 1
  VectorX(i) = (TensorXX(i) - TensorXX(i - 1))*invdx + (TensorXY(i) - TensorXY(i - 1))*invdx + (TensorXZ(i) - TensorXZ(i - 1))*invdx
  VectorY(i) = (TensorXY(i) - TensorXY(i - 1))*invdx + (TensorYY(i) - TensorYY(i - 1))*invdx + (TensorYZ(i) - TensorYZ(i - 1))*invdx
  VectorZ(i) = (TensorXZ(i) - TensorXZ(i - 1))*invdx + (TensorYZ(i) - TensorYZ(i - 1))*invdx + (TensorZZ(i) - TensorXZ(i - 1))*invdx
         End do
         VectorX(1) = 2.d0*VectorX(2) - VectorX(3)
         VectorY(1) = 2.d0*VectorY(2) - VectorY(3)
         VectorZ(1) = 2.d0*VectorZ(2) - VectorZ(3)
         VectorX(No) = 2.d0*VectorX(No - 1) - VectorX(No - 2)
         VectorY(No) = 2.d0*VectorY(No - 1) - VectorY(No - 2)
         VectorZ(No) = 2.d0*VectorZ(No - 1) - VectorZ(No - 2)
      Else
         WRite (*, *) "Error Curl!!"
         Stop
      End If
      Return
   End subroutine DivSymmTensorN2C

   !

End Module ModuleVector1DX
