Module Numrical
     Implicit None
     contains
     subroutine Interpolation1D(Ni,Input,No,Output)
        implicit none
        Integer(4),Intent(in) :: Ni,No
        Real(8),Intent(in) :: Input(Ni)
        Real(8),intent(inout) :: Output(No)  
        Integer(4) :: i,NTemp
        Real(8) :: dx1,dx2,S1,S2,XTemp
        dx1=1.d0/dble(Ni-1)
        dx2=1.d0/dble(No-1)
        do i=1,No
              Xtemp=dble(i-1)/dble(No-1)
              NTemp=Ceiling(XTemp/dx1)
             If(NTemp==0) then
                  Output(i)=(Xtemp+0.5d0)*Input(1)+(0.5-Xtemp)*Input(Ni)
             Else 
                 S1=XTemp-dble(NTemp)*dx1
                 S2=1.d0-S1 
                 Output(i)=S1*Input(NTemp)+S2*Input(NTemp+1)
             End if 
         end do
         return
     end subroutine Interpolation1D

  Subroutine AnyAvg( rN , rM)
    Real (8), Intent( IN ) :: rN(:)
    Real (8), Intent( OUT ) :: rM(:)
    integer :: n , m , i , j , e , k
    real :: rs , re , rl , rtmp , add
    n = size( rN )
    m = size( rM )
    if ( n < m ) Return
    if ( n == m ) then
      rM = rN
      return
    end if
    If ( mod( n , m ) == 0 ) then !// 常规整个数平緛E
      e = n / m
      Do i = 1 , m
        rs = 0.0
        Do j = (i-1)*e+1 , i*e
          rs = rs + rN(j)/e
        End Do
        rM(i) = rs
      End Do
    Else !// 小数个数平緛E
      re = n*1.0 / m*1.0
      rl = 0.0
      k = 1
      Do i = 1 , m
        rs = 0.0
        add = re
        rtmp = rl - int(rl)
        if ( rtmp > 0.0001 ) then !// 如果上次有剩觼E
          rs = rs + rN(k)*(1.0-rtmp)/re
          add = add - (1.0-rtmp)
          k = k + 1
        end if
        Do while ( add >= 1.0 )
          rs = rs + rN(k)/re
          k = k + 1
          add = add - 1.0
        End Do
        if ( add > 0.0001 ) then !//如果本次有剩觼E
          rs = rs + rN(k)*add / re
        end if
        rM(i) = rs
        rl = rl + re
      End Do
    End If
  End Subroutine AnyAvg
  
  SUBROUTINE tridag(a,b,c,r,u,n)
    integer,PARAMETER ::nmax=10000
    Integer,intent(in) :: n
    real(8), intent(in) :: a(n),b(n),c(n),r(n)
    real(8), intent(inout) :: u(n)
    INTEGER j
    REAL(8) gam(nmax)
    real(8) bet
        !A=1.d0
        !B=-2.d0
        !C=1.d0
        bet=b(1)
        u(1)=r(1)/bet
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        u(j)=(r(j)-a(j)*u(j-1))/bet
    end do
    do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    end do
    return
   END SUBROUTINE tridag

     
            SUBROUTINE locate(xx,n,x,j)
                    INTEGER(4) j,n
                    REAL(8) aaa,x,xx(n)
                    INTEGER(4) jl,jm,ju
                    jl=0
                    ju=n+1
                    do
                      aaa=0.
                      if(ju-jl>1) then
                        aaa=2.
                        jm=(ju+jl)/2
                        if((xx(n)>=xx(1)).eqv.(x>=xx(jm))) then
                          jl=jm
                        else
                          ju=jm
                        endif
                        if(.not.aaa>1.) exit
                      endif
                      if(.not.aaa>1.) exit
                    end do
                    if(x==xx(1)) then
                      j=1
                    else if(x==xx(n)) then
                      j=n-1
                    else
                      j=jl
                    endif
            END SUBROUTINE locate
           
          
          
            ! nr.f90
            ! Some subroutines come from Numerical Recipes in FORTRAN 
            ! C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
            ! Only tiny tune for performance (Wang Hongyu,2007)
            ! Subroutines: 
            ! FOUR -----calculate the Fast Fourier Transform 
            ! REALFT -------calculate the real Fast Fourier Transform 
            ! SINFT ------calculate the sine Fast Fourier Transform 
            ! TRIDIAG ----- solve the tridiagonal linear equations

            ! subroutine realft
            ! Calculates the Fourier transform of a set of n real-valued data points. Replaces this data
            ! (which is stored in array data(1:n)) by the positive frequency half of its complex Fourier
            ! transform. The real-valued rst and last components of the complex transform are returned
            ! as elements data(1) and data(2), respectively. n must be a power of 2. This routine
            ! also calculates the inverse transform of a complex data array if it is the transform of real
            ! data. (Result in this case must be multiplied by 2/n.)
            ! subroutine sinft
            ! Calculates the sine transform of a set of n real-valued data points stored in array y(1:n).
            ! The number n must be a power of 2. On exit y is replaced by its transform. This program,
            ! without changes, also calculates the inverse sine transform, but in this case the output array
            ! should be multiplied by 2/n. 
             
                  SUBROUTINE sinft(y,n)
                  implicit none
                  INTEGER n
                  double precision y(n)
                  INTEGER j
                  double precision sum,y1,y2
                  DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
                  integer:: sinftinited=0,localn=0
                  save sinftinited,localn,wpr,wpi
                  theta=3.1415926535897932385d0/dble(n)
                  if (sinftinited==0 .or. localn/=n) then
                  sinftinited=1
                  localn=n
                  wpr=-2.0d0*dsin(0.5d0*theta)**2
                  wpi=dsin(theta)
                  endif
                  wr=1.0d0
                  wi=0.0d0
                  y(1)=0.0
                  do j=1,n/2
                    wtemp=wr
                    wr=wr*wpr-wi*wpi+wr
                    wi=wi*wpr+wtemp*wpi+wi
                    y1=wi*(y(j+1)+y(n-j+1))
                    y2=0.5d0*(y(j+1)-y(n-j+1))
                    y(j+1)=y1+y2
                    y(n-j+1)=y1-y2
                 enddo
                  call realft(y,n,+1)
                  sum=0.0
                  y(1)=0.5d0*y(1)
                  y(2)=0.0
                  do j=1,n-1,2
                    sum=sum+y(j)
                    y(j)=y(j+1)
                    y(j+1)=sum
                 enddo
                  return
                  END SUBROUTINE

               SUBROUTINE realft(data,n,isign)
                   implicit none
                  INTEGER isign,n
                  double precision data(n)
                  INTEGER i,i1,i2,i3,i4,n2p3
                  double precision c1,c2,h1i,h1r,h2i,h2r,wis,wrs
                  DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
                  integer :: rftinited=0,rftln=0
                  save rftinited,rftln,wpi,wpr
                  theta=3.1415926535897932385d0/dble(n/2)
                  c1=0.5
                  if (isign.eq.1) then
                    c2=-0.5
                    call four1(data,n/2,+1)
                  else
                    c2=0.5
                    theta=-theta
                  endif
                  if(rftinited==0 .or.rftln/=n) then
                  wpr=-2.0d0*dsin(0.5d0*theta)**2
                  wpi=dsin(theta)
                  endif
                  wr=1.0d0+wpr
                  wi=wpi
                  n2p3=n+3
                  do i=2,n/4
                    i1=2*i-1
                    i2=i1+1
                    i3=n2p3-i2
                    i4=i3+1
                    wrs=dble(wr)
                    wis=dble(wi)
                    h1r=c1*(data(i1)+data(i3))
                    h1i=c1*(data(i2)-data(i4))
                    h2r=-c2*(data(i2)+data(i4))
                    h2i=c2*(data(i1)-data(i3))
                    data(i1)=h1r+wrs*h2r-wis*h2i
                    data(i2)=h1i+wrs*h2i+wis*h2r
                    data(i3)=h1r-wrs*h2r+wis*h2i
                    data(i4)=-h1i+wrs*h2i+wis*h2r
                    wtemp=wr
                    wr=wr*wpr-wi*wpi+wr
                    wi=wi*wpr+wtemp*wpi+wi
                 enddo
                  if (isign.eq.1) then
                    h1r=data(1)
                    data(1)=h1r+data(2)
                    data(2)=h1r-data(2)
                  else
                    h1r=data(1)
                    data(1)=c1*(h1r+data(2))
                    data(2)=c1*(h1r-data(2))
                    call four1(data,n/2,-1)
                  endif
                  return
                  END SUBROUTINE
            ! subroutine four1
            ! Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces
            ! data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as ?1.
            ! data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn
            !  MUST be an integer power of 2 (this is not checked).

            SUBROUTINE four1(data1,nn,isign)
            double precision  wsinpr(40),wsinpi(40)
            INTEGER isign,nn
            double precision data1(2*nn)

            INTEGER i,istep,j,m,mmax,n,ilog
            double precision tempi,tempr
            DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
            integer,save :: inited=0,localnn=0
            save wsinpr,wsinpi
            n=2*nn
            j=1
            if(inited==0 .or. localnn/=nn) then
            inited=1
            localnn=nn
            mmax=2
            ilog=1
            do while(n>mmax)
                istep=2*mmax
                theta=6.28318530717959d0/(isign*mmax)
                wsinpr(ilog)=-2.d0*dsin(0.5d0*theta)**2
                wsinpi(ilog)=dsin(theta)
                ilog=ilog+1
                mmax=istep
            enddo    
            endif        

            do i=1,n,2
              if(j>i)then
                tempr=data1(j)
                tempi=data1(j+1)
                data1(j)=data1(i)
                data1(j+1)=data1(i+1)
                data1(i)=tempr
                data1(i+1)=tempi
              endif
              m=n/2
              do while((m>=2).and.(j>m)) 
                j=j-m
                m=m/2
              end do
              j=j+m
            end do
            mmax=2
            ilog=1
            do while(n>mmax) 
              istep=2*mmax
              theta=6.28318530717959d0/(isign*mmax)
              
              wpr=wsinpr(ilog)
              wpi=wsinpi(ilog)
              wr=1.d0
              wi=0.d0
              do m=1,mmax,2
                do i=m,n,istep
                  j=i+mmax
                  tempr=dble(wr)*data1(j)-dble(wi)*data1(j+1)
                  tempi=dble(wr)*data1(j+1)+dble(wi)*data1(j)
                  data1(j)=data1(i)-tempr
                  data1(j+1)=data1(i+1)-tempi
                  data1(i)=data1(i)+tempr
                  data1(i+1)=data1(i+1)+tempi
                end do
                wtemp=wr
                wr=wr*wpr-wi*wpi+wr
                wi=wi*wpr+wtemp*wpi+wi
              end do
              mmax=istep
              ilog=ilog+1
            end do
            END SUBROUTINE four1

            ! subroutine tridiag
            ! Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
            ! a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modied.
            ! Parameter: NMAX is the maximum expected value of n.
         subroutine mytdginit(a,b,c,gam,bet,n)
            implicit none
            integer n
            double precision a(n),b(n),c(n),gam(n)
            double precision bet(n)

            integer j
                bet(1)=1.d0/b(1)
                  do j=2,n
                    gam(j)=c(j-1)*bet(j-1)
                    bet(j)=1.d0/(b(j)-a(j)*gam(j))
                     enddo
            end subroutine mytdginit
            
            subroutine mytdrun(a,b,c,gam,bet,r,u,n)
            implicit none
            integer n
            double precision a(n),b(n),c(n),r(n),u(n)
            double precision gam(n),bet(n)
            integer j
               u(1)=r(1)*bet(1)
              do j=2,n
                    u(j)=(r(j)-a(j)*u(j-1))*bet(j)
                 enddo
              
            do j=n-1,1,-1
                    u(j)=u(j)-gam(j+1)*u(j+1)
                 enddo   
            end subroutine mytdrun
            
            SUBROUTINE tridagopt(a,b,c,r,u,n)

            implicit none
            INTEGER n,NMAX
            double precision a(n),b(n),c(n),r(n),u(n)
            PARAMETER (NMAX=10000)
            INTEGER j
            double precision bet(n),gam(NMAX)
                bet(1)=1.d0/b(1)
                u(1)=r(1)*bet(1)
                  do j=2,n
                    gam(j)=c(j-1)*bet(j-1)
                    bet(j)=1.d0/(b(j)-a(j)*gam(j))
                    u(j)=(r(j)-a(j)*u(j-1))*bet(j)
                 enddo
                  do j=n-1,1,-1
                    u(j)=u(j)-gam(j+1)*u(j+1)
                 enddo   
                  return
                  END SUBROUTINE tridagopt
 End Module Numrical




        
