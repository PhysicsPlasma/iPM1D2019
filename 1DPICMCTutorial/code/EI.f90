module ModuleEI
    use ModuleParticleBundle
    use ModuleParticleBoundary
    use Constants
    implicit none
    
contains
    subroutine EI(PB,dx,dt,Xmin)
    implicit none
        type(ParticleBundle),INTENT(INOUT)::PB
        type(ParticleOne) :: PTemp
        REAL(8),INTENT(IN)::dx,dt,Xmin
        INTEGER(4)::i,Min=0
        real(8):: Qeb=0,Veb=0,XFactor,VFactor,Residue=0,XRandom=0.01,Neb=0
        real(8):: Eeb=20.d0,Jeb=0.1d0,Reb=0.1d0
        Veb=dsqrt(2.0*Eeb*JtoeV/ElectronMass)
        Neb=Jeb/(ElectronCharge*Veb*Reb**2*PI)/PB%Weight

        Min=int(Neb)
        Residue=Neb-dble(Min)
        call random_number(R)
        if ( R<Residue ) Min=Min+1
        !write(*,*) '+Ne:',Min

        VFactor=1.d0/PB%VFactor
        XFactor=1.d0/PB%XFactor
        do i = 1, Min
            call random_number(R)
            ! PTemp%X=Xmin+R*XRandom
            PTemp%X=Xmin/XFactor+Veb*dt*R
            CALL PTemp%PosRes(XFactor)
            call PTemp%VelInpInit(Veb,0.d0,0.d0)
            call PTemp%AccInpInit()
            
            call PTemp%VelRes(VFactor)
            PTemp%Vx=abs(PTemp%Vx)
            call PB%AddOne(PTemp)
        end do
        
    end subroutine EI
end module ModuleEI