Module ModuleDiagOneStep
      Use ModuleControlFlow
      Use DiagnosticsCollisionRate
      Use DiagnosticsEEPF
      Use DiagnosticsMomentum
      Use ModuleOneStepField
      Use ModuleOneStep
      Use DiagnosticsTestParticle
      
      !Use DiagnosticsTestParticle
      Implicit none
      Type(Grid1D(Nx=NxMax,Ns=11)),save :: G1DDiagParticleField
      Type(Grid2D(Nx=NxMax,Ny=100,Ns=11)),save :: G2DDiagParticleField
      
      Type(Grid1D(Nx=NxMax,Ns=3)),save :: G1DDiagParticleCR
      Type(Grid2D(Nx=NxMax,Ny=100,Ns=3)),save :: G2DDiagParticleCR

      Type(Grid1D(Nx=500,Ns=2)),save :: G1DDiagParticleEDF
      Type(Grid2D(Nx=500,Ny=100,Ns=2)),save :: G2DDiagParticleEDF
      
      !Type(ParticleElectrode),save :: ParticleElectrodeGlobal(0:NsMax)
      
      contains
      Subroutine DiagInitilalization(CF)
            Implicit none
            Class(ControlFlow), intent(in) :: CF
            !Integer(4),intent(in) ::  Ns
            !Type(ParticleBundle),intent(in) :: PB(0:Ns)
                 Call G1DDiagParticleField%Init(CF)
                 Call G2DDiagParticleField%Init(CF)
                 Call G1DDiagParticleCR%Init(CF)
                 Call G2DDiagParticleCR%Init(CF)
                 Call G1DDiagParticleEDF%Init(CF)
                 Call G2DDiagParticleEDF%Init(CF)
                 
          return  
      End Subroutine DiagInitilalization
      
     Subroutine DiagInitilalization2()
            Implicit none
            !Class(ControlFlow), intent(in) :: CF
            !Integer(4),intent(in) ::  Ns
            !Type(ParticleBundle),intent(in) :: PB(0:Ns)
            Call DiagParticleTestParticle(ParticleGlobal(0),ParticleBDOneGlobal(0),FieldGlobal,-1)              
          return  
      End Subroutine DiagInitilalization2
      
      
      Subroutine  DiagOneStep()
          Implicit none   
             Call DiagParticleFieldPeriod(G1DDiagParticleField,1,ParticleGlobal,FieldGlobal,0)
             Call DiagParticleFieldPeriod(G2DDiagParticleField,1,ParticleGlobal,FieldGlobal,0)
             Call DiagParticleEDFOne(G1DDiagParticleEDF,ParticleGlobal(0),0)
             Call DiagParticleEDFOne(G2DDiagParticleEDF,ParticleGlobal(0),0)
             Call DiagParticleCollisionRateOne(G1DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),0)
             Call DiagParticleCollisionRateOne(G2DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),0)
        return  
      End Subroutine DiagOneStep
      
      Subroutine  DiagOneStep2()
          Implicit none   
           Call DiagParticleTestParticle(ParticleGlobal(0),ParticleBDOneGlobal(0),FieldGlobal,0)              
           return  
      End Subroutine DiagOneStep2
       !
       Subroutine DiagOneStepFinal()
       Implicit none  
             Call DiagParticleFieldPeriod(G1DDiagParticleField,1,ParticleGlobal,FieldGlobal,1)
             Call DiagParticleFieldPeriod(G2DDiagParticleField,1,ParticleGlobal,FieldGlobal,1)
             Call DiagParticleEDFOne(G1DDiagParticleEDF,ParticleGlobal(0),1)
             Call DiagParticleEDFOne(G2DDiagParticleEDF,ParticleGlobal(0),1)
             Call DiagParticleCollisionRateOne(G1DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),1)
             Call DiagParticleCollisionRateOne(G2DDiagParticleCR,ParticleGlobal(0),MCCBundleGlobal(0,1),1)
          return  
       End Subroutine DiagOneStepFinal
       
      Subroutine DiagOneStepFinal2()
       Implicit none  
          Call DiagParticleTestParticle(ParticleGlobal(0),ParticleBDOneGlobal(0),FieldGlobal,1)              
          return  
        End Subroutine DiagOneStepFinal2

  End Module ModuleDiagOneStep