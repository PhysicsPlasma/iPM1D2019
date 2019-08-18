Module MCCEnergyKai
     Use Constants
     Implicit none
     contains
     Function IsotropicCosKai(Energy) 
       Implicit none
       Real(8) :: IsotropicCosKai,Energy
       CALL RANDOM_NUMBER(R)
       IsotropicCosKai=1.d0-2.d0*R
       return
    end Function IsotropicCosKai
    
     Function CosKaiVahedi(energy)
       !!!!!Warinig: ArcosTheta cant not be equal to +-1.0.
	    implicit none
    	real(8) :: energy
        real(8) :: CosKaiVahedi
        if(Energy < Minreal) then
		     CosKaiVahedi= 1.d0
          else
		     CALL RANDOM_NUMBER(R)
		     CosKaiVahedi= (energy +2.0d0 -2.0d0*(energy+1.0d0)**R)/energy
		end if
 		return 
     end function CosKaiVahedi
     
    Function CosKaiOkhrimovskyy(Chi)
       !!!!!Warinig: ArcosTheta cant not be equal to +-1.0.
	    implicit none
    	real(8) :: Chi
        real(8) :: CosKaiOkhrimovskyy
        CALL RANDOM_NUMBER(R)  
        CosKaiOkhrimovskyy= 1.d0-2.0d0*R*(1-Chi)/(1+Chi*(1-2.d0*R))
 		return 
     end function CosKaiOkhrimovskyy
   
    Function IsotropicEnergy(Energy)
        implicit none
    	Real(8) :: Energy
        Real(8) :: IsotropicEnergy
        CALL RANDOM_NUMBER(R)
        IsotropicEnergy=R*Energy
        If(abs(IsotropicEnergy)<MinReal) then
                   IsotropicEnergy=MinReal
         end  if   
        return
    end Function  IsotropicEnergy
   
     Function CreatEnergyVehadi(Energy)
        implicit none
    	Real(8) :: Energy
        Real(8) :: CreatEnergyVehadi
        CALL RANDOM_NUMBER(R)
        CreatEnergyVehadi=10.d0*Dtan(R*DAtan(Energy)/20.d0)
        If(abs(CreatEnergyVehadi)<MinReal) then
                   CreatEnergyVehadi=MinReal
         end  if   
        return
    end Function  CreatEnergyVehadi
End Module  MCCEnergyKai