  Í\  ä   k820309              19.1        äí_                                                                                                          
       D:\Source\iPM1D2019\1DPICMCTutorial\code\mc\MCCElectron.F90 MODULEMCCELECTRON                                                     
                                                           
                         @               A                'è                    #REACTIONINDEX    #PARTICLEANNIHILATION    #PARTICLECREATION    #POI    #MIU i   #MASSRATIO j   #ENERGY k   #BETA l   #GX m   #GY n   #GZ o   #GPER p   #G q   #POT r   #NPONEW s   #PON t   #SELECT    #UPDATER    #VELOCITYUPDATER                                                                                                                                              0                                                                                                                                                                                                                                                                                                                                                     8                    #PARTICLEONE                                          y#PARTICLEONE                                                                     @                                '8                    #X 	   #VX 
   #VY    #VZ    #AX    #AY    #AZ    #POSINIT    #VELINPINIT    #VELMAXINIT    #VELRANINIT $   #ACCINPINIT (   #POSRES .   #VELRES 2   #ACCRES 6   #ENERGY :   #COPY ?   #SWAP C   #WEIGHTP2C G   #WEIGHTC2PES L   #MOVEES R   #WEIGHTC2PEM V   #MOVEEM `                                              	                
                                              
               
                                                             
                                                             
                                                              
                                                   (          
                                                   0          
   1         À                                                  #POSITIONRANDOMINITIALIZATIONPARTICLEONE    #         @                                                      #PO    #XL    #XU    #YL    #YU    #ZL    #ZU              
                                     8               #PARTICLEONE              
                                     
                
                                     
                
                                     
                
                                     
                
                                     
                
                                     
      1         À                                             	     #VELOCITYINPUTINITIALIZATIONPARTICLEONE    #         @                                                       #PO    #VX    #VY    #VZ              
                                     8               #PARTICLEONE              
                                     
                
                                     
                
                                     
      1         À                                             
     #VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE     #         @                                                        #PO !   #MASS "   #TEMPERATURE #             
                                !     8               #PARTICLEONE              
                                 "     
                
                                 #     
      1         À                                $                  #VELOCITYRANDOMINITIALIZATIONPARTICLEONE %   #         @                                   %                    #PO &   #V '             
                                &     8               #PARTICLEONE              
                                 '     
      1         À                                (                  #ACCELERATIONINPUTINITIALIZATIONPARTICLEONE )   #         @                                   )                    #PO *   #AX +   #AY ,   #AZ -             
                                *     8               #PARTICLEONE              
                                +     
                
                                ,     
                
                                -     
      1         À                                .                  #POSITIONRESCALEPARTICLEONE /   #         @                                   /                    #PO 0   #XFACTOR 1             
                                0     8               #PARTICLEONE              
                                 1     
      1         À                                2                  #VELOCITYRESCALEPARTICLEONE 3   #         @                                   3                    #PO 4   #VFACTOR 5             
                                4     8               #PARTICLEONE              
                                 5     
      1         À                                6                  #ACCELERATIONRESCALEPARTICLEONE 7   #         @                                   7                    #PO 8   #AFACTOR 9             
                                8     8               #PARTICLEONE              
                                 9     
      1         À                               :              	    #ENERGYPARTICLEONE ;   %         @                                ;                    
       #PO <   #MASS =   #VFACTOR >             
                                 <     8              #PARTICLEONE              
                                 =     
                
                                >     
      1         À                                ?              
    #COPYPARTICLEONE @   #         @                                   @                    #POD A   #POC B             
                                A     8               #PARTICLEONE              
                                 B     8              #PARTICLEONE    1         À                                C                  #SWAPPARTICLEONE D   #         @                                   D                    #POD E   #POC F             
                                E     8               #PARTICLEONE              
                                F     8               #PARTICLEONE    1         À                                G                  #WEIGHTP2CPARTICLEONE H   #         @                                   H                    #PO I   #FO J             
                                 I     8              #PARTICLEONE              
                                 J     (              #FIELDONE K   1         À                                L                  #WEIGHTC2PELECTROSTATICPARTICLEONE M   #         @                                   M                    #PO N   #FG O   #EX Q             
                                N     8               #PARTICLEONE              
                                  O     `Z             #FIELD P                                             Q     
       1         À                                R                  #MOVEELECTROSTATICPARTICLEONE S   #         @                                   S                    #PO T   #EX U             
                                T     8               #PARTICLEONE              
                                U     
      1         À                                V                  #WEIGHTC2PELECTROMAGNETICPARTICLEONE W   #         @                                   W                    #PO X   #FG Y   #EX Z   #EY [   #EZ \   #BX ]   #BY ^   #BZ _             
                                X     8               #PARTICLEONE              
                                  Y     `Z             #FIELD P                                             Z     
                                                 [     
                                                 \     
                                                 ]     
                                                 ^     
                                                 _     
       1         À                                `                  #MOVEELECTROMAGNETICPARTICLEONE a   #         @                                   a                    #PO b   #EX c   #EY d   #EZ e   #BX f   #BY g   #BZ h             
                                b     8               #PARTICLEONE              
                                 c     
                
                                 d     
                
                                 e     
                
                                 f     
                
                                 g     
                
                                 h     
                                                 i               
                                              j                
                                              k     (          
                                             l     0         
                                                
                                 0.D0                                               m     8       	   
                                              n     @       
   
                                              o     H          
                                              p     P          
                                              q     X          
                                               r     8       `              #PARTICLEONE                                              s                                                                                               0                                              t                    @             #PARTICLEONEINDEX u             &                                                             @                          u     '@                    #PARTICLEONE v   #INDEX w   #INDEXINIT x   #COPY |   #SWAP                                                v     8                      #PARTICLEONE                                               w     8             1         À                               x                  #INDEXINITIALIZATIONPARTICLEONEINDEX y   #         @                                  y                    #POI z   #INDEX {             
                                z     @               #PARTICLEONEINDEX u             
                                 {           1         À                               |              
    #COPYPARTICLEONEINDEX }   #         @                                   }                    #POD ~   #POC              
                                ~     @               #PARTICLEONEINDEX u             
                                      8              #PARTICLEONE    1         À                                                 #SWAPPARTICLEONEINDEX    #         @                                                       #POD    #POC              
                                     @               #PARTICLEONEINDEX u             
                                     8               #PARTICLEONE    1         À                                                  #SELECTMCCPARTICLEONE    #         @                                                       #MCPO    #PO              
                                     è               #MCCPARTICLEONE              
                                       8              #PARTICLEONE    1         À                                                  #UPDATEMCCPARTICLEONE    #         @                                                       #MCPO    #SO    #GO              
                                     è               #MCCPARTICLEONE              
                                       °              #SPECYONE              
                                       °              #GASONE    1         À                                                 #UPDATEVELOCITYMCCPARTICLEONE    #         @                                                      #MCPO    #VFACTOR              
                                     è               #MCCPARTICLEONE              
                                      
                        @                               '°                    #NAME    #SPECYINDEX    #GASINDEX    #CHARGE    #MASS    #RADIUS    #NATOM    #INITDENSITY    #DENSITY    #INITTEMPERATURE    #TEMPERATURE                                                    c                                                                     d                                                                             ÿÿÿÿÿÿÿÿ                                                           h                                                                             ÿÿÿÿÿÿÿÿ                                                            p          
                                                   x          
                                                             
                                                             
                                                             
                                                          	   
                                                           
   
                                                   ¨          
                     @                               '°                    #NAME    #MCMODEL    #NS     #INDEXSTART ¡   #MASS ¢   #RADIUS £   #NATOM ¤   #BETAMAX ¥   #INITDENSITY ¦   #DENSITY §   #INITTEMPERATURE ¨   #TEMPERATURE ©                                                   c                                                                     d                                                                             ÿÿÿÿÿÿÿÿ                                                             h                                                        ¡     l                                                        ¢     p          
                                              £     x          
                                              ¤               
                                              ¥               
                                              ¦            	   
                                              §            
   
                                              ¨                
                                              ©     ¨          
                     @                          ª     '                    #REACTANT «   #REACTIONTYPE ¬   #RESULTANT ­   #THRESHOLD ®                                              «                                                              ¬                                                             ­                                                             ®               
                     @                          K     '(                   #NX ¯   #DX °   #DT ±   #RHOONE ²   #CHIONE ³                                             ¯                                                                                 A              321                                              °              
                                                  
                  {®Gáz4?                                                      ±              
                                                 
                 »½×Ùß|Û=        1.D-10                                               ²     A                       
  p          & p        p A          p A                                                                   ³     A       
                
  p          & p        p A          p A                                          @                          P     '`Z                   #NX ´   #DX µ   #DT ¶   #EX ·   #EY ¸   #EZ ¹   #BX º   #BY »   #BZ ¼   #RHO ½   #PHI ¾   #CHI ¿   #DUMP À   #LOAD Ä                                             ´                                                                                 A              321                                              µ              
                                                  
                  {®Gáz4?                                                      ¶              
                                                 
                 »½×Ùß|Û=        1.D-10                                              ·     A                      
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          ¸     A       
               
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          ¹     A      (               
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          º     A      0               
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          »     A      8(               
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          ¼     A      @2             	  
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                           ½     A      H<             
   
  p          & p        p A          p A                                                                   ¾     A      PF                
  p          & p        p A          p A                                                                   ¿     A      XP                
  p          & p        p A          p A                        1         À                                À                  #DUMPFIELD Á   #         @                                   Á                    #FG Â   #MODE Ã             
                                Â     `Z              #FIELD P                                             Ã            1         À                                Ä                  #LOADFIELD Å   #         @                                   Å                    #FG Æ   #STATUS Ç             
                                Æ     `Z              #FIELD P             
                                 Ç                              @               A           È     'ð                    #MODEL É   #NREACTION Ê   #NSIGMA Ë   #ENERGYMIN Ì   #ENERGYINTERVAL Í   #ENERGYMAX Î   #COLLISIONRATIO Ï   #SIGMAMAX Ð   #SO Ñ   #GO Ò   #REACTION Ó   #PROBILITY Ô   #DUMP Õ                                              É                                                             Ê                                                                                               0                                              Ë                                                                                               0                                               Ì               
                                              Í               
                                              Î                
                                              Ï     (          
                                              Ð     0          
                                             Ñ     °       8       	      #SPECYONE                                          y#SPECYONE                                                                                             Ò     °       @       
      #GASONE                                          y#GASONE                                                                                             Ó            H                    #REACTIONONE ª             &                                                                                    Ô                             
            &                   &                                           1         À                                Õ                  #DUMPMCCBUNDLE Ö   #         @                                   Ö                    #MCB ×             
                                ×     ð               #MCCBUNDLE È   #         @                                  Ø                    #MCPO Ù   #COSKAI Ú             
                                 Ù     è               #MCCPARTICLEONE    "                                        Ú       
     %         @  @                            Û                    
       #ENERGY Ü                                             Ü     
       %         @  @                            Ý                    
       #ENERGY Þ                                             Þ     
       #         @                                   ß                    #MCPO à   #SO á   #GO â   #RO ã             
D @                               à     è               #MCCPARTICLEONE              
  @                               á     °              #SPECYONE              
  @                               â     °              #GASONE              
  @                               ã                   #REACTIONONE ª          V      fn#fn    ö   @   J   MODULETYPEMCC    6  @   J   MCCENERGYKAI -   v  ?      MCCPARTICLEONE+MODULETYPEMCC ;   µ  ¥   a   MCCPARTICLEONE%REACTIONINDEX+MODULETYPEMCC B   Z  ¤   a   MCCPARTICLEONE%PARTICLEANNIHILATION+MODULETYPEMCC >   þ  ¤   a   MCCPARTICLEONE%PARTICLECREATION+MODULETYPEMCC 1   ¢  Ò   a   MCCPARTICLEONE%POI+MODULETYPEMCC .   t  a      PARTICLEONE+MODULEPARTICLEONE 0   Õ  H   a   PARTICLEONE%X+MODULEPARTICLEONE 1     H   a   PARTICLEONE%VX+MODULEPARTICLEONE 1   e  H   a   PARTICLEONE%VY+MODULEPARTICLEONE 1   ­  H   a   PARTICLEONE%VZ+MODULEPARTICLEONE 1   õ  H   a   PARTICLEONE%AX+MODULEPARTICLEONE 1   =  H   a   PARTICLEONE%AY+MODULEPARTICLEONE 1     H   a   PARTICLEONE%AZ+MODULEPARTICLEONE 6   Í  u   a   PARTICLEONE%POSINIT+MODULEPARTICLEONE J   B	         POSITIONRANDOMINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE M   Â	  Y   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE M   
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%XL+MODULEPARTICLEONE M   [
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%XU+MODULEPARTICLEONE M   
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%YL+MODULEPARTICLEONE M   Û
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%YU+MODULEPARTICLEONE M     @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%ZL+MODULEPARTICLEONE M   [  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%ZU+MODULEPARTICLEONE 9     t   a   PARTICLEONE%VELINPINIT+MODULEPARTICLEONE I     h       VELOCITYINPUTINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE L   w  Y   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE L   Ð  @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VX+MODULEPARTICLEONE L     @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VY+MODULEPARTICLEONE L   P  @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VZ+MODULEPARTICLEONE 9     y   a   PARTICLEONE%VELMAXINIT+MODULEPARTICLEONE N   	  k       VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE Q   t  Y   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE S   Í  @   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%MASS+MODULEPARTICLEONE Z     @   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%TEMPERATURE+MODULEPARTICLEONE 9   M  u   a   PARTICLEONE%VELRANINIT+MODULEPARTICLEONE J   Â  W       VELOCITYRANDOMINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE M     Y   a   VELOCITYRANDOMINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE L   r  @   a   VELOCITYRANDOMINITIALIZATIONPARTICLEONE%V+MODULEPARTICLEONE 9   ²  x   a   PARTICLEONE%ACCINPINIT+MODULEPARTICLEONE M   *  h       ACCELERATIONINPUTINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE P     Y   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE P   ë  @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AX+MODULEPARTICLEONE P   +  @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AY+MODULEPARTICLEONE P   k  @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AZ+MODULEPARTICLEONE 5   «  h   a   PARTICLEONE%POSRES+MODULEPARTICLEONE =     ]       POSITIONRESCALEPARTICLEONE+MODULEPARTICLEONE @   p  Y   a   POSITIONRESCALEPARTICLEONE%PO+MODULEPARTICLEONE E   É  @   a   POSITIONRESCALEPARTICLEONE%XFACTOR+MODULEPARTICLEONE 5   	  h   a   PARTICLEONE%VELRES+MODULEPARTICLEONE =   q  ]       VELOCITYRESCALEPARTICLEONE+MODULEPARTICLEONE @   Î  Y   a   VELOCITYRESCALEPARTICLEONE%PO+MODULEPARTICLEONE E   '  @   a   VELOCITYRESCALEPARTICLEONE%VFACTOR+MODULEPARTICLEONE 5   g  l   a   PARTICLEONE%ACCRES+MODULEPARTICLEONE A   Ó  ]       ACCELERATIONRESCALEPARTICLEONE+MODULEPARTICLEONE D   0  Y   a   ACCELERATIONRESCALEPARTICLEONE%PO+MODULEPARTICLEONE I     @   a   ACCELERATIONRESCALEPARTICLEONE%AFACTOR+MODULEPARTICLEONE 5   É  _   a   PARTICLEONE%ENERGY+MODULEPARTICLEONE 4   (  o       ENERGYPARTICLEONE+MODULEPARTICLEONE 7     Y   a   ENERGYPARTICLEONE%PO+MODULEPARTICLEONE 9   ð  @   a   ENERGYPARTICLEONE%MASS+MODULEPARTICLEONE <   0  @   a   ENERGYPARTICLEONE%VFACTOR+MODULEPARTICLEONE 3   p  ]   a   PARTICLEONE%COPY+MODULEPARTICLEONE 2   Í  Z       COPYPARTICLEONE+MODULEPARTICLEONE 6   '  Y   a   COPYPARTICLEONE%POD+MODULEPARTICLEONE 6     Y   a   COPYPARTICLEONE%POC+MODULEPARTICLEONE 3   Ù  ]   a   PARTICLEONE%SWAP+MODULEPARTICLEONE 2   6  Z       SWAPPARTICLEONE+MODULEPARTICLEONE 6     Y   a   SWAPPARTICLEONE%POD+MODULEPARTICLEONE 6   é  Y   a   SWAPPARTICLEONE%POC+MODULEPARTICLEONE 8   B  b   a   PARTICLEONE%WEIGHTP2C+MODULEPARTICLEONE 7   ¤  X       WEIGHTP2CPARTICLEONE+MODULEPARTICLEONE :   ü  Y   a   WEIGHTP2CPARTICLEONE%PO+MODULEPARTICLEONE :   U  V   a   WEIGHTP2CPARTICLEONE%FO+MODULEPARTICLEONE :   «  o   a   PARTICLEONE%WEIGHTC2PES+MODULEPARTICLEONE D     `       WEIGHTC2PELECTROSTATICPARTICLEONE+MODULEPARTICLEONE G   z  Y   a   WEIGHTC2PELECTROSTATICPARTICLEONE%PO+MODULEPARTICLEONE G   Ó  S   a   WEIGHTC2PELECTROSTATICPARTICLEONE%FG+MODULEPARTICLEONE G   &  @   a   WEIGHTC2PELECTROSTATICPARTICLEONE%EX+MODULEPARTICLEONE 5   f  j   a   PARTICLEONE%MOVEES+MODULEPARTICLEONE ?   Ð  X       MOVEELECTROSTATICPARTICLEONE+MODULEPARTICLEONE B   (  Y   a   MOVEELECTROSTATICPARTICLEONE%PO+MODULEPARTICLEONE B     @   a   MOVEELECTROSTATICPARTICLEONE%EX+MODULEPARTICLEONE :   Á  q   a   PARTICLEONE%WEIGHTC2PEM+MODULEPARTICLEONE F   2          WEIGHTC2PELECTROMAGNETICPARTICLEONE+MODULEPARTICLEONE I   º   Y   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%PO+MODULEPARTICLEONE I   !  S   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%FG+MODULEPARTICLEONE I   f!  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EX+MODULEPARTICLEONE I   ¦!  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EY+MODULEPARTICLEONE I   æ!  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EZ+MODULEPARTICLEONE I   &"  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BX+MODULEPARTICLEONE I   f"  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BY+MODULEPARTICLEONE I   ¦"  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BZ+MODULEPARTICLEONE 5   æ"  l   a   PARTICLEONE%MOVEEM+MODULEPARTICLEONE A   R#         MOVEELECTROMAGNETICPARTICLEONE+MODULEPARTICLEONE D   Ò#  Y   a   MOVEELECTROMAGNETICPARTICLEONE%PO+MODULEPARTICLEONE D   +$  @   a   MOVEELECTROMAGNETICPARTICLEONE%EX+MODULEPARTICLEONE D   k$  @   a   MOVEELECTROMAGNETICPARTICLEONE%EY+MODULEPARTICLEONE D   «$  @   a   MOVEELECTROMAGNETICPARTICLEONE%EZ+MODULEPARTICLEONE D   ë$  @   a   MOVEELECTROMAGNETICPARTICLEONE%BX+MODULEPARTICLEONE D   +%  @   a   MOVEELECTROMAGNETICPARTICLEONE%BY+MODULEPARTICLEONE D   k%  @   a   MOVEELECTROMAGNETICPARTICLEONE%BZ+MODULEPARTICLEONE 1   «%  H   a   MCCPARTICLEONE%MIU+MODULETYPEMCC 7   ó%  H   a   MCCPARTICLEONE%MASSRATIO+MODULETYPEMCC 4   ;&  H   a   MCCPARTICLEONE%ENERGY+MODULETYPEMCC 2   &  ¨   a   MCCPARTICLEONE%BETA+MODULETYPEMCC 0   +'  H   a   MCCPARTICLEONE%GX+MODULETYPEMCC 0   s'  H   a   MCCPARTICLEONE%GY+MODULETYPEMCC 0   »'  H   a   MCCPARTICLEONE%GZ+MODULETYPEMCC 2   (  H   a   MCCPARTICLEONE%GPER+MODULETYPEMCC /   K(  H   a   MCCPARTICLEONE%G+MODULETYPEMCC 1   (  a   a   MCCPARTICLEONE%POT+MODULETYPEMCC 4   ô(  ¥   a   MCCPARTICLEONE%NPONEW+MODULETYPEMCC 1   )  ª   a   MCCPARTICLEONE%PON+MODULETYPEMCC 3   C*         PARTICLEONEINDEX+MODULEPARTICLEONE ?   Ò*  a   a   PARTICLEONEINDEX%PARTICLEONE+MODULEPARTICLEONE 9   3+  H   a   PARTICLEONEINDEX%INDEX+MODULEPARTICLEONE =   {+  q   a   PARTICLEONEINDEX%INDEXINIT+MODULEPARTICLEONE F   ì+  \       INDEXINITIALIZATIONPARTICLEONEINDEX+MODULEPARTICLEONE J   H,  ^   a   INDEXINITIALIZATIONPARTICLEONEINDEX%POI+MODULEPARTICLEONE L   ¦,  @   a   INDEXINITIALIZATIONPARTICLEONEINDEX%INDEX+MODULEPARTICLEONE 8   æ,  b   a   PARTICLEONEINDEX%COPY+MODULEPARTICLEONE 7   H-  Z       COPYPARTICLEONEINDEX+MODULEPARTICLEONE ;   ¢-  ^   a   COPYPARTICLEONEINDEX%POD+MODULEPARTICLEONE ;    .  Y   a   COPYPARTICLEONEINDEX%POC+MODULEPARTICLEONE 8   Y.  b   a   PARTICLEONEINDEX%SWAP+MODULEPARTICLEONE 7   ».  Z       SWAPPARTICLEONEINDEX+MODULEPARTICLEONE ;   /  ^   a   SWAPPARTICLEONEINDEX%POD+MODULEPARTICLEONE ;   s/  Y   a   SWAPPARTICLEONEINDEX%POC+MODULEPARTICLEONE 4   Ì/  b   a   MCCPARTICLEONE%SELECT+MODULETYPEMCC 3   .0  Z       SELECTMCCPARTICLEONE+MODULETYPEMCC 8   0  \   a   SELECTMCCPARTICLEONE%MCPO+MODULETYPEMCC 6   ä0  Y   a   SELECTMCCPARTICLEONE%PO+MODULETYPEMCC 5   =1  b   a   MCCPARTICLEONE%UPDATER+MODULETYPEMCC 3   1  b       UPDATEMCCPARTICLEONE+MODULETYPEMCC 8   2  \   a   UPDATEMCCPARTICLEONE%MCPO+MODULETYPEMCC 6   ]2  V   a   UPDATEMCCPARTICLEONE%SO+MODULETYPEMCC 6   ³2  T   a   UPDATEMCCPARTICLEONE%GO+MODULETYPEMCC =   3  j   a   MCCPARTICLEONE%VELOCITYUPDATER+MODULETYPEMCC ;   q3  _       UPDATEVELOCITYMCCPARTICLEONE+MODULETYPEMCC @   Ð3  \   a   UPDATEVELOCITYMCCPARTICLEONE%MCPO+MODULETYPEMCC C   ,4  @   a   UPDATEVELOCITYMCCPARTICLEONE%VFACTOR+MODULETYPEMCC (   l4  é       SPECYONE+MODULESPECYONE -   U5  P   a   SPECYONE%NAME+MODULESPECYONE 3   ¥5  ¤   a   SPECYONE%SPECYINDEX+MODULESPECYONE 1   I6  ¤   a   SPECYONE%GASINDEX+MODULESPECYONE /   í6  H   a   SPECYONE%CHARGE+MODULESPECYONE -   57  H   a   SPECYONE%MASS+MODULESPECYONE /   }7  H   a   SPECYONE%RADIUS+MODULESPECYONE .   Å7  H   a   SPECYONE%NATOM+MODULESPECYONE 4   8  H   a   SPECYONE%INITDENSITY+MODULESPECYONE 0   U8  H   a   SPECYONE%DENSITY+MODULESPECYONE 8   8  H   a   SPECYONE%INITTEMPERATURE+MODULESPECYONE 4   å8  H   a   SPECYONE%TEMPERATURE+MODULESPECYONE &   -9  ñ       GASONE+MODULESPECYONE +   :  P   a   GASONE%NAME+MODULESPECYONE .   n:  ¤   a   GASONE%MCMODEL+MODULESPECYONE )   ;  H   a   GASONE%NS+MODULESPECYONE 1   Z;  H   a   GASONE%INDEXSTART+MODULESPECYONE +   ¢;  H   a   GASONE%MASS+MODULESPECYONE -   ê;  H   a   GASONE%RADIUS+MODULESPECYONE ,   2<  H   a   GASONE%NATOM+MODULESPECYONE .   z<  H   a   GASONE%BETAMAX+MODULESPECYONE 2   Â<  H   a   GASONE%INITDENSITY+MODULESPECYONE .   
=  H   a   GASONE%DENSITY+MODULESPECYONE 6   R=  H   a   GASONE%INITTEMPERATURE+MODULESPECYONE 2   =  H   a   GASONE%TEMPERATURE+MODULESPECYONE *   â=         REACTIONONE+MODULETYPEMCC 3   p>  H   a   REACTIONONE%REACTANT+MODULETYPEMCC 7   ¸>  H   a   REACTIONONE%REACTIONTYPE+MODULETYPEMCC 4    ?  H   a   REACTIONONE%RESULTANT+MODULETYPEMCC 4   H?  H   a   REACTIONONE%THRESHOLD+MODULETYPEMCC %   ?         FIELDONE+MODULEFIELD (   @  §   a   FIELDONE%NX+MODULEFIELD (   ·@  ¤   a   FIELDONE%DX+MODULEFIELD (   [A  ª   a   FIELDONE%DT+MODULEFIELD ,   B  ¬   a   FIELDONE%RHOONE+MODULEFIELD ,   ±B  ¬   a   FIELDONE%CHIONE+MODULEFIELD "   ]C  Ç       FIELD+MODULEFIELD %   $D  §   a   FIELD%NX+MODULEFIELD %   ËD  ¤   a   FIELD%DX+MODULEFIELD %   oE  ª   a   FIELD%DT+MODULEFIELD %   F    a   FIELD%EX+MODULEFIELD %   %G    a   FIELD%EY+MODULEFIELD %   1H    a   FIELD%EZ+MODULEFIELD %   =I    a   FIELD%BX+MODULEFIELD %   IJ    a   FIELD%BY+MODULEFIELD %   UK    a   FIELD%BZ+MODULEFIELD &   aL  ¬   a   FIELD%RHO+MODULEFIELD &   M  ¬   a   FIELD%PHI+MODULEFIELD &   ¹M  ¬   a   FIELD%CHI+MODULEFIELD '   eN  W   a   FIELD%DUMP+MODULEFIELD &   ¼N  Z       DUMPFIELD+MODULEFIELD )   O  S   a   DUMPFIELD%FG+MODULEFIELD +   iO  @   a   DUMPFIELD%MODE+MODULEFIELD '   ©O  W   a   FIELD%LOAD+MODULEFIELD &    P  \       LOADFIELD+MODULEFIELD )   \P  S   a   LOADFIELD%FG+MODULEFIELD -   ¯P  @   a   LOADFIELD%STATUS+MODULEFIELD (   ïP        MCCBUNDLE+MODULETYPEMCC .   ðQ  H   a   MCCBUNDLE%MODEL+MODULETYPEMCC 2   8R  ¥   a   MCCBUNDLE%NREACTION+MODULETYPEMCC /   ÝR  ¥   a   MCCBUNDLE%NSIGMA+MODULETYPEMCC 2   S  H   a   MCCBUNDLE%ENERGYMIN+MODULETYPEMCC 7   ÊS  H   a   MCCBUNDLE%ENERGYINTERVAL+MODULETYPEMCC 2   T  H   a   MCCBUNDLE%ENERGYMAX+MODULETYPEMCC 7   ZT  H   a   MCCBUNDLE%COLLISIONRATIO+MODULETYPEMCC 1   ¢T  H   a   MCCBUNDLE%SIGMAMAX+MODULETYPEMCC +   êT  Ì   a   MCCBUNDLE%SO+MODULETYPEMCC +   ¶U  È   a   MCCBUNDLE%GO+MODULETYPEMCC 1   ~V  ¥   a   MCCBUNDLE%REACTION+MODULETYPEMCC 2   #W  ¬   a   MCCBUNDLE%PROBILITY+MODULETYPEMCC -   ÏW  [   a   MCCBUNDLE%DUMP+MODULETYPEMCC ,   *X  Q       DUMPMCCBUNDLE+MODULETYPEMCC 0   {X  W   a   DUMPMCCBUNDLE%MCB+MODULETYPEMCC 4   ÒX  ^       POSTCOLLISIONVELOCITY+MODULETYPEMCC 9   0Y  \   a   POSTCOLLISIONVELOCITY%MCPO+MODULETYPEMCC ;   Y  @      POSTCOLLISIONVELOCITY%COSKAI+MODULETYPEMCC *   ÌY  \       COSKAIVAHEDI+MCCENERGYKAI 1   (Z  @   a   COSKAIVAHEDI%ENERGY+MCCENERGYKAI /   hZ  \       CREATENERGYVEHADI+MCCENERGYKAI 6   ÄZ  @   a   CREATENERGYVEHADI%ENERGY+MCCENERGYKAI (   [  j       SELECTCOLLISIONELECTRON -   n[  \   a   SELECTCOLLISIONELECTRON%MCPO +   Ê[  V   a   SELECTCOLLISIONELECTRON%SO +    \  T   a   SELECTCOLLISIONELECTRON%GO +   t\  Y   a   SELECTCOLLISIONELECTRON%RO 