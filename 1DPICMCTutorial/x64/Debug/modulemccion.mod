  '\  â   k820309              19.1        äí_                                                                                                          
       D:\Source\iPM1D2019\1DPICMCTutorial\code\mc\MCCIon.F90 MODULEMCCION                                                     
                                                           
                         @               A                'è                    #REACTIONINDEX    #PARTICLEANNIHILATION    #PARTICLECREATION    #POI    #MIU i   #MASSRATIO j   #ENERGY k   #BETA l   #GX m   #GY n   #GZ o   #GPER p   #G q   #POT r   #NPONEW s   #PON t   #SELECT    #UPDATER    #VELOCITYUPDATER                                                                                                                                              0                                                                                                                                                                                                                                                                                                                                                     8                    #PARTICLEONE                                          y#PARTICLEONE                                                                     @                                '8                    #X 	   #VX 
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
      1         À                               ?              
    #COPYPARTICLEONE @   #         @                                  @                    #POD A   #POC B             
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
                                                
                                 0.D0                                               m     8       	   
                                              n     @       
   
                                              o     H          
                                              p     P          
                                              q     X          
                                               r     8       `              #PARTICLEONE                                              s                                                                                               0                                              t                    @             #PARTICLEONEINDEX u             &                                                             @                          u     '@                    #PARTICLEONE v   #INDEX w   #INDEXINIT x   #COPY |   #SWAP                                                v     8                      #PARTICLEONE                                               w     8             1         À                               x                  #INDEXINITIALIZATIONPARTICLEONEINDEX y   #         @                                  y                    #POI z   #INDEX {             
                                z     @               #PARTICLEONEINDEX u             
                                 {           1         À                               |              
    #COPYPARTICLEONEINDEX }   #         @                                  }                    #POD ~   #POC              
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
                        @                               '°                    #NAME    #SPECYINDEX    #GASINDEX    #CHARGE    #MASS    #RADIUS    #NATOM    #INITDENSITY    #DENSITY    #INITTEMPERATURE    #TEMPERATURE                                                    c                                                                     d                                                                             ÿÿÿÿÿÿÿÿ                                                           h                                                                             ÿÿÿÿÿÿÿÿ                                                            p          
                                                   x          
                                                             
                                                             
                                                             
                                                          	   
                                                           
   
                                                   ¨          
                     @                               '°                    #NAME    #MCMODEL    #NS     #INDEXSTART ¡   #MASS ¢   #RADIUS £   #NATOM ¤   #BETAMAX ¥   #INITDENSITY ¦   #DENSITY §   #INITTEMPERATURE ¨   #TEMPERATURE ©                                                   c                                                                     d                                                                             ÿÿÿÿÿÿÿÿ                                                             h                                                        ¡     l                                                        ¢     p          
                                              £     x          
                                              ¤               
                                              ¥               
                                              ¦            	   
                                              §            
   
                                              ¨                
                                              ©     ¨          
                     @                          ª     '                    #REACTANT «   #REACTIONTYPE ¬   #RESULTANT ­   #THRESHOLD ®                                              «                                                              ¬                                                             ­                                                             ®               
                     @                          K     '(                   #NX ¯   #DX °   #DT ±   #RHOONE ²   #CHIONE ³                                             ¯                                                                                 A              321                                              °              
                                                  
                  {®Gáz4?                                                      ±              
                                                 
                 »½×Ùß|Û=        1.D-10                                               ²     A                       
  p          & p        p A          p A                                                                   ³     A       
                
  p          & p        p A          p A                                          @                          P     '`Z                   #NX ´   #DX µ   #DT ¶   #EX ·   #EY ¸   #EZ ¹   #BX º   #BY »   #BZ ¼   #RHO ½   #PHI ¾   #CHI ¿   #DUMP À   #LOAD Ä                                             ´                                                                                 A              321                                              µ              
                                                  
                  {®Gáz4?                                                      ¶              
                                                 
                 »½×Ùß|Û=        1.D-10                                              ·     A                      
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          ¸     A       
               
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          ¹     A      (               
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          º     A      0               
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          »     A      8(               
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                          ¼     A      @2             	  
  p          & p        p A          p A                                       A            A                    
                                 0.D0                                           ½     A      H<             
   
  p          & p        p A          p A                                                                   ¾     A      PF                
  p          & p        p A          p A                                                                   ¿     A      XP                
  p          & p        p A          p A                        1         À                                À                  #DUMPFIELD Á   #         @                                   Á                    #FG Â   #MODE Ã             
                                Â     `Z              #FIELD P                                             Ã            1         À                                Ä                  #LOADFIELD Å   #         @                                   Å                    #FG Æ   #STATUS Ç             
                                Æ     `Z              #FIELD P             
                                 Ç                              @               A           È     'ð                    #MODEL É   #NREACTION Ê   #NSIGMA Ë   #ENERGYMIN Ì   #ENERGYINTERVAL Í   #ENERGYMAX Î   #COLLISIONRATIO Ï   #SIGMAMAX Ð   #SO Ñ   #GO Ò   #REACTION Ó   #PROBILITY Ô   #DUMP Õ                                              É                                                             Ê                                                                                               0                                              Ë                                                                                               0                                               Ì               
                                              Í               
                                              Î                
                                              Ï     (          
                                              Ð     0          
                                             Ñ     °       8       	      #SPECYONE                                          y#SPECYONE                                                                                             Ò     °       @       
      #GASONE                                          y#GASONE                                                                                             Ó            H                    #REACTIONONE ª             &                                                                                    Ô                             
            &                   &                                           1         À                                Õ                  #DUMPMCCBUNDLE Ö   #         @                                   Ö                    #MCB ×             
                                ×     ð               #MCCBUNDLE È   #         @                                  Ø                    #MCPO Ù   #COSKAI Ú             
                                 Ù     è               #MCCPARTICLEONE    "                                        Ú       
     %         @  @                            Û                    
       #ENERGY Ü                                             Ü     
       #         @                                   Ý                    #MCPO Þ   #SO ß   #GO à   #RO á             
D @                               Þ     è               #MCCPARTICLEONE              
                                  ß     °              #SPECYONE              
                                  à     °              #GASONE              
  @                               á                   #REACTIONONE ª          L      fn#fn    ì   @   J   MODULETYPEMCC    ,  @   J   MCCENERGYKAI -   l  ?      MCCPARTICLEONE+MODULETYPEMCC ;   «  ¥   a   MCCPARTICLEONE%REACTIONINDEX+MODULETYPEMCC B   P  ¤   a   MCCPARTICLEONE%PARTICLEANNIHILATION+MODULETYPEMCC >   ô  ¤   a   MCCPARTICLEONE%PARTICLECREATION+MODULETYPEMCC 1     Ò   a   MCCPARTICLEONE%POI+MODULETYPEMCC .   j  a      PARTICLEONE+MODULEPARTICLEONE 0   Ë  H   a   PARTICLEONE%X+MODULEPARTICLEONE 1     H   a   PARTICLEONE%VX+MODULEPARTICLEONE 1   [  H   a   PARTICLEONE%VY+MODULEPARTICLEONE 1   £  H   a   PARTICLEONE%VZ+MODULEPARTICLEONE 1   ë  H   a   PARTICLEONE%AX+MODULEPARTICLEONE 1   3  H   a   PARTICLEONE%AY+MODULEPARTICLEONE 1   {  H   a   PARTICLEONE%AZ+MODULEPARTICLEONE 6   Ã  u   a   PARTICLEONE%POSINIT+MODULEPARTICLEONE J   8	         POSITIONRANDOMINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE M   ¸	  Y   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE M   
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%XL+MODULEPARTICLEONE M   Q
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%XU+MODULEPARTICLEONE M   
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%YL+MODULEPARTICLEONE M   Ñ
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%YU+MODULEPARTICLEONE M     @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%ZL+MODULEPARTICLEONE M   Q  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%ZU+MODULEPARTICLEONE 9     t   a   PARTICLEONE%VELINPINIT+MODULEPARTICLEONE I     h       VELOCITYINPUTINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE L   m  Y   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE L   Æ  @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VX+MODULEPARTICLEONE L     @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VY+MODULEPARTICLEONE L   F  @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VZ+MODULEPARTICLEONE 9     y   a   PARTICLEONE%VELMAXINIT+MODULEPARTICLEONE N   ÿ  k       VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE Q   j  Y   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE S   Ã  @   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%MASS+MODULEPARTICLEONE Z     @   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%TEMPERATURE+MODULEPARTICLEONE 9   C  u   a   PARTICLEONE%VELRANINIT+MODULEPARTICLEONE J   ¸  W       VELOCITYRANDOMINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE M     Y   a   VELOCITYRANDOMINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE L   h  @   a   VELOCITYRANDOMINITIALIZATIONPARTICLEONE%V+MODULEPARTICLEONE 9   ¨  x   a   PARTICLEONE%ACCINPINIT+MODULEPARTICLEONE M      h       ACCELERATIONINPUTINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE P     Y   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE P   á  @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AX+MODULEPARTICLEONE P   !  @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AY+MODULEPARTICLEONE P   a  @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AZ+MODULEPARTICLEONE 5   ¡  h   a   PARTICLEONE%POSRES+MODULEPARTICLEONE =   	  ]       POSITIONRESCALEPARTICLEONE+MODULEPARTICLEONE @   f  Y   a   POSITIONRESCALEPARTICLEONE%PO+MODULEPARTICLEONE E   ¿  @   a   POSITIONRESCALEPARTICLEONE%XFACTOR+MODULEPARTICLEONE 5   ÿ  h   a   PARTICLEONE%VELRES+MODULEPARTICLEONE =   g  ]       VELOCITYRESCALEPARTICLEONE+MODULEPARTICLEONE @   Ä  Y   a   VELOCITYRESCALEPARTICLEONE%PO+MODULEPARTICLEONE E     @   a   VELOCITYRESCALEPARTICLEONE%VFACTOR+MODULEPARTICLEONE 5   ]  l   a   PARTICLEONE%ACCRES+MODULEPARTICLEONE A   É  ]       ACCELERATIONRESCALEPARTICLEONE+MODULEPARTICLEONE D   &  Y   a   ACCELERATIONRESCALEPARTICLEONE%PO+MODULEPARTICLEONE I     @   a   ACCELERATIONRESCALEPARTICLEONE%AFACTOR+MODULEPARTICLEONE 5   ¿  _   a   PARTICLEONE%ENERGY+MODULEPARTICLEONE 4     o       ENERGYPARTICLEONE+MODULEPARTICLEONE 7     Y   a   ENERGYPARTICLEONE%PO+MODULEPARTICLEONE 9   æ  @   a   ENERGYPARTICLEONE%MASS+MODULEPARTICLEONE <   &  @   a   ENERGYPARTICLEONE%VFACTOR+MODULEPARTICLEONE 3   f  ]   a   PARTICLEONE%COPY+MODULEPARTICLEONE 2   Ã  Z       COPYPARTICLEONE+MODULEPARTICLEONE 6     Y   a   COPYPARTICLEONE%POD+MODULEPARTICLEONE 6   v  Y   a   COPYPARTICLEONE%POC+MODULEPARTICLEONE 3   Ï  ]   a   PARTICLEONE%SWAP+MODULEPARTICLEONE 2   ,  Z       SWAPPARTICLEONE+MODULEPARTICLEONE 6     Y   a   SWAPPARTICLEONE%POD+MODULEPARTICLEONE 6   ß  Y   a   SWAPPARTICLEONE%POC+MODULEPARTICLEONE 8   8  b   a   PARTICLEONE%WEIGHTP2C+MODULEPARTICLEONE 7     X       WEIGHTP2CPARTICLEONE+MODULEPARTICLEONE :   ò  Y   a   WEIGHTP2CPARTICLEONE%PO+MODULEPARTICLEONE :   K  V   a   WEIGHTP2CPARTICLEONE%FO+MODULEPARTICLEONE :   ¡  o   a   PARTICLEONE%WEIGHTC2PES+MODULEPARTICLEONE D     `       WEIGHTC2PELECTROSTATICPARTICLEONE+MODULEPARTICLEONE G   p  Y   a   WEIGHTC2PELECTROSTATICPARTICLEONE%PO+MODULEPARTICLEONE G   É  S   a   WEIGHTC2PELECTROSTATICPARTICLEONE%FG+MODULEPARTICLEONE G     @   a   WEIGHTC2PELECTROSTATICPARTICLEONE%EX+MODULEPARTICLEONE 5   \  j   a   PARTICLEONE%MOVEES+MODULEPARTICLEONE ?   Æ  X       MOVEELECTROSTATICPARTICLEONE+MODULEPARTICLEONE B     Y   a   MOVEELECTROSTATICPARTICLEONE%PO+MODULEPARTICLEONE B   w  @   a   MOVEELECTROSTATICPARTICLEONE%EX+MODULEPARTICLEONE :   ·  q   a   PARTICLEONE%WEIGHTC2PEM+MODULEPARTICLEONE F   (          WEIGHTC2PELECTROMAGNETICPARTICLEONE+MODULEPARTICLEONE I   °   Y   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%PO+MODULEPARTICLEONE I   	!  S   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%FG+MODULEPARTICLEONE I   \!  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EX+MODULEPARTICLEONE I   !  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EY+MODULEPARTICLEONE I   Ü!  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EZ+MODULEPARTICLEONE I   "  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BX+MODULEPARTICLEONE I   \"  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BY+MODULEPARTICLEONE I   "  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BZ+MODULEPARTICLEONE 5   Ü"  l   a   PARTICLEONE%MOVEEM+MODULEPARTICLEONE A   H#         MOVEELECTROMAGNETICPARTICLEONE+MODULEPARTICLEONE D   È#  Y   a   MOVEELECTROMAGNETICPARTICLEONE%PO+MODULEPARTICLEONE D   !$  @   a   MOVEELECTROMAGNETICPARTICLEONE%EX+MODULEPARTICLEONE D   a$  @   a   MOVEELECTROMAGNETICPARTICLEONE%EY+MODULEPARTICLEONE D   ¡$  @   a   MOVEELECTROMAGNETICPARTICLEONE%EZ+MODULEPARTICLEONE D   á$  @   a   MOVEELECTROMAGNETICPARTICLEONE%BX+MODULEPARTICLEONE D   !%  @   a   MOVEELECTROMAGNETICPARTICLEONE%BY+MODULEPARTICLEONE D   a%  @   a   MOVEELECTROMAGNETICPARTICLEONE%BZ+MODULEPARTICLEONE 1   ¡%  H   a   MCCPARTICLEONE%MIU+MODULETYPEMCC 7   é%  H   a   MCCPARTICLEONE%MASSRATIO+MODULETYPEMCC 4   1&  H   a   MCCPARTICLEONE%ENERGY+MODULETYPEMCC 2   y&  ¨   a   MCCPARTICLEONE%BETA+MODULETYPEMCC 0   !'  H   a   MCCPARTICLEONE%GX+MODULETYPEMCC 0   i'  H   a   MCCPARTICLEONE%GY+MODULETYPEMCC 0   ±'  H   a   MCCPARTICLEONE%GZ+MODULETYPEMCC 2   ù'  H   a   MCCPARTICLEONE%GPER+MODULETYPEMCC /   A(  H   a   MCCPARTICLEONE%G+MODULETYPEMCC 1   (  a   a   MCCPARTICLEONE%POT+MODULETYPEMCC 4   ê(  ¥   a   MCCPARTICLEONE%NPONEW+MODULETYPEMCC 1   )  ª   a   MCCPARTICLEONE%PON+MODULETYPEMCC 3   9*         PARTICLEONEINDEX+MODULEPARTICLEONE ?   È*  a   a   PARTICLEONEINDEX%PARTICLEONE+MODULEPARTICLEONE 9   )+  H   a   PARTICLEONEINDEX%INDEX+MODULEPARTICLEONE =   q+  q   a   PARTICLEONEINDEX%INDEXINIT+MODULEPARTICLEONE F   â+  \       INDEXINITIALIZATIONPARTICLEONEINDEX+MODULEPARTICLEONE J   >,  ^   a   INDEXINITIALIZATIONPARTICLEONEINDEX%POI+MODULEPARTICLEONE L   ,  @   a   INDEXINITIALIZATIONPARTICLEONEINDEX%INDEX+MODULEPARTICLEONE 8   Ü,  b   a   PARTICLEONEINDEX%COPY+MODULEPARTICLEONE 7   >-  Z       COPYPARTICLEONEINDEX+MODULEPARTICLEONE ;   -  ^   a   COPYPARTICLEONEINDEX%POD+MODULEPARTICLEONE ;   ö-  Y   a   COPYPARTICLEONEINDEX%POC+MODULEPARTICLEONE 8   O.  b   a   PARTICLEONEINDEX%SWAP+MODULEPARTICLEONE 7   ±.  Z       SWAPPARTICLEONEINDEX+MODULEPARTICLEONE ;   /  ^   a   SWAPPARTICLEONEINDEX%POD+MODULEPARTICLEONE ;   i/  Y   a   SWAPPARTICLEONEINDEX%POC+MODULEPARTICLEONE 4   Â/  b   a   MCCPARTICLEONE%SELECT+MODULETYPEMCC 3   $0  Z       SELECTMCCPARTICLEONE+MODULETYPEMCC 8   ~0  \   a   SELECTMCCPARTICLEONE%MCPO+MODULETYPEMCC 6   Ú0  Y   a   SELECTMCCPARTICLEONE%PO+MODULETYPEMCC 5   31  b   a   MCCPARTICLEONE%UPDATER+MODULETYPEMCC 3   1  b       UPDATEMCCPARTICLEONE+MODULETYPEMCC 8   ÷1  \   a   UPDATEMCCPARTICLEONE%MCPO+MODULETYPEMCC 6   S2  V   a   UPDATEMCCPARTICLEONE%SO+MODULETYPEMCC 6   ©2  T   a   UPDATEMCCPARTICLEONE%GO+MODULETYPEMCC =   ý2  j   a   MCCPARTICLEONE%VELOCITYUPDATER+MODULETYPEMCC ;   g3  _       UPDATEVELOCITYMCCPARTICLEONE+MODULETYPEMCC @   Æ3  \   a   UPDATEVELOCITYMCCPARTICLEONE%MCPO+MODULETYPEMCC C   "4  @   a   UPDATEVELOCITYMCCPARTICLEONE%VFACTOR+MODULETYPEMCC (   b4  é       SPECYONE+MODULESPECYONE -   K5  P   a   SPECYONE%NAME+MODULESPECYONE 3   5  ¤   a   SPECYONE%SPECYINDEX+MODULESPECYONE 1   ?6  ¤   a   SPECYONE%GASINDEX+MODULESPECYONE /   ã6  H   a   SPECYONE%CHARGE+MODULESPECYONE -   +7  H   a   SPECYONE%MASS+MODULESPECYONE /   s7  H   a   SPECYONE%RADIUS+MODULESPECYONE .   »7  H   a   SPECYONE%NATOM+MODULESPECYONE 4   8  H   a   SPECYONE%INITDENSITY+MODULESPECYONE 0   K8  H   a   SPECYONE%DENSITY+MODULESPECYONE 8   8  H   a   SPECYONE%INITTEMPERATURE+MODULESPECYONE 4   Û8  H   a   SPECYONE%TEMPERATURE+MODULESPECYONE &   #9  ñ       GASONE+MODULESPECYONE +   :  P   a   GASONE%NAME+MODULESPECYONE .   d:  ¤   a   GASONE%MCMODEL+MODULESPECYONE )   ;  H   a   GASONE%NS+MODULESPECYONE 1   P;  H   a   GASONE%INDEXSTART+MODULESPECYONE +   ;  H   a   GASONE%MASS+MODULESPECYONE -   à;  H   a   GASONE%RADIUS+MODULESPECYONE ,   (<  H   a   GASONE%NATOM+MODULESPECYONE .   p<  H   a   GASONE%BETAMAX+MODULESPECYONE 2   ¸<  H   a   GASONE%INITDENSITY+MODULESPECYONE .    =  H   a   GASONE%DENSITY+MODULESPECYONE 6   H=  H   a   GASONE%INITTEMPERATURE+MODULESPECYONE 2   =  H   a   GASONE%TEMPERATURE+MODULESPECYONE *   Ø=         REACTIONONE+MODULETYPEMCC 3   f>  H   a   REACTIONONE%REACTANT+MODULETYPEMCC 7   ®>  H   a   REACTIONONE%REACTIONTYPE+MODULETYPEMCC 4   ö>  H   a   REACTIONONE%RESULTANT+MODULETYPEMCC 4   >?  H   a   REACTIONONE%THRESHOLD+MODULETYPEMCC %   ?         FIELDONE+MODULEFIELD (   @  §   a   FIELDONE%NX+MODULEFIELD (   ­@  ¤   a   FIELDONE%DX+MODULEFIELD (   QA  ª   a   FIELDONE%DT+MODULEFIELD ,   ûA  ¬   a   FIELDONE%RHOONE+MODULEFIELD ,   §B  ¬   a   FIELDONE%CHIONE+MODULEFIELD "   SC  Ç       FIELD+MODULEFIELD %   D  §   a   FIELD%NX+MODULEFIELD %   ÁD  ¤   a   FIELD%DX+MODULEFIELD %   eE  ª   a   FIELD%DT+MODULEFIELD %   F    a   FIELD%EX+MODULEFIELD %   G    a   FIELD%EY+MODULEFIELD %   'H    a   FIELD%EZ+MODULEFIELD %   3I    a   FIELD%BX+MODULEFIELD %   ?J    a   FIELD%BY+MODULEFIELD %   KK    a   FIELD%BZ+MODULEFIELD &   WL  ¬   a   FIELD%RHO+MODULEFIELD &   M  ¬   a   FIELD%PHI+MODULEFIELD &   ¯M  ¬   a   FIELD%CHI+MODULEFIELD '   [N  W   a   FIELD%DUMP+MODULEFIELD &   ²N  Z       DUMPFIELD+MODULEFIELD )   O  S   a   DUMPFIELD%FG+MODULEFIELD +   _O  @   a   DUMPFIELD%MODE+MODULEFIELD '   O  W   a   FIELD%LOAD+MODULEFIELD &   öO  \       LOADFIELD+MODULEFIELD )   RP  S   a   LOADFIELD%FG+MODULEFIELD -   ¥P  @   a   LOADFIELD%STATUS+MODULEFIELD (   åP        MCCBUNDLE+MODULETYPEMCC .   æQ  H   a   MCCBUNDLE%MODEL+MODULETYPEMCC 2   .R  ¥   a   MCCBUNDLE%NREACTION+MODULETYPEMCC /   ÓR  ¥   a   MCCBUNDLE%NSIGMA+MODULETYPEMCC 2   xS  H   a   MCCBUNDLE%ENERGYMIN+MODULETYPEMCC 7   ÀS  H   a   MCCBUNDLE%ENERGYINTERVAL+MODULETYPEMCC 2   T  H   a   MCCBUNDLE%ENERGYMAX+MODULETYPEMCC 7   PT  H   a   MCCBUNDLE%COLLISIONRATIO+MODULETYPEMCC 1   T  H   a   MCCBUNDLE%SIGMAMAX+MODULETYPEMCC +   àT  Ì   a   MCCBUNDLE%SO+MODULETYPEMCC +   ¬U  È   a   MCCBUNDLE%GO+MODULETYPEMCC 1   tV  ¥   a   MCCBUNDLE%REACTION+MODULETYPEMCC 2   W  ¬   a   MCCBUNDLE%PROBILITY+MODULETYPEMCC -   ÅW  [   a   MCCBUNDLE%DUMP+MODULETYPEMCC ,    X  Q       DUMPMCCBUNDLE+MODULETYPEMCC 0   qX  W   a   DUMPMCCBUNDLE%MCB+MODULETYPEMCC 4   ÈX  ^       POSTCOLLISIONVELOCITY+MODULETYPEMCC 9   &Y  \   a   POSTCOLLISIONVELOCITY%MCPO+MODULETYPEMCC ;   Y  @      POSTCOLLISIONVELOCITY%COSKAI+MODULETYPEMCC -   ÂY  \       ISOTROPICCOSKAI+MCCENERGYKAI 4   Z  @   a   ISOTROPICCOSKAI%ENERGY+MCCENERGYKAI #   ^Z  j       SELECTCOLLISIONION (   ÈZ  \   a   SELECTCOLLISIONION%MCPO &   $[  V   a   SELECTCOLLISIONION%SO &   z[  T   a   SELECTCOLLISIONION%GO &   Î[  Y   a   SELECTCOLLISIONION%RO 