  G\  â   k820309    9          19.0        ñ³^                                                                                                          
       C:\Users\physi\Source\Repos\PhysicsPlasma\iPM1D2019\1DPICMCTutorial\code\mc\MCCIon.F90 MODULEMCCION                                                     
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
  @                               á                   #REACTIONONE ª          l      fn#fn      @   J   MODULETYPEMCC    L  @   J   MCCENERGYKAI -     ?      MCCPARTICLEONE+MODULETYPEMCC ;   Ë  ¥   a   MCCPARTICLEONE%REACTIONINDEX+MODULETYPEMCC B   p  ¤   a   MCCPARTICLEONE%PARTICLEANNIHILATION+MODULETYPEMCC >     ¤   a   MCCPARTICLEONE%PARTICLECREATION+MODULETYPEMCC 1   ¸  Ò   a   MCCPARTICLEONE%POI+MODULETYPEMCC .     a      PARTICLEONE+MODULEPARTICLEONE 0   ë  H   a   PARTICLEONE%X+MODULEPARTICLEONE 1   3  H   a   PARTICLEONE%VX+MODULEPARTICLEONE 1   {  H   a   PARTICLEONE%VY+MODULEPARTICLEONE 1   Ã  H   a   PARTICLEONE%VZ+MODULEPARTICLEONE 1     H   a   PARTICLEONE%AX+MODULEPARTICLEONE 1   S  H   a   PARTICLEONE%AY+MODULEPARTICLEONE 1     H   a   PARTICLEONE%AZ+MODULEPARTICLEONE 6   ã  u   a   PARTICLEONE%POSINIT+MODULEPARTICLEONE J   X	         POSITIONRANDOMINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE M   Ø	  Y   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE M   1
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%XL+MODULEPARTICLEONE M   q
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%XU+MODULEPARTICLEONE M   ±
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%YL+MODULEPARTICLEONE M   ñ
  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%YU+MODULEPARTICLEONE M   1  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%ZL+MODULEPARTICLEONE M   q  @   a   POSITIONRANDOMINITIALIZATIONPARTICLEONE%ZU+MODULEPARTICLEONE 9   ±  t   a   PARTICLEONE%VELINPINIT+MODULEPARTICLEONE I   %  h       VELOCITYINPUTINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE L     Y   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE L   æ  @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VX+MODULEPARTICLEONE L   &  @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VY+MODULEPARTICLEONE L   f  @   a   VELOCITYINPUTINITIALIZATIONPARTICLEONE%VZ+MODULEPARTICLEONE 9   ¦  y   a   PARTICLEONE%VELMAXINIT+MODULEPARTICLEONE N     k       VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE Q     Y   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE S   ã  @   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%MASS+MODULEPARTICLEONE Z   #  @   a   VELOCITYMAXWELLIANINITIALIZATIONPARTICLEONE%TEMPERATURE+MODULEPARTICLEONE 9   c  u   a   PARTICLEONE%VELRANINIT+MODULEPARTICLEONE J   Ø  W       VELOCITYRANDOMINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE M   /  Y   a   VELOCITYRANDOMINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE L     @   a   VELOCITYRANDOMINITIALIZATIONPARTICLEONE%V+MODULEPARTICLEONE 9   È  x   a   PARTICLEONE%ACCINPINIT+MODULEPARTICLEONE M   @  h       ACCELERATIONINPUTINITIALIZATIONPARTICLEONE+MODULEPARTICLEONE P   ¨  Y   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%PO+MODULEPARTICLEONE P     @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AX+MODULEPARTICLEONE P   A  @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AY+MODULEPARTICLEONE P     @   a   ACCELERATIONINPUTINITIALIZATIONPARTICLEONE%AZ+MODULEPARTICLEONE 5   Á  h   a   PARTICLEONE%POSRES+MODULEPARTICLEONE =   )  ]       POSITIONRESCALEPARTICLEONE+MODULEPARTICLEONE @     Y   a   POSITIONRESCALEPARTICLEONE%PO+MODULEPARTICLEONE E   ß  @   a   POSITIONRESCALEPARTICLEONE%XFACTOR+MODULEPARTICLEONE 5     h   a   PARTICLEONE%VELRES+MODULEPARTICLEONE =     ]       VELOCITYRESCALEPARTICLEONE+MODULEPARTICLEONE @   ä  Y   a   VELOCITYRESCALEPARTICLEONE%PO+MODULEPARTICLEONE E   =  @   a   VELOCITYRESCALEPARTICLEONE%VFACTOR+MODULEPARTICLEONE 5   }  l   a   PARTICLEONE%ACCRES+MODULEPARTICLEONE A   é  ]       ACCELERATIONRESCALEPARTICLEONE+MODULEPARTICLEONE D   F  Y   a   ACCELERATIONRESCALEPARTICLEONE%PO+MODULEPARTICLEONE I     @   a   ACCELERATIONRESCALEPARTICLEONE%AFACTOR+MODULEPARTICLEONE 5   ß  _   a   PARTICLEONE%ENERGY+MODULEPARTICLEONE 4   >  o       ENERGYPARTICLEONE+MODULEPARTICLEONE 7   ­  Y   a   ENERGYPARTICLEONE%PO+MODULEPARTICLEONE 9     @   a   ENERGYPARTICLEONE%MASS+MODULEPARTICLEONE <   F  @   a   ENERGYPARTICLEONE%VFACTOR+MODULEPARTICLEONE 3     ]   a   PARTICLEONE%COPY+MODULEPARTICLEONE 2   ã  Z       COPYPARTICLEONE+MODULEPARTICLEONE 6   =  Y   a   COPYPARTICLEONE%POD+MODULEPARTICLEONE 6     Y   a   COPYPARTICLEONE%POC+MODULEPARTICLEONE 3   ï  ]   a   PARTICLEONE%SWAP+MODULEPARTICLEONE 2   L  Z       SWAPPARTICLEONE+MODULEPARTICLEONE 6   ¦  Y   a   SWAPPARTICLEONE%POD+MODULEPARTICLEONE 6   ÿ  Y   a   SWAPPARTICLEONE%POC+MODULEPARTICLEONE 8   X  b   a   PARTICLEONE%WEIGHTP2C+MODULEPARTICLEONE 7   º  X       WEIGHTP2CPARTICLEONE+MODULEPARTICLEONE :     Y   a   WEIGHTP2CPARTICLEONE%PO+MODULEPARTICLEONE :   k  V   a   WEIGHTP2CPARTICLEONE%FO+MODULEPARTICLEONE :   Á  o   a   PARTICLEONE%WEIGHTC2PES+MODULEPARTICLEONE D   0  `       WEIGHTC2PELECTROSTATICPARTICLEONE+MODULEPARTICLEONE G     Y   a   WEIGHTC2PELECTROSTATICPARTICLEONE%PO+MODULEPARTICLEONE G   é  S   a   WEIGHTC2PELECTROSTATICPARTICLEONE%FG+MODULEPARTICLEONE G   <  @   a   WEIGHTC2PELECTROSTATICPARTICLEONE%EX+MODULEPARTICLEONE 5   |  j   a   PARTICLEONE%MOVEES+MODULEPARTICLEONE ?   æ  X       MOVEELECTROSTATICPARTICLEONE+MODULEPARTICLEONE B   >  Y   a   MOVEELECTROSTATICPARTICLEONE%PO+MODULEPARTICLEONE B     @   a   MOVEELECTROSTATICPARTICLEONE%EX+MODULEPARTICLEONE :   ×  q   a   PARTICLEONE%WEIGHTC2PEM+MODULEPARTICLEONE F   H          WEIGHTC2PELECTROMAGNETICPARTICLEONE+MODULEPARTICLEONE I   Ð   Y   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%PO+MODULEPARTICLEONE I   )!  S   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%FG+MODULEPARTICLEONE I   |!  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EX+MODULEPARTICLEONE I   ¼!  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EY+MODULEPARTICLEONE I   ü!  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%EZ+MODULEPARTICLEONE I   <"  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BX+MODULEPARTICLEONE I   |"  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BY+MODULEPARTICLEONE I   ¼"  @   a   WEIGHTC2PELECTROMAGNETICPARTICLEONE%BZ+MODULEPARTICLEONE 5   ü"  l   a   PARTICLEONE%MOVEEM+MODULEPARTICLEONE A   h#         MOVEELECTROMAGNETICPARTICLEONE+MODULEPARTICLEONE D   è#  Y   a   MOVEELECTROMAGNETICPARTICLEONE%PO+MODULEPARTICLEONE D   A$  @   a   MOVEELECTROMAGNETICPARTICLEONE%EX+MODULEPARTICLEONE D   $  @   a   MOVEELECTROMAGNETICPARTICLEONE%EY+MODULEPARTICLEONE D   Á$  @   a   MOVEELECTROMAGNETICPARTICLEONE%EZ+MODULEPARTICLEONE D   %  @   a   MOVEELECTROMAGNETICPARTICLEONE%BX+MODULEPARTICLEONE D   A%  @   a   MOVEELECTROMAGNETICPARTICLEONE%BY+MODULEPARTICLEONE D   %  @   a   MOVEELECTROMAGNETICPARTICLEONE%BZ+MODULEPARTICLEONE 1   Á%  H   a   MCCPARTICLEONE%MIU+MODULETYPEMCC 7   	&  H   a   MCCPARTICLEONE%MASSRATIO+MODULETYPEMCC 4   Q&  H   a   MCCPARTICLEONE%ENERGY+MODULETYPEMCC 2   &  ¨   a   MCCPARTICLEONE%BETA+MODULETYPEMCC 0   A'  H   a   MCCPARTICLEONE%GX+MODULETYPEMCC 0   '  H   a   MCCPARTICLEONE%GY+MODULETYPEMCC 0   Ñ'  H   a   MCCPARTICLEONE%GZ+MODULETYPEMCC 2   (  H   a   MCCPARTICLEONE%GPER+MODULETYPEMCC /   a(  H   a   MCCPARTICLEONE%G+MODULETYPEMCC 1   ©(  a   a   MCCPARTICLEONE%POT+MODULETYPEMCC 4   
)  ¥   a   MCCPARTICLEONE%NPONEW+MODULETYPEMCC 1   ¯)  ª   a   MCCPARTICLEONE%PON+MODULETYPEMCC 3   Y*         PARTICLEONEINDEX+MODULEPARTICLEONE ?   è*  a   a   PARTICLEONEINDEX%PARTICLEONE+MODULEPARTICLEONE 9   I+  H   a   PARTICLEONEINDEX%INDEX+MODULEPARTICLEONE =   +  q   a   PARTICLEONEINDEX%INDEXINIT+MODULEPARTICLEONE F   ,  \       INDEXINITIALIZATIONPARTICLEONEINDEX+MODULEPARTICLEONE J   ^,  ^   a   INDEXINITIALIZATIONPARTICLEONEINDEX%POI+MODULEPARTICLEONE L   ¼,  @   a   INDEXINITIALIZATIONPARTICLEONEINDEX%INDEX+MODULEPARTICLEONE 8   ü,  b   a   PARTICLEONEINDEX%COPY+MODULEPARTICLEONE 7   ^-  Z       COPYPARTICLEONEINDEX+MODULEPARTICLEONE ;   ¸-  ^   a   COPYPARTICLEONEINDEX%POD+MODULEPARTICLEONE ;   .  Y   a   COPYPARTICLEONEINDEX%POC+MODULEPARTICLEONE 8   o.  b   a   PARTICLEONEINDEX%SWAP+MODULEPARTICLEONE 7   Ñ.  Z       SWAPPARTICLEONEINDEX+MODULEPARTICLEONE ;   +/  ^   a   SWAPPARTICLEONEINDEX%POD+MODULEPARTICLEONE ;   /  Y   a   SWAPPARTICLEONEINDEX%POC+MODULEPARTICLEONE 4   â/  b   a   MCCPARTICLEONE%SELECT+MODULETYPEMCC 3   D0  Z       SELECTMCCPARTICLEONE+MODULETYPEMCC 8   0  \   a   SELECTMCCPARTICLEONE%MCPO+MODULETYPEMCC 6   ú0  Y   a   SELECTMCCPARTICLEONE%PO+MODULETYPEMCC 5   S1  b   a   MCCPARTICLEONE%UPDATER+MODULETYPEMCC 3   µ1  b       UPDATEMCCPARTICLEONE+MODULETYPEMCC 8   2  \   a   UPDATEMCCPARTICLEONE%MCPO+MODULETYPEMCC 6   s2  V   a   UPDATEMCCPARTICLEONE%SO+MODULETYPEMCC 6   É2  T   a   UPDATEMCCPARTICLEONE%GO+MODULETYPEMCC =   3  j   a   MCCPARTICLEONE%VELOCITYUPDATER+MODULETYPEMCC ;   3  _       UPDATEVELOCITYMCCPARTICLEONE+MODULETYPEMCC @   æ3  \   a   UPDATEVELOCITYMCCPARTICLEONE%MCPO+MODULETYPEMCC C   B4  @   a   UPDATEVELOCITYMCCPARTICLEONE%VFACTOR+MODULETYPEMCC (   4  é       SPECYONE+MODULESPECYONE -   k5  P   a   SPECYONE%NAME+MODULESPECYONE 3   »5  ¤   a   SPECYONE%SPECYINDEX+MODULESPECYONE 1   _6  ¤   a   SPECYONE%GASINDEX+MODULESPECYONE /   7  H   a   SPECYONE%CHARGE+MODULESPECYONE -   K7  H   a   SPECYONE%MASS+MODULESPECYONE /   7  H   a   SPECYONE%RADIUS+MODULESPECYONE .   Û7  H   a   SPECYONE%NATOM+MODULESPECYONE 4   #8  H   a   SPECYONE%INITDENSITY+MODULESPECYONE 0   k8  H   a   SPECYONE%DENSITY+MODULESPECYONE 8   ³8  H   a   SPECYONE%INITTEMPERATURE+MODULESPECYONE 4   û8  H   a   SPECYONE%TEMPERATURE+MODULESPECYONE &   C9  ñ       GASONE+MODULESPECYONE +   4:  P   a   GASONE%NAME+MODULESPECYONE .   :  ¤   a   GASONE%MCMODEL+MODULESPECYONE )   (;  H   a   GASONE%NS+MODULESPECYONE 1   p;  H   a   GASONE%INDEXSTART+MODULESPECYONE +   ¸;  H   a   GASONE%MASS+MODULESPECYONE -    <  H   a   GASONE%RADIUS+MODULESPECYONE ,   H<  H   a   GASONE%NATOM+MODULESPECYONE .   <  H   a   GASONE%BETAMAX+MODULESPECYONE 2   Ø<  H   a   GASONE%INITDENSITY+MODULESPECYONE .    =  H   a   GASONE%DENSITY+MODULESPECYONE 6   h=  H   a   GASONE%INITTEMPERATURE+MODULESPECYONE 2   °=  H   a   GASONE%TEMPERATURE+MODULESPECYONE *   ø=         REACTIONONE+MODULETYPEMCC 3   >  H   a   REACTIONONE%REACTANT+MODULETYPEMCC 7   Î>  H   a   REACTIONONE%REACTIONTYPE+MODULETYPEMCC 4   ?  H   a   REACTIONONE%RESULTANT+MODULETYPEMCC 4   ^?  H   a   REACTIONONE%THRESHOLD+MODULETYPEMCC %   ¦?         FIELDONE+MODULEFIELD (   &@  §   a   FIELDONE%NX+MODULEFIELD (   Í@  ¤   a   FIELDONE%DX+MODULEFIELD (   qA  ª   a   FIELDONE%DT+MODULEFIELD ,   B  ¬   a   FIELDONE%RHOONE+MODULEFIELD ,   ÇB  ¬   a   FIELDONE%CHIONE+MODULEFIELD "   sC  Ç       FIELD+MODULEFIELD %   :D  §   a   FIELD%NX+MODULEFIELD %   áD  ¤   a   FIELD%DX+MODULEFIELD %   E  ª   a   FIELD%DT+MODULEFIELD %   /F    a   FIELD%EX+MODULEFIELD %   ;G    a   FIELD%EY+MODULEFIELD %   GH    a   FIELD%EZ+MODULEFIELD %   SI    a   FIELD%BX+MODULEFIELD %   _J    a   FIELD%BY+MODULEFIELD %   kK    a   FIELD%BZ+MODULEFIELD &   wL  ¬   a   FIELD%RHO+MODULEFIELD &   #M  ¬   a   FIELD%PHI+MODULEFIELD &   ÏM  ¬   a   FIELD%CHI+MODULEFIELD '   {N  W   a   FIELD%DUMP+MODULEFIELD &   ÒN  Z       DUMPFIELD+MODULEFIELD )   ,O  S   a   DUMPFIELD%FG+MODULEFIELD +   O  @   a   DUMPFIELD%MODE+MODULEFIELD '   ¿O  W   a   FIELD%LOAD+MODULEFIELD &   P  \       LOADFIELD+MODULEFIELD )   rP  S   a   LOADFIELD%FG+MODULEFIELD -   ÅP  @   a   LOADFIELD%STATUS+MODULEFIELD (   Q        MCCBUNDLE+MODULETYPEMCC .   R  H   a   MCCBUNDLE%MODEL+MODULETYPEMCC 2   NR  ¥   a   MCCBUNDLE%NREACTION+MODULETYPEMCC /   óR  ¥   a   MCCBUNDLE%NSIGMA+MODULETYPEMCC 2   S  H   a   MCCBUNDLE%ENERGYMIN+MODULETYPEMCC 7   àS  H   a   MCCBUNDLE%ENERGYINTERVAL+MODULETYPEMCC 2   (T  H   a   MCCBUNDLE%ENERGYMAX+MODULETYPEMCC 7   pT  H   a   MCCBUNDLE%COLLISIONRATIO+MODULETYPEMCC 1   ¸T  H   a   MCCBUNDLE%SIGMAMAX+MODULETYPEMCC +    U  Ì   a   MCCBUNDLE%SO+MODULETYPEMCC +   ÌU  È   a   MCCBUNDLE%GO+MODULETYPEMCC 1   V  ¥   a   MCCBUNDLE%REACTION+MODULETYPEMCC 2   9W  ¬   a   MCCBUNDLE%PROBILITY+MODULETYPEMCC -   åW  [   a   MCCBUNDLE%DUMP+MODULETYPEMCC ,   @X  Q       DUMPMCCBUNDLE+MODULETYPEMCC 0   X  W   a   DUMPMCCBUNDLE%MCB+MODULETYPEMCC 4   èX  ^       POSTCOLLISIONVELOCITY+MODULETYPEMCC 9   FY  \   a   POSTCOLLISIONVELOCITY%MCPO+MODULETYPEMCC ;   ¢Y  @      POSTCOLLISIONVELOCITY%COSKAI+MODULETYPEMCC -   âY  \       ISOTROPICCOSKAI+MCCENERGYKAI 4   >Z  @   a   ISOTROPICCOSKAI%ENERGY+MCCENERGYKAI #   ~Z  j       SELECTCOLLISIONION (   èZ  \   a   SELECTCOLLISIONION%MCPO &   D[  V   a   SELECTCOLLISIONION%SO &   [  T   a   SELECTCOLLISIONION%GO &   î[  Y   a   SELECTCOLLISIONION%RO 