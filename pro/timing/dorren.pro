pro dorren, x, A, f, pder,ldstr=ldstr,ldspt=ldspt,rdns=rdns,$
photT=photT,spotT=spotT,wvl=wvl,_extra=e

;+
;function  dorren 
;
; computes spot modulated stellar light curve using Dorren(1987)
;
;syntax 
;      curve =  x, A, f, pder,ldstr=ldstr,ldspt=ldspt,rdns=rdns,$
;      photT=photT,spotT=spotT,wvl=wvl
; 
;parameters 
;       x   [INPUT;required] absissca points at which to compute spot
;           model. These are angular displacements with respect to the 
;           initial spot longitude specified in A, in units of degrees. 
;       A   [INPUT;required] spot parameters for dorrens model in
;                            radians where: 
;                            a(0) = stellar inclination to LOS
;                            a(1) = spot latitude 
;                            a(2) = spot size 
;                            a(3) = initial spot longitude 
;       pder[OUTPUT]         partial derivative with respect to each 
;                            parameter in A 
;keywords 
;       ldstr [INPUT] limb-darkening coefficient for star 
;                          default = 0.32 Van Hamme(1994) 
;       ldspt [INPUT] limb-darkening coefficient for spot (umbra) 
;       rdns  [INPUT] set this parameter if spot parameters and absissca are input 
;                     in radians rather than degrees
;       photT [INPUT] photospheric temperature in degrees Kelvin
;       spotT [INPUT] spot temperature in degrees Kelvin 
;       wvl   [INPUT] wavelength in nanometers
;subroutines 
;       MOD
;       
;history
;       liwei lin (Sep 03) better documentation/functionality than
;                          mod.pro. switch to modular
;                 (Jul 06) removed loop over longitudes and add pder
;                 (Aug 06) BUGFIX: note that T should be defined as: 
;                           Pi - ATAN((-SIN(h(tmp)))*(TAN(b1(tmp)))) when B1>=Pi/2           
;                           ATAN((SIN(h(tmp)))*(TAN(b1(tmp)))) when B1<Pi/2 
;                          to avoid dicontinuities. This differs from Dorren (1987)
;                          formalism only in the '=' sign placement 
;-            

nx = n_elements(x) & f  = fltarr(nx)  & pi = !pi 
if not keyword_set(ldstr) then ldstr = 0.306D ;limb dark coeiff star  
if not keyword_set(ldspt) then ldspt = 0.36D  ;limb dark coeff spot using umbra
if not keyword_set(photT) then phott = 4000D  ;photospheric temperature 
if not keyword_set(spotT) then spott = 3300D  ;spot temperature 
if not keyword_set(wvl)   then wvl   = 809    ;wavelength 
if not keyword_set(radns) then begin 
    a = !dtor*a & x = !dtor*x
endif
; cgs constants 
c = 2.9979246e+10
h = 6.6261760e-27 
k = 1.3806620e-16

; spot flux / photospheric flux  
ff = exp(((h*c)/(k*wvl*1d-7))*((1/double(photT))-(1/double(spotT))))
  ; longitudinal displacement with respect to chosen reference
  lng = double(x - a(3))

  ; angular distance between LOS and surface normal at spot center
  b= ACOS((COS(a(0))*SIN(a(1)))+(SIN(a(0))*COS(a(1))*COS(Lng)))
 
  ; initialize AAA and BBB
  AAA = dblarr(nx) & BBB = dblarr(nx) ; T = dblarr(nx)

  ;1. the spot is partially out of view
      oo1 = where((((pi/2.)-a(2)) lt b) and (((pi/2.)+a(2)) gt b),noo1)
      if oo1(0) gt -1 then begin
          b1 = b(oo1) 
          h  = ACOS(COS(a(2))/SIN(b1))
          d  = ACOS((COS(a(2))/SIN(a(2)))*(COS(b1)/SIN(b1)))	
          T  = findgen(noo1) 
          ;if (b1 gt Pi/2) then begin 
              tmp = where(b1 ge Pi/2) 
              if tmp(0) gt -1 then T(tmp)   = Pi - ATAN((-SIN(h(tmp)))*(TAN(b1(tmp))))           
              tmp = where(b1 lt Pi/2) 
              if tmp(0) gt -1 then T(tmp)   = ATAN((SIN(h(tmp)))*(TAN(b1(tmp))))
          AAA(oo1) = h + (Pi - d)*COS(b1)*(SIN(a(2)))^2 - SIN(h)*SIN(b1)*COS(a(2))
          BBB(oo1) = ((1/3.)*(Pi - d)*((-2.0*(COS(a(2)))^3)  - $
                (3.0*((SIN(b1))^2)*COS(a(2))*((SIN(a(2)))^2))))+ $
                ((2/3.)*(Pi-T))+((1/6.)*(SIN(h)*SIN(2*b1)*(2-3*(COS(a(2)))^2)))
      endif

  ;2. the spot is completely in view
      oo2 = where(((b+a(2)) LE Pi/2.))
      ;b2 = b(oo2) 
      if oo2(0) gt -1 then begin
      b2 = b(oo2) 
          d=0 & h=0 & T=0
          AAA(oo2) = h + (Pi - d)*COS(b2)*(SIN(a(2)))^2 - SIN(h)*SIN(b2)*COS(a(2))
          BBB(oo2) = ((1/3.)*(Pi - d)*((-2.0*(COS(a(2)))^3)- $
            (3.0*((SIN(b2))^2)*COS(a(2))*((SIN(a(2)))^2))))+ $ 
            ((2/3.)*(Pi-T))+((1/6.)*(SIN(h)*SIN(2*b2)*(2-3*(COS(a(2)))^2)))
      endif

  ;3.  the spot is completely out of view
      oo3 = where(b-(a(2)) ge Pi/2) 
      ;b3 = b(oo3)
      if oo3(0) gt -1 then begin
          b3 = b(oo3)
          AAA(oo3) = 0
          BBB(oo3) = 0
      endif
  aa = (1 - ldstr) - (1 - ldspt)*ff
  bb = ldstr - ldspt*(0.23D)
  diffmag =Double (-2.5*alog10(1-(aa*AAA+bb*BBB)/(Pi*(1-(ldstr)/3.))) )
  f = diffmag
  pder= fltarr(nx,4) 

  ;1. the spot is partially out of view

  if oo1(0) gt -1 then begin
       tmp = where(b1 ge Pi/2 )
       ; note: no difference between T's for PDERs 
       h  = ACOS(COS(a(2))/SIN(b1))
       d  = ACOS((COS(a(2))/SIN(a(2)))*(COS(b1)/SIN(b1)))	
        pder(oo1,0) =  (0.345600567762181*((1 - 0.23*D*(1 - ldspt) - ldstr)* $
    (-((Cos(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
        (Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1))))/  $
       ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)* $
        Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
             2)))) + (AAA*Cos(a(2))^3*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
        Cos(a(0))*Sin(a(1)))*(Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1))))/ $
      ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))) + (AAA*Cos(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
       (Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1)))* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2)))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2) +  $
     (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))* $
      (Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1)))*Sin(a(2))^2 -  $
     ((Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
       (-((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2* $
           (Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1))))/ $
          (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)) -  $
        (1/Tan(a(2))*(Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))*Sin(a(2))^2)/ $
      Sqrt(1 - (1/Tan(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)/ $
         (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))) +  $
   (-0.23*D*ldspt + ldstr)*((0.3333333333333333*(2 - 3*Cos(a(2))^2)* $
       Cos(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))* $
       (-(Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1))) + Sin(a(0))*Sin(a(1)))* $ 
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^  $
            2)))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2) -  $
     (0.6666666666666666*(-((Cos(a(2))^2*(Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) -  $
            Sin(a(0))*Sin(a(1))))/((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))* $
                Sin(a(1)))^2)^(3/2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))* $
                  Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))) -  $
        ((Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1)))* $
          Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
               2)))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2) - ((Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1)))* $
          Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
          Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
               2)))/(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))/ $
      (1 + ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
         (1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))/ $
        (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2) +  $
     2.*(Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))*Cos(a(2))* $
      (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
      (Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1)))*Sin(a(2))^2 -  $
     (0.3333333333333333* $
       (-((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2* $
           (Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1))))/ $
          (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)) -  $
        (1/Tan(a(2))*(Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))*  $ 
       (-2.*Cos(a(2))^3 - 3.*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $  
            Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2))/ $
      Sqrt(1 - (1/Tan(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)/ $
         (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)) -  $
     (0.16666666666666666*Cos(a(2))^2*(2 - 3*Cos(a(2))^2)* $
       (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
       (Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo1)) - Sin(a(0))*Sin(a(1)))* $
       Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))))/ $
      ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^2* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))))))/((1 - ldstr/3)* $
  (1 - ((1 - 0.23*D*(1 - ldspt) - ldstr)* $
      (ACos(Cos(a(2))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2)) - AAA*Cos(a(2))*Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1)))^2)* $
        Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
             2)) + (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1))))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))* $
                Sin(a(1)))^2)))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
        Sin(a(2))^2) + (-0.23*D*ldspt + ldstr)* $
      (-0.6666666666666666*ATan( $
         (Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
           Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
                 Cos(a(0))*Sin(a(1)))^2)))/(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
           Cos(a(0))*Sin(a(1)))) + 0.3333333333333333* $
        (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
           Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))* $
        (-2.*Cos(a(2))^3 - 3.*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2) + 0.16666666666666666* $
        (2 - 3*Cos(a(2))^2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
              Cos(a(0))*Sin(a(1)))^2))*Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1))))))/((1 - ldstr/3)*Pi))) 

    pder(oo1,1) =  (0.345600567762181*((1 - 0.23*D*(1 - ldspt) - ldstr)* $
    (-((Cos(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
        (Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1))))/ $
       ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)* $
        Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
             2)))) + (AAA*Cos(a(2))^3*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
        Cos(a(0))*Sin(a(1)))*(Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1))))/ $
      ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))) + (AAA*Cos(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
       (Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1)))* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2)))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2) +  $
     (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))* $
      (Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1)))*Sin(a(2))^2 -  $
     ((Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
       (-((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2* $
           (Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1))))/ $
          (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)) -  $
        (1/Tan(a(2))*(Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))*Sin(a(2))^2) / $
      Sqrt(1 - (1/Tan(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)/ $
         (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))) +  $
   (-0.23*D*ldspt + ldstr)*((0.3333333333333333*(2 - 3*Cos(a(2))^2)* $
       Cos(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))* $
       (-(Cos(a(0))*Cos(a(1))) + Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1)))* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2)))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2) -  $
     (0.6666666666666666*(-((Cos(a(2))^2*(Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))* $
             Sin(a(1))))/((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^ $
            (3/2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
                 Cos(a(0))*Sin(a(1)))^2)))) -  $
        ((Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1)))* $
          Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
               2)))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2) - ((Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1)))* $
          Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
          Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
               2)))/(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))/ $
      (1 + ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
         (1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))/ $
        (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2) +  $
     2.*(Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))*Cos(a(2))* $
      (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))*(Cos(a(0))*Cos(a(1)) -  $
       Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1)))*Sin(a(2))^2 -  $
     (0.3333333333333333* $
       (-((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2* $
           (Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1))))/ $
          (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)) -  $
        (1/Tan(a(2))*(Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))* $
       (-2.*Cos(a(2))^3 - 3.*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2))/ $
      Sqrt(1 - (1/Tan(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)/ $
         (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)) -  $
     (0.16666666666666666*Cos(a(2))^2*(2 - 3*Cos(a(2))^2)* $
       (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))*(Cos(a(0))*Cos(a(1)) -  $
        Cos((a(3)) - x(oo1))*Sin(a(0))*Sin(a(1)))* $
       Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))))/ $ 
      ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^2* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))))))/((1 - ldstr/3)* $
  (1 - ((1 - 0.23*D*(1 - ldspt) - ldstr)* $
      (ACos(Cos(a(2))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2)) - AAA*Cos(a(2))*Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1)))^2)* $
        Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
             2)) + (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1))))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))* $
                Sin(a(1)))^2)))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
        Sin(a(2))^2) + (-0.23*D*ldspt + ldstr)* $
      (-0.6666666666666666*ATan( $
         (Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
           Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
                 Cos(a(0))*Sin(a(1)))^2)))/(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
           Cos(a(0))*Sin(a(1)))) + 0.3333333333333333* $
        (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
           Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))* $
        (-2.*Cos(a(2))^3 - 3.*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2) + 0.16666666666666666* $
        (2 - 3*Cos(a(2))^2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
              Cos(a(0))*Sin(a(1)))^2))*Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1))))))/((1 - ldstr/3)*Pi))) 

    pder(oo1,2) =   (0.345600567762181*((1 - 0.23*D*(1 - ldspt) - ldstr)* $
    (-((Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2/ $
       (Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
        Sqrt(1 - (1/Tan(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)/ $
           (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))) +  $
     2*(Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))*Cos(a(2))* $
      (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))*Sin(a(2)) +  $
     Sin(a(2))/(Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))) - (AAA*Cos(a(2))^2*Sin(a(2)))/ $
      (Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))) + AAA*Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
         2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1)))^2))*Sin(a(2))) + (-0.23*D*ldspt + ldstr)* $
    ((-0.6666666666666666*Cos(a(2))*Sin(a(2)))/((Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
        Cos(a(0))*Sin(a(1)))*Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
          2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2))* $
       (1 + ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
          (1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
              2)))/(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)) -  $
     (0.3333333333333333*1/Sin(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
        Cos(a(0))*Sin(a(1)))*(-2.*Cos(a(2))^3 - 3.*Cos(a(2))* $
         (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2))/ $
      (Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
       Sqrt(1 - (1/Tan(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)/ $
          (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))) +  $
     0.3333333333333333*(Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
           Cos(a(0))*Sin(a(1))))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2)))*(6.*Cos(a(2))^2*Sin(a(2)) -  $
       6.*Cos(a(2))^2*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
        Sin(a(2)) + 3.*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
        Sin(a(2))^3) + (0.16666666666666666*Cos(a(2))*(2 - 3*Cos(a(2))^2)*Sin(a(2))* $
       Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))))/ $ 
      ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))) + 1.*Cos(a(2))* $
      Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))* $
      Sin(a(2))*Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))))))/ $
 ((1 - ldstr/3)*(1 - ((1 - 0.23*D*(1 - ldspt) - ldstr)* $
      (ACos(Cos(a(2))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2)) - AAA*Cos(a(2))*Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1)))^2)* $
        Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
             2)) + (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1))))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))* $
                Sin(a(1)))^2)))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
        Sin(a(2))^2) + (-0.23*D*ldspt + ldstr)* $
      (-0.6666666666666666*ATan( $
         (Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
           Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
                 Cos(a(0))*Sin(a(1)))^2)))/(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
           Cos(a(0))*Sin(a(1)))) + 0.3333333333333333* $
        (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
           Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))* $
        (-2.*Cos(a(2))^3 - 3.*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2) + 0.16666666666666666* $
        (2 - 3*Cos(a(2))^2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
              Cos(a(0))*Sin(a(1)))^2))*Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1))))))/((1 - ldstr/3)*Pi))) 

    pder(oo1,3) = (0.345600567762181*((1 - 0.23*D*(1 - ldspt) - ldstr)* $
    ((Cos(a(1))*Cos(a(2))*Sin(a(0))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
       Sin((a(3)) - x(oo1)))/((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^ $
        (3/2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2))) - (AAA*Cos(a(1))*Cos(a(2))^3*Sin(a(0))* $
       (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))*Sin((a(3)) - x(oo1)))/ $
      ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))) - (AAA*Cos(a(1))*Cos(a(2))*Sin(a(0))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
        Cos(a(0))*Sin(a(1)))*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2))*Sin((a(3)) - x(oo1)))/ $
      Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2) -  $
     (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))*Cos(a(1))* $
      Sin(a(0))*Sin(a(2))^2*Sin((a(3)) - x(oo1)) -  $
     ((Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))*Sin(a(2))^2* $
       ((Cos(a(1))*1/Tan(a(2))*Sin(a(0))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2* $
          Sin((a(3)) - x(oo1)))/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^ $
          (3/2) + (Cos(a(1))*1/Tan(a(2))*Sin(a(0))*Sin((a(3)) - x(oo1)))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))/ $
      Sqrt(1 - (1/Tan(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)/ $
         (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2))) +  $
   (-0.23*D*ldspt + ldstr)*((0.3333333333333333*Cos(a(1))*(2 - 3*Cos(a(2))^2)* $
       Cos(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))*Sin(a(0))* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2))*Sin((a(3)) - x(oo1)))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
          Cos(a(0))*Sin(a(1)))^2) -  $
     2.*(Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))*Cos(a(1))* $
      Cos(a(2))*Sin(a(0))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))*Sin(a(2))^2* $
      Sin((a(3)) - x(oo1)) - (0.3333333333333333*(-2.*Cos(a(2))^3 -  $
        3.*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
         Sin(a(2))^2)*((Cos(a(1))*1/Tan(a(2))*Sin(a(0))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1)))^2*Sin((a(3)) - x(oo1)))/ $
         (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2) +  $
        (Cos(a(1))*1/Tan(a(2))*Sin(a(0))*Sin((a(3)) - x(oo1)))/ $
         Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))/ $
      Sqrt(1 - (1/Tan(a(2))^2*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)/ $
         (1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)) -  $
     (0.6666666666666666*((Cos(a(1))*Cos(a(2))^2*Sin(a(0))*Sin((a(3)) - x(oo1)))/ $
         ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^(3/2)* $
          Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
               2))) + (Cos(a(1))*Sin(a(0))* $
          Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
               2))*Sin((a(3)) - x(oo1)))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2) + (Cos(a(1))*Sin(a(0))* $
          Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
          Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
               2))*Sin((a(3)) - x(oo1)))/(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
          2))/(1 + ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
         (1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))/ $
        (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2) +  $
     (0.16666666666666666*Cos(a(1))*Cos(a(2))^2*(2 - 3*Cos(a(2))^2)*Sin(a(0))* $
       (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))*Sin((a(3)) - x(oo1))* $
       Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))))/ $
      ((1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)^2* $
       Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^  $
            2))))))/((1 - ldstr/3)* $ 
  (1 - ((1 - 0.23*D*(1 - ldspt) - ldstr)* $
      (ACos(Cos(a(2))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
            2)) - AAA*Cos(a(2))*Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1)))^2)* $ 
        Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^ $
             2)) + (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1))))/Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))* $
                Sin(a(1)))^2)))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
        Sin(a(2))^2) + (-0.23*D*ldspt + ldstr)* $
      (-0.6666666666666666*ATan( $
         (Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)* $
           Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
                 Cos(a(0))*Sin(a(1)))^2)))/(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
           Cos(a(0))*Sin(a(1)))) + 0.3333333333333333* $
        (Pi - ACos((1/Tan(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1))))/ $
           Sqrt(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)))* $
        (-2.*Cos(a(2))^3 - 3.*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
             Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2) + 0.16666666666666666* $
        (2 - 3*Cos(a(2))^2)*Sqrt(1 - Cos(a(2))^2/(1 - (Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
              Cos(a(0))*Sin(a(1)))^2))*Sin(2*ACos(Cos(a(1))*Cos((a(3)) - x(oo1))*Sin(a(0)) +  $
            Cos(a(0))*Sin(a(1))))))/((1 - ldstr/3)*Pi)))
endif 
if oo2(0) gt -1 then begin
        d=0 & h=0 & T=0
    pder(oo2,0) =  (0.345600567762181*(Pi*aa*(Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo2)) - Sin(a(0))*Sin(a(1)))* $
            Sin(a(2))^2 + 2*(-0.23*D*ldspt + ldstr)*Pi*Cos(a(2))* $
            (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
            (Cos(a(0))*Cos(a(1))*Cos((a(3)) - x(oo2)) - Sin(a(0))*Sin(a(1)))*Sin(a(2))^2))/ $
            ((1 - ldstr/3)*(1 - (Pi*aa*(Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
            Sin(a(2))^2 + (-0.23*D*ldspt + ldstr)*((2*Pi)/3 + $
            (Pi*(-2*Cos(a(2))^3 - 3*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + $ 
            Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2))/3))/((1 - ldstr/3)*Pi))) 
    pder(oo2,1) = (0.345600567762181*(Pi*aa*(Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo2))*Sin(a(0))*Sin(a(1)))* $
            Sin(a(2))^2 + 2*(-0.23*D*ldspt + ldstr)*Pi*Cos(a(2))* $
            (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
            (Cos(a(0))*Cos(a(1)) - Cos((a(3)) - x(oo2))*Sin(a(0))*Sin(a(1)))*Sin(a(2))^2))/ $ 
            ((1 - ldstr/3)*(1 - (Pi*aa*(Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
            Sin(a(2))^2 + (-0.23*D*ldspt + ldstr)*((2*Pi)/3 + $
            (Pi*(-2*Cos(a(2))^3 - 3*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + $
            Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2))/3))/((1 - ldstr/3)*Pi))) 
    pder(oo2,2) =  (0.345600567762181*(2*Pi*aa*Cos(a(2))*(Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + $
            Cos(a(0))*Sin(a(1)))*Sin(a(2)) + ((-0.23*D*ldspt + ldstr)*Pi* $
            (6*Cos(a(2))^2*Sin(a(2)) - 6*Cos(a(2))^2*(1 - (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + $
            Cos(a(0))*Sin(a(1)))^2)*Sin(a(2)) + $
            3*(1 - (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^3))/3))/ $
            ((1 - ldstr/3)*(1 - (Pi*aa*(Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
            Sin(a(2))^2 + (-0.23*D*ldspt + ldstr)*((2*Pi)/3 + $
            (Pi*(-2*Cos(a(2))^3 - 3*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + $
            Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2))/3))/((1 - ldstr/3)*Pi))) 
    pder(oo2,3) =   (0.345600567762181*(-(Pi*aa*Cos(a(1))*Sin(a(0))*Sin(a(2))^2*Sin((a(3)) - x(oo2))) - $
            2*(-0.23*D*ldspt + ldstr)*Pi*Cos(a(1))*Cos(a(2))*Sin(a(0))* $
            (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))*Sin(a(2))^2*Sin((a(3)) - x(oo2))))/ $
            ((1 - ldstr/3)*(1 - (Pi*aa*(Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + Cos(a(0))*Sin(a(1)))* $
            Sin(a(2))^2 + (-0.23*D*ldspt + ldstr)*((2*Pi)/3 + $
            (Pi*(-2*Cos(a(2))^3 - 3*Cos(a(2))*(1 - (Cos(a(1))*Cos((a(3)) - x(oo2))*Sin(a(0)) + $
            Cos(a(0))*Sin(a(1)))^2)*Sin(a(2))^2))/3))/((1 - ldstr/3)*Pi))) 
endif
if oo3(0) gt -1 then begin
    pder(oo3,0) = 0
    pder(oo3,1) = 0
    pder(oo3,2) = 0
    pder(oo3,3) = 0
endif
if not keyword_set(radn) then begin 
    a = !radeg*a & x = !radeg*x
endif
f= (f-mean(f))
end

