function hrcs_geom,x,offaxis=offaxis,rowd=rowd,Sdelta=Sdelta,$
	tilt1=tilt1,tilt2=tilt2,tilt3=tilt3,w1=w1,w3=w3,wl=wl,wr=wr,$
	e12=e12,e21=e21,e23=e23,e32=e32, order=order,dper=dper,$
	ang2mm=ang2mm, gaps=gaps, _extra=e
;+
;function	hrcs_geom
;	incorporates the HRC-S/LETG geometry and converts to/from
;	dispersion coordinate/dispersion angle
;
;	given source location and position on detector in mm, returns
;	the appropriate wavelength, and vice versa.
;
;syntax
;	y=hrcs_geom(x,offaxis=offaxis,rowd=rowd,Sdelta=Sdelta,$
;	tilt1=tilt1,tilt2=tilt2,tilt3=tilt3,w1=w1,w3=w3,wl=wl,wr=wr,$
;	e12=e12,e21=e21,e23=e23,e32=e32, order=order,dper=dper,$
;	/ang2mm, gaps=gaps)
;
;parameters
;	x	[INPUT; required] dispersion coordinate
;		(or wavelength if ANG2MM is set)
;
;keywords
;	offaxis	[INPUT; 0 arcmin] offset of 0th order source
;		relative to nominal
;	rowd	[INPUT; 8637 mm] diameter of Rowland circle
;	Sdelta	[INPUT; 0.1 mm] offset of plate S2 from best
;		focus at nominal position
;	tilt1	[INPUT; 178.5679 deg] tilt of plate S1
;	tilt2	[INPUT; 0.0 deg] tilt of plate S2
;	tilt3	[INPUT; 1.2202 deg] tilt of plate S3
;	w1	[INPUT; 103 mm] length of plate S1
;	w3	[INPUT; 103 mm] length of plate S3
;	wl	[INPUT; 48 mm] length of plate S1 left of nominal 0th order
;	wr	[INPUT; 55 mm] length of plate S1 right of nominal 0th order
;	e12	[INPUT; 0 mm] from edge of S1 to intersection of plane of S2
;	e21	[INPUT; 0 mm] from edge of S2 to intersection of plane of S1
;	e23	[INPUT; 0 mm] from edge of S2 to intersection of plane of S3
;	e32	[INPUT; 0 mm] from edge of S3 to intersection of plane of S2
;	order	[INPUT; 1] in order to convert from angle to angstrom for
;		higher orders
;	dper	[INPUT; 9912.5 Ang] grating period
;	ang2mm	[INPUT; 0] if set, assumes that X are in wavelengths in
;		[Ang] and converts to [mm]
;	gaps	[OUTPUT] wavelengths corresponding to the edges of the plates
;		in the order [lS1,rS1, lS2,rS2, lS3,rS3]
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	file://localhost/data/toad1/kashyap/Olivia/hrcs_geom_fig.ps
;	file://localhost/data/toad1/kashyap/Olivia/hrcs_geom_dvi.ps
;
;history
;	vinay kashyap (JunMM)
;	ROWD changed (VK; FebMMI)
;-

;	Usage
nx=n_elements(x)
if nx eq 0 then begin
  print,'Usage: l=hrcs_geom(x,offaxis=offaxis,rowd=rowd,Sdelta=Sdelta,$'
  print,'       tilt1=tilt1,tilt2=tilt2,tilt3=tilt3,w1=w1,w3=w3,wl=wl,wr=wr,$'
  print,'       e12=e12,e21=e21,e23=e23,e32=e32, order=order,dper=dper,$'
  print,'       /ang2mm, gaps=gaps)'
  print,'  convert dispersion coordinate to wavelength and vice versa'
  return,-1L
endif

;	keywords
if not keyword_set(offaxis) then offaxis=0.	;[arcmin]
;if not keyword_set(rowd) then rowd=8632.15	;[mm]
;if not keyword_set(rowd) then rowd=8632.31	;[mm]
if not keyword_set(rowd) then rowd=8637.0	;[mm]	in POG 3.
if not keyword_set(Sdelta) then Sdelta=0.1	;[mm]
if not keyword_set(tilt1) then tilt1=178.5679	;[deg]
if not keyword_set(tilt2) then tilt2=0.0	;[deg]
if not keyword_set(tilt3) then tilt3=1.2202	;[deg]
if not keyword_set(w1) then w1=103.0	;[mm]
if not keyword_set(w3) then w3=103.0	;[mm]
if not keyword_set(wl) then wl=48.0	;[mm]
if not keyword_set(wr) then wr=55.0	;[mm]
if not keyword_set(e12) then e12=0.0	;[mm]
if not keyword_set(e21) then e21=0.0	;[mm]
if not keyword_set(e23) then e23=0.0	;[mm]
if not keyword_set(e32) then e32=0.0	;[mm]
if not keyword_set(ang2mm) then ang2mm=0
ord=1 & if keyword_set(order) then ord=abs(fix(order)) > 1
;if not keyword_set(pltscl) then pltscl=1.148	;[Ang/mm]
if not keyword_set(dper) then dper=9912.5	;[Ang]

;	initialize
xx=[x[*]]
R=rowd/2.	;radius of Rowland circle [mm]
;	convert all angles to [radian]
omega=(offaxis/60.)*!pi/180.
alpha1=tilt1*!pi/180. & alpha2=tilt2*!pi/180. & alpha3=tilt3*!pi/180.
;
;	some special items
U_G = (sin(omega)/cos(alpha2-omega)) * (R-Sdelta)
H_G = R + (cos(alpha2)/cos(alpha2-omega))*(R-Sdelta)
;
G_C = wl + U_G
H_C = sqrt(G_C^2+H_G^2+2.*G_C*H_G*sin(alpha2-omega))
theta_C = asin((G_C/H_C)*cos(alpha2-omega))
G_D = wr - U_G
H_D = sqrt(G_D^2+H_G^2+2.*G_D*H_G*sin(alpha2-omega))
theta_D = asin((G_D/H_D)*cos(alpha2-omega))
;
G_M = e21 + G_C
H_M = sqrt(G_M^2+H_G^2+2.*G_M*H_G*sin(alpha2-omega))
theta_M = asin((G_M/H_M)*cos(alpha2-omega))
G_N = e23 + G_D
H_N = sqrt(G_N^2+H_G^2+2.*G_N*H_G*sin(alpha2-omega))
theta_N = asin((G_N/H_N)*cos(alpha2-omega))
;
H_E = sqrt(H_N^2+e32^2-2*e32*H_N*sin(alpha3-omega-theta_N))
theta_E = theta_N + asin((e32/H_E)*cos(alpha3-omega-theta_N))
H_F = sqrt(H_N^2+(w3+e32)^2-2*(w3+e32)*H_N*sin(alpha3-omega-theta_N))
theta_F = theta_N + asin(((w3+e32)/H_F)*cos(alpha3-omega-theta_N))
;
H_B = sqrt(H_M^2+e12^2-2*e12*H_M*sin(alpha1-omega+theta_M))
theta_B = theta_M - asin((e12/H_B)*cos(alpha1-omega+theta_M))
H_A = sqrt(H_M^2+(w1+e12)^2-2*(w1+e12)*H_M*sin(alpha1-omega+theta_M))
theta_A = theta_M - asin(((w1+e12)/H_A)*cos(alpha1-omega+theta_M))

;	get plate edge wavelengths
theta_A=-theta_A & theta_B=-theta_B & theta_C=-theta_C & theta_M=-theta_M
theta=[theta_A, theta_B, theta_C, theta_D, theta_E, theta_F]
nlam=dper*sin(theta) & gaps=nlam/ord

;	convert from wavelength to dispersion coordinate if ANG2MM is set
if keyword_set(ang2mm) then begin
  th=asin(ord*xx/dper) & y=0*th
  ;o1=where(th ge theta_A and th lt theta_M,mo1)
  ;o3=where(th ge theta_N and th lt theta_F,mo3)
  o1=where(th lt theta_M,mo1)
  o2=where(th ge theta_M and th lt theta_N,mo2)
  o3=where(th ge theta_N,mo3)
  if mo1 gt 0 then begin
    theta=abs(th[o1])
    G_M=U_G+wl+e21
    G_Tpp=(sin(theta)/cos(alpha2-omega+theta))*H_G
    M_Tpp=G_Tpp-G_M
    z_pp=-(cos(alpha2-omega+theta)/cos(alpha1-omega+theta)) * M_Tpp
    y[o1]=-(z_pp[*]+G_M)
  endif
  if mo2 gt 0 then begin
    theta=th[o2]
    z=(sin(theta)/cos(alpha2-omega-theta))*H_G
    y[o2]=z[*]
  endif
  if mo3 gt 0 then begin
    theta=th[o3]
    G_N=U_G+wr+e23
    G_Tp=(sin(theta)/cos(alpha2-omega-theta))*H_G
    N_Tp=G_Tp-G_N
    z_p=(cos(alpha2-omega-theta)/cos(alpha3-omega-theta)) * N_Tp
    y[o3]=z_p[*]+G_N
  endif
endif

;	convert from dispersion coordinate to wavelength in ANG2MM is not set
if not keyword_set(ang2mm) then begin
  y=0*xx
  o1=where(xx lt -G_M,mo1)
  o2=where(xx ge -G_M and xx lt G_N,mo2)
  o3=where(xx ge G_N,mo3)
  if mo1 gt 0 then begin
    M_Spp = abs(xx[o1] + G_M)
    H_Spp = sqrt(H_M^2+M_Spp^2-2*M_Spp*H_M*sin(alpha1-omega+abs(theta_M)))
    theta = abs(theta_M) - asin((M_Spp/H_Spp)*cos(alpha1-omega+abs(theta_M)))
    y[o1]=dper*sin(-theta[*])/ord
  endif
  if mo2 gt 0 then begin
    G_S = xx[o2]
    H_S = sqrt(G_S^2+H_G^2-2*G_S*H_G*sin(alpha2-omega))
    theta=asin((G_S/H_S)*cos(alpha2-omega))
    y[o2]=dper*sin(theta[*])/ord
  endif
  if mo3 gt 0 then begin
    N_Sp = abs(xx[o3] - G_N)
    H_Sp = sqrt(H_N^2+N_Sp^2-2*N_Sp*H_N*sin(alpha3-omega-theta_N))
    theta = theta_N + asin((N_Sp/H_Sp)*cos(alpha3-omega-theta_N))
    y[o3]=dper*sin(theta[*])/ord
  endif
endif

return,y
end
