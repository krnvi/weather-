load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

  in     = addfile("/user1/wrfarw/WRFV3/run/wrfout_d04_2009-05-01_00:00:00.nc","r")



times   = wrf_user_list_times(in)         
ntimes  = dimsizes(times)  
          
TH      = new(ntimes,float)
T       = new(ntimes,float)               
P       = new(ntimes,float)
PB      = new(ntimes,float)      
T2      = new(ntimes,float)
TH2     = new(ntimes,float)
PSFC    = new(ntimes,float)
P_l     = new(ntimes,float)
RH      = new(ntimes,float)
Q       = new(ntimes,float)
UST     = new(ntimes,float)
QV      = new(ntimes,float)
wind_10m = new(ntimes,float)
wind_speed= new(ntimes,float)
THGT    = new(ntimes,float)
HGT     = new(ntimes,float)
rh      = new(ntimes,float)
UC      = new(ntimes,float)
VC      = new(ntimes,float)
PH      = new(ntimes,float)               
PHB     = new(ntimes,float)               
ghgt    = new(ntimes,float)
Gl      = new(ntimes,float)

;j= in->ZNW
;k= sizeof(j)
;Levels=k-1
Levels=27


do il = 0,Levels-1 ;Levels-1
do it = 0,ntimes-1 ;ntimes-1                  ;Loop for the time: it= starting time
time = it

Latitude = 26.799
Longitude = 71.066

res = True
res@returnInt = True                       ; False : return real values, True: return interger values
point = wrf_user_ll_to_ij(in,Longitude,Latitude,res)       ; wrf_user_ll_to_ij(nc_file,lon,lat,opt)

x = point(0)
y = point(1)

;print("X location is: " + x)            ; print the value of X at the screen
;print("Y location is: " + y)            ; print the value of Y at the screen



;wind spd
  U     = wrf_user_getvar(in, "U", time)
  ua_in = wrf_user_unstagger(U,U@stagger)
  ua    = ua_in(il:il,x,y)             ;ua(bottom_up,grid_lon,grid_lat)
  UC(it) = ua 

  V     = wrf_user_getvar(in, "V", time)
  va_in = wrf_user_unstagger(V,V@stagger)
  va    = va_in(il:il,x,y)            ;va(bottom_up,grid_lon,grid_lat)
  VC(it)= va  

  wind_speed(it) = sqrt(ua^2+va^2)
  copy_VarCoords(ua,wind_speed(it))                         ; copy coord vars to speed
 wind_speed@Units= "m/s"
; print(wind_speed(it))

;terrain hgt

 thgt     = wrf_user_getvar(in,"HGT",it)
 THGT(it) = thgt(x,y)

;level height

  hgt       = wrf_user_getvar(in,"height",it) ;
  kht       = hgt(il:il,x,y)
  HGT(it)   = kht

;pert. temp at eta level

  TH_in     = wrf_user_getvar(in, "T", time)
  th        = TH_in(il:il,x,y)
  TH(it)    = th 

; Pot. temp @ eta level

  T_in     = wrf_user_getvar(in, "tk", time)
  t        = T_in(il:il,x,y)             
  T(it)    = t

;print(T)   

;Qvapor
 QV_in     = wrf_user_getvar(in, "QVAPOR", time)
 qv        = QV_in(il:il,x,y)
 QV(it)    = qv


;pert pressure @ eta level

  P_in     = wrf_user_getvar(in, "P", time)
  p        = P_in(il:il,x,y) 
  P(it)    = p           
 ;print(P)

; base state pressure @ eta level

  PB_in     = wrf_user_getvar(in, "PB", time)
  pb        = PB_in(il:il,x,y) 
  PB(it)    = pb           
 ;print(PB)

;total pressure @ eta level 
 
  P_l(it)       = (P(it)+PB(it))*0.01

  rh_in         = wrf_user_getvar(in, "rh", time)
  rhl           = rh_in(il:il,x,y)
  rh(it)        = rhl
  
  Ph_in     = wrf_user_getvar(in, "PH", time)
  ph        = Ph_in(il:il,x,y)
  PH(it)    = ph

  PHB_in     = wrf_user_getvar(in, "PHB", time)
  phb        = PHB_in(il:il,x,y)
  PHB(it)    = phb

  ghgt(it)     =  (PHB(it)+PH(it))/9.81
; Ghgt(it) = ghgt

   gl       = ghgt(it)-THGT(it)
   Gl(it)   = gl

end do

print("  Time  level     Ter.hgt        Lv.hgt         Ucmpnt        Vcmpnt        wspd         Pert.Temp     pot.temp      pert.pres     base.pres       pres@level        WVMR         Rel.Hmdt       Pgeop       bgeop    geop.hgt   agl")

do it = 0,ntimes-1 ;ntimes-1

    print (sprintf("%5.0f",it)    +" " \
               +sprintf("%5.0f",il)    +" " \         
                           +sprintf("%13.3f",THGT(it))    +" " \
                                  +sprintf("%13.3f",HGT(it))    +" " \
                                          +sprintf("%13.3f",UC(it))    +" " \
                                                   +sprintf("%13.3f",VC(it))    +" " \
                                                      +sprintf("%13.3f",wind_speed(it))    +" " \
                                                                    + sprintf("%13.3f", TH(it)) +"  " \
                                                                           + sprintf("%13.3f", T(it)) +"  " \
                                                                                   +sprintf("%13.3f", P(it)) +"  " \
                                                                                           +sprintf("%13.3f", PB(it)) +"  " \
                                                                                                   +sprintf("%13.3f", P_l(it)) +"  " \
                                                                                                            +sprintf("%13.3f", QV(it)) +"  " \
                                                                                                                    +sprintf("%13.3f", rh(it)) +"  " \
                                                                                                                            + sprintf("%13.2f", PH(it)) +"  " \  
                                                                                                                                     + sprintf("%13.2f", PB(it)) +"  " \
                                                                                                                                                + sprintf("%13.2f", ghgt(it)) +"  " \
                                                                                                                                                             + sprintf("%13.2f", Gl(it)) +"  " )   


end do

end do

do it= 0,ntimes-1 ;ntimes-1
time= it

;surface pressure               

PSFC_in = wrf_user_getvar(in, "PSFC", time)   ; Extract surface temperature (Pa)
PSFC(it) = PSFC_in(x,y)
P_rh= PSFC_in(x,y)

;print(P(it))

; Temperature at 2m

T2_in = wrf_user_getvar(in, "T2", time)      ; Extract Temperature (K) at 2m
T2(it) = T2_in(x,y)            ; Temperature at 2m in Kelvin
T2_rh = T2_in(x,y)

;print(T(it))

;potential Temp at 2m

th2 = wrf_user_getvar(in,"TH2",it) ; Get U10
TH2(it)    = th2(x,y)
;print(TH2(it))

Q_in = wrf_user_getvar(in, "Q2", time)      ; Extract Mixing ratio at 2m
q    = Q_in(x,y)
Q(it) = q
Q_rh=Q_in(x,y)

RH(it) = wrf_rh(Q_rh, P_rh, T2_rh)            ; Calculate relative humidity
;print(RH(it))

;wind @ 10m

U10 = wrf_user_getvar(in,"U10",it) ; Get U10
u10    = U10(x,y)

V10 = wrf_user_getvar(in,"V10",it) ; get V10
v10    = V10(x,y)
wind_10m(it) = sqrt(u10^2+v10^2)
wind_10m@Units= "m/s"

;fric vel
ust= wrf_user_getvar(in,"UST",it)
UST(it)= ust(x,y)

end do

print("  Time       temp@2m      pot.temp@2m     sur.pres          mixR            RH           wind10         fric.vel")

do it = 0,ntimes-1 ;ntimes-1

    print (sprintf("%5.0f",it)    +" " \
         + sprintf("%13.2f", T2(it)) +"  " \
                   + sprintf("%13.2f", TH2(it)) +"  " \
                            +sprintf("%13.2f", PSFC(it)) +"  " \
					+sprintf("%13.2f", Q(it)) +"  " \
                                                 +sprintf("%13.2f", RH(it)) +"  " \
                                                               +sprintf("%13.2f", wind_10m(it)) +"  " \ 
                                                                            +sprintf("%13.2f", UST(it)) +"  " )
												 
end do

end
