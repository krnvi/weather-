load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

  infile     = addfile("wrfout/wrfout_d01_2014-10-10_00:00:00.nc","r")  
  times   = wrf_user_list_times(infile)
  ntimes  = dimsizes(times)
  ;print (ntimes) 
  
  type   = "X11"
  ;type   = "png"
  wks=gsn_open_wks(type,"plt_Precip")

  res = True
  res@MainTitle = "WRF PRECIPITATION"

; Following are needed for control of frame advance:
  pltres = True
  pltres@gsnFrame = False

  mpres = True
  mpres@mpGeophysicalLineColor = "Black"
  mpres@mpNationalLineColor    = "Black"
  mpres@mpUSStateLineColor     = "Black"
  mpres@mpGridLineColor        = "Black"
  mpres@mpLimbLineColor        = "Black"
  mpres@mpPerimLineColor       = "Black"

  FirstTime = True
  do it=0, ntimes-1
  if (FirstTime) then            
      times_sav = times(it)
  end if
  res@TimeLabel = times(it)
 
; First get the variables we will need        

    slp = wrf_user_getvar(infile,"slp",it)  ; slp
    wrf_smooth_2d( slp, 3 )            ; smooth slp

  ; Get non-convective, convective and total precipitation
  ; Calculate tendency values                               
    rain_exp = wrf_user_getvar(infile,"RAINNC",it)
    rain_con = wrf_user_getvar(infile,"RAINC",it)
    rain_tot = rain_exp + rain_con
    rain_tot@description = "Total Precipitation"

   if( FirstTime ) then
      if ( it .eq. 0 ) then
        rain_exp_save = rain_exp
        rain_con_save = rain_con
        rain_tot_save = rain_tot
      else
        rain_exp_save = wrf_user_getvar(a,"RAINNC",it-1)
        rain_con_save = wrf_user_getvar(a,"RAINC",it-1)
        rain_tot_save = rain_exp_save + rain_con_save
        FirstTime = False
        times_sav = times(it-1)
      end if
    end if

    rain_exp_tend = rain_exp - rain_exp_save
    rain_con_tend = rain_con - rain_con_save
    rain_tot_tend = rain_tot - rain_tot_save
    rain_exp_tend@description = "Explicit Precipitation Tendency"
    rain_con_tend@description = "Param  Precipitation Tendency"
    rain_tot_tend@description = "Precipitation Tendency"

  ; Bookkeeping, just to allow the tendency at the next time step
    rain_exp_save = rain_exp
    rain_con_save = rain_con
    rain_tot_save = rain_tot

    ; Plotting options for Sea Level Pressure
        opts_psl = res          
        opts_psl@ContourParameters = (/ 900., 1100., 2. /)
        opts_psl@cnLineColor       = "Blue"
        opts_psl@cnInfoLabelOn     = False
        opts_psl@cnLineLabelFontHeightF = 0.01
        opts_psl@cnLineLabelPerimOn = False
        opts_psl@gsnContourLineThicknessesScale = 1.5
        contour_psl = wrf_contour(infile,wks,slp,opts_psl)
        delete(opts_psl)
    

      ; Plotting options for Precipitation
;        opts_r = res                        
;        opts_r@UnitLabel            = "mm"
;        opts_r@cnLevelSelectionMode = "ExplicitLevels"
;        opts_r@cnLevels             = (/ .1, .2, .4, .8, 1.6, 3.2, 6.4, \
;                                        12.8, 25.6, 51.2, 102.4/)
;        opts_r@cnFillColors         = (/"White","White","DarkOliveGreen1", \
;                                        "DarkOliveGreen3","Chartreuse", \
;                                        "Chartreuse3","Green","ForestGreen", \
;                                        "Yellow","Orange","Red","Violet"/)
;        opts_r@cnInfoLabelOn        = False
;        opts_r@cnConstFLabelOn      = False
;        opts_r@cnFillOn             = True
    

      ; Total Precipitation (color fill)
;        contour_tot = wrf_contour(infile,wks, rain_tot, opts_r)

      ; Precipitation Tendencies 
;        opts_r@SubFieldTitle = "from " + times_sav + " to " + times(it)
    
;        contour_tend = wrf_contour(infile,wks, rain_tot_tend,opts_r) ; total (color)
;        contour_res = wrf_contour(infile,wks,rain_exp_tend,opts_r)  ; exp (color)
;        opts_r@cnFillOn = False
;        opts_r@cnLineColor = "Red4"
	  ; con (red lines)
;        contour_prm = wrf_contour(infile,wks,rain_con_tend,opts_r)   
;        delete(opts_r)

      ; MAKE PLOTS                                       

       ; Sea Level Pressure
       plot=wrf_map_overlays(infile,wks,contour_psl,pltres,mpres)

        ; Total Precipitation:  Don't plot 
;         plot = wrf_map_overlays(infile,wks,contour_tot,pltres,mpres)

        ; Total Precipitation Tendency + SLP : Do plot
;          plot = wrf_map_overlays(infile,wks,(/contour_tend,contour_psl/),pltres,mpres)

 
       ; Non-Convective and Convective Precipiation Tendencies: Don't plot
;          plot = wrf_map_overlays(infile,wks,(/contour_res,contour_prm/),pltres,mpres)

    times_sav = times(it)
    FirstTime = False
  end do        ; END OF TIME LOOP

end  
quit
