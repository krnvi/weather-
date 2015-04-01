load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

shpfile   = "subdiv/India_subdiv.shp" 
nfiles    = dimsizes(shpfile)

infile  = addfile("wrfout/wrfout_d01_2014-10-10_00:00:00.nc","r")
times   = wrf_user_list_times(infile)
ntimes  = dimsizes(times)

tmp = wrf_user_getvar(infile,"SST",-1)
tmp=tmp-273.15

lndmsk= wrf_user_getvar(infile,"XLAND",-1)
tmp_l  = mask(tmp,lndmsk,2)

tmp_1 = dim_avg_n_Wrap(tmp_l(0:11,:,:),0)
tmp_2 = dim_avg_n_Wrap(tmp_l(12:23,:,:),0)
tmp_3 = dim_avg_n_Wrap(tmp_l(24:35,:,:),0)
tmp_4 = dim_avg_n_Wrap(tmp_l(36:47,:,:),0)

tmp@description = "Temperature degC"
  wks_type = "png"
  wks_type@wkWidth = 1000
  wks_type@wkHeight = 1000
  wks  = gsn_open_wks(wks_type,"plot_sst")
  ;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  ;gsn_define_colormap(wks,"sunshine_9lev")
  
  plot=new(4,graphic)
  res                             =True
  res@gsnDraw                     = False
  res@gsnFrame                    = False
  res@gsnLeftString               = ""
  res@gsnRightString              = ""
  res@tfDoNDCOverlay              = True
  res@gsnMaximize                = True
; Plotting options for Sea Level Pressure

  levels      = ispan(1,100,1)
  res_rf                                = res
  res_rf@cnLevelSelectionMode           = "ExplicitLevels"
  res_rf@cnLevels                       = (/25,27,29,31,33/)
  ;res_rf@cnFillColors                   = (/0,88,80,56,32,246,241,237,231,228,222,220/)  
  res_rf@cnFillColors                   = (/0,252,249,246,242,239,236,232,229,225,222,220/)  
  ;res_rf@cnFillColors                   = (/0,108,103,98,93,88,80,72,64,66,58,50/)  
  ;res_rf@cnFillColors                   = (/"white","gold1","gold2","goldenrod1","goldenrod2","goldenrod3"/)  
  res_rf@cnLineColor                    = "Black"
  res_rf@cnInfoLabelOn                  = False
  res_rf@cnLineLabelFontHeightF         = 0.01
  res_rf@cnLineLabelPerimOn             = False
  res_rf@gsnContourLineThicknessesScale = 1.00
  res_rf@cnLineLabelBackgroundColor     = -1    ; transparent
  res_rf@cnLineThicknessF               = 2.5
  res_rf@gsnDraw                        = False
  res_rf@gsnFrame                       = False
  ;res_rf@cnLevels                       = levels
  res_rf@cnFillOn                       = True

  contour_tmp1 = gsn_csm_contour(wks,tmp_1,res_rf)
  contour_tmp2 = gsn_csm_contour(wks,tmp_2,res_rf)
  contour_tmp3 = gsn_csm_contour(wks,tmp_3,res_rf)
  contour_tmp4 = gsn_csm_contour(wks,tmp_4,res_rf)
  delete(res_rf)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  res_rf                                = res
  res_rf@cnLevelSelectionMode           = "ExplicitLevels"
  res_rf@cnLevels                       = (/25,27,29,31,33,35/)
  res_rf@cnLineColor                    = "Black"
  res_rf@cnInfoLabelOn                  = False
  res_rf@cnLineLabelFontHeightF         = 0.01
  res_rf@cnLineLabelPerimOn             = False
  res_rf@gsnContourLineThicknessesScale = 1.00
  res_rf@cnLineLabelBackgroundColor     = -1    ; transparent
  res_rf@cnLineThicknessF               = 2.5
  res_rf@gsnDraw                        = False
  res_rf@gsnFrame                       = False
  res_rf@cnFillOn                       = False

  contour_utmp1 = gsn_csm_contour(wks,tmp_1,res_rf)
  contour_utmp2 = gsn_csm_contour(wks,tmp_2,res_rf)
  contour_utmp3 = gsn_csm_contour(wks,tmp_3,res_rf)
  contour_utmp4 = gsn_csm_contour(wks,tmp_4,res_rf)
  delete(res_rf)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  map_res                      = True
  map_res@gsnFrame             = False
  map_res@gsnDraw              = False
  map_res@gsnPaperOrientation  = "auto"
  map_res@gsnMaximize          = True
  map_res@tiMainFont           = 25
  map_res@tiMainFontHeightF    = 0.020
  ;map_res@tiMainJust           = "CenterLeft"
  map_res@tiMainPosition       = "Center"
  map_res@tiYAxisAngleF        = 0.0
  map_res@tiYAxisFont          = 25
  map_res@tiYAxisFontAspectF   = 1.5
  map_res@tiYAxisFontHeightF   = 0.015
  map_res@tiYAxisFuncCode      = "\"
  map_res@tiYAxisJust          = "BottomRight"

; Resources for X and Y axes. The top, left, and right axes
; will be turned off.

  map_res@tmXTOn               = False      ; Turn off top tickmarks
  map_res@tmYLOn               = True      ; Turn off left tickmarks
  map_res@tmYROn               = False      ; Turn off bottom tickmarks
  map_res@tmXTBorderOn         = True      ; Turn off top border
  map_res@tmYLBorderOn         = True      ; Turn off left border
  map_res@tmYRBorderOn         = True      ; Turn off right border;

; Resources for the tickmarks on the bottom axis.

  map_res@tmXBMode             = "Manual"
  map_res@tmXBTickSpacingF     = 2.5
  map_res@tmXBLabelFont        = 25
  map_res@tmXBLabelFontHeightF = 0.015
  map_res@tmXBMajorLengthF     = 0.02
  map_res@tmXBMinorLengthF     = 0.01
  map_res@tmXBMinorPerMajor    = 4


; Control range of X and Y axis.
;  map_res@trXMaxF              = 20.0
;  map_res@trXMinF              = 0.0
;  map_res@trYMinF              = 0.0


; Width, height, and position of X axis. Every time the
; axis is redrawn, the vpYF resource will be changed to change
; the position.

;  map_res@vpXF                 = 0.2
;  map_res@vpYF                 = 0.9
;  map_res@vpHeightF            = 0.02
;  map_res@vpWidthF             = 0.7


  map_res@tiMainString         = "Sea Surface Temperature (Deg C)"
  map_res=wrf_map_resources(infile,map_res)
  map_res@mpNationalLineColor  = "Red"
  map_res@mpFillColors         = (/0,8,0,20/)
  map_res@mpOceanFillColor     = "White"
  map1=gsn_csm_map(wks,map_res)
  map2=gsn_csm_map(wks,map_res)
  map3=gsn_csm_map(wks,map_res)
  map4=gsn_csm_map(wks,map_res)

  lnres             = True
  lnres@gsLineColor = "black"
  poly1=gsn_add_shapefile_polylines(wks,map1,shpfile,lnres)
  poly2=gsn_add_shapefile_polylines(wks,map2,shpfile,lnres)
  poly3=gsn_add_shapefile_polylines(wks,map3,shpfile,lnres)
  poly4=gsn_add_shapefile_polylines(wks,map4,shpfile,lnres)

  overlay(map1,contour_tmp1)
  overlay(map2,contour_tmp2)
  overlay(map3,contour_tmp3)
  overlay(map4,contour_tmp4)
 
  overlay(map1,contour_utmp1)
  overlay(map2,contour_utmp2)
  overlay(map3,contour_utmp3)
  overlay(map4,contour_utmp4)
 

  
  plot(0)=map1
  plot(1)=map2
  plot(2)=map3
  plot(3)=map4

  pres             = True
  pres@gsnMaximize = True
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5

  gsn_panel(wks,(/plot(0),plot(1),plot(2),plot(3)/),(/2,2/),pres)


end
