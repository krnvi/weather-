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


slp=wrf_user_getvar(infile,"slp",-1)
wrf_smooth_2d( slp, 3 )
slp_avg1=dim_avg_n_Wrap(slp(0:12,:,:),0)
slp_avg2=dim_avg_n_Wrap(slp(13:24,:,:),0)
slp_avg3=dim_avg_n_Wrap(slp(25:36,:,:),0)
slp_avg4=dim_avg_n_Wrap(slp(37:48,:,:),0)

U10   = wrf_user_getvar(infile,"U10",-1)
V10   = wrf_user_getvar(infile,"V10",-1)

u10_1   = dim_avg_n_Wrap(U10(0:12,:,:),0)
v10_1   = dim_avg_n_Wrap(V10(0:12,:,:),0)
u10_2   = dim_avg_n_Wrap(U10(13:24,:,:),0)
v10_2   = dim_avg_n_Wrap(V10(13:24,:,:),0)
u10_3   = dim_avg_n_Wrap(U10(25:36,:,:),0)
v10_3   = dim_avg_n_Wrap(V10(25:36,:,:),0)
u10_4   = dim_avg_n_Wrap(U10(37:48,:,:),0)
v10_4   = dim_avg_n_Wrap(V10(37:48,:,:),0)

ws      = sqrt(U10^2+V10^2)
ws_1    = dim_avg_n_Wrap(ws(0:12,:,:),0)

slp_avg1@description  = "Sea Level Pressure"
slp_avg1@units        = "hpa"
u10_1@units           = "m/s"
v10_1@units           = "m/s"

  wks_type = "png"
  wks_type@wkWidth = 1000
  wks_type@wkHeight = 1000

  wks  = gsn_open_wks(wks_type,"plot_slp")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  plot=new(4,graphic)
  res                             =True
  res@gsnDraw                     = False
  res@gsnFrame                    = False
  res@gsnLeftString               = ""
  res@gsnRightString              = ""
  res@tfDoNDCOverlay              = True
  res@gsnMaximize                = True
; Plotting options for Sea Level Pressure

  levels      = ispan(970,1006,1)
  msl_psl                          = res
  msl_psl@cnLevelSelectionMode     = "ExplicitLevels"
  msl_psl@cnLineColor              = "Blue"
  msl_psl@cnInfoLabelOn            = False
  msl_psl@cnLineLabelFontHeightF   = 0.01
  msl_psl@cnLineLabelPerimOn       = False
  msl_psl@gsnContourLineThicknessesScale = 1.00
  msl_psl@cnLineLabelBackgroundColor     = -1    ; transparent
  msl_psl@cnLineThicknessF               = 2.5
  msl_psl@gsnDraw                   = False
  msl_psl@gsnFrame                  = False
  msl_psl@cnLevels                  = levels

  contour_psl1 = gsn_csm_contour(wks,slp_avg1,msl_psl)
  contour_psl2 = gsn_csm_contour(wks,slp_avg2,msl_psl)
  contour_psl3 = gsn_csm_contour(wks,slp_avg3,msl_psl)
  contour_psl4 = gsn_csm_contour(wks,slp_avg4,msl_psl)
  delete(msl_psl)

;---Wind vector plot
  vec_res                 	  = res
  vec_res@vcMinDistanceF  	  = 0.02
  vec_res@vcRefLengthF    	  = 0.02
  vec_res@vcMinFracLengthF	  = 0.2
  vec_res@vcGlyphStyle   	  = "CurlyVector"
  vec_res@vcRefAnnoOn     	  = False
  vec_res@vcRefMagnitudeF         = 5.0
  vec_res@vcMinMagnitudeF         = 1.0
  vec_res@vcLineArrowColor        = "black"           ; change vector color
  vec_res@gsnSpreadColors         = True
  vec_res@vcLineArrowThicknessF   = 2.0               ; change vector thickness
  vec_res@vcVectorDrawOrder       = "PostDraw"
  vec_res@vcMonoLineArrowColor    = True
  vec_res@gsnLeftString           = ""
  vec_res@gsnRightString          = ""
  vec_res@gsnDraw                 = False
  vec_res@gsnFrame                = False
  vector_1                        = gsn_csm_vector(wks,u10_1,v10_1,vec_res)
  vector_2                        = gsn_csm_vector(wks,u10_2,v10_2,vec_res)
  vector_3                        = gsn_csm_vector(wks,u10_3,v10_3,vec_res)
  vector_4                        = gsn_csm_vector(wks,u10_4,v10_4,vec_res)

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


  map_res@tiMainString         = "WRF Sea Level Pressure (Contour) & Wind Vectors"
  map_res=wrf_map_resources(infile,map_res)
  map_res@mpNationalLineColor  = "Red"
  map_res@mpFillColors         = (/16,8,0,20/)
  map_res@mpOceanFillColor     = "SkyBlue"
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

  overlay(map1,contour_psl1)
  overlay(map1,vector_1)

  overlay(map2,contour_psl2)
  overlay(map2,vector_2)

  overlay(map3,contour_psl3)
  overlay(map3,vector_3)

  overlay(map4,contour_psl4)
  overlay(map4,vector_4)

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
