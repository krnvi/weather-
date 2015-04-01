
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

  wks  = gsn_open_wks("ps","plot_slp")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  res                             =True
  ;res@MainTitle                  = "WRF Sea Level Pressure (Contour)"
  res@tiMainFont                  = "helvetica"
  res@tiMainString                = ""
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
  msl_psl@gsnContourLineThicknessesScale = 1.50
  msl_psl@cnLineLabelBackgroundColor     = -1    ; transparent
  msl_psl@cnLineThicknessF               = 2.5
  msl_psl@gsnDraw                   = False
  msl_psl@gsnFrame                  = False
  msl_psl@cnLevels                  = levels
  
  contour_psl = gsn_csm_contour(wks,slp_avg1,msl_psl)
  delete(msl_psl)

;---Wind vector plot
  vec_res                  = res
  vec_res@vcMinDistanceF   = 0.045
  vec_res@vcRefLengthF     = 0.02  
  vec_res@vcMinFracLengthF = 0.2
  vec_res@vcGlyphStyle     = "LineArrow"
  vec_res@vcRefAnnoOn      = False
  vec_res@vcRefMagnitudeF  = 5.0
  vec_res@vcMinMagnitudeF  = 1.0
  vec_res@vcLineArrowColor        = "Red"           ; change vector color
  vec_res@gsnSpreadColors          =True
  vec_res@vcLineArrowThicknessF    = 2.0               ; change vector thickness
  vec_res@vcVectorDrawOrder        = "PostDraw"    
  vec_res@vcMonoLineArrowColor     = True 
  vec_res@gsnLeftString               = ""
  vec_res@gsnRightString              = ""
  vec_res@gsnDraw       = False
  vec_res@gsnFrame      = False

  vector_1                   = gsn_csm_vector(wks,u10_1,v10_1,vec_res)
  

  map_res                      = True
  map_res@gsnFrame             = False
  map_res@gsnDraw              = False 
  map_res@gsnPaperOrientation  = "auto"
  map_res@gsnMaximize          = True
  map_res@tiMainString         = "WRF Sea Level Pressure (Contour) & Wind Vectors" 
  map_res=wrf_map_resources(infile,map_res)
  map_res@mpNationalLineColor  = "Red"
  map_res@mpFillColors         = (/16,8,0,20/)
  map_res@mpOceanFillColor     = "SkyBlue"
  map=gsn_csm_map(wks,map_res)
  
  lnres             = True
  lnres@gsLineColor = "black"
  poly=gsn_add_shapefile_polylines(wks,map,shpfile,lnres)  

  overlay(map,contour_psl)
  overlay(map,vector_1)
  draw(map)
  frame(wks)
end
