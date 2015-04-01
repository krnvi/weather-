load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
  infile     = addfile("wrfout/wrfout_d01_2014-10-10_00:00:00.nc","r")  
  times   = wrf_user_list_times(infile)
  ntimes  = dimsizes(times)
  ;print (ntimes) 
  type   = "ps"
  
  res=True
  res@MainTitle                   = "WRF Sea Level Pressure (Contour)"

;  res@mpDataBaseVersion           = "Ncarg4_1"
;  res@mpMaxLatF                   = 38           ; choose subregion to plot          
;  res@mpMinLatF                   = 6
;  res@mpMaxLonF                   = 100
;  res@mpMinLonF                   = 65

  pltres = True
  pltres@gsnFrame = False

  mpres = True
  mpres@mpGeophysicalLineColor = "Black"
  mpres@mpNationalLineColor    = "Black"
  mpres@mpUSStateLineColor     = "Black"
  mpres@mpGridLineColor        = "Black"
  mpres@mpLimbLineColor        = "Black"
  mpres@mpPerimLineColor       = "Black"


  slp=wrf_user_getvar(infile,"slp",-1)

  do it=0, 70,6
  slp_avg=dim_avg_n_Wrap(slp(it:it+6,:,:),0)
  flnm="ncl_slp"+it 
  wks=gsn_open_wks(type,flnm)
  gsn_merge_colormaps(wks,"rainbow", "lightgray")
  ;gsn_draw_colormap(wks)
   
  ; Plotting options for Sea Level Pressure
  opts_psl = res
  opts_psl@ContourParameters = (/ 980., 1000., 2. /)
  opts_psl@cnLineColor       = "Black"
  opts_psl@cnInfoLabelOn     = False
  opts_psl@cnLineLabelFontHeightF = 0.01
  opts_psl@cnLineLabelPerimOn = False
  opts_psl@gsnContourLineThicknessesScale = 1.5
  contour_psl = wrf_contour(infile,wks,slp_avg,opts_psl)
  delete(opts_psl)

  plot=wrf_map_overlays(infile,wks,contour_psl,pltres,mpres)


  end do        ; END OF TIME LOOP


end
 
quit
