pro debug

xSize=73
ySize=104


B_r=fltarr(xSize,ySize) 
B_phi=fltarr(xSize,ySize) 
B_z=fltarr(xSize,ySize) 


openr,1,'debug.bin'

readu,1, B_r,B_phi,B_z

close,1


B_mag=B_r^2+B_phi^2+B_z^2


B_prof_x=reform(B_mag[*,0])
B_prof_y=reform(B_mag[50,*])



plot,B_prof_y



;window,0,xsize=xSize,ysize=ySize
;tvscl,B_z

end

