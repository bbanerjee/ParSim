# Requires:
#  install.packages('cda')
#  install.packages('rgl')
#  install.packages('orientlib')

require('cda')
require('rgl')
require('orientlib')

rm(list = ls())

#------------------------------------------------------------------
# 1) Plot the bounding boxe
#------------------------------------------------------------------
open3d()
cube <- cube3d(color = "green", alpha = 0.3)

xmax = 0.0015
ymax = 0.0015
zmax = 0.0015
cube_scaled <- rgl::scale3d(cube, xmax, ymax, zmax)
cube_shift <- rgl::translate3d(cube_scaled, 0, 0, zmax)
shade3d(cube_shift)
for (i in 1:6)
  lines3d(t(cube_shift$vb)[cube_shift$ib[,i],])

#------------------------------------------------------------------
# 5) Plot arrows
#------------------------------------------------------------------
arrow3d(c(xmax, 0, zmax), c(1.5*xmax, 0, zmax), s = 1/3, width = 1/3, n = 10, 
        type="rotation", col="red")
arrow3d(c(0, -ymax, zmax), c(0, -1.5*ymax, zmax), s = 1/3, width = 1/3, n = 10, 
        type="rotation", col="red")
arrow3d(c(0, ymax, zmax), c(0, 1.5*ymax, zmax), s = 1/3, width = 1/3, n = 10, 
        type="rotation", col="red")
arrow3d(c(0, 0, 2*zmax), c(0, 0, 2.5*zmax), s = 1/3, width = 1/3, n = 10, 
        type="rotation", col="red")
arrow3d(c(0, 0, 0), c(0, 0, -0.5*zmax), s = 1/3, width = 1/3, n = 10, 
        type="rotation", col="red")

