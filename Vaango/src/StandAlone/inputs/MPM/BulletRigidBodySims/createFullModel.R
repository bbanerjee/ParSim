setwd(".")

#----------------------------------------------------------------------------------
data = read.csv("Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1", 
                header = FALSE, sep = "")

reflect_x = data.frame()
reflect_x = rbind(reflect_x, data)
reflect_x$V2 = -reflect_x$V2
reflect_x$V5 = -reflect_x$V5

data_upd = rbind(data, reflect_x)
reflect_y = data.frame()
reflect_y = rbind(reflect_y, data_upd)
reflect_y$V3 = -reflect_y$V3
reflect_y$V6 = -reflect_y$V6

data_full = rbind(data_upd, reflect_y)

write.table(data_full, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.full.pos_1", 
            sep = " ", row.names = FALSE, col.names = FALSE)

#----------------------------------------------------------------------------------
data = read.csv("Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra1", 
                header = FALSE, sep = "")

orig = data.frame()
orig = rbind(orig, data)
orig$V8 = 0
write.table(orig, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra1.0", 
            sep = " ", row.names = FALSE, col.names = FALSE)

reflect_x = data.frame()
reflect_x = rbind(reflect_x, orig)
reflect_x$V2 = -reflect_x$V2
reflect_x$V5 = -reflect_x$V5
write.table(reflect_x, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra1.1", 
            sep = " ", row.names = FALSE, col.names = FALSE)

reflect_y = data.frame()
reflect_y = rbind(reflect_y, orig)
reflect_y$V3 = -reflect_y$V3
reflect_y$V6 = -reflect_y$V6
write.table(reflect_y, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra1.2", 
            sep = " ", row.names = FALSE, col.names = FALSE)

reflect_xy = data.frame()
reflect_xy = rbind(reflect_y, orig)
reflect_xy$V2 = -reflect_xy$V2
reflect_xy$V5 = -reflect_xy$V5
reflect_xy$V3 = -reflect_xy$V3
reflect_xy$V6 = -reflect_xy$V6
write.table(reflect_xy, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra1.3", 
            sep = " ", row.names = FALSE, col.names = FALSE)

#----------------------------------------------------------------------------------
data = read.csv("Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra2", 
                header = FALSE, sep = "")

orig = data.frame()
orig = rbind(orig, data)
orig$V8 = 0
write.table(orig, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra2.0", 
            sep = " ", row.names = FALSE, col.names = FALSE)

reflect_x = data.frame()
reflect_x = rbind(reflect_x, orig)
reflect_x$V2 = -reflect_x$V2
reflect_x$V5 = -reflect_x$V5
write.table(reflect_x, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra2.1", 
            sep = " ", row.names = FALSE, col.names = FALSE)

reflect_y = data.frame()
reflect_y = rbind(reflect_y, orig)
reflect_y$V3 = -reflect_y$V3
reflect_y$V6 = -reflect_y$V6
write.table(reflect_y, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra2.2", 
            sep = " ", row.names = FALSE, col.names = FALSE)

reflect_xy = data.frame()
reflect_xy = rbind(reflect_y, orig)
reflect_xy$V2 = -reflect_xy$V2
reflect_xy$V5 = -reflect_xy$V5
reflect_xy$V3 = -reflect_xy$V3
reflect_xy$V6 = -reflect_xy$V6
write.table(reflect_xy, 
            file = "Centrifuge_FlatHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1_extra2.3", 
            sep = " ", row.names = FALSE, col.names = FALSE)

#----------------------------------------------------------------------------------
data = read.csv("Centrifuge_RoundHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1", 
                header = FALSE, sep = "")

reflect_x = data.frame()
reflect_x = rbind(reflect_x, data)
reflect_x$V2 = -reflect_x$V2
reflect_x$V5 = -reflect_x$V5

data_upd = rbind(data, reflect_x)
reflect_y = data.frame()
reflect_y = rbind(reflect_y, data_upd)
reflect_y$V3 = -reflect_y$V3
reflect_y$V6 = -reflect_y$V6

data_full = rbind(data_upd, reflect_y)

write.table(data_full, 
            file = "Centrifuge_RoundHull_BoulderClay_20g_13ww_midPBC.uda.000.full.pos_1", 
            sep = " ", row.names = FALSE, col.names = FALSE)

#----------------------------------------------------------------------------------
data = read.csv("Centrifuge_VHull_BoulderClay_20g_13ww_midPBC.uda.000.pos_1", 
                header = FALSE, sep = "")

reflect_x = data.frame()
reflect_x = rbind(reflect_x, data)
reflect_x$V2 = -reflect_x$V2
reflect_x$V5 = -reflect_x$V5

data_upd = rbind(data, reflect_x)
reflect_y = data.frame()
reflect_y = rbind(reflect_y, data_upd)
reflect_y$V3 = -reflect_y$V3
reflect_y$V6 = -reflect_y$V6

data_full = rbind(data_upd, reflect_y)

write.table(data_full, 
            file = "Centrifuge_VHull_BoulderClay_20g_13ww_midPBC.uda.000.full.pos_1", 
            sep = " ", row.names = FALSE, col.names = FALSE)


