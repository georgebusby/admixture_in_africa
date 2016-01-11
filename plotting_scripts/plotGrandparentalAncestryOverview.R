

### PLOT CHROMOSOMES 

library(plotrix) ## for ellipses

# EMPTY PLOT


chromlength <- 1
chromheight <- 0.5
dist2box <- 0.25
dist2boxhap <- 0.1
indlwd <- 2
left1 <- 0.5
left2 <- 2.5
bottom1 <- 9
bottom2 <- 8
bottom3 <- 7
bottom4 <- 6
bottom5 <- 5    
bottom6 <- 4    

par(mar=c(1,1,1,1))
plot(0,0,xlim=c(0,4),ylim=c(0,10),type='n')

#### FATHER ####
rect(left1-dist2box,bottom2-dist2box,
     left1+chromlength+dist2box,bottom1+chromheight+dist2box,
     lwd=indlwd)
# father father
rect(left1,bottom1,
     left1+chromlength,bottom1+chromheight,
     col="lightgreen")
# father mother
rect(left1,bottom2,
     left1+chromlength,bottom2+chromheight,
     col="hotpink")

#### MOTHER ####
rect(left2-dist2box,bottom2-dist2box,
     left2+chromlength+dist2box,bottom1+chromheight+dist2box,
     lwd=indlwd)
# mother father
rect(left2,bottom1,
     left2+chromlength,bottom1+chromheight,
     col="darkgreen")
# mother mother
rect(left2,bottom2,
     left2+chromlength,bottom2+chromheight,
     col="brown1")


#### GAMETES #####
adm1 <- 0.75
adm2 <- 0.6
    
#### FATHER GAMETES
# gamete 1
rect(left1-dist2box,bottom3-dist2boxhap,
     left1+chromlength+dist2box,bottom3+chromheight+dist2boxhap,
     lwd=indlwd)

rect(left1,bottom3,
     left1+chromlength,bottom3+chromheight,
     col="lightgreen")
rect(left1+(chromlength*adm1),bottom3,
     left1+chromlength,bottom3+chromheight,
     col="hotpink")

########################
# gamete 2
rect(left1-dist2box,bottom4-dist2boxhap,
     left1+chromlength+dist2box,bottom4+chromheight+dist2boxhap,
     lwd=indlwd)

rect(left1,bottom4,
     left1+chromlength,bottom4+chromheight,
     col="hotpink")
rect(left1+(chromlength*adm1),bottom4,
     left1+chromlength,bottom4+chromheight,
     col="lightgreen")

#### MOTHER GAMETES
rect(left2-dist2box,bottom3-dist2boxhap,
     left2+chromlength+dist2box,bottom3+chromheight+dist2boxhap,
     lwd=indlwd)
rect(left2,bottom3,
     left2+chromlength,bottom3+chromheight,
     col="darkgreen")
rect(left2+(chromlength*adm2),bottom3,
     left2+chromlength,bottom3+chromheight,
     col="brown1")
# gamete 2
rect(left2-dist2box,bottom4-dist2boxhap,
     left2+chromlength+dist2box,bottom4+chromheight+dist2boxhap,
     lwd=indlwd)
rect(left2,bottom4,
     left2+chromlength,bottom4+chromheight,
     col="brown1")
rect(left2+(chromlength*adm2),bottom4,
     left2+chromlength,bottom4+chromheight,
     col="darkgreen")


##### ME#
mex <- mean(c(left1,left2))
rect(mex-dist2box,bottom6-dist2boxhap,
     mex+chromlength+dist2box,bottom5+chromheight+dist2boxhap,
     lwd=indlwd)

## first paternal gamete
rect(mex,bottom5,
     mex+chromlength,bottom5+chromheight,
     col="lightgreen")
rect(mex+(chromlength*adm1),bottom5,
     mex+chromlength,bottom5+chromheight,
     col="hotpink")
## second maternal gamete
rect(mex,bottom6,
     mex+chromlength,bottom6+chromheight,
     col="brown1")
rect(mex+(chromlength*adm2),bottom6,
     mex+chromlength,bottom6+chromheight,
     col="darkgreen")

##########################################################################################################
## SPACING
## ARROWS
## RECOMBINATION



