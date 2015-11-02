## This code will make the maps of the distribution of green crab SNP data. 

## Load Libraries --------------
library(maps) # tool for maps
library(mapdata) # all your basemaps are here
library(marmap) # for bathymetry if needed
library(mapplots) # for add.pie
library(gplots) # for colour range

## Load in data for station coordinates -----------
coords=read.csv("Crab Coordinates.csv")

## Load in data for Haplotypes 
Haplo=read.csv("Crab_haplotypes.csv")
Haplo[is.na(Haplo)] <- 0 # convert NA's to 0

#Merge Haplo with coords
hcoords=merge(coords,Haplo,by="Code")
#add alternative points for the pie charts which are clustered together near Cape Breton
hcoords$lat=hcoords$Latitude
hcoords$long=hcoords$Longitude
##BE CAREFUL TO ONLY RUN THIS ONCE this will offset the coordinates for overlapping pies.
hcoords[which(hcoords$Code=="SYH"),"long"]=hcoords[which(hcoords$Code=="SYH"),"long"]+1.5
hcoords[which(hcoords$Code=="SYH"),"lat"]=hcoords[which(hcoords$Code=="SYH"),"lat"]+0.75
hcoords[which(hcoords$Code=="MBO"),"lat"]=hcoords[which(hcoords$Code=="MBO"),"lat"]-2
hcoords[which(hcoords$Code=="MBO"),"long"]=hcoords[which(hcoords$Code=="MBO"),"long"]+0.75


## Sampling range ####
        Sample.Lat.lim=c(38.5,52)
        Sample.Long.lim=c(-75.2,-51.9)

            png(filename = "GreenCrabSampleRange.png", 
                width = 2400, height = 2400, res = 300, bg="transparent")
            
            map("worldHires", xlim=Sample.Long.lim, ylim=Sample.Lat.lim, col="white", fill=TRUE, resolution=0);map.axes();map.scale(ratio=FALSE)
            points(coords$Longitude,  coords$Latitude,pch=19,cex=1.5) #Add points
            
            dev.off()

## Zoomed out east coast ####
        EC.Lat.lim=c(25,55)
        EC.Long.lim=c(-83,-50)

            png(filename = "EastCoast.png", 
                width = 2400, height = 2400, res = 300, bg="transparent")
            
            map("worldHires", xlim=EC.Long.lim, ylim=EC.Lat.lim, col="white", fill=TRUE, resolution=0);map.axes()
            map.scale(-65,27,ratio=FALSE,relwidth=0.2)
           
            dev.off()

## Zoomed out east coast with sample box ####

        BoxCoords=data.frame(x=c(Sample.Long.lim[1],Sample.Long.lim[1],Sample.Long.lim[2],Sample.Long.lim[2]),
                             y=c(Sample.Lat.lim[1],Sample.Lat.lim[2],Sample.Lat.lim[2],Sample.Lat.lim[1]))

        png(filename = "EastCoast_RangeBox.png", 
            width = 2400, height = 2400, res = 300, bg="transparent")
        
        map("worldHires", xlim=EC.Long.lim, ylim=EC.Lat.lim, col="white", fill=TRUE, resolution=0);map.axes()
        map.scale(-65,27,ratio=FALSE,relwidth=0.2)
        polygon(BoxCoords$x,BoxCoords$y,lty=2)
        
        dev.off()


## Sampling range ####
        
        Pie.Lat.lim=c(38,52)
        Pie.Long.lim=c(-75.9,-51.9)
        
        #specify colours smooth panel from red - blue
        PieCols=colorpanel(length(hcoords[i,c("H01","H02","H04","H05","H06","H08","ABL01","ABL02")]),
                           low="yellow",mid="red",high="darkblue")
        
        png(filename = "GreenCrabSampleRange_PieSamples.png", 
            width = 2400, height = 2400, res = 300, bg="transparent")
        
            map("worldHires", xlim=Pie.Long.lim, ylim=Pie.Lat.lim, col="white", fill=TRUE, resolution=0);map.axes();map.scale(ratio=FALSE)
           
            segments(hcoords$long[which(hcoords$Code=="SYH")],hcoords$lat[which(hcoords$Code=="SYH")],
                     hcoords$Longitude[which(hcoords$Code=="SYH")],
                     hcoords$Latitude[which(hcoords$Code=="SYH")],lty=1,lwd=1.25)
            
            segments(hcoords$long[which(hcoords$Code=="MBO")],hcoords$lat[which(hcoords$Code=="MBO")],
                     hcoords$Longitude[which(hcoords$Code=="MBO")],
                     hcoords$Latitude[which(hcoords$Code=="MBO")],lty=1,lwd=1.25)
            
            for(i in 1:nrow(hcoords))
            {
              add.pie(as.integer(hcoords[i,c("H01","H02","H04","H05","H06","H08","ABL01","ABL02")]),
                      x=hcoords$long[i],y=hcoords$lat[i],labels="",radius = 0.5,
                      col=PieCols)
            }
            
            #Add legend
            legend("bottomright",
                   c("H01","H02","H04","H05","H06","H08","ABL01","ABL02"), 
                   pch=rep(15,length(c("H01","H02","H04","H05","H06","H08","ABL01","ABL02"))),
                   col=PieCols)
               
        
        dev.off()#complete plot
        
    
