############################################
####Plotting Treemix Trees and Residuals####
############################################
########Nick Jeffery December 2015##########
############################################

#Set wd and get the 'plot_tree' function from Treemix
setwd("~/treemix-1.12")
source("src/plotting_funcs.R")

#Change wd again to where your output is located
setwd("~/treemix-1.12/src/")

#plot with whatever you called your output - the function will use all the zipped folders so you don't need to include the extension (e.g. .gz not needed)
#Save the image
png("CrabMigrationAllSNPSmig3.png",width=2400,height=2000,res = 300)
plot_tree("crabmigout3")
dev.off()

#Now plot and save the residuals
png("CrabAllSNPSmig3Resid.png",width=2400,height=2000,res=300)
plot_resid("crabmigout3","poporder") #poporder is a text file with populations in the order you want them
dev.off()

#Now move your way through all the migration events you tested (for me it was 1 to 10) and see which are best based on the residuals.