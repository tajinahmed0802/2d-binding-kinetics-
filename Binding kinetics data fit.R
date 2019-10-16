#The purpose of this code is to estimate equilibrium binding values from experimentally derived binding probability data.
#The first part fits the binding probability data to a known function and show the fitted parameters and the last portion of code uses those parameters to fit the data to a binding probability vs. time graph as well as residuals


######################################### LOAD AND PROCESS DATA #######################################################
#Set working directory 
setwd("C:/Users/Tajin/Documents/CPSC 441/HW and keys")

#call probability data for Desmoglein 3 Fc tagged protein , saved as"binding.data"
binding.data<- read.csv("Probability Data.csv", header=TRUE)

#set binding.data as list
as.list(binding.data)

#check out what the data are listed in "binding.data"
head(binding.data)

#assign "time" to variable "x"
x<- binding.data$time

#assign "probability" to variable "y"
y<- binding.data$probability

#create "mr" and set experimentally determined receptor density to 64 molecules/um^2
mr <- 64

#create "ml" and set set experimentally determined ligand density to 64 molecules/um^2
ml <- 64 

# create "Ac" for area of contact (um^2) and set it to experimentally determined value of 6.2
Ac <- 6.2

################################### CREATE FIT FUNCTION ###############################################################

#create function to estimate 2D equilibrium constant (Ka) and off-rate constant(kr) values from binding probability data
fit.func<- function(x,Ka,kr){1-exp(-mr*ml*Ac*Ka*(1-exp(-kr*x)))}

#Use non-linear least square fit to estimate Ka and kr  
params.fit <- nls(y~fit.func(x,Ka,kr),data=binding.data,start=list(Ka=0.00001,kr=1))

#show fitted parameters (Ka, kr)
params.fit

#compare data and the model by using fitted parameters (Ka,kr)
plot(x,y,xlab="Time (s)", ylab="Binding frequency",main="Binding Probability of Dsg-3 Fc on RBC")
curve(fit.func(x,0.0000385,0.628),add=TRUE)


################################# GENERATE PLOT WITH FIT #############################################################

#create a publication quality plot of the fit using ggplot2 package 
#install.packages("ggplot2")        #turn this on if package is not installed yet
library(ggplot2)

bindingplot<- ggplot(binding.data, aes(x, y)) + 
  geom_point(shape = 16, size = 3, color = 'black', show.legend = FALSE) +
  labs(title = "Binding Frequency vs. Time", x = "Time(s)", y = "Binding Probability") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) +
  stat_smooth(colour = "black",         #non linear modeling
              linetype = 1,
              se = TRUE,                #set se to false if desired
              method = "glm",           #Log-link Gaussian fitting
              formula=y~(exp(-mr*ml*Ac*Ka*(1-exp(-kr*x)))))

print(bindingplot)


#######################FIT BINDING PROBABILITY DATA WITH GENERALIZED LINEAR MODEL AND DISPLAY HOW TO PLOT RESIDUALS OF FIT################
#fit data and obtain residual plots
fit<- glm(y~x,data=binding.data)
fit
#install.packages("ggfortify")  #TURN ON IF NOT INSTALLED ALREADY 
library(ggfortify)              #call ggfortify
fortify(fit)
myfits<-autoplot(fit)           #obtain residual vs fitted, normal Q-Q, Scale-Location and Residual vs Leverage plots
print(myfits)


################################ OUTPUT THE PLOTS AS .PNG ####################################################################

#save plot as "binding plot" png files
ggsave("bindingplot.png", plot = bindingplot)

#save generated residual plots as .png file
ggsave("residual.png", plot=myfits)

############################# END OF CODE##################






