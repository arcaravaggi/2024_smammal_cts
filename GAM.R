setwd("C:/Users/phoeb/OneDrive/Desktop/UNI/RESEARCH PROJECT/Data")

library(lubridate)
library(mgcv)
library(mgcViz)
library(reshape2)
library(multcomp)
library(tidyverse)
library(GGally)
library(overlap)
library(scales)
library(ggimage)
library(activity)
library(nlme)

# Initial steps ####
#import data 
camdata <- read.csv("Camera Trap Data.csv")
names(camdata)

# Exclude records that aren't small mammals
camdata <- camdata[!(camdata$animal %in% 
                        c("bird", "invertebrate", "mammal", 
                          "none", "reptile")),]

class(camdata$date) # Check class of Date column; should be 'Date'

camdata$date.time <- as.POSIXct(paste0(camdata$time," ", camdata$date), format = "%H:%M:%S %d/%m/%Y")
class(camdata$date.time)

# Add columns for day and month
camdata$day <- day(camdata$date.time)
camdata$month <- month(camdata$date.time)

# Create a data frame for relevant weeks
camdata$week <- NA

# Add week data based on values in day and month
camdata$week[camdata$day %in% c(11:17) & camdata$month == 7] <- 1

camdata$week[camdata$day %in% c(18:24) & camdata$month == 7] <- 2

camdata$week[camdata$day %in% c(25:31) & camdata$month == 7] <- 3

camdata$week[camdata$day %in% c(1:7) & camdata$month == 8] <- 4

camdata$week[camdata$day %in% c(8:14) & camdata$month == 8] <- 5

camdata$week[camdata$day %in% c(15:21) & camdata$month == 8] <- 6

camdata$week[camdata$day %in% c(22:28) & camdata$month == 8] <- 7

camdata$week[camdata$day %in% c(29:31) & camdata$month == 8] <- 8

camdata$week[camdata$day %in% c(1:4) & camdata$month == 9] <- 8

camdata$week[camdata$day %in% c(5:6) & camdata$month == 9] <- 9

# extract rows with NA
camdata[is.na(camdata$week),]

# Add count column
camdata$count <- 1
names(camdata)
camdata.2 <- camdata[,c(3, 5, 13, 15)]
names(camdata.2)

# Create cumulative counts per week and species
camdata.3 <- camdata.2 %>%
  group_by(species, camera, week) %>%
  mutate(f = 1) %>%
  summarize(count = sum(f)) %>%
  mutate(sp_group = species) %>%
  mutate(sp_group = replace(sp_group, sp_group %in% c("bank vole", "field vole", "vole sp."), "vole")) %>%
  mutate(sp_group = replace(sp_group, sp_group %in% c("common shrew", "pygmy shrew", "water shrew", "shrew sp."), "shrew"))

# Remove unidentified records
camdata.4 <- camdata.3[!(camdata.3$sp_group == "unidentified"),]

m1 <- gamm(count ~ s(week, bs = "tp", k = 9) + 
             as.factor(camera) * as.factor(sp_group), 
           family=poisson, 
           data = camdata.4)
summary(m1$gam)
plot(m1$gam)

# Plot models ####
b.m1 <- getViz(m1$gam)
plot(sm(b.m1, 1)) + l_fitLine(colour = "red") + 
  l_rug(mapping = aes(y=y), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = 5, colour = "blue", linetype = 2) +
  l_points(shape = 19, size = 1, alpha = 0.1) +
  scale_x_continuous(name="Week", breaks=c(1,2,3,4,5,6,7,8,9), 
                     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                     limits=c(1, 9))+ 
  geom_hline(yintercept = 0, lty=2) +
  ylab("Smoothed number of detections")  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))


# Plot of counts by species group
camdata.4 %>% 
  ggplot(aes(x = sp_group, y = after_stat(count), fill = camera)) +
  geom_bar(position = position_dodge(), stat = "prop") 

# Activity analyses
# Recode month
camdata$month <- recode(camdata$month, "1" = "01_Jan", "2" = "02_Feb", "3" = "03_Mar", "4" = "04_Apr", 
                       "5" = "05_May", "6" = "06_Jun", "7" = "07_Jul", "8" = "08_Aug", "9" = "09_Sep",
                       "10" = "10_Oct", "11" = "11_Nov", "12" = "12_Dec")

# Extract hour from time
camdata$hour <- hour(hms(as.character(camdata$time)))

camdata <- camdata %>%
  mutate(sp_group = species) %>%
  mutate(sp_group = replace(sp_group, sp_group %in% c("bank vole", "field vole", "vole sp."), "vole")) %>%
  mutate(sp_group = replace(sp_group, sp_group %in% c("common shrew", "pygmy shrew", "water shrew", "shrew sp."), "shrew"))

camdata <- camdata[!(camdata$sp_group == "unidentified"),]

# Overlap plots ####
# Convert time to decimal and rescale to between 0 and 1
camdata$time_adj <- sapply(strsplit(camdata$time,":"), 
                      function(x) {
                        x <- as.numeric(x)
                        x[1]+x[2]/60
                      }
)
camdata$time_adj <- scales:::rescale(camdata$time_adj, to = c(0, 1))

# Convert to radians
timeRad.fr <- camdata$time_adj  * 2 * pi
mouse.act <- timeRad.fr[camdata$sp_group == 'mouse']
vole.act <- timeRad.fr[camdata$sp_group == 'vole']
shrew.act <- timeRad.fr[camdata$sp_group == 'shrew']

# Calculate overlap
mousevole.over <- overlapEst(mouse.act, vole.act, type="Dhat4")
mouseshrew.over <- overlapEst(mouse.act, shrew.act, type="Dhat4")
voleshrew.over <- overlapEst(vole.act, shrew.act, type="Dhat4")

# Bootstrap data
mouseboot <- resample(mouse.act, 1000) # 1000 resamples
voleboot <- resample(vole.act, 1000) # 1000 resamples
shrewboot <- resample(shrew.act, 1000) # 1000 resamples

mousevole<- bootEst(mouseboot, voleboot, type="Dhat4") # takes a few seconds
(BSmean.mousevole <- mean(mousevole))
mv.cis <- bootCI(mousevole.over, mousevole)
mv.cis [3,]
       
# Plot overlap
dev.control(displaylist="enable")
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)

par(family = 'sans')
overlapPlot(mouse.act, vole.act, main=NULL, font.lab=2 )
legend("topright", inset=c(-0.2,0), c("Mouse", "Vole"), lty=c(1,2), col=c(1,4), bty='n', title = "Species")
text(12, 0.10, bquote(Delta == .(round(mousevole.over, 2)) ~ "(" * .(round(mv.cis[3,1], 2)) * "-" * 
                      .(round(mv.cis[3,2],2)) * ")"))
p2 <- recordPlot()


# Group activity by species group ####
actsp <- camdata %>%  
  group_by(sp_group, camera)  %>% 
  count(hour) %>% #
  mutate(freq = n / sum(n)) %>%
  mutate(freq_scale = rescale(freq, to = c(0.1, .9))) %>%
  na.omit(camera)

# Create hourly data and merge to fill in missing hours
sp <- unique(camdata$sp_group)
y <- expand.grid(hour = 0:23, sp_group = sp, camera = unique(camdata$camera))
actsp <- merge(actsp, y, all = TRUE) %>% replace_na(list(n = 0, freq = 0))


# Mouse data - activity and overlap calculations ####
mact.u <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'mouse' & camdata$camera == "no box"], #convert to radians
                           sample = "model", reps = 1000)
mact.o <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'mouse' & camdata$camera == "open box"], 
                           sample = "model", reps = 1000)
mact.t <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'mouse' & camdata$camera == "tube box"], 
                           sample = "model", reps = 1000)

# Compare activity patterns
mouse.comp <- compareAct(list(mact.u, mact.o, mact.t))

# Calculate overlap
mouse.uo <- compareCkern(mact.u, mact.o, reps = 1000)
mouse.ut <- compareCkern(mact.u, mact.t, reps = 1000)
mouse.ot <- compareCkern(mact.o, mact.t, reps = 1000)

# Rescale for plotting
mact.u@pdf[,2] <- rescale(mact.u@pdf[,2], to = c(0, 1))
mact.o@pdf[,2] <- rescale(mact.o@pdf[,2], to = c(0, 1))
mact.t@pdf[,2] <- rescale(mact.t@pdf[,2], to = c(0, 1))

# Mouse data -  Plotting overlap  ####
dev.control(displaylist="enable")
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
par(family = 'sans')

plot(mact.u, yunit="density", data="none", las=1, lwd=2, lty = 1,
     tline=list(lwd=2, lty = 1), # Thick line 
     cline=list(lty=0)) # Supress confidence intervals

plot(mact.o, yunit="density", data="none", add=TRUE, 
     tline=list(col="red", lwd=2, lty = 2),
     cline=list(lty=0))

plot(mact.t, yunit="density", data="none", add=TRUE, 
     tline=list(col="blue", lwd=2, lty = 3),
     cline=list(lty=0))

legend("top", 
       inset=c(-0.26,-0.16), 
       c("Unboxed (U)", "Open box (O)", "Tube box (T)"), 
       col=c(1,2,4), 
       lty=c(1,2,3), 
       lwd=2,
       bty="n", 
       horiz = T)

mocorners <- par("usr")

text(x = mocorners[2], y = mocorners[4]-0.0123, pos = 4,
     bquote(bold("Overlap (" ~ Delta ~ ")")), xpd = T)
text(x = mocorners[2], y = mocorners[4]-0.0243, pos = 4,
     bquote("U-O" == .(round(mouse.uo[1], 2)) ~ "(P = " * .(round(1-mouse.uo[4], 3)) * ")"), xpd=T)
text(x = mocorners[2], y = mocorners[4]-0.0363, pos = 4,
     bquote("U-T" == .(round(mouse.ut[1], 2)) ~ "(P = " * .(round(1-mouse.ut[4], 3)) * ")"), xpd=T)
text(x = mocorners[2], y = mocorners[4]-0.0483, pos = 4,
     bquote("O-T" == .(round(mouse.ot[1], 2)) ~ "(P = 0.952)"), xpd=T)

text(x = mocorners[2], y = mocorners[4]-0.0723, pos = 4, 
     bquote(bold("Sample size (n)")), xpd = T)
text(x = mocorners[2], y = mocorners[4]-0.0843, pos = 4,
     bquote("U" == .(length(mact.u@data))))
text(x = mocorners[2], y = mocorners[4]-0.0963, pos = 4, 
     bquote("O" == .(length(mact.o@data))))
text(x = mocorners[2], y = mocorners[4]-0.1083, pos = 4,
     bquote("T" == .(length(mact.t@data))))

grid::grid.raster(png::readPNG('mouse.png'), x = .87, y=0.25, width = .15) 

text(x = mocorners[1]-2.75, y = mocorners[4]+.015, pos = 4, 
     bquote(bold("a)")), xpd = T, cex = 1.2)

mouse.activity.plot <- recordPlot()
invisible(dev.off())

png("mouse_activity.png", width = 8.5, height = 4.5, units = 'in', res = 300)
mouse.activity.plot
dev.off()

# Vole data - activity and overlap calculations ####
voct.u <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'vole' & camdata$camera == "no box"], #convert to radians
                 sample = "model", reps = 1000)
voct.o <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'vole' & camdata$camera == "open box"], 
                 sample = "model", reps = 1000)
voct.t <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'vole' & camdata$camera == "tube box"], 
                 sample = "model", reps = 1000)

# Compare activity patterns
vole.comp <- compareAct(list(voct.u, voct.o, voct.t))

# Calculate overlap
vole.uo <- compareCkern(voct.u, voct.o, reps = 1000)
vole.ut <- compareCkern(voct.u, voct.t, reps = 1000)
vole.ot <- compareCkern(voct.o, voct.t, reps = 1000)

# Rescale for plotting
voct.u@pdf[,2] <- rescale(voct.u@pdf[,2], to = c(0, 1))
voct.o@pdf[,2] <- rescale(voct.o@pdf[,2], to = c(0, 1))
voct.t@pdf[,2] <- rescale(voct.t@pdf[,2], to = c(0, 1))

# Vole data - Plotting overlap ####
dev.control(displaylist="enable")
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
par(family = 'sans')

plot(voct.u, yunit="density", data="none", las=1, lwd=2, lty = 1,
     tline=list(lwd=2, lty = 1), # Thick line 
     cline=list(lty=0))

plot(voct.o, yunit="density", data="none", add=TRUE, 
     tline=list(col="red", lwd=2, lty = 2),
     cline=list(lty=0))

plot(voct.t, yunit="density", data="none", add=TRUE, 
     tline=list(col="blue", lwd=2, lty = 3),
     cline=list(lty=0))

legend("top", 
       inset=c(-0.26,-0.16), 
       c("Unboxed (U)", "Open box (O)", "Tube box (T)"), 
       col=c(1,2,4), 
       lty=c(1,2,3), 
       lwd=2,
       bty="n", 
       horiz = T)

vocorners <- par("usr")

text(x = vocorners[2], y = vocorners[4]-0.0223, pos = 4,
     bquote(bold("Overlap (" ~ Delta ~ ")")), xpd = T)
text(x = vocorners[2], y = vocorners[4]-0.0423, pos = 4, 
     bquote("U-O" == .(round(vole.uo[1], 2)) ~ "(P = " * .(round(vole.uo[4], 3)) * ")"), xpd=T)
text(x = vocorners[2], y = vocorners[4]-0.0623, pos = 4, 
     bquote("U-T" == .(format(round(vole.ut[1], digits = 2), nsmall = 2)) ~ "(P = " * .(round(vole.ut[4], 3)) * ")"), xpd=T)
text(x = vocorners[2], y = vocorners[4]-0.0823, pos = 4,
     bquote("O-T" == .(format(round(vole.ot[1], digits = 2), nsmall = 2)) ~ "(P = " * .(round(vole.ot[4], 3)) * ")"), xpd=T)

text(x = vocorners[2], y = vocorners[4]-0.1223, pos = 4,
     bquote(bold("Sample size (n)")), xpd = T)
text(x = vocorners[2], y = vocorners[4]-0.1423, pos = 4,
     bquote("U" == .(length(voct.u@data))))
text(x = vocorners[2], y = vocorners[4]-0.1623, pos = 4, 
     bquote("O" == .(length(voct.o@data))))
text(x = vocorners[2], y = vocorners[4]-0.18223, pos = 4,
     bquote("T" == .(length(voct.t@data))))

grid::grid.raster(png::readPNG('vole.png'), x = .87, y=0.25, width = .15) 

text(x = mocorners[1]-2.75, y = mocorners[4]+.135, pos = 4, 
     bquote(bold("c)")), xpd = T, cex = 1.2)

vole.activity.plot <- recordPlot()
invisible(dev.off())

png("vole_activity.png", width = 8.5, height = 4.5, units = 'in', res = 300)
vole.activity.plot
dev.off()

# Shrew data - activity and overlap calculations ####
shct.u <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'shrew' & camdata$camera == "no box"], #convert to radians
                 sample = "model", reps = 1000)
shct.o <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'shrew' & camdata$camera == "open box"], 
                 sample = "model", reps = 1000)
shct.t <- fitact((camdata$time_adj  * 2 * pi)[camdata$sp_group == 'shrew' & camdata$camera == "tube box"], 
                 sample = "model", reps = 1000)

# Compare activity patterns
shrew.comp <- compareAct(list(shct.u, shct.o, shct.t))

# Calculate overlap
shrew.uo <- compareCkern(shct.u, voct.o, reps = 1000)
shrew.ut <- compareCkern(shct.u, voct.t, reps = 1000)
shrew.ot <- compareCkern(shct.o, voct.t, reps = 1000)

# Rescale for plotting
shct.u@pdf[,2] <- rescale(shct.u@pdf[,2], to = c(0, 1))
shct.o@pdf[,2] <- rescale(shct.o@pdf[,2], to = c(0, 1))
shct.t@pdf[,2] <- rescale(shct.t@pdf[,2], to = c(0, 1))

# Shrew data - Plotting overlaps #### 
dev.control(displaylist="enable")
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
par(family = 'sans')

plot(shct.u, yunit="density", data="none", las=1, lwd=2, lty = 1,
     tline=list(lwd=2, lty = 1), # Thick line 
     cline=list(lty=0)) # Supress confidence intervals

plot(shct.o, yunit="density", data="none", add=TRUE, 
     tline=list(col="red", lwd=2, lty = 2),
     cline=list(lty=0))

plot(shct.t, yunit="density", data="none", add=TRUE, 
     tline=list(col="blue", lwd=2, lty = 3),
     cline=list(lty=0))

legend("top", 
       inset=c(-0.26,-0.16), 
       c("Unboxed (U)", "Open box (O)", "Tube box (T)"), 
       col=c(1,2,4), 
       lty=c(1,2,3), 
       lwd=2,
       bty="n", 
       horiz = T)

shcorners <- par("usr")

text(x = shcorners[2], y = shcorners[4]-0.0223, pos = 4,
     bquote(bold("Overlap (" ~ Delta ~ ")")), xpd = T)
text(x = shcorners[2], y = shcorners[4]-0.0423, pos = 4,
     bquote("U-O" == .(round(shrew.uo[1], 2)) ~ "(P < 0.001)"), xpd=T)
text(x = shcorners[2], y = shcorners[4]-0.0623, pos = 4,
     bquote("U-T" == .(format(round(shrew.ut[1], digits = 2), nsmall = 2)) ~ "(P < 0.001)"), xpd=T)
text(x = shcorners[2], y = shcorners[4]-0.0823, pos = 4, 
     bquote("O-T" == .(format(round(shrew.ot[1], digits = 2), nsmall = 2)) ~ "(P < 0.001)"), xpd=T)

text(x = shcorners[2], y = shcorners[4]-0.1223, pos = 4,
     bquote(bold("Sample size (n)")), xpd = T)
text(x = shcorners[2], y = shcorners[4]-0.1423, pos = 4,
     bquote("U" == .(length(shct.u@data))))
text(x = shcorners[2], y = shcorners[4]-0.1623, pos = 4,
     bquote("O" == .(length(shct.o@data))))
text(x = shcorners[2], y = shcorners[4]-0.1823, pos = 4,
     bquote("T" == .(length(shct.t@data))))

grid::grid.raster(png::readPNG('shrew.png'), x = .87, y=0.25, width = .15) 

text(x = mocorners[1]-2.75, y = mocorners[4]+.135, pos = 4, 
     bquote(bold("b)")), xpd = T, cex = 1.2)

shrew.activity.plot <- recordPlot()
invisible(dev.off())

png("shrew_activity.png", width = 8.5, height = 4.5, units = 'in', res = 300)
shrew.activity.plot
dev.off()

# Save/load image #####
#save.image(file='myEnvironment.RData')
load(file='myEnvironment.RData')

# Latency to first detection ####
# extract min and max dates for each species group and station
min_max <- camdata %>%
  group_by(station, camera) %>%
  summarize(min_date = min(date.time, na.rm=T),
            max_date = max(date.time, na.rm=T))

# create columns for set date and retrieval date
min_max$set_date <- as.POSIXct(paste(as.character(date(min_max$min_date)-1), "18:00:00", sep =  " "), format = "%Y-%m-%d %H:%M:%S")
min_max$col_date <- as.POSIXct(paste(as.character(date(min_max$max_date)+1), "18:00:00", sep =  " "), format = "%Y-%m-%d %H:%M:%S")

# Create all sequences based on min_max dates
seq_list <- list()
for(i in 1:nrow(min_max)) {
  seq_list[[i]] <- data.frame(date = seq(min_max$set_date[[i]], min_max$col_date[[i]], by = "hour"))
  # do stuff with row
}

b <- c(paste0("m", a, sep = ""), paste0("s", a, sep = ""), paste0("v", a, sep = ""))
names(seq_list) <- paste(substr(min_max$camera, 1,1), min_max$station, sep = "")
seq_cams <- bind_rows(seq_list, .id = "cams") # Collapse list to single data frame

# create camera and station columns for merging
names(seq_cams)[names(seq_cams) == "date"] <- "date.time"
seq_cams$station <- as.vector(unlist(lapply(strsplit(seq_cams$cams, ""), `[`, 2)))
seq_cams$camera <- as.vector(unlist(lapply(strsplit(seq_cams$cams, ""), `[`, 1)))
seq_cams$camera[seq_cams$camera  == "n"] <- "no box" 
seq_cams$camera[seq_cams$camera  == "o"] <- "open box"
seq_cams$camera[seq_cams$camera  == "t"] <- "tube box"
  

# combine with camdata, keeping all rows
sm.cams <- merge(camdata, seq_cams, by = c("camera", "station", "date.time"), 
                 all = T) 

# Change NA presence to 0
sm.cams$count[is.na(sm.cams$count)] <- 0

# Split by camera type
cam.split <- sm.cams %>% 
  group_by(camera, station) %>%
  arrange(date,time) %>% 
  ungroup() %>%
  group_split(camera, station) %>%
  set_names("n1", "n2", "n4", "n6", "o1", "o2", "o3", "o4", "o5", "o6", "t2", "t3", "t4", "t5", "t6")
  
cam.split <- lapply(cam.split, function(df) df[order(df$date.time),])

# Function to randomly sample data according to null followed by detection
smoot <- function(df){
  df$time2 <- strftime(df$date.time, '%H:%M:%S')
  d <- df[sample(which(df$count == 0 & 
                                   (as.POSIXct(df$time2, format = "%H:%M:%S") <= format(as.POSIXct("18:00:00", format = "%H:%M:%S")) |
                                    as.POSIXct(df$time2, format = "%H:%M:%S") >= format(as.POSIXct("07:00:00", format = "%H:%M:%S")))), 
                           1), ] # randomly sample row with count == 0 and between 07:00 and 18:00
  e <- df[df$date.time > d$date.time & 
                                df$count == 1, ] # find next instance of species detection
  
  if (nrow(e) == 0) {
    # No matching e found
    return(NULL)
  }
  
  e <- e[1,] # take the first row
  f <- as.numeric(difftime(e$date.time, d$date.time), "hours") # Calculate time between detections
  
  if (is.na(f) || length(f) == 0) {
    return(NULL)
  }
  
  return(data.frame(set = d$date.time, detection = e$date.time, diff = f))
  }

# Test function
p <- smoot(cam.split$n1)

# Repeat N times
n_replicates <- 1000  # Number of replicates
set.seed(1507)
final_df <- map_df(1:n_replicates, ~ smoot(cam.split$n1)) 
final_df2 <- final_df %>% distinct(set, detection) # remove duplicates

# Run across all data frames
set.seed(1308)
stat.boot <- lapply(cam.split, function(x) map_df(1:n_replicates, ~ smoot(x)))
  
#  Remove duplicates across list of data frames based on specific columns
columns_to_check <- c("set", "detection")

stat.boot2 <- stat.boot %>%
  bind_rows(.id = "source") %>%    # Combine all data frames with an identifier column for splitting
  distinct(across(all_of(columns_to_check)), .keep_all = TRUE) %>%  # Remove duplicates based on specified columns
  group_split(source)  # Split back into the list of data frames by the original source

# Convert to data frame and prepare for plotting
stat.df <- do.call(rbind.data.frame, stat.boot2)
stat.df$cam_type[substr(stat.df$source, 1, 1) == "n"] <- "Unboxed"
stat.df$cam_type[substr(stat.df$source, 1, 1) == "o"] <- "Open box"
stat.df$cam_type[substr(stat.df$source, 1, 1) == "t"] <- "Tube box"
stat.df$station <- as.factor(substr(stat.df$source, 2, 2))

latency.data <- stat.df %>% group_by(cam_type) %>%
  summarise(m = mean(diff, na.rm = T),
            sd = sd(diff, na.rm = T),
            se = sd(diff, na.rm = T)/sqrt(length(diff)))

# ggplot theme
theme_ac1 <- function(base_family = "serif", base_size_a = 12, base_size_t = 12){
  theme_bw(base_family = base_family) %+replace%
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),   
      axis.text = element_text(size = base_size_a),
      axis.title = element_text(size=base_size_t,face="bold"),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = NA)
    )
}

# Plot
latency.plot <- ggplot(latency.data, aes(cam_type, m)) +
  geom_point() +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2,
                position=position_dodge(.9)) +
  theme_ac1() +
  ylab("Hours before detection") +
  xlab("Camera trap build")

# Save plot
png("latency_plot.png", width = 4.5, height = 4.5, units = 'in', res = 300)
latency.plot
dev.off()

# Kruskal-Wallis test (can't normalise residuals to an acceptable level)
kruskal.test(diff ~ as.factor(cam_type), data = stat.df)

# Pairwise Wilcox test
pairwise.wilcox.test(stat.df$diff, stat.df$cam_type,
                     p.adjust.method = "BH")
