setwd(file.path("~", "Downloads"))

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
library(magick)

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

# Standardised activity-overlap panel plot function
# 
# Draws:
# main plot, rug plot
# Right margin, top to bottom: legend, overlap stats, sample sized
# Bottom right: species silhouette
plot_activity_panel <- function(models,                     # list(u=, o=, t=) overlap::activity fits
                                dhat,                       # list(uo=, ut=, ot=), each vector where [1] = Dhat
                                panel_letter,                # "a)", "b)", "c)" ...
                                line_col = c("black", "red", "blue"),
                                line_lty = c(1, 2, 3),
                                legend_labels = c("Unboxed (U)", "Open box (O)", "Tube box (T)"),
                                comparison_labels = c("U-O", "U-T", "O-T"),
                                line_spacing_mult = 1.8,     # line pitch, as a multiple of one text line's natural height
                                letter_x_offset = -2.75,     # in x-axis units (hours), same scale for all panels
                                letter_y_frac = 0.03,        # letter offset above plot, as fraction of y-range
                                legend_y_intersp = 1.4) {    # vertical spacing between legend rows
  
  fmt_d <- function(d) format(round(d, 2), nsmall = 2)
  
  layout(matrix(c(1, 2), ncol = 1), heights = c(4.5, 1))
  par(mar = c(4.1, 4.1, 2, 10.1), xpd = TRUE, bg = "white", family = "sans")
  
  plot(models$u, yunit = "density", data = "none", las = 1, lwd = 2, lty = line_lty[1],
       tline = list(lwd = 2, lty = line_lty[1]),
       cline = list(lty = 0))
  plot(models$o, yunit = "density", data = "none", add = TRUE,
       tline = list(col = line_col[2], lwd = 2, lty = line_lty[2]),
       cline = list(lty = 0))
  plot(models$t, yunit = "density", data = "none", add = TRUE,
       tline = list(col = line_col[3], lwd = 2, lty = line_lty[3]),
       cline = list(lty = 0))
  
  corners <- par("usr") 
  
  # Legend
  leg <- legend(x = corners[2], y = corners[4], xjust = 0, yjust = 1,
                legend_labels, col = c(1, 2, 4), lty = line_lty, lwd = 2, bty = "n",
                y.intersp = legend_y_intersp)
  legend_bottom <- corners[4] - leg$rect$h
  
  line_height_actual <- graphics::strheight("Xg", cex = par("cex")) * line_spacing_mult
  
  overlap_top <- legend_bottom - line_height_actual   # small gap below the legend
  y_at <- function(i) overlap_top - (i - 1) * line_height_actual
  
  # Overlap block
  text(x = corners[2], y = y_at(1), pos = 4, bquote(bold("Overlap (" ~ Delta ~ ")")), xpd = TRUE)
  text(x = corners[2], y = y_at(2), pos = 4,
       bquote(.(comparison_labels[1]) == .(fmt_d(dhat$uo[1]))), xpd = TRUE)
  text(x = corners[2], y = y_at(3), pos = 4,
       bquote(.(comparison_labels[2]) == .(fmt_d(dhat$ut[1]))), xpd = TRUE)
  text(x = corners[2], y = y_at(4), pos = 4,
       bquote(.(comparison_labels[3]) == .(fmt_d(dhat$ot[1]))), xpd = TRUE)
  
  # line 5 left blank as a section gap
  text(x = corners[2], y = y_at(6), pos = 4, bquote(bold("Sample size (n)")), xpd = TRUE)
  text(x = corners[2], y = y_at(7), pos = 4, bquote("U" == .(length(models$u@data))), xpd = TRUE)
  text(x = corners[2], y = y_at(8), pos = 4, bquote("O" == .(length(models$o@data))), xpd = TRUE)
  text(x = corners[2], y = y_at(9), pos = 4, bquote("T" == .(length(models$t@data))), xpd = TRUE)
  
  text(x = corners[1] + letter_x_offset,
       y = corners[4] + letter_y_frac * (corners[4] - corners[3]),
       pos = 4, bquote(bold(.(panel_letter))), xpd = TRUE, cex = 1.2)
}

add_bottom_right_silhouette <- function(silhouette_png, width_frac = 0.9, x_pad_frac = 0) {
  # Anchors a silhouette to the bottom-right corner of the most recent panel
  #
  # x_pad_frac = fraction of the right-margin width to leave as empty space
  # between the silhouette and the true right-hand edge (0 = flush against
  # the edge; increase to shunt it left off the edge, e.g. 0.05).
  corners <- par("usr")             # xmin, xmax, ymin, ymax of the active panel
  pin <- par("pin")                 # plot region size, inches (w, h)
  right_margin_in <- par("mai")[4]  # right margin width of the active panel, inches
  
  img <- png::readPNG(silhouette_png)
  aspect <- dim(img)[1] / dim(img)[2]   # height/width in pixels, to preserve on resize
  
  width_in  <- right_margin_in * width_frac
  height_in <- width_in * aspect
  
  # convert inches to this panel's user units
  in_to_user_x <- (corners[2] - corners[1]) / pin[1]
  in_to_user_y <- (corners[4] - corners[3]) / pin[2]
  
  width_user  <- width_in  * in_to_user_x
  height_user <- height_in * in_to_user_y
  
  margin_user   <- right_margin_in * in_to_user_x
  right_edge    <- corners[2] + margin_user
  xright        <- right_edge - (right_margin_in * x_pad_frac * in_to_user_x)
  xleft         <- xright - width_user
  
  rasterImage(img,
              xleft   = xleft,
              ybottom = corners[3],
              xright  = xright,
              ytop    = corners[3] + height_user,
              xpd = TRUE)
}


save_activity_panel <- function(panel_expr, rug_expr, silhouette_png, filename,
                                width = 8.5, height = 4.5, res = 300,
                                silhouette_width_frac = 0.75) {
  png(filename, width = width, height = height, units = "in", res = res, bg = "white")
  panel_expr()
  rug_expr()
  add_bottom_right_silhouette(silhouette_png, width_frac = silhouette_width_frac)
  dev.off()
}

# Mouse panel
save_activity_panel(
  panel_expr = function() plot_activity_panel(
    models = list(u = mact.u, o = mact.o, t = mact.t),
    dhat   = list(uo = mouse.uo, ut = mouse.ut, ot = mouse.ot),
    panel_letter = "a)"
  ),
  rug_expr = function() plot_rugs(mact.u, mact.o, mact.t),
  silhouette_png = "mouse.png",
  filename = "mouse_activity_proof.png"
)

# Vole panel
save_activity_panel(
  panel_expr = function() plot_activity_panel(
    models = list(u = voct.u, o = voct.o, t = voct.t),
    dhat   = list(uo = vole.uo, ut = vole.ut, ot = vole.ot),
    panel_letter = "c)"
  ),
  rug_expr = function() plot_rugs(voct.u, voct.o, voct.t),
  silhouette_width_frac = 0.55,
  silhouette_png = "vole.png",
  filename = "vole_activity_proof.png"
)

# Shrew panel
save_activity_panel(
  panel_expr = function() plot_activity_panel(
    models = list(u = shct.u, o = shct.o, t = shct.t),
    dhat   = list(uo = shrew.uo, ut = shrew.ut, ot = shrew.ot),
    panel_letter = "b)"
  ),
  rug_expr = function() plot_rugs(shct.u, shct.o, shct.t),
  silhouette_width_frac = 0.55,
  silhouette_png = "shrew.png",
  filename = "shrew_activity_proof.png"
)


# Combine into a single a-b-c stacked figure
panels <- image_read(c(
  "mouse_activity_proof.png",   # a)
  "shrew_activity.png",         # b)
  "vole_activity.png"           # c)
))

combined <- image_append(panels, stack = TRUE)
image_write(combined, path = "activity_panels_combined.png")

