require(plyr)
require(dplyr)
require(ggplot2)
require(lomb)
require(RobPer)
require(reshape2)

diel_phys_dat = read.table("~/martiny_diel_phys_dat/moisture.tep.res.rin.txt", header = T, sep = "\t")

diel_phys_dat$fixtime = strptime(as.character(diel_phys_dat$Time.bead.beat), '%m/%d/%Y %H:%M')
diel_phys_dat$temptime = strptime(as.character(diel_phys_dat$Temp.time), '%m/%d/%Y %H:%M')


#########################################
###########Lomb-Scargle###########

#data.frames col1 is time col 2 is obs

diel_temp.lomb_dat = data.frame(hrs = diel_phys_dat$hours.cumulative.temp, temp = diel_phys_dat$Temperature)
diel_moist.lomb_dat = data.frame(hrs = diel_phys_dat$hours.cumulative.temp, moist = diel_phys_dat$Moisture)
diel_resp.lomb_dat = data.frame(hrs = diel_phys_dat$hours.cumulative.temp, moist = diel_phys_dat$resp.per.dry.weight)
diel_cell.lomb_dat = data.frame(hrs = diel_phys_dat$hours.cumulative.temp, cell = diel_phys_dat$cell.counts)#cell counts is per dry weight
diel_resp_cell.lomb_dat = data.frame(hrs = diel_phys_dat$hours.cumulative.temp, cell = diel_phys_dat$resp.per.cell)

#By plots
diel_temp.P1.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P1")$hours.cumulative.temp, temp = subset(diel_phys_dat, Plot == "P1")$Temperature)
diel_moist.P1.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P1")$hours.cumulative.temp, moist = subset(diel_phys_dat, Plot == "P1")$Moisture)
diel_resp.P1.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P1")$hours.cumulative.temp, moist = subset(diel_phys_dat, Plot == "P1")$resp.per.dry.weight)
diel_cell.P1.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P1")$hours.cumulative.temp, cell = subset(diel_phys_dat, Plot == "P1")$cell.counts)
diel_resp_cell.P1.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P1")$hours.cumulative.temp, cell = subset(diel_phys_dat, Plot == "P1")$resp.per.cell)

diel_temp.P2.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P2")$hours.cumulative.temp, temp = subset(diel_phys_dat, Plot == "P2")$Temperature)
diel_moist.P2.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P2")$hours.cumulative.temp, moist = subset(diel_phys_dat, Plot == "P2")$Moisture)
diel_resp.P2.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P2")$hours.cumulative.temp, moist = subset(diel_phys_dat, Plot == "P2")$resp.per.dry.weight)
diel_cell.P2.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P2")$hours.cumulative.temp, cell = subset(diel_phys_dat, Plot == "P2")$cell.counts)
diel_resp_cell.P2.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P2")$hours.cumulative.temp, cell = subset(diel_phys_dat, Plot == "P2")$resp.per.cell)

diel_temp.P3.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P3")$hours.cumulative.temp, temp = subset(diel_phys_dat, Plot == "P3")$Temperature)
diel_moist.P3.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P3")$hours.cumulative.temp, moist = subset(diel_phys_dat, Plot == "P3")$Moisture)
diel_resp.P3.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P3")$hours.cumulative.temp, moist = subset(diel_phys_dat, Plot == "P3")$resp.per.dry.weight)
diel_cell.P3.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P3")$hours.cumulative.temp, cell = subset(diel_phys_dat, Plot == "P3")$cell.counts)
diel_resp_cell.P3.lomb_dat = data.frame(hrs = subset(diel_phys_dat, Plot == "P3")$hours.cumulative.temp, cell = subset(diel_phys_dat, Plot == "P3")$resp.per.cell)

#overall
diel_temp.lsp = lsp(diel_temp.lomb_dat, type = "period", from = 22, to = 26)
diel_moist.lsp = lsp(diel_moist.lomb_dat, type = "period", from = 22, to = 26)
diel_resp.lsp = lsp(diel_resp.lomb_dat, type = "period", from = 23, to = 25)
diel_cell.lsp = lsp(diel_cell.lomb_dat, type = "period", from = 22, to = 26)
diel_resp_cell.lsp = lsp(diel_cell.lomb_dat, type = "period", from = 22, to = 26)

#by plots
diel_temp.P1.lsp = lsp(diel_temp.P1.lomb_dat, type = "period", from = 22, to = 26)
diel_moist.P1.lsp = lsp(diel_moist.P1.lomb_dat, type = "period", from = 22, to = 26)
diel_resp.P1.lsp = lsp(diel_resp.P1.lomb_dat, type = "period", from = 22, to = 26)
diel_cell.P1.lsp = lsp(diel_cell.P1.lomb_dat, type = "period", from = 22, to = 26)
diel_resp_cell.P1.lsp = lsp(diel_cell.P1.lomb_dat, type = "period", from = 22, to = 26)

diel_temp.P2.lsp = lsp(diel_temp.P2.lomb_dat, type = "period", from = 22, to = 26)
diel_moist.P2.lsp = lsp(diel_moist.P2.lomb_dat, type = "period", from = 22, to = 26)
diel_resp.P2.lsp = lsp(diel_resp.P2.lomb_dat, type = "period", from = 22, to = 26)
diel_cell.P2.lsp = lsp(diel_cell.P2.lomb_dat, type = "period", from = 22, to = 26)
diel_resp_cell.P2.lsp = lsp(diel_cell.P2.lomb_dat, type = "period", from = 22, to = 26)

diel_temp.P3.lsp = lsp(diel_temp.P3.lomb_dat, type = "period", from = 22, to = 26)
diel_moist.P3.lsp = lsp(diel_moist.P3.lomb_dat, type = "period", from = 22, to = 26)
diel_resp.P3.lsp = lsp(diel_resp.P3.lomb_dat, type = "period", from = 22, to = 26)
diel_cell.P3.lsp = lsp(diel_cell.P3.lomb_dat, type = "period", from = 22, to = 26)
diel_resp_cell.P3.lsp = lsp(diel_cell.P3.lomb_dat, type = "period", from = 22, to = 26)

#overall
diel_temp.RobPer = RobPer(ts = diel_temp.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_moist.RobPer = RobPer(ts = diel_moist.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_resp.RobPer = RobPer(ts = diel_resp.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_cell.RobPer = RobPer(ts = diel_cell.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_resp_cell.RobPer = RobPer(ts = diel_resp_cell.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

#by plots

diel_temp.P1.RobPer = RobPer(ts = diel_temp.P1.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_moist.P1.RobPer = RobPer(ts = diel_moist.P1.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_resp.P1.RobPer = RobPer(ts = diel_resp.P1.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_cell.P1.RobPer = RobPer(ts = diel_cell.P1.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_resp_cell.P1.RobPer = RobPer(ts = diel_resp_cell.P1.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_temp.P2.RobPer = RobPer(ts = diel_temp.P2.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_moist.P2.RobPer = RobPer(ts = diel_moist.P2.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_resp.P2.RobPer = RobPer(ts = diel_resp.P2.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_cell.P2.RobPer = RobPer(ts = diel_cell.P2.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_resp_cell.P2.RobPer = RobPer(ts = diel_resp_cell.P2.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_temp.P3.RobPer = RobPer(ts = diel_temp.P3.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_moist.P3.RobPer = RobPer(ts = diel_moist.P3.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_resp.P3.RobPer = RobPer(ts = diel_resp.P3.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_cell.P3.RobPer = RobPer(ts = diel_cell.P3.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

diel_resp_cell.P3.RobPer = RobPer(ts = diel_resp_cell.P3.lomb_dat,
weighting = F,
periods = seq(from = 12, to = 36, by = 0.5),
regression = "L2",
model = "sine")

temp_moist_rob_per_df = data.frame(period = seq(from = 12, to = 36, by = 0.5),
    temp_overall = diel_temp.RobPer,
    moisture_overall = diel_moist.RobPer,
    resp_overall = diel_resp.RobPer,
cell_overall = diel_cell.RobPer,
    temp_P1 = diel_temp.P1.RobPer,
    moisture_P1 = diel_moist.P1.RobPer,
resp_P1 = diel_resp.P1.RobPer,
cell_P1 = diel_cell.P1.RobPer,
   temp_P2 = diel_temp.P2.RobPer,
    moisture_P2 = diel_moist.P2.RobPer,
resp_P2 = diel_resp.P2.RobPer,
cell_P2 = diel_cell.P2.RobPer,
   temp_P3 = diel_temp.P3.RobPer,
    moisture_P3 = diel_moist.P3.RobPer,
resp_P3 = diel_resp.P3.RobPer,
cell_P3 = diel_cell.P3.RobPer
)

temp_moist_robper_melt = melt(temp_moist_rob_per_df, id.vars = "period")

pdf("martiny_diel_phys_dat/plots/robper_analysis.pdf")
ggplot(temp_moist_robper_melt, aes(period, value)) +
geom_point() +
facet_wrap(~variable, ncol = 4) +
ylab(expression(paste("r"^2))) +
xlab("Period (hrs)") +
my_gg_theme
dev.off()




############################################################
#DIRECT regression of respiration against temp/moist
############################################################

moist_resp_lm = lme(fixed = log10(resp.per.dry.weight+0.1) ~ Moisture, random = ~1|Plot/Temp.time, data = diel_phys_dat,
na.action = na.omit) #needs log
moist_resp_lm.poly = lme(fixed = log10(resp.per.dry.weight+0.1) ~ poly(Moisture,2), random = ~1|Plot/Temp.time, data = diel_phys_dat, na.action = na.omit) #needs log trans
anova(moist_resp_lm,moist_resp_lm.poly)


plot(residuals(moist_resp_lm.poly))
qqnorm(residuals(moist_resp_lm.poly))
summary(moist_resp_lm.poly)

fit = lm(log10(resp.per.dry.weight+0.1) ~ poly(Moisture,2), data = diel_phys_dat)


pdf("diel_phys_data_moisture_vs_respiration.pdf")
ggplot(diel_phys_dat, aes(Moisture, resp.per.dry.weight+0.1)) +
geom_point() +
geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
scale_y_log10() +
my_gg_theme +
labs(y = expression(paste("Respiration (",mu,"g CO"[2], "-C minute"^-1," g"^-1," litter)")) ,
title = paste("R2 = ", signif(summary(fit)$adj.r.squared, 2),
" P < 0.001")
)
dev.off()

temp_resp_lm = lme(fixed = log10(resp.per.dry.weight+0.1) ~ Temperature, random = ~1|Plot/Temp.time, data = diel_phys_dat,
na.action = na.omit) #needs log trans
plot(residuals(temp_resp_lm))
qqnorm(residuals(temp_resp_lm))
summary(temp_resp_lm)
anova(temp_resp_lm)


fit = lm(log10(resp.per.dry.weight+0.1) ~ Temperature, data = diel_phys_dat)

pdf("diel_phys_data_temp_vs_respiration.pdf")
ggplot(diel_phys_dat, aes(Temperature, resp.per.dry.weight+0.1)) +
geom_point() +
geom_smooth(method = "lm") +
scale_y_log10() +
my_gg_theme +
labs(y = expression(paste("Respiration (",mu,"g CO"[2], "-C minute"^-1," g"^-1," litter)")) ,
title = paste("R2 = ", signif(summary(fit)$adj.r.squared, 2),
" P = 0 .042")
)
dev.off()

#####################
#Temp moisture vs log10(cell.counts)


moist_resp_lm = lme(fixed = log10(cell.counts) ~ Moisture, random = ~1|Plot/Temp.time, data = diel_phys_dat,
na.action = na.omit) #needs log
moist_resp_lm.poly = lme(fixed = log10(cell.counts) ~ poly(Moisture,2), random = ~1|Plot/Temp.time, data = diel_phys_dat, na.action = na.omit) #needs log trans
anova(moist_resp_lm,moist_resp_lm.poly)


plot(residuals(moist_resp_lm))
qqnorm(residuals(moist_resp_lm))
summary(moist_resp_lm)
anova(moist_resp_lm)


fit = lm(log10(cell.counts) ~ Moisture, data = diel_phys_dat)
fit.poly = lm(log10(cell.counts) ~ poly(Moisture,2), data = diel_phys_dat)
anova(fit, fit.poly)

pdf("diel_phys_data_moisture_vs_cells.pdf")
ggplot(diel_phys_dat, aes(Moisture, cell.counts)) +
geom_point() +
geom_smooth(method = "lm") +
#scale_y_log10() +
my_gg_theme +
labs(y = expression(paste("cells g"^-1," litter")) ,
title = paste("R2 = ", signif(summary(fit)$adj.r.squared, 2),
" P < 0.001")
)
dev.off()

temp_resp_lm = lme(fixed = log10(cell.counts) ~ Temperature, random = ~1|Plot/Temp.time, data = diel_phys_dat,
na.action = na.omit) #needs log trans
plot(residuals(temp_resp_lm))
qqnorm(residuals(temp_resp_lm))
summary(temp_resp_lm)
anova(temp_resp_lm)


fit = lm(log10(cell.counts) ~ Temperature, data = diel_phys_dat)

pdf("diel_phys_data_temp_vs_cells.pdf")
ggplot(diel_phys_dat, aes(Temperature, cell.counts)) +
geom_point() +
#geom_smooth(method = "lm") +
#scale_y_log10() +
my_gg_theme +
labs(y = expression(paste("cells g"^-1," litter")) ,
title = paste("R2 = ", signif(summary(fit)$r.squared, 2),
" P = 0 .47")
)
dev.off()

################
#Cells vs resp#
###############

cells_resp_lm = lme(fixed = log10(resp.per.dry.weight+0.1) ~ log10(cell.counts), random = ~1|Plot/Temp.time, data = diel_phys_dat,
na.action = na.omit) #needs log
cells_resp_lm.poly = lme(fixed = log10(resp.per.dry.weight+0.1) ~ poly(log10(cell.counts),2), random = ~1|Plot/Temp.time, data = diel_phys_dat, na.action = na.omit) #needs log trans
anova(cells_resp_lm,cells_resp_lm.poly)


cells_resp_lm = lme(fixed = log10(resp.per.dry.weight+0.1) ~ (cell.counts), random = ~1|Plot/Temp.time, data = diel_phys_dat,
na.action = na.omit) #needs log

plot(residuals(cells_resp_lm))
qqnorm(residuals(cells_resp_lm))
summary(cells_resp_lm)
anova(cells_resp_lm)

fit = lm(log10(resp.per.dry.weight+0.1) ~ log10(cell.counts), data = diel_phys_dat)

fit = lm(log10(resp.per.dry.weight+0.1) ~ (cell.counts), data = diel_phys_dat)


pdf("diel_phys_data_resp_vs_cells.pdf")
ggplot(diel_phys_dat, aes(cell.counts, resp.per.dry.weight+0.1)) +
geom_point() +
geom_smooth(method = "lm") +
scale_y_log10() +
scale_x_log10() +
my_gg_theme +
labs(x = expression(paste("cells g"^-1," litter")) ,
y = expression(paste("Respiration (",mu,"g CO"[2], "-C minute"^-1," g"^-1," litter)")) ,
title = paste("R2 = ", signif(summary(fit)$adj.r.squared, 2),
" P < 0.001")
)
dev.off()

model.tables(aov((cell.counts) ~ Plot, data =diel_phys_dat), "means")
summary(aov((cell.counts) ~ Plot, data =diel_phys_dat))
TukeyHSD(aov((cell.counts) ~ Plot, data =diel_phys_dat))


model.tables(aov(Moisture ~ Plot, data =diel_phys_dat), "means")
summary(aov(Moisture ~ Plot, data =diel_phys_dat))
TukeyHSD(aov(Moisture ~ Plot, data =diel_phys_dat))


temp_resp_gt_50_moist_lm = lm(resp.per.dry.weight+0.1 ~ Temperature, data = subset(diel_phys_dat, Moisture > 0.5)) #needs log trans
plot(residuals(temp_resp_gt_50_moist_lm))
qqnorm(residuals(temp_resp_gt_50_moist_lm))
summary(temp_resp_gt_50_moist_lm)

moist_temp_resp_lm = lm(log10(resp.per.dry.weight+0.1) ~ Moisture*Temperature, data = diel_phys_dat) #needs log trans
plot(residuals(moist_temp_resp_lm))
qqnorm(residuals(moist_temp_resp_lm))
summary(moist_temp_resp_lm)

moist_temp_resp_plot_lm_1 = lm(log10(resp.per.dry.weight+0.1) ~ Moisture+Temperature+Temperature:Plot+Moisture:Plot+Plot + hours.cumulative.temp, data = diel_phys_dat) #needs log trans
moist_temp_resp_plot_lm_2 = lm(log10(resp.per.dry.weight+0.1) ~ Moisture+Temperature+Temperature:Plot+Plot + hours.cumulative.temp, data = diel_phys_dat) #no sig dif from 1. use this BEST MODEL
moist_temp_resp_plot_lm_3 = lm(log10(resp.per.dry.weight+0.1) ~ Moisture+Temperature+Moisture:Plot+Plot + hours.cumulative.temp, data = diel_phys_dat) #1 sig better than 3
moist_temp_resp_plot_lm_4 = lm(log10(resp.per.dry.weight+0.1) ~ Moisture+Temperature+Temperature:Plot + hours.cumulative.temp, data = diel_phys_dat) #no sig dif from 1. use this
moist_temp_resp_plot_lm_5 = lm(log10(resp.per.dry.weight+0.1) ~ Moisture+I(Moisture^2)+Temperature+Temperature:Plot + hours.cumulative.temp, data = diel_phys_dat) #ns than 2


anova(moist_resp_lm, moist_temp_resp_plot_lm_2) #sig
anova(temp_resp_lm,moist_temp_resp_plot_lm_2) #sig
anova(moist_temp_resp_lm,moist_temp_resp_plot_lm_2) #sig
anova(moist_temp_resp_plot_lm_4,moist_temp_resp_plot_lm_2) #sig
anova(moist_temp_resp_plot_lm_2, moist_temp_resp_plot_lm_5) #2 better

anova(moist_temp_resp_plot_lm_2)

moist_temp_resp_plot_aov = aov(log10(resp.per.dry.weight+0.1) ~ Moisture+Temperature+Temperature:Plot+Plot + hours.cumulative.temp, data = diel_phys_dat, na.action = na.exclude)

plot(log10(resp.per.dry.weight+0.1) ~ Moisture, data = diel_phys_dat)
lines(fitted(moist_sin_plot_time_plot_lm) ~Moisture, data= diel_phys_dat)





############################################################
#Harmonic regressions LM

#TEMPERATURE
#sin fit overall temp

temp_sin_lm = lm(Temperature ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat) # plot is not sig
plot(residuals(temp_sin_lm))
qqnorm(residuals(temp_sin_lm))
summary(temp_sin_lm)

#sin fit + plot nesting
temp_sin_plot_lm = lm(Temperature ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat)

#sin fit + time
temp_sin_time_lm = lm(Temperature ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24)  + hours.cumulative.temp, data = diel_phys_dat)
plot(residuals(temp_sin_time_lm))
qqnorm(residuals(temp_sin_time_lm))
#sin fit + plot
temp_sin_plotEffect_lm = lm(Temperature ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24)  + Plot, data = diel_phys_dat)
plot(residuals(temp_sin_time_lm))
qqnorm(residuals(temp_sin_time_lm))

anova(temp_sin_lm,temp_sin_plot_lm)#marg, P = 0.059
anova(temp_sin_lm,temp_sin_time_lm) #sig
anova(temp_sin_lm,temp_sin_plotEffect_lm) #ns

anova(temp_sin_time_lm)

summary(temp_sin_time_lm)


pdf("martiny_diel_phys_dat/plots/temp_sine.pdf")
plot(Temperature ~ hours.cumulative.temp, data = diel_phys_dat)
lines(fitted(temp_sin_time_lm) ~ hours.cumulative.temp, data = diel_phys_dat, col = "blue")
dev.off()

#To find amplitude according to (https://stats.stackexchange.com/questions/77543/how-do-i-get-the-amplitude-and-phase-for-sine-wave-from-lm-summary)
b0 <- coef(temp_sin_lm)[1]
alpha <- coef(temp_sin_lm)[2]
beta <- coef(temp_sin_lm)[3]

r <- sqrt(alpha^2 + beta^2)
phi <- atan2(beta, alpha)


#MOISTURE

#sin fit overall moist
#basic harmonic model

moist_sin_lm = lm(log10(Moisture*100) ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
plot(residuals(moist_sin_lm))
qqnorm(residuals(moist_sin_lm))
summary(moist_sin_lm)

#To find amplitude according to (https://stats.stackexchange.com/questions/77543/how-do-i-get-the-amplitude-and-phase-for-sine-wave-from-lm-summary)
b0 <- coef(moist_sin_lm)[1]
alpha <- coef(moist_sin_lm)[2]
beta <- coef(moist_sin_lm)[3]

#r is amplitude and "offset" (or phase) is phi
r <- sqrt(alpha^2 + beta^2)
phi <- atan2(beta, alpha)


#sin fit moist by plot
moist_sin_plot_lm = lm(log10(Moisture*100) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
#adding time
moist_sin_plot_time_lm = lm(log10(Moisture*100) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude) #sig better model with time included
#adding plot
moist_sin_plot_plot_lm = lm(log10(Moisture*100) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + Plot, data = diel_phys_dat[-35,], na.action = na.exclude) #sig better model with time included
#adding plot and time
moist_sin_plot_time_plot_lm = lm(log10(Moisture*100) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp + Plot, data = diel_phys_dat[-35,], na.action = na.exclude) #sig better model with time

anova(moist_sin_lm,moist_sin_plot_lm) #sig
anova(moist_sin_lm,moist_sin_plot_time_lm) #sig
anova(moist_sin_lm,moist_sin_plot_plot_lm) #sig
anova(moist_sin_lm,moist_sin_plot_time_plot_lm) #sig

anova(moist_sin_plot_time_lm,moist_sin_plot_plot_lm) #sig
anova(moist_sin_plot_plot_lm, moist_sin_plot_time_plot_lm) #sig
anova(moist_sin_plot_time_plot_lm)

moist_sin_plot_time_plot_aov = aov(log10(Moisture*100) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp + Plot, data = diel_phys_dat[-35,], na.action = na.exclude)
summary(moist_sin_plot_time_plot_aov)
model.tables(moist_sin_plot_time_plot_aov, "means")

moist_sin_plot_time_plot_lm

pdf("martiny_diel_phys_dat/plots/moist_sine.pdf")
plot(log10(Moisture) ~ hours.cumulative.temp, data = diel_phys_dat[-35,])
lines(fitted(moist_sin_plot_time_plot_lm) ~ hours.cumulative.temp, data = diel_phys_dat[-35,], col = "blue")
dev.off()

plot(residuals(moist_sin_plot_time_lm))
qqnorm(residuals(moist_sin_plot_time_lm))
summary(moist_sin_plot_time_lm)
anova(moist_sin_plot_time_lm)

anova(moist_sin_plot_time_lm)

b0.p1 <- coef(moist_sin_plot_time_lm)[1]
b0.p2 <- coef(moist_sin_plot_time_lm)[2]
b0.p3 <- coef(moist_sin_plot_time_lm)[3]
alpha.p1 <- coef(moist_sin_plot_time_lm)[5]
alpha.p2 <- coef(moist_sin_plot_time_lm)[6]
alpha.p3 <- coef(moist_sin_plot_time_lm)[7]
beta.p1 <- coef(moist_sin_plot_time_lm)[8]
beta.p2 <- coef(moist_sin_plot_time_lm)[9]
beta.p3 <- coef(moist_sin_plot_time_lm)[10]


r.p1 <- sqrt(alpha.p1^2 + beta.p1^2)
phi.p1 <- atan2(beta.p1, alpha.p1)
r.p2 <- sqrt(alpha.p2^2 + beta.p2^2)
phi.p2 <- atan2(beta.p2, alpha.p2)
r.p3 <- sqrt(alpha.p3^2 + beta.p3^2)
phi.p3 <- atan2(beta.p3, alpha.p3)


#RESPIRATION

#sin fit resp
resp_sin_lm = lm(resp.per.dry.weight ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
plot(residuals(resp_sin_lm))
qqnorm(residuals(resp_sin_lm))
#sin fit resp by plot
resp_sin_plot_lm = lm(resp.per.dry.weight ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
#adding time
resp_sin_plot_time_lm = lm(resp.per.dry.weight ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude)
resp_sin_noPlot_time_lm = lm(resp.per.dry.weight ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude)
resp_plot_lm = lm(resp.per.dry.weight ~ Plot, data = diel_phys_dat[-35,], na.action = na.exclude)

anova(resp_sin_lm, resp_sin_plot_lm)
plot(residuals(resp_sin_plot_lm))
qqnorm(residuals(resp_sin_plot_lm))

anova(resp_sin_lm, resp_sin_plot_lm)
anova(resp_sin_lm, resp_sin_noPlot_time_lm)
anova(resp_sin_plot_lm,resp_sin_plot_time_lm)
anova(resp_sin_plot_lm,resp_sin_noPlot_time_lm)
anova(resp_sin_lm, resp_plot_lm)
anova(resp_sin_plot_lm, resp_plot_lm)

#Time is not adding to model, model with harmonic plus plot is better than plot only

pdf("martiny_diel_phys_dat/plots/resp_sine.pdf")
plot(log10(resp.per.dry.weight+1) ~ hours.cumulative.temp, data = diel_phys_dat[-35,])
lines(fitted(resp_sin_plot_lm) ~ hours.cumulative.temp, data = diel_phys_dat[-35,], col = "blue")
dev.off()


#std error of r = r*sqrt(0.5 * ( (2*alpha.err/alpha)^2 + (2*beta.err/beta)^2 ) )
#errors
alpha.p1.err <- summary(moist_sin_plot_time_lm)$coef[5,2]
alpha.p2.err <- summary(moist_sin_plot_time_lm)$coef[6,2]
alpha.p3.err <- summary(moist_sin_plot_time_lm)$coef[7,2]
beta.p1.err <- summary(moist_sin_plot_time_lm)$coef[8,2]
beta.p2.err <- summary(moist_sin_plot_time_lm)$coef[9,2]
beta.p3.err <- summary(moist_sin_plot_time_lm)$coef[10,2]

r.p1.err <- r.p1 * sqrt(0.5 * ( ((2*alpha.p1.err/alpha)^2) + ((2*beta.p1.err/beta)^2) ))
r.p2.err <- r.p2 * sqrt(0.5 * ( ((2*alpha.p2.err/alpha)^2) + ((2*beta.p2.err/beta)^2) ))
r.p3.err <- r.p3 * sqrt(0.5 * ( ((2*alpha.p3.err/alpha)^2) + ((2*beta.p3.err/beta)^2) ))


t.test2(r.p1, r.p2, r.p1.err, r.p2.err, 35, 35)
t.test2(r.p1, r.p3, r.p1.err, r.p3.err, 35, 35)
t.test2(r.p2, r.p3, r.p2.err, r.p3.err, 35, 35)


#CELL COUNTS

#sin fit cells
cells_sin_lm = lm(log10(cell.counts) ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
plot(residuals(cells_sin_lm))
qqnorm(residuals(cells_sin_lm))
#sin fit cells by plot
cells_sin_plot_lm = lm(log10(cell.counts) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
#adding time
cells_sin_plot_time_lm = lm(log10(cell.counts) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude)
cells_sin_noPlot_time_lm = lm(log10(cell.counts) ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude)
cells_plot_no_harmonic_lm = lm(log10(cell.counts) ~ Plot, data = diel_phys_dat[-35,], na.action = na.exclude)
cells_sin_plot_no_harmonic_aov = aov(log10(cell.counts) ~ Plot, data = diel_phys_dat[-35,], na.action = na.exclude)

anova(cells_sin_lm, cells_sin_plot_lm)
plot(residuals(cells_sin_plot_lm))
qqnorm(residuals(cells_sin_plot_lm))

anova(cells_sin_lm, cells_sin_plot_lm)
anova(cells_sin_lm, cells_sin_noPlot_time_lm)
anova(cells_sin_plot_lm,cells_sin_plot_time_lm)
anova(cells_sin_plot_lm,cells_sin_noPlot_time_lm)
anova(cells_sin_lm, cells_sin_plot_no_harmonic_lm)
anova(cells_sin_plot_lm, cells_plot_no_harmonic_lm)


#Time is not adding to model, plot only fit is best (better than basic harmonic)

anova(cells_sin_plot_lm)
summary(cells_sin_plot_lm)

pdf("martiny_diel_phys_dat/plots/resp_sine.pdf")
plot(log10(resp.per.dry.weight+1) ~ hours.cumulative.temp, data = diel_phys_dat[-35,])
lines(fitted(resp_sin_plot_lm) ~ hours.cumulative.temp, data = diel_phys_dat[-35,], col = "blue")
dev.off()

#RESP PER CELL COUNTS

#sin fit cells
resp_p_cells_sin_lm = lm(resp.per.cell ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
plot(residuals(resp_p_cells_sin_lm))
qqnorm(residuals(resp_p_cells_sin_lm))
#sin fit cells by plot
resp_p_cells_sin_plot_lm = lm(resp.per.cell ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
#adding time
resp_p_cells_sin_plot_time_lm = lm(resp.per.cell ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude)
resp_p_cells_sin_noPlot_time_lm = lm(resp.per.cell ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude)
resp_p_cells_plot_no_harmonic_lm = lm(resp.per.cell ~ Plot, data = diel_phys_dat[-35,], na.action = na.exclude)
resp_p_cells_sin_plot_no_harmonic_aov = aov(resp.per.cell ~ Plot, data = diel_phys_dat[-35,], na.action = na.exclude)

anova(resp_p_cells_sin_lm, resp_p_cells_sin_plot_lm)
plot(residuals(resp_p_cells_sin_plot_lm))
qqnorm(residuals(resp_p_cells_sin_plot_lm))

anova(resp_p_cells_sin_lm, resp_p_cells_sin_plot_lm)
anova(resp_p_cells_sin_lm, resp_p_cells_sin_noPlot_time_lm)
anova(resp_p_cells_sin_plot_lm,resp_p_cells_sin_plot_time_lm)
anova(resp_p_cells_sin_plot_lm,resp_p_cells_sin_noPlot_time_lm)#no time better
anova(resp_p_cells_sin_lm, resp_p_cells_plot_no_harmonic_lm)
anova(resp_p_cells_sin_plot_lm, resp_p_cells_plot_no_harmonic_lm)

summary(resp_p_cells_sin_plot_lm)
resp_p_cells_sin_plot_lm = lm(resp.per.cell ~ sin(2*pi*hours.cumulative.temp/24):Plot + cos(2*pi*hours.cumulative.temp/24):Plot, data = diel_phys_dat[-35,], na.action = na.exclude)

################################################
#TESTING RESIDUALS corelations from lms

#resp_sin_plot_lm
#temp_sin_lm
#moist_sin_plot_time_lm

temp_moist_resid_lm = lm( residuals(moist_sin_plot_time_plot_lm) ~ residuals(temp_sin_lm)[-35] )
plot(residuals(temp_moist_resid_lm))
qqnorm(residuals(temp_moist_resid_lm))
summary(temp_moist_resid_lm)
summary(temp_moist_resid_lm)


temp_resp_resid_lm = lm( log10( residuals(resp_sin_plot_lm) + 0.5) ~ residuals(temp_sin_lm)[-35] )
qqnorm(residuals(temp_resp_resid_lm))
plot(residuals(temp_resp_resid_lm))
summary(temp_resp_resid_lm)

moist_resp_resid_lm = lm(residuals(resp_sin_plot_lm) ~ residuals(moist_sin_plot_time_lm) )
qqnorm(residuals(moist_resp_resid_lm))
plot(residuals(moist_resp_resid_lm))
summary(moist_resp_resid_lm)

moist_temp_resid_lm = lm( residuals(resp_sin_plot_lm) ~ residuals(moist_sin_plot_time_lm) + residuals(moist_sin_plot_time_lm):residuals(temp_sin_lm)[-35])
anova(moist_resp_resid_lm, moist_temp_resid_lm)


##With quadratic
moist_resp_resid_sq_lm = lm(residuals(resp_sin_plot_lm) ~ residuals(moist_sin_plot_time_lm) + I(residuals(moist_sin_plot_time_lm)^2))
qqnorm(residuals(moist_resp_resid_lm))
plot(residuals(moist_resp_resid_lm))
summary(moist_resp_resid_sq_lm)
anova(moist_resp_resid_lm,moist_resp_resid_sq_lm)

residuals.df = data.frame(Moisture.residuals = residuals(moist_sin_plot_time_lm), Respiration.residuals = residuals(resp_sin_plot_lm), Temperature.residuals = residuals(temp_sin_lm)[-35] )

pdf("martiny_diel_phys_dat/plots/moisture_resp_residuals_fit.pdf")
ggplot(residuals.df, aes(Moisture.residuals, Respiration.residuals)) +
geom_point() +
stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
my_gg_theme +
labs(x = "Moisture periodic model residuals", y = "Respiration periodic model residuals")
dev.off()









###############
#####PLOTS#####
###############
pdf("martiny_diel_phys_dat/plots/moisture_v_time.pdf")
ggplot(diel_phys_dat, aes(temptime, Moisture*100, color = Plot)) +
geom_point() +
scale_color_manual(values = cbPalette) +
ylab("Gravimetric moisture content (%)") +
xlab("Time") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/moisture_v_time_w_fit_log10.pdf", width = 5, height = 4)
ggplot(diel_phys_dat, aes(hours.cumulative.temp, (Moisture*100), color = Plot)) +
geom_point() +
stat_smooth(method = "lm", se = F, formula = y~sin(2*pi*x/24)+ cos(2*pi*x/24) + x) + #
scale_color_manual(values = cbPalette, labels = c("P1" = "Plot 1", "P2" = "Plot 2", "P3" = "Plot 3")) +
scale_y_log10(breaks = c(1,2,5,10,20,50,100,200), limits = c(1,200)) +
ylab("Moisture content (%)") +
xlab("Time (hrs)") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/moisture_v_time_w_fit_no_log10.pdf")
ggplot(diel_phys_dat, aes(hours.cumulative.temp, (Moisture*100), color = Plot)) +
geom_point() +
stat_smooth(method = "lm", se = F, formula = y~sin(2*pi*x/24)+ cos(2*pi*x/24) + x) + #
scale_color_manual(values = cbPalette, labels = c("P1" = "Plot 1", "P2" = "Plot 2", "P3" = "Plot 3")) +
ylab("Gravimetric moisture content") +
xlab("Time") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/temp_v_time.pdf")
ggplot(diel_phys_dat, aes(temptime, Temperature, color = Plot)) +
geom_point() +
scale_color_manual(values = cbPalette) +
xlab("Time") +
ylab(expression("Temperature " ( degree*C))) +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/temp_v_time_w_fit.pdf")
ggplot(diel_phys_dat, aes(hours.cumulative.temp, Temperature, color = Plot)) +
geom_point() +
stat_smooth(method = "lm", se = F, formula = y~sin(2*pi*x/24) + cos(2*pi*x/24)) +
scale_color_manual(values = cbPalette) +
xlab("Time") +
ylab(expression("Temperature " ( degree*C))) +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/temp_v_time_w_fit_constant_phase.pdf", width = 5, height = 4)
ggplot() +
geom_point(data = diel_phys_dat, aes(x = hours.cumulative.temp, y = Temperature, color = Plot)) +
stat_smooth(data = subset(diel_phys_dat, Plot == "P1"),
    aes(x = hours.cumulative.temp+0.6, y = Temperature),
    method = "lm",
    se = F,
    color = "#E69F00",
    formula = y~sin(2*pi*x/24) + cos(2*pi*x/24)) +
stat_smooth(data = subset(diel_phys_dat, Plot == "P2"),
aes(x = hours.cumulative.temp, y = Temperature),
method = "lm",
se = F,
color = "#56B4E9",
formula = y~sin(2*pi*x/24) + cos(2*pi*x/24)) +
stat_smooth(data = subset(diel_phys_dat, Plot == "P3"),
aes(x = hours.cumulative.temp-0.25, y = Temperature),
method = "lm",
se = F,
color = "#009E73",
formula = y~sin(2*pi*x/24) + cos(2*pi*x/24)) +
scale_color_manual(values = cbPalette, labels = c("P1" = "Plot 1", "P2" = "Plot 2", "P3" = "Plot 3")) +
xlab("Time (hrs)") +
ylab(expression("Temperature " ( degree*C))) +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/temp_v_time_w_fit_no_plot.pdf")
ggplot(diel_phys_dat, aes(hours.cumulative.temp, Temperature)) +
geom_point(aes(color = Plot)) +
stat_smooth(method = "lm", color = "black", se = F, formula = y~sin(2*pi*x/24) + cos(2*pi*x/24)) +
scale_color_manual(values = cbPalette) +
xlab("Time") +
ylab(expression("Temperature " ( degree*C))) +
my_gg_theme
dev.off()


pdf("martiny_diel_phys_dat/plots/resp_v_time.pdf", width = 5, height = 4)
ggplot(diel_phys_dat, aes(temptime, resp.per.dry.weight, color = Plot)) +
geom_point() +
scale_color_manual(values = cbPalette, labels = c("P1" = "Plot 1", "P2" = "Plot 2", "P3" = "Plot 3")) +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " g"^-1, "litter" ))) +
xlab("Time (hrs)") +
my_gg_theme
dev.off()

ggplot(diel_phys_dat, aes(temptime, resp.per.log10.cell, color = Plot)) +
geom_point() +
scale_color_manual(values = cbPalette) +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " bacterial cell"^-1 ))) +
xlab("Time") +
my_gg_theme

ggplot(diel_phys_dat, aes(temptime, resp.per.cell, color = Plot)) +
geom_point() +
scale_color_manual(values = cbPalette) +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " g"^-1, " bacterial cell" ))) +
xlab("Time") +
my_gg_theme

pdf("martiny_diel_phys_dat/plots/resp_v_time_w_fit.pdf", width = 5, height = 4)
ggplot(diel_phys_dat, aes(hours.cumulative.temp, (resp.per.dry.weight), color = Plot)) + #+0.0729325
geom_point() +
stat_smooth(method = "lm", se = F, formula = y~sin(2*pi*x/24) + cos(2*pi*x/24)) +
scale_color_manual(values = cbPalette, labels = c("P1" = "Plot 1", "P2" = "Plot 2", "P3" = "Plot 3")) +
#scale_y_log10(breaks = c(0,0.01,0.1,1,3)) +
xlab("Time (hrs)") +
ylab(expression(paste(mu, "g CO"[2],"-C min"^-1, " g"^-1, " litter" ))) +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/cells_v_time.pdf")
ggplot(diel_phys_dat, aes(temptime, log10(cell.counts), color = Plot)) +
geom_point() +
scale_color_manual(values = cbPalette) +
ylab("Bacterial cell count") +
xlab("Time") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/resp_per_cell_v_time_w_fit.pdf", height = 4, width = 5)
ggplot(diel_phys_dat, aes(hours.cumulative.temp, resp.per.cell, color = Plot)) +
geom_point() +
stat_smooth(method = "lm", se = F, formula = y~sin(2*pi*x/24) + cos(2*pi*x/24)) +
scale_color_manual(values = cbPalette, labels = c("P1" = "Plot 1", "P2" = "Plot 2", "P3" = "Plot 3")) +
xlab("Time (hrs)") +
ylab(expression(paste(mu, "g CO"[2],"-C min"^-1, " bacterial cell"^-1))) +
scale_y_continuous(labels = fancy_scientific)+
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/cells_v_time_w_fit.pdf")
ggplot(diel_phys_dat, aes(Plot, log10(cell.counts), color = Plot)) +
geom_boxplot() +
scale_color_manual(values = cbPalette, labels = c("Plot 1", "Plot 2", "Plot 3")) +
xlab("") +
ylab("Bacterial cell count") +
my_gg_theme +
theme(
#legend.position = "none"
)
dev.off()

########################################################
pdf("martiny_diel_phys_dat/plots/resp_v_moisture.pdf")
ggplot(diel_phys_dat, aes(Moisture, log10(resp.per.dry.weight+0.1))) +
geom_point(aes(color = Plot)) +
scale_color_manual(values = cbPalette) +
stat_smooth(method = "lm", se = F, color = "black", formula = y~poly(x,2)) +
ylab(expression(paste("log10 ", mu, "g CO"[2],"-C minute"^-1, " g"^-1, " dry litter" ))) +
xlab("Gravimetric moisture content") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/resp_v_moisture.pdf")
ggplot(diel_phys_dat, aes(Moisture, resp.per.dry.weight)) +
geom_point(aes(color = Plot)) +
scale_color_manual(values = cbPalette) +
stat_smooth(method = "lm", se = F, color = "black", formula = y~poly(x,2)) +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " g"^-1, " dry litter" ))) +
xlab("Gravimetric moisture content") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/resp_v_moisture_temp.pdf")
ggplot(diel_phys_dat, aes(Moisture*100, log10(resp.per.dry.weight+0.1), color = Temperature)) +
geom_point() +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
stat_smooth(method = "lm", se = F, color = "black") +
ylab(expression(paste("log10 ", mu, "g CO"[2],"-C minute"^-1, " g"^-1, " dry litter" ))) +
xlab("Gravimetric moisture content (%)") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/resp_v_temp.pdf")
ggplot(diel_phys_dat, aes(Temperature, log10(resp.per.dry.weight+0.1), color = Plot)) +
geom_point() +
scale_color_manual(values = cbPalette) +
stat_smooth(method = "lm", se = F) +
ylab(expression(paste("log10 ", mu, "g CO"[2],"-C minute"^-1, " g"^-1, " dry litter" ))) +
xlab(expression("Temperature " ( degree*C))) +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/resp_v_temp_moisture_gt_50.pdf")
#dplyr filter requires POSIXct instread of lt
ggplot(subset(diel_phys_dat, Moisture > 0.5), aes(Temperature, resp.per.dry.weight, color = Plot)) +
geom_point(size = 5) +
scale_color_manual(values = cbPalette) +
stat_smooth(method = "lm", se = F, color = "black") +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " g"^-1, " dry litter" ))) +
xlab(expression("Temperature " ( degree*C))) +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/plots/moisture_v_temp.pdf")
ggplot(diel_phys_dat, aes(Temperature, Moisture, color = Plot)) +
geom_point() +
scale_color_manual(values = cbPalette) +
ylab("Gravimetric moisture content") +
xlab("Temperature") +
my_gg_theme
dev.off()



ggplot(diel_phys_dat, aes(log10(cell.counts), log10(resp.per.dry.weight+1), color = Plot)) +
geom_point() +
stat_smooth(method = "lm", se = F, color = "black") +
scale_color_manual(values = cbPalette) +
xlab("Bacterial cell count") +
ylab("resp") +
my_gg_theme

