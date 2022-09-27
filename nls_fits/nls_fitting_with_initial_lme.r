require(tidyverse)

diel_phys_dat = read.table("~/martiny_diel_phys_dat/moisture.tep.res.rin.txt", header = T, sep = "\t")

diel_phys_dat$fixtime = strptime(as.character(diel_phys_dat$Time.bead.beat), '%m/%d/%Y %H:%M')
diel_phys_dat$temptime = strptime(as.character(diel_phys_dat$Temp.time), '%m/%d/%Y %H:%M')

#Get quartiles for T and moisture

temp_moisture_quants =
    diel_phys_dat %>%
        group_by(Plot) %>%
        summarize(
            Moisture_25 = quantile(Moisture, probs = 0.25, na.rm = T),
            Moisture_75 = quantile(Moisture, probs = 0.75, na.rm = T),
            Temperature_25 = quantile(Temperature, probs = 0.25, na.rm = T),
            Temperature_75 = quantile(Temperature, probs = 0.75, na.rm = T)
        ) %>%
    pivot_longer(-Plot, names_to = c("var", "prob"), names_sep = "_")

diel_phys_dat$moisture_quant = vector(mode = "character", length = nrow(diel_phys_dat))
diel_phys_dat$temp_quant = vector(mode = "character", length = nrow(diel_phys_dat))


for(i in 1:nrow(diel_phys_dat)){
    plot.temp = diel_phys_dat$Plot[i]
    mois_25 = ( filter(temp_moisture_quants, Plot == plot.temp & var == "Moisture" & prob == "25") )$value
    mois_75 = ( filter(temp_moisture_quants, Plot == plot.temp & var == "Moisture" & prob == "75") )$value
    temp_25 = ( filter(temp_moisture_quants, Plot == plot.temp & var == "Temperature" & prob == "25") )$value
    temp_75 = ( filter(temp_moisture_quants, Plot == plot.temp & var == "Temperature" & prob == "75") )$value
    
    if(is.na(diel_phys_dat$Moisture[i]) == T){
        diel_phys_dat$moisture_quant[i] = NA
    }else{
    
        if(diel_phys_dat$Moisture[i] < mois_25){
            diel_phys_dat$moisture_quant[i] = "lower"
        }
        else if(diel_phys_dat$Moisture[i] > mois_75){
            diel_phys_dat$moisture_quant[i] = "upper"
        }
        else{
            diel_phys_dat$moisture_quant[i] = "inner"
        }
    }
    
    if(is.na(diel_phys_dat$Temperature[i]) == T){
        diel_phys_dat$temp_quant[i] = NA
    }else{
        
        if(diel_phys_dat$Temperature[i] < temp_25){
            diel_phys_dat$temp_quant[i] = "lower"
        }
        else if(diel_phys_dat$Temperature[i] > temp_75){
            diel_phys_dat$temp_quant[i] = "upper"
        }
        else{
            diel_phys_dat$temp_quant[i] = "inner"
        }
    }
}

#########
#LME fits are the best fit of multiple models compared in `diel_phys_dat_code.r`
#Routines are repeated here

#LME TEMP

temp_sin_lm = lm(Temperature ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat) # plot is not sig
plot(residuals(temp_sin_lm))
qqnorm(residuals(temp_sin_lm))
summary(temp_sin_lm)

temp_sin_plot_lm = lm(Temperature ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat)
temp_sin_time_lm = lm(Temperature ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24)  + hours.cumulative.temp, data = diel_phys_dat)

summary(temp_sin_plot_lm)
anova(temp_sin_lm,temp_sin_plot_lm)
anova(temp_sin_lm,temp_sin_time_lm)

############
#NLS temp

y = diel_phys_dat$Temperature[-35]
x.time = diel_phys_dat$hours.cumulative.temp[-35]

#coefs from estimated linear model
b0 <- coef(temp_sin_lm)[[1]]
alpha <- coef(temp_sin_lm)[[2]]
beta <- coef(temp_sin_lm)[[3]]

#r is amplitude and "offset" (or phase) is phi
r <- sqrt(alpha^2 + beta^2)
phi <- atan2(beta, alpha)

#omega = 2*pi/x where x results in unit of 1 period (e.g. 24 for 24 hr period)
#or x can be estimated by lomb-scargle OR fft https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
#omega = 2*pi/24

nls.mod = nls(y ~ amp*sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))

nls.coefs = coef(nls.mod)


f = function(x, amp, phase, period, C){amp*sin(2*pi*x/24+phase)+C}

plot(y ~ x.time)
curve(f(x, amp = nls.coefs["amp"], phase = nls.coefs["phase"], C = nls.coefs["C"]), add=TRUE ,lwd=2, )


summary(nls.mod)


##########
#Moisture#
##########

diel_phys_dat$Mositure[35] = NA

moist_sin_lm = lm(log10(Moisture*100) ~ sin(2*pi*(hours.cumulative.temp/24)) + cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
#sin fit moist by plot
moist_sin_plot_lm = lm(log10(Moisture*100) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
#adding time
moist_sin_plot_time_lm = lm(log10(Moisture*100) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude) #sig better model with time included
moist_sin_plot_time_plot_lm = lm(log10(Moisture*100) ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24) + Plot*hours.cumulative.temp, data = diel_phys_dat[-35,], na.action = na.exclude)

anova(moist_sin_lm,moist_sin_plot_lm) #sig
anova(moist_sin_lm,moist_sin_plot_time_lm) #sig
anova(moist_sin_plot_lm,moist_sin_plot_time_lm) #sig
anova(moist_sin_plot_lm,moist_sin_plot_time_plot_lm)
anova(moist_sin_plot_time_lm,moist_sin_plot_time_plot_lm)

summary(moist_sin_lm)
summary(moist_sin_plot_lm)
summary(moist_sin_plot_time_lm)
summary(moist_sin_plot_time_plot_lm)

#y = diel_phys_dat$Moisture[-35]
#x.time = diel_phys_dat$hours.cumulative.temp[-35]

#coefs from estimated linear model - using simple sin model
#b0 <- coef(moist_sin_lm)[[1]]
#alpha <- coef(moist_sin_lm)[[2]]
#beta <- coef(moist_sin_lm)[[3]]


#r is amplitude and "offset" (or phase) is phi
#r <- sqrt(alpha^2 + beta^2)
#phi <- atan2(beta, alpha)


moist_sin_plot_time_lm

#omega = 2*pi/x where x results in unit of 1 period (e.g. 24 for 24 hr period)
#or x can be estimated by lomb-scargle OR fft https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
#omega = 2*pi/24

#nls.mod = nls(y ~ amp*sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))


#nls.mod.moist = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r, phase = phi, C = 1), data = diel_phys_dat, na.action = na.omit)

#Fit models separately by plot

#Using best fit lm model for moisture wihtin plots
#
b0.p1 <- coef(moist_sin_plot_lm)[1]
b0.p2 <- coef(moist_sin_plot_lm)[2]
b0.p3 <- coef(moist_sin_plot_lm)[3]
alpha.p1 <- coef(moist_sin_plot_lm)[4]
alpha.p2 <- coef(moist_sin_plot_lm)[5]
alpha.p3 <- coef(moist_sin_plot_lm)[6]
beta.p1 <- coef(moist_sin_plot_lm)[7]
beta.p2 <- coef(moist_sin_plot_lm)[8]
beta.p3 <- coef(moist_sin_plot_lm)[9]


r.p1 <- sqrt(alpha.p1^2 + beta.p1^2)
phi.p1 <- atan2(beta.p1, alpha.p1)
r.p2 <- sqrt(alpha.p2^2 + beta.p2^2)
phi.p2 <- atan2(beta.p2, alpha.p2)
r.p3 <- sqrt(alpha.p3^2 + beta.p3^2)
phi.p3 <- atan2(beta.p3, alpha.p3)

#With linear effect of time
nls.mod.moist.P1.time = nls(Moisture ~ amp*sin(2*pi*(hours.cumulative.temp)/24+phase)+C *hours.cumulative.temp, start = list(amp = r.p1, phase = phi.p1, C = 1), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
#including time(*C) results in non-sig fit of amplitude, changes C and amp so near 0
nls.mod.moist.P2.time = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C*hours.cumulative.temp, start = list(amp = r.p2, phase = phi.p2, C = 1), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.moist.P3.time = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C*hours.cumulative.temp, start = list(amp = r.p3, phase = phi.p3, C = 1), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)

#timeDiff = -24 * (phase/2*pi)
#Above looks incorrect time diff =
#period*phase/(2*pi)

summary(nls.mod.moist.P1.time)
summary(nls.mod.moist.P2.time)
summary(nls.mod.moist.P3.time)

plot(residuals(nls.mod.moist.P2))
plot(residuals(nls.mod.moist.P2.time))
plot(nls.mod.moist.P2)
plot(nls.mod.moist.P2)
#Without linear effect of time
nls.mod.moist.P1 = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C , start = list(amp = r.p1, phase = phi.p1, C = 1), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
#including time(*C) results in non-sig fit of amplitude, changes C and amp so near 0
nls.mod.moist.P2 = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p2, phase = phi.p2, C = 1), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.moist.P3 = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p3, phase = phi.p3, C = 1), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)


nls.mod.moist.P1.log = nls(log10(Moisture*100) ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C , start = list(amp = r.p1, phase = phi.p1, C = 1), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
#including time(*C) results in non-sig fit of amplitude, changes C and amp so near 0
nls.mod.moist.P2.log = nls(log10(Moisture*100) ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p2, phase = phi.p2, C = 1), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.moist.P3.log = nls(log10(Moisture*100) ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p3, phase = phi.p3, C = 1), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)


plot(residuals(nls.mod.moist.P2))
plot(residuals(nls.mod.moist.P2.time))
plot(nls.mod.moist.P1)
plot(nls.mod.moist.P1.time)
plot(nls.mod.moist.P2)
plot(nls.mod.moist.P2.time)
plot(nls.mod.moist.P3)
plot(nls.mod.moist.P3.time)

summary(nls.mod.moist.P1)
summary(nls.mod.moist.P2)
summary(nls.mod.moist.P3)

plot(subset(diel_phys_dat, Plot == "P1")$temptime, predict(nls.mod.moist.P1))

plot(Moisture ~ hours.cumulative.temp, data = subset(diel_phys_dat, Plot == "P1"))
lines(subset(diel_phys_dat, Plot == "P1")$hours.cumulative.temp, predict(nls.mod.moist.P1))

saveRDS(predict(nls.mod.moist.P1), file = "moisture.nls_fit.P1.rds")
saveRDS(predict(nls.mod.moist.P2), file = "moisture.nls_fit.P2.rds")
saveRDS(predict(nls.mod.moist.P3), file = "moisture.nls_fit.P3.rds")

#############
#Temperature#
#############

#coefs from estimated linear model - using simple sin model
b0 <- coef(temp_sin_lm)[[1]]
alpha <- coef(temp_sin_lm)[[2]]
beta <- coef(temp_sin_lm)[[3]]


#r is amplitude and "offset" (or phase) is phi
r <- sqrt(alpha^2 + beta^2)
phi <- atan2(beta, alpha)

#omega = 2*pi/x where x results in unit of 1 period (e.g. 24 for 24 hr period)
#or x can be estimated by lomb-scargle OR fft https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
#omega = 2*pi/24


#Phase help constant
nls.mod.temp = nls(Temperature ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r, C = 1, phase = phi), data = diel_phys_dat, na.action = na.omit)

#Fit models separately by plot
#Plot subset is not sig in lm fit so just use parameters from basic model

b0.p1 <- coef(temp_sin_plot_lm)[1]
b0.p2 <- coef(temp_sin_plot_lm)[2]
b0.p3 <- coef(temp_sin_plot_lm)[3]
alpha.p1 <- coef(temp_sin_plot_lm)[4]
alpha.p2 <- coef(temp_sin_plot_lm)[5]
alpha.p3 <- coef(temp_sin_plot_lm)[6]
beta.p1 <- coef(temp_sin_plot_lm)[7]
beta.p2 <- coef(temp_sin_plot_lm)[8]
beta.p3 <- coef(temp_sin_plot_lm)[9]


r.p1 <- sqrt(alpha.p1^2 + beta.p1^2)
phi.p1 <- atan2(beta.p1, alpha.p1)
r.p2 <- sqrt(alpha.p2^2 + beta.p2^2)
phi.p2 <- atan2(beta.p2, alpha.p2)
r.p3 <- sqrt(alpha.p3^2 + beta.p3^2)
phi.p3 <- atan2(beta.p3, alpha.p3)


nls.mod.temp.P1 = nls(Temperature ~ amp*sin(2*pi*hours.cumulative.temp/24+phi)+C, start = list(amp = r.p1, C = 8), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
nls.mod.temp.P2 = nls(Temperature ~ amp*sin(2*pi*hours.cumulative.temp/24+phi)+C, start = list(amp = r.p2, C = 8), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.temp.P3 = nls(Temperature ~ amp*sin(2*pi*hours.cumulative.temp/24+phi)+C, start = list(amp = r.p3, C = 8), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)

predict(nls.mod.temp.P1) #y coordinates for line function


########
#Respiration by litter dry weight
########
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

#coefs from estimated linear model - using simple sin model
b0 <- coef(resp_sin_lm)[[1]]
alpha <- coef(resp_sin_lm)[[2]]
beta <- coef(resp_sin_lm)[[3]]


#r is amplitude and "offset" (or phase) is phi
r <- sqrt(alpha^2 + beta^2)
phi <- atan2(beta, alpha)

#omega = 2*pi/x where x results in unit of 1 period (e.g. 24 for 24 hr period)
#or x can be estimated by lomb-scargle OR fft https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
#omega = 2*pi/24

C.start = sum(diel_phys_dat$resp.per.dry.weight, na.rm = T)/84


nls.mod.resp = nls(resp.per.dry.weight ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r, phase = phi, C = C.start), data = diel_phys_dat, na.action = na.omit)

#Fit models separately by plot

b0.p1 <- coef(resp_sin_plot_lm)[1]
b0.p2 <- coef(resp_sin_plot_lm)[2]
b0.p3 <- coef(resp_sin_plot_lm)[3]
alpha.p1 <- coef(resp_sin_plot_lm)[4]
alpha.p2 <- coef(resp_sin_plot_lm)[5]
alpha.p3 <- coef(resp_sin_plot_lm)[6]
beta.p1 <- coef(resp_sin_plot_lm)[7]
beta.p2 <- coef(resp_sin_plot_lm)[8]
beta.p3 <- coef(resp_sin_plot_lm)[9]


r.p1 <- sqrt(alpha.p1^2 + beta.p1^2)
phi.p1 <- atan2(beta.p1, alpha.p1)
r.p2 <- sqrt(alpha.p2^2 + beta.p2^2)
phi.p2 <- atan2(beta.p2, alpha.p2)
r.p3 <- sqrt(alpha.p3^2 + beta.p3^2)
phi.p3 <- atan2(beta.p3, alpha.p3)



nls.mod.resp.P1 = nls(resp.per.dry.weight ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p1, phase = phi.p1, C = C.start), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
nls.mod.resp.P2 = nls(resp.per.dry.weight ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p2, phase = phi.p2, C = C.start), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.resp.P3 = nls(resp.per.dry.weight ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p3, phase = phi.p3, C = C.start), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)

predict(nls.mod.resp.P3) #y coordinates for line function
plot(resp.per.dry.weight ~ hours.cumulative.temp, data = subset(diel_phys_dat, Plot == "P3"))
lines(na.omit(subset(diel_phys_dat, Plot == "P3")$hours.cumulative.temp), predict(nls.mod.resp.P3))

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
anova(cells_sin_lm)
summary(cells_sin_lm)

anova(cells_sin_plot_lm)
summary(cells_sin_plot_lm)

#Fit models separately by plot

b0.p1 <- coef(cells_sin_plot_lm)[1]
b0.p2 <- coef(cells_sin_plot_lm)[2]
b0.p3 <- coef(cells_sin_plot_lm)[3]
alpha.p1 <- coef(cells_sin_plot_lm)[4]
alpha.p2 <- coef(cells_sin_plot_lm)[5]
alpha.p3 <- coef(cells_sin_plot_lm)[6]
beta.p1 <- coef(cells_sin_plot_lm)[7]
beta.p2 <- coef(cells_sin_plot_lm)[8]
beta.p3 <- coef(cells_sin_plot_lm)[9]


r.p1 <- sqrt(alpha.p1^2 + beta.p1^2)
phi.p1 <- atan2(beta.p1, alpha.p1)
r.p2 <- sqrt(alpha.p2^2 + beta.p2^2)
phi.p2 <- atan2(beta.p2, alpha.p2)
r.p3 <- sqrt(alpha.p3^2 + beta.p3^2)
phi.p3 <- atan2(beta.p3, alpha.p3)



nls.mod.cells.P1 = nls(log10(cell.counts) ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p1, phase = phi.p1, C = C.start), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
nls.mod.cells.P2 = nls(log10(cell.counts) ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p2, phase = phi.p2, C = C.start), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.cells.P3 = nls(log10(cell.counts) ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p3, phase = phi.p3, C = C.start), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)

summary(nls.mod.cells.P1)
summary(nls.mod.cells.P2)
summary(nls.mod.cells.P3)
####################################
#Plot for pub

temperature_plot_dat = data.frame(
Plot = paste0("P", c(
subset(diel_phys_dat, Plot == "P1")$Plot,
subset(diel_phys_dat, Plot == "P2")$Plot,
subset(diel_phys_dat, Plot == "P3")$Plot
) ),
Temperature = c(
subset(diel_phys_dat, Plot == "P1")$Temperature,
subset(diel_phys_dat, Plot == "P2")$Temperature,
subset(diel_phys_dat, Plot == "P3")$Temperature
),
temptime = c(
subset(diel_phys_dat, Plot == "P1")$temptime,
subset(diel_phys_dat, Plot == "P2")$temptime,
subset(diel_phys_dat, Plot == "P3")$temptime
),
NLS_fits = c(
predict(nls.mod.temp.P1),
predict(nls.mod.temp.P2),
predict(nls.mod.temp.P3)
),
quants = c(
subset(diel_phys_dat, Plot == "P1")$temp_quant,
subset(diel_phys_dat, Plot == "P2")$temp_quant,
subset(diel_phys_dat, Plot == "P3")$temp_quant
)
)

#with quants and faceted
sample_exclude = c("P1 T0", "P1 T1", "P1 T2", "P1 T3", "P1 T4", "P1 T5",
"P2 T0", "P2 T1", "P2 T2", "P2 T3", "P2 T4", "P2 T5",
"P3 T0", "P3 T1", "P3 T2", "P3 T3", "P3 T4", "P3 T5"
)


p1.quant = ggplot(temperature_plot_dat, aes(x = as.POSIXct(temptime), y = Temperature, fill = quants)) +
facet_wrap(~Plot) +
geom_point(size = 3, shape = 21) +
geom_line(aes(y = NLS_fits, group = Plot))  +
scale_fill_manual(
labels = c("lower" = "lower quantile", "inner"="inner quantile", "upper"="upper quantile"),
values = c("lower" = "black", "inner"="dark grey", "upper"="white"),
breaks = c("lower", "inner", "upper"),
) +
xlab("Date time") +
ylab(expression(paste("Temperature ("~degree, "C)"))) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
labs(title = "A") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.08),
axis.title.x = element_blank(),
axis.text.x = element_blank())

p1 = ggplot(temperature_plot_dat, aes(x = as.POSIXct(temptime), y = Temperature, fill = Plot, shape = Plot, linetype = Plot)) +
geom_point(size = 3) +
geom_line(aes(y = NLS_fits))  +
#scale_y_continuous(labels = fancy_scientific, trans = "log10")+
scale_fill_manual(
labels = c("P1" = "Plot 1", "P2"="Plot 2", "P3"="Plot 3"),
values = c("P1" = "black", "P2"="white", "P3"="dark grey"),
breaks = c("P1", "P2", "P3"),
) +
scale_linetype_manual(
labels = c("P1" = "Plot 1", "P2"="Plot 2", "P3"="Plot 3"),
values = c("P1" = 1, "P2"=2, "P3"=3),
breaks = c("P1", "P2", "P3")
) +
scale_shape_manual(
labels = c("P1" = "Plot 1", "P2"="Plot 2", "P3"="Plot 3"),
values = c("P1" = 21, "P2"=22, "P3"=24),
breaks = c("P1", "P2", "P3")
)  +
xlab("Date time") +
ylab(expression(paste("Temperature ("~degree, "C)"))) +
#scale_x_datetime( date_breaks = "6 hours", date_labels = "%m/%d\n%H:%M") +
scale_x_datetime( date_breaks = "2 hours", date_labels = "%H:%M") +
labs(title = "A") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.08),
axis.title.x = element_blank(),
axis.text.x = element_text(angle = 90, size = 5)
)


moisture_plot_dat = data.frame(
Plot = paste0("P", c(
subset(diel_phys_dat, Plot == "P1")$Plot,
subset(diel_phys_dat, Plot == "P2" & Moisture >= 0)$Plot,
subset(diel_phys_dat, Plot == "P3")$Plot
) ),
Moisture = c(
subset(diel_phys_dat, Plot == "P1")$Moisture,
subset(diel_phys_dat, Plot == "P2" & Moisture >= 0)$Moisture,
subset(diel_phys_dat, Plot == "P3")$Moisture
),
temptime = c(
subset(diel_phys_dat, Plot == "P1")$temptime,
subset(diel_phys_dat, Plot == "P2" & Moisture >= 0)$temptime,
subset(diel_phys_dat, Plot == "P3")$temptime
),
NLS_fits = c(
predict(nls.mod.moist.P1),
predict(nls.mod.moist.P2),
predict(nls.mod.moist.P3.time)
),
quants = c(
subset(diel_phys_dat, Plot == "P1")$moisture_quant,
subset(diel_phys_dat, Plot == "P2" & Moisture >= 0)$moisture_quant,
subset(diel_phys_dat, Plot == "P3")$moisture_quant
)
)


moisture_plot_dat[72:107,"NLS_fits"] = NA

p2.quant = ggplot(moisture_plot_dat, aes(x = as.POSIXct(temptime), y = Moisture, fill = quants)) +
facet_wrap(~Plot) +
geom_point(size = 3, shape = 21) +
geom_line(aes(y = NLS_fits, group = Plot))  +
scale_fill_manual(
labels = c("lower" = "lower quantile", "inner"="inner quantile", "upper"="upper quantile"),
values = c("lower" = "black", "inner"="dark grey", "upper"="white"),
breaks = c("lower", "inner", "upper"),
) +
xlab("Date time") +
ylab("Moisture content (%)") +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
labs(title = "B") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.08),
axis.title.x = element_blank(),
axis.text.x = element_blank())

pdf("moisture_and_temperature_w_quantiles.pdf", width = 8, height = 8)
grid.arrange(p1.quant, p2.quant, nrow = 2)
dev.off()

p2 = ggplot(moisture_plot_dat, aes(x = as.POSIXct(temptime), y = Moisture, fill = Plot, shape = Plot, linetype = Plot)) +
geom_point(size = 3) +
geom_line(aes(y = NLS_fits))  +
#scale_y_continuous(labels = fancy_scientific, trans = "log10")+
scale_fill_manual(
labels = c("P1" = "Plot 1", "P2"="Plot 2", "P3"="Plot 3"),
values = c("P1" = "black", "P2"="white", "P3"="dark grey"),
breaks = c("P1", "P2", "P3"),
) +
scale_linetype_manual(
labels = c("P1" = "Plot 1", "P2"="Plot 2", "P3"="Plot 3"),
values = c("P1" = 1, "P2"=2, "P3"=3),
breaks = c("P1", "P2", "P3")
) +
scale_shape_manual(
labels = c("P1" = "Plot 1", "P2"="Plot 2", "P3"="Plot 3"),
values = c("P1" = 21, "P2"=22, "P3"=24),
breaks = c("P1", "P2", "P3")
)  +
xlab("Date time") +
ylab("Moisture content (%)") +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
labs(title = "B") +
my_gg_theme #+
#theme(plot.title = element_text(hjust = -0.08),
#axis.title.x = element_blank(),
#axis.text.x = element_blank())

pdf("martiny_diel_phys_dat/resp_per_lit_nls_date_time.pdf", width = 8, height = 4)
ggplot() +
geom_point(data = subset(diel_phys_dat, Plot == "P1"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P1"), aes(y = c(predict(nls.mod.resp.P1), rep(NA,8)), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P2"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[3]), size = 2) +
#geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P2")), aes(y = predict(nls.mod.resp.P2), x = as.POSIXct(temptime)), color = cbPalette[3], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P3"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P3"), aes(y = c(predict(nls.mod.resp.P3), rep(NA,8)), x = as.POSIXct(temptime)), color = cbPalette[4], size = 1.5) +
xlab("Date time") +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " g"^-1, "litter" ))) +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "#56B4E9"="Plot 2", "#009E73"="Plot 3"),
    values = c("#E69F00" = cbPalette[2], "#56B4E9"=cbPalette[3], "#009E73"=cbPalette[4]),
    breaks = c("#E69F00", "#56B4E9", "#009E73")) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
my_gg_theme
dev.off()

p3 = ggplot() +
geom_point(data = subset(diel_phys_dat, Plot == "P1"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P1"), aes(y = c(predict(nls.mod.resp.P1), rep(NA,8)), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P2"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = "black"), size = 2) +
#geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P2")), aes(y = predict(nls.mod.resp.P2), x = as.POSIXct(temptime)), color = cbPalette[3], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P3"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P3"), aes(y = c(predict(nls.mod.resp.P3), rep(NA,8)), x = as.POSIXct(temptime)), color = cbPalette[4], size = 1.5) +
xlab("Date time") +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " g"^-1, "litter" ))) +
labs(title = "C") +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "black"="Plot 2", "#009E73"="Plot 3"),
values = c("#E69F00" = cbPalette[2], "black"="black", "#009E73"=cbPalette[4]),
breaks = c("#E69F00", "black", "#009E73")) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.0775, margin=margin(b = 30)),
legend.position = c(0.9,0.5))

###################
#Plot for pub
require(gridExtra)

p1 = ggplot() +
geom_point(data = subset(diel_phys_dat, Plot == "P1"), aes(y = Temperature, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P1"), aes(y = predict(nls.mod.temp.P1), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P2"), aes(y = Temperature, x = as.POSIXct(temptime), color = "black"), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P2"), aes(y = predict(nls.mod.temp.P2), x = as.POSIXct(temptime)), color = "black", size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P3"), aes(y = Temperature, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P3"), aes(y = predict(nls.mod.temp.P3), x = as.POSIXct(temptime)), color = cbPalette[4], size = 1.5) +
xlab("Date time") +
ylab(expression("Temperature " ( degree*C))) +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "black"="Plot 2", "#009E73"="Plot 3"),
values = c("#E69F00" = cbPalette[2], "black"="black", "#009E73"=cbPalette[4]),
breaks = c("#E69F00", "black", "#009E73"),
guide = F) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
#labs(title = "A") +
my_gg_theme +
theme(
plot.title = element_text(hjust = -0.08),
axis.title.x = element_blank(),
axis.text.x = element_blank()
)

p2 = ggplot() +
geom_point(data = subset(diel_phys_dat, Plot == "P1"), aes(y = Moisture, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P1"), aes(y = predict(nls.mod.moist.P1.time), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = subset(diel_phys_dat[-35,], Plot == "P2"), aes(y = Moisture, x = as.POSIXct(temptime), color = "black"), size = 2) +
geom_line(data = subset(diel_phys_dat[-35,], Plot == "P2"), aes(y = predict(nls.mod.moist.P2.time), x = as.POSIXct(temptime)), color = "black", size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P3"), aes(y = Moisture, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
#geom_line(data = subset(diel_phys_dat, Plot == "P3"), aes(y = predict(nls.mod.moist.P3.time), x = as.POSIXct(datetime)), color = cbPalette[4], size = 1.5) +
xlab("Date time") +
ylab("Moisture content (%)") +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "black"="Plot 2", "#009E73"="Plot 3"),
values = c("#E69F00" = cbPalette[2], "black"="black", "#009E73"=cbPalette[4]),
breaks = c("#E69F00", "black", "#009E73"),
guide = F) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
#labs(title = "B") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.08),
axis.title.x = element_blank(),
axis.text.x = element_blank())

p3 = ggplot() +
geom_point(data = subset(diel_phys_dat, Plot == "P1"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P1"), aes(y = c(predict(nls.mod.resp.P1), rep(NA,8)), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P2"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = "black"), size = 2) +
#geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P2")), aes(y = predict(nls.mod.resp.P2), x = as.POSIXct(temptime)), color = cbPalette[3], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P3"), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P3"), aes(y = c(predict(nls.mod.resp.P3), rep(NA,8)), x = as.POSIXct(temptime)), color = cbPalette[4], size = 1.5) +
xlab("Date time") +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " g"^-1, "litter" ))) +
#labs(title = "C") +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "black"="Plot 2", "#009E73"="Plot 3"),
values = c("#E69F00" = cbPalette[2], "black"="black", "#009E73"=cbPalette[4]),
breaks = c("#E69F00", "black", "#009E73")) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.0775, margin=margin(b = 30)),
legend.position = c(0.9,0.5))

pdf(file = "~/martiny_diel_phys_dat/Fig_1_temp_moisture_resp_nls.pdf", width = 10, height = 12)
grid.arrange(p1,p2,p3, nrow = 3)
dev.off()

pdf(file = "~/martiny_diel_phys_dat/Fig_1_temp_moisture_resp_nls.smaller.pdf", width = 7, height = 10)
grid.arrange(p1,p2,p3, nrow = 3, heights = c(0.3,0.3,0.4))
dev.off()

########
#Respiration by cell_abd
########
resp_p_cells_sin_plot_lm = lm(resp.per.cell ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)
resp_p_cells_sin_lm = lm(resp.per.cell ~ sin(2*pi*hours.cumulative.temp/24) + cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)

#coefs from estimated linear model - using simple sin model
b0 <- coef(resp_p_cells_sin_lm)[[1]]
alpha <- coef(resp_p_cells_sin_lm)[[2]]
beta <- coef(resp_p_cells_sin_lm)[[3]]


#r is amplitude and "offset" (or phase) is phi
r <- sqrt(alpha^2 + beta^2)
phi <- atan2(beta, alpha)

#omega = 2*pi/x where x results in unit of 1 period (e.g. 24 for 24 hr period)
#or x can be estimated by lomb-scargle OR fft https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
#omega = 2*pi/24

C.start = sum(diel_phys_dat$resp.per.cell, na.rm = T)/84


nls.mod.resp.cell = nls(resp.per.cell ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r, phase = phi, C = C.start), data = diel_phys_dat, na.action = na.omit)

#Fit models separately by plot

b0.p1 <- coef(resp_p_cells_sin_plot_lm)[1]
b0.p2 <- coef(resp_p_cells_sin_plot_lm)[2]
b0.p3 <- coef(resp_p_cells_sin_plot_lm)[3]
alpha.p1 <- coef(resp_p_cells_sin_plot_lm)[4]
alpha.p2 <- coef(resp_p_cells_sin_plot_lm)[5]
alpha.p3 <- coef(resp_p_cells_sin_plot_lm)[6]
beta.p1 <- coef(resp_p_cells_sin_plot_lm)[7]
beta.p2 <- coef(resp_p_cells_sin_plot_lm)[8]
beta.p3 <- coef(resp_p_cells_sin_plot_lm)[9]


r.p1 <- sqrt(alpha.p1^2 + beta.p1^2)
phi.p1 <- atan2(beta.p1, alpha.p1)
r.p2 <- sqrt(alpha.p2^2 + beta.p2^2)
phi.p2 <- atan2(beta.p2, alpha.p2)
r.p3 <- sqrt(alpha.p3^2 + beta.p3^2)
phi.p3 <- atan2(beta.p3, alpha.p3)



nls.mod.resp.cell.P1 = nls(resp.per.cell ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p1, phase = phi.p1, C = C.start), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
nls.mod.resp.cell.P2 = nls(resp.per.cell ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p2, phase = phi.p2, C = C.start), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.resp.cell.P3 = nls(resp.per.cell ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p3, phase = phi.p3, C = C.start), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)

predict(nls.mod.temp.P1) #y coordinates for line function

summary(nls.mod.resp.cell.P1)
summary(nls.mod.resp.cell.P2)
summary(nls.mod.resp.cell.P3)

pdf("martiny_diel_phys_dat/resp_per_cell_nls.pdf", width = 5, height = 4)
ggplot() +
geom_point(data = na.omit(subset(diel_phys_dat, Plot == "P1")), aes(y = resp.per.cell, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P1")), aes(y = predict(nls.mod.resp.cell.P1), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = na.omit(subset(diel_phys_dat, Plot == "P2")), aes(y = resp.per.cell, x = as.POSIXct(temptime), color = cbPalette[3]), size = 2) +
geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P2")), aes(y = predict(nls.mod.resp.cell.P2), x = as.POSIXct(temptime)), color = cbPalette[3], size = 1.5) +
geom_point(data = na.omit(subset(diel_phys_dat, Plot == "P3")), aes(y = resp.per.cell, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P3")), aes(y = predict(nls.mod.resp.cell.P3), x = as.POSIXct(temptime)), color = cbPalette[4], size = 1.5) +
xlab("Time (hrs)") +
ylab(expression(paste(mu, "g CO"[2],"-C min"^-1, " bacterial cell"^-1))) +
scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "#56B4E9"="Plot 2", "#009E73"="Plot 3"),
    values = c("#E69F00" = cbPalette[2], "#56B4E9"=cbPalette[3], "#009E73"=cbPalette[4]),
    breaks = c("#E69F00", "#56B4E9", "#009E73")) +
my_gg_theme
dev.off()





##########################
#Plotting of temp and moisture with quantiles
#this routine excludes the first six timepoints because the quantile data will be used for DE expression analysis

#with quants and faceted
sample_exclude = c("P1 T0", "P1 T1", "P1 T2", "P1 T3", "P1 T4", "P1 T5",
"P2 T0", "P2 T1", "P2 T2", "P2 T3", "P2 T4", "P2 T5",
"P3 T0", "P3 T1", "P3 T2", "P3 T3", "P3 T4", "P3 T5"
)
#sample_exclude = c(1,2,3,4,5)


moisture_plot_dat.quant = data.frame(
Plot = paste0("P", c(
subset(diel_phys_dat, Plot == "P1")$Plot,
subset(diel_phys_dat, Plot == "P2" & Moisture >= 0)$Plot,
subset(diel_phys_dat, Plot == "P3")$Plot
) ),
Moisture = c(
subset(diel_phys_dat, Plot == "P1")$Moisture,
subset(diel_phys_dat, Plot == "P2" & Moisture >= 0)$Moisture,
subset(diel_phys_dat, Plot == "P3")$Moisture
),
temptime = c(
subset(diel_phys_dat, Plot == "P1")$temptime,
subset(diel_phys_dat, Plot == "P2" & Moisture >= 0)$temptime,
subset(diel_phys_dat, Plot == "P3")$temptime
),
NLS_fits = c(
predict(nls.mod.moist.P1),
predict(nls.mod.moist.P2),
predict(nls.mod.moist.P3.time)
),
Sample.ID = rbind(
subset(diel_phys_dat, Plot == "P1")[,"Sample.ID", drop = F],
subset(diel_phys_dat, Plot == "P2" & Moisture >= 0)[,"Sample.ID", drop = F],
subset(diel_phys_dat, Plot == "P3")[,"Sample.ID", drop = F]
)
)

moisture_plot_dat.quant = moisture_plot_dat.quant %>%
filter(!Sample.ID %in% sample_exclude)

moisture_quants =
moisture_plot_dat.quant %>%
group_by(Plot) %>%
summarize(
Moisture_25 = quantile(Moisture, probs = 0.25, na.rm = T),
Moisture_75 = quantile(Moisture, probs = 0.75, na.rm = T)
) %>%
pivot_longer(-Plot, names_to = c("var", "prob"), names_sep = "_")

moisture_plot_dat.quant$moisture_quant = vector(mode = "character", length = nrow(moisture_plot_dat.quant))

for(i in 1:nrow(moisture_plot_dat.quant)){
    plot.temp = moisture_plot_dat.quant$Plot[i]
    moisture_25 = ( filter(moisture_quants, Plot == plot.temp & prob == "25") )$value
    moisture_75 = ( filter(moisture_quants, Plot == plot.temp & prob == "75") )$value
    
    if(is.na(moisture_plot_dat.quant$Moisture[i]) == T){
        moisture_plot_dat.quant$temp_quant[i] = NA
    }else{
        
        if(moisture_plot_dat.quant$Moisture[i] < moisture_25){
            moisture_plot_dat.quant$moisture_quant[i] = "lower"
        }
        else if(moisture_plot_dat.quant$Moisture[i] > moisture_75){
            moisture_plot_dat.quant$moisture_quant[i] = "upper"
        }
        else{
            moisture_plot_dat.quant$moisture_quant[i] = "inner"
        }
    }
}

p2.quant = ggplot(moisture_plot_dat.quant, aes(x = as.POSIXct(temptime), y = Moisture, fill = moisture_quant)) +
facet_wrap(~Plot) +
geom_point(size = 3, shape = 21) +
geom_line(aes(y = NLS_fits, group = Plot))  +
scale_fill_manual(
labels = c("lower" = "lower quantile", "inner"="inner quantile", "upper"="upper quantile"),
values = c("lower" = "black", "inner"="dark grey", "upper"="white"),
breaks = c("lower", "inner", "upper"),
) +
xlab("Date time") +
ylab("Moisture content (%)") +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
labs(title = "B") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.08),
axis.title.x = element_blank(),
axis.text.x = element_blank())

#Temp
temperature_plot_dat.quants = data.frame(
Plot = paste0("P", c(
subset(diel_phys_dat, Plot == "P1")$Plot,
subset(diel_phys_dat, Plot == "P2")$Plot,
subset(diel_phys_dat, Plot == "P3")$Plot
) ),
Temperature = c(
subset(diel_phys_dat, Plot == "P1")$Temperature,
subset(diel_phys_dat, Plot == "P2")$Temperature,
subset(diel_phys_dat, Plot == "P3")$Temperature
),
temptime = c(
subset(diel_phys_dat, Plot == "P1")$temptime,
subset(diel_phys_dat, Plot == "P2")$temptime,
subset(diel_phys_dat, Plot == "P3")$temptime
),
NLS_fits = c(
predict(nls.mod.temp.P1),
predict(nls.mod.temp.P2),
predict(nls.mod.temp.P3)
),
Sample.ID = rbind(
subset(diel_phys_dat, Plot == "P1")[,"Sample.ID", drop = F],
subset(diel_phys_dat, Plot == "P2")[,"Sample.ID", drop = F],
subset(diel_phys_dat, Plot == "P3")[,"Sample.ID", drop = F]
)
)

temperature_plot_dat.quants = temperature_plot_dat.quants %>%
    filter(!Sample.ID %in% sample_exclude)


temp_quants =
temperature_plot_dat.quants %>%
group_by(Plot) %>%
summarize(
Temperature_25 = quantile(Temperature, probs = 0.25, na.rm = T),
Temperature_75 = quantile(Temperature, probs = 0.75, na.rm = T)
) %>%
pivot_longer(-Plot, names_to = c("var", "prob"), names_sep = "_")

temperature_plot_dat.quants$temp_quant = vector(mode = "character", length = nrow(temperature_plot_dat.quants))


for(i in 1:nrow(temperature_plot_dat.quants)){
    plot.temp = temperature_plot_dat.quants$Plot[i]
    temp_25 = ( filter(temp_quants, Plot == plot.temp & prob == "25") )$value
    temp_75 = ( filter(temp_quants, Plot == plot.temp & prob == "75") )$value
    
    if(is.na(temperature_plot_dat.quants$Temperature[i]) == T){
        temperature_plot_dat.quants$temp_quant[i] = NA
    }else{
        
        if(temperature_plot_dat.quants$Temperature[i] < temp_25){
            temperature_plot_dat.quants$temp_quant[i] = "lower"
        }
        else if(temperature_plot_dat.quants$Temperature[i] > temp_75){
            temperature_plot_dat.quants$temp_quant[i] = "upper"
        }
        else{
            temperature_plot_dat.quants$temp_quant[i] = "inner"
        }
    }
}


p1.quant = ggplot(temperature_plot_dat.quants, aes(x = as.POSIXct(temptime), y = Temperature, fill = temp_quant)) +
facet_wrap(~Plot) +
geom_point(size = 3, shape = 21) +
geom_line(aes(y = NLS_fits, group = Plot))  +
scale_fill_manual(
labels = c("lower" = "lower quantile", "inner"="inner quantile", "upper"="upper quantile"),
values = c("lower" = "black", "inner"="dark grey", "upper"="white"),
breaks = c("lower", "inner", "upper"),
) +
xlab("Date time") +
ylab(expression(paste("Temperature ("~degree, "C)"))) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
labs(title = "A") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.08),
axis.title.x = element_blank(),
axis.text.x = element_blank())



pdf("moisture_and_temperature_w_quantiles.pdf", width = 8, height = 8)
grid.arrange(p1.quant, p2.quant, nrow = 2)
dev.off()
