

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


#y = diel_phys_dat$Moisture[-35]
#x.time = diel_phys_dat$hours.cumulative.temp[-35]

#coefs from estimated linear model - using simple sin model
b0 <- coef(moist_sin_lm)[[1]]
alpha <- coef(moist_sin_lm)[[2]]
beta <- coef(moist_sin_lm)[[3]]


#r is amplitude and "offset" (or phase) is phi
r <- sqrt(alpha^2 + beta^2)
phi <- atan2(beta, alpha)


moist_sin_plot_time_lm

#omega = 2*pi/x where x results in unit of 1 period (e.g. 24 for 24 hr period)
#or x can be estimated by lomb-scargle OR fft https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
#omega = 2*pi/24

#nls.mod = nls(y ~ amp*sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))


nls.mod.moist = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r, phase = phi, C = 1), data = diel_phys_dat, na.action = na.omit)

#Fit models separately by plot

#Using best fit lm model for moisture wihtin plots #4 is skipped because this is the moisture parameter estimate
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


nls.mod.moist.P1.time = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C *hours.cumulative.temp, start = list(amp = r.p1, phase = phi.p1, C = 1), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
#including time(*C) results in non-sig fit of amplitude, changes C and amp so near 0
nls.mod.moist.P2.time = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C*hours.cumulative.temp, start = list(amp = r.p2, phase = phi.p2, C = 1), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.moist.P3.time = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C*hours.cumulative.temp, start = list(amp = r.p3, phase = phi.p3, C = 1), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)


nls.mod.moist.P1 = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C , start = list(amp = r.p1, phase = phi.p1, C = 1), data = subset(diel_phys_dat, Plot == "P1"), na.action = na.omit)
#including time(*C) results in non-sig fit of amplitude, changes C and amp so near 0
nls.mod.moist.P2 = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p2, phase = phi.p2, C = 1), data = subset(diel_phys_dat, Plot == "P2"), na.action = na.omit)
nls.mod.moist.P3 = nls(Moisture ~ amp*sin(2*pi*hours.cumulative.temp/24+phase)+C, start = list(amp = r.p3, phase = phi.p3, C = 1), data = subset(diel_phys_dat, Plot == "P3"), na.action = na.omit)

plot(Moisture ~ hours.cumulative.temp, data = subset(diel_phys_dat, Plot == "P1"))
lines(subset(diel_phys_dat, Plot == "P1")$hours.cumulative.temp, predict(nls.mod.moist.P1))

pdf("martiny_diel_phys_dat/moisture_nls.pdf", width = 5, height = 4)
ggplot() +
geom_point(data = subset(diel_phys_dat, Plot == "P1"), aes(y = Moisture*100, x = hours.cumulative.temp, color = cbPalette[1]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P1"), aes(y = predict(nls.mod.moist.P1)*100, x = hours.cumulative.temp), color = cbPalette[1], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P2"), aes(y = Moisture*100, x = hours.cumulative.temp, color = cbPalette[2]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P2"), aes(y = predict(nls.mod.moist.P2)*100, x = hours.cumulative.temp), color = cbPalette[2], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P3"), aes(y = Moisture*100, x = hours.cumulative.temp, color = cbPalette[3]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P3"), aes(y = predict(nls.mod.moist.P3)*100, x = hours.cumulative.temp), color = cbPalette[3], size = 1.5) +
xlab("Time (hrs)") +
ylab("Moisture content (%)") +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "#56B4E9"="Plot 2", "#009E73"="Plot 3"),
values = c("#E69F00" = cbPalette[1], "#56B4E9"=cbPalette[2], "#009E73"=cbPalette[3]),
breaks = c("#E69F00", "#56B4E9", "#009E73")) +
my_gg_theme
dev.off()


pdf("martiny_diel_phys_dat/moisture_time_nls_date_time.pdf", width = 10, height = 4)
ggplot() +
geom_point(data = subset(diel_phys_dat, Plot == "P1"), aes(y = Moisture, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P1"), aes(y = predict(nls.mod.moist.P1.time), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P2"), aes(y = Moisture, x = as.POSIXct(temptime), color = cbPalette[3]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P2"), aes(y = predict(nls.mod.moist.P2.time), x = as.POSIXct(temptime)), color = cbPalette[3], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P3"), aes(y = Moisture, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
#geom_line(data = subset(diel_phys_dat, Plot == "P3"), aes(y = predict(nls.mod.moist.P3.time), x = as.POSIXct(datetime)), color = cbPalette[4], size = 1.5) +
xlab("Date time") +
ylab("Moisture content (%)") +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "#56B4E9"="Plot 2", "#009E73"="Plot 3"),
values = c("#E69F00" = cbPalette[2], "#56B4E9"=cbPalette[3], "#009E73"=cbPalette[4]),
breaks = c("#E69F00", "#56B4E9", "#009E73")) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
my_gg_theme
dev.off()

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
nls.mod.temp = nls(Temperature ~ amp*sin(2*pi*hours.cumulative.temp/24+phi)+C, start = list(amp = r, C = 1), data = diel_phys_dat, na.action = na.omit)

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



pdf("martiny_diel_phys_dat/temperature_nls_date.pdf", width = 10, height = 4)
ggplot() +
geom_point(data = subset(diel_phys_dat, Plot == "P1"), aes(y = Temperature, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P1"), aes(y = predict(nls.mod.temp.P1), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P2"), aes(y = Temperature, x = as.POSIXct(temptime), color = cbPalette[3]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P2"), aes(y = predict(nls.mod.temp.P2), x = as.POSIXct(temptime)), color = cbPalette[3], size = 1.5) +
geom_point(data = subset(diel_phys_dat, Plot == "P3"), aes(y = Temperature, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
geom_line(data = subset(diel_phys_dat, Plot == "P3"), aes(y = predict(nls.mod.temp.P3), x = as.POSIXct(temptime)), color = cbPalette[4], size = 1.5) +
xlab("Date time") +
ylab(expression("Temperature " ( degree*C))) +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "#56B4E9"="Plot 2", "#009E73"="Plot 3"),
values = c("#E69F00" = cbPalette[2], "#56B4E9"=cbPalette[3], "#009E73"=cbPalette[4]),
breaks = c("#E69F00", "#56B4E9", "#009E73")) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
my_gg_theme
dev.off()

########
#Respiration by litter dry weight
########

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


pdf("martiny_diel_phys_dat/resp_per_lit_nls_date_time.pdf", width = 8, height = 4)
ggplot() +
geom_point(data = na.omit(subset(diel_phys_dat, Plot == "P1")), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[2]), size = 2) +
geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P1")), aes(y = predict(nls.mod.resp.P1), x = as.POSIXct(temptime)), color = cbPalette[2], size = 1.5) +
geom_point(data = na.omit(subset(diel_phys_dat, Plot == "P2")), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[3]), size = 2) +
#geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P2")), aes(y = predict(nls.mod.resp.P2), x = as.POSIXct(temptime)), color = cbPalette[3], size = 1.5) +
geom_point(data = na.omit(subset(diel_phys_dat, Plot == "P3")), aes(y = resp.per.dry.weight, x = as.POSIXct(temptime), color = cbPalette[4]), size = 2) +
geom_line(data = na.omit(subset(diel_phys_dat, Plot == "P3")), aes(y = predict(nls.mod.resp.P3), x = as.POSIXct(temptime)), color = cbPalette[4], size = 1.5) +
xlab("Date time") +
ylab(expression(paste(mu, "g CO"[2],"-C minute"^-1, " g"^-1, "litter" ))) +
#scale_y_continuous(labels = fancy_scientific)+
scale_colour_manual(labels = c("#E69F00" = "Plot 1", "#56B4E9"="Plot 2", "#009E73"="Plot 3"),
    values = c("#E69F00" = cbPalette[2], "#56B4E9"=cbPalette[3], "#009E73"=cbPalette[4]),
    breaks = c("#E69F00", "#56B4E9", "#009E73")) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
my_gg_theme
dev.off()



########
#Respiration by cell_abd
########
resp_p_cells_sin_plot_lm = lm(resp.per.cell ~ Plot/sin(2*pi*hours.cumulative.temp/24) + Plot/cos(2*pi*hours.cumulative.temp/24), data = diel_phys_dat[-35,], na.action = na.exclude)

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










#Models resids comparisons


 summary(lm(residuals(nls.mod.resp) ~ residuals(nls.mod.moist)[1:84]))

plot(residuals(nls.mod.resp) ~ residuals(nls.mod.moist)[1:84])






#df = data.frame(moist.estimate = summary(nls.mod.moist)$coef[,1], moist.std.error = summary(nls.mod.moist)$coef[,2], moisture.p.value = summary(nls.mod.moist)$coef[,4])

###############################
#use this loop to extract data#
###############################

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.moist)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.moist)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.moist)$coef)[p]
        df$value[ii] = summary(nls.mod.moist)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.moist)$coef)[i]
        df$variable[ii] = "Moisture"
        df$plot[ii] = "Overall"
        ii=ii+1
    }
}
everything.df = data.frame(df)

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.moist.P1)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.moist.P1)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.moist.P1)$coef)[p]
        df$value[ii] = summary(nls.mod.moist.P1)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.moist.P1)$coef)[i]
        df$variable[ii] = "Moisture"
        df$plot[ii] = "P1"
        ii=ii+1
    }
}
#moisture.P1.df = data.frame(df)
everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.moist.P2)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.moist.P2)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.moist.P2)$coef)[p]
        df$value[ii] = summary(nls.mod.moist.P2)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.moist.P2)$coef)[i]
        df$variable[ii] = "Moisture"
        df$plot[ii] = "P2"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.moist.P3)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.moist.P3)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.moist.P3)$coef)[p]
        df$value[ii] = summary(nls.mod.moist.P3)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.moist.P3)$coef)[i]
        df$variable[ii] = "Moisture"
        df$plot[ii] = "P3"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.temp)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.temp)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.temp)$coef)[p]
        df$value[ii] = summary(nls.mod.temp)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.temp)$coef)[i]
        df$variable[ii] = "Temperature"
        df$plot[ii] = "Overall"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.temp.P1)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.temp.P1)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.temp.P1)$coef)[p]
        df$value[ii] = summary(nls.mod.temp.P1)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.temp.P1)$coef)[i]
        df$variable[ii] = "Temperature"
        df$plot[ii] = "P1"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.temp.P2)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.temp.P2)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.temp.P2)$coef)[p]
        df$value[ii] = summary(nls.mod.temp.P2)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.temp.P2)$coef)[i]
        df$variable[ii] = "Temperature"
        df$plot[ii] = "P2"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.temp.P3)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.temp.P3)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.temp.P3)$coef)[p]
        df$value[ii] = summary(nls.mod.temp.P3)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.temp.P3)$coef)[i]
        df$variable[ii] = "Temperature"
        df$plot[ii] = "P3"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.resp)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.resp)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.resp)$coef)[p]
        df$value[ii] = summary(nls.mod.resp)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.resp)$coef)[i]
        df$variable[ii] = "Respiration"
        df$plot[ii] = "Overall"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.resp.P1)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.resp.P1)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.resp.P1)$coef)[p]
        df$value[ii] = summary(nls.mod.resp.P1)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.resp.P1)$coef)[i]
        df$variable[ii] = "Respiration"
        df$plot[ii] = "P1"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.resp.P2)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.resp.P2)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.resp.P2)$coef)[p]
        df$value[ii] = summary(nls.mod.resp.P2)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.resp.P2)$coef)[i]
        df$variable[ii] = "Respiration"
        df$plot[ii] = "P2"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.resp.P3)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.resp.P3)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.resp.P3)$coef)[p]
        df$value[ii] = summary(nls.mod.resp.P3)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.resp.P3)$coef)[i]
        df$variable[ii] = "Respiration"
        df$plot[ii] = "P3"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.resp.cell)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.resp.cell)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.resp.cell)$coef)[p]
        df$value[ii] = summary(nls.mod.resp.cell)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.resp.cell)$coef)[i]
        df$variable[ii] = "Respiration per cell"
        df$plot[ii] = "Overall"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.resp.cell.P1)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.resp.cell.P1)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.resp.cell.P1)$coef)[p]
        df$value[ii] = summary(nls.mod.resp.cell.P1)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.resp.cell.P1)$coef)[i]
        df$variable[ii] = "Respiration per cell"
        df$plot[ii] = "P1"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.resp.cell.P2)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.resp.cell.P2)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.resp.cell.P2)$coef)[p]
        df$value[ii] = summary(nls.mod.resp.cell.P2)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.resp.cell.P2)$coef)[i]
        df$variable[ii] = "Respiration per cell"
        df$plot[ii] = "P2"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

df = NULL#data.frame(colnames(summary(nls.mod.moist)$coef))
ii= 1
for(i in 1:length(rownames(summary(nls.mod.resp.cell.P3)$coef)) )
{
    for(p in 1:length(colnames(summary(nls.mod.resp.cell.P3)$coef)) )
    {
        df$coef[ii] = colnames(summary(nls.mod.resp.cell.P3)$coef)[p]
        df$value[ii] = summary(nls.mod.resp.cell.P3)$coef[i,p]
        df$parameter[ii] = rownames(summary(nls.mod.resp.cell.P3)$coef)[i]
        df$variable[ii] = "Respiration per cell"
        df$plot[ii] = "P3"
        ii=ii+1
    }
}

everything.df = rbind(everything.df, data.frame(df))

#cast to wide form on coefs
everything.wide = dcast(everything.df,  variable + plot + parameter ~ coef, value.var = "value")
colnames(everything.wide)[5] = "p.value"
colnames(everything.wide)[6] = "Std.err"
for(i in 1:length(everything.wide$p.value))
{
    if(everything.wide$p.value[i] < 0.05)
    {
        everything.wide$sig[i] = "P < 0.05"
    }else{
        everything.wide$sig[i] = "n.s."
    }
}

#Plot coefs& variable == "Moisture"
pdf("martiny_diel_phys_dat/moisture_nls_estimates.pdf", width = 8, height = 4)
ggplot() +
geom_point(data = filter(everything.wide, plot != "Overall" & variable == "Moisture"),
    aes(x = plot, y = Estimate, color = sig),
    stat = "identity",
    position = position_dodge(width = 0.9),
    size = 3) +
geom_errorbar(data = filter(everything.wide, plot != "Overall" & variable == "Moisture"),
    aes(ymin = Estimate-Std.err, ymax = Estimate+Std.err, x = plot),
    stat = "identity",
    position = position_dodge(width = 0.9),
    width = 0.2) +
scale_color_manual(values = c("P < 0.05" = "blue", "n.s." = "dark grey"), breaks = c("P < 0.05", "n.s."), labels = c("P < 0.05" = expression(paste(italic("P"), "< 0.05")))) +
facet_wrap(~parameter, scales = "free_y" ) +
labs(x = "", title = "Moisture") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/temperature_nls_estimates.pdf", width = 6, height = 4)
ggplot() +
geom_point(data = filter(everything.wide, plot != "Overall" & variable == "Temperature"),
aes(x = plot, y = Estimate, color = sig),
stat = "identity",
position = position_dodge(width = 0.9),
size = 3) +
geom_errorbar(data = filter(everything.wide, plot != "Overall" & variable == "Temperature"),
aes(ymin = Estimate-Std.err, ymax = Estimate+Std.err, x = plot),
stat = "identity",
position = position_dodge(width = 0.9),
width = 0.2) +
scale_color_manual(values = c("P < 0.05" = "blue", "n.s." = "dark grey"), breaks = c("P < 0.05", "n.s."), labels = c("P < 0.05" = expression(paste(italic("P"), "< 0.05")))) +
facet_wrap(~parameter, scales = "free_y" ) +
labs(x = "", title = "Temperature") +
my_gg_theme
dev.off()

pdf("martiny_diel_phys_dat/respiration_nls_estimates.pdf", width = 8, height = 4)
ggplot() +
geom_point(data = filter(everything.wide, plot != "Overall" & variable == "Respiration"),
aes(x = plot, y = Estimate, color = sig),
stat = "identity",
position = position_dodge(width = 0.9),
size = 3) +
geom_errorbar(data = filter(everything.wide, plot != "Overall" & variable == "Respiration"),
aes(ymin = Estimate-Std.err, ymax = Estimate+Std.err, x = plot),
stat = "identity",
position = position_dodge(width = 0.9),
width = 0.2) +
scale_color_manual(values = c("P < 0.05" = "blue", "n.s." = "dark grey"), breaks = c("P < 0.05", "n.s."), labels = c("P < 0.05" = expression(paste(italic("P"), "< 0.05")))) +
facet_wrap(~parameter, scales = "free_y" ) +
labs(x = "", title = "Respiration per litter dry weight") +
my_gg_theme
dev.off()


pdf("martiny_diel_phys_dat/resp_per_cell_nls_estimates.pdf", width = 8, height = 4)
ggplot() +
geom_point(data = filter(everything.wide, plot != "Overall" & variable == "Respiration per cell"),
aes(x = plot, y = Estimate, color = sig),
stat = "identity",
position = position_dodge(width = 0.9),
size = 3) +
geom_errorbar(data = filter(everything.wide, plot != "Overall" & variable == "Respiration per cell"),
aes(ymin = Estimate-Std.err, ymax = Estimate+Std.err, x = plot),
stat = "identity",
position = position_dodge(width = 0.9),
width = 0.2) +
scale_color_manual(values = c("P < 0.05" = "blue", "n.s." = "dark grey"), breaks = c("P < 0.05", "n.s."), labels = c("P < 0.05" = expression(paste(italic("P"), "< 0.05")))) +
facet_wrap(~parameter, scales = "free_y" ) +
labs(x = "", title = "Respiration per cell") +
my_gg_theme
dev.off()



##################################################################
#Old stuff by maybe useful for fitting single models subset by plot

gnls.mod = gnls(Temperature ~ amp*sin(2*pi* hours.cumulative.temp /24 + phase) + C, groups = Plot, start = list(amp = r, phase = phi, C = 8), data = diel_phys_dat[-35], na.action = na.omit)

anova(nls.mod, nls.mod.2)

nls.coefs = coef(nls.mod.2)


f = function(x, amp, phase, period, C){amp*sin(2*pi*x/24+phase)+C+x}

plot(y ~ x.time)
curve(f(x, amp = nls.coefs["amp"], phase = nls.coefs["phase"], C = nls.coefs["C"]), add=TRUE ,lwd=2, )


summary(nls.mod)


####Estimating different levels of coefficients by dummy vars
#amp is amplitude, phase is offset, period is harcoded to 24 hours, C is harmonic mean

diel_phys_dat$P1 = as.numeric(diel_phys_dat$Plot=="P1")
diel_phys_dat$P2 = as.numeric(diel_phys_dat$Plot=="P2")
diel_phys_dat$P3 = as.numeric(diel_phys_dat$Plot=="P3")

nls.mod.2 = nls(Temperature ~ amp*sin(2*pi* hours.cumulative.temp /24 + phase) + C1*P1 + C2*P2 + C3*P3, start = list(amp = r, phase = phi, C1 = 8, C2 = 8, C3 = 8), data = diel_phys_dat[-35], na.action = na.omit)

#No phase dif
nls.mod.2 = nls(Temperature ~ amp1*C1*sin(2*pi* hours.cumulative.temp /24 + phase) + amp2*C2*sin(2*pi* hours.cumulative.temp /24 + phase) + amp3*C3*sin(2*pi* hours.cumulative.temp /24 + phase) + C1*P1 + C2*P2 + C3*P3, start = list(amp1 = r+2, amp2 = r, amp3 = r-2,phase = phi, C1 = 8, C2 = 8, C3 = 8), data = diel_phys_dat[-35], na.action = na.omit)

nls.mod.2 = nls(Temperature ~ I(amp1*C1)*I(amp2*C2)*I(amp3*C3)*sin(2*pi* hours.cumulative.temp /24 + phase)  + C1*P1 + C2*P2 + C3*P3, start = list(amp1 = r+1, amp2 = r, amp3 = r-1,phase = phi, C1 = 8, C2 = 8, C3 = 8), data = diel_phys_dat[-35], na.action = na.omit)



nls.mod.2 = nls(Temperature ~ amp*sin(2*pi* hours.cumulative.temp /24 + phase) + C1*P1 + C2*P2 + C3*P3, start = list(amp = r, phase = phi, C1 = 8, C2 = 8, C3 = 8), data = diel_phys_dat[-35], na.action = na.omit)

nlme(

