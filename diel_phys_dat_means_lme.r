require(plyr)
require(dplyr)
require(ggplot2)
require(nlme)

diel_phys_dat = read.table("~/martiny_diel_phys_dat/moisture.tep.res.rin.txt", header = T, sep = "\t")

diel_phys_dat$fixtime = strptime(as.character(diel_phys_dat$Time.bead.beat), '%m/%d/%Y %H:%M')
diel_phys_dat$temptime = strptime(as.character(diel_phys_dat$Temp.time), '%m/%d/%Y %H:%M')

mod1 = lme(fixed = Moisture ~ Plot, random = ~1|hours.cumulative.temp, data = diel_phys_dat)
anova(mod1)
mod2 = aov(Moisture ~ Plot, data = diel_phys_dat)
model.tables(mod2, "means")

mod1 = lme(fixed = Temperature ~ Plot, random = ~1|hours.cumulative.temp, data = diel_phys_dat)
anova(mod1)
mod2 = aov(Temperature ~ Plot, data = diel_phys_dat)
model.tables(mod2, "means")

mod1 = lme(fixed = cell.counts ~ Plot, random = ~1|hours.cumulative.temp, data = diel_phys_dat, na.action = "na.exclude")
anova(mod1)
mod2 = aov(cell.counts ~ Plot, data = diel_phys_dat)
model.tables(mod2, "means")


mod1 = lme(fixed = resp.per.dry.weight ~ Plot, random = ~1|hours.cumulative.temp, data = diel_phys_dat, na.action = "na.exclude")
anova(mod1)
mod2 = aov(resp.per.dry.weight ~ Plot, data = diel_phys_dat)
model.tables(mod2, "means")

mod1 = lme(fixed = resp.per.cell ~ Plot, random = ~1|hours.cumulative.temp, data = diel_phys_dat, na.action = "na.exclude")
anova(mod1)
mod2 = aov(resp.per.cell ~ Plot, data = diel_phys_dat)
model.tables(mod2, "means")
