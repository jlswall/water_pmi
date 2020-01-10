#mixed effects modeling
environment = read.table("/Volumes/LaCie/2019 Dissertation Data/2020_01_08_Combined_HL_RRC/2020_Dissertation_Combined/HL_RRC_Environmental.txt", header=TRUE, sep="\t")
require(reshape2)
environment2 = melt(environment, id.vars = c("Site","Collection","ADD","Days","Season"), variable.name = "Measure", value.name = "Parameter")
environment3 = na.omit(environment2)
environment4 = na.omit(environment)
require(nlme)
install.packages("multcomp")

#Attempt with lm 
env.lm = lm(Parameter ~ Site*Measure*ADD,
   data = environment3)
summary(env.lm)
anova(env.lm)



#Attempt with lme 
?lme
env.lme <-  lme(Parameter ~ Site*Measure,
                 data = environment3, random = ~1 | Collection)
anova(env.lme) #Site <.0001, Measure <.0001 and Site:Measure <.0001

#significant interaction between environmental measures and site 
#test each seperately 
summary(aov(Temp ~ Site,data=environment4)) #significant
summary(aov(SpCond ~ Site,data=environment4)) #significant
summary(aov(pH ~ Site,data=environment)) #significant
summary(aov(Sal ~ Site,data=environment)) #significant



