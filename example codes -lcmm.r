library(lcmm)
library(tidyr)
library(ggplot2)
library(dplyr)


set.seed(5)
sim_data <- read.csv(file.choose())  # import long_data.csv from RxPharmStat/adherence-sim-data repository at GitHub
View(sim_data)

# View data in plot. 
ggplot(sim_data, aes(time, Adherence, color=factor(PatientID))) + geom_point() + geom_line() + 
  theme(legend.position = "none")

# unconditional 1-4 class models
m1.l <- hlme(Adherence~time, subject = "PatientID", ng=1, data=sim_data)
m2.l <- hlme(Adherence~time, mixture= ~time, subject="PatientID", ng=2, data=sim_data, B=m1.l) 
m3.l <- hlme(Adherence~time, mixture= ~time, subject="PatientID", ng=3, data=sim_data, B=m1.l) 
m4.l <- hlme(Adherence~time, mixture= ~time, subject="PatientID", ng=4, data=sim_data, B=m1.l) 

summarytable(m1.l ,m2.l, m3.l, m4.l)  # show fit statistics
summarytable(m4.l, m4.q,  m4.c)
summaryplot(m1.l ,m2.l, m3.l, m4.l)

# show estimated trajectories from 1-4 class models
par(mfrow=c(1,4))
plot(m1.l, which="fit", var.time="time", shades=TRUE, font.lab = 2, legend=NULL)
plot(m2.l, which="fit", var.time="time", shades=TRUE, font.lab = 2, legend=NULL)
plot(m3.l, which="fit", var.time="time", shades=TRUE, font.lab = 2, legend=NULL)
plot(m4.l, which="fit", var.time="time", shades=TRUE, font.lab = 2, legend=NULL)

# add qudratic term  
m1.q <- hlme(Adherence~time+I(time^2), subject = "PatientID", ng=1, data=sim_data)
m4.q <- hlme(Adherence~time+I(time^2), mixture= ~time+I(time^2), subject="PatientID", ng=4, data=sim_data, B=m1.q) #converged

# add cubic term  
m1.c <- hlme(Adherence~time+I(time^2)+I(time^3), subject = "PatientID", ng=1, data=sim_data)
m4.c <- hlme(Adherence~time+I(time^2)+I(time^3), mixture= ~time+I(time^2)+I(time^3), subject="PatientID", ng=4, data=sim_data, B=m1.c) # not converged
m4.c <- hlme(Adherence~time+I(time^2)+I(time^3), mixture= ~time+I(time^2)+I(time^3), subject="PatientID", ng=4, data=sim_data, B=random(m1.c)) #converged after changing the initial values
summary(m4.c) # output 4-class model with cubic term

# show estimated trajectories from 4 class model with cubic term
plot(m4.c, which="fit", var.time="time", shades=TRUE, font.lab = 2)

# fit conditional model by adding 1 binary covaiate 
# linear 
m1.con.l <- hlme(Adherence~time+X1, subject = "PatientID", ng=1, data=sim_data)
m4.con.l <- hlme(Adherence~time+X1, mixture= ~time+X1, classmb= ~X1, subject="PatientID", ng=4, data=sim_data, B=m1.con.l) 

# quadratic 
m1.con.q <- hlme(Adherence~time+I(time^2)+X1, subject = "PatientID", ng=1, data=sim_data)
m4.con.q <- hlme(Adherence~time+I(time^2)+X1, mixture= ~time+I(time^2), classmb= ~X1,  subject="PatientID", ng=4, data=sim_data, B=m1.con.q)

# cubic
m1.con.c <- hlme(Adherence~time+I(time^2)+I(time^3)+X1, subject = "PatientID", ng=1, data=sim_data)
m4.con.c <- hlme(Adherence~time+I(time^2)+I(time^3)+X1, mixture= ~time+I(time^2)+I(time^3)+X1, classmb= ~X1, subject="PatientID", ng=4, data=sim_data, B=random(m1.con.c)) #converged
summary(m4.con.c)
summarytable(m4.con.c)

# add random effects in conditional 4-class model
m1.con.ran.c <- hlme(Adherence~time+I(time^2)+I(time^3)+X1,random=~time, subject = "PatientID", ng=1, data=sim_data)
m4.con.ran.c <- hlme(Adherence~time+I(time^2)+I(time^3)+X1, mixture= ~time+I(time^2)+I(time^3)+X1, classmb= ~X1,random=~time, subject="PatientID", ng=4, data=sim_data, B=random(m1.con.ran.c)) #converged
summary(m4.con.ran.c)

# add random effects in unconditional 4-class model
m4.ran.c <- hlme(Adherence~time+I(time^2)+I(time^3)+X1, mixture= ~time+I(time^2)+I(time^3)+X1, random=~time, subject="PatientID", ng=4, data=sim_data, B=random(m1.con.ran.c))

# compare 4-class models: unconditional without random effects vs. conditional without random effects vs. conditional with random effects
summarytable(m4.c, m4.con.c, m4.con.ran.c)


