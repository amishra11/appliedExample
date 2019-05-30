### Anu Mishra ###
### HIV Predictors in young pop ###
### 6/1/17 ####

library(survival)
library(epiR)

### Purpose of Code: #####
## 1) Baseline characteristics for MTN-020
## 2) Predictors by Regions (and chi-sq)
## 3) Incidence rates and HR for MTN-020

######## Updated: 3-16-18 Redo MTN South Africa analyses for other age groups #####
## 4) Other Age analyses for MTN-020

mtn020 <- read.csv("/Volumes/amishra$/HIV RS - Young Pop - predictors/adata/mtn020hiv.csv")

table(mtn020$censor)
table(mtn020$censor12mon)

################## 1: Baseline characteristics for MTN-020 ######################
#age
mtn.age <-  cbind(table(mtn020$age),round(100*prop.table(table(mtn020$age)),2))
c(median(mtn020$age,na.rm=T),summary(mtn020$age)[c(2,5)])
table(mtn020$age,exclude=NULL)
sum(is.na(mtn020$age))

#education
mtn020.edu <- cbind(table(mtn020$educ),round(100*prop.table(table(mtn020$educ)),2))
sum(is.na(mtn020$educ))

#alcohol
mtn020$drinkcat <- NULL
mtn020$drinkcat[mtn020$drink==0] <- 0
mtn020$drinkcat[mtn020$drink>0 & mtn020$drink <= 3] <- 1
mtn020$drinkcat[mtn020$drink>3 & mtn020$drink <= 6] <- 2
mtn020$drinkcat[mtn020$drink>6] <- 3
mtn020$drinkcat <- factor(mtn020$drinkcat)

mtn020.drink <- cbind(table(mtn020$drinkcat),round(100*prop.table(table(mtn020$drinkcat)),2))
sum(is.na(mtn020$drink))

#income
mtn020.income <- cbind(sum(mtn020$demincom=="Yes"),round(100*mean(ifelse(mtn020$demincom=="Yes",1,0)),2))
sum(is.na(mtn020$demincom))

#married
mtn020.married <- cbind(sum(mtn020$demmerrd==1,na.rm=T),round(100*mean(mtn020$demmerrd,na.rm=T),2))
sum(is.na(mtn020$demmerrd))

#primary partner
mtn020.part <-  cbind(sum(mtn020$primpart==1,na.rm=T),round(100*mean(mtn020$primpart,na.rm=T),2))
sum(is.na(mtn020$primpart))

#baseline sti
mtn020.trr <- cbind(sum(mtn020$trr==1,na.rm=T),round(100*mean(mtn020$trr,na.rm=T),2))
sum(is.na(mtn020$trr))
mtn020.ct <- cbind(sum(mtn020$ct==1,na.rm=T),round(100*mean(mtn020$ct,na.rm=T),2))
sum(is.na(mtn020$ct))
mtn020.gon <- cbind(sum(mtn020$gon==1,na.rm=T),round(100*mean(mtn020$gon,na.rm=T),2))
sum(is.na(mtn020$gon))
mtn020.syph <- cbind(sum(mtn020$syph==1,na.rm=T),round(100*mean(mtn020$syph,na.rm=T),2))
sum(is.na(mtn020$syph))

#anysti
mtn020$anysti <- ifelse(mtn020$trr==1 | mtn020$ct==1 | mtn020$gon==1 | mtn020$syph==1, 1,0)
mtn020.anysti <- cbind(sum(mtn020$anysti==1,na.rm=T),round(100*mean(mtn020$anysti,na.rm=T),2))
sum(is.na(mtn020$anysti))

#vaginal ph
mtn020$vagPH_cat <- NULL
mtn020$vagPH_cat[mtn020$vagPH< 4.5] <- 1
mtn020$vagPH_cat[mtn020$vagPH>= 4.5 & mtn020$vagPH < 5.5] <- 2
mtn020$vagPH_cat[mtn020$vagPH>= 5.5] <- 3
mtn020.ph <-cbind(table(mtn020$vagPH_cat),round(100*prop.table(table(mtn020$vagPH_cat)),2))
sum(is.na(mtn020$vagPH_cat))
mtn020$vagPH_cat <- factor(mtn020$vagPH_cat)

#reversible contraceptive
mtn020.contra <- cbind(sum(mtn020$revContra==1,na.rm=T),round(100*mean(mtn020$revContra,na.rm=T),2))
sum(is.na(mtn020$revContra))
  
#parity
table(mtn020$BFPprgct)
100*prop.table(table(mtn020$BFPprgct))
mtn020$pregcat <- NULL
mtn020$pregcat[mtn020$BFPprgct==0] <- 0
mtn020$pregcat[mtn020$BFPprgct>0 & mtn020$BFPprgct <= 3] <- 1
mtn020$pregcat[mtn020$BFPprgct>3] <- 2
mtn020.preg <- cbind(table(mtn020$pregcat),round(100*prop.table(table(mtn020$pregcat)),2))
sum(is.na(mtn020$BFPprgct))
mtn020$pregcat <- factor(mtn020$pregcat)

#live births
table(mtn020$liveBirth)
100*prop.table(table(mtn020$liveBirth))
mtn020$liveBirth_relev[mtn020$liveBirth==0] <- 0
mtn020$liveBirth_relev[mtn020$liveBirth>0] <- 1
mtn020.livebirths <- cbind(table(mtn020$liveBirth_relev),round(100*prop.table(table(mtn020$liveBirth_relev)),2))
sum(is.na(mtn020$liveBirth_relev))
mtn020$liveBirth_relev <- factor(mtn020$liveBirth_relev)

#hiv partner
mtn020.parthiv <- cbind(table(mtn020$parthiv),round(100*prop.table(table(mtn020$parthiv)),2))
sum(is.na(mtn020$parthiv))

#circumcision
mtn020.parcirc <- cbind(table(mtn020$partcirc),round(100*prop.table(table(mtn020$partcirc)),2))
sum(is.na(mtn020$partcirc))

#same partner for last 3 mons
mtn020.samepart <- cbind(sum(mtn020$part3mon==1,na.rm=T),round(100*mean(mtn020$part3mon,na.rm=T),2))
sum(is.na(mtn020$part3mon))

#number partners last 3 months
mtn020$BBAnopsp <- as.numeric(mtn020$BBAnopsp)
mtn020$numPartcat <- NULL
mtn020$numPartcat[mtn020$BBAnopsp==0] <- 1
mtn020$numPartcat[mtn020$BBAnopsp>0 & mtn020$BBAnopsp <= 1] <- 2
mtn020$numPartcat[mtn020$BBAnopsp>1] <- 3
mtn020.numPartcat <- cbind(table(mtn020$numPartcat),round(100*prop.table(table(mtn020$numPartcat)),2))
sum(is.na(mtn020$BBAnopsp))
mtn020$BBAnopsp <- factor(mtn020$BBAnopsp)

#sex acts in last 3 months
mtn020$BBA3mvsn[mtn020$BBA3mvsn==99] <- NA
mtn020.vagSex <- c(median(mtn020$BBA3mvsn,na.rm=T),summary(mtn020$BBA3mvsn)[c(2,5)])
sum(is.na(mtn020$BBA3mvsn))

mtn020$analSex_cat <- NULL
mtn020$analSex_cat[mtn020$BBA3masn==0] <- 1
mtn020$analSex_cat[mtn020$BBA3masn!=0] <- 2
mtn020.analSex <- cbind(table(mtn020$analSex_cat),round(100*prop.table(table(mtn020$analSex_cat)),2))
sum(is.na(mtn020$analSex_cat))
mtn020$BBA3masn <- factor(mtn020$BBA3masn)
  
#worried about ring use
mtn020.ringworry <- cbind(table(mtn020$ringworry),round(100*prop.table(table(mtn020$ringworry)),2))
sum(is.na(mtn020$ringworry))

#vaginal washing
mtn020.vagwash <- cbind(sum(mtn020$vagwash==1,na.rm=T),round(100*mean(mtn020$vagwash,na.rm=T),2))
sum(is.na(mtn020$vagwash))

#transactional sex
mtn020.trans <- cbind(sum(mtn020$qexch==1,na.rm=T),round(100*mean(mtn020$qexch,na.rm=T),2))
sum(is.na(mtn020$qexch))

#other partners
mtn020.otherpart <- cbind(table(mtn020$otherpart),round(100*prop.table(table(mtn020$otherpart)),2))
sum(is.na(mtn020$otherpart))

#worried about HIV
mtn020.hivworry <- cbind(table(mtn020$hivworry),round(100*prop.table(table(mtn020$hivworry)),2))[1:3,]
sum(is.na(mtn020$hivworry))

#region 
mtn020$region1 <- relevel(mtn020$region1,ref="Durban")
mtn020.region1 <- cbind(table(mtn020$region1),round(100*prop.table(table(mtn020$region1)),2))

levels(mtn020$region2) <-  c("North Durban","West Durban","South Durban","Cape Town","Johannesburg")
mtn020.region2 <- cbind(table(mtn020$region2),round(100*prop.table(table(mtn020$region2)),2))



tab1 <- rbind(mtn.age,mtn020.edu,mtn020.drink,mtn020.income,mtn020.married,mtn020.part,mtn020.trr,mtn020.ct,mtn020.gon,
      mtn020.syph,mtn020.anysti,mtn020.ph,mtn020.contra,mtn020.preg,mtn020.livebirths,mtn020.parthiv,mtn020.parcirc,mtn020.samepart,
      mtn020.numPartcat,mtn020.analSex,mtn020.ringworry,mtn020.vagwash,mtn020.trans,mtn020.otherpart,
      mtn020.hivworry,mtn020.region1,mtn020.region2)

#write.csv(tab1,"/Users/anumishra/Dropbox/RA/HIV predictors - Young SA pop/raw results/tab1.csv")

########################## 2: Descriptives by Region ###########################
age <- rbind(c(median(mtn020$age[mtn020$region1=="Durban"]), summary(mtn020$age[mtn020$region1=="Durban"])[c(2,5)]),
             c(median(mtn020$age[mtn020$region1=="Cape Town"]), summary(mtn020$age[mtn020$region1=="Cape Town"])[c(2,5)]),
             c(median(mtn020$age[mtn020$region1=="Johannesburg"]), summary(mtn020$age[mtn020$region1=="Johannesburg"])[c(2,5)]))

eduRegion <- cbind(table(mtn020$educ[mtn020$region1=="Durban"]),
                   round(100*prop.table(table(mtn020$educ[mtn020$region1=="Durban"])),2),
                   table(mtn020$educ[mtn020$region1=="Cape Town"]),
                   round(100*prop.table(table(mtn020$educ[mtn020$region1=="Cape Town"])),2),
                   table(mtn020$educ[mtn020$region1=="Johannesburg"]),
                   round(100*prop.table(table(mtn020$educ[mtn020$region1=="Johannesburg"])),2))


drinkRegion <- cbind(table(mtn020$drinkcat[mtn020$region1=="Durban"]),
                     round(100*prop.table(table(mtn020$drinkcat[mtn020$region1=="Durban"])),2),
                     table(mtn020$drinkcat[mtn020$region1=="Cape Town"]),
                     round(100*prop.table(table(mtn020$drinkcat[mtn020$region1=="Cape Town"])),2),
                     table(mtn020$drinkcat[mtn020$region1=="Johannesburg"]),
                     round(100*prop.table(table(mtn020$drinkcat[mtn020$region1=="Johannesburg"])),2))

incomeRegion <- cbind(table(mtn020$demincom[mtn020$region1=="Durban"]),
                      round(100*prop.table(table(mtn020$demincom[mtn020$region1=="Durban"])),2),
                      table(mtn020$demincom[mtn020$region1=="Cape Town"]),
                      round(100*prop.table(table(mtn020$demincom[mtn020$region1=="Cape Town"])),2),
                      table(mtn020$demincom[mtn020$region1=="Johannesburg"]),
                      round(100*prop.table(table(mtn020$demincom[mtn020$region1=="Johannesburg"])),2))


marriedRegion <- cbind(table(mtn020$demmerrd[mtn020$region1=="Durban"]),
                       round(100*prop.table(table(mtn020$demmerrd[mtn020$region1=="Durban"])),2),
                       table(mtn020$demmerrd[mtn020$region1=="Cape Town"]),
                       round(100*prop.table(table(mtn020$demmerrd[mtn020$region1=="Cape Town"])),2),
                       table(mtn020$demmerrd[mtn020$region1=="Johannesburg"]),
                       round(100*prop.table(table(mtn020$demmerrd[mtn020$region1=="Johannesburg"])),2))

pspRegion <- cbind(table(mtn020$primpart[mtn020$region1=="Durban"]),
                   round(100*prop.table(table(mtn020$primpart[mtn020$region1=="Durban"])),2),
                   table(mtn020$primpart[mtn020$region1=="Cape Town"]),
                   round(100*prop.table(table(mtn020$primpart[mtn020$region1=="Cape Town"])),2),
                   table(mtn020$primpart[mtn020$region1=="Johannesburg"]),
                   round(100*prop.table(table(mtn020$primpart[mtn020$region1=="Johannesburg"])),2))

trrRegion <- cbind(table(mtn020$trr[mtn020$region1=="Durban"]),
                   round(100*prop.table(table(mtn020$trr[mtn020$region1=="Durban"])),2),
                   table(mtn020$trr[mtn020$region1=="Cape Town"]),
                   round(100*prop.table(table(mtn020$trr[mtn020$region1=="Cape Town"])),2),
                   table(mtn020$trr[mtn020$region1=="Johannesburg"]),
                   round(100*prop.table(table(mtn020$trr[mtn020$region1=="Johannesburg"])),2))

ctRegion <- cbind(table(mtn020$ct[mtn020$region1=="Durban"]),
                  round(100*prop.table(table(mtn020$ct[mtn020$region1=="Durban"])),2),
                  table(mtn020$ct[mtn020$region1=="Cape Town"]),
                  round(100*prop.table(table(mtn020$ct[mtn020$region1=="Cape Town"])),2),
                  table(mtn020$ct[mtn020$region1=="Johannesburg"]),
                  round(100*prop.table(table(mtn020$ct[mtn020$region1=="Johannesburg"])),2))

gonRegion <- cbind(table(mtn020$gon[mtn020$region1=="Durban"]),
                   round(100*prop.table(table(mtn020$gon[mtn020$region1=="Durban"])),2),
                   table(mtn020$gon[mtn020$region1=="Cape Town"]),
                   round(100*prop.table(table(mtn020$gon[mtn020$region1=="Cape Town"])),2),
                   table(mtn020$gon[mtn020$region1=="Johannesburg"]),
                   round(100*prop.table(table(mtn020$gon[mtn020$region1=="Johannesburg"])),2))

syphRegion <- cbind(table(mtn020$syph[mtn020$region1=="Durban"]),
                    round(100*prop.table(table(mtn020$syph[mtn020$region1=="Durban"])),2),
                    table(mtn020$syph[mtn020$region1=="Cape Town"]),
                    round(100*prop.table(table(mtn020$syph[mtn020$region1=="Cape Town"])),2),
                    table(mtn020$syph[mtn020$region1=="Johannesburg"]),
                    round(100*prop.table(table(mtn020$syph[mtn020$region1=="Johannesburg"])),2))

anystiRegion <- cbind(table(mtn020$anysti[mtn020$region1=="Durban"]),
                      round(100*prop.table(table(mtn020$anysti[mtn020$region1=="Durban"])),2),
                      table(mtn020$anysti[mtn020$region1=="Cape Town"]),
                      round(100*prop.table(table(mtn020$anysti[mtn020$region1=="Cape Town"])),2),
                      table(mtn020$anysti[mtn020$region1=="Johannesburg"]),
                      round(100*prop.table(table(mtn020$anysti[mtn020$region1=="Johannesburg"])),2))


phRegion <- cbind(table(mtn020$vagPH_cat[mtn020$region1=="Durban"]),
                 round(100*prop.table(table(mtn020$vagPH_cat[mtn020$region1=="Durban"])),2),
                 table(mtn020$vagPH_cat[mtn020$region1=="Cape Town"]),
                 round(100*prop.table(table(mtn020$vagPH_cat[mtn020$region1=="Cape Town"])),2),
                 table(mtn020$vagPH_cat[mtn020$region1=="Johannesburg"]),
                 round(100*prop.table(table(mtn020$vagPH_cat[mtn020$region1=="Johannesburg"])),2))

contraRegion <- cbind(table(mtn020$revContra[mtn020$region1=="Durban"]),
                      round(100*prop.table(table(mtn020$revContra[mtn020$region1=="Durban"])),2),
                      table(mtn020$revContra[mtn020$region1=="Cape Town"]),
                      round(100*prop.table(table(mtn020$revContra[mtn020$region1=="Cape Town"])),2),
                      table(mtn020$revContra[mtn020$region1=="Johannesburg"]),
                      round(100*prop.table(table(mtn020$revContra[mtn020$region1=="Johannesburg"])),2))

pregRegion <- cbind(table(mtn020$pregcat[mtn020$region1=="Durban"]),
                    round(100*prop.table(table(mtn020$pregcat[mtn020$region1=="Durban"])),2),
                    table(mtn020$pregcat[mtn020$region1=="Cape Town"]),
                    round(100*prop.table(table(mtn020$pregcat[mtn020$region1=="Cape Town"])),2),
                    table(mtn020$pregcat[mtn020$region1=="Johannesburg"]),
                    round(100*prop.table(table(mtn020$pregcat[mtn020$region1=="Johannesburg"])),2))

livebirthRegion <- cbind(table(mtn020$liveBirth_relev[mtn020$region1=="Durban"]),
                         round(100*prop.table(table(mtn020$liveBirth_relev[mtn020$region1=="Durban"])),2),
                         table(mtn020$liveBirth_relev[mtn020$region1=="Cape Town"]),
                         round(100*prop.table(table(mtn020$liveBirth_relev[mtn020$region1=="Cape Town"])),2),
                         table(mtn020$liveBirth_relev[mtn020$region1=="Johannesburg"]),
                         round(100*prop.table(table(mtn020$liveBirth_relev[mtn020$region1=="Johannesburg"])),2))


hivpartRegion <- cbind(table(mtn020$parthiv[mtn020$region1=="Durban"]),
                       round(100*prop.table(table(mtn020$parthiv[mtn020$region1=="Durban"])),2),
                       table(mtn020$parthiv[mtn020$region1=="Cape Town"]),
                       round(100*prop.table(table(mtn020$parthiv[mtn020$region1=="Cape Town"])),2),
                       table(mtn020$parthiv[mtn020$region1=="Johannesburg"]),
                       round(100*prop.table(table(mtn020$parthiv[mtn020$region1=="Johannesburg"])),2))

partcircRegion <- cbind(table(mtn020$partcirc[mtn020$region1=="Durban"]),
                        round(100*prop.table(table(mtn020$partcirc[mtn020$region1=="Durban"])),2),
                        c(table(mtn020$partcirc[mtn020$region1=="Cape Town"]),0),
                        c(round(100*prop.table(table(mtn020$partcirc[mtn020$region1=="Cape Town"])),2),0),
                        table(mtn020$partcirc[mtn020$region1=="Johannesburg"]),
                        round(100*prop.table(table(mtn020$partcirc[mtn020$region1=="Johannesburg"])),2))
                        
samepartRegion <- cbind(table(mtn020$part3mon[mtn020$region1=="Durban"]),
                        round(100*prop.table(table(mtn020$part3mon[mtn020$region1=="Durban"])),2),
                        table(mtn020$part3mon[mtn020$region1=="Cape Town"]),
                        round(100*prop.table(table(mtn020$part3mon[mtn020$region1=="Cape Town"])),2),
                        table(mtn020$part3mon[mtn020$region1=="Johannesburg"]),
                        round(100*prop.table(table(mtn020$part3mon[mtn020$region1=="Johannesburg"])),2))

part3monRegion <- cbind(table(mtn020$numPartcat[mtn020$region1=="Durban"]),
                        round(100*prop.table(table(mtn020$numPartcat[mtn020$region1=="Durban"])),2),
                        table(mtn020$numPartcat[mtn020$region1=="Cape Town"]),
                        round(100*prop.table(table(mtn020$numPartcat[mtn020$region1=="Cape Town"])),2),
                        c(table(mtn020$numPartcat[mtn020$region1=="Johannesburg"]),0),
                        c(round(100*prop.table(table(mtn020$numPartcat[mtn020$region1=="Johannesburg"])),2),0))

sexRegion <- rbind(c(median(mtn020$BBA3mvsn[mtn020$region1=="Durban"]), summary(mtn020$BBA3mvsn[mtn020$region1=="Durban"])[c(2,5)]),
             c(median(mtn020$BBA3mvsn[mtn020$region1=="Cape Town"]), summary(mtn020$BBA3mvsn[mtn020$region1=="Cape Town"])[c(2,5)]),
             c(median(mtn020$BBA3mvsn[mtn020$region1=="Johannesburg"]), summary(mtn020$BBA3mvsn[mtn020$region1=="Johannesburg"])[c(2,5)]))


analRegion <- cbind(table(mtn020$analSex_cat[mtn020$region1=="Durban"]),
                    round(100*prop.table(table(mtn020$analSex_cat[mtn020$region1=="Durban"])),2),
                    table(mtn020$analSex_cat[mtn020$region1=="Cape Town"]),
                    round(100*prop.table(table(mtn020$analSex_cat[mtn020$region1=="Cape Town"])),2),
                    table(mtn020$analSex_cat[mtn020$region1=="Johannesburg"]),
                    round(100*prop.table(table(mtn020$analSex_cat[mtn020$region1=="Johannesburg"])),2))

ringworryRegion <- cbind(table(mtn020$ringworry[mtn020$region1=="Durban"]),
                         round(100*prop.table(table(mtn020$ringworry[mtn020$region1=="Durban"])),2),
                         c(table(mtn020$ringworry[mtn020$region1=="Cape Town"]),0),
                         c(round(100*prop.table(table(mtn020$ringworry[mtn020$region1=="Cape Town"])),2),0),
                         c(table(mtn020$ringworry[mtn020$region1=="Johannesburg"]),0),
                         c(round(100*prop.table(table(mtn020$ringworry[mtn020$region1=="Johannesburg"])),2),0))

vagwashRegion <- cbind(table(mtn020$vagwash[mtn020$region1=="Durban"]),
                       round(100*prop.table(table(mtn020$vagwash[mtn020$region1=="Durban"])),2),
                       table(mtn020$vagwash[mtn020$region1=="Cape Town"]),
                       round(100*prop.table(table(mtn020$vagwash[mtn020$region1=="Cape Town"])),2),
                       table(mtn020$vagwash[mtn020$region1=="Johannesburg"]),
                       round(100*prop.table(table(mtn020$vagwash[mtn020$region1=="Johannesburg"])),2))

transRegion <- cbind(table(mtn020$qexch[mtn020$region1=="Durban"]),
                     round(100*prop.table(table(mtn020$qexch[mtn020$region1=="Durban"])),2),
                     table(mtn020$qexch[mtn020$region1=="Cape Town"]),
                     round(100*prop.table(table(mtn020$qexch[mtn020$region1=="Cape Town"])),2),
                     table(mtn020$qexch[mtn020$region1=="Johannesburg"])[1:2],
                     round(100*prop.table(table(mtn020$qexch[mtn020$region1=="Johannesburg"])),2)[1:2])

otherRegion <- cbind(table(mtn020$otherpart[mtn020$region1=="Durban"]),
                     round(100*prop.table(table(mtn020$otherpart[mtn020$region1=="Durban"])),2),
                     table(mtn020$otherpart[mtn020$region1=="Cape Town"]),
                     round(100*prop.table(table(mtn020$otherpart[mtn020$region1=="Cape Town"])),2),
                     table(mtn020$otherpart[mtn020$region1=="Johannesburg"]),
                     round(100*prop.table(table(mtn020$otherpart[mtn020$region1=="Johannesburg"])),2))

hivworryRegion <- cbind(table(mtn020$hivworry[mtn020$region1=="Durban"]),
                        round(100*prop.table(table(mtn020$hivworry[mtn020$region1=="Durban"])),2),
                        table(mtn020$hivworry[mtn020$region1=="Cape Town"]),
                        round(100*prop.table(table(mtn020$hivworry[mtn020$region1=="Cape Town"])),2),
                        table(mtn020$hivworry[mtn020$region1=="Johannesburg"]),
                        round(100*prop.table(table(mtn020$hivworry[mtn020$region1=="Johannesburg"])),2))

tab2 <- rbind(eduRegion,drinkRegion,incomeRegion,marriedRegion,pspRegion,phRegion,contraRegion,pregRegion,hivpartRegion,partcircRegion,
              samepartRegion,part3monRegion,analRegion,ringworryRegion,vagwashRegion,transRegion,
              otherRegion,hivworryRegion)

#write.csv(tab2,"/Users/anumishra/Dropbox/RA/HIV predictors - Young SA pop/raw results/tab2.csv")


eduFish <- fisher.test(table(mtn020$educ,mtn020$region1))$p.value
drinkFish <- fisher.test(table(mtn020$drinkcat,mtn020$region1))$p.value
incomFish <- fisher.test(table(mtn020$demincom,mtn020$region1))$p.value
marriedFish <- fisher.test(table(mtn020$demmerrd,mtn020$region1))$p.value
primpartFish <- fisher.test(table(mtn020$primpart,mtn020$region1))$p.value
stiFish <- fisher.test(table(mtn020$anysti,mtn020$region1))$p.value
phFish <- fisher.test(table(mtn020$vagPH_cat,mtn020$region1))$p.value
contraFish <- fisher.test(table(mtn020$revContra,mtn020$region1))$p.value
pregFish <- fisher.test(table(mtn020$pregcat,mtn020$region1))$p.value
liveFish <- fisher.test(table(mtn020$liveBirth_relev,mtn020$region1))$p.value
parthivFish <- fisher.test(table(mtn020$parthiv,mtn020$region1))$p.value
partcircFish <- fisher.test(table(mtn020$partcirc,mtn020$region1))$p.value
sameFish <- fisher.test(table(mtn020$part3mon,mtn020$region1))$p.value
numPartFish <- fisher.test(table(mtn020$numPartcat,mtn020$region1))$p.value
analFish <- fisher.test(table(mtn020$analSex_cat,mtn020$region1))$p.value
ringFish <- fisher.test(table(mtn020$ringworry,mtn020$region1))$p.value
vagWashFish <- fisher.test(table(mtn020$vagwash,mtn020$region1))$p.value
transFish <- fisher.test(table(mtn020$qexch,mtn020$region1))$p.value
otherFish <- fisher.test(table(mtn020$otherpart,mtn020$region1))$p.value
hivWorryFish <- fisher.test(table(mtn020$hivworry,mtn020$region1))$p.value

tab2fish <- round(cbind(eduFish,drinkFish,incomFish,marriedFish,primpartFish,phFish,contraFish,pregFish,parthivFish,partcircFish,sameFish,numPartFish,
                  analFish,ringFish,vagWashFish,transFish,otherFish,hivWorryFish),4)
#write.csv(tab2fish,"/Users/anumishra/Dropbox/RA/HIV predictors - Young SA pop/raw results/tab2fish.csv")



########################## 3: Incidence rates and HR for MTN-020 ###########################
ir.fun <- function(x,event,fuTime,type){
  
  if(type=="number"){
    mod <- coxph(Surv(fuTime, event)~x)
    low <- round(exp(summary(mod)$coef[1] - (1.96*summary(mod)$coef[3])),2)
    upp <- round(exp(summary(mod)$coef[1] + (1.96*summary(mod)$coef[3])),2)
    
    hr <- as.matrix(paste(round(summary(mod)$coef[2],2)," (",low,",",upp,")",sep=""))
    
    res <- c("-","-",hr)
  }
    
  else if (type=="cat"){
    nevent <- table(x,event)[,2]
    py <- round(as.matrix(by(fuTime,x,sum))/12,2)
    
    events.py <- as.matrix(paste(nevent,"/",py,sep=""))
    irTmp <- round(100*epi.conf(cbind(nevent,py),ctype = "inc.rate",method="exact"),2)
    ir <- as.matrix(paste(irTmp[,1]," (",irTmp[,2],",",irTmp[,3],")",sep=""))
    
    if(0 %in% nevent){
      xtmp <- factor(x)
      
      if(length(levels(xtmp))<= 2){
        res <- cbind(events.py,ir,"-")
      }
      else{
        keepLev <- levels(xtmp)[which(nevent!=0)]
        xtmp[!(xtmp %in% keepLev)] <- NA
        xtmp2 <- factor(xtmp)
      
        mod <- coxph(Surv(fuTime, event)~xtmp2)
        low <- round(exp(summary(mod)$coef[,1] - (1.96*summary(mod)$coef[,3])),2)
        upp <- round(exp(summary(mod)$coef[,1] + (1.96*summary(mod)$coef[,3])),2)
      
        hrTmp <- rbind("-",as.matrix(paste(round(summary(mod)$coef[,2],2)," (",low,",",upp,")",sep="")))
        point <- which(table(x,event)[,2]==0)
        hr <- append(hrTmp,"-",(point-1))
        res <- cbind(events.py, ir, hr)
      }
    }
    else{
      x <- factor(x)
      mod <- coxph(Surv(fuTime, event)~x)
      low <- round(exp(summary(mod)$coef[,1] - (1.96*summary(mod)$coef[,3])),2)
      upp <- round(exp(summary(mod)$coef[,1] + (1.96*summary(mod)$coef[,3])),2)
      hr <- rbind("-",as.matrix(paste(round(summary(mod)$coef[,2],2)," (",low,",",upp,")",sep="")))
      res <- cbind(events.py, ir, hr)
      
    }
    
  }
  return(res)
}

#overall
nevent <- table(mtn020$censor)[2]
py <- sum(mtn020$fu_mos)/12
ir.allFU <- round(100*epi.conf(cbind(nevent,py),ctype = "inc.rate",method="exact"),2)



age.allFU <- ir.fun(x=mtn020$age,event=mtn020$censor,fuTime=mtn020$fu_mos,type="number")
age.cat.allFU <- ir.fun(x=mtn020$age,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")

mtn020$educ_relev <- ifelse(mtn020$educ=="post-secondary",1,0)
educ.allFU <- ir.fun(x=mtn020$educ_relev,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")

mtn020$drink_relv <- ifelse(mtn020$drinkcat==0,0,1)
alc.allFU <- ir.fun(x=mtn020$drink_relv,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")

income.allFU <- ir.fun(x=mtn020$demmerrd,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")

tv.allFU <- ir.fun(x=mtn020$trr,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
ct.allFU <- ir.fun(x=mtn020$ct,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
gon.allFU <- ir.fun(x=mtn020$gon,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
syph.allFU <- ir.fun(x=mtn020$syph,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
any.allFU <- ir.fun(x=mtn020$anysti,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")

ph.allFU <- ir.fun(x=mtn020$vagPH_cat,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")


mtn020$revContraNo <- ifelse(mtn020$revContra==0,1,0)
contra.allFU <- ir.fun(x=mtn020$revContraNo,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
preg.allFU <- ir.fun(x=mtn020$pregcat,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
live.allFU <- ir.fun(x=mtn020$liveBirth_relev,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
parthiv.allFU <- ir.fun(x=mtn020$parthiv,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
partcirc.allFU <- ir.fun(x=mtn020$partcirc,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
part3mon.allFU <- ir.fun(x=mtn020$part3mon,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")

mtn020$numpart_relev <- ifelse(mtn020$numPartcat>=2,1,0)
numpart.allFU <- ir.fun(x=mtn020$numpart_relev,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
vagacts.allFU <- ir.fun(x=mtn020$BBA3mvsn,event=mtn020$censor,fuTime=mtn020$fu_mos,type="number")
analacts.allFU <- ir.fun(x=mtn020$analSex_cat,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
ringworry.allFU <- ir.fun(x=mtn020$ringworry,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
vagwash.allFU <- ir.fun(x=mtn020$vagwash,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
transsex.allFU <- ir.fun(x=mtn020$qexch,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
otherPart.allFU <- ir.fun(x=mtn020$otherpart,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")


hivworry.allFU <- ir.fun(x=mtn020$hivworry,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
mtn020$region1 <- relevel(mtn020$region1, ref="Durban")
region1.allFU <- ir.fun(x=mtn020$region1,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")

mtn020$region2 <- relevel(mtn020$region2, ref="North Durban")
region2.allFU <- ir.fun(x=mtn020$region2,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")

tab3a <- rbind(age.allFU,educ.allFU,alc.allFU,income.allFU,
               tv.allFU,ct.allFU,gon.allFU,syph.allFU,any.allFU,
               ph.allFU,contra.allFU,preg.allFU,live.allFU,
               parthiv.allFU,
               partcirc.allFU,part3mon.allFU,numpart.allFU,vagacts.allFU,analacts.allFU,ringworry.allFU,
               vagwash.allFU,transsex.allFU,otherPart.allFU,hivworry.allFU,region1.allFU,region2.allFU)
#write.csv(tab3a,"/Users/anumishra/Dropbox/RA/HIV predictors - Young SA pop/raw results/tab3a.csv")

#overall
nevent <- table(mtn020$censor12)[2]
py <- sum(mtn020$fu_mon12)/12
ir <- round(100*epi.conf(cbind(nevent,py),ctype = "inc.rate",method="exact"),2)


age.12mos <- ir.fun(x=mtn020$age,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="number")
age.cat.12mos <- ir.fun(x=mtn020$age,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")


educ.12mos <- ir.fun(x=mtn020$educ_relev,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
alc.12mos <- ir.fun(x=mtn020$drink_relv,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
income.12mos <- ir.fun(x=mtn020$demmerrd,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
tv.12mos <- ir.fun(x=mtn020$trr,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
ct.12mos <- ir.fun(x=mtn020$ct,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
gon.12mos <- ir.fun(x=mtn020$gon,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
syph.12mos <- ir.fun(x=mtn020$syph,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
any.12mos  <- ir.fun(x=mtn020$anysti,event=mtn020$censor12mon,fuTime=mtn020$fu_mos,type="cat")
ph.12mos <- ir.fun(x=mtn020$vagPH_cat,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
contra.12mos <- ir.fun(x=mtn020$revContraNo,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
preg.12mos <- ir.fun(x=mtn020$pregcat,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
live.12mos <- ir.fun(x=mtn020$liveBirth_relev,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
parthiv.12mos <- ir.fun(x=mtn020$parthiv,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
partcirc.12mos <- ir.fun(x=mtn020$partcirc,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
part3mon.12mos <- ir.fun(x=mtn020$part3mon,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
numpart.12mos <- ir.fun(x=mtn020$numpart_relev,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
vagacts.12mos <- ir.fun(x=mtn020$BBA3mvsn,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="number")
analacts.12mos <- ir.fun(x=mtn020$analSex_cat,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
ringworry.12mos <- ir.fun(x=mtn020$ringworry,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
vagwash.12mos <- ir.fun(x=mtn020$vagwash,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
transsex.12mos <- ir.fun(x=mtn020$qexch,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
otherPart.12mos <- ir.fun(x=mtn020$otherpart,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
hivworry.12mos <- ir.fun(x=mtn020$hivworry,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
region1.12mos <- ir.fun(x=mtn020$region1,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")
region2.12mos <- ir.fun(x=mtn020$region2,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")

tab3b <- rbind(age.12mos,educ.12mos,alc.12mos,income.12mos,
               tv.12mos,ct.12mos,gon.12mos,syph.12mos,any.12mos,ph.12mos,contra.12mos,preg.12mos,parthiv.12mos,
               partcirc.12mos,part3mon.12mos,numpart.12mos,vagacts.12mos,analacts.12mos,ringworry.12mos,
               vagwash.12mos,transsex.12mos,otherPart.12mos,hivworry.12mos,region1.12mos,region2.12mos)
#write.csv(tab3b,"/Users/anumishra/Dropbox/RA/HIV predictors - Young SA pop/raw results/tab3b.csv")



########## 4: partner charateristics ############
addmargins(table(mtn020$otherpart,mtn020$numpart_relev,exclude=NULL))
round(100*prop.table(table(mtn020$otherpart,mtn020$numPart>1),margin=2),2)

ba <- read.csv("/Volumes/trials/mtn/p020/analysis/manuscripts/tertiary/sti_contraceptive/adata/ba_fmt.csv") #bevahior dataset
ba <- ba[,c("ptid","visit","BApsp","BA3mpsp")]
ba$ptid <- as.integer(gsub("-","",ba$ptid))
#change in partner 
u <- unique(mtn020$ptid)
partChg <- ba[ba$ptid %in% u,]
partChg$BApsp[partChg$BApsp==""] <- NA
partChg$BApsp[partChg$BApsp=="Not Available"] <- NA
partChg$primPart <- ifelse(partChg$BApsp=="Yes" ,1,0)

partChg$partChg[partChg$BA3mpsp=="" | partChg$BA3mpsp=="Not Available"] <- NA
partChg$partChg[is.na(partChg$BApsp) | partChg$BApsp=="No"] <- 99
partChg$partChg[partChg$BA3mpsp=="Yes"] <- 0
partChg$partChg[partChg$BA3mpsp=="No"] <- 1

#censoring partChg variable so FU is the same as mtn020 dataset
u <- unique(partChg$ptid)
for (i in 1:length(u)){
  pc.tmp <- partChg[partChg$ptid==u[i],]
  mtn.tmp <- mtn020[mtn020$ptid==u[i],]
  
  #any partner change during FU (before last FU Visi)
  partChgVis <- pc.tmp$partChg[pc.tmp$visit <= mtn.tmp$lastfuvisit]
  partChgVis12 <- pc.tmp$partChg[pc.tmp$visit <= min(mtn.tmp$lastfuvisit,1201)]
  
  mtn020$anyPartChg[mtn020$ptid==u[i]] <- ifelse(1 %in% partChgVis,1,0)
  mtn020$anyPartChg12[mtn020$ptid==u[i]] <- ifelse(1 %in% partChgVis12,1,0)
  
  #number of partner changes
  mtn020$numPartChg[mtn020$ptid==u[i]] <- sum(partChgVis==1)
  mtn020$numPartChg12[mtn020$ptid==u[i]] <- sum(partChgVis12==1)
  
  #if seroconverted was detected at visit where question asked did partner change occur?
  mtn020$seroPart[mtn020$ptid==u[i]] <- ifelse(mtn.tmp$lastfuvisit %in% pc.tmp$visit & mtn.tmp$censor==1,
                     pc.tmp$partChg[pc.tmp$visit==mtn.tmp$lastfuvisit],
                     99)
  lastfuVist12 <- min(mtn.tmp$lastfuvisit,1201)
  mtn020$seroPart12[mtn020$ptid==u[i]] <- ifelse(lastfuVist12 %in% pc.tmp$visit & mtn.tmp$censor12==1,
                     pc.tmp$partChg[pc.tmp$visit==mtn.tmp$lastfuvisit],
                     99)
  
  #change in partner 3 months after baseline
  mtn020$postBaseChg[mtn020$ptid==u[i]] <- ifelse(1 %in% pc.tmp$partChg[pc.tmp$visit <= 301],1,0)
}


#all FU
partChgAny <- ir.fun(x=mtn020$anyPartChg,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
partChg12 <- ir.fun(x=mtn020$anyPartChg12,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
#12 mon
partChg12.12 <- ir.fun(x=mtn020$anyPartChg12,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")


## Change in partner over time descriptives
table(mtn020$anyPartChg,exclude=NULL)
table(mtn020$anyPartChg12,exclude=NULL)
table(mtn020$postBaseChg,exclude=NULL)

table(mtn020$numPartChg,exclude=NULL)
table(mtn020$numPartChg12,exclude=NULL)

#women who had only one partner at baseline 
addmargins(table(mtn020$anyPartChg,mtn020$numpart_relev,exclude=NULL))
round(100*prop.table(table(mtn020$anyPartChg,mtn020$numpart_relev),margin=2),2)
addmargins(table(mtn020$anyPartChg12,mtn020$numpart_relev,exclude=NULL))
round(100*prop.table(table(mtn020$anyPartChg12,mtn020$numpart_relev),margin=2),2)
addmargins(table(mtn020$postBaseChg,mtn020$numpart_relev,exclude=NULL))
round(100*prop.table(table(mtn020$postBaseChg,mtn020$numpart_relev),margin=2),2)


######### 6: PAR for MTN-020 ONLY GROUP #########
## sti ##
#PAR
sti.iu.mtn <- round(100*epi.conf(cbind(table(mtn020$anysti,mtn020$censor)[,2],
                                       round(as.matrix(by(mtn020$fu_mos,mtn020$anysti,sum))/12,2)),
                                 ctype = "inc.rate",method="exact"),2)[1,1]
100*(ir.allFU[1]-sti.iu.mtn)/ir.allFU[1]

#PAF
sti.ir <- 100*epi.conf(cbind(table(mtn020$anysti,mtn020$censor)[,2],
                             round(as.matrix(by(mtn020$fu_mos,mtn020$anysti,sum))/12,2)),
                       ctype = "inc.rate",method="exact")[,1]
sti.rr <- sti.ir[2]/sti.ir[1]
p.sti <- mean(mtn020$anysti,na.rm=T)
100*(p.sti * ((sti.rr - 1)/(1 + p.sti*(sti.rr-1))))

## live births ##
#PAR
birth.iu.mtn <- round(100*epi.conf(cbind(table(mtn020$liveBirth_relev,mtn020$censor)[,2],
                                         round(as.matrix(by(mtn020$fu_mos,mtn020$liveBirth_relev,sum))/12,2)),
                                   ctype = "inc.rate",method="exact"),2)[1,1]
100*(ir.allFU[1]-birth.iu.mtn)/ir.allFU[1]

#PAF
birth.ir <- 100*epi.conf(cbind(table(mtn020$liveBirth_relev,mtn020$censor)[,2],
                               round(as.matrix(by(mtn020$fu_mos,mtn020$liveBirth_relev,sum))/12,2)),
                         ctype = "inc.rate",method="exact")[,1]
birth.rr <- birth.ir[2]/birth.ir[1]
p.birth <- mean(as.integer(as.character(mtn020$liveBirth_relev)),na.rm=T)
100*(p.birth * ((birth.rr - 1)/(1 + p.birth*(birth.rr-1))))

## other partners ##
mtn020$otherPartIND <- ifelse(mtn020$otherpart==0,0,1)
#PAR
otherPart.comb.iu <- round(100*epi.conf(cbind(table(mtn020$otherPartIND,mtn020$censor)[,2],
                                              round(as.matrix(by(mtn020$fu_mos,mtn020$otherPartIND,sum))/12,2)),
                                        ctype = "inc.rate",method="exact"),2)[1,1]
100*(ir.allFU[1]-otherPart.comb.iu)/ir.allFU[1]

#PAF
partInd.ir <- 100*epi.conf(cbind(table(mtn020$otherPartIND,mtn020$censor)[,2],
                                 round(as.matrix(by(mtn020$fu_mos,mtn020$otherPartIND,sum))/12,2)),
                           ctype = "inc.rate",method="exact")[,1]
partInd.rr <- partInd.ir[2]/partInd.ir[1]
p.partInd <- mean(mtn020$otherPartIND,na.rm=T)
100*(p.partInd * ((partInd.rr - 1)/(1 + p.partInd*(partInd.rr-1))))

## other partners PAF for individual levels ##
part.ir <- 100*epi.conf(cbind(table(mtn020$otherpart,mtn020$censor)[,2],
                                 round(as.matrix(by(mtn020$fu_mos,mtn020$otherpart,sum))/12,2)),
                           ctype = "inc.rate",method="exact")[,1]
part.rr <- part.ir/partInd.ir[1]
p.part <- prop.table(table(mtn020$otherpart))

#Yes
100*(p.part[2] * (part.rr[2] - 1) / (1 + (p.part[2]*(part.rr[2]-1)) ))
#Don't know
100*(p.part[3] * (part.rr[3] - 1) / (1 + (p.part[3]*(part.rr[3]-1)) ))


## >1 partner ##
#PAR
numpart.iu.mtn <- round(100*epi.conf(cbind(table(mtn020$numpart_relev,mtn020$censor)[,2],
                                           round(as.matrix(by(mtn020$fu_mos,mtn020$numpart_relev,sum))/12,2)),
                                     ctype = "inc.rate",method="exact"),2)[1,1]
100*(ir.allFU[1]-numpart.iu.mtn)/ir.allFU[1]

#PAF
numpart_relev.ir <- 100*epi.conf(cbind(table(mtn020$numpart_relev,mtn020$censor)[,2],
                                 round(as.matrix(by(mtn020$fu_mos,mtn020$numpart_relev,sum))/12,2)),
                           ctype = "inc.rate",method="exact")[,1]
numpart_relev.rr <- numpart_relev.ir[2]/numpart_relev.ir[1]
p.numpart_relev <- mean(mtn020$numpart_relev,na.rm=T)
100*(p.numpart_relev * ((numpart_relev.rr - 1)/(1 + p.numpart_relev*(numpart_relev.rr-1))))



## any partner change ##
#PAR
partchg.iu.mtn <- round(100*epi.conf(cbind(table(mtn020$anyPartChg,mtn020$censor)[,2],
                                           round(as.matrix(by(mtn020$fu_mos,mtn020$anyPartChg,sum))/12,2)),
                                     ctype = "inc.rate",method="exact"),2)[,1]
100*(ir.allFU[1]-partchg.iu.mtn)/ir.allFU[1]

#PAF
partChg.ir <- 100*epi.conf(cbind(table(mtn020$anyPartChg,mtn020$censor)[,2],
                                 round(as.matrix(by(mtn020$fu_mos,mtn020$anyPartChg,sum))/12,2)),
                           ctype = "inc.rate",method="exact")[,1]
partChg.rr <- partChg.ir[2]/partChg.ir[1]
p.partChg <- mean(mtn020$anyPartChg,na.rm=T)
100*(p.partChg * ((partChg.rr - 1)/(1 + p.partChg*(partChg.rr-1))))



## par/paf calculation functions
par.fun <- function(x,censor,futime,ir.all,levels){
  if(levels==1){
    #par
    iu <- round(100*epi.conf(cbind(table(x,censor)[,2],
                                   round(as.matrix(by(futime,x,sum))/12,2)),
                             ctype = "inc.rate",method="exact"),2)[,1]
    par <- 100*(ir.all-iu[1])/ir.all
    
    #paf
    rr <- iu[2]/iu[1]
    p <- prop.table(table(x))[2]
    paf <- 100*(p * ((rr - 1)/(1 + p*(rr-1))))
    
    res <- as.matrix(c(par,paf))
    return(res)
  }
  else if (levels==2){
    ir <- 100*epi.conf(cbind(table(x,censor)[,2],
                            round(as.matrix(by(futime,x,sum))/12,2)),
                            ctype = "inc.rate",method="exact")[,1]
    rr <- ir/ir[1]
    p <- prop.table(table(x))
    
    #Yes
    paf1 <- 100*(p[2] * (rr[2] - 1) / (1 + (p[2]*(rr[2]-1)) ))
    #Don't know
    paf2 <- 100*(p[3] * (rr[3] - 1) / (1 + (p[3]*(rr[3]-1)) ))
    
    res <- as.matrix(c(paf1,paf2))
    return(res)
  }
  
}


######## 7: Missing data and multiple imputation #######
#subsetting to only include predictors fitting model too (not original versions, etc)
mtn020$trans <- mtn020$qexch
mtn020$trans[mtn020$qexch==88] <- NA


mtn020.imp <- mtn020[c("ptid","age","educ_relev","drink_relv",
                        "anysti","vagPH_cat","revContraNo","pregcat",
                         "liveBirth_relev","parthiv","partcirc","part3mon",
                         "numpart_relev","BBA3mvsn","ringworry",
                         "vagwash","trans","otherpart","hivworry",
                         "region1","region2","anyPartChg","censor","censor12mon","fu_mos",
                         "fu_mon12")]

#cleaning up data for mice package
#binary factors need to be coded as 0/1 
mtn020.imp$parthiv[mtn020.imp$parthiv==0] <- 1 #no
mtn020.imp$parthiv[mtn020.imp$parthiv==88] <- 2 #don't know

mtn020.imp$partcirc[mtn020.imp$partcirc==88] <- 3 #don't know

mtn020.imp$hivworry[mtn020.imp$hivworry==3] <- 2 

#turning everything into factor for clean imputation
factorNames <- names(mtn020.imp) [!(names(mtn020.imp) %in% c("ptid","BBA3mvsn","censor","censor12mon","fu_mos","fu_mon12"))]
mtn020.imp[,factorNames] <- lapply(mtn020.imp[,factorNames],factor)
apply(mtn020.imp,2,function(x){sum(is.na(x))})


library(mice)
set.seed(604)
mi.dat <- mice(mtn020.imp,m=10,meth=c(rep("",9),"logreg","polyreg","logreg",
                                      "","","","logreg","logreg","polyreg",
                                      "polyreg","","","logreg","","","",""))


list.beta <- list("age.beta" <- NULL, "educ.beta" <- NULL, "drink.beta" <- NULL,
                  "anysti.beta" <- NULL, "vagPH.beta" <- NULL, 
                  "revContra.beta" <- NULL, "pregcat.beta" <- NULL,
                  "liveBirth.beta" <- NULL, "parthiv.beta"<- NULL,
                  "partcirc.beta" <- NULL, "part3mon.beta" <- NULL,
                  "numpart.beta" <- NULL, "numSex"<- NULL, "ringworry.beta" <-NULL, 
                  "vagwash.beta" <- NULL, "trans.beta"<-NULL, 
                  "otherpart.beta" <- NULL, "hivworry.beta" <- NULL,
                  "region1.beta" <- NULL, "region2.beta" <- NULL, "partChg" <- NULL)
names(list.beta) <- names(mtn020.imp[2:22])


list.beta12 <- list("age.beta" <- NULL, "educ.beta" <- NULL, "drink.beta" <- NULL,
                  "anysti.beta" <- NULL, "vagPH.beta" <- NULL, 
                  "revContra.beta" <- NULL, "pregcat.beta" <- NULL,
                  "liveBirth.beta" <- NULL, "parthiv.beta"<- NULL,
                  "partcirc.beta" <- NULL,
                  "numpart.beta" <- NULL, "numSex"<- NULL, "ringworry.beta" <-NULL, 
                  "vagwash.beta" <- NULL, "trans.beta"<-NULL, 
                  "otherpart.beta" <- NULL, "hivworry.beta" <- NULL,
                  "region1.beta" <- NULL, "region2.beta" <- NULL, "partChg" <- NULL)
names(list.beta12) <- names(mtn020.imp[c(2:11,13:22)])


#averaging across imputations
for(i in 1:10){
  dat <- complete(mi.dat,i)
  for(j in 2:22){
    mod <- coxph(Surv(dat$fu_mos,dat$censor)~dat[,j])
    list.beta[[j-1]] <- cbind(list.beta[[j-1]],exp(mod$coefficients))
  }
}
  

round(unlist(lapply(list.beta,function(x){apply(x,1,mean)})),2)

for(i in 1:10){
  dat <- complete(mi.dat,i)
  dat <- dat[,names(dat)!="part3mon"]
  for(j in c(2:21)){
    mod <- coxph(Surv(dat$fu_mon12,dat$censor12mon)~dat[,j])
    list.beta12[[j-1]] <- cbind(list.beta12[[j-1]],exp(mod$coefficients))
  }
}
round(unlist(lapply(list.beta12,function(x){apply(x,1,mean)})),2)


par.list <- list("sti" <- NULL, "livebirth" <- NULL, "numpart" <- NULL,
                  "partChg" <- NULL)
names(par.list) <- c("sti","livebirth","numpart","partChg")

for(i in 1:10){
  dat <- complete(mi.dat,i)
  dat <- dat[,c("anysti","liveBirth_relev","numpart_relev","anyPartChg","censor","fu_mos")]
  nevent <- table(dat$censor)[2]
  py <- sum(dat$fu_mos)/12
  ir.all <- round(100*epi.conf(cbind(nevent,py),ctype = "inc.rate",method="exact"),2)[1]
  for(j in 1:4){
    res <- as.matrix(par.fun(x = dat[,j],censor = dat$censor,futime = dat$fu_mos,ir.all = ir.all,levels = 1))
    par.list[[j]] <-cbind(par.list[[j]],res)
  }
}
lapply(par.list$sti[,1],mean)
lapply(par.list$livebirth[,1],mean)
lapply(par.list$numpart[,1],mean)
lapply(par.list$partChg[,1],mean)


otherPart.mat <- matrix(NA,nrow=2,ncol=10)
for(i in 1:10){
  dat <- complete(mi.dat,i)
  dat <- dat[,c("otherpart","censor","fu_mos")]
  nevent <- table(dat$censor)[2]
  py <- sum(dat$fu_mos)/12
  ir.all <- round(100*epi.conf(cbind(nevent,py),ctype = "inc.rate",method="exact"),2)[1]
  
  otherPart.mat[,i] <- par.fun(x = dat$otherpart,censor = dat$censor,futime = dat$fu_mos,ir.all = ir.all,levels = 2)
}
apply(otherPart.mat,1,mean)

##### Missing data table ########
miss <- is.na(mtn020$age) | is.na(mtn020$educ_relev) | is.na(mtn020$anysti) | is.na(mtn020$vagPH_cat) |
  is.na(mtn020$revContra) | is.na(mtn020$pregcat) | is.na(mtn020$liveBirth_relev) |
  is.na(mtn020$partcirc) | is.na(mtn020$part3mon) | is.na(mtn020$numpart_relev) | 
  is.na(mtn020$analSex_cat) | is.na(mtn020$ringworry) | is.na(mtn020$vagwash) |
  is.na(mtn020$qexch) | is.na(mtn020$otherpart) | is.na(mtn020$hivworry) | 
  is.na(mtn020$anyPartChg)
missDat <- mtn020[miss,]
dim(missDat)
table(missDat$region1)
round(100*(table(missDat$region1)/17),2)
sum(is.na(missDat$region1))

table(missDat$region2)
round(100*(table(missDat$region2)/17),2)
sum(is.na(missDat$region2))



###### Additional Analysis ##############
## Chi-squared of sti and live births 

addmargins(table(mtn020$liveBirth_relev,mtn020$anysti))
chisq.test(table(mtn020$liveBirth_relev,mtn020$anysti))

## adjusted model with live birth and sti
adj <- coxph(Surv(fu_mos, censor)~anysti+liveBirth_relev, data = mtn020)
summary(adj)

adj12 <- coxph(Surv(fu_mon12, censor12mon)~anysti+liveBirth_relev, data = mtn020)
summary(adj12)


## association between any pregnancy and HIV
pr <- read.csv("/Volumes/amishra$/HIV RS - Young Pop - predictors/adata/MTN_Data/pr.csv") #pregnancy outcomes
pr$ptid <- noquote(gsub("'",'',pr$Ptid))
pr$ptid <- as.integer(gsub("-","",pr$ptid))
pr$visit <- as.integer(noquote(gsub("'",'',pr$Visit)))
pr <- pr[,c("ptid","visit")]
#row 107 is missing so removing it
pr <- pr[-107,]


# matrix of each ptid that had pregancy report with first visit of pregnancy
firstPreg <- aggregate(pr$visit,by=list(pr$ptid),min)
names(firstPreg) <- c("ptid","visit")
mtn020$pregFU <- ifelse(mtn020$ptid %in% firstPreg$ptid,1,0)

#only keep those in firstPreg that are in mtn020
firstPregKeep <- firstPreg$ptid %in% mtn020$ptid
firstPregMTN <- firstPreg[firstPregKeep,]
mtn020$pregVis[mtn020$pregFU==1] <- firstPregMTN[,2]
mtn020$pregVis[mtn020$pregFU==0] <- 9999

#censoring partChg variable so FU is the same as mtn020 dataset
u <- unique(mtn020$ptid)
for (i in 1:length(u)){
  mtn.tmp <- mtn020[mtn020$ptid==u[i],]
  
  #any pregnancy during FU (before seroconversion)
  mtn020$pregFUsero[mtn020$ptid==u[i]] <- ifelse(mtn.tmp$pregVis <= mtn.tmp$lastfuvisit,1,0)
  mtn020$preg12sero[mtn020$ptid==u[i]] <- ifelse(mtn.tmp$pregVis <= min(mtn.tmp$lastfuvisit,1201),1,0)
}

#all FU
pregAny <- ir.fun(x=mtn020$pregFUsero,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
preg12 <- ir.fun(x=mtn020$preg12sero,event=mtn020$censor,fuTime=mtn020$fu_mos,type="cat")
#12 mon
preg12.12 <- ir.fun(x=mtn020$preg12sero,event=mtn020$censor12mon,fuTime=mtn020$fu_mon12,type="cat")



########### 8) Other Age Analyses ##########
mtn020AllAge <- read.csv("/Volumes/amishra$/HIV RS - Young Pop - predictors/adata/MTN_data/mtn020AllAge.csv")

## reformatting variables ##

mtn020AllAge$educ_relev <- ifelse(mtn020AllAge$educ=="post-secondary",1,0)

mtn020AllAge$drinkcat <- NULL
mtn020AllAge$drinkcat[mtn020AllAge$drink==0] <- 0
mtn020AllAge$drinkcat[mtn020AllAge$drink>0 & mtn020AllAge$drink <= 3] <- 1
mtn020AllAge$drinkcat[mtn020AllAge$drink>3 & mtn020AllAge$drink <= 6] <- 2
mtn020AllAge$drinkcat[mtn020AllAge$drink>6] <- 3
mtn020AllAge$drink_relv <- ifelse(mtn020AllAge$drinkcat==0,0,1)

mtn020AllAge$anysti <- ifelse(mtn020AllAge$trr==1 | mtn020AllAge$ct==1 | mtn020AllAge$gon==1 | mtn020AllAge$syph==1, 1,0)

mtn020AllAge$vagPH_cat <- NULL
mtn020AllAge$vagPH_cat[mtn020AllAge$vagPH< 4.5] <- 1
mtn020AllAge$vagPH_cat[mtn020AllAge$vagPH>= 4.5 & mtn020AllAge$vagPH < 5.5] <- 2
mtn020AllAge$vagPH_cat[mtn020AllAge$vagPH>= 5.5] <- 3

mtn020AllAge$revContraNo <- ifelse(mtn020AllAge$revContra==0,1,0)

mtn020AllAge$numpart_relev <- ifelse(mtn020AllAge$numPartcat>=2,1,0)

mtn020AllAge$region1 <- relevel(mtn020AllAge$region1, ref="Durban")

mtn020AllAge$region2 <- relevel(mtn020AllAge$region2, ref="North Durban")

mtn020AllAge$pregcat <- NULL
mtn020AllAge$pregcat[mtn020AllAge$BFPprgct==0] <- 0
mtn020AllAge$pregcat[mtn020AllAge$BFPprgct>0 ] <- 1

mtn020AllAge$liveBirth_relev[mtn020AllAge$liveBirth==0] <- 0
mtn020AllAge$liveBirth_relev[mtn020AllAge$liveBirth>0] <- 1

mtn020AllAge$BBAnopsp <- as.numeric(mtn020AllAge$BBAnopsp)
mtn020AllAge$numPartcat <- NULL
mtn020AllAge$numPartcat[mtn020AllAge$BBAnopsp==0] <- 1
mtn020AllAge$numPartcat[mtn020AllAge$BBAnopsp>0 & mtn020AllAge$BBAnopsp <= 1] <- 2
mtn020AllAge$numPartcat[mtn020AllAge$BBAnopsp>1] <- 3

mtn020AllAge$BBA3mvsn[mtn020AllAge$BBA3mvsn==99] <- NA

mtn020AllAge$analSex_cat <- NULL
mtn020AllAge$analSex_cat[mtn020AllAge$BBA3masn==0] <- 1
mtn020AllAge$analSex_cat[mtn020AllAge$BBA3masn!=0] <- 2
mtn020AllAge$BBA3masn <- factor(mtn020AllAge$BBA3masn)

mtn020AllAge$region1 <- relevel(mtn020AllAge$region1,ref="Durban")
levels(mtn020AllAge$region2) <-  c("North Durban","West Durban","South Durban","Cape Town","Johannesburg")

#### 22-26 #######
mtn020age2 <- mtn020AllAge[mtn020AllAge$age >= 22 & mtn020AllAge$age <= 26,]

nrow(mtn020age2)
length(unique(mtn020age2$ptid))
table(mtn020age2$censor)
table(mtn020age2$censor12mon)

#overall
nevent <- table(mtn020age2$censor)[2]
py <- sum(mtn020age2$fu_mos)/12
ir.allFU <- round(100*epi.conf(cbind(nevent,py),ctype = "inc.rate",method="exact"),2)

#covariates
age.allFU <- ir.fun(x=mtn020age2$age,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="number")
age.cat.allFU <- ir.fun(x=mtn020age2$age,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

educ.allFU <- ir.fun(x=mtn020age2$educ_relev,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

alc.allFU <- ir.fun(x=mtn020age2$drink_relv,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

income.allFU <- ir.fun(x=mtn020age2$demmerrd,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

tv.allFU <- ir.fun(x=mtn020age2$trr,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
ct.allFU <- ir.fun(x=mtn020age2$ct,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
gon.allFU <- ir.fun(x=mtn020age2$gon,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
syph.allFU <- ir.fun(x=mtn020age2$syph,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
any.allFU <- ir.fun(x=mtn020age2$anysti,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

ph.allFU <- ir.fun(x=mtn020age2$vagPH_cat,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

contra.allFU <- ir.fun(x=mtn020age2$revContraNo,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
preg.allFU <- ir.fun(x=mtn020age2$pregcat,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
live.allFU <- ir.fun(x=mtn020age2$liveBirth_relev,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
parthiv.allFU <- ir.fun(x=mtn020age2$parthiv,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
partcirc.allFU <- ir.fun(x=mtn020age2$partcirc,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
part3mon.allFU <- ir.fun(x=mtn020age2$part3mon,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

numpart.allFU <- ir.fun(x=mtn020age2$numpart_relev,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
vagacts.allFU <- ir.fun(x=mtn020age2$BBA3mvsn,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="number")
analacts.allFU <- ir.fun(x=mtn020age2$analSex_cat,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
ringworry.allFU <- ir.fun(x=mtn020age2$ringworry,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
vagwash.allFU <- ir.fun(x=mtn020age2$vagwash,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
transsex.allFU <- ir.fun(x=mtn020age2$qexch,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
otherPart.allFU <- ir.fun(x=mtn020age2$otherpart,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")


hivworry.allFU <- ir.fun(x=mtn020age2$hivworry,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")
region1.allFU <- ir.fun(x=mtn020age2$region1,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

region2.allFU <- ir.fun(x=mtn020age2$region2,event=mtn020age2$censor,fuTime=mtn020age2$fu_mos,type="cat")

tabAge2 <- rbind(age.allFU,age.cat.allFU,educ.allFU,alc.allFU,income.allFU,
               tv.allFU,ct.allFU,gon.allFU,syph.allFU,any.allFU,
               ph.allFU,contra.allFU,preg.allFU,live.allFU,
               parthiv.allFU,
               partcirc.allFU,part3mon.allFU,numpart.allFU,vagacts.allFU,analacts.allFU,ringworry.allFU,
               vagwash.allFU,transsex.allFU,otherPart.allFU,hivworry.allFU,region1.allFU,region2.allFU)
write.csv(tabAge2,"/Users/anumishra/Dropbox/RA/HIV predictors - Young SA pop/raw results/tabAge2.csv")

## PAF table for abstract ##

#STI 
table(mtn020age2$anysti)
round(100*prop.table(table(mtn020age2$anysti)),3)
par.fun(mtn020age2$anysti,censor = mtn020age2$censor,futime=mtn020age2$fu_mos,ir.all = ir.allFU[1], levels = 1)

#live births
table(mtn020age2$liveBirth_relev)
round(100*prop.table(table(mtn020age2$liveBirth_relev)),3)
par.fun(mtn020age2$liveBirth_relev,censor = mtn020age2$censor,futime=mtn020age2$fu_mos,ir.all = ir.allFU[1], levels = 1)

#other partners
table(mtn020age2$otherpart)
round(100*prop.table(table(mtn020age2$otherpart)),3)
par.fun(x = mtn020age2$otherpart,censor = mtn020age2$censor,futime=mtn020age2$fu_mos,ir.all = ir.allFU[1], levels = 2)

#>1 sex partner
table(mtn020age2$numpart_relev)
round(100*prop.table(table(mtn020age2$numpart_relev)),3)
par.fun(x = mtn020age2$numpart_relev,censor = mtn020age2$censor,futime=mtn020age2$fu_mos,ir.all = ir.allFU[1], levels = 1)

#### 27 + #######
mtn020age3 <- mtn020AllAge[mtn020AllAge$age >= 27,]

nrow(mtn020age3)
length(unique(mtn020age3$ptid))
table(mtn020age3$censor)
table(mtn020age3$censor12mon)


#overall
nevent <- table(mtn020age3$censor)[2]
py <- sum(mtn020age3$fu_mos)/12
ir.allFU <- round(100*epi.conf(cbind(nevent,py),ctype = "inc.rate",method="exact"),2)

#covariates
age.allFU <- ir.fun(x=mtn020age3$age,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="number")
mtn020age3$age_cat <- NULL
mtn020age3$age_cat[mtn020age3$age>=27 & mtn020age3$age<=30] <- 1
mtn020age3$age_cat[mtn020age3$age>=31 & mtn020age3$age<=35] <- 2
mtn020age3$age_cat[mtn020age3$age>=36] <- 3
age.cat.allFU <- ir.fun(x=mtn020age3$age_cat,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

educ.allFU <- ir.fun(x=mtn020age3$educ_relev,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

alc.allFU <- ir.fun(x=mtn020age3$drink_relv,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

income.allFU <- ir.fun(x=mtn020age3$demmerrd,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

tv.allFU <- ir.fun(x=mtn020age3$trr,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
ct.allFU <- ir.fun(x=mtn020age3$ct,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
gon.allFU <- ir.fun(x=mtn020age3$gon,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
syph.allFU <- ir.fun(x=mtn020age3$syph,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
any.allFU <- ir.fun(x=mtn020age3$anysti,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

ph.allFU <- ir.fun(x=mtn020age3$vagPH_cat,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

contra.allFU <- ir.fun(x=mtn020age3$revContraNo,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
preg.allFU <- ir.fun(x=mtn020age3$pregcat,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
live.allFU <- ir.fun(x=mtn020age3$liveBirth_relev,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
parthiv.allFU <- ir.fun(x=mtn020age3$parthiv,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
partcirc.allFU <- ir.fun(x=mtn020age3$partcirc,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
part3mon.allFU <- ir.fun(x=mtn020age3$part3mon,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

numpart.allFU <- ir.fun(x=mtn020age3$numpart_relev,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
vagacts.allFU <- ir.fun(x=mtn020age3$BBA3mvsn,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="number")
analacts.allFU <- ir.fun(x=mtn020age3$analSex_cat,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
ringworry.allFU <- ir.fun(x=mtn020age3$ringworry,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
vagwash.allFU <- ir.fun(x=mtn020age3$vagwash,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
transsex.allFU <- ir.fun(x=mtn020age3$qexch,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
otherPart.allFU <- ir.fun(x=mtn020age3$otherpart,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")


hivworry.allFU <- ir.fun(x=mtn020age3$hivworry,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")
region1.allFU <- ir.fun(x=mtn020age3$region1,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

region2.allFU <- ir.fun(x=mtn020age3$region2,event=mtn020age3$censor,fuTime=mtn020age3$fu_mos,type="cat")

tabAge3 <- rbind(age.allFU,age.cat.allFU,educ.allFU,alc.allFU,income.allFU,
                 tv.allFU,ct.allFU,gon.allFU,syph.allFU,any.allFU,
                 ph.allFU,contra.allFU,preg.allFU,live.allFU,
                 parthiv.allFU,
                 partcirc.allFU,part3mon.allFU,numpart.allFU,vagacts.allFU,analacts.allFU,ringworry.allFU,
                 vagwash.allFU,otherPart.allFU,hivworry.allFU,region1.allFU,region2.allFU)
write.csv(tabAge3,"/Users/anumishra/Dropbox/RA/HIV predictors - Young SA pop/raw results/tabAge3.csv")



## PAF table for abstract ##

#STI 
table(mtn020age3$anysti)
round(100*prop.table(table(mtn020age3$anysti)),3)
par.fun(mtn020age3$anysti,censor = mtn020age3$censor,futime=mtn020age3$fu_mos,ir.all = ir.allFU[1], levels = 1)

#live births
table(mtn020age3$liveBirth_relev)
round(100*prop.table(table(mtn020age3$liveBirth_relev)),3)
par.fun(x = mtn020age3$liveBirth_relev,censor = mtn020age3$censor,futime=mtn020age3$fu_mos,ir.all = ir.allFU[1], levels = 1)

#other partners
table(mtn020age3$otherpart)
round(100*prop.table(table(mtn020age3$otherpart)),3)
par.fun(x = mtn020age3$otherpart,censor = mtn020age3$censor,futime=mtn020age3$fu_mos,ir.all = ir.allFU[1], levels = 2)

#>1 sex partner
table(mtn020age3$numpart_relev)
round(100*prop.table(table(mtn020age3$numpart_relev)),3)
par.fun(x = mtn020age3$numpart_relev,censor = mtn020age3$censor,futime=mtn020age3$fu_mos,ir.all = ir.allFU[1], levels = 1)

