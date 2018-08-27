
####Set working directory under your home directory ('R'is my home directory)
setwd("~/R")

##Plot the results from ecoPCR
##tutorial http://metabarcoding.org/IMG/html/primerdesign.html


library(devtools)

#install rtool first https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows


# Sys.getenv('PATH')#make sure add the path to Rtool into thr system PATH

#download ROBITools, ROBITaxonomy, ROBIBarcodes packages from https://git.metabarcoding.org/obitools save to src directory
#install the packages use install.packages() as below

#install.packages('../src/ROBITaxonomy-master.tar.gz', repos = NULL, type="source")

library(ROBITaxonomy)


#install.packages('../src/ROBITools-master.tar.gz', repos = NULL, type="source")

library(ROBITools)


#install.packages('../src/ROBIBarcodes-master.tar.gz', repos = NULL, type="source")


library(ROBIBarcodes)


taxo = read.taxonomy('ecoPCR/ncbi20180701/ncbi20180701')


#The four files ncbi20180701.adx, ncbi20180701.ndx, ncbi20180701.rdx, ncbi20180701.tdx constitue the reformated taxonomy.
#need to prepare in Linux see 'Formating the taxonomy for OBITools' in http://metabarcoding.org/IMG/html/primerdesign.html


###12S_01_Tele02####

#read the ecoPCR results
Tele02_12S <- read.ecopcr.result('ecoPCR/Tele02-12S-mismatch_3/EcoPCR_report.txt')



#EXculde the title infromation from the EcoPCR_report.txt

Tele02_12S_sub <- Tele02_12S [17:nrow(Tele02_12S),]


head(Tele02_12S_sub,n=2)


#Testing the conservation of the priming sites
Tele02_12S_forward = ecopcr.forward.shanon(ecopcr = Tele02_12S_sub)

Tele02_12S_reverse = ecopcr.reverse.shanon(ecopcr = Tele02_12S_sub)

par(mfcol=c(2,1))


dnalogoplot(Tele02_12S_forward,
            primer = "AAACTCGTGCCAGCCACC",
            main='Tele02_F',xlab="Forward primer",ylab='bits')

dnalogoplot(Tele02_12S_reverse,
            primer = "GGGTATCTAATCCCAGTTTG",
            main='Tele02_R',xlab="Reverse primer",ylab='bits')


#How mismatches influence taxonomical selection in primers


mismatchplot(Tele02_12S_sub)+title(xlab="Number of mismatches with forward primer", 
                                   ylab="Number of mismatches with reverse primer",
                                   main = 'Distribution of the number of mismatches with primer')


#check the mismatch in primers is not 0
Tele02_12S_mis <- Tele02_12S_sub [(Tele02_12S_sub$forward_mismatch > 1 |
                                     Tele02_12S_sub$reverse_mismatch > 1 ),]

#What is the taxonomic resolution of our marker


res = resolution(taxo,Tele02_12S_sub)

resolution = with(Tele02_12S_sub,
                  unique(data.frame(species_name,taxid,rank=res)))

t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


round (table(res)/sum(table(res))*100) ##sequences resolution


#check the resolution results with species level


speciesname <- unique(Tele02_12S_sub[res=='species',]$species_name)
attributes(speciesname)

#check the different between two species list to check which species can not be distinguished 
setdiff(levels(speciesname),levels(droplevels(speciesname)))



unique(Tele02_12S_sub[res=='genus',]$genus_name)

levels(unique(Tele02_12S_sub[res=='genus',]$genus_name))

unique(Tele02_12S_sub[res=='family',]$family_name)

resolution[which(resolution$rank =='genus'),]

resolution[which(resolution$rank =='family'),]


#amplicom distribution
min(Tele02_12S_sub$amplicon_length)
max(Tele02_12S_sub$amplicon_length)
mean(Tele02_12S_sub$amplicon_length)



# Make plot


png("./ecoPCR/Figures/Tele02_12S.png", width = 7, height = 8, units = 'in', res = 500, pointsize = 12)




layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE),
       heights=c(2,2,3)) #arrange the plots 

dnalogoplot(Tele02_12S_forward,
            primer = "AAACTCGTGCCAGCCACC",
            main='Tele02_F',xlab="Forward primer",ylab='bits')
mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1) 

dnalogoplot(Tele02_12S_reverse,
            primer = "GGGTATCTAATCCCAGTTTG",
            main='Tele02_R',xlab="Reverse primer",ylab='bits')
mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1) 


mismatchplot(Tele02_12S_sub,col='grey')+title(xlab="Number of mismatches with forward primer", 
                                              ylab="Number of mismatches with reverse primer")
mtext('(c)', side = 3, line = 1, adj = 0.0, cex = 1) 



hist(Tele02_12S_sub$amplicon_length,xlab ="Metabarcode length in bp", ylab= "Frequency",
     main=NULL, xlim = range(172),xaxt="n",col = 'grey')
axis(1, at=seq(10,210,by=10), labels=seq(10,210,by=10))

mtext('(d)', side = 3, line = 1, adj = 0.0, cex = 1) 

dev.off()




###12S_02_Mifish####

#read the ecoPCR results
Mifisha_12S <- read.ecopcr.result('ecoPCR/Mifish_adpt-12S-mismatch_3/EcoPCR_report.txt')





#EXculde the title infromation from the EcoPCR_report.txt

Mifisha_12S_sub <- Mifisha_12S [17:nrow(Mifisha_12S),]


head(Mifisha_12S_sub,n=2)


#Testing the conservation of the priming sites
Mifisha_12S_forward = ecopcr.forward.shanon(ecopcr = Mifisha_12S_sub)

Mifisha_12S_reverse = ecopcr.reverse.shanon(ecopcr = Mifisha_12S_sub)

par(mfcol=c(2,1))


dnalogoplot(Mifisha_12S_forward,
            primer = "GCCGGTAAAACTCGTGCCAGC",
            main='MiFish-E-Fa',xlab="Forward primer",ylab='bits')

dnalogoplot(Mifisha_12S_reverse,
            primer = "CATAGTGGGGTATCTAATCCCAGTTTG",
            main='MiFish-E-R',xlab="Reverse primer",ylab='bits')


#How mismatches influence taxonomical selection in primers


mismatchplot(Mifisha_12S_sub)+title(xlab="Number of mismatches with forward primer", 
                                  ylab="Number of mismatches with reverse primer",
                                  main = 'Distribution of the number of mismatches with primer')


#check the mismatch in primers is not 0
Mifisha_12S_mis <- Mifisha_12S_sub [(Mifisha_12S_sub$forward_mismatch > 1 |
                                   Mifisha_12S_sub$reverse_mismatch > 1 ),]

#What is the taxonomic resolution of our marker


res = resolution(taxo,Mifisha_12S_sub)

resolution = with(Mifisha_12S_sub,
                  unique(data.frame(species_name,taxid,rank=res)))

t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


round (table(res)/sum(table(res))*100) ##sequences resolution


#check the resolution results with species level


speciesname <- unique(Mifisha_12S_sub[res=='species',]$species_name)
attributes(speciesname)

#check the different between two species list to check which species can not be distinguished 
setdiff(levels(speciesname),levels(droplevels(speciesname)))



unique(Mifisha_12S_sub[res=='genus',]$genus_name)

levels(unique(Mifisha_12S_sub[res=='genus',]$genus_name))

unique(Mifisha_12S_sub[res=='family',]$family_name)

resolution[which(resolution$rank =='genus'),]

resolution[which(resolution$rank =='family'),]


#amplicom distribution
min(Mifisha_12S_sub$amplicon_length)
max(Mifisha_12S_sub$amplicon_length)
mean(Mifisha_12S_sub$amplicon_length)



# Make plot


png("./ecoPCR/Figures/Mifisha_12S.png", width = 7, height = 8, units = 'in', res = 500, pointsize = 12)


layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE),
       heights=c(2,2,3)) #arrange the plots 

dnalogoplot(Mifisha_12S_forward,
            primer = "GCCGGTAAAACTCGTGCCAGC",
            main='MiFish-E-Fa',xlab="Forward primer",ylab='bits')
mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1) 

dnalogoplot(Mifisha_12S_reverse,
            primer = "CATAGTGGGGTATCTAATCCCAGTTTG",
            main='MiFish-E-R',xlab="Reverse primer",ylab='bits')

mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1) 


mismatchplot(Mifisha_12S_sub,col='grey')+title(xlab="Number of mismatches with forward primer", 
                                             ylab="Number of mismatches with reverse primer")
mtext('(c)', side = 3, line = 1, adj = 0.0, cex = 1) 



hist(Mifisha_12S_sub$amplicon_length,xlab ="Metabarcode length in bp", ylab= "Frequency",
     main=NULL, xlim = range(172),xaxt="n",col = 'grey')
axis(1, at=seq(10,210,by=10), labels=seq(10,210,by=10))

mtext('(d)', side = 3, line = 1, adj = 0.0, cex = 1) 

dev.off()


###12S_03_Kelly####

#read the ecoPCR results
Kelly_12S <- read.ecopcr.result('ecoPCR/kelly-12S-mismatch_3/EcoPCR_report.txt')




#EXculde the title infromation from the EcoPCR_report.txt

Kelly_12S_sub <- Kelly_12S [17:nrow(Kelly_12S),]


head(Kelly_12S_sub,n=2)


#Testing the conservation of the priming sites
Kelly_12S_forward = ecopcr.forward.shanon(ecopcr = Kelly_12S_sub)

Kelly_12S_reverse = ecopcr.reverse.shanon(ecopcr = Kelly_12S_sub)

par(mfcol=c(2,1))


dnalogoplot(Kelly_12S_forward,
            primer = "ACTGGGATTAGATACCCC",
            main='12S_V5_F',xlab="Forward primers",ylab='bits')

dnalogoplot(Kelly_12S_reverse,
            primer = "TAGAACAGGCTCCTCTAG",
            main='12S_V5_R',xlab="Reverse primers",ylab='bits')


#How mismatches influence taxonomical selection in primers


mismatchplot(Kelly_12S_sub)+title(xlab="Number of mismatches with forward primer", 
                                  ylab="Number of mismatches with reverse primer",
                                  main = 'Distribution of the number of mismatches with primer')


#check the mismatch in primers is not 0
Kelly_12S_mis <- Kelly_12S_sub [(Kelly_12S_sub$forward_mismatch > 1 |
                                   Kelly_12S_sub$reverse_mismatch > 1 ),]

#What is the taxonomic resolution of our marker


res = resolution(taxo,Kelly_12S_sub)

resolution = with(Kelly_12S_sub,
                  unique(data.frame(species_name,taxid,rank=res)))

t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


round (table(res)/sum(table(res))*100) ##sequences resolution


#check the resolution results with species level


speciesname <- unique(Kelly_12S_sub[res=='species',]$species_name)
attributes(speciesname)

#check the different between two species list to check which species can not be distinguished 
setdiff(levels(speciesname),levels(droplevels(speciesname)))



unique(Kelly_12S_sub[res=='genus',]$genus_name)

levels(unique(Kelly_12S_sub[res=='genus',]$genus_name))

unique(Kelly_12S_sub[res=='family',]$family_name)

resolution[which(resolution$rank =='genus'),]

resolution[which(resolution$rank =='family'),]


#amplicom distribution
min(Kelly_12S_sub$amplicon_length)
max(Kelly_12S_sub$amplicon_length)
mean(Kelly_12S_sub$amplicon_length)




# Make plot

png("./ecoPCR/Figures/Kelly_12S.png", width = 7, height = 8, units = 'in', res = 500, pointsize = 12)



layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE),
       heights=c(2,2,3)) #arrange the plots 

dnalogoplot(Kelly_12S_forward,
            primer = "ACTGGGATTAGATACCCC",
            main='12S_V5_F', xlab="Forward primers",ylab='bits')
mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1) 

dnalogoplot(Kelly_12S_reverse,
            primer = "TAGAACAGGCTCCTCTAG",
            main='12S_V5_R', xlab="Reverse primers",ylab='bits')

mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1) 


mismatchplot(Kelly_12S_sub,col='grey')+title(xlab="Number of mismatches with forward primer", 
                                             ylab="Number of mismatches with reverse primer")
mtext('(c)', side = 3, line = 1, adj = 0.0, cex = 1) 



hist(Kelly_12S_sub$amplicon_length,xlab ="Metabarcode length in bp", ylab= "Frequency",
     main=NULL, xlim = range(80),xaxt="n",col = 'grey')
axis(1, at=seq(10,120,by=2), labels=seq(10,120,by=2))

mtext('(d)', side = 3, line = 1, adj = 0.0, cex = 1) 

dev.off()



# Kelly_length<-ggplot(Kelly_12S_sub,aes(amplicon_length))+
#   geom_bar()+
#   scale_y_continuous(breaks=seq(0, 150, 20))+
#   scale_x_continuous(limits = c(80, 120),breaks=seq(0, 150, 10))+
#   labs(x="Metabarcode length in bp", y= "Frequency",title="Amplicom length distribution")+
#   theme_bw()+theme(text=element_text(size=15),
#                    panel.grid.minor=element_blank(),panel.grid.major =element_blank(), 
#                    plot.title = element_text(size = rel(0.8),face = "bold",hjust=0.5))




###12S_04_kellyR2####

#read the ecoPCR results
KellyR2_12S <- read.ecopcr.result('ecoPCR/V5_R2-12S-mismatch_3/EcoPCR_report.txt')





#EXculde the title infromation from the EcoPCR_report.txt

KellyR2_12S_sub <- KellyR2_12S [17:nrow(KellyR2_12S),]


head(KellyR2_12S_sub,n=2)


#Testing the conservation of the priming sites
KellyR2_12S_forward = ecopcr.forward.shanon(ecopcr = KellyR2_12S_sub)

KellyR2_12S_reverse = ecopcr.reverse.shanon(ecopcr = KellyR2_12S_sub)

par(mfcol=c(2,1))


dnalogoplot(KellyR2_12S_forward,
            primer = "AAACTCGTGCCAGCCACC",
            main='12S_V5_F',xlab="Forward primer",ylab='bits')

dnalogoplot(KellyR2_12S_reverse,
            primer = "CTACACCTCGACCTGACG",
            main='12S_V5_R2',xlab="Reverse primer",ylab='bits')


#How mismatches influence taxonomical selection in primers


mismatchplot(KellyR2_12S_sub)+title(xlab="Number of mismatches with forward primer", 
                                   ylab="Number of mismatches with reverse primer",
                                   main = 'Distribution of the number of mismatches with primer')


#check the mismatch in primers is not 0
KellyR2_12S_mis <- KellyR2_12S_sub [(KellyR2_12S_sub$forward_mismatch > 1 |
                                     KellyR2_12S_sub$reverse_mismatch > 1 ),]

#What is the taxonomic resolution of our marker


res = resolution(taxo,KellyR2_12S_sub)

resolution = with(KellyR2_12S_sub,
                  unique(data.frame(species_name,taxid,rank=res)))

t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


round (table(res)/sum(table(res))*100) ##sequences resolution


#check the resolution results with species level


speciesname <- unique(KellyR2_12S_sub[res=='species',]$species_name)
attributes(speciesname)

#check the different between two species list to check which species can not be distinguished 
setdiff(levels(speciesname),levels(droplevels(speciesname)))

levels(unique(KellyR2_12S_sub[res=='species',]$species_name))

unique(KellyR2_12S_sub[res=='genus',]$genus_name)

levels(unique(KellyR2_12S_sub[res=='genus',]$genus_name))

unique(KellyR2_12S_sub[res=='family',]$family_name)

resolution[which(resolution$rank =='genus'),]

resolution[which(resolution$rank =='family'),]


#amplicom distribution
min(KellyR2_12S_sub$amplicon_length)
max(KellyR2_12S_sub$amplicon_length)
mean(KellyR2_12S_sub$amplicon_length)



# Make plot


png("./ecoPCR/Figures/KellyR2_12S.png", width = 7, height = 8, units = 'in', res = 500, pointsize = 12)




layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE),
       heights=c(2,2,3)) #arrange the plots 


dnalogoplot(KellyR2_12S_forward,
            primer = "AAACTCGTGCCAGCCACC",
            main='12S_V5_F',xlab="Forward primer",ylab='bits')
mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1) 

dnalogoplot(KellyR2_12S_reverse,
            primer = "CTACACCTCGACCTGACG",
            main='12S_V5_R2',xlab="Reverse primer",ylab='bits')
mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1) 


mismatchplot(KellyR2_12S_sub,col='grey')+title(xlab="Number of mismatches with forward primer", 
                                              ylab="Number of mismatches with reverse primer")
mtext('(c)', side = 3, line = 1, adj = 0.0, cex = 1) 



hist(KellyR2_12S_sub$amplicon_length,xlab ="Metabarcode length in bp", ylab= "Frequency",
     main=NULL, xlim = range(200),xaxt="n",col = 'grey')
axis(1, at=seq(10,250,by=10), labels=seq(10,250,by=10))

mtext('(d)', side = 3, line = 1, adj = 0.0, cex = 1) 

dev.off()


###12S_05_kellyR3####

#read the ecoPCR results
KellyR3_12S <- read.ecopcr.result('ecoPCR/V5_R3-12S-mismatch_3/EcoPCR_report.txt')





#EXculde the title infromation from the EcoPCR_report.txt

KellyR3_12S_sub <- KellyR3_12S [17:nrow(KellyR3_12S),]


head(KellyR3_12S_sub,n=2)


#Testing the conservation of the priming sites
KellyR3_12S_forward = ecopcr.forward.shanon(ecopcr = KellyR3_12S_sub)

KellyR3_12S_reverse = ecopcr.reverse.shanon(ecopcr = KellyR3_12S_sub)

par(mfcol=c(2,1))


dnalogoplot(KellyR3_12S_forward,
            primer = "AAACTCGTGCCAGCCACC",
            main='12S_V5_F',xlab="Forward primer",ylab='bits')

dnalogoplot(KellyR3_12S_reverse,
            primer = "GAGAGTGACGGGCGGTGT",
            main='12S_V5_R3',xlab="Reverse primer",ylab='bits')


#How mismatches influence taxonomical selection in primers


mismatchplot(KellyR3_12S_sub)+title(xlab="Number of mismatches with forward primer", 
                                    ylab="Number of mismatches with reverse primer",
                                    main = 'Distribution of the number of mismatches with primer')


#check the mismatch in primers is not 0
KellyR3_12S_mis <- KellyR3_12S_sub [(KellyR3_12S_sub$forward_mismatch > 1 |
                                       KellyR3_12S_sub$reverse_mismatch > 1 ),]

#What is the taxonomic resolution of our marker


res = resolution(taxo,KellyR3_12S_sub)

resolution = with(KellyR3_12S_sub,
                  unique(data.frame(species_name,taxid,rank=res)))

t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


round (table(res)/sum(table(res))*100) ##sequences resolution


#check the resolution results with species level


speciesname <- unique(KellyR3_12S_sub[res=='species',]$species_name)
attributes(speciesname)

#check the different between two species list to check which species can not be distinguished 
setdiff(levels(speciesname),levels(droplevels(speciesname)))

levels(unique(KellyR3_12S_sub[res=='species',]$species_name))

unique(KellyR3_12S_sub[res=='genus',]$genus_name)

levels(unique(KellyR3_12S_sub[res=='genus',]$genus_name))

unique(KellyR3_12S_sub[res=='family',]$family_name)

resolution[which(resolution$rank =='genus'),]

resolution[which(resolution$rank =='family'),]


#amplicom distribution
min(KellyR3_12S_sub$amplicon_length)
max(KellyR3_12S_sub$amplicon_length)
mean(KellyR3_12S_sub$amplicon_length)



# Make plot


png("./ecoPCR/Figures/KellyR3_12S.png", width = 7, height = 8, units = 'in', res = 500, pointsize = 12)




layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE),
       heights=c(2,2,3)) #arrange the plots 


dnalogoplot(KellyR3_12S_forward,
            primer = "AAACTCGTGCCAGCCACC",
            main='12S_V5_F',xlab="Forward primer",ylab='bits')
mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1) 

dnalogoplot(KellyR3_12S_reverse,
            primer = "GAGAGTGACGGGCGGTGT",
            main='Ac12s_R',xlab="Reverse primer",ylab='bits')

mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1) 


mismatchplot(KellyR3_12S_sub,col='grey')+title(xlab="Number of mismatches with forward primer", 
                                               ylab="Number of mismatches with reverse primer")
mtext('(c)', side = 3, line = 1, adj = 0.0, cex = 1) 



hist(KellyR3_12S_sub$amplicon_length,xlab ="Metabarcode length in bp", ylab= "Frequency",
     main=NULL, xlim = range(350),xaxt="n",col = 'grey')
axis(1, at=seq(100,400,by=20), labels=seq(100,400,by=20))

mtext('(d)', side = 3, line = 1, adj = 0.0, cex = 1) 

dev.off()


###12S_06_Teleo####

#read the ecoPCR results
Teleo_12S <- read.ecopcr.result('ecoPCR/teleo-12S-mismatch_3/EcoPCR_report.txt')




#EXculde the title infromation from the EcoPCR_report.txt

Teleo_12S_sub <- Teleo_12S [17:nrow(Teleo_12S),]


head(Teleo_12S_sub,n=2)


#Testing the conservation of the priming sites
Teleo_12S_forward = ecopcr.forward.shanon(ecopcr = Teleo_12S_sub)

Teleo_12S_reverse = ecopcr.reverse.shanon(ecopcr = Teleo_12S_sub)

par(mfcol=c(2,1))


dnalogoplot(Teleo_12S_forward,
            primer = "ACACCGCCCGTCACTCT",
            main='Teleo_F',xlab="Forward primers",ylab='bits')

dnalogoplot(Teleo_12S_reverse,
            primer = "CTTCCGGTACACTTACCATG",
            main='Teleo_R',xlab="Reverse primers",ylab='bits')


#How mismatches influence taxonomical selection in primers


mismatchplot(Teleo_12S_sub)+title(xlab="Number of mismatches with forward primer", 
                                  ylab="Number of mismatches with reverse primer",
                                  main = 'Distribution of the number of mismatches with primer')


#check the mismatch in primers is not 0
Teleo_12S_mis <- Teleo_12S_sub [(Teleo_12S_sub$forward_mismatch > 1 |
                                   Teleo_12S_sub$reverse_mismatch > 1 ),]

#What is the taxonomic resolution of our marker


res = resolution(taxo,Teleo_12S_sub)

resolution = with(Teleo_12S_sub,
                  unique(data.frame(species_name,taxid,rank=res)))

t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


round (table(res)/sum(table(res))*100) ##sequences resolution


#check the resolution results with species level


speciesname <- unique(Teleo_12S_sub[res=='species',]$species_name)
attributes(speciesname)

#check the different between two species list to check which species can not be distinguished 
setdiff(levels(speciesname),levels(droplevels(speciesname)))



unique(Teleo_12S_sub[res=='genus',]$genus_name)

levels(unique(Teleo_12S_sub[res=='genus',]$genus_name))

unique(Teleo_12S_sub[res=='family',]$family_name)

resolution[which(resolution$rank =='genus'),]

resolution[which(resolution$rank =='family'),]


#amplicom distribution
min(Teleo_12S_sub$amplicon_length)
max(Teleo_12S_sub$amplicon_length)
mean(Teleo_12S_sub$amplicon_length)






# Make plot

png("./ecoPCR/Figures/Teleo_12S.png", width = 7, height = 8, units = 'in', res = 500, pointsize = 12)



layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE),
       heights=c(2,2,3)) #arrange the plots 

dnalogoplot(Teleo_12S_forward,
            primer = "ACACCGCCCGTCACTCT",
            main='Teleo_F',xlab="Forward primers",ylab='bits')

mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1) 

dnalogoplot(Teleo_12S_reverse,
            primer = "CTTCCGGTACACTTACCATG",
            main='Teleo_R',xlab="Reverse primers",ylab='bits')

mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1) 


mismatchplot(Teleo_12S_sub,col='grey')+title(xlab="Number of mismatches with forward primer", 
                                             ylab="Number of mismatches with reverse primer")
mtext('(c)', side = 3, line = 1, adj = 0.0, cex = 1) 



hist(Teleo_12S_sub$amplicon_length,xlab ="Metabarcode length in bp", ylab= "Frequency",
     main=NULL, xlim = range(80),xaxt="n",col = 'grey')
axis(1, at=seq(10,100,by=2), labels=seq(10,100,by=2))

mtext('(d)', side = 3, line = 1, adj = 0.0, cex = 1) 

dev.off()



###Cytb02####

#read the ecoPCR results
Cytb02 <- read.ecopcr.result('ecoPCR/Fish2CBL-Cytb-mismatch_3/EcoPCR_report.txt')




#EXculde the title infromation from the EcoPCR_report.txt

Cytb02_sub <- Cytb02 [17:nrow(Cytb02),]


head(Cytb02_sub,n=2)


#Testing the conservation of the priming sites
Cytb02_forward = ecopcr.forward.shanon(ecopcr = Cytb02_sub)

Cytb02_reverse = ecopcr.reverse.shanon(ecopcr = Cytb02_sub)

par(mfcol=c(2,1))


dnalogoplot(Cytb02_forward,
            primer = "ACAACTTCACCCCTGCAAAC",
            main='Fish2CBL',xlab="Forward primers",ylab='bits')

dnalogoplot(Cytb02_reverse,
            primer = "GATGGCGTAGGCAAACAAGA",
            main='Fish2bCBR',xlab="Reverse primers",ylab='bits')


#How mismatches influence taxonomical selection in primers


mismatchplot(Cytb02_sub)+title(xlab="Number of mismatches with forward primer", 
                                  ylab="Number of mismatches with reverse primer",
                                  main = 'Distribution of the number of mismatches with primer')


#check the mismatch in primers is not 0
Cytb02_mis <- Cytb02_sub [(Cytb02_sub$forward_mismatch > 1 |
                                   Cytb02_sub$reverse_mismatch > 1 ),]

#What is the taxonomic resolution of our marker


res = resolution(taxo,Cytb02_sub)

resolution = with(Cytb02_sub,
                  unique(data.frame(species_name,taxid,rank=res)))

t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


round (table(res)/sum(table(res))*100) ##sequences resolution


#check the resolution results with species level


speciesname <- unique(Cytb02_sub[res=='species',]$species_name)
attributes(speciesname)

#check the different between two species list to check which species can not be distinguished 
setdiff(levels(speciesname),levels(droplevels(speciesname)))



unique(Cytb02_sub[res=='genus',]$genus_name)

levels(unique(Cytb02_sub[res=='genus',]$genus_name))

unique(Cytb02_sub[res=='family',]$family_name)

resolution[which(resolution$rank =='genus'),]

resolution[which(resolution$rank =='family'),]


#amplicom distribution
min(Cytb02_sub$amplicon_length)
max(Cytb02_sub$amplicon_length)
mean(Cytb02_sub$amplicon_length)







# Make plot

png("./ecoPCR/Figures/Cytb02.png", width = 7, height = 8, units = 'in', res = 500, pointsize = 12)



layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE),
       heights=c(2,2,3)) #arrange the plots 

dnalogoplot(Cytb02_forward,
            primer = "ACAACTTCACCCCTGCAAAC",
            main='Fish2CBL',xlab="Forward primers",ylab='bits')

mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1) 

dnalogoplot(Cytb02_reverse,
            primer = "GATGGCGTAGGCAAACAAGA",
            main='Fish2bCBR',xlab="Reverse primers",ylab='bits')


mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1) 


mismatchplot(Cytb02_sub,col='grey')+title(xlab="Number of mismatches with forward primer", 
                                             ylab="Number of mismatches with reverse primer")
mtext('(c)', side = 3, line = 1, adj = 0.0, cex = 1) 


counts_Cytb02 <- table(Cytb02_sub$amplicon_length)

barplot(counts_Cytb02,xlab ="Metabarcode length in bp", ylab= "Frequency",
        main=NULL,col = 'grey')


mtext('(d)', side = 3, line = 1, adj = 0.0, cex = 1) 

dev.off()


###Cytb03####

#read the ecoPCR results
Cytb03 <- read.ecopcr.result('ecoPCR/Fish2degCBL-Cytb-mismatch_3/EcoPCR_report.txt')




#EXculde the title infromation from the EcoPCR_report.txt

Cytb03_sub <- Cytb03 [17:nrow(Cytb03),]


head(Cytb03_sub,n=2)


#Testing the conservation of the priming sites
Cytb03_forward = ecopcr.forward.shanon(ecopcr = Cytb03_sub)

Cytb03_reverse = ecopcr.reverse.shanon(ecopcr = Cytb03_sub)

par(mfcol=c(2,1))


dnalogoplot(Cytb03_forward,
            primer = "ACAACTTCACCCCTGCRAAY",
            main='Fish2degCBL',xlab="Forward primers",ylab='bits')

dnalogoplot(Cytb03_reverse,
            primer = "GATGGCGTAGGCAAATAGGA",
            main='Fish2CBR',xlab="Reverse primers",ylab='bits')


#How mismatches influence taxonomical selection in primers


mismatchplot(Cytb03_sub)+title(xlab="Number of mismatches with forward primer", 
                               ylab="Number of mismatches with reverse primer",
                               main = 'Distribution of the number of mismatches with primer')


#check the mismatch in primers is not 0
Cytb03_mis <- Cytb03_sub [(Cytb03_sub$forward_mismatch > 1 |
                             Cytb03_sub$reverse_mismatch > 1 ),]

#What is the taxonomic resolution of our marker


res = resolution(taxo,Cytb03_sub)

resolution = with(Cytb03_sub,
                  unique(data.frame(species_name,taxid,rank=res)))

t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


round (table(res)/sum(table(res))*100) ##sequences resolution


#check the resolution results with species level


speciesname <- unique(Cytb03_sub[res=='species',]$species_name)
attributes(speciesname)

#check the different between two species list to check which species can not be distinguished 
setdiff(levels(speciesname),levels(droplevels(speciesname)))



unique(Cytb03_sub[res=='genus',]$genus_name)

levels(unique(Cytb03_sub[res=='genus',]$genus_name))

unique(Cytb03_sub[res=='family',]$family_name)

resolution[which(resolution$rank =='genus'),]

resolution[which(resolution$rank =='family'),]


#amplicom distribution
min(Cytb03_sub$amplicon_length)
max(Cytb03_sub$amplicon_length)
mean(Cytb03_sub$amplicon_length)







# Make plot

png("./ecoPCR/Figures/Cytb03.png", width = 7, height = 8, units = 'in', res = 500, pointsize = 12)



layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE),
       heights=c(2,2,3)) #arrange the plots 

dnalogoplot(Cytb03_forward,
            primer = "ACAACTTCACCCCTGCRAAY",
            main='Fish2degCBL',xlab="Forward primers",ylab='bits')

mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1) 


dnalogoplot(Cytb03_reverse,
            primer = "GATGGCGTAGGCAAATAGGA",
            main='Fish2CBR',xlab="Reverse primers",ylab='bits')


mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1) 


mismatchplot(Cytb03_sub,col='grey')+title(xlab="Number of mismatches with forward primer", 
                                          ylab="Number of mismatches with reverse primer")
mtext('(c)', side = 3, line = 1, adj = 0.0, cex = 1) 


counts_Cytb03 <- table(Cytb03_sub$amplicon_length)

barplot(counts_Cytb03,xlab ="Metabarcode length in bp", ylab= "Frequency",
        main=NULL,col = 'grey')


mtext('(d)', side = 3, line = 1, adj = 0.0, cex = 1) 

dev.off()




