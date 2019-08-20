setwd("~/Stanford/writings/Thesis/finaldataforpaper/Ch2_dataprocessing/randomforests")

# predictions

rFsite_predictions_summary<-read.csv('rFsite_predictions_summary.csv', header=T)
rFsite_predictions_table<-table(rFsite_predictions_summary[,2:3])
rFsite_predictions_proportions<-t(rFsite_predictions_table/rowSums(rFsite_predictions_table))

rFSal_predictions_summary<-read.csv('rFSal_predictions_summary.csv', header=T)
rFSal_predictions_summary$color<-as.character(rFSal_predictions_summary$color)

rFNtot_predictions_summary<-read.csv('rFNtot_predictions_summary.csv', header=T)
rFNtot_predictions_summary$color<-as.character(rFNtot_predictions_summary$color)

rFTemp_predictions_summary<-read.csv('rFTemp_predictions_summary.csv', header=T)
rFTemp_predictions_summary$color<-as.character(rFTemp_predictions_summary$color)


#tiff(filename='~/lab/writings/Thesis/finaldataforpaper/Ch2_dataprocessing/newfigures/rFpredictions.tiff',
#     width=7,height=8,units='in',res=200)
#tiff(filename='~/lab/Stanford/writings/Thesis/forpub/Paper2/Paper2_finalfigs/rf_predictions.tif',
#     width=7,height=8,units='in',res=200)
setEPS()
postscript("~/Stanford/writings/Thesis/forpub/Paper2/EM_second_revisions/Fig4/rFpredictions.eps", width=7, height=8)
#pdf("~/Stanford/writings/Thesis/forpub/Paper2/EM_second_revisions/Fig4/rFpredictions.pdf", width=7, height=8)
par(mfrow=c(2,2))

barplot(rFsite_predictions_proportions, 
     col=c("#663300","#FF0000", "#FFCC00","#66CC00","#0000CC"),
     xlab=expression(bold(Sample~Site)), ylab=expression(bold(Fraction~of~predicted~sites)))

plot(rFSal_predictions_summary$prediction~rFSal_predictions_summary$Sal, pch=16, 
     xlim=c(0, 35), ylim=c(0, 35),
     col=rFSal_predictions_summary$color, xlab=expression(bold(Sample~Salinity~(PSU))), ylab=expression(bold(Predicted~Salinity~(PSU))))
abline(lm(rFSal_predictions_summary$prediction~rFSal_predictions_summary$Sal))

plot(rFNtot_predictions_summary$prediction~rFNtot_predictions_summary$Ntot, pch=16, 
     xlim=c(0, 0.12), ylim=c(0,0.12),
     col=rFNtot_predictions_summary$color, xlab=expression(paste(bold(Sample~Total~N),bold(' (%)'))), ylab=expression(paste(bold(Predicted~Total~N),bold(' (%)'))))
abline(lm(rFNtot_predictions_summary$prediction~rFNtot_predictions_summary$Ntot))

plot(rFTemp_predictions_summary$prediction~rFTemp_predictions_summary$Temp, pch=16,
     xlim=c(8,20), ylim=c(8,20),
     col=rFTemp_predictions_summary$color, xlab=expression(bold(Sample~Temperature~(ºC))), ylab=expression(bold(Predicted~Temperature~(ºC))))
abline(lm(rFTemp_predictions_summary$prediction~rFTemp_predictions_summary$Temp))

dev.off()

summary(lm(rFSal_predictions_summary$prediction~rFSal_predictions_summary$Sal))
anova(lm(rFSal_predictions_summary$prediction~rFSal_predictions_summary$Site+rFSal_predictions_summary$Sal))
summary(lm(rFNtot_predictions_summary$prediction~rFNtot_predictions_summary$Ntot))
anova(lm(rFNtot_predictions_summary$prediction~rFNtot_predictions_summary$Site*rFNtot_predictions_summary$Ntot))
anova(lm(rFNtot_predictions_summary$prediction~rFNtot_predictions_summary$Site))
summary(lm(rFTemp_predictions_summary$prediction~rFTemp_predictions_summary$Temp))
anova(lm(rFTemp_predictions_summary$prediction~rFTemp_predictions_summary$Site*rFTemp_predictions_summary$Temp))


##################################
# variable importance
uniqueOTUs_site<-as.vector(read.csv('rF_site_topimportantOTUs.csv', row.names=1)[,1])
uniqueOTUs_Sal<-as.vector(read.csv('rF_Sal_topimportantOTUs.csv', row.names=1)[,1])
uniqueOTUs_Ntot<-as.vector(read.csv('rF_Ntot_topimportantOTUs.csv', row.names=1)[,1])

##################################
#opar<-par()
## barplots of OTUs that are variable by Site
# to get nirSotumat.site, see nirSresamp_randomForest.R
# a new way to re-order the sites, for plotting:
nirSotumat.site2<-nirSotumat.site[order(nirSotumat.site$Site, rownames(nirSotumat.site)),]
# now subset just to get the important OTUs
Site_important<-subset(nirSotumat.site2, select=uniqueOTUs_site[1:10])
Site_interesting<-subset(nirSotumat.site2, select=)
# and color scheme 2:
site_colors<-rep(c("#663300","#FF0000","#FFCC00","#66CC00","#0000CC"),each=7)

for (i in 1:length(Site_important)){
  barplot(Site_important[,i], main=colnames(Site_important)[i], 
          names.arg=rownames(Site_important), cex.names=.6, las=2,
          col=site_colors,
          ylab='abundance')
}

# to plot OTUs that are not important to site, eg 488, 533, 1849, 1460, 61, 463
barplot(nirSotumat.site2$OTU_1140, main='OTU_1140', 
        names.arg=rownames(nirSotumat.site2), cex.names=.6, las=2,
        col=site_colors,
        ylab='abundance')

#par(opar)

site_colors_norep<-rep(c("#663300","#FF0000","#FFCC00","#66CC00","#0000CC"),7)
## plots of OTUs that are variable by Sal
Sal_important<-subset(otumat.sal, select=c('Sal',uniqueOTUs_Sal[1:10]))
for (i in 2:length(Sal_important)){
  plot(Sal_important[,i]~Sal_important$Sal, main=colnames(Sal_important)[i], 
       xlab='Salinity (PSU)', ylab='abundance', pch=16, col=site_colors_norep)
}

# and unimportant ones: 1666, 67, 1519, 300, 775, 985
head(rFsal_train1$importance[order(rFsal_train1$importance[,1]),])
plot(otumat.sal$OTU_985~otumat.sal$Sal, main='OTU_985', 
     xlab='Salinity (PSU)', ylab='abundance', pch=16, col=site_colors_norep)



## plots of OTUs that are variable by Ntot
Ntot_important<-subset(otumat.Ntot, select=c('Ntot',uniqueOTUs_Ntot[1:10]))
for (i in 2:length(Ntot_important)){
  plot(Ntot_important[,i]~Ntot_important$Ntot, main=colnames(Ntot_important)[i], 
       xlab='total sediment N', ylab='abundance', pch=16, col=site_colors_norep)
}
# and unimportant ones: 1623, 139, 1575, 1519, 1044, 811
head(rFNtot_train1$importance[order(rFNtot_train1$importance[,1]),])
plot(otumat.Ntot$OTU_1044~otumat.Ntot$Ntot, main='OTU_1044', 
     xlab='total sediment N (%)', ylab='abundance', pch=16, col=site_colors_norep)



## plots of OTUs that are variable by Temp: 391, 985, 1584, 397, 1609, 277
tail(rFTemp_train1$importance[order(rFTemp_train1$importance[,1]),])
plot(otumat.Temp$OTU_1609~otumat.Temp$Temp, main='OTU_1609', 
     xlab='Temperature (C)', ylab='abundance', pch=16, col=site_colors_norep)


############################################

# more phyloseq
# from nirS88_bigtree.R: bigtree_PGMimportant<-subset_taxa(bigtree, rF %in% c('Ntot','Sal','Site'))

bigtree_PGMsite<-subset_taxa(bigtree, rF=='Site')
sample_data(bigtree_PGMsite)<-sample_data(read.csv('~/lab/writings/Thesis/finaldataforpaper/Ch2_dataprocessing/nirS_resamp_sampledata.csv', header=T, row.names=1))
sample_data(bigtree_PGMsite)$Month<-factor(sample_data(bigtree_PGMsite)$Month, levels=c("Jul","Sep","Oct","Nov","Jan","Mar","May",'cultured'))
sample_data(bigtree_PGMsite)$Site<-factor(sample_data(bigtree_PGMsite)$Site, levels=c('4.1','8.1','13','21','24','cultured'))
Polaris_color_tree<-c("#FFCC00","#66CC00","#0000CC","#663300","#FF0000")
plot_tree(bigtree_PGMsite, color='Site', size='abundance', ladderize='left') + scale_color_manual(values=Polaris_color_tree)

bigtree_PGMsite_refs<-subset_taxa(bigtree, rF=='Site' | sampletype=='cultured')
sample_data(bigtree_PGMsite_refs)<-sample_data(read.csv('~/lab/writings/Thesis/finaldataforpaper/Ch2_dataprocessing/bigtree/nirS_bigtree2_sampledata.csv', header=T, row.names=1))
sample_data(bigtree_PGMsite_refs)$Month<-factor(sample_data(bigtree_PGMsite_refs)$Month, levels=c("Jul","Sep","Oct","Nov","Jan","Mar","May",'cultured'))
sample_data(bigtree_PGMsite_refs)$Site<-factor(sample_data(bigtree_PGMsite_refs)$Site, levels=c('4.1','8.1','13','21','24','cultured'))
bigtree_PGMsite_plot<-plot_tree(bigtree_PGMsite_refs, color='Site', size='abundance', ladderize='left', label.tips='Species') + scale_color_manual(values=c(Polaris_color_tree, 'darkgray'))
tiff('rFtree_site.tiff', width=6, height=10, units='in',res=200)
bigtree_PGMsite_plot
dev.off()

unifracdist_PGMsite = phyloseq::distance(bigtree_PGMsite, method="unifrac", weighted=TRUE)
unifracord_PGMsite<-ordinate(bigtree_PGMsite, "PCoA", unifracdist_PGMsite)
unifracplot_PGMsite<-plot_ordination(bigtree_PGMsite, unifracord_PGMsite, color="Site", shape='Month') 
unifracplot_PGMsite<-unifracplot_PGMsite + geom_point(size = 4, fill='white') 
unifracplot_PGMsite<-unifracplot_PGMsite + scale_colour_manual(values=c(Polaris_color)) + scale_shape_manual(values=Polaris_shape) + ggtitle("PCoA, weighted Unifrac\n 51 Site-associated OTUs")
unifracplot_PGMsite

adonis(unifracdist_PGMsite ~ Site, as(sample_data(bigtree_PGMsite), "data.frame")) 
adonis(unifracdist_PGMsite ~ NO3+Temp+Ntot+NH4+Sal+Fe, as(sample_data(bigtree_PGMsite), "data.frame")) 

##
bigtree_PGMSal<-subset_taxa(bigtree, rF=='Sal')
bigtree_PGMSal_refs<-subset_taxa(bigtree, rF=='Sal' | sampletype=='cultured')
sample_data(bigtree_PGMSal_refs)<-sample_data(read.csv('~/lab/writings/Thesis/finaldataforpaper/Ch2_dataprocessing/bigtree/nirS_bigtree2_sampledata.csv', header=T, row.names=1))
sample_data(bigtree_PGMSal_refs)$Month<-factor(sample_data(bigtree_PGMSal_refs)$Month, levels=c("Jul","Sep","Oct","Nov","Jan","Mar","May",'cultured'))
sample_data(bigtree_PGMSal_refs)$Site<-factor(sample_data(bigtree_PGMSal_refs)$Site, levels=c('4.1','8.1','13','21','24','cultured'))
plot_tree(bigtree_PGMSal, color='Sal', size='abundance', ladderize='left') + scale_color_gradient(limits=c(0, 32), low="#CC3300", high="#0066CC")
plot_tree(bigtree_PGMSal, color='Site', size='abundance', ladderize='left') + scale_color_manual(values=Polaris_color_tree)
bigtree_PGMSalplot<-plot_tree(bigtree_PGMSal_refs, color='Sal', size='abundance', ladderize='left', label.tips='Species') 
tiff('rFtree_Sal.tiff', width=6, height=10, units='in',res=200)
bigtree_PGMSalplot
dev.off()


unifracdist_PGMSal = phyloseq::distance(bigtree_PGMSal, method="unifrac", weighted=TRUE)
unifracord_PGMSal<-ordinate(bigtree_PGMSal, "PCoA", unifracdist_PGMSal)
unifracwtplot_PGMSal<-plot_ordination(bigtree_PGMSal, unifracord_PGMSal, color="Sal") + scale_color_gradient(limits=c(0, 32), low="#CC3300", high="#0066CC") + geom_point(size = 4) + ggtitle("PCoA, weighted Unifrac\n 46 Salinity-associated OTUs")
unifracwtplot_PGMSal

adonis(unifracdist_PGMSal ~ Site, as(sample_data(bigtree_PGMsite), "data.frame")) 
adonis(unifracdist_PGMSal ~ NO3+Temp+Ntot+NH4+Sal+Fe, as(sample_data(bigtree_PGMsite), "data.frame")) 
row.names(sample_data(bigtree_PGMSal))==row.names(unifracord_PGMSal$vectors) #just to make sure
summary(lm(unifracord_PGMSal$vectors[,1]~sample_data(bigtree_PGMSal)$Sal))

##
bigtree_PGMNtot<-subset_taxa(bigtree, rF=='Ntot')
sample_data(bigtree_PGMNtot)<-sample_data(read.csv('~/lab/writings/Thesis/finaldataforpaper/Ch2_dataprocessing/nirS_resamp_sampledata.csv', header=T, row.names=1))
sample_data(bigtree_PGMNtot)$Month<-factor(sample_data(bigtree_PGMNtot)$Month, levels=c("Jul","Sep","Oct","Nov","Jan","Mar","May",'cultured'))
sample_data(bigtree_PGMNtot)$Site<-factor(sample_data(bigtree_PGMNtot)$Site, levels=c('4.1','8.1','13','21','24','cultured'))

bigtree_PGMNtot_refs<-subset_taxa(bigtree, rF=='Ntot' | sampletype=='cultured')
sample_data(bigtree_PGMNtot_refs)<-sample_data(read.csv('~/lab/writings/Thesis/finaldataforpaper/Ch2_dataprocessing/bigtree/nirS_bigtree2_sampledata.csv', header=T, row.names=1))
sample_data(bigtree_PGMNtot_refs)$Month<-factor(sample_data(bigtree_PGMNtot_refs)$Month, levels=c("Jul","Sep","Oct","Nov","Jan","Mar","May",'cultured'))
sample_data(bigtree_PGMNtot_refs)$Site<-factor(sample_data(bigtree_PGMNtot_refs)$Site, levels=c('4.1','8.1','13','21','24','cultured'))
bigtree_PGMNtot_plot<-plot_tree(bigtree_PGMNtot_refs, color='Ntot', size='abundance', ladderize='left', label.tips='Species') + scale_color_gradient(low="#CCFF66", high="#336600") 
tiff('rFtree_Ntot.tiff', width=6, height=10, units='in',res=200)
bigtree_PGMNtot_plot
dev.off()

unifracdist_PGMNtot = phyloseq::distance(bigtree_PGMNtot, method="unifrac", weighted=TRUE)
unifracord_PGMNtot<-ordinate(bigtree_PGMNtot, "PCoA", unifracdist_PGMNtot)
unifracwtplot_PGMNtot<-plot_ordination(bigtree_PGMNtot, unifracord_PGMNtot, color="Ntot") + scale_color_gradient(low="#CCFF66", high="#336600")  + geom_point(size = 4) + ggtitle("PCoA, weighted Unifrac\n 38 Nitrogen-associated OTUs")
unifracwtplot_PGMNtot

adonis(unifracdist_PGMNtot ~ Site, as(sample_data(bigtree_PGMsite), "data.frame")) 
adonis(unifracdist_PGMNtot ~ NO3+Temp+Ntot+NH4+Sal+Fe, as(sample_data(bigtree_PGMsite), "data.frame")) 
row.names(sample_data(bigtree_PGMNtot))==row.names(unifracord_PGMNtot$vectors) #just to make sure
summary(lm(unifracord_PGMNtot$vectors[,1]~sample_data(bigtree_PGMNtot)$Ntot))
