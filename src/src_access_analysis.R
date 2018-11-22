### ANALYSIS PART

# load packages, data, and options
library(ggplot2)
library(raster)
library(geosphere)
source("./src/geo_functions.R")
# LOAD data
poly_df = raster::shapefile("./input/polygons")
parkrun_marker = raster::shapefile("./input/marker_england")
parkrun_marker = get_event_res(event = parkrun_marker,
                               poly_pop = poly_df,
                               pop_n = poly_df$pop,
                               objective = "dist^2 *pop")



### ANALYSIS PART
# tbc
options(scipen = 999)


# plotting and analyses
hist_dist_per_lsoa  = ggplot() + 
  geom_histogram(aes(poly_df$mn_dstn),bins=100,fill = "cyan", col ="cyan", alpha = 0.7) + theme_minimal() +
  scale_x_continuous(labels = function(x) format(round(x),big.mark = ","),
                     name = "Distance from LSOA to the nearest event (in km)",
                     limits = c(0,40)) +
  ylab(paste("Frequency (N =",length(poly_df$mn_dstn),")",sep="")) +
  geom_vline(xintercept = median(poly_df$mn_dstn), color="orange", alpha = 0.6) +
  geom_text(aes(x=median(poly_df$mn_dstn)+5,y=2300),label=paste("Median =",median(poly_df$mn_dstn),"km"),color="orange")


temp_dist = poly_df$mn_dstn[order(poly_df$mn_dstn,decreasing = F)]
cum.n = (1:length(poly_df$mn_dstn))/length(poly_df$mn_dstn)

cumsum_dist_per_lsoa =
  ggplot() +
  geom_line(aes(x=temp_dist, y=cum.n )) +
  scale_y_continuous(minor_breaks=seq(0,1,by=0.1),
                     breaks = c(0,0.5,0.9,1)) +
  scale_x_continuous(minor_breaks=seq(0,30,by=2.5),
                     breaks = c(0,5,10,20,30),
                     limits = c(0,30)) +
  theme_minimal() +
  ylab("Cumulative proportion of LSOA") +
  xlab("Distance to the nearest event") 


# population per event
hist_srvd_pop_per_event = ggplot() + geom_histogram(aes(parkrun_marker$srvd_pop),bins=50,fill = "cyan", col ="cyan", alpha = 0.7) + theme_minimal() +
  scale_x_continuous(labels = function(x) paste(format(round(x/1000),big.mark = ","),"k")) +
  geom_vline(xintercept = median(parkrun_marker$srvd_pop ), size=2,color="orange", alpha = 0.6)  +
  ylab("Frequency") + 
  xlab("Served population per parkrun event (in thousand)")





steps1 = seq(min(poly_df$a),max(poly_df$a),length.out = 100)
temp.cut = cut(poly_df$a,steps1)
temp.group.means1 = by(data=poly_df$mn_dstn,INDICES = temp.cut,function(x) mean(x,na.rm=T))


lm_dist_depri = lm(mn_dstn ~ a, data=poly_df)

ma_dist_depri = 
  ggplot() +
  theme_minimal() +
  geom_point(aes(x=poly_df$a,y=poly_df$mn_dstn), size = 0.5, alpha = 0.9,color="pink") +
  geom_smooth(aes(x=poly_df$a,y=poly_df$mn_dstn, col ="Model"),method = 'lm' , formula= 'y ~ x', se = F, size = 2) +
  geom_line(aes(x=steps1[-1],y=as.numeric(temp.group.means1), col = "Moving average"),size=2)+
  geom_point(aes(x=poly_df$a,y=poly_df$mn_dstn), size = 0.5, alpha = 0.9, color ="lightblue") +
  geom_smooth(aes(x=poly_df$a,y=poly_df$mn_dstn, col ="Model"),method = 'lm' , formula= 'y ~ x', se = F) +
  geom_line(aes(x=steps1[-1],y=as.numeric(temp.group.means1), col = "Moving average"))+
  theme(legend.position="bottom") +
  guides(col=guide_legend(title="")) +
  ylab("Distance to the nearest parkrun event location") +
  xlab("Multiple deprivation score") 


steps2 = seq(min(log(poly_df$pp_dnst)),max(log(poly_df$pp_dnst)),length.out = 100)
temp.cut = cut(log(poly_df$pp_dnst),steps2)
temp.group.means2 = by(data=(poly_df$mn_dstn),INDICES = temp.cut,function(x) mean(x,na.rm=T))

lm_dist_density = lm(mn_dstn ~ poly(log(pp_dnst),4),data=poly_df)

ma_dist_density = 
  ggplot() +
  theme_minimal() +
  geom_point(aes(x=log(poly_df$pp_dnst),y=poly_df$mn_dstn), size = 0.5, alpha = 0.9, color ="lightgreen") +
  geom_line(aes(x=log(poly_df$pp_dnst),y=predict(lm_dist_density),col ="Model"),size =2) +
  geom_line(aes(x=steps2[-1],y=as.numeric(temp.group.means2), col = "Moving average"), size = 2)+
  geom_point(aes(x=log(poly_df$pp_dnst),y=poly_df$mn_dstn), size = 0.5, alpha = 0.7, color ="lightgreen") +
  geom_line(aes(x=log(poly_df$pp_dnst),y=predict(lm_dist_density),col ="Model")) +
  geom_line(aes(x=steps2[-1],y=as.numeric(temp.group.means2), col = "Moving average"))+
  theme(legend.position="bottom") +
  guides(col=guide_legend(title="")) +
  ylab("Distance to the nearest parkrun event location") +
  xlab("Population density (log scale)") +
  ggtitle("Association between pop. density and access to parkrun")



steps3 = seq(min(log(poly_df$pp_dnst)),max(log(poly_df$pp_dnst)),length.out = 100)
temp.cut = cut(log(poly_df$pp_dnst),steps3)
temp.group.means3 = by(data=(poly_df$a),INDICES = temp.cut,function(x) mean(x,na.rm=T))

lm_a_density = lm(a ~ poly(log(pp_dnst),2),data=poly_df)

ma_depri_density = 
  ggplot() +
  geom_point(aes(x=log(poly_df$pp_dnst),y=poly_df$a), size = 0.5, alpha = 0.9, color ="pink") +
  geom_line(aes(x=log(poly_df$pp_dnst),y=predict(lm_a_density),col ="Model")) +
  geom_line(aes(x=steps3[-1],y=as.numeric(temp.group.means3), col = "Moving average"))+
  theme(legend.position="bottom") +
  theme_minimal() +
  guides(col=guide_legend(title="")) +
  ylab("Multiple deprivation score") +
  xlab("Population density (log scale)") +
  ggtitle("Association between pop. density and deprivation")


lm_dist_depri_dens = lm(mn_dstn ~ a * poly(log(pp_dnst),4),data=poly_df)
new.df = data.frame(a = rep(seq(from=min(poly_df$a),to=max(poly_df$a),length.out = 10),
                            times = 3),
                    pp_dnst = rep(quantile(poly_df$pp_dnst,c(0.2,0.5,0.8)),
                                  each = 10),
                    group = rep(c("low density",
                                  "medium density",
                                  "high density"),
                                each = 10)
)
new.preds = predict.lm(lm_dist_depri_dens,newdata=new.df)


model_est_a_density_dist = 
  ggplot() +
  theme_minimal() +
  geom_line(aes(x=new.df$a,y=new.preds, col =new.df$group)) +
  theme(legend.position="bottom") +
  guides(col=guide_legend(title="")) +
  ylab("Distance to nearest parkrun event") +
  xlab("Multiple deprivation score") 

