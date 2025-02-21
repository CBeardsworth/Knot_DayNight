################################################
# Figure 1 #####
##############################################
library(ggplot2)
library(sf)
library(ggspatial)
library(data.table)
library(dplyr)
library(cowplot)
library(scales)
library(stars)

bath <- read_sf("C:/Users/BESCBEAR/OneDrive - NIOZ/2_Reference Data/Basemaps/utm31 Bathy/DutchWaddenBathymetry-144_LAT.shp")#srows 6,17,21,25 seem to be corrupted.
coast_NL <- read_sf("C:/Users/BESCBEAR/OneDrive - NIOZ/2_Reference Data/Basemaps/WaddenCoast.gpkg")#%>% #https://www.eea.europa.eu/data-and-maps/data/eea-coastline-for-analysis-2/gis-data/eea-coastline-polygon 
coast <- read_sf("C:/Users/BESCBEAR/OneDrive - NIOZ/2_Reference Data/Basemaps/utm31 Bathy/DutchWaddenBathymetry100.shp")


griend <- st_as_sfc(st_bbox(c(xmin = 645464, xmax = 656867, ymin = 5900456, ymax = 5908406), crs=32631)) #for box on fig 1
box <- st_as_sfc(st_bbox(c(xmin = 665000, xmax = 670000, ymin = 5912923, ymax = 5922923), crs=32631)) #for box on fig 1

wadden <- ggplot() +  
  geom_sf(data = bath[c(1:5,7:16,18:20,22:24,26:nrow(bath)),], aes(fill="Mudflat"),col="rosybrown", lwd=0.3) +
  geom_sf(data = coast_NL, aes(fill="Land"), col="grey22", lwd=0.1) +
  #geom_sf(data=receivers, fill="red", aes(shape="Receivers"), size = 3)+
  #geom_sf(data=griend, col="blue", fill=NA)+
  coord_sf(xlim= c(605107,670000), ylim=c(5860000,5922923))+
  annotate(geom = "text", x = 615200, y = 5915100, label = "North Sea", 
           fontface = "italic", color = "black", size = 4) +
  annotate(geom = "text", x = 635000, y = 5880100, label = "Wadden Sea", 
           fontface = "italic", color = "black", size = 4) +
  annotate(geom = "text", x = 647900, y = 5907800, label = "Griend", 
           fontface = "italic", color = "black", size = 3) +
  annotate(geom = "text", x = 639200, y = 5911100, label = "Richel", 
           fontface = "italic", color = "black", size = 3) +
  geom_segment(aes(x = 647900, y = 5906800, xend = 649700, yend = 5903400), linewidth=1.2, #Griend arrow
               arrow = arrow(type="closed", length= unit(0.15, "cm")))+
  geom_segment(aes(x = 639200, y = 5910100, xend = 641500, yend = 5908000), linewidth = 1.2,#Richel Arrow
                            arrow = arrow(type="closed", length= unit(0.15, "cm")))+
  #geom_sf_text(data=griend, label="Griend Mudflat", size=2.5, nudge_x=-2000, nudge_y=3200, col="blue")+
  xlab("Longitude")+
  ylab("Latitude")+
  annotation_north_arrow(location="br",height= unit(0.8,"cm"), width= unit(0.6, "cm"), pad_y=unit(0.8, "cm"),
                         style=north_arrow_orienteering(text_size=8))+
  annotation_scale(location="br", width_hint=0.2)+
  scale_fill_manual(values = c("Land" = "grey65", "Mudflat" = "papayawhip"), 
                    name = NULL) +
  scale_shape_manual(values = c("Receivers" = 24), 
                     name = NULL) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill=NA),
        legend.spacing.y = unit(-0.15, "cm"),
        legend.key=element_blank())
wadden
#ggsave("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/figures/Fig_1.eps", wadden,units="cm", dpi=1200, height=16, width =16, device="eps")

cairo_ps("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/figures/Fig_1.eps", wadden, width = 7, height = 7, pointsize= 12, fallback_resolution = 1200)
print(wadden)
dev.off()

################################################################################
######################### Figure 2 #############################################
################################################################################

bird2117_535 <- readRDS("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/islandica-tag_2117-tide_2020535-date_2020-10-03.rds")
bird2117_536 <- readRDS("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/islandica-tag_2117-tide_2020536-date_2020-10-04.rds")

b535xy_res <- rbindlist(bird2117_535$residence_patches$patchdata) %>% st_as_sf(coords = c("x","y"), crs = 32631, remove = F) 
b535xy_com <- bird2117_535$travel_segments %>% st_as_sf(coords = c("x","y"), crs = 32631, remove = F) 
b535xy <- rbind(b535xy_res[,c("time","geometry")], b535xy_com[,c("time","geometry")]) %>%  arrange(time)
b535res <- do.call(c,bird2117_535$residence_patches$polygons) %>%  st_as_sf(.,crs = 32631) %>% mutate(patch_id = seq(1:nrow(.)), duration = bird2117_535$residence_patches$duration/60)
b535line <- b535xy %>% summarize(do_union=F) %>% st_cast("LINESTRING")

b536xy_res <- rbindlist(bird2117_536$residence_patches$patchdata) %>% st_as_sf(coords = c("x","y"), crs = 32631, remove = F) 
b536xy_com <- bird2117_536$travel_segments %>% st_as_sf(coords = c("x","y"), crs = 32631, remove = F) 
b536xy <- rbind(b536xy_res[,c("time","geometry")], b536xy_com[,c("time","geometry")]) %>%  arrange(time)
b536res <- do.call(c,bird2117_536$residence_patches$polygons) %>%  st_as_sf(.,crs = 32631) %>% mutate(patch_id = seq(1:nrow(.)), duration = bird2117_536$residence_patches$duration/60)
b536line <- b536xy %>% summarize(do_union=F) %>% st_cast("LINESTRING")

dayresplot <- ggplot() +  
  geom_sf(data = bath, fill="floralwhite",col="grey50", lwd=0.1) +
  geom_sf(data = coast,fill="grey", col="grey40", lwd=0.1) +
  geom_sf(data=b535res, aes(fill=duration), col="grey20",alpha = 0.9) +
  geom_sf(data=b535line, col="grey")+
  geom_sf(data=b535xy, col="grey1",pch = 18,size =0.8) +
  coord_sf(xlim= c(648864,654867), ylim=c(5901456,5905706))+ #griend only
  xlab("Longitude")+
  ylab("Latitude")+
  scale_fill_gradient2(limits = c(4, 250), low="white", high="blue")+
  annotation_north_arrow(location="br",height= unit(0.6,"cm"), width= unit(0.5, "cm"), pad_y=unit(0.8, "cm"),
                         style=north_arrow_orienteering(text_size=6))+
  annotation_scale(location="br", width_hint=0.2)+
  scale_y_continuous(breaks=c(53.240, 53.250, 53.260,53.270, 53.280))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill="NA"),
        #legend.key=element_blank(), 
        legend.position = "none",
        #legend.box.background = element_rect(colour="black"),
        legend.spacing.y = unit(-0.1, "cm"))
dayresplot

nightresplot <- ggplot() +  
  geom_sf(data = bath, fill="floralwhite",col="grey50", lwd=0.1) +
  geom_sf(data = coast,fill="grey", col="grey40", lwd=0.1) +
  geom_sf(data=b536res, aes(fill=duration), col="grey20",alpha = 0.9) +
  geom_sf(data=b536line, col="grey")+
  geom_sf(data=b536xy, col="grey2", pch = 18,size =0.8) +
  coord_sf(xlim= c(648864,654867), ylim=c(5901456,5905706))+ #griend only
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="duration \n (mins)")+
  scale_fill_gradient2(limits = c(4, 250), low="white", high="blue")+
  annotation_north_arrow(location="br",height= unit(0.6,"cm"), width= unit(0.5, "cm"), pad_y=unit(0.8, "cm"),
                         style=north_arrow_orienteering(text_size=6))+
  annotation_scale(location="br", width_hint=0.2)+
  scale_y_continuous(breaks=c(53.240, 53.250, 53.260,53.270, 53.280))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill="NA"),
        #legend.key=element_blank(), 
        legend.position = "right",
        legend.frame = element_rect(colour="grey40"),
        #legend.box.background = element_rect(colour="black"),
        legend.spacing.y = unit(-0.1, "cm"))#+

legend <- get_plot_component(nightresplot, 'guide-box-right', return_all = TRUE)

nightresplot<- nightresplot + theme(legend.position = "none")
res_patch_plots <- plot_grid(plot_grid(dayresplot, nightresplot, ncol=1,labels="auto"), legend, ncol=2,rel_widths = c(1,0.2), background=)
res_patch_plots

#ggsave("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/figures/Fig_2_DayNight_Examplebirds.eps", res_patch_plots, width =15.5, height = 17, units="cm",dpi=1200, device=cairo_eps)

cairo_ps("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/figures/Fig_2_DayNight_Examplebirds.eps", res_patch_plots, width = 7, height = 7, pointsize= 12, fallback_resolution = 1200)
print(res_patch_plots)
dev.off()
####################################################################################
############################ Create figure 4
###################################################################################

# Figure 4 shows the mean revisiting rate of all birds (per year and per day/night) to 75mx75m cells across a raster of the griend mudflat. revisiting rates are calculated first per bird by checking overlap with 75x75 cells (>50% = visit). More than 2 visits across multiple tides counts as a revisited cell. Visits per tide are calculated per bird and a proportion is calculated based on how many tides they visited. We then calculate the mean across all birds. Cells that are not revisited by a bird don't count towards a mean.

griend_poly <- read_sf("C:/Users/BESCBEAR/OneDrive - NIOZ/2_Reference Data/Griend/GriendHighTide_2015elevation_2018Tides.gpkg")

tides <- read.csv("C:/Users/BESCBEAR/OneDrive - NIOZ/2_Reference Data/Water/Griend/Wide/HLH/allYears-tidalPattern-HLH-griend_Wide-UTC.csv")%>%
  mutate(year = substring(tideID,1,4),
         yearly_id = as.numeric(substring(tideID, 6,8)),
         tideID = as.numeric(sub("_", "", tideID)),
         high_start_time = as.POSIXct(high_start_time, format = "%Y-%m-%d %H:%M:%S"),
         date = as.Date(high_start_time, tz="UTC"),
         month = month(high_start_time)) %>%     #allert saved dates on the files in CET (I think) so need to change this.  
  left_join(SA, by = c("low_level"="waterlevel")) %>% #space available in the tide
  left_join(waterdn, by="tideID") %>% 
  filter(month %in% c(9,10), tod_summary!="twilight" & start_tod == tod_summary) %>% #choose only after September and when start time of day is the same as majority (not twilight)
  rename(tide_tod = start_tod)

b_tides <- st_read("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/2018-2022_bird_tides_multipolygons_min5tides.gpkg") %>% 
  mutate(bird_tides = paste(tideID, id, sep="_")) %>% 
  filter(tideID %in% tides$tideID)

fig_data <- b_tides %>% 
  group_by(tide_tod, id, year) %>% 
  summarise(n_tides = n(), 
            geom=st_union(geom)) %>% 
  ungroup()

summary_fig <- fig_data %>% 
  group_by(tide_tod, year) %>% 
  summarise(n_birds = n(),
            min_tides = min(n_tides),
            mean_tides = mean(n_tides),
            median_tides = median(n_tides),
            max_tides = max(n_tides),
            geom=st_union(geom))
summary_fig$geom <- NULL             

summary_fig
bath <- read_sf("C:/Users/BESCBEAR/OneDrive - NIOZ/2_Reference Data/Basemaps/utm31 Bathy/DutchWaddenBathymetry-144_LAT.shp") #srows 6,17,21,25 seem to be corrupted.

ras <- read_sf("C:/Users/BESCBEAR/OneDrive - NIOZ/2_Reference Data/Wadden Regions/griendonly_bufferapx1km.gpkg")  %>% 
  st_rasterize(dx=75,dy=75, inside = T)
template <- ras-1

for(y in 2018:2021){
  
  day_ras <- template
  night_ras <- template
  r_day = list()
  r_night = list()
  print(y)
  for(tod in c("day", "night")){
    print(tod)
    for(ind in unique(b_tides[b_tides$year==y,][["id"]])){
      sub <- b_tides %>%
        filter(year == y, tide_tod == tod, id == ind)
      n_tides <- nrow(sub)
      ras_sub <- st_rasterize(sub[1,"geom"], template, options = "MERGE_ALG=ADD")
      
      for(row in 2:nrow(sub)){ # add all visits to an area from all birds
        ras_sub2 <- st_rasterize(sub[row,"geom"], template, options = "MERGE_ALG=ADD")
        ras_sub <- ras_sub + ras_sub2
      }
      ras_sub[ras_sub<=1]<- NA # remove non-revisited
      ras_sub <- ras_sub/n_tides
    
      if(tod=="day"){r_day[[ind]] <- ras_sub}
      if(tod=="night"){r_night[[ind]] <- ras_sub}
    }
  }
  
  r_day = do.call("c", r_day) 
  r_day = st_redimension(r_day)
  r_day = st_set_dimensions(r_day, names = c("x", "y", "id"))
  r_day = st_apply(r_day, c("x", "y"), mean,na.rm=T)
  
  r_night = do.call("c", r_night) 
  r_night = st_redimension(r_night)
  r_night = st_set_dimensions(r_night, names = c("x", "y", "id"))
  r_night = st_apply(r_night, c("x", "y"), mean, na.rm=T)
  
  print(paste("day: ",max(r_night$mean, na.rm=T)))
  print(paste("night: ", max(r_day$mean, na.rm=T)))
  print(paste("day: ",min(r_night$mean, na.rm=T)))
  print(paste("night: ", min(r_day$mean, na.rm=T)))
  
  day_plot <- ggplot() +  
    geom_sf(data = bath[c(22,30:nrow(bath)),], fill="floralwhite",col="grey50", lwd=0.1) +
    geom_stars(data = r_day, na.action=na.omit) +
    geom_sf(data = griend_poly, fill="grey65", col="grey40", lwd=0.1) +
    viridis::scale_fill_viridis(limits= c(0,0.8), na.value=NA)+
    # scale_fill_gradientn(colours = c("#440154FF","#1F968BFF","#FDE725FF"), 
    #                   values = rescale(c(0, 0.2,0.8)),
    #                 guide = "colorbar", limits=c(0,0.8))+
    xlab("Longitude")+
    ylab(paste0(y, "\n \n Latitude"))+
    coord_sf(xlim= c(645864,656867), ylim=c(5899956,5908406)) + #griend only
    annotation_north_arrow(location="br",height= unit(0.5,"cm"), width= unit(0.3, "cm"), pad_y=unit(0.8, "cm"),
                           style=north_arrow_orienteering(text_size=5)) +
    annotation_scale(location="br", width_hint=0.2)+
    annotate("text",  x=-Inf, y = Inf, label = paste("n birds = ", length(unique(b_tides[b_tides$year==y,][["id"]]))), vjust=2, hjust=-0.1, size = 2)+
    annotate("text",  x=-Inf, y = Inf, label = paste("n tides = ", length(unique(b_tides[b_tides$year==y & b_tides$tide_tod=="day",][["tideID"]]))), vjust=3.5, hjust=-0.1, size = 2)+
    annotate("text",  x=-Inf, y = Inf, label = paste("n bird-tides = ", length(unique(b_tides[b_tides$year==y & b_tides$tide_tod=="day",][[c("bird_tides")]]))), vjust=5, hjust=-0.065, size = 2)+
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
          panel.border= element_rect(colour="black", fill="NA"),
          legend.position="none", 
          axis.text.x=element_blank(), 
          axis.title.x=element_blank())#+
  #ggtitle(paste(y, "day", sep=" "))
  
  night_plot <-  ggplot() +  
    geom_sf(data = bath[c(22,30:nrow(bath)),], fill="floralwhite",col="grey50", lwd=0.1) +
    geom_stars(data = r_night, na.action=na.omit) +
    geom_sf(data = griend_poly, fill="grey65", col="grey40", lwd=0.1) +
    viridis::scale_fill_viridis(limits= c(0,0.8), na.value=NA)+
    #scale_fill_gradientn(colours = c("#440154FF","#1F968BFF","#FDE725FF"), 
    #                    values = rescale(c(0, 0.2,0.8)),
    #                   guide = "colorbar", limits=c(0,0.8))+
    xlab("Longitude")+
    #ylab("Latitude")+
    coord_sf(xlim= c(645864,656867), ylim=c(5899956,5908406)) + #griend only
    annotation_north_arrow(location="br",height= unit(0.5,"cm"), width= unit(0.3, "cm"), pad_y=unit(0.8, "cm"),
                           style=north_arrow_orienteering(text_size=5)) +
    annotation_scale(location="br", width_hint=0.2)+
    annotate("text",  x=-Inf, y = Inf, label = paste("n birds = ", length(unique(b_tides[b_tides$year==y,][["id"]]))), vjust=2, hjust=-0.1, size = 2)+
    annotate("text",  x=-Inf, y = Inf, label = paste("n tides = ", length(unique(b_tides[b_tides$year==y & b_tides$tide_tod=="night",][["tideID"]]))), vjust=3.5, hjust=-0.1, size = 2)+
    annotate("text",  x=-Inf, y = Inf, label = paste("n bird-tides = ", length(unique(b_tides[b_tides$year==y & b_tides$tide_tod=="night",][[c("bird_tides")]]))), vjust=5, hjust=-0.065, size = 2)+
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
          panel.border= element_rect(colour="black", fill="NA"),
          legend.position="none", 
          axis.text=element_blank(), 
          axis.title=element_blank())#+
  #ggtitle(paste(y,"night", sep=" "))
  
  
  if(y==2018){p_2018 <- cowplot::plot_grid(day_plot, night_plot, rel_widths = c(1,0.67))}
  if(y==2019){p_2019 <- cowplot::plot_grid(day_plot, night_plot, rel_widths = c(1,0.67))}
  if(y==2020){p_2020 <- cowplot::plot_grid(day_plot, night_plot, rel_widths = c(1,0.67))}
  if(y==2021){p_2021 <- cowplot::plot_grid(day_plot, night_plot, rel_widths = c(1,0.67))}
  
  
}

y=2022
day_ras <- template
night_ras <- template
r_day = list()
r_night = list()
print(y)
for(tod in c("day", "night")){
  print(tod)
  for(ind in unique(b_tides[b_tides$year==y,][["id"]])){
    sub <- b_tides %>%
      filter(year == y, tide_tod == tod, id == ind)
    n_tides <- nrow(sub)
    ras_sub <- st_rasterize(sub[1,"geom"], template, options = "MERGE_ALG=ADD")
    
    for(row in 2:nrow(sub)){ # add all visits to an area from all birds
      ras_sub2 <- st_rasterize(sub[row,"geom"], template, options = "MERGE_ALG=ADD")
      ras_sub <- ras_sub + ras_sub2
    }
    ras_sub[ras_sub<=1]<- NA # remove non-revisited
    ras_sub <- ras_sub/n_tides
    
    if(tod=="day"){r_day[[ind]] <- ras_sub}
    if(tod=="night"){r_night[[ind]] <- ras_sub}
  }
}

r_day = do.call("c", r_day) 
r_day = st_redimension(r_day)
r_day = st_set_dimensions(r_day, names = c("x", "y", "id"))
r_day = st_apply(r_day, c("x", "y"), mean,na.rm=T)

r_night = do.call("c", r_night) 
r_night = st_redimension(r_night)
r_night = st_set_dimensions(r_night, names = c("x", "y", "id"))
r_night = st_apply(r_night, c("x", "y"), mean, na.rm=T)

print(paste("day: ",max(r_night$mean, na.rm=T)))
print(paste("night: ", max(r_day$mean, na.rm=T)))
print(paste("day: ",min(r_night$mean, na.rm=T)))
print(paste("night: ", min(r_day$mean, na.rm=T)))

day_plot <-   ggplot() +  
  geom_sf(data = bath[c(22,30:nrow(bath)),], fill="floralwhite",col="grey50", lwd=0.1) +
  geom_stars(data = r_day, na.action=na.omit) +
  geom_sf(data = griend_poly, fill="grey65", col="grey40", lwd=0.1) +
  viridis::scale_fill_viridis(limits= c(0,0.8), na.value=NA)+
  #scale_fill_gradientn(colours = c("#440154FF","#1F968BFF","#FDE725FF"), 
  #                    values = rescale(c(0, 0.2,0.8)),
  #                   guide = "colorbar", limits=c(0,0.8))+
  xlab("Longitude \n \n Day")+
  ylab(paste0(y, "\n \n Latitude"))+
  coord_sf(xlim= c(645864,656867), ylim=c(5899956,5908406)) + #griend only
  annotation_north_arrow(location="br",height= unit(0.5,"cm"), width= unit(0.3, "cm"), pad_y=unit(0.8, "cm"),
                         style=north_arrow_orienteering(text_size=5)) +
  annotation_scale(location="br", width_hint=0.2)+
  annotate("text",  x=-Inf, y = Inf, label = paste("n birds = ", length(unique(b_tides[b_tides$year==y,][["id"]]))), vjust=2, hjust=-0.1, size = 2)+
  annotate("text",  x=-Inf, y = Inf, label = paste("n tides = ", length(unique(b_tides[b_tides$year==y & b_tides$tide_tod=="day",][["tideID"]]))), vjust=3.5, hjust=-0.1, size = 2)+
  annotate("text",  x=-Inf, y = Inf, label = paste("n bird-tides = ", length(unique(b_tides[b_tides$year==y & b_tides$tide_tod=="day",][[c("bird_tides")]]))), vjust=5, hjust=-0.065, size = 2)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill="NA"),
        legend.title = element_text("mean \n revisiting\ rate"),
        legend.position="none")+
  labs(fill = "mean %\nrevisits")
# ggtitle(paste(y, "day", sep=" "))

night_plot <-  ggplot() +  
  geom_sf(data = bath[c(22,30:nrow(bath)),], fill="floralwhite",col="grey50", lwd=0.1) +
  geom_stars(data = r_night, na.action=na.omit) +
  geom_sf(data = griend_poly, fill="grey65", col="grey40", lwd=0.1) +
  viridis::scale_fill_viridis(limits= c(0,0.8), na.value=NA)+
  #scale_fill_gradientn(colours = c("#440154FF","#1F968BFF","#FDE725FF"), 
  #                    values = rescale(c(0, 0.2,0.8)),
  #                   guide = "colorbar", limits=c(0,0.8))+
  xlab("Longitude \n \n Night")+
  #ylab("Latitude")+
  coord_sf(xlim= c(645864,656867), ylim=c(5899956,5908406)) + #griend only
  annotation_north_arrow(location="br",height= unit(0.5,"cm"), width= unit(0.3, "cm"), pad_y=unit(0.8, "cm"),
                         style=north_arrow_orienteering(text_size=5)) +
  annotation_scale(location="br", width_hint=0.2)+
  annotate("text",  x=-Inf, y = Inf, label = paste("n birds = ", length(unique(b_tides[b_tides$year==y,][["id"]]))), vjust=2, hjust=-0.1, size = 2)+
  annotate("text",  x=-Inf, y = Inf, label = paste("n tides = ", length(unique(b_tides[b_tides$year==y & b_tides$tide_tod=="night",][["tideID"]]))), vjust=3.5, hjust=-0.1, size = 2)+
  annotate("text",  x=-Inf, y = Inf, label = paste("n bird-tides = ", length(unique(b_tides[b_tides$year==y & b_tides$tide_tod=="night",][[c("bird_tides")]]))), vjust=5, hjust=-0.065, size = 2)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill="NA"),
        legend.position="none", 
        axis.text.y=element_blank(), 
        axis.title.y=element_blank())#+
#ggtitle(paste(y,"night", sep=" "))


if(y==2022){p_2022 <- cowplot::plot_grid(day_plot, night_plot, rel_widths = c(1,0.67))}

leg <- get_legend(day_plot + theme(legend.position="right", legend.key.height = unit(1, "cm")))

#cowplot::plot_grid(day_plot, leg)
p_all <- cowplot::plot_grid(plot_grid(p_2018,p_2019, p_2020, p_2021, p_2022, ncol=1, rel_heights = c(0.695,0.695,0.695,0.695,1)),leg, ncol=2, rel_widths=c(1,0.2))

#ggsave("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/figures/Fig_3.eps", p_all, dpi=1200, height = 8.5, width = 6, unit = "in", device="eps")

cairo_ps("C:/Users/BESCBEAR/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/figures/Fig_4.eps",  width = 7, height = 9, pointsize=12, fallback_resolution = 1200)
print(p_all)
dev.off()

