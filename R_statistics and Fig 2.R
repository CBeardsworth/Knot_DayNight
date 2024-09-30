# Aim: make raster per tide per individual where polygons are imported and made into a raster (polygon must overlap over 50% with cell) The data given to the raster is duration. 

#----------------------------
#import libraries
library(ggspatial)
library(cowplot)
library(stringr)
library(dplyr)
library(raster)
library(sf)
library(maptools)
library(stars)
library(lubridate)
library(data.table)
library(ggplot2)
library(purrr)
library(lme4)
library(effects)
library(introdataviz) # split violin plots
library(extrafont)

#--------------------------
Sys.setenv(TZ='UTC') # change time zone for session
setwd("D:/OneDrive - NIOZ/")
#setwd("C:/Users/cbeardsworth/OneDrive - NIOZ")

ele <- raster("2_Reference Data/Elevation/griend_bathy_2015-2020.tif")
SA <- read.csv("2_Reference Data/Water/spacebywaterlevel.csv")
griend_poly <- read_sf("2_Reference Data/Griend/Griend_2015-2020Bathy_100NAP.gpkg")
meta <- read.csv("2_Reference Data/Bird Metadata/2018-2022_islandica.csv") %>%
            dplyr::select(release_ts, tag, captivity, age)

waterdn <- read.csv("2_Reference Data/Water/2017-2022_waterlevel_under50_propdaynight.csv") %>% # proportions of low tide at each part of solar cycle. 
    mutate(total_light = prop_day + prop_civil, 
           total_dark = prop_naut + prop_astro + prop_night, 
           tod_summary = case_when(total_light >= 0.70 ~ "day", #tide majority in day
                                   total_light <0.70 & total_dark <0.70 ~ "twilight", 
                                   total_dark >=0.70 ~ "night")) #tide majority at night

tides <- read.csv("2_Reference Data/Water/Griend/Wide/HLH/allYears-tidalPattern-HLH-griend_Wide-UTC.csv")%>%
    mutate(year = substring(tideID,1,4),
           yearly_id = as.numeric(substring(tideID, 6,8)),
           tideID = as.numeric(sub("_", "", tideID)),
           high_start_time = as.POSIXct(high_start_time),
           date = as.Date(high_start_time, tz="UTC"),
           month = month(high_start_time)) %>%     #allert saved dates on the files in CET (I think) so need to change this.  
    left_join(SA, by = c("low_level"="waterlevel")) %>% #space available in the tide
    left_join(waterdn, by="tideID") %>% 
    filter(month>=9, tod_summary!="twilight" & start_tod == tod_summary) %>% #choose only after September and when start time of day is the same as majority (not twilight)
    rename(tide_tod = start_tod)

#choose residence patches that are not roosting (<50NAP), adults only, not our birds that were in captivity and not 2 other birds that were identified as acting very strange. 
df <- fread("1_WATLAS/Data/res_patches/2018-2022_patches.csv") %>% 
    as.data.frame() %>% 
    filter(!tag  %in% c(911,923), captivity == "no", age==3) %>%   #choose birds for analysis
    mutate(tod_start = ifelse(tod_start %in% c("day", "civil twilight"), "day", "night"),
           datetime_start_UTC = as.POSIXct(datetime_start_UTC, tz = "UTC"),
           month_num = str_pad(lubridate::month(datetime_start_UTC),2,pad="0"),
           month_name = lubridate::month(datetime_start_UTC, label=T),
           month = paste(month_num, month_name, sep="_"),
           tag= str_pad(tag, 4, pad = "0"),
           bird_elevation_start = raster::extract(ele, cbind(x_start, y_start)), # get elevation for where bird actually is to determine whether roosting or not. 
           roost_start = ifelse(bird_elevation_start>=50, 1, 0),
           bird_elevation_end = raster::extract(ele, cbind(x_end, y_end)),
           roost_end = ifelse(bird_elevation_end>=50, 1, 0),
           roost_site = ifelse(roost_start ==1 | roost_end == 1, 1,0)) %>% 
    filter(roost_site == 0,
           tideID %in% tides$tideID) %>%  # only choose tides that started at the same time of day as the majority proportion (night or day)
    left_join(., tides[,c("tideID", "tide_tod")], by = c("tideID"= "tideID")) %>% 
    st_as_sf(coords = c("x_start","y_start"), crs = 32631, remove = F) %>% 
    mutate(dist_from_griendpoly = as.numeric(st_distance(., griend_poly)))

df$geometry <- NULL

# Calculate duration tracked per ind, per tide. 
bird_tides <- df %>%
    group_by(year, tideID, tag, tide_tod) %>%
    summarise(sum_duration_min = sum(duration_s)/60, #low tide durations, we calculate the sum of the residence patch durations in case they leave the area in the middle of the tide and return at the end.
              furthest_patch = max(dist_from_griendpoly)) %>% #calculate furthest distance from griend 
    filter(sum_duration_min >=240) %>%  # which birds on which tides should be analysed. 
    left_join(tides[,c("tideID","low_level")], by = "tideID") #add low level for analysis on distance from griend. 

min_tides <- 5 ########## choose min tides
tides_table <- as.data.frame.matrix(table(bird_tides$tag,bird_tides$tide_tod)) %>% #which birds have at least 5 day tides and 5 night tides?
    filter(day>= min_tides, night >= min_tides)

bird_tides <- bird_tides %>% 
    filter(tag %in% row.names(tides_table)) %>%  #choose only birds that had minimum tides tracked for min durations. 
    mutate(year = as.factor(year),
           tag = as.factor(tag))

#sample sizes 

table(bird_tides$year,bird_tides$tide_tod) #total tides tracked in day/night
tide_sample_size <- unique(as.data.frame(bird_tides[, c("tag", "tideID", "tide_tod")])) #number of tides studied in each year
tide_sample_size <- as.data.frame.matrix(table(tide_sample_size$tag, tide_sample_size$tide_tod))
tide_sample_size$dif <- tide_sample_size$day - tide_sample_size$night
summary(tide_sample_size$day) # how many day tides were birds tracked for
summary(tide_sample_size$night) # how many night tides were they tracked for
summary(tide_sample_size$dif) #what were the differences within individuals for tracking day or night tides. 

sample_size <- unique(as.data.frame(bird_tides[, c("year", "tag")]))
table(sample_size$year)

#------------------------------------------------------------------------------------------------------------------------------------------------------
#### Question 1: Do knots stay longer in residence patches at night? #####

model_dur <- lmer(data=df, log_duration_min ~ tod_start*inout + (1|year) + (1|tag) + (1|tideID))
plot(model_dur)
#anova(model_dur, type="III")
summary(model_dur)
plot(effect("tod_start", model_dur))
plot(effect("tod_start*inout", model_dur))

duration_df <- df[,c("tod_start", "tag", "year", "tideID","inout")]
duration_df$log_duration_min <- predict(model_dur, newdata=duration_df,type="response")
duration_df$duration_min <- 10^duration_df$log_duration_min
duration_df[duration_df$inout == "incoming",]$inout <- "flood"
duration_df[duration_df$inout == "outgoing",]$inout <- "ebb"

p1 <- ggplot(data= duration_df)+
  geom_split_violin(aes(y = duration_min, x = inout, fill=tod_start))+
  scale_fill_manual(values=c("orange","grey15"))+
  scale_y_continuous(expand=c(0,0), limits = c(0,80), name = "duration (mins)")+
  scale_x_discrete(name = "tidal direction")+
  theme_classic()+
  theme(legend.position="none",
        text = element_text(size=22)) #poster size
p1

#----------------------------------------------------------------------------------------------------------------------------------------------------------


##### Question 2: Do red knots stay closer to Griend at night ####
# If birds are more site faithful at night, are they going further out during the day? Base this on elevation or distance to the centre of griend????

bird_tides$low_level_scaled <- scale(bird_tides$low_level, center = T, scale = F) # center subtracts the mean, scale subtracts SD. Do this to make the intercept meaningful (the mean low tide) for the table
scale_c <- attributes(bird_tides$low_level_scaled)$`scaled:center`

model_furth <- lmer(data=bird_tides, furthest_patch ~ tide_tod + low_level_scaled + (1|year) + (1|tag)) 
anova(model_furth, type="III")
summary(model_furth)
coef(model_furth)

intercept_day <- mean(c(coef(model_furth)$tag[,1],coef(model_furth)$year[,1])) - scale_c*coef(model_furth)$tag[1,"low_level_scaled"] #calculate intercept at mean low level rather than 0
intercept_night <- intercept_day + coef(model_furth)$tag[1,"tide_todnight"]
confint_day <- quantile(c(coef(model_furth)$tag[,1],coef(model_furth)$year[,1]), c(0.025, 0.975))

slope <- coef(model_furth)$tag[1,"low_level_scaled"]

# Predict
furth_df <- unique(as.data.frame(bird_tides[,c("tide_tod", "tag", "year", "tideID","low_level_scaled")]
                    ))
furth_df$furthest_patch <- predict(model_furth, newdata=furth_df,type="response")
furth_df$low_level <- furth_df$low_level_scaled + mean(bird_tides$low_level)

p2 <- ggplot(data = bird_tides, aes(x=low_level, y = furthest_patch/1000, col = tide_tod))+
    geom_point(alpha = 0.9, pch = 4)+
    geom_abline(intercept = intercept_day/1000, slope = slope/1000, col = "orange", linewidth = 1)+
    geom_abline(intercept = intercept_night/1000, slope = slope/1000, col = "grey15", linewidth = 1)+
    scale_fill_manual(values=c("orange","grey15"))+
    scale_color_manual(values=c("orange","grey15"))+
    scale_y_continuous(expand=c(0,0), limits = c(0,6.2), name = "furthest patch from Griend (km)")+
    scale_x_continuous(expand=c(0,0), limits = c(-145,25), name = "lowest tide level (NAP)")+
    theme_classic()+
  theme(legend.position="none",
        text = element_text(size=22)) #poster size
p2

#--------------------------------------------------------------------------------------------------------------------------------------------------#


#### Merge polygons into one multipolygon per bird, per tide #### ---------------------------------------------------------------------------------------
# NO NEED TO RUN THIS - DATA IS SAVED  as 2018-2022_bird_tides_multipolygons_min5tides.gpkg
# d_df <- NULL
# 
# for(i in 1:nrow(bird_tides)){
#     
#     id <- bird_tides$tag[i]
#     tide_id <- bird_tides$tideID[i]
#     tod <- bird_tides$tide_tod[i]
#     year <- as.character(bird_tides$year[i])
# 
#     print(paste0("reading in bird...", id))
#     
#     filename <- str_subset(string = get(paste("files", year, sep = "_")), pattern = paste0("//zeus/cos/birds/bijleveld/fieldwork/WATLAS/residence_patches/", year, "/islandica/batch1/underlying_data/islandica-tag_",id, "-tide_", tide_id)) #files for this bird
# 
# # import polygons
#         
#         data <- readRDS(filename)$residence_patches %>% 
#             filter(patch %in% df[c(df$tag==id & df$tideID ==tide_id), "patch"]) %>%  # select only non-roosting patches
#                 mutate(duration_min = (time_end-time_start)/60)%>% 
#                 dplyr::select(duration_min, polygons) %>% 
#                 mutate(polygons = st_sfc(unlist(polygons, recursive = F))) %>% #not sure why we have to unlist this first. this was veeeeeeeery annoying to figure out as it would not recognise as an sf column. 
#                 st_set_geometry("polygons") %>% 
#                 st_set_crs(32631)
#         
#         # get polygons to merge per bird per tide and create data frame (sf object with multipolygon per row).  
# 
#             d <- st_union(data)
#             duration_tracked <- sum(data$duration_min)
#             area_d <- as.numeric(st_area(d))/10^6 #square km
#             d_df_new <- data.frame(tideID = tide_id, id = id, duration_tracked = duration_tracked,  tide_tod = tod, area_km2 = area_d, polygon = d)
#             d_df <- rbind(d_df, d_df_new)
# }
# 
# d_df2 <- merge(tides[,c("tideID","date")], d_df)
# d_df2$year <- year(d_df2$date)
# #st_write(d_df2,"F:/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/2018-2022_bird_tides_multipolygons_min5tides.gpkg")   



########################################################################################

# plot_grid(`2018_day`, `2019_night`)
# 
# 
# 
# do.call(grid.arrange,p)
#         ras_df <- data.frame(ras_sub) %>%
#             filter(!is.na(ID)) %>%
#             rename(sum = "ID")
# 
#         sum_df <- data.frame(year = first(sub$year),
#                                  id = bird_id,
#                                  tide_tod = tod,
#                                  n_tides = nrow(sub),
#                                  visited1_cells = length(ras_df[ras_df$sum>=1,"sum"]),
#                                  visited2_cells = length(ras_df[ras_df$sum>=2,"sum"]),
#                                  visited3_cells = length(ras_df[ras_df$sum>=3,"sum"]),
#                                  visited4_cells = length(ras_df[ras_df$sum>=4,"sum"]),
#                                  visited5_cells = length(ras_df[ras_df$sum>=5,"sum"]))
# 
# 
#         summary_df <- rbind(summary_df, sum_df)
# 
#         if(tod == "day"){  r_day[[bird_id]] <- ras_sub  }
#         if(tod == "night"){  r_night[[bird_id]] <- ras_sub  }
#     }
# }


#merge all night locs and all day locs for each bird
#--------------------------------------------------------------------------------------------------------------------------------------------------#
#### Question 3: Do knots use smaller areas at night?
b_tides <- st_read("F:/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/2018-2022_bird_tides_multipolygons_min5tides.gpkg") %>% 
  mutate(bird_tides = paste(tideID, id, sep="_"))

area_ptide <- b_tides %>% 
    arrange(id) %>% 
    mutate(area_ptide = (area_km2/duration_tracked)*6*60) # area covered in a 6h low tide period. 

area_ptide$geom <- NULL

model_area <- lmer(data=area_ptide, area_ptide ~ tide_tod + (1|year) + (1|id))
anova(model_area, type="III")
summary(model_area)
plot(effect("tide_tod", model_area))


area_df <- area_ptide[,c("tide_tod", "id", "year")]
area_df$area_pred <- predict(model_area, newdata=area_df,type="response")

#mean predicted area for night
mean(area_df[area_df$tide_tod=="night",]$area_pred)
sd(area_df[area_df$tide_tod=="night",]$area_pred)
#mean predicted area during the day
mean(area_df[area_df$tide_tod=="day",]$area_pred)
sd(area_df[area_df$tide_tod=="day",]$area_pred)


#raw data
# p3 <- ggplot(data = area_ptide, aes(x=tide_tod, y = area_ptide, fill = tide_tod))+
#     geom_boxplot(notch=T)+
#     scale_fill_manual(values = c("orange","grey15"))+
#     scale_x_discrete(name = "time of day")+
#     scale_y_continuous(limits = c(0,1.05), expand = c(0,0), name = "area covered per tide (km2)")+
#     theme_classic()+
#     theme(legend.position="none")

#predicted data
p3 <- ggplot(data = area_df, aes(x=tide_tod, y = area_pred, fill = tide_tod))+
    #geom_point(data = area_ptide, aes(x=tide_tod, y = area_ptide, fill = tide_tod), position = position_jitter(0.25), alpha = 0.1)+
    geom_boxplot(notch=T)+
    scale_fill_manual(values = c("orange","grey15"))+
    scale_x_discrete(name = "time")+
    scale_y_continuous(limits = c(0,0.55), expand = c(0,0), name = expression(paste("area covered per tide (km"^~~2, ")")))+
    theme_classic()+
  theme(legend.position="none",
        text = element_text(size=22)) #poster size


p3 


#--------------------------------------------------------------------------------------------------------------------------------------------------#
#### Question 4: Do day and night areas overlap? ####

# merge polygons for each bird (overlap) and calculate area used per tide

overlap <- st_read("F:/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/2018-2022_bird_tides_multipolygons_min5tides.gpkg") %>% 
    group_by(year, id, tide_tod) %>% 
    summarise(geom= st_union(geom),
              duration = sum(duration_tracked),
              num_tides = n()) %>% 
    ungroup() %>% 
    arrange(id) %>% 
    tidyr::pivot_wider(names_from = "tide_tod", values_from = c("duration", "num_tides", "geom")) %>% 
    st_set_geometry("geom_day")  %>% 
    group_by(id) %>% 
    mutate(intersection = st_union(st_intersection(geom_day, geom_night)),
           area_day = as.numeric(st_area(geom_day))/10^6,
           area_night = as.numeric(st_area(geom_night))/10^6,
           area_intersect = as.numeric(st_area(intersection))/10^6, 
           pc_overlap_day = as.numeric((area_intersect/area_day)*100), 
           pc_overlap_night = as.numeric((area_intersect/area_night)*100), 
           geom_day = NULL,
           geom_night = NULL, 
           intersection = NULL) %>% 
    as.data.frame()

#mean total space use (night)
mean(overlap$area_night)
sd(overlap$area_night)
#summary(overlap$area_night)
mean(overlap$num_tides_night)
sd(overlap$num_tides_night)

#mean total space use (day)
mean(overlap$area_day)
sd(overlap$area_day)
#summary(overlap$area_day)
mean(overlap$num_tides_day)
sd(overlap$num_tides_day)

mean(overlap$pc_overlap_day)
sd(overlap$pc_overlap_day)
sd(overlap$pc_overlap_night)

#plot for individuals
overlap_long <- overlap %>% 
    tidyr::pivot_longer(cols = c("pc_overlap_night","pc_overlap_day"), 
                        names_to = "time",
                        values_to = "nonoverlap")



p4 <- ggplot(overlap_long, aes(x = nonoverlap, fill = time))+
    geom_density(alpha = 0.6)+
    scale_fill_manual(values=c("orange","grey15"))+
    scale_x_continuous(limits = c(0,100), name = "overlap (%)")+
    scale_y_continuous(expand=c(0,0), limits = c(0,0.065))+
    theme_classic()+
  theme(legend.position="none",
        text = element_text(size=22)) #poster size

p4

#UNIQUE AREAS ONLY 

#--------------------------------------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------------------------------------#

### NO NEED TO RUN THIS - DATA IS SAVED AS 

### Site fidelity

# The below code creates the raster per individual for each of the selected tides (start in day and 70% of day, or night equivalent)
# ONE RASTER PER INDIVIDUAL OVER ALL TIDES

# polygons <- st_read("F:/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/2018-2022_bird_tides_multipolygons_min5tides.gpkg") 
# ras <- read_sf("2_Reference Data/Wadden Regions/griendonly_bufferapx1km.gpkg")  %>% 
#     st_rasterize(dx=75,dy=75, inside = T)
# template <- ras-1
# 
# summary_df <- data.frame()
# r_day = list()
# r_night = list()
# 
# for(bird_id in unique(polygons$id)){
#     
#     for(tod in c("day", "night")){
#         sub <- polygons %>% 
#             filter(id == bird_id, tide_tod ==tod)
#         ras_sub <- st_rasterize(sub[1,"geom"], template, options = "MERGE_ALG=ADD")
#         
#             for(row in 2:nrow(sub)){
#                 ras_sub2 <- st_rasterize(sub[row,"geom"], template, options = "MERGE_ALG=ADD")
#                 ras_sub <- ras_sub + ras_sub2
#             }
#         
#         ras_df <- data.frame(ras_sub) %>% 
#             filter(!is.na(ID)) %>% 
#             rename(sum = "ID")
#         
#         sum_df <- data.frame(year = first(sub$year),
#                                  id = bird_id, 
#                                  tide_tod = tod,
#                                  n_tides = nrow(sub),
#                                  visited1_cells = length(ras_df[ras_df$sum>=1,"sum"]),
#                                  visited2_cells = length(ras_df[ras_df$sum>=2,"sum"]), 
#                                  visited3_cells = length(ras_df[ras_df$sum>=3,"sum"]),
#                                  visited4_cells = length(ras_df[ras_df$sum>=4,"sum"]),
#                                  visited5_cells = length(ras_df[ras_df$sum>=5,"sum"]))
#         
#         
#         summary_df <- rbind(summary_df, sum_df)
#         
#         if(tod == "day"){  r_day[[bird_id]] <- ras_sub  }
#         if(tod == "night"){  r_night[[bird_id]] <- ras_sub  }
#     }
# }
#
#write.csv(summary_df, "F:/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/2018-2022_postSept_SiteFidelityData.csv")

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#### Question 5: Are birds more site faithful during the day or night? ####

summary_df <- read.csv("F:/OneDrive - Liverpool John Moores University/Synced files/Projects/1_NIOZ projects/2_DayNight/data/2018-2022_postSept_SiteFidelityData.csv")
m_revisits1 <- glmer(cbind(visited2_cells, visited1_cells - visited2_cells) ~ tide_tod + offset(log(n_tides)) + (1 | year) + (1|id), data=summary_df, family = "binomial") 
summary(m_revisits1)
plot(m_revisits1)

plot(effect("tide_tod", m_revisits1))

newdf1 <- data.frame(revisits = 1, 
                     summary_df[,c("tide_tod", "id", "year")],
                     n_tides = mean(summary_df$n_tides))
newdf1$revisits_prob <- predict(m_revisits1,newdata=newdf1,type="response")

summary_df$revisits_pc <- summary_df$visited2_cells/summary_df$visited1_cells

p5 <- ggplot()+
    geom_boxplot(data= newdf1, aes(x= tide_tod, y = revisits_prob, fill = tide_tod), col="black", notch=T)+
    #geom_point(data= summary_df, aes(x= tide_tod, y = revisits_pc, fill = tide_tod), alpha=0.9, col="black",position = position_jitter(0.3))+
    scale_fill_manual(values=c("orange","grey15"))+
    scale_y_continuous(limits = c(0,0.5), expand = c(0,0), name = "revisit probability")+
    scale_x_discrete(name = "time")+
    theme_classic()+
    theme(legend.position="none",
          text = element_text(size=22)) #poster size

p5

legend <- get_legend(p1 + theme(legend.position="top", legend.title = element_blank()))


r1 <- plot_grid(p1,p2, rel_widths = c(1,2), labels = c("a","b"), label_size=20)
r2<- plot_grid(p3,p4,p5, labels=c("c","d", "e"), nrow=1, label_size=20)
plot_all<- plot_grid(legend, r1, r2, nrow=3, rel_heights = c(0.1,1,1))

plot_all
#ggsave("F:\\OneDrive - Liverpool John Moores University\\Synced files\\Projects\\1_NIOZ projects\\2_DayNight\\figures\\fig_2.png", plot_all,units="cm", height=30, width =40, dpi = 600)

#### For supplementary ####
#different thresholds for site fidelity

m_revisits2 <- glmer(cbind(visited3_cells, visited1_cells - visited3_cells) ~ tide_tod + offset(log(n_tides)) + (1 | year/id), data=summary_df, family = "binomial")  
m_revisits3 <- glmer(cbind(visited4_cells, visited1_cells - visited4_cells) ~ tide_tod + offset(log(n_tides)) + (1 | year/id), data=summary_df, family = "binomial")
m_revisits4 <- glmer(cbind(visited5_cells, visited1_cells - visited5_cells) ~ tide_tod + offset(log(n_tides)) + (1 | year/id), data=summary_df, family = "binomial")  



newdf2 <- data.frame(revisits = 2, 
                     summary_df[,c("tide_tod", "id", "year")],
                     n_tides = mean(summary_df$n_tides))
newdf2$revisits_log <- predict(m_revisits2,newdata=newdf2,type="response")

newdf3 <- data.frame(revisits = 3, 
                     summary_df[,c("tide_tod", "id", "year")],
                     n_tides = mean(summary_df$n_tides))
newdf3$revisits_log <- predict(m_revisits3,newdata=newdf3,type="response")

newdf4 <- data.frame(revisits = 4, 
                     summary_df[,c("tide_tod", "id", "year")],
                     n_tides = mean(summary_df$n_tides))
newdf4$revisits_log <- predict(m_revisits4,newdata=newdf4,type="response")


newdf <- rbind(newdf1, newdf2, newdf3,newdf4)


realdf <- rbind(data.frame(revisits = 1, summary_df[,c("tide_tod", "id", "year")],  revisits_pc = summary_df[,"visited2_cells"]/summary_df[,"visited1_cells"]),
                data.frame(revisits = 2,summary_df[,c("tide_tod", "id", "year")],  revisits_pc = summary_df[,"visited3_cells"]/summary_df[,"visited1_cells"]),
                data.frame(revisits = 3, summary_df[,c("tide_tod", "id", "year")],  revisits_pc = summary_df[,"visited4_cells"]/summary_df[,"visited1_cells"]),
                data.frame(revisits = 4, summary_df[,c("tide_tod", "id", "year")],  revisits_pc = summary_df[,"visited5_cells"]/summary_df[,"visited1_cells"]))

newdf2 <- merge(newdf, realdf, by = c("revisits", "tide_tod", "id", "year"))
newdf2$revisits = factor(newdf$revisits)



ggplot()+
  #  geom_split_violin(data= newdf2, aes(x= revisits, y = revisits_log, fill = tide_tod, col = tide_tod))+
    geom_boxplot(data= newdf1, aes(x= tide_tod, y = revisits_prob, fill = tide_tod), col="black", notch=T)+
    #geom_point(data= newdf2, aes(x= revisits, y = revisits_log, fill = tide_tod), pch = 19, alpha = 0.9, position = position_jitterdodge(dodge.width = 0.75))+
    scale_fill_manual(values=c("orange","grey15"))+
    scale_color_manual(values=c("black", "black"))+
    theme_bw()#+
   # facet_grid(cols = vars(year))

#--------------------------------------------------------------------------------------------------------------------------------------------------#

