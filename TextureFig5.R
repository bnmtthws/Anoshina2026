# Mosquito texture paper – Figure 5 pipeline updated 2026.03.18
# streamlined version of TextureFig5.R
# individual panels end up as ggplot2 objects
# Fig5.B
# Fig5.C.texture
# Fig5.C.velocity
# Fig5.D
# Fig5.E
# Fig5.F

### data files can be downloaded along with code here: https://github.com/bnmtthws/Anoshina2026


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(ggridges)
library(Hmisc)

############################
# LOAD DATA
############################

allSessionsDF  <- fread("smoothedTrackedEggDFWithVelocity.csv")
bufferedEggsDF <- fread("bufferedEggsLaidRaw.csv")
bufferedEggsDF2 <- fread("bufferedEggsLaid2.csv")
midline <- fread("2choiceTextureMidlineCoords-mar2024.csv")

setnames(allSessionsDF,"V1","rowID")

############################
# ORDER MOSQUITOES BY FIRST EGG
############################

firstEggTimeDF <- allSessionsDF %>%
  filter(eggLaidEvent == 1) %>%
  group_by(session) %>%
  summarise(firstEggTime = min(timeSeconds), .groups="drop") %>%
  arrange(firstEggTime) %>%
  mutate(mosquitoID = row_number())

allMosquitoesDF <- allSessionsDF %>%
  left_join(firstEggTimeDF[,c("session","mosquitoID")], by="session") %>%
  mutate(timeMinutes = timeSeconds/60) %>%
  arrange(mosquitoID,rowID)

############################
# MIDLINE COORDINATES
############################

pixelsPerMM <- 32.47597

midline <- midline %>%
  mutate(
    xCoord = xCoord/pixelsPerMM,
    yCoord = yCoord/pixelsPerMM
  )

midline <- midline %>%
  group_by(session) %>%
  mutate(topOrBottom = c("top","bottom")) %>%
  ungroup()

############################
# FIG 5B EXAMPLE SESSION
############################

sessionID <- 3

midline_s <- midline %>% filter(session==sessionID)
eggs_s <- allMosquitoesDF %>%
  filter(session==sessionID, eggLaidEvent==1)

top <- midline_s %>% filter(topOrBottom=="top")
bottom <- midline_s %>% filter(topOrBottom=="bottom")

rotationAngle <- atan(
  (bottom$xCoord-top$xCoord)/
    (top$yCoord-bottom$yCoord)
)

Fig5.B <- ggplot(eggs_s,
                 aes(xCoordMeanK21,yCoordMeanK21)) +
  geom_point(aes(colour=timeMinutes),
             size=4,alpha=0.6) +
  geom_path(aes(colour=timeMinutes),
            alpha=0.15) +
  geom_line(data=midline_s,
            aes(xCoord,yCoord),
            linetype="dotted",
            linewidth=1) +
  scale_color_viridis(option="cividis") +
  coord_fixed() +
  xlim(0,65)+ylim(0,65)+
  theme_classic()

############################
# SURFACE OCCUPANCY SUMMARY
############################

surfaceTimingSummaryDF <- allMosquitoesDF %>%
  group_by(mosquitoID) %>%
  summarise(
    numberSmoothFrames =
      sum(agaroseSurface=="smooth"),
    numberTexturedFrames =
      sum(agaroseSurface=="textured"),
    totalFrames = n(),
    
    numberSmoothEggs =
      sum(agaroseSurface=="smooth" &
            eggLaidEvent==1),
    
    numberTexturedEggs =
      sum(agaroseSurface=="textured" &
            eggLaidEvent==1),
    
    .groups="drop"
  ) %>%
  mutate(
    proportionSmoothFrames =
      numberSmoothFrames/totalFrames,
    proportionTexturedFrames =
      numberTexturedFrames/totalFrames,
    totalEggs =
      numberSmoothEggs+numberTexturedEggs
  )

surfaceTimingSummaryDFLong <-
  surfaceTimingSummaryDF %>%
  select(mosquitoID,
         proportionSmoothFrames,
         proportionTexturedFrames) %>%
  pivot_longer(
    -mosquitoID,
    names_to="smoothVsTextured",
    values_to="proportionOfFrames"
  ) %>%
  mutate(
    smoothVsTextured =
      recode(smoothVsTextured,
             proportionSmoothFrames="smooth",
             proportionTexturedFrames="textured")
  )

############################
# FIG 5D
############################

Fig5.D <- ggplot(
  surfaceTimingSummaryDFLong,
  aes(x = smoothVsTextured,
      y = proportionOfFrames,
      group = mosquitoID)
) +
  
  geom_line(alpha = 0.3, colour = "black") +
  
  geom_point(
    aes(colour = smoothVsTextured),
    alpha = 0.5,
    size = 5
  ) +
  
  stat_summary(
    aes(
      x = smoothVsTextured,
      y = proportionOfFrames,
      group = smoothVsTextured,      # key fix
      colour = smoothVsTextured
    ),
    fun.data = mean_cl_boot,
    geom = "pointrange",
    linewidth = 2,
    size = 1.35,
    alpha = 0.7,
    position = position_nudge(x = c(-0.25, 0.25))
  ) +
  
  scale_color_manual(
    values = c("darkorange", "darkgreen")
  ) +
  
  scale_y_continuous(
    breaks = seq(0,1,0.25),
    limits = c(0,1)
  ) +
  
  labs(
    x = "",
    y = "proportion of time on agarose substrate"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18)
  )

############################
# PREFERENCE INDICES
############################

# determine first and last egg frame for each mosquito
egg_bounds <- allMosquitoesDF %>%
  filter(eggLaidEvent == 1) %>%
  group_by(mosquitoID) %>%
  summarise(
    firstEggRowID = min(rowID),
    lastEggRowID  = max(rowID),
    .groups = "drop"
  )

# keep only frames during egg laying
framesWhileEggsLaidDF <- allMosquitoesDF %>%
  inner_join(egg_bounds, by = "mosquitoID") %>%
  filter(rowID >= firstEggRowID & rowID <= lastEggRowID)


frame_counts <- framesWhileEggsLaidDF %>%
  group_by(mosquitoID, agaroseSurface) %>%
  summarise(nFrames = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = agaroseSurface,
    values_from = nFrames,
    values_fill = 0
  ) %>%
  rename(
    numberSmoothFrames = smooth,
    numberTexturedFrames = textured
  )



preferenceIndexSummaryDF <- surfaceTimingSummaryDF %>%
  select(mosquitoID,
         numberSmoothEggs,
         numberTexturedEggs) %>%
  left_join(frame_counts, by = "mosquitoID")

calculatePreferenceIndex <- function(treatment, control){
  (treatment - control) / (treatment + control)
}


preferenceIndexSummaryDF <- preferenceIndexSummaryDF %>%
  mutate(
    
    ovipositionPreferenceIndex =
      calculatePreferenceIndex(
        numberTexturedEggs,
        numberSmoothEggs
      ),
    
    agaroseSurfacePreferenceIndex =
      calculatePreferenceIndex(
        numberTexturedFrames,
        numberSmoothFrames
      )
  )

preferenceIndexSummaryDFLong <-
  preferenceIndexSummaryDF %>%
  pivot_longer(
    cols = c(
      ovipositionPreferenceIndex,
      agaroseSurfacePreferenceIndex
    ),
    names_to = "eggsLaidVStimeSpent",
    values_to = "preferenceIndex"
  ) %>%
  mutate(
    eggsLaidVStimeSpent = recode(
      eggsLaidVStimeSpent,
      ovipositionPreferenceIndex = "eggsLaid",
      agaroseSurfacePreferenceIndex = "timeSpent"
    )
  )

############################
# FIG 5E
############################

Fig5.E <- ggplot(
  preferenceIndexSummaryDFLong,
  aes(
    x = eggsLaidVStimeSpent,
    y = preferenceIndex,
    group = mosquitoID
  )
) +
  geom_line(alpha = 0.3, colour = "black") +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    aes(colour = eggsLaidVStimeSpent),
    size = 5,
    alpha = 0.5
  ) +
  
  stat_summary(
    aes(
      colour = eggsLaidVStimeSpent,
      group = eggsLaidVStimeSpent
    ),
    fun.data = mean_cl_boot,
    geom = "pointrange",
    linewidth = 2,
    size = 1.35,
    alpha = 0.5,
    position = position_nudge(x = c(-0.25, 0.25))
  ) +
  
  scale_y_continuous(
    breaks = seq(-1,1,0.25),
    limits = c(-1,1)
  ) +
  
  scale_color_manual(
    values = c("#56B4E9", "#009E73")
  ) +
  
  labs(
    x = "",
    y = "preference index"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0.2,"cm"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18)
  )

############################
# FIG 5F – MOVEMENT RIBBON
############################

framesPerSecond <- 70

rollingFramesHalfSecond <- round(0.5 * framesPerSecond)

############################
# FIG 5F – velocity distributions by surface
############################

movingVelocityDF <- allMosquitoesDF %>%
  arrange(mosquitoID, frames) %>%
  group_by(mosquitoID) %>%
  mutate(
    rollingMeanVelocityHalfSecond = zoo::rollmean(
      velocityMMPerSecond,
      k = rollingFramesHalfSecond,
      fill = NA,
      align = "center"
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(rollingMeanVelocityHalfSecond)) %>%
  mutate(
    logVelocity = log10(rollingMeanVelocityHalfSecond+.01),
    surface = factor(agaroseSurface,
                     levels = c("smooth","textured"))
  )

Fig5.F <- ggplot(movingVelocityDF,
                 aes(x = logVelocity,
                     fill = surface,
                     colour = surface)) +
  
  geom_density(alpha = 0.35,
               size = 1) +
  
  scale_fill_manual(values = c(
    smooth = "darkorange",
    textured = "darkgreen"
  )) +
  
  scale_colour_manual(values = c(
    smooth = "darkorange",
    textured = "darkgreen"
  )) +
  
  coord_cartesian(xlim = c(log10(0.001), log10(15))) +
  
  labs(
    x = "log10 (rolling average velocity mm/s)",
    y = "density"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 18)
  )

############################
# FIG 5G – VELOCITY AROUND EGG
############################

eggEvents <- bufferedEggsDF %>%
  filter(focalEggEvent==1)

smoothEggEvents <-
  bufferedEggsDF %>%
  filter(eggID %in%
           eggEvents$eggID[
             eggEvents$agaroseSurface=="smooth"])

texturedEggEvents <-
  bufferedEggsDF %>%
  filter(eggID %in%
           eggEvents$eggID[
             eggEvents$agaroseSurface=="textured"])

Fig5.G <- ggplot()+
  geom_smooth(
    data=texturedEggEvents,
    aes(timeRelativeToEgg,
        velocityMMPerSecond,
        colour="textured"))+
  geom_smooth(
    data=smoothEggEvents,
    aes(timeRelativeToEgg,
        velocityMMPerSecond,
        colour="smooth"))+
  geom_vline(xintercept=0,
             linetype="dashed")+
  scale_color_manual(
    values=c("darkorange","darkgreen")
  )+
  coord_cartesian(ylim=c(0,1.6))+
  theme_classic()


############################
# FIG 5C – ETHOGRAM
############################

numberSessions <- n_distinct(allMosquitoesDF$mosquitoID)
orderMosquitoIDs <- sort(unique(allMosquitoesDF$mosquitoID), decreasing = TRUE)

# Use the full rolling-velocity dataframe here so stationary periods are still present
rollingVelocityEthogramDF <- movingVelocityDF %>%
  mutate(
    agaroseSurface = factor(
      agaroseSurface,
      levels = c("smooth", "textured")
    )
  )

# Determine stationary threshold from representative stationary samples
stationarySample1 <- rollingVelocityEthogramDF[
  rollingVelocityEthogramDF$mosquitoID == 7 &
    rollingVelocityEthogramDF$frames > 75000 &
    rollingVelocityEthogramDF$frames < 145000,
  "rollingMeanVelocityHalfSecond"
][[1]]

stationarySample2 <- rollingVelocityEthogramDF[
  rollingVelocityEthogramDF$mosquitoID == 7 &
    rollingVelocityEthogramDF$frames > 395000 &
    rollingVelocityEthogramDF$frames < 414000,
  "rollingMeanVelocityHalfSecond"
][[1]]

stationarySample3 <- rollingVelocityEthogramDF[
  rollingVelocityEthogramDF$mosquitoID == 8 &
    rollingVelocityEthogramDF$frames > 20000 &
    rollingVelocityEthogramDF$frames < 22450,
  "rollingMeanVelocityHalfSecond"
][[1]]

fastest100S1 <- head(sort(stationarySample1, decreasing = TRUE), 70)
fastest100S2 <- head(sort(stationarySample2, decreasing = TRUE), 70)
fastest100S3 <- head(sort(stationarySample3, decreasing = TRUE), 70)

meanFastestS1 <- mean(fastest100S1, na.rm = TRUE)
meanFastestS2 <- mean(fastest100S2, na.rm = TRUE)
meanFastestS3 <- mean(fastest100S3, na.rm = TRUE)

meanFastestVelocityWhileStationary <- mean(
  c(meanFastestS1, meanFastestS2, meanFastestS3),
  na.rm = TRUE
)

stationaryVelocityLimit <- meanFastestVelocityWhileStationary

# Complete missing frames for each mosquito so geom_raster has a full regular grid
filteredStationaryVelocityEthogramDF <- rollingVelocityEthogramDF %>%
  arrange(mosquitoID, frames) %>%
  group_by(mosquitoID) %>%
  tidyr::complete(frames = seq(min(frames, na.rm = TRUE), max(frames, na.rm = TRUE), by = 1)) %>%
  tidyr::fill(session, mosquitoID, agaroseSurface, .direction = "downup") %>%
  mutate(
    timeMinutes = frames / framesPerSecond / 60,
    rollingMeanVelocityHalfSecond = dplyr::if_else(
      is.na(rollingMeanVelocityHalfSecond) |
        rollingMeanVelocityHalfSecond < stationaryVelocityLimit,
      0,
      rollingMeanVelocityHalfSecond
    ),
    log10RollingMeanVelocityHalfSecond = log10(rollingMeanVelocityHalfSecond+0.01),
    mosquitoFacet = factor(mosquitoID, levels = orderMosquitoIDs),
    rasterY = 1
  ) %>%
  ungroup()

ethogramEggEventsDF <- allMosquitoesDF %>%
  filter(eggLaidEvent == 1) %>%
  distinct(mosquitoID, frames, .keep_all = TRUE) %>%
  mutate(
    mosquitoFacet = factor(mosquitoID, levels = orderMosquitoIDs)
  )

frameBreaks <- seq(0, 120, by = 15) * 60 * framesPerSecond
frameLimits <- c(0, 120 * 60 * framesPerSecond)

Fig5.C.texture <- ggplot(
  filteredStationaryVelocityEthogramDF,
  aes(x = frames, y = rasterY)
) +
  geom_raster(aes(fill = agaroseSurface)) +
  scale_fill_manual(
    values = c(
      smooth = "darkorange",
      textured = "darkgreen"
    )
  ) +
  facet_wrap(
    vars(mosquitoFacet),
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = frameBreaks,
    labels = seq(0, 120, by = 15),
    limits = frameLimits
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    line = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

Fig5.C.velocity <- ggplot(
  filteredStationaryVelocityEthogramDF,
  aes(x = frames, y = rasterY)
) +
  geom_raster(aes(fill = log10RollingMeanVelocityHalfSecond),alpha=1) +
  geom_point(
    data = ethogramEggEventsDF,
    aes(
      x = frames,
      y = 1   # same as rasterY
    ),
    shape = 10,          # "X"
    colour = "white",
    stroke = 1.2,       # thickness of the X
    size = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_fill_viridis_c(option = "magma") +
  facet_wrap(
    vars(mosquitoFacet),
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = frameBreaks,
    labels = seq(0, 120, by = 15),
    limits = frameLimits
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "time (minutes)",
    y = NULL,
    fill = "log10 rolling average velocity\n(mm/s, 0.5 s window)"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

Fig5.C <- Fig5.C.texture / Fig5.C.velocity +
  plot_layout(heights = c(1, 6))



#### Fig5.F #####
Fig5.F2 <- ggplot(movingVelocityDF,
                 aes(x = logVelocity,
                     fill = surface,
                     colour = surface)) +
  
  geom_density(alpha = 0.35,
               size = 1) +
  
  scale_fill_manual(values = c(
    smooth = "darkorange",
    textured = "darkgreen"
  )) +
  
  scale_colour_manual(values = c(
    smooth = "darkorange",
    textured = "darkgreen"
  )) +
  
  coord_cartesian(xlim = c(log10(stationaryVelocityLimit - 0.1), log10(15)),ylim=c(0,.8)) +
  
  labs(
    x = "log10 (rolling average velocity mm/s)",
    y = "density"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 14),a
    axis.title = element_text(size = 18)
  ) + 
  
  geom_vline(xintercept=log10(stationaryVelocityLimit+0.01),linetype = "dashed")

