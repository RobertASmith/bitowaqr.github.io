##########################
#     
### TITLE
#
#
#########################

time1 = Sys.time()

## Load the requires packages, functions, and data
    source("./src/geo_functions.R")
    install_n_load(c("ggplot2","cowplot","shiny","leaflet","sp","raster","rgeos","geosphere"))

# LOAD DATA
    # LOAD LSOA with population and deprivation info
    poly_df = raster::shapefile("./input/polygons")
    
    # England polygon
    england_poly <- maptools::unionSpatialPolygons(poly_df, IDs = rep("England",times=dim(poly_df)[1]))
    
    # LOAD PARKRUN EVENTS
      parkrun_marker = raster::shapefile("./input/marker_england")
      # in_england = gIntersects(england_poly,parkrun_marker, byid = T)
      # parkrun_marker = subset(parkrun_marker,in_england[,1])
      # false_data = parkrun_marker$Club %in% c("newboroughforest" ,"standrews") # there is an error in the locations, for some rason
      # parkrun_marker = subset(parkrun_marker,!false_data)
      # compute event - lsoa data
      parkrun_marker = get_event_res(event = parkrun_marker,
                                     poly_pop = poly_df,
                                     pop_n = poly_df$pop,
                                     objective = "dist^2 *pop")
      
      
    # LOAD PARKS & GARDENS
      greenspaces = shapefile("./input/greenspaces")
      # 1. select subset of parks: large enough
      greenspaces$area = raster::area(greenspaces)
      greenspaces$area = greenspaces$area/1000^2 # in square km
      enough_space_to_run = greenspaces$area >= 0.05
      greenspaces = subset(greenspaces,enough_space_to_run)
      # 2. select subset of parks: only English
      in_england = gIntersects(england_poly,greenspaces, byid = T)
      greenspaces = subset(greenspaces,in_england[,1])
      rm("in_england","enough_space_to_run")


### FIRST PCIKS: FIND BEST CANDIDATES FOR A SINGLE NEW EVENT
    candidates_1st = greenspaces
    first.picks = distanc0r2(candidates = candidates_1st,
                             pop_poly = poly_df,
                             pop_n = poly_df$pop,
                             objective =  "dist^2 * pop", 
                             events = parkrun_marker,
                             showstatus = F)
    
    candidates_1st$d_sum_raw = with(first.picks,round((candidates.res$sum_raw - baseline.res$sum_0_raw)/baseline.res$sum_0_raw,5)*100)
    candidates_1st$d_mean_raw = with(first.picks,round((candidates.res$mean_raw - baseline.res$mean_0_raw)/baseline.res$mean_0_raw,5)*100)
    candidates_1st$d_mdn_raw = with(first.picks,round((candidates.res$median_raw - baseline.res$median_0_raw)/baseline.res$median_0_raw,5)*100)
    candidates_1st$d_sum_obj = with(first.picks,round((candidates.res$sum_obj - baseline.res$sum_0_obj)/baseline.res$sum_0_obj,5)*100)
    candidates_1st$d_mean_obj = with(first.picks,round((candidates.res$mean_obj - baseline.res$mean_0_obj)/baseline.res$mean_0_obj,5)*100)
    candidates_1st$d_mdn_obj = with(first.picks,round((candidates.res$median_obj - baseline.res$median_0_obj)/baseline.res$median_0_obj,5)*100)
    candidates_1st$srvd_lsoa = first.picks$candidates.res$served_lsoa_i
    candidates_1st$pos = rank(candidates_1st$d_sum_obj,ties.method="random") 
    # top 50 first candidates
    top_n = 1:50
    top_first_parks = subset(candidates_1st, candidates_1st$pos %in% top_n)
    top_first_parks$rel_size = (abs(top_first_parks$d_sum_obj)/max(abs(top_first_parks$d_sum_obj))) + 1/4
    top_first_parks.coord = coordinates(top_first_parks)
    top_first_parks$lon = top_first_parks.coord[,1]
    top_first_parks$lat = top_first_parks.coord[,2]
    # save
    shapefile(top_first_parks,filename = "./output/top_first_parks", overwrite =T)
    write.csv(first.picks$candidates.res,file = "./output/first_candidate_res.csv",row.names = F)
    write.csv(first.picks$baseline.res,file = "./output/first_baseline_res.csv",row.names = F)
    rm("candidates_1st","first.picks","top_n","top_first_parks.coord") # free memory
    

### CONSECUTIVE PCIKS: FIND BEST CANDIDATES FOR CONSECUTIVE NEW EVENTS 
    candidates_consecutive = distance_l00per(candidates = greenspaces,
                                             pop_poly = poly_df,
                                             pop_n =  poly_df$pop,
                                             events = parkrun_marker,
                                             top_n = 25,
                                             showstatus = F) #top 25
    
    top_consecutive_parks = greenspaces
    top_consecutive_parks$pos = match(1:length(top_consecutive_parks[,1]),candidates_consecutive$index[-1])
    top_consecutive_parks_select = 1:length(greenspaces[,1]) %in% candidates_consecutive$index
    top_consecutive_parks = subset(top_consecutive_parks,top_consecutive_parks_select)
    top_consecutive_parks = top_consecutive_parks[order(top_consecutive_parks$pos),]
    top_consecutive_parks$objective = candidates_consecutive$objective[-1]/1e+9 # in million
    top_consecutive_parks$change = candidates_consecutive$change[-1]
    top_consecutive_parks$rel_size = (abs(top_consecutive_parks$change)/max(abs(top_consecutive_parks$change))) + 1/4
    top_consecutive_parks.coord = coordinates(top_consecutive_parks)
    top_consecutive_parks$lon = top_consecutive_parks.coord[,1] 
    top_consecutive_parks$lat = top_consecutive_parks.coord[,2]
    # save
    shapefile(top_consecutive_parks,filename = "./output/top_consecutive_parks", overwrite =T)
    write.csv(candidates_consecutive,"./output/candidates_cons_res.csv",row.names = F)
    rm("top_consecutive_parks.coord","top_consecutive_parks_select")


# # #   BUILDING THE MAP   # # #
    # color palette
    q_dists = as.numeric(quantile(poly_df$mn_dstn,probs = c(seq(0,1,by=0.1))))
    pal_dist <- colorBin("RdYlGn", domain = poly_df$mn_dstn, bins = q_dists,reverse =T)
    # pal_dist_new <- colorBin("RdYlGn", domain = poly_df$new_mn_dstn, bins = q_dists,reverse =T)
    
    
    map_distance_squared = 
      leaflet() %>%
      # addTiles(group = "OSM Map") %>%
      addProviderTiles("CartoDB.Positron",group= "Carto Map", 
                       options= providerTileOptions(opacity = 0.99)
      ) %>%
      # # DISTANCES POLYGONS
      addPolygons(data = poly_df, group = "Show baseline distances",
                  color = "gray", smoothFactor = 1,stroke = T,opacity = 0.5,
                  weight = 0.1, fillOpacity = 0.5,fillColor = ~pal_dist(mn_dstn),
                  highlight = highlightOptions(
                    weight = 1,color = "cyan",opacity = 1,bringToFront = FALSE,sendToBack = TRUE),
                  popup = paste(
                    poly_df$name,"<br>",
                    "Nearest event:", poly_df$nrst_vn,"<br>",
                    "Distance: ",poly_df$mn_dstn," km <br>",
                    "IMD score:", poly_df$a,"<br>",
                    "Population:", poly_df$pop)
      ) %>%
      # # add distance legend
      addLegend(group = "Show legend",
                position = "bottomright", 
                pal = pal_dist, values = poly_df$mn_dstn,
                opacity = 0.7, title="Distance deciles",
                labFormat = labelFormat(suffix  = "km") 
      ) %>%
      # # ESTABLISHED EVENTS
      addCircleMarkers(
        group = "Parkrun events",
        data = parkrun_marker,
        lng = ~lon,lat = ~lat,
        radius = 3,fillColor = "blue",
        stroke = FALSE, fillOpacity = 0.9,
        popup = paste("Course:",parkrun_marker$Club,"<br>",
                      "Established:", parkrun_marker$Estblsh,"<br>",
                      "Served pop:", round(parkrun_marker$srvd_pop),"<br>",
                      "Served LSOA:", round(parkrun_marker$srvd_lsoa),"<br>",
                      "Mean distance:", round(parkrun_marker$mean_dist/1000,1),"km <br>",
                      "Objective:", parkrun_marker$objective,"<br>",
                      "Mean participants:", round(parkrun_marker$Mn_prtc),"<br>",
                      "Mean volunteers:", round(parkrun_marker$Mn_vlnt))
      ) %>%
      # # ANY PARK
      addPolygons(
        data = greenspaces, group = "Parks considered",
        smoothFactor = 2,stroke = F, fillOpacity = 0.8,fillColor = "lightgreen",
        highlight = highlightOptions(
          weight = 1,color = "cyan",opacity = 1,bringToFront = T,sendToBack = TRUE),
        popup = paste(greenspaces$distName1,"<br>",
                      "Area:",round(greenspaces$area,2),"km2")
      ) %>%
      # # PARKS FOR FIRST PICK
      addCircleMarkers(
        data = top_first_parks, group = "New candidates (first)", 
        lng = ~lon, lat = ~lat,
        radius = top_first_parks$rel_size*15,fillColor = "pink",stroke = FALSE, fillOpacity = 0.9,
        popup = paste("Park:", top_first_parks$distName1,"<br>",
                      "Pick no:", top_first_parks$pos,"<br>",
                      "% change in objective:",top_first_parks$d_sum_obj,"% <br>",
                      "Served LSOA:", top_first_parks$lsoa_served,"<br>"),
        label = paste(top_first_parks$pos),
        labelOptions = labelOptions(noHide = T, textsize = "9px",textOnly = T, direction ="center") 
      ) %>%
      # # PARKS FOR CONSECUTIVE PICKS
      addCircleMarkers(
        data = top_consecutive_parks, group = "New candidates (consecutive)", 
        lng = ~lon, lat = ~lat,
        radius = top_consecutive_parks$rel_size*15,fillColor = "orange",stroke = FALSE, fillOpacity = 0.9,
        popup = paste("Park:", top_consecutive_parks$distName1,"<br>",
                      "Pick no:", 1:length(top_consecutive_parks$pos),"<br>",
                      "% change in objective:",top_consecutive_parks$change,"% <br>"),
        label = paste(top_consecutive_parks$pos),
        labelOptions = labelOptions(noHide = T, textsize = "9px",textOnly = T, direction ="center")
      ) %>%
      # # LAYER CONTROL
      addLayersControl(
        baseGroups = c("Hide baseline distances","Show baseline distances"),
        overlayGroups = c(
          "Parkrun events",
          "Parks considered",
          "New candidates (first)",
          "New candidates (consecutive)",
          "Show legend"),
        options = layersControlOptions(collapsed = F,autoZIndex=T)
      ) %>% 
      setView(lng = -1.43, lat = 53.36, zoom = 7
      ) %>%
      # # hide all layers except parkrun events
      hideGroup("Show baseline distances") %>%
      hideGroup("New candidates (first)") %>%
      hideGroup("New candidates (consecutive)") %>%
      hideGroup("Parks considered") %>%
      hideGroup("Show legend") 
    

  ## Save map as stand-alone html widget
    # htmlwidgets::saveWidget(map_distance_squared, file = "distance_squared_map.html", selfcontained = F)
    

# take time
  time2 = Sys.time()
  print("Duration:")
  print(time2-time1)