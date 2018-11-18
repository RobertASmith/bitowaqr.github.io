
install_n_load <- function(package){
  for(i in 1:length(package)){
    if(eval(parse(text=paste("require(",package[i],")")))==0) {
      install.packages(package)
    }
  }
  return (eval(parse(text=paste("require(",package,")"))))
}

get_event_res = function(event = parkrun_marker,
                         poly_pop = poly_df,
                         pop_n = poly_df$pop,
                         objective = "dist^2 *pop"
                         ){
    ## Baseline distances and objective
    objective_function = function(dist,pop) {eval(parse(text=objective))}
    baseline_distances = distm(coordinates(event),coordinates(poly_pop))
    # get distance, pop, objective info for events
    dist.M.temp = obj.M.temp = pop.M.temp = 
      matrix(nrow = dim(baseline_distances)[1],
             ncol= dim(baseline_distances)[2],
             data =0)
    for(col in 1:dim(baseline_distances)[2]){
      min.dist.temp = min(baseline_distances[,col])
      row.index = which(baseline_distances[,col] == min.dist.temp)
      pop_i = pop_n[col]
      
      dist.M.temp[row.index,col]= min.dist.temp
      pop.M.temp[row.index,col] = pop_i
      obj.M.temp[row.index,col] = objective_function(dist = min.dist.temp, pop = pop_i)
    }
    
    event$srvd_pop =  rowSums(pop.M.temp)
    event$objective = round(rowSums(obj.M.temp)/1e+9,0) # in millions
    event$srvd_lsoa = apply(dist.M.temp,1,function(x) sum(x != 0))
    event$mean_dist = apply(dist.M.temp,1,function(x) mean(x[x!=0]))
    
    return(event)
    }



# function to find a new park to minimize the objective = distance^2 
distanc0r2 = function(candidates = greenspaces, # new park coorindates matrix
                      pop_poly = poly_df, # lsoa location
                      pop_n = poly_df$pop, # lsoa population
                      objective = "dist^2 * pop", # can be set to x^2 or x < 2
                      events = parkrun_marker, #parkrunevent coord matrix
                      showstatus = F
                      ){ 
      # transform distance (not a real objectve function)
      objective_function = function(dist,pop) {eval(parse(text=objective))}
      
      cat("\n Computing baseline distances")
      # distances between population coord and events
      population.coord = coordinates(pop_poly)
      events.coord = coordinates(events)
      dist_0_raw = distm(population.coord,events.coord)
      
      # for each LSOA there is a min distance to a parkrun event
      dist_0_raw.min = apply(dist_0_raw,1,function(x) min(x))
      dist_0_obj = objective_function(pop = pop_n,dist = dist_0_raw.min)
      
      sum_0_raw = sum(dist_0_raw.min)
      mean_0_raw = mean(dist_0_raw.min)
      median_0_raw = median(dist_0_raw.min)
      sum_0_obj = sum(dist_0_obj)  
      mean_0_obj = mean(dist_0_obj)
      median_0_obj = median(dist_0_obj)
      
      
      baseline_res = data.frame(sum_0_raw=sum_0_raw,
                                mean_0_raw = mean_0_raw,
                                median_0_raw = median_0_raw,
                                sum_0_obj = sum_0_obj,
                                mean_0_obj = mean_0_obj,
                                median_0_obj = median_0_obj)
      
      # distances between population coord and candidate parks
      cat("\n Computing new distances")
      candidates.coord = coordinates(candidates)
      dist_new_raw = distm(population.coord,candidates.coord)
      dist_new_obj = objective_function(pop = pop_n,dist = dist_new_raw)
      
      len = length(candidates[,1])
      candidates.res = data.frame(id = NA, 
                                  lon = NA,
                                  lat = NA,
                                  sum_raw=NA,
                                  mean_raw = NA,
                                  median_raw = NA,
                                  sum_obj = NA,
                                  mean_obj = NA,
                                  median_obj = NA,
                                  served_lsoa_i = rep(NA,times=len)) # define size in advance
      
      # check new total distances for each candidate
      timecheck = 0
      cat("\n Checking changes in total distances")
      for(i in 1:len){
        
        
        if(showstatus == T){
          # Report status
          progress = round((i/len),4)*100
          if(abs(progress - timecheck)>1){
            timecheck = progress
            cat("\n",progress,"%")
          }
        }
        
        id = i
        reassign_population = dist_new_raw[,i] < dist_0_raw.min
        served_lsoa_i = sum(reassign_population)
        
        dist_i_raw = c(dist_new_raw[reassign_population,i] ,dist_0_raw.min[!reassign_population])
        sum_i_raw = sum(dist_i_raw)
        mean_i_raw = mean(dist_i_raw)
        median_i_raw = median(dist_i_raw)
        
        dist_i_obj = c(dist_new_obj[reassign_population,i] ,dist_0_obj[!reassign_population])
        sum_i_obj = sum(dist_i_obj)
        mean_i_obj = mean(dist_i_obj)
        median_i_obj = median(dist_i_obj)
        
        temp.df = data.frame(id = id, 
                             lon = candidates.coord[i,1],
                             lat = candidates.coord[i,2],
                             sum_raw= sum_i_raw,
                             mean_raw = mean_i_raw,
                             median_raw = median_i_raw,
                             sum_obj = sum_i_obj,
                             mean_obj = mean_i_obj,
                             median_obj = median_i_obj,
                             served_lsoa_i = served_lsoa_i)
        
        candidates.res[i,] = temp.df
        
      }
      
      return(list("candidates.res" = candidates.res,
                  "baseline.res" = baseline_res))
      }


### function to find the top consecutive parkrun event locations
distance_l00per = function(candidates = greenspaces, # new park coorindates matrix
                           events = parkrun_marker, #parkrunevent coord matrix
                           pop_poly = poly_df, # lsoa location
                           pop_n = poly_df$pop, # lsoa population
                           objective = "dist^2 * pop", # can be set to x^2 or x < 2
                           showstatus = T,
                           top_n = 3
                           ){ 
      # transform distance (not a real objectve function)
      objective_function = function(dist,pop) {eval(parse(text=objective))}
      
      cat("\n Computing baseline baseline distances")
      # distances between population coord and events
      population.coord = coordinates(pop_poly)
      events.coord = coordinates(events)
      candidates.coord = coordinates(candidates)
      candidate_names = 1:length(candidates.coord[,1])
      
      cat("\n Computing distances between lsoa and parkrun events")
      dist_pop_events = distm(population.coord,events.coord)
      min_dist_pop_events = apply(dist_pop_events,1,function(x) min(x))
      trans_min_dist_pop_events = objective_function(pop = pop_n, dist= min_dist_pop_events)
      objective_state.append = sum(trans_min_dist_pop_events)  
      
      cat("\n Computing distances between lsoa and candidates")
      dist_pop_candidates = distm(population.coord,candidates.coord)
      trans_dist_pop_candidates = objective_function(pop = pop_n, dist= dist_pop_candidates)
      
      save = data.frame(index=NA,objective=NA,change=rep(NA,times=top_n+1))
      save[1,] = c(0,objective_state.append,0)
      
      for(k in 1:top_n){
        timecheck = 0
        len = length(candidate_names)
        for(i in 1:len){
          
          if(showstatus == T){
            # Report status
            progress = round((i/len),4)*100
            if(abs(progress - timecheck)>5){
              timecheck = progress
              cat("\n Outer loop: ",k, "/",top_n,"; inner loop: ",progress,"%", sep="")
            }
          }
          reassign_population = dist_pop_candidates[,i] < min_dist_pop_events
          if(sum(reassign_population)>0){
            new_trans_min_dist_pop_combi = rep(NA,times=length(reassign_population))
            new_trans_min_dist_pop_combi[which(reassign_population)] = trans_dist_pop_candidates[reassign_population,i]
            new_trans_min_dist_pop_combi[-which(reassign_population)] = trans_min_dist_pop_events[!reassign_population]
            prop_objective_state = sum(new_trans_min_dist_pop_combi)
            
            if(prop_objective_state < objective_state.append){
              cat("\n new objective at",i,": ",prop_objective_state - objective_state.append)
              objective_state.append = prop_objective_state
              index = i
              candidate_i = candidate_names[i]
              append.min_dist_pop_event = rep(NA,times=length(reassign_population))
              append.min_dist_pop_event[which(reassign_population)] = dist_pop_candidates[reassign_population,i]
              append.min_dist_pop_event[-which(reassign_population)] = min_dist_pop_events[!reassign_population]
              append.trans_min_dist_pop_event = new_trans_min_dist_pop_combi
            }
          }
          
        }
        
        save[k+1,] = c(candidate_i,
                       objective_state.append,
                       round((objective_state.append-save$objective[k])/save$objective[k],5)*100)
        
        
        min_dist_pop_events = append.min_dist_pop_event
        trans_min_dist_pop_events = append.trans_min_dist_pop_event
        
        dist_pop_candidates = dist_pop_candidates[,-index]
        trans_dist_pop_candidates = trans_dist_pop_candidates[,-index]
        candidate_names = candidate_names[-index]
        
      }
      
      return(save)
      }



# Place to store unused functions
#### SENSITIVITY ANALYSIS : 
# #    Check if results depend on where parks are located.
# #    Using centroids of LSOA to check this.
# 
# # Analysis requires too much memory to use the function used earlier...
# distanc0r_for_heavy_compuation = function(LSOAcandidates = poly_df,
#                                         population =  poly_df,
#                                         objective = "x^2 * population$pop",
#                                         events = parkrun_marker){
#       # transform distance (not a real objectve function)
#       objective_function = function(x) {eval(parse(text=objective))}
#       
#       cat("\n Computing baseline distances")
#       # distances between population coord and events
#       population.coord = coordinates(population)
#       events.coord = coordinates(events)
#       dist_0_raw = distm(population.coord,events.coord)
#       
#       dist_0_raw.min = apply(dist_0_raw,1,function(x) min(x))
#       dist_0_obj = objective_function(dist_0_raw.min)
#       sum_0_obj = sum(dist_0_obj)  
#       
#       # distances between population coord and candidate parks
#       LSOAcandidates.coord = coordinates(LSOAcandidates)
#       
#       len = length(LSOAcandidates.coord[,1])
#       save.df = data.frame(index = 0:len,
#                            objective = NA)
#       save.df$objective[1] = sum_0_obj
#       
#       draw_seq = seq(0,len,length.out = 50)
#       draw_seq = floor(draw_seq)
#       seq_len = length(draw_seq)
#       for(i in 2:seq_len){
#         candidates = (draw_seq[i-1]+1):draw_seq[i]
#         dist_i_raw = distm(population.coord,LSOAcandidates.coord[candidates,]) 
#         dist_obj_j = objective_function(dist_i_raw)
#         len_j = length(candidates)
#         for(j in 1:len_j){
#           cat("\n i:",i-1,"/",seq_len-1,"  - j: ",j,"/",len_j,sep="")
#           reassign = dist_i_raw[,j] <  dist_0_raw.min
#           sum_j_obj = sum(c(dist_obj_j[reassign,j],dist_0_obj[!reassign]))  
#           save.df$objective[candidates[j]+1] = sum_j_obj
#           }
#         }
#       return(save.df)
#       }
#     
#     
# Identity_parks = distanc0r_for_heavy_compuation(LSOAcandidates = poly_df,
#                                                 population =  poly_df,
#                                                 objective = "x^2 * population$pop",
#                                                 events = parkrun_marker)
# 
# 
# 
# 
# top_first_identity_parks = poly_df
# top_first_identity_parks$rank = rank(Identity_parks$objective[-1])
# top_first_identity_parks$diff = (Identity_parks$objective[-1] - Identity_parks$objective[1])/1e+9 # in millions
# top_first_identity_parks$change = round(top_first_identity_parks$diff/Identity_parks$objective[1],5)*100
# top_n = 1:500
# top_first_identity_parks = subset(top_first_identity_parks,top_first_identity_parks$rank %in% top_n)
# 
# shapefile(top_first_identity_parks,filename = "./output/top_first_identity_parks", overwrite =T)
# write.csv(Identity_parks,file = "./output/Identity_parks.csv",row.names = F)
#  rm("Identity_parks","top_n") # free memory
#
# # Spoiler: Park candidate results are robust!

