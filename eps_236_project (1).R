
library(sp)
library(rgdal)
library(maps)
library(dplyr)
library(geosphere)
library(IsoplotR)
library(ncdf4)
path = "/Users/jackbruno/Documents/EPS_236/filtered_data/"
setwd(path)

trop = read.csv('total_output_2019_11_22_july_only.csv')
trop_weekend = read.csv("total_output_2019_11_29_weekend.csv")
trop_weekday = read.csv("total_output_2019_11_29_weekday.csv")

pan = read.csv("pandora_timed_2019_11_16.csv")

panlon = -71.1046
panlat = 42.3502


#This function filters for pixel distance based off of calculated geometric distance
filter_distance = function(dat1, latref, lonref, thresh) {
  lats = dat1['lats'][[1]]
  lons = dat1['lons'][[1]]
  dist = 1:length(lats)
  for (j in 1:length(lats)) {
    dist[j] = distm(c(latref, lonref), c(lats[[j]], lons[[j]]), fun = distHaversine)
  }
  dist_mask = dist < thresh
  return(dist_mask)
}

#This function filters points based on pixel corner values
filter_distance_inside = function(dat1, latref, lonref) {
  lats_1 = dat1$cornerlat1
  lats_2 = dat1$cornerlat2
  lats_3 = dat1$cornerlat3
  lats_4 = dat1$cornerlat4
  lons_1 = dat1$cornerlon1
  lons_2 = dat1$cornerlon2
  lons_3 = dat1$cornerlon3
  lons_4 = dat1$cornerlon4
  mask = matrix(0, 1, length(lats_1))
  
  for (j in 1:length(lats_1)) {
    mask[j] = point.in.polygon(
      latref,
      lonref,
      c(lats_1[j], lats_2[j], lats_3[j], lats_4[j], lats_1[j]),
      c(lons_1[j], lons_2[j], lons_3[j], lons_4[j], lons_1[j])
    )
    
  }
  
  return(mask)
}


#This function finds the closest pandora time to a given tropomi input timestamp
closest_pan_time = function(pan, trop_time) {
  pan_time = pan['timestamp']
  answer = which(abs(pan_time - trop_time) == min(abs(pan_time - trop_time)))
  return(answer)
}

#This function time averages the pandora data and returns a mean and error on the mean
time_avg_pandora = function(time, times, pan, window_mins) {
  window_s = window_mins * 60
  time_mask = (times > (time - window_s)) &
    (times < (time + window_s))
  no2 = pan[time_mask, 'VCD.DU.']
  er = sum(pan[time_mask, 'UNCERTAINTY.DU.'])**(1/2)/(length(no2))
  out = mean(no2)
  return(cbind(out,er))
  
}

#This function fits and bootstraps the data. This is the MAIN EVENT
match_pan_trop = function(pan,
                          trop,
                          thresh,
                          panlat,
                          panlon,
                          avg_time = 10) {
  #mask = filter_distance(trop, panlat, panlon, thresh)
  mask = as.logical(filter_distance_inside(trop,panlat,panlon))
  # print(mask&masktest)
  trop_time = trop[mask, 'timestamp']
  trop_no2 = trop[mask, 'no2'] * (6.02214e19)
  trop_error = trop[mask, 'no2_error'] * (6.02214e19)
  pan_indices = 1:length(trop_time)
  pan_no2_avg = 1:length(trop_time)
  pan_no2_avg_errors = 1:length(trop_time)
  for (j in 1:length(trop_time)) {
    pan_indices[j] = closest_pan_time(pan, trop_time[j])
    panavg = time_avg_pandora(pan$timestamp[pan_indices[j]], pan$timestamp, pan, avg_time)
    pan_no2_avg[j] = panavg[1]
    pan_no2_avg_errors[j] = panavg[2]
  }
  
  pan_no2 = pan$VCD.DU.[pan_indices] * 4.462e-8 * 6.02214e23
  pan_no2_error = pan$UNCERTAINTY.DU.[pan_indices] * 4.462e-8 * 6.02214e23
  pan_no2_avg = pan_no2_avg * 4.462e-8 * 6.02214e23
  pan_no2_avg_errors = pan_no2_avg_errors * 4.462e-8 * 6.02214e23
  
  plot(
    pan_no2,
    trop_no2,
    xlab = "Pandora NO2 molecules/cm2",
    ylab = "TROPOMI NO2 molecules/cm2",
    xlim = c(0, 2e16),
    ylim = c(0, 2e16)
  )
  
  points(
    pan_no2_avg,
    trop_no2,
    xlab = "Pandora NO2 molecules/cm2",
    ylab = "TROPOMI NO2 molecules/cm2",
    xlim = c(0, 2e16),
    ylim = c(0, 2e16),
    col = 'violet'
  )
  
  
  boot_slopes = 1:1000
  #bootstrap
  for (j in 1:1000) {
    dummy = 1:length(pan_no2_avg)
    boot_i = sample(dummy, length(dummy), replace = TRUE)
#    boot_pan = pan_indices[boot_i]
    
#    boot_pan_no2 = pan$VCD.DU.[boot_pan] * 4.462e-8 * 6.02214e23
#    boot_pan_no2_error = pan$UNCERTAINTY.DU.[boot_pan] * 4.462e-8 * 6.02214e23
    
    boot_pan_no2 = pan_no2_avg[boot_i]
    boot_pan_no2_error = pan_no2_avg_errors[boot_i]
    
    boot_trop_no2 = trop_no2[boot_i]
    boot_trop_no2_error = trop_error[boot_i]
    
    
    yorkfit = york(cbind(
      boot_pan_no2,
      boot_pan_no2_error,
      boot_trop_no2,
      boot_trop_no2_error
    ))
    
    slopey = yorkfit$b[1]
    
    intercepty = yorkfit$a[1]
    
    x = c(0,boot_pan_no2,2e16)
    
    yy = slopey * x + intercepty
    lines(x,  yy)
    
    
    #  print(fit)
    boot_slopes[j] = slopey
  }

  print(max(boot_slopes))
  print(min(boot_slopes))
  #bootstrap
  
  
  
  
  
  
  
  #BAD FITTING
  # fit = lm(f = trop_no2 ~ pan_no2)
  #
  # slope = fit$coefficients[2]
  # intercept = fit$coefficients[1]
  # x = pan_no2
  # y = slope * pan_no2 + intercept
  #
  #
  #
  # fita = lm(f = trop_no2 ~ pan_no2_avg)
  # slopea = fita$coefficients[2]
  # intercepta = fita$coefficients[1]
  # ya = slopea * x + intercepta
  #BAD FITTING
  
  
  yorkfit_avg = york(cbind(pan_no2_avg, pan_no2_avg_errors, trop_no2, trop_error))
  
  slopey_avg = yorkfit_avg$b[1]
  print(slopey_avg)
  
  intercepty_avg = yorkfit_avg$a[1]
  
  yy_avg = slopey_avg * x + intercepty_avg
  
  
  yorkfit = york(cbind(pan_no2, pan_no2_error, trop_no2, trop_error))
  
  slopey = yorkfit$b[1]
  print(slopey)
  
  intercepty = yorkfit$a[1]
  
  yy = slopey * x + intercepty
  

  lines(x, yy, col = 'green')
  lines(x,  yy_avg, col = 'red')
  lines(c(0, 2e16), c(0, 2e16), col = 'blue')
  
  
  return(cbind(yorkfit, yorkfit_avg))
}




fit_all = match_pan_trop(pan, trop, 3000, panlat, panlon,10)
fit_weekday = match_pan_trop(pan, trop_weekday, 3000, panlat, panlon,10)
fit_weekend = match_pan_trop(pan, trop_weekend, 3000, panlat, panlon,10)






