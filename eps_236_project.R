
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
pan = read.csv("pandora_timed_2019_11_16.csv")

panlon = -71.1046
panlat = 42.3502

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

filter_distance_inside = function(dat1, latref, lonref){
  lats_1 = dat1$cornerlat1
  lats_2 = dat1$cornerlat2
  lats_3 = dat1$cornerlat3
  lats_4 = dat1$cornerlat4
  lons_1 = dat1$cornerlon1
  lons_2 = dat1$cornerlon2
  lons_3 = dat1$cornerlon3
  lons_4 = dat1$cornerlon4
  mask = matrix(0,1,length(lats_1))

  for (j in length(lats_1)){
    mask[j] = point.in.polygon(latref,lonref,c(lats_1[j],lats_2[j],lats_3[j],lats_4[j],lats_1[j]), c(lons_1[j],lons_2[j],lons_3[j],lons_4[j],lons_1[j]))

  }

  return(mask)
}



closest_pan_time = function(pan, trop_time) {
  pan_time = pan['timestamp']
  answer = which(abs(pan_time - trop_time) == min(abs(pan_time - trop_time)))
  return(answer)
}

time_avg_pandora = function(time, times, pan, window_mins) {
  window_s = window_mins * 60
  time_mask = (times > (time - window_s)) & (times < (time + window_s))
  no2 = pan[time_mask, 'VCD.DU.']
  out = mean(no2)
  return(out)
  
}
time_avg_pandora_errors = function(time, times, pan, window_mins) {
  window_s = window_mins * 60
  time_mask = (times > (time - window_s)) & (times < (time + window_s))
  un = pan[time_mask, 'UNCERTAINTY.DU.']
  out = mean(un)/(length(un))^(1/2)
  return(out)
  
}

#functiontesting
times = pan$timestamp
time = times[400]
win = 20
#for (j in 400:500){
print(time_avg_pandora(time, times, pan, win))
#}
#functiontesting




match_pan_trop = function(pan,
                          trop,
                          thresh,
                          panlat,
                          panlon,
                          avg_time = 10) {
  mask = filter_distance(trop, panlat, panlon, thresh)
#  mask = filter_distance_inside(trop,panlat,panlon)
 # print(mask&masktest)
  trop_time = trop[mask, 'timestamp']
  trop_no2 = trop[mask, 'no2']*(6.02214e19)
  trop_error = trop[mask,'no2_error']*(6.02214e19)
  pan_indices = 1:length(trop_time)
  pan_no2_avg = 1:length(trop_time)
  pan_no2_avg_errors = 1:length(trop_time)
  for (j in 1:length(trop_time)) {
    pan_indices[j] = closest_pan_time(pan, trop_time[j])
    pan_no2_avg[j] = time_avg_pandora(pan$timestamp[pan_indices[j]], pan$timestamp, pan, avg_time)
    pan_no2_avg_errors[j]= time_avg_pandora_errors(pan$timestamp[pan_indices[j]], pan$timestamp, pan, avg_time)
  }

  pan_no2 = pan$VCD.DU.[pan_indices] * 4.462e-8*6.02214e23
  pan_no2_error = pan$UNCERTAINTY.DU.[pan_indices] * 4.462e-8*6.02214e23
  pan_no2_avg = pan_no2_avg * 4.462e-8*6.02214e23
  pan_no2_avg_errors = pan_no2_avg_errors*4.462e-8*6.02214e23
  
  plot(
    pan_no2,
    trop_no2,
    xlab = "Pandora NO2 molecules/cm2",
    ylab = "TROPOMI NO2 molecules/cm2",
    xlim = c(0, 2e16),
    ylim = c(0, 2e16)
  )
  
  boot_slopes = 1:100
  #bootstrap
  for (j in 1:1000){
  dummy = 1:length(pan_indices)
  boot_i = sample(dummy,length(dummy),replace = TRUE)
  boot_pan = pan_indices[boot_i]
  
  boot_pan_no2 = pan$VCD.DU.[boot_pan] * 4.462e-8*6.02214e23
  boot_pan_no2_error = pan$UNCERTAINTY.DU.[boot_pan] * 4.462e-8*6.02214e23
  
  boot_trop_no2 = trop_no2[boot_i]
  boot_trop_no2_error = trop_error[boot_i]
  
  fit = lm(f = boot_trop_no2 ~ boot_pan_no2)
  
  yorkfit = york(cbind(boot_pan_no2, boot_pan_no2_error,boot_trop_no2,boot_trop_no2_error))
  
  slopey = yorkfit$b[1]
  
  intercepty = yorkfit$a[1]
  
  x = boot_pan_no2
  
  yy = slopey*x+intercepty
  lines(c(0,x),c(intercepty,yy))
  
  
#  print(fit)
  boot_slopes[j]=slopey
  }
#  hist(boot_slopes)
  print(max(boot_slopes))
  print(min(boot_slopes))
  
  
  
  
  #bootstrap
  
  
  
  
  
  
  
  
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

  
  yorkfit_avg = york(cbind(pan_no2_avg, pan_no2_avg_errors,trop_no2,trop_error))
  
  slopey_avg = yorkfit_avg$b[1]
  print(slopey_avg)
  
  intercepty_avg = yorkfit_avg$a[1]
  
  yy_avg = slopey_avg*x+intercepty_avg
  
  
  yorkfit = york(cbind(pan_no2, pan_no2_error,trop_no2,trop_error))
  
  slopey = yorkfit$b[1]
  print(slopey)

  intercepty = yorkfit$a[1]

  yy = slopey*x+intercepty
  
  
#  lines(c(0,x), c(intercept,y), col = 'green')
#  lines(c(0,x), c(intercepta,ya), col = 'blue')
  lines(c(0,x),c(intercepty,yy), col='green')
  

  lines(c(0,x),c(intercepty_avg,yy_avg), col='red')
  lines(c(0, 2e16), c(0, 2e16), col = 'blue')
  
#  lines(c(0, 5), c(intercepta, (intercepta + 5)))
  

  #  fake_trop_error = 0*1:length(pan_no2)
  
  #  york_vector = t(matrix(c(pan_no2,pan_no2_error,trop_no2,fake_trop_error)), nrow=1,ncol=4)
  #  print(york_vector)
  
  #  print(dim(york_vector))
  #  print(york(york_vector))
  
  
  return(cbind(yorkfit,yorkfit_avg))
}

fit2 = match_pan_trop(pan,trop,3000,panlat,panlon)

#fit3 = match_pan_trop(pan,trop,3000,panlat,panlon)

#fit4 = match_pan_trop(pan, trop, 4000, panlat, panlon)

#fit5 = match_pan_trop(pan,trop,5000,panlat,panlon)

#fitbig = match_pan_trop(pan,trop,15000,panlat,panlon)



