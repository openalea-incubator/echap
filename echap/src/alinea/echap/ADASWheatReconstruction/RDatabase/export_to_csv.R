# Permits to export the data from '.RData' to '../data/' directory.
# 
# Author: cchambon
###############################################################################


# import R data base
rm(meteodb)
attach('.RData')
meteodb <- get('meteodb')
detach()

# export to csv files
meteodb$HMj <- data.frame(date=meteodb$HMj$Date, Tmin=meteodb$HMj$Minimum, Tmax=meteodb$HMj$Maximum)

completeDates <- as.POSIXct(strptime(paste(meteodb$RMj$year, '-', meteodb$RMj$month, '-', meteodb$RMj$day, sep=''), '%Y-%m-%d'))
meteodb$RMj$date <- as.character(format(completeDates,'%Y-%m-%d'))
meteodb$RMj <- data.frame(date=meteodb$RMj$date, Tmin=meteodb$RMj$min, Tmax=meteodb$RMj$max)

tableList = c('HMj', 'RMh', 'RMj')
sapply(tableList, function(x) write.csv(meteodb[[x]], 
                                        file = paste('../data/', x, '.csv', sep=''), 
                                        row.names = FALSE))
