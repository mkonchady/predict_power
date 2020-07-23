#
# predict power generation for a given location and configuration
#
from sun_functions import *
from datetime import *
import calendar
import csv
from collections import defaultdict

# get monthly cloud percentages for the city
def get_cloud_data():
    INFILE = "cloud.csv"
    with open(INFILE) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        cloud_arr = []
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                cloud_arr.append(float(row[1]))
    return cloud_arr

# return the month and day for the date
def daynum_to_date(year : int, daynum : int) -> datetime.date:
    month = 1
    day = daynum
    while month < 13:
        month_days = calendar.monthrange(year, month)[1]
        if day <= month_days:
            return (day, month)
        day -= month_days
        month += 1
    return (day, month)

# calculate the power for the given day and location
def calc_power(LAT, LON, year, month, mday, peak, tilt, rotation, cloud):
    sund = get_sun_parms(LAT, LON, year, month, mday, 12.5)

    sunrise = int(sund["sunrise_time"])
    sunset = int(sund['sunset_time'])

    C = (5.0 / 6.0) * peak  # approximate
    M = -0.009 * peak
    tot_power = 0
    #for (minutes = sunrise; minutes < sunset; minutes = minutes + 1):
    PERIOD = 5
    for minutes in range (sunrise, sunset, PERIOD):
        theta = get_incident_angle(LAT, LON, year, month, mday, minutes / 60.0, rotation, tilt )
        if (theta > 0 and theta < 90):
            #x = (0.23 * cloud + 0.77 * theta)
            x = (0.25 * cloud + 0.75 * theta)
            current_power = M * x + C 
            current_power = 0.86 * current_power
            tot_power += current_power / (MINS_PER_HOUR / PERIOD)
    return tot_power            

# main section
# Constant parameters
MINS_PER_HOUR = 60.0
LAT = 33.17    
LON = -96.82 
year = 2020
peak = 10000
tilt = int(LAT)
rotation = 180

tot_power = 0
cloud_arr = get_cloud_data()
month_power = defaultdict(int)
for day in range (1, 366):
    (mday, month) = daynum_to_date(year, day)
    cloud = cloud_arr[month-1]
    pred_power = calc_power(LAT, LON, year, month, mday, peak, tilt, rotation, cloud)
    tot_power = tot_power + pred_power
    month_power[month] = month_power[month] + pred_power
    #print (mday, month, cloud, pred_power)

# dump the power by month
f = open ("city.dat", "w")
for key in range (1,13):
    outline = str(key) + "," + str(month_power[key])
    f.write(outline + "\n")
f.close()

print ("Annual power: " + str(tot_power))