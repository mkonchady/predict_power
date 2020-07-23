#
# solar time functions 
# Translation of NOAA's code from
# https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
#  based on equations from Astronomical Algorithms, by Jean Meeus
#
import time
from datetime import datetime, date, timedelta
import math
from timezonefinder import TimezoneFinder
import tz
import numexpr
import pytz


monthName = [ "January", "February", "March", "April", "May", "June", "July", "August", \
    "September", "October", "November", "December"]
monthDays = [31,28,31,30,31,30,31,31,30,31,30,31]
monthAbbrev = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# is the string a floating point number
def isfloat(value):
  try:
    float(str(value))
    return True
  except ValueError:
    return False

# pad the number of lead zeroes
def zeroPad(n, leads):
  n = str(n)
  while (len(n) < leads):
    n = '0' + n
  return n

# take the mod of a floating point 360
def mod360(number):
    return number - math.floor(number / 360.0) * 360

# fix the longitude to show negative values
def fix_lon(lon):
    if (lon <= 180):
        return lon
    return -(360 - lon)

def fix_bearing(angle):
    if (angle <= 180):
        return angle
    return (360 - angle)

def sin(degrees):
    return math.sin(math.radians(degrees))

def cos(degrees):
    return math.cos(math.radians(degrees))

def tan(degrees):
    return math.tan(math.radians(degrees))

def acos(num):
    if (num > 1):
        num = 1
    if (num < -1):
        num = -1
    if (num == 1):
        return 0
    return math.degrees(1.570796327 - math.atan(num / math.sqrt(-1 * num * num + 1)))

def is_dst(zonename):
    tz = pytz.timezone(zonename)
    now = pytz.utc.localize(datetime.utcnow())
    return now.astimezone(tz).dst() != timedelta(0)

def get_timezone(LAT, LON):
    tf = TimezoneFinder()
    return tf.timezone_at(lng=LON, lat=LAT) 

# get the time difference between the local time zone and gmt in hours
# lon is expressed from -180 to +180
def get_gmt_time_diff(LAT, LON):
    timezone = get_timezone(LAT, LON)
    if (timezone is None):  # could not get the timezone and so use the longitude
        diff = LON * 12 / 180.0
        frac, numb = math.modf(diff)
        if (frac < 0.25):
            return math.floor(diff)
        if (frac < 0.75):
            return math.floor(numb + 0.5)
        return math.ceil(numb)
    
    diff = tz.get_tz_hours(timezone)
    if (is_dst(timezone)):
        diff = diff + 1
    return diff

# return the current time as current year, month, day, and hour as a fraction of hours
def get_current_time(LAT, LON):
    unix_epoch = time.time()
    tz_diff = math.ceil(get_gmt_time_diff(LAT, LON) * 3600)
    ltime = time.gmtime(unix_epoch + tz_diff )
    return ltime.tm_year, ltime.tm_mon, ltime.tm_mday, \
        ltime.tm_hour + (ltime.tm_min / 60.0) + (ltime.tm_sec / 3600.0)

# return the julian day for the date
# check https://idlastro.gsfc.nasa.gov/ftp/pro/astro/jdcnv.pro
def get_julian_day(year, month, day):
    if (month <= 2):
        year -= 1
        month += 12
    century = math.floor( year / 100)
    leap_day = 2 - century + math.floor( century / 4)
    p1 = math.floor(365.25 * (year + 4716))
    p2 = math.floor(30.6001 * (month + 1))
    p3 = day + leap_day - 1524.5
    return p1 + p2 + p3

# return the julian day for the date and time in minutes 
def get_julian_day1(year, month, day, ctime):
    jd = get_julian_day(year, month, day)
    return jd + (ctime / 60.0) / 24.0
    #return jd + (hour + (min / 60.0) + (sec / 3600.0) ) / 24.0

# return the julian century fraction
def get_julian_century(julian_day):
    return (julian_day - 2451545.0) / 36525.0

# return the julian century fraction
def get_julian_century1(year, month, mday, ctime):
    julian_day = get_julian_day1(year, month, mday, ctime)
    return (julian_day - 2451545.0) / 36525.0

def isLeapYear(yr): 
    return ((yr % 4 == 0 and yr % 100 != 0) or yr % 400 == 0)
  
# return a string from minutes
def timeString(minutes, flag):
    if (flag == 3):
        return time.strftime("%H:%M:%S", time.gmtime(60 * minutes))
    return time.strftime("%H:%M", time.gmtime(60 * minutes))

#  calculate the Geometric Mean Longitude of the Sun
def get_geom_mean_long_sun(julian_century):
    longitude = 280.46646 + julian_century * (36000.76983 + julian_century * (0.0003032))
    while(longitude > 360.0):
        longitude -= 360.0
    while(longitude < 0.0):
        longitude += 360.0
    return longitude		# in degrees

# calculate the Geometric Mean Anomaly of the Sun
def get_geom_mean_anom_sun(julian_century):
  return 357.52911 + julian_century * (35999.05029 - 0.0001537 * julian_century)

# calculate the eccentricity of earth's orbit
def get_eccent_earth_orbit(julian_century):
    return 0.016708634 - julian_century * (0.000042037 + 0.0000001267 * julian_century)

# calculate the equation of center for the sun
def get_sun_eq_of_center(julian_century):
  m = get_geom_mean_anom_sun(julian_century)
  sinm = sin(m)
  sin2m = sin(m+m)
  sin3m = sin(m+m+m)
  return sinm * (1.914602 - julian_century * (0.004817 + 0.000014 * julian_century)) + \
    sin2m * (0.019993 - 0.000101 * julian_century) + sin3m * 0.000289  # in degrees

# calculate the true longitude of the sun	
def get_sun_true_long(julian_century):
    return get_geom_mean_long_sun(julian_century) + get_sun_eq_of_center(julian_century)

# calculate the true anomaly of the sun
def get_sun_true_anom(julian_century):
    return get_geom_mean_anom_sun(julian_century) + get_sun_eq_of_center(julian_century)

# calculate the distance to the sun in AU
def get_sun_rad_vector(julian_century):
    anom = get_sun_true_anom(julian_century)
    eccent = get_eccent_earth_orbit(julian_century)
    return (1.000001018 * (1 - eccent * eccent)) / (1 + eccent * cos(anom)) # in aus

# calculate the apparent longitude of the sun
def get_sun_app_long(julian_century):
    true_long = get_sun_true_long(julian_century)
    omega = 125.04 - 1934.136 * julian_century
    return true_long - 0.00569 - 0.00478 * sin(omega)

# calculate the mean obliquity of the ecliptic
def get_mean_obliq_ecliptic(julian_century):
    seconds = 21.448 - julian_century * (46.8150 + julian_century*(0.00059 - julian_century*(0.001813)))
    return 23.0 + (26.0 + (seconds/60.0))/60.0 # in degrees

# calculate the corrected obliquity of the ecliptic
def get_obliq_corr(julian_century):
  e0 = get_mean_obliq_ecliptic(julian_century)
  omega = 125.04 - 1934.136 * julian_century
  return e0 + 0.00256 * cos(omega) # in degrees

# calculate the right ascension of the sun	in degrees
def get_sun_rt_ascen(julian_century):
    epsilon_d = get_obliq_corr(julian_century)
    lamb_d = get_sun_app_long(julian_century)
    num = (cos(epsilon_d) * sin(lamb_d))
    denom = cos(lamb_d)
    return math.degrees(math.atan2(num, denom))

# calculate the declination of the sun
def get_sun_decln(julian_century):
    epsilon_d = get_obliq_corr(julian_century)
    lambda_d = get_sun_app_long(julian_century)
    sint = sin(epsilon_d) * sin(lambda_d)
    return math.degrees(math.asin(sint))

# calculate the difference between true solar time and mean	solar time
def get_equation_of_time(julian_century):
    obliq_corr = get_obliq_corr(julian_century)
    geom_mean_anom_sun = get_geom_mean_anom_sun(julian_century)
    geom_mean_long_sun = get_geom_mean_long_sun(julian_century)
    eccent_earth_orbit = get_eccent_earth_orbit(julian_century)
    y = tan(obliq_corr / 2.0) * tan(obliq_corr / 2.0)
    sin2 = sin(2.0 * geom_mean_long_sun)
    sinm = sin(geom_mean_anom_sun)
    cos2 = cos(2.0 * geom_mean_long_sun)
    sin4 = sin(4.0 * geom_mean_long_sun)
    sin2m = sin(2.0 * geom_mean_anom_sun)
    etime = y * sin2 - (2.0 * eccent_earth_orbit * sinm) + \
        (4.0 * eccent_earth_orbit * y * sinm * cos2) - \
        (0.5 * y * y * sin4) - \
        (1.25 * eccent_earth_orbit * eccent_earth_orbit * sin2m)
    return 4.0 * math.degrees(etime)	# in minutes of time

# get the sunrise hour angle for the LAT and declination
def get_ha_sunrise(LAT, sun_decln):
    if (LAT > 66 or LAT < -68):
        return 0
    cos_arg = (cos(90.833) / (cos(LAT)*cos(sun_decln))-tan(LAT) * tan(sun_decln))
    return math.degrees(math.acos(cos_arg))

# return the solar noon when sun crosses the given longitude
def get_solar_noon(julian_day, LON, tz_diff, formatted, dst = False):
    tnoon = get_julian_century(julian_day - (LON / 360.0))
    eqTime = get_equation_of_time(tnoon)
    solNoonOffset = 720.0 - (LON * 4) - eqTime # in minutes
    newt = get_julian_century(julian_day + (solNoonOffset / 1440.0))
    eqTime = get_equation_of_time(newt)
    solNoonLocal = 720 - (LON * 4) - eqTime + (tz_diff * 60.0) # in minutes
    if(dst):
        solNoonLocal += 60.0
    while (solNoonLocal < 0.0):
        solNoonLocal += 1440.0
    while (solNoonLocal >= 1440.0):
        solNoonLocal -= 1440.0
    if (formatted):
        return timeString(solNoonLocal, 3)
    return solNoonLocal

# get sunrise / set in UTC
def get_sunrise_set_UTC(rise, julian_day, LAT, LON):
    julian_century = get_julian_century(julian_day)
    eqTime = get_equation_of_time(julian_century)
    sun_decln = get_sun_decln(julian_century)
    hourAngle = get_ha_sunrise(LAT, sun_decln)
    if (not rise):
        hourAngle = -hourAngle
    delta = LON + hourAngle
    return 720 - (4.0 * delta) - eqTime	# in minutes

#rise = 1 for sunrise, 0 for sunset
def get_sunriseset(rise, julian_day, LAT, LON, timezone, dst):
    timeUTC = get_sunrise_set_UTC(rise, julian_day, LAT, LON)
    newTimeUTC = get_sunrise_set_UTC(rise, julian_day + timeUTC/1440.0, LAT, LON)
    if (isfloat(newTimeUTC)):
        timeLocal = newTimeUTC + (timezone * 60.0)
        if (dst):
            timeLocal += 60.0
        if ( (timeLocal < 0.0) or (timeLocal > 1440.0) ):
            if (timeLocal < 0):
                increment = 1
            else:
                increment = -1
            while ((timeLocal < 0.0) or (timeLocal >= 1440.0)):
                timeLocal += increment * 1440.0
    return timeLocal

# get the azimuth, elevation, and zenith angle for the time, local time in minutes
def get_az_el_ze(julian_century, LAT, LON, tz_diff, localtime):
    eqTime = get_equation_of_time(julian_century)
    theta = get_sun_decln(julian_century)
    solarTimeFix = eqTime + 4.0 * LON - 60.0 * tz_diff
    trueSolarTime = localtime + solarTimeFix
    while (trueSolarTime > 1440):
        trueSolarTime -= 1440
    hourAngle = trueSolarTime / 4.0 - 180.0
   
    if (hourAngle < -180): 
        hourAngle += 360.0

    csz = sin(LAT) * sin(theta) + cos(LAT) * cos(theta) * cos(hourAngle)
    if (csz > 1.0):
        csz = 1.0
    elif (csz < -1.0): 
        csz = -1.0
  
    zenith = math.degrees(math.acos(csz))
    azDenom = cos(LAT) * sin(zenith)
    if (abs(azDenom) > 0.001):
        azRad = ( (sin(LAT) * cos(zenith)) - sin(theta) ) / azDenom
        if (abs(azRad) > 1.0):
            if (azRad < 0):
	            azRad = -1.0
            else:
	            azRad = 1.0
    
        azimuth = 180.0 - math.degrees(math.acos(azRad))
        if (hourAngle > 0.0):
            azimuth = -azimuth
    else:
        if (LAT > 0.0):
            azimuth = 180.0
        else: 
            azimuth = 0.0
    
    if (azimuth < 0.0):
        azimuth += 360.0
  
    exoatmElevation = 90.0 - zenith

    #Atmospheric Refraction correction
    if (exoatmElevation > 85.0):
        refractionCorrection = 0.0
    else: 
        te = tan (exoatmElevation)
        if (exoatmElevation > 5.0):
            refractionCorrection = 58.1 / te - 0.07 / (te*te*te) + 0.000086 / (te*te*te*te*te)
        elif (exoatmElevation > -0.575):
            refractionCorrection = 1735.0 + exoatmElevation * (-518.2 + exoatmElevation * (103.4 + exoatmElevation * (-12.79 + exoatmElevation * 0.711) ) )
        else:
            refractionCorrection = -20.774 / te
    
        refractionCorrection = refractionCorrection / 3600.0

    #refractionCorrection = 0
    solarZen = zenith - refractionCorrection
    elevation = math.floor((90.0 - solarZen) * 100 + 0.5) / 100.0
    azimuth = math.floor(azimuth * 100 + 0.5) / 100.0
    return (azimuth, elevation, solarZen)

#
# get the earth longitude of the sun
# https://bit.ly/3a52bdI  from Mathworks site
# Keplerian Elements for the Sun (geocentric)
def get_sun_lon(julian_day):
    d = julian_day - 2451543.5                  # get the julian date from Jan 1, 2000 
    lst = get_lst(julian_day)                   # get local sidereal time
    frac_date = d - math.floor(d)               # fraction of julian date
    ra = get_sun_rt_ascen(get_julian_century(julian_day))
    return fix_lon(mod360(ra - lst - frac_date * 360))

# get local sidereal time
def get_lst(julian_day):
    d = julian_day - 2451543.5                  # get the julian date from Jan 1, 2000 
    w = 282.9404 + 4.70935e-5 * d               #    (longitude of perihelion degrees)
    m = mod360(356.0470 + 0.9856002585 * d)     # (mean anomaly degrees)
    mean_long = mod360(w + m)                   # (Sun's mean longitude degrees)
    lst = (mean_long + 180)                     # local sidereal time
    return lst

# return a dictinary of solar parameters for a given lat/lon in decimal 
# and time in fraction of hours (1:30 pm = 13.5)
# if year is provided, then ALL other optional parms must be provided
def get_sun_parms(LAT, LON, year = -1, month = -1, mday = -1, time = -1):
    if (year == -1):
        year, month, mday, time = get_current_time(LAT, LON)

    #tz_diff = get_gmt_time_diff(LON) # in hours
    tz_diff = get_gmt_time_diff(LAT, LON)
    julian_day = get_julian_day(year, month, mday) + (time - tz_diff) / 24.0
    julian_century = get_julian_century(julian_day)
    geom_mean_long_sun = get_geom_mean_long_sun(julian_century)
    geom_mean_anom_sun = get_geom_mean_anom_sun(julian_century)
    eccent_earth_orbit = get_eccent_earth_orbit(julian_century)

    sun_eq_of_ctr = get_sun_eq_of_center(julian_century)
    sun_true_long = get_sun_true_long(julian_century)
    sun_true_anom = get_sun_true_anom(julian_century)
    sun_rad_vector = get_sun_rad_vector(julian_century)
    sun_app_long = get_sun_app_long(julian_century)

    mean_obliq_ecliptic = get_mean_obliq_ecliptic(julian_century)
    obliq_corr = get_obliq_corr(julian_century)
    sun_rt_ascen = get_sun_rt_ascen(julian_century)
    sun_decln = get_sun_decln(julian_century)
    eq_of_time = get_equation_of_time(julian_century)
    ha_sunrise = get_ha_sunrise(LAT, sun_decln)
    solar_noon = get_solar_noon(get_julian_day(year, month, mday), LON, tz_diff, True)

    #sunrise_time = timeString(get_sunriseset(True, julian_day, LAT, LON, tz_diff, False), 2) 
    #sunset_time = timeString(get_sunriseset(False, julian_day, LAT, LON, tz_diff, False), 2)
    sunrise_time = get_sunriseset(True, julian_day, LAT, LON, tz_diff, False)
    sunset_time = get_sunriseset(False, julian_day, LAT, LON, tz_diff, False)
    sunlight_duration = timeString(8.0 * ha_sunrise, 3)    # minutes
    sd = sunset_time - sunrise_time
    sunlight_duration = 8.0 * ha_sunrise;   
    (solar_azimuth_angle, solar_elevation_angle, solar_zenith_angle) = \
        get_az_el_ze(julian_century, LAT, LON, tz_diff, time * 60)
    solar_lon = fix_lon(get_sun_lon(julian_day)) # get the location of the sun on earth
    solar_lat = sun_decln

    sund = {}
    sund["tz_diff"] = tz_diff
    sund["julian_day"] = julian_day
    sund["julian_century"] = julian_century
    sund['geom_mean_long_sun'] = geom_mean_long_sun
    sund["geom_mean_anom_sun"] = geom_mean_anom_sun
    sund["eccent_earth_orbit"] = eccent_earth_orbit
    sund["sun_eq_of_ctr"] = sun_eq_of_ctr 
    sund["sun_true_long"] = sun_true_long
    sund["sun_true_anom"] = sun_true_anom
    sund["sun_rad_vector"] =  sun_rad_vector
    sund["sun_app_long"] = sun_app_long
    sund["mean_obliq_ecliptic"] = mean_obliq_ecliptic
    sund["obliq_corr"] = obliq_corr
    sund["sun_rt_ascen"] = sun_rt_ascen
    sund["sun_decln"] = sun_decln
    sund["eq_of_time"] = eq_of_time
    sund["ha_sunrise"] = ha_sunrise
    sund["solar_noon"] = solar_noon
    sund["sunrise_time"] = sunrise_time
    sund["sunset_time"] = sunset_time
    sund["sunlight_duration"] = sunlight_duration
    sund["solar_zenith_angle"] = solar_zenith_angle
    sund["solar_elevation_angle"] = solar_elevation_angle
    sund["solar_azimuth_angle"] = solar_azimuth_angle
    sund['solar_lat'] = solar_lat
    sund['solar_lon'] = solar_lon

    return sund

# get the incident angle  and panel rotation/tilt angle
def get_incident_angle(LAT, LON, year, month, mday, time, panel_azimuth_angle, panel_tilt_angle):
    tz_diff = get_gmt_time_diff(LAT, LON)
    julian_day = get_julian_day(year, month, mday) + (time - tz_diff) / 24.0
    julian_century = get_julian_century(julian_day)
    (solar_azimuth_angle, solar_elevation_angle, solar_zenith_angle) = \
        get_az_el_ze(julian_century, LAT, LON, tz_diff, time * 60)
    if (solar_elevation_angle < 0):
        return -1
    solar_noon = get_solar_noon(julian_day, LON, tz_diff, False) / 60.0 # in hour fractions

    if (LAT > 0): # northern hemisphere -- panel facing south
        psi = abs(180 - panel_azimuth_angle)
        phi = abs(180 - solar_azimuth_angle)
        # gamma is angular difference between solar and panel azimuth angles
        if (panel_azimuth_angle < 180): # East of South  
            if time < solar_noon:       # morning
                gamma = phi - psi
            else:
                gamma = phi + psi
        else:                           # West of South
            if time < solar_noon:       # morning
                gamma = phi + psi
            else:
                gamma = phi - psi
    else:   # southern hemisphere -- panel facing north
        psi = fix_bearing(panel_azimuth_angle)
        phi = fix_bearing(solar_azimuth_angle)
        if (panel_azimuth_angle < 180): # east of north
            if time < solar_noon:       # morning
                gamma = phi - psi
            else:
                gamma = phi + psi
        else:                           # west of north
            if time < solar_noon:       # morning
                gamma = phi + psi
            else:
                gamma = phi - psi

    p1 = cos(solar_elevation_angle) * cos(gamma) * sin(panel_tilt_angle)
    p2 = sin(solar_elevation_angle) * cos(panel_tilt_angle)
    return math.degrees(math.acos(p1 + p2))

# Convert julian day to Local Sidereal Time
# check    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ct2lst.pro
# tz_diff takes into account DST
# julian_date calculated using local time
def ct2lst(julian_date, LAT, LON):
    
    c1 = 280.46061837
    c2 = 360.98564736629
    c3 = 0.000387933
    c4 = 38710000.0
    tz_diff = get_gmt_time_diff(LAT, LON)
    corrected_julian_date = julian_date - (tz_diff / 24.0)
    t0 = corrected_julian_date - 2451545
    t = t0 / 36525
    
    # Compute GST in seconds
    theta = c1 + (c2 * t0) + (t**2) * (c3 - (t / c4))
    
    # Compute LST in hours
    lst = (theta + LON) / 15.0
    
    # Deal with LST out of bounds
    if lst < 0.0:
        lst = 24.0 + lst % 24
    lst = lst % 24.0
    return lst

def mod(x, y):
    return x%y

# taken from http://www.stargazing.net/mas/al_az.htm
def azel2radec(julian_date, azimuth, elevation, LAT, LON):
    # Get the local Sidereal Time of the date
    LON = fix_lon(LON)
    lst = ct2lst(julian_date, LAT, LON)

    # calculate declination
    dec = math.degrees(math.asin(sin(elevation) * sin(LAT) + cos(elevation) * cos(azimuth) * cos(LAT)))

    # calculate right ascension
    numer = (sin(elevation) - sin(LAT) * sin(dec)) / (cos(LAT) * cos(dec))
    s = acos(numer)
    j = sin(azimuth)
    if (j > 0):
        s = 360 - s
    ra = lst - (s / 15)
    if (ra < 0):
        ra = ra + 24 
    return ra, dec