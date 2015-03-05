#!/usr/bin/env python

'''
timeutils.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Sept 2013

Contains various useful tools for dealing with time in astronomical contexts.

'''

import logging
import ConfigParser
import time
import os.path
import os

import numpy as np

# import the SPICE library
# this uses PySPICE: https://github.com/rca/PySPICE
# get CSPICE from here:
# http://naif.jpl.nasa.gov/naif/toolkit_C_PC_Linux_GCC_64bit.html
import spice

# import the SPICE kernels we need for calculating our values
CONF = ConfigParser.ConfigParser()
CONF.read('lcserver.conf')

spice_ephem_dir = os.path.abspath(CONF.get('data','spice_ephemerides'))
spice_planetdata_dir = os.path.abspath(CONF.get('data','spice_planetdata'))
spice_leapsecond_dir = os.path.abspath(CONF.get('data','spice_leapsecond'))

spice.furnsh(CONF.get('data','spice_ephemerides'))
spice.furnsh(CONF.get('data','spice_planetdata'))
spice.furnsh(CONF.get('data','spice_leapsecond'))


####################
## GENERAL CONFIG ##
####################

# setup a logger
LOGGER = logging.getLogger('timeutils')
LOGGER.addHandler(logging.NullHandler())


######################
## USEFUL CONSTANTS ##
######################

# physical constants
CLIGHT_KPS = 299792.458

# various JDs
JD1800 = 2378495.0
JD2000 = 2451545.0
JD2000INT = 2451545
JD2050 = 2469807.5

# conversion factors
MAS_P_YR_TO_RAD_P_DAY = 1.3273475e-11
ARCSEC_TO_RADIANS = 4.84813681109536e-6
KM_P_AU = 1.49597870691e8
SEC_P_DAY = 86400


#######################
## UTILITY FUNCTIONS ##
#######################

def precess_coordinates(ra, dec,
                        epoch_one, epoch_two,
                        jd=None,
                        mu_ra=0.0,
                        mu_dec=0.0,
                        outscalar=False):
    '''
    Precesses target coordinates ra, dec from epoch_one to epoch_two, given the
    jd of the observations, as well as the proper motion of the target mu_ra,
    mu_dec. Adapted from hatpipe/source/vartools/converttime.c [coordprecess].

    epoch_one, epoch_two = epochs (e.g. 1985.0, 2013.0, etc.)

    jd = Julian date (full JD, not reduced JD)

    ra = right ascension in decimal degrees
    dec = declination in decimal degrees

    mu_ra = proper motion in RA (mas/yr)
    mu_dec = proper motion in Dec (mas/yr)

    '''

    raproc, decproc = np.radians(ra), np.radians(dec)

    if ((mu_ra != 0.0) and (mu_dec != 0.0) and jd):

        jd_epoch_one = JD2000 + (epoch_one - epoch_two)*365.25
        raproc = (
            raproc +
            (jd - jd_epoch_one)*mu_ra*MAS_P_YR_TO_RAD_P_DAY/np.cos(decproc)
            )
        decproc = decproc + (jd - jd_epoch_one)*mu_dec*MAS_P_YR_TO_RAD_P_DAY

    ca = np.cos(raproc)
    cd = np.cos(decproc)
    sa = np.sin(raproc)
    sd = np.sin(decproc)

    if epoch_one != epoch_two:

        t1 = 1.0e-3 * (epoch_two - epoch_one)
        t2 = 1.0e-3 * (epoch_one - 2000.0)

        a = ( t1*ARCSEC_TO_RADIANS * (23062.181 + t2*(139.656 + 0.0139*t2) +
                                      t1*(30.188 - 0.344*t2+17.998*t1)) )
        b = t1*t1*ARCSEC_TO_RADIANS*(79.280 + 0.410*t2 + 0.205*t1) + a
        c = (
            ARCSEC_TO_RADIANS*t1*(20043.109 - t2*(85.33 + 0.217*t2) +
                                  t1*(-42.665 - 0.217*t2 - 41.833*t2))
            )
        sina, sinb, sinc = np.sin(a), np.sin(b), np.sin(c)
        cosa, cosb, cosc = np.cos(a), np.cos(b), np.cos(c)

        precmatrix = np.matrix([[cosa*cosb*cosc - sina*sinb,
                                 sina*cosb + cosa*sinb*cosc,
                                 cosa*sinc],
                                [-cosa*sinb - sina*cosb*cosc,
                                  cosa*cosb - sina*sinb*cosc,
                                  -sina*sinc],
                                [-cosb*sinc,
                                  -sinb*sinc,
                                  cosc]])

        precmatrix = precmatrix.transpose()

        x = (np.matrix([cd*ca, cd*sa, sd])).transpose()

        x2 = precmatrix * x

        outra = np.arctan2(x2[1],x2[0])
        outdec = np.arcsin(x2[2])


        outradeg = np.rad2deg(outra)
        outdecdeg = np.rad2deg(outdec)

        if outradeg < 0.0:
            outradeg = outradeg + 360.0

        if outscalar:
            return float(outradeg), float(outdecdeg)
        else:
            return outradeg, outdecdeg

    else:

        # if the epochs are the same and no proper motion, this will be the same
        # as the input values. if the epochs are the same, but there IS proper
        # motion (and a given JD), then these will be perturbed from the input
        # values of ra, dec by the appropriate amount of motion
        return np.degrees(raproc), np.degrees(decproc)


def earthmoon_orbital_elements(jd):
    '''
    Returns the approximate Kepler elements and rates for the Earth-Moon
    Barycenter at jd. Adapted from hatpipe/source/vartools/converttime.c
    [get_EM_orbital_elements].

    '''

    if (JD1800 < jd < JD2050):

        return (1.00000261, 0.00000562,       # a, adot
                0.01671123, -0.00004392,      # e, edot
                -0.00001531, -0.01294668,     # i, idot
                100.46457166, 35999.37244981, # L, Ldot
                102.93768193, 0.32327364,     # om, omdot
                0.0, 0.0)                     # Om, Omdot

    else:

        return (1.00000018, -0.00000003,      # a, adot
                0.01673163, -0.00003661,      # e, edot
                -0.00054346, -0.01337178,     # i, idot
                100.46691572, 35999.37306329, # L, Ldot
                102.93005885, 0.31795260,     # om, omdot
                -5.11260389, -0.24123856)     # Om, Omdot


def eccentric_anomaly(M, e):
    '''
    Returns the eccentric anomaly for an object of mean anomaly M, and
    eccentricity e. Adapted from hatpipe/source/vartools/transit.c
    [eccentricAnomaly]

    '''

    BISECTION_MAX_ITER = 1000
    ECC_ANOMALY_MAX_ERR = 1.0e-4

    Emin = M - 1.0
    Emax = M + 1.0
    E = M

    f = M + (e * np.sin(E)) - E
    dfm = 1.0 - (e  * np.cos(E))
    dE = dfm * ECC_ANOMALY_MAX_ERR
    counter = 0

    while ((dE > ECC_ANOMALY_MAX_ERR) and (counter < BISECTION_MAX_ITER)):

        f = M + (e * np.sin(E)) - E
        dfm = 1.0 - (e  * np.cos(E))
        dE = dfm * ECC_ANOMALY_MAX_ERR

        if f > dE:
            Emin = E
            dE = Emax - Emin

            if ((dfm * dE) > f):
                E = E + f/dfm
            else:
                E = 0.5 * (Emax + Emin)

        elif (f < (-dE)):

            Emax = E
            dE = Emax - Emin

            if ((dfm * dE) > -f):
                E = E + f/dfm
            else:
                E = 0.5 * (Emax + Emin)

        else:
            break

        counter = counter + 1

    return E


###########################
## JULIAN DATE FUNCTIONS ##
###########################

def unixtime_to_jd(unix_time):
    '''
    This converts UNIX time in seconds to a julian date.

    '''
    timestr = time.gmtime(unix_time)
    year, month, day, hour, minute, second = (timestr.tm_year,
                                              timestr.tm_mon,
                                              timestr.tm_mday,
                                              timestr.tm_hour,
                                              timestr.tm_min,
                                              timestr.tm_sec)

    a = (14 - month)/12
    y = year + 4800 - a
    m = month + 12*a - 3

    jdn = day + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400 - 32045

    jd = jdn + (hour-12.)/24. + minute/1440. + second/86400.

    return jd


def datetime_to_jd(dt):
    '''
    This converts a Python datetime object (naive, time in UT) to JD.

    '''

    year, month, day, hour, minute, second = (dt.year,
                                              dt.month,
                                              dt.day,
                                              dt.hour,
                                              dt.minute,
                                              (dt.second +
                                               dt.microsecond*(10.0**(-6))))

    a = (14 - month)/12
    y = year + 4800 - a
    m = month + 12*a - 3

    jdn = day + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400 - 32045

    jd = jdn + (hour-12.)/24. + minute/1440. + second/86400.

    return jd


def jd_now():
    '''
    Returns the JD at the current time.

    '''
    return unixtime_to_jd(time.time())


def jd_to_lmst(jd, longitude):
    '''
    Returns LMST (local mean sidereal time) in decimal hours for a given jd and
    longitude.

    Mostly stolen from:
    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ct2lst.pro

    '''

    # constants also stolen from:
    # http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ct2lst.pro
    c = [280.46061837, 360.98564736629, 0.000387933, 38710000.0]
    t0 = jd - JD2000
    t = t0/36525.0

    theta = c[0] + (c[1] * t0) + (t**2.0)*(c[2] - t/c[3])
    lst = (theta + longitude)/15.0

    if lst < 0:
        lst = 24.0 + (lst % 24.0)

    return lst % 24.0


def jd_to_mjd(jd):
    '''
    Converts Julian Date to Modified Julian Date.

    MJD = JD - 2400000.5

    '''

    return jd - 2400000.5


def mjd_to_jd(mjd):
    '''
    Converts Julian Date to Modified Julian Date.

    JD = MJD + 2400000.5

    '''

    return mjd + 2400000.5


def jd_to_hjd(jd, ra, dec,
              epoch=2000.0,
              mu_ra=0.0,
              mu_dec=0.0,
              timediff=False):
    '''
    Converts Julian Date to Heliocentric Julian Date for object at ra, dec
    observed at jd.

    epoch = epoch of coordinates, default is J2000.0

    mu_ra = proper motion in RA (mas/yr)

    mu_dec = proper_motion in Dec (mas/yr)

    timediff = True -> return difference in seconds between JD and HJD

    '''

    # if the epoch isn't J2000.0, or if there's any proper motion, get the
    # precessed coordinates for the object
    if (epoch != 2000.0) or (mu_ra != 0.0) or (mu_dec != 0.0):

        raproc, decproc = precess_coordinates(ra, dec, epoch, 2000.0,
                                              mu_ra=mu_ra, mu_dec=mu_dec)

        # FIXME: this makes this function work on scalars only, add vector
        # support later
        raproc, decproc = (float(np.radians(raproc[0])),
                           float(np.radians(decproc[0])))


    else:

        raproc, decproc = np.radians(ra), np.radians(dec)


    # get the Earth-Moon barycenter coordinates for jd
    (aEM, adotEM, eEM, edotEM, iEM, idotEM, LEM, LdotEM,
     omEM, omdotEM, OmEM, OmdotEM) = earthmoon_orbital_elements(jd)

    cd, sd, ca, sa  = (np.cos(decproc), np.sin(decproc),
                       np.cos(raproc), np.sin(raproc))

    # get the orbital elements for the EM system barycenter at time of
    # observation
    Tcent = (jd - JD2000)/36525.0
    aEM, eEM, iEM = (aEM + Tcent * adotEM,
                     eEM + Tcent*edotEM,
                     iEM + Tcent*idotEM)
    LEM, omEM, OmEM = (LEM + Tcent*LdotEM,
                       omEM + Tcent*omdotEM,
                       OmEM + Tcent*OmdotEM)

    # get the eccentric anomaly
    omperi = omEM - OmEM
    M = LEM - omEM
    M = M - 360.0*np.floor(M/360.0)

    if M > 180.0:
        M = M - 360.0

    M = np.radians(M)
    omperi = np.radians(omperi)
    OmEM = np.radians(OmEM)
    iEM = np.radians(iEM)
    E = eccentric_anomaly(M, eEM)

    # get the heliocentric coordinates of the EM barycenter in the EM orbital
    # plane
    xprime = np.array([aEM*(np.cos(E) - eEM),
                       aEM*np.sqrt(1.0 - eEM*eEM)*np.sin(E),
                       0.0])

    # get the coordinates in the J2000 ecliptic plane
    co, so = np.cos(omperi), np.sin(omperi)
    cO, sO = np.cos(OmEM), np.sin(OmEM)
    ci, si = np.cos(iEM), np.sin(iEM)

    xecl = np.array([(co*cO - so*sO*ci)*xprime[0]
                     + (-so*cO-co*sO*ci)*xprime[1],
                     (co*sO+so*cO*ci)*xprime[0]
                     + (-so*sO+co*cO*ci)*xprime[1],
                     so*si*xprime[0] + co*si*xprime[1]])

    # coordinates in the J2000 frame
    xeq = np.array([xecl[0],
                    0.917482139208287*xecl[1] - 0.397776978008764*xecl[2],
                    0.397776978008764*xecl[1] + 0.917482139208287*xecl[2]])

    # the Romer delay calculation
    delta_romer = ( (xeq[0]*cd*ca + xeq[1]*cd*sa + xeq[2]*sd) *
                    KM_P_AU/CLIGHT_KPS/SEC_P_DAY )

    if timediff:
        return delta_romer*SEC_P_DAY
    else:
        return jd + delta_romer


def jd_to_bjd(jd, ra, dec,
              obslat, obslon, obsalt,
              epoch=2000.0,
              mu_ra=0.0,
              mu_dec=0.0,
              timediff=False):
    '''
    Converts Julian Date to Baryocentric Julian Date for object at ra, dec
    observed at jd, at an observatory located at obslon, obslat, and at an
    altitude of obsalt.

    epoch = epoch of coordinates, default is J2000.0

    mu_ra = proper motion in RA (mas/yr)

    mu_dec = proper_motion in Dec (mas/yr)

    timediff = True -> return difference in seconds between JD and HJD

    '''

    # if the epoch isn't J2000.0, or if there's any proper motion, get the
    # precessed coordinates for the object
    if (epoch != 2000.0) or (mu_ra != 0.0) or (mu_dec != 0.0):

        raproc, decproc = precess_coordinates(ra, dec, epoch, 2000.0,
                                              mu_ra=mu_ra, mu_dec=mu_dec)

        # FIXME: this makes this function work on scalars only, add vector
        # support later
        raproc, decproc = (float(np.radians(raproc[0])),
                           float(np.radians(decproc[0])))


    else:

        raproc, decproc = np.radians(ra), np.radians(dec)

    cd, sd, ca, sa  = (np.cos(decproc), np.sin(decproc),
                       np.cos(raproc), np.sin(raproc))


    # calculate the seconds elapsed since J2000.0 for jd
    elapsed_seconds = (jd - JD2000)*SEC_P_DAY

    # use SPICE to get the J2000 frame x, y, z coords of the Earth relative to
    # the Solar System Barycenter "SSB".
    (x,y,z,vx,vy,vz), light_time = spice.spkezr('Earth',
                                                elapsed_seconds,
                                                'J2000',
                                                'None',
                                                'SSB')
    # get the J2000 frame xyz coordinates of the observatory relative to the
    # center of the Earth
    if obsalt > -9998:

        _, radii = spice.bodvcd(399,'RADII',3)
        equator, polar = radii[0], radii[2]
        flatcoeff = (equator - polar)/equator
        earthposition = spice.georec(obslon, obslat, obsalt/1000.0,
                                     equator, flatcoeff)
        transformmatrix = spice.pxform('IAU_EARTH',
                                       'J2000',
                                       elapsed_seconds)
        framex, framey, framez = spice.mxv(transformmatrix, earthposition)

        # get the vector from the SSB to observer position
        x, y, z = x + framex, y + framey, z + framez

    # now calculate the Romer delay
    delta_romer = (x*cd*ca + y*cd*sa + z*sd)/CLIGHT_KPS/SEC_P_DAY

    if timediff:
        return delta_romer*SEC_P_DAY
    else:
        return jd + delta_romer
