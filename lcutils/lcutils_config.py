#!/usr/bin/env python

'''
lcutils_settings.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Nov 2013

Contains config variables for lcutils.py.

'''

import ConfigParser
import os.path

LCCACHE = "/Users/jah5/Documents/Fall2014_Gaspar/rrlyr_classification/lccache"

'''
#########################
## DATABASES AND HOSTS ##
#########################

 # parse the configuration file to get the database credentials
CONF_FILE = 'lcserver.conf'

CONF = ConfigParser.ConfigParser()
CONF.read(CONF_FILE)

# database config
DBUSER = CONF.get('database','user')
DBPASS = CONF.get('database','password')
DBDATA = CONF.get('database','database')

# lightcurve cache directory
LCCACHE = CONF.get('paths','lccache')
LCCACHE = os.path.abspath(LCCACHE)

## MYSQL DATABASES ##

# DB config for the HATMASTER DB
MYSQL_USER_HN = CONF.get('mysql','hatmaster_user')
MYSQL_PASS_HN = CONF.get('mysql','hatmaster_password')
MYSQL_DATA_HN = CONF.get('mysql','hatmaster_database')
MYSQL_HOST_HN = CONF.get('mysql','hatmaster_host')

# DB config for the HATSOUTH DB
MYSQL_USER_HS = CONF.get('mysql','hatsouth_user')
MYSQL_PASS_HS = CONF.get('mysql','hatsouth_password')
MYSQL_DATA_HS = CONF.get('mysql','hatsouth_database')
MYSQL_HOST_HS = CONF.get('mysql','hatsouth_host')

# DB config for the HATRED DB
MYSQL_USER_HNCAND = CONF.get('mysql','hatred_user')
MYSQL_PASS_HNCAND = CONF.get('mysql','hatred_password')
MYSQL_DATA_HNCAND = CONF.get('mysql','hatred_database')
MYSQL_HOST_HNCAND = CONF.get('mysql','hatred_host')

# DB config for the HSCAND DB
MYSQL_USER_HSCAND = CONF.get('mysql','hscand_user')
MYSQL_PASS_HSCAND = CONF.get('mysql','hscand_password')
MYSQL_DATA_HSCAND = CONF.get('mysql','hscand_database')
MYSQL_HOST_HSCAND = CONF.get('mysql','hscand_host')


# HOST CONFIG for the HATSouth lightcurve fetching functions
HSLC_HOST = CONF.get('hosts','hatsouth_lcs')
HSLC_HOST = HSLC_HOST.split(',')
HSLC_HOST = [x.strip() for x in HSLC_HOST]

'''
# FALLBACK path for HN LCs to be used when field location isn't found in the
# LCDB
FALLBACK_HN_LCPATH = '/nfs/phn11/ar1/H/BIGPROJ/hatuser/2007_hatnet_phot'
#FALLBACK_SSHUSER = CONF.get('hosts','sshuser')


##########################
## TEXTLC GENERAL CONFIG ##
##########################

TEXTLC_OUTPUT_COLUMNS = {
    'BJD':['time in Baryocentric Julian Date',
           '%20.7f','D',float],
    'MJD':['time in Modified Julian Date',
           '%20.7f','D',float],
    'HJD':['time in Heliocentric Julian Date',
           '%20.7f','D',float],
    'RJD':['time in Reduced Julian Date',
           '%20.7f','D',float],
    'FJD':['time in Full Julian Date',
           '%20.7f','D',float],
    'XCC':['x coordinate on CCD',
           '%.1f','E',float],
    'YCC':['y coordinate on CCD',
           '%.1f','E',float],
    'IM1':['instrumental lightcurve magnitude (aperture 1)',
           '%12.5f','D',float],
    'IE1':['instrumental lightcurve measurement error (aperture 1)',
           '%12.5f','D',float],
    'IQ1':['instrumental lightcurve quality flag (aperture 1)',
           '%s','1A',str],
    'IM2':['instrumental lightcurve magnitude (aperture 2)',
           '%12.5f','D',float],
    'IE2':['instrumental lightcurve measurement error (aperture 2)',
           '%12.5f','D',float],
    'IQ2':['instrumental lightcurve quality flag (aperture 2)',
           '%s','1A',str],
    'IM3':['instrumental lightcurve magnitude (aperture 3)',
           '%12.5f','D',float],
    'IE3':['instrumental lightcurve measurement error (aperture 3)',
           '%12.5f','D',float],
    'IQ3':['instrumental lightcurve quality flag (aperture 3)',
           '%s','1A',str],
    'RM1':['reduced lightcurve magnitude (aperture 1)',
           '%12.5f','D',float],
    'RM2':['reduced lightcurve magnitude (aperture 2)',
           '%12.5f','D',float],
    'RM3':['reduced lightcurve magnitude (aperture 3)',
           '%12.5f','D',float],
    'EP1':['EPD lightcurve magnitude (aperture 1)',
           '%12.5f','D',float],
    'EP2':['EPD lightcurve magnitude (aperture 2)',
           '%12.5f','D',float],
    'EP3':['EPD lightcurve magnitude (aperture 3)',
           '%12.5f','D',float],
    'TF1':['TFA lightcurve magnitude (aperture 1)',
           '%12.5f','D',float],
    'TF2':['TFA lightcurve magnitude (aperture 2)',
           '%12.5f','D',float],
    'TF3':['TFA lightcurve magnitude (aperture 3)',
           '%12.5f','D',float],
    'RSTF':['HAT station and frame number of this LC point',
            '%s','10A',str],
    'RSTFC':['HAT station, frame number, and CCD of this LC point',
             '%s','10A',str],
    'ESTF':['HAT station and frame number of this LC point',
            '%s','10A',str],
    'ESTFC':['HAT station, frame number, and CCD of this LC point',
             '%s','10A',str],
    'TSTF':['HAT station and frame number of this LC point',
            '%s','10A',str],
    'TSTFC':['HAT station, frame number, and CCD of this LC point',
             '%s','10A',str],
    'FSV':['PSF fit S value',
           '%12.5f','E',float],
    'FDV':['PSF fit D value',
           '%12.5f','E',float],
    'FKV':['PSF fit K value',
           '%12.5f','E',float],
    'FLT':['filter used for this LC point',
           '%s','1A',str],
    'FLD':['observed HAT field',
           '%s','15A',str],
    'CCD':['CCD taking this LC point ',
           '%s','I',int],
    'CFN':['CCD frame number',
           '%s','I',int],
    'STF':['HAT station taking this LC point',
           '%s','I',int],
    'BGV':['Background value',
           '%12.5f','E',float],
    'BGE':['Background error',
           '%12.5f','E',float],
    'IHA':['Hour angle of object [hr]',
           '%12.5f','E',float],
    'IZD':['Zenith distance of object [deg]',
           '%12.5f','E',float],
    'NET':['HAT network responsible for this LC point',
           '%s','2A',str],
    'EXP':['exposure time for this LC point [seconds]',
           '%12.3f','E',float],
    'CAM':['camera taking the exposure for this LC point',
           '%s','2A',str],
    'TEL':['telescope taking the exposure for this LC point',
           '%s','2A',str],
    'XIC':['image-subtraction X coordinate on CCD',
           '%.1f','E',float],
    'YIC':['image-subtraction Y coordinate on CCD',
           '%.1f','E',float],
    'IRM1':['image-subtraction lightcurve reduced magnitude (aperture 1)',
           '%12.5f','D',float],
    'IRE1':['image-subtraction lightcurve measurement error (aperture 1)',
           '%12.5f','D',float],
    'IRQ1':['image-subtraction lightcurve quality flag (aperture 1)',
           '%s','1A',str],
    'IRM2':['image-subtraction lightcurve reduced magnitude (aperture 2)',
           '%12.5f','D',float],
    'IRE2':['image-subtraction lightcurve measurement error (aperture 2)',
           '%12.5f','D',float],
    'IRQ2':['image-subtraction lightcurve quality flag (aperture 2)',
           '%s','1A',str],
    'IRM3':['image-subtraction lightcurve reduced magnitude (aperture 3)',
           '%12.5f','D',float],
    'IRE3':['image-subtraction lightcurve measurement error (aperture 3)',
           '%12.5f','D',float],
    'IRQ3':['image-subtraction lightcurve quality flag (aperture 3)',
           '%s','1A',str],
    'IEP1':['image-subtraction EPD lightcurve magnitude (aperture 1)',
           '%12.5f','D',float],
    'IEP2':['image-subtraction EPD lightcurve magnitude (aperture 2)',
           '%12.5f','D',float],
    'IEP3':['image-subtraction EPD lightcurve magnitude (aperture 3)',
           '%12.5f','D',float],
    'ITF1':['image-subtraction TFA lightcurve magnitude (aperture 1)',
           '%12.5f','D',float],
    'ITF2':['image-subtraction TFA lightcurve magnitude (aperture 2)',
           '%12.5f','D',float],
    'ITF3':['image-subtraction TFA lightcurve magnitude (aperture 3)',
           '%12.5f','D',float],
    }

# this notes all the columns in each filetype in order
HATLC_COL_DEFS = {'hn':{'rlc':['HAT',
                               'RSTF',
                               'XCC',
                               'YCC',
                               'BGV',
                               'BGE',
                               'IM1',
                               'IE1',
                               'IQ1',
                               'IM2',
                               'IE2',
                               'IQ2',
                               'IM3',
                               'IE3',
                               'IQ3',
                               'RM1',
                               'RM2',
                               'RM3',
                               'FSV',
                               'FDV',
                               'FKV'],
                        'epdlc':['ESTF',
                                 'RJD',
                                 'IM1',
                                 'IE1',
                                 'IQ1',
                                 'IM2',
                                 'IE2',
                                 'IQ2',
                                 'IM3',
                                 'IE3',
                                 'IQ3',
                                 'RM1',
                                 'RM2',
                                 'RM3',
                                 'EP1',
                                 'EP2',
                                 'EP3'],
                        'tfalc':['TSTF',
                                 'RJD',
                                 'IM1',
                                 'IE1',
                                 'IQ1',
                                 'IM2',
                                 'IE2',
                                 'IQ2',
                                 'IM3',
                                 'IE3',
                                 'IQ3',
                                 'RM1',
                                 'RM2',
                                 'RM3',
                                 'EP1',
                                 'EP2',
                                 'EP3',
                                 'TF1',
                                 'TF2',
                                 'TF3'],
                        'ilc':['RSTF',
                               'RJD',
                               'HAT',
                               'XIC',
                               'YIC',
                               'XCC',
                               'YCC',
                               'FSV',
                               'FDV',
                               'FKV',
                               'BJD',
                               None,
                               None,
                               None,
                               None,
                               None,
                               'IRQ1',
                               'IRM1',
                               'IRE1',
                               'IEP1',
                               'ITF1',
                               'IRQ2',
                               'IRM2',
                               'IRE2',
                               'IEP2',
                               'ITF2',
                               'IRQ3',
                               'IRM3',
                               'IRE3',
                               'IEP3',
                               'ITF3'],
                        'oilc':['RSTF',
                                'HAT',
                                'XIC',
                                'YIC',
                                'XCC',
                                'YCC',
                                'FSV',
                                'FDV',
                                'FKV',
                                'IRQ1',
                                'IRM1',
                                'IRE1',
                                'IRQ2',
                                'IRM2',
                                'IRE2',
                                'IRQ3',
                                'IRM3',
                                'IRE3'],
                        'iepdlc':['ESTF',
                                  'RJD',
                                  'HAT',
                                  'XIC',
                                  'YIC',
                                  'XCC',
                                  'YCC',
                                  'FSV',
                                  'FDV',
                                  'FKV',
                                  'BJD',
                                  None,
                                  None,
                                  None,
                                  None,
                                  None,
                                  'IRQ1',
                                  'IRM1',
                                  'IRE1',
                                  'IEP1',
                                  'ITF1',
                                  'IRQ2',
                                  'IRM2',
                                  'IRE2',
                                  'IEP2',
                                  'ITF2',
                                  'IRQ3',
                                  'IRM3',
                                  'IRE3',
                                  'IEP3',
                                  'ITF3'],
                        'itfalc':['TSTF',
                                  'RJD',
                                  'HAT',
                                  'XIC',
                                  'YIC',
                                  'XCC',
                                  'YCC',
                                  'FSV',
                                  'FDV',
                                  'FKV',
                                  'BJD',
                                  None,
                                  None,
                                  None,
                                  None,
                                  None,
                                  'IRQ1',
                                  'IRM1',
                                  'IRE1',
                                  'IEP1',
                                  'ITF1',
                                  'IRQ2',
                                  'IRM2',
                                  'IRE2',
                                  'IEP2',
                                  'ITF2',
                                  'IRQ3',
                                  'IRM3',
                                  'IRE3',
                                  'IEP3',
                                  'ITF3']},
                  'hs':{'rawlc':['STF', # binary format rawlc
                                 'CFN',
                                 'CCD',
                                 'FLD',
                                 'BJD',
                                 'IM1',
                                 'IE1',
                                 'IQ1',
                                 'IM2',
                                 'IE2',
                                 'IQ2',
                                 'IM3',
                                 'IE3',
                                 'IQ3',
                                 'RM1',
                                 'RM2',
                                 'RM3',
                                 'EP1',
                                 'EP2',
                                 'EP3',
                                 'TF1',
                                 'TF2',
                                 'TF3',
                                 'XCC',
                                 'YCC',
                                 'BGV',
                                 'BGE',
                                 'FSV',
                                 'FDV',
                                 'FKV',
                                 'IHA',
                                 'IZD',
                                 'RJD'],
                        'rlc':['RSTFC', # text format rawlc
                               'FLD',
                               'BJD',
                               'IM1',
                               'IE1',
                               'IQ1',
                               'IM2',
                               'IE2',
                               'IQ2',
                               'IM3',
                               'IE3',
                               'IQ3',
                               'RM1',
                               'RM2',
                               'RM3',
                               'EP1',
                               'EP2',
                               'EP3',
                               'TF1',
                               'TF2',
                               'TF3',
                               'XCC',
                               'YCC',
                               'BGV',
                               'BGE',
                               'FSV',
                               'FDV',
                               'FKV',
                               'IHA',
                               'IZD',
                               'RJD'],
                        'epdlc':['HAT', # text format epdlc
                                 'ESTFC',
                                 'FLD',
                                 'BJD',
                                 'IM1',
                                 'IE1',
                                 'IQ1',
                                 'IM2',
                                 'IE2',
                                 'IQ2',
                                 'IM3',
                                 'IE3',
                                 'IQ3',
                                 'RM1',
                                 'RM2',
                                 'RM3',
                                 'EP1',
                                 'EP2',
                                 'EP3',
                                 'TF1',
                                 'TF2',
                                 'TF3',
                                 'XCC',
                                 'YCC',
                                 'BGV',
                                 'BGE',
                                 'FSV',
                                 'FDV',
                                 'FKV',
                                 'IHA',
                                 'IZD',
                                 'RJD'],
                        'tfalc':['HAT', # text format tfalc
                                 'TSTFC',
                                 'FLD',
                                 'BJD',
                                 'IM1',
                                 'IE1',
                                 'IQ1',
                                 'IM2',
                                 'IE2',
                                 'IQ2',
                                 'IM3',
                                 'IE3',
                                 'IQ3',
                                 'RM1',
                                 'RM2',
                                 'RM3',
                                 'EP1',
                                 'EP2',
                                 'EP3',
                                 'TF1',
                                 'TF2',
                                 'TF3',
                                 'XCC',
                                 'YCC',
                                 'BGV',
                                 'BGE',
                                 'FSV',
                                 'FDV',
                                 'FKV',
                                 'IHA',
                                 'IZD',
                                 'RJD']},
                  }


# maximum line lengths for each type oF LC, required to filter out weird reduced
# LCs
TEXTLC_LINE_LENGTHS = {'rlc':180,
                       'epdlc':150,
                       'tfalc':190,
                       'ilc':280,
                       'iepdlc':280,
                       'itfalc':280,
                       'hsrlc':400,
                       'hsepd':400,
                       'hstfa':400}


TEXTLC_CONSOLIDATION_SOURCES = {'hn':{'XCC':['rlc'],
                                      'YCC':['rlc'],
                                      'IM1':['rlc','epdlc','tfalc'],
                                      'IE1':['rlc','epdlc','tfalc'],
                                      'IQ1':['rlc','epdlc','tfalc'],
                                      'IM2':['rlc','epdlc','tfalc'],
                                      'IE2':['rlc','epdlc','tfalc'],
                                      'IQ2':['rlc','epdlc','tfalc'],
                                      'IM3':['rlc','epdlc','tfalc'],
                                      'IE3':['rlc','epdlc','tfalc'],
                                      'IQ3':['rlc','epdlc','tfalc'],
                                      'RM1':['rlc','epdlc','tfalc'],
                                      'RM2':['rlc','epdlc','tfalc'],
                                      'RM3':['rlc','epdlc','tfalc'],
                                      'EP1':['epdlc','tfalc'],
                                      'EP2':['epdlc','tfalc'],
                                      'EP3':['epdlc','tfalc'],
                                      'TF1':['tfalc'],
                                      'TF2':['tfalc'],
                                      'TF3':['tfalc'],
                                      'RSTF':['rlc'],
                                      'ESTF':['epdlc'],
                                      'TSTF':['tfalc'],
                                      'FSV':['rlc'],
                                      'FDV':['rlc'],
                                      'FKV':['rlc'],
                                      'BGV':['rlc'],
                                      'BGE':['rlc']},
                                'hs':{'XCC':['rawlc'],
                                      'YCC':['rawlc'],
                                      'IM1':['rawlc','epdlc','tfalc'],
                                      'IE1':['rawlc','epdlc','tfalc'],
                                      'IQ1':['rawlc','epdlc','tfalc'],
                                      'IM2':['rawlc','epdlc','tfalc'],
                                      'IE2':['rawlc','epdlc','tfalc'],
                                      'IQ2':['rawlc','epdlc','tfalc'],
                                      'IM3':['rawlc','epdlc','tfalc'],
                                      'IE3':['rawlc','epdlc','tfalc'],
                                      'IQ3':['rawlc','epdlc','tfalc'],
                                      'RM1':['rawlc','epdlc','tfalc'],
                                      'RM2':['rawlc','epdlc','tfalc'],
                                      'RM3':['rawlc','epdlc','tfalc'],
                                      'EP1':['epdlc','tfalc'],
                                      'EP2':['epdlc','tfalc'],
                                      'EP3':['epdlc','tfalc'],
                                      'TF1':['tfalc'],
                                      'TF2':['tfalc'],
                                      'TF3':['tfalc'],
                                      'ESTFC':['epdlc'],
                                      'TSTFC':['tfalc'],
                                      'FSV':['rawlc'],
                                      'FDV':['rawlc'],
                                      'FKV':['rawlc'],
                                      'BGV':['rawlc'],
                                      'BGE':['rawlc'],
                                      'STF':['rawlc'],
                                      'CFN':['rawlc'],
                                      'CCD':['rawlc'],
                                      'FLD':['rawlc'],
                                      'IHA':['rawlc'],
                                      'IZD':['rawlc'],
                                      'RJD':['rawlc']}
                                }

TEXTLC_DATATABLE_COLUMNS = ['RJD',
                            'NET',
                            'STF',
                            'CFN',
                            'CCD',
                            'FLT',
                            'FLD',
                            'EXP',
                            'XCC',
                            'YCC',
                            'BGV',
                            'BGE',
                            'FSV',
                            'FDV',
                            'FKV',
                            'IHA',
                            'IZD',
                            'IM1','IE1','IQ1',
                            'IM2','IE2','IQ2',
                            'IM3','IE3','IQ3',
                            'RM1','RM2','RM3',
                            'EP1','EP2','EP3',
                            'TF1','TF2','TF3',
                            'IRM1','IRE1','IRQ1',
                            'IRM2','IRE2','IRQ2',
                            'IRM3','IRE3','IRQ3',
                            'IEP1','IEP2','IEP3',
                            'ITF1','ITF2','ITF3']


TEXTLC_DEFAULT_OUTPUTCOLS = ['BJD',
                             'RJD',
                             'NET',
                             'STF',
                             'CCD',
                             'FLT',
                             'EXP',
                             'IM1',
                             'IE1',
                             'IQ1',
                             'IM2',
                             'IE2',
                             'IQ2',
                             'IM3',
                             'IE3',
                             'IQ3',
                             'RM1',
                             'RM2',
                             'RM3',
                             'EP1',
                             'EP2',
                             'EP3',
                             'TF1',
                             'TF2',
                             'TF3']


HATLC_OUTPUT_COLUMNS = {'default': TEXTLC_DEFAULT_OUTPUTCOLS,
                        'full': ['BJD','HJD'] + TEXTLC_DATATABLE_COLUMNS,
                        'epdlc':['BJD',
                                 'RJD',
                                 'NET',
                                 'STF',
                                 'CCD',
                                 'FLT',
                                 'EXP',
                                 'RM1',
                                 'RM2',
                                 'RM3',
                                 'EP1',
                                 'EP2',
                                 'EP3'],
                        'tfalc':['BJD',
                                 'RJD',
                                 'NET',
                                 'STF',
                                 'CCD',
                                 'FLT',
                                 'EXP',
                                 'RM1',
                                 'RM2',
                                 'RM3',
                                 'EP1',
                                 'EP2',
                                 'EP3',
                                 'TF1',
                                 'TF2',
                                 'TF3'],
                        'redlc':['BJD',
                                 'RJD',
                                 'NET',
                                 'STF',
                                 'CCD',
                                 'FLT',
                                 'EXP',
                                 'IM1',
                                 'IE1',
                                 'IQ1',
                                 'IM2',
                                 'IE2',
                                 'IQ2',
                                 'IM3',
                                 'IE3',
                                 'IQ3',
                                 'RM1',
                                 'RM2',
                                 'RM3'],
                        'iepdlc':['BJD',
                                  'RJD',
                                  'NET',
                                  'STF',
                                  'CCD',
                                  'FLT',
                                  'EXP',
                                  'RM1',
                                  'RM2',
                                  'RM3',
                                  'EP1',
                                  'EP2',
                                  'EP3',
                                  'IRM1',
                                  'IRM2',
                                  'IRM3',
                                  'IEP1',
                                  'IEP2',
                                  'IEP3'],
                        'itfalc':['BJD',
                                  'RJD',
                                  'NET',
                                  'STF',
                                  'CCD',
                                  'FLT',
                                  'EXP',
                                  'RM1',
                                  'RM2',
                                  'RM3',
                                  'EP1',
                                  'EP2',
                                  'EP3',
                                  'TF1',
                                  'TF2',
                                  'TF3',
                                  'IRM1',
                                  'IRM2',
                                  'IRM3',
                                  'IEP1',
                                  'IEP2',
                                  'IEP3',
                                  'ITF1',
                                  'ITF2',
                                  'ITF3'],
                        'iredlc':['BJD',
                                  'RJD',
                                  'NET',
                                  'STF',
                                  'CCD',
                                  'FLT',
                                  'EXP',
                                  'IM1',
                                  'IE1',
                                  'IQ1',
                                  'IM2',
                                  'IE2',
                                  'IQ2',
                                  'IM3',
                                  'IE3',
                                  'IQ3',
                                  'RM1',
                                  'RM2',
                                  'RM3',
                                  'IRM1',
                                  'IRE1',
                                  'IRQ1',
                                  'IRM2',
                                  'IRE2',
                                  'IRQ2',
                                  'IRM3',
                                  'IRE3',
                                  'IRQ3'],
                        'fullred':['BJD',
                                   'RJD',
                                   'NET',
                                   'RSTF',
                                   'CCD',
                                   'FLT',
                                   'FLD',
                                   'EXP',
                                   'XCC',
                                   'YCC',
                                   'BGV',
                                   'BGE',
                                   'FSV',
                                   'FDV',
                                   'FKV',
                                   'IHA',
                                   'IZD',
                                   'IM1',
                                   'IE1',
                                   'IQ1',
                                   'IM2',
                                   'IE2',
                                   'IQ2',
                                   'IM3',
                                   'IE3',
                                   'IQ3',
                                   'RM1',
                                   'RM2',
                                   'RM3']}



HATNET_FILTER_DESCRIPTIONS = {
    0:['none','no filter'],
    1:['none','no filter'],
    2:['V','Bessel V (Omega Optical)'],
    3:['I','Bessel I (Omega Optical)'],
    4:['I','Bessel I (Omega Optical)'],
    5:['V','Bessel V (Omega Optical)'],
    6:['I','Bessel I (Omega Optical)'],
    7:['I','Bessel I (Omega Optical 0303)'],
    8:['I','Bessel I (Omega Optical)'],
    9:['I','Bessel I (Omega Optical)'],
    10:['V','Bessel V (Omega Optical)'],
    11:['I','Bessel I (Omega Optical)'],
    12:['V','Bessel V (Omega Optical)'],
    13:['I','Bessel I (Omega Optical)'],
    14:['V','Bessel V (Omega Optical)'],
    15:['I','Bessel I (Omega Optical)'],
    16:['V','Bessel V (Omega Optical)'],
    17:['R','Bessel R (Omega Optical)'],
    18:['I','Bessel I (Omega Optical)'],
    19:['R','Bessel R (Omega Optical)'],
    20:['R','Bessel R (Omega Optical)'],
    21:['R','Bessel R (Omega Optical)'],
    22:['R','Bessel R (Omega Optical)'],
    23:['R','Bessel R (Omega Optical)'],
    24:['r','SDSS r SN399 AST0166'],
    25:['r','SDSS r SN403 AST0166'],
    26:['r','SDSS r SN401 AST0166'],
    27:['r','SDSS r SN400 AST0166'],
    28:['r','SDSS r SN402 AST0166'],
    29:['r','SDSS r SN404 AST0166'],
    30:['r','SDSS r SN400 AST0166'],
    31:['i',"Sloan i Astrodon Gen1"],
    32:['z_s2','Sloan z_s2 Astrodon Gen1'],
    }

HATSOUTH_FILTER_DESCRIPTIONS = {
    0:['none','no filter'],
    1:['r','Sloan r AST0285'],
    2:['r','Sloan r AST0285'],
    3:['r','Sloan r AST0285'],
    4:['r','Sloan r AST0285'],
    5:['r','Sloan r AST0285'],
    6:['r','Sloan r AST0285'],
    7:['r','Sloan r AST0285'],
    8:['r','Sloan r AST0285'],
    9:['r','Sloan r SN072 AST0285'],
    10:['r','Sloan r SN068 AST0285'],
    11:['r','Sloan r SN069 AST0285'],
    12:['r','Sloan r SN070 AST0285'],
    13:['r','Sloan r SN071 AST0285'],
    14:['r','Sloan r SN072 AST0285'],
    15:['r','Sloan r SN073 AST0285'],
    16:['r','Sloan r SN074 AST0285'],
    17:['r','Sloan r SN075 AST0285'],
    18:['r','Sloan r SN076 AST0285'],
    19:['r','Sloan r SN077 AST0285'],
    20:['r','Sloan r SN078 AST0285'],
    21:['r','Sloan r SN080 AST0285'],
    22:['r','Sloan r SN079 AST0285'],
    23:['r','Sloan r SN081 AST0285'],
    24:['r','Sloan r SN082 AST0285'],
    25:['r','Sloan r SN083 AST0285'],
    26:['r','Sloan r AST0285'],
    28:['i',"Sloan i' Astrodon Gen2 Round1 AST0285"],
    29:['z_s2','Sloan z_s2 Astrodon Gen2 Square1'],
    101:['r','no description'],
    102:['r','no description'],
    103:['r','no description'],
    104:['r','no description'],
    105:['r','no description'],
    106:['r','no description'],
    107:['r','no description'],
    108:['r','no description'],
    109:['r','no description'],
    }


TEXTLC_HEADER_TEMPLATE = '''\
# {hatid} - 2MASS J{twomassid}
# RA = {ra:.3f} deg, DEC = {dec:.3f} deg
# V = {vmag:.2f}, R = {rmag:.2f}, I = {imag:.2f}, \
J = {jmag:.2f}, H = {hmag:.2f}, K = {kmag:.2f}
#
# Total LC detections: {ndet}
# HAT stations: {hatstations}
#
# Filters used:
#
{filterlist}
#
# Columns:
#
{columnlist}

'''


###########################
## HAT STATION LOCATIONS ##
###########################


#                           lat      lon     alt
HAT_LOCATIONS = {'HN':{5:[31.6811,-110.8783,2320.0],
                       6:[31.6811,-110.8783,2320.0],
                       7:[31.6811,-110.8783,2320.0],
                       8:[19.8244,-155.4733,4160.0],
                       9:[19.8244,-155.4733,4160.0],
                       10:[31.6811,-110.8783,2320.0],
                       11:[30.59583,34.76333,875.0], # WISE HAT
                       12:[31.6811,-110.8783,2320.0]},
                 'HS':{1:[-29.0146,-70.6926,2282.0],
                       2:[-29.0146,-70.6926,2282.0],
                       3:[-23.2716,16.5,1800.0],
                       4:[-23.2716,16.5,1800.0],
                       5:[-31.273333,149.064444,1149.0],
                       6:[-31.273333,149.064444,1149.0]}
                 }




#######################################
## CONFIG FOR FILTER_CONSOLIDATED_LC ##
#######################################

CONDITION_OPERATORS = {'=':"(lc['{col}'] == {operand1})",
                       '!=':"(lc['{col}'] != {operand1})",
                       '>':"(lc['{col}'] > {operand1})",
                       '<':"(lc['{col}'] < {operand1})",
                       '<=':"(lc['{col}'] <= {operand1})",
                       '>=':"(lc['{col}'] >= {operand1})",
                       'between':("((lc['{col}'] > {operand1}) & "
                                  "(lc['{col}'] < {operand2}))")}