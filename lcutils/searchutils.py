#!/usr/bin/env python

'''
searchutils.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Aug 2013

Contains various useful tools for finding and processing HAT LCs.


'''

import logging
import ConfigParser
import os
import os.path
import subprocess
import shlex
import glob
import json
import time
from datetime import datetime

import re
from tornado.escape import squeeze
from lcutils_config import LCCACHE

import numpy.random
from numpy import nan

from textwrap import TextWrapper
TW = TextWrapper()
TW.initial_indent = '#   '
TW.subsequent_indent = '#   '
TW.width = 80
TW.break_long_words = False
os.environ['DYLD_LIBRARY_PATH'] = '/opt/local/lib'
import psycopg2 as pg

import lcdb
import twomass

from coordutils import hms_to_decimal, dms_to_decimal, \
    hms_str_to_tuple, dms_str_to_tuple

from lcaccessdb import auth_hatlcs
from lcauthdb import check_user, check_apikey, get_apikey_username

from zmqutils import datetime_to_jsondate

# setup a logger
LOGGER = logging.getLogger('searchutils')
LOGGER.addHandler(logging.NullHandler())

'''
##########################
### USEFUL DEFINITIONS ###
##########################

# parse the configuration file to get the database credentials

CONF_FILE = 'lcserver.conf'

CONF = ConfigParser.ConfigParser()
CONF.read(CONF_FILE)

# database config
DBUSER = CONF.get('database','user')
DBPASS = CONF.get('database','password')
DBDATA = CONF.get('database','database')
DBHOST = CONF.get('database','host')

# lightcurve temporary cache directory
LCCACHE = CONF.get('paths','lccache')
'''
LCCACHE = os.path.abspath(LCCACHE)

TWOMASS_COLS = {'ra':'%8.3f',
                'decl':'%8.3f',
                'err_maj':'%8.3f','err_min':'%8.3f','err_ang':'%8.3f',
                'twomass_id':'2MASSJ%s',
                'jmag':'%8.3f','jmag_sig':'%8.3f',
                'jmag_sigcom':'%8.3f','jmag_snr':'%8.3f',
                'hmag':'%8.3f','hmag_sig':'%8.3f',
                'hmag_sigcom':'%8.3f','hmag_snr':'%8.3f',
                'kmag':'%8.3f','kmag_sig':'%8.3f',
                'kmag_sigcom':'%8.3f','kmag_snr':'%8.3f',
                'qual_flag':'%3s','rd_flag':'%8s','bl_flag':'%8s',
                'cc_flag':'%8s','ndet_flag':'%8s',
                'nbr_prox':'%8.3f','nbr_pxpa':'%8.3f','nbr_ptskey':'%8s',
                'gal_contam':'%8s','mp_flag':'%8s',
                'ptskey':'%8s','hemis':'%8s','obsdate':'%8s','scan':'%8s',
                'glon':'%8.3f','glat':'%8.3f',
                'x_scan':'%8s','obs_juliandate':'%12.3f',
                'jmag_psfchi':'%8.3f','jmag_stdap':'%8.3f',
                'jmag_sig_stdap':'%8.3f',
                'hmag_psfchi':'%8.3f','hmag_stdap':'%8.3f',
                'hmag_sig_stdap':'%8.3f',
                'kmag_psfchi':'%8.3f','kmag_stdap':'%8.3f',
                'kmag_sig_stdap':'%8.3f',
                'dist_edge_ns':'%8s','dist_edge_ew':'%8s','dist_edge_flg':'%8s',
                'dup_src':'%8s','use_src':'%8s','a':'%8s',
                'dist_opt':'%8s','phi_opt':'%8s','b_m_opt':'%8s',
                'vr_m_opt':'%8s','nopt_mchs':'%8s',
                'ext_key':'%8s','scan_key':'%8s','coadd_key':'%8s','coadd':'%8s',
                'hat_field':'%3i','hat_field_objid':'%07i',
                'hat_id':'%15s',
                'bmag':'%8.3f',
                'vmag':'%8.3f',
                'rmag':'%8.3f',
                'imag':'%8.3f',
                'sdssu':'%8.3f',
                'sdssg':'%8.3f',
                'sdssr':'%8.3f',
                'sdssi':'%8.3f',
                'sdssz':'%8.3f',
                'jh_color':'%8.3f',
                'jk_color':'%8.3f',
                'hk_color':'%8.3f',
                'dist_arcsec':'%8.3f',
                }

COND_OPERATORS = {'=':'({col} = %s)',
                  '!=':'({col} != %s)',
                  '>':'({col} > %s)',
                  '<':'({col} < %s)',
                  '<=':'({col} <= %s)',
                  '>=':'({col} >= %s)',
                  'between':'({col} between %s and %s)'}



############################################
## REGEXES FOR PARSING QUICKSEARCH PARAMS ##
############################################

SDSSID_REGEX = re.compile(
    r'^([SDSS J|SDSSJ|J])*(\d{6}\.\d{2})([\+\-]\d{6}\.\d{1})$'
    )
TWOMASSID_REGEX = re.compile(r'^([2MASS J|2MASSJ|J])*(\d{8})([\+\-]\d{7})$')

COORD_SEXAGESIMAL_REGEX = re.compile(
    r'^(\d{1,2}[: ]\d{2}[: ]\d{2}\.{0,1}\d*) '
    '([+\-]\d{1,2}[: ]\d{2}[: ]\d{2}\.{0,1}\d*)$'
    )
COORD_DECIMAL_REGEX = re.compile(
    r'^(\d{1,3}\.{0,1}\d*) ([+\-]\d{1,2}\.{0,1}\d*)$'
    )

HATID_REGEX = re.compile(r'^(HAT{1}\-\d{3}\-\d{7})$')

HATFIELD_REGEX = re.compile(r'^(\d{3})$')
HNCAND_REGEX = re.compile(r'^(HTR{1}\d{3}\-\d{3})$')
HNPLANET_REGEX = re.compile(r'^(HAT\-P{1}\-\d{1,3}[b-z]{1})$')
HSCAND_REGEX = re.compile(r'^(HATS{1}\d{3}\-\d{3})$')
HSPLANET_REGEX = re.compile(r'^(HATS{1}\-\d{1,3}[b-z]{1})$')

HATIDLIST_REGEX = re.compile(r'HAT{1}\-\d{3}\-\d{7}')
COORD_SEARCH_REGEX = re.compile(
    r'^(\d{1,3}\.{0,1}\d*) ([+\-]\d{1,2}\.{0,1}\d*) (\d{1,2}\.{0,1}\d*)$'
    )


####################################
## DATASET COLLECTIONS FOR HATLCS ##
####################################

COLLECTIONS = {
    'hnplanets':{
        'dataset_name':'Exoplanets discovered by the HATNet survey',
        'selector_required':False,
        'selector_regex':HNPLANET_REGEX,
        'selector_default':'HAT-P-',
        'search_query':(
            "select ('HAT-' || to_char(a.hat_field,'FM000') "
            "|| '-' || to_char(a.hat_field_objid, 'FM0000000')) as hat_id, "
            "c.ndet, a.object_type, a.candidate_id, a.final_id, "
            "c.access_groups, "
            "b.ucac4_ra, b.ucac4_decl, b.jmag, b.hmag, b.kmag, "
            "a.ads_ref, a.arxiv_ref, "
            "substring(final_id from '[0-9]{1,3}')::integer as planetnum "
            "from interesting_objects a join twomass_ucac4 b "
            "on ((a.hat_field = b.hat_field) and "
            "(a.hat_field_objid = b.hat_field_objid)) "
            "left join hat_lightcurves c on ((a.hat_field = c.hat_field) and "
            "(a.hat_field_objid = c.hat_field_objid)) "
            "where a.final_id like %s and a.access_groups "
            "like %s order by planetnum desc"
        ),
        'query_params':('selector','access_group'),
        'query_results':('hatid','ndet',
                         'object_type','candidate_id','final_id',
                         'access_groups',
                         'ra','decl','jmag','hmag','kmag',
                         'ads_ref','arxiv_ref'),
    },
    'hsplanets':{
        'dataset_name':'Exoplanets discovered by the HATSouth survey',
        'selector_required':False,
        'selector_regex':HSPLANET_REGEX,
        'selector_default':'HATS-',
        'search_query':(
            "select ('HAT-' || to_char(a.hat_field,'FM000') "
            "|| '-' || to_char(a.hat_field_objid, 'FM0000000')) as hat_id, "
            "c.ndet, a.object_type, a.candidate_id, a.final_id, "
            "c.access_groups, "
            "b.ucac4_ra, b.ucac4_decl, b.jmag, b.hmag, b.kmag, "
            "a.ads_ref, a.arxiv_ref, "
            "substring(final_id from '[0-9]{1,3}')::integer as planetnum "
            "from interesting_objects a join twomass_ucac4 b "
            "on ((a.hat_field = b.hat_field) and "
            "(a.hat_field_objid = b.hat_field_objid)) "
            "left join hat_lightcurves c on ((a.hat_field = c.hat_field) and "
            "(a.hat_field_objid = c.hat_field_objid)) "
            "where a.final_id like %s and a.access_groups "
            "like %s order by planetnum desc"
        ),
        'query_params':('selector','access_group'),
        'query_results':('hatid','ndet',
                         'object_type','candidate_id','final_id',
                         'access_groups',
                         'ra','decl','jmag','hmag','kmag',
                         'ads_ref','arxiv_ref'),
    },
    'hncandidates':{
        'dataset_name':'HATNet exoplanet transit candidates',
        'selector_required':False,
        'selector_regex':HNCAND_REGEX,
        'selector_default':'HTR',
        'search_query':(
            "select ('HAT-' || to_char(a.hat_field,'FM000') "
            "|| '-' || to_char(a.hat_field_objid, 'FM0000000')) as hat_id, "
            "c.ndet, a.object_type, a.candidate_id, a.final_id, "
            "c.access_groups, "
            "b.ucac4_ra, b.ucac4_decl, b.jmag, b.hmag, b.kmag, "
            "a.ads_ref, a.arxiv_ref "
            "from interesting_objects a join twomass_ucac4 b "
            "on ((a.hat_field = b.hat_field) and "
            "(a.hat_field_objid = b.hat_field_objid)) "
            "left join hat_lightcurves c on ((a.hat_field = c.hat_field) and "
            "(a.hat_field_objid = c.hat_field_objid)) "
            "where a.candidate_id like %s and a.access_groups "
            "like %s order by a.candidate_id asc"
        ),
        'query_params':('selector','access_group'),
        'query_results':('hatid','ndet',
                         'object_type','candidate_id','final_id',
                         'access_groups',
                         'ra','decl','jmag','hmag','kmag',
                         'ads_ref','arxiv_ref'),
    },
    'hscandidates':{
        'dataset_name':'HATSouth exoplanet transit candidates',
        'selector_required':False,
        'selector_regex':HSCAND_REGEX,
        'selector_default':'HATS',
        'search_query':(
            "select ('HAT-' || to_char(a.hat_field,'FM000') "
            "|| '-' || to_char(a.hat_field_objid, 'FM0000000')) as hat_id, "
            "c.ndet, a.object_type, a.candidate_id, a.final_id, "
            "c.access_groups, "
            "b.ucac4_ra, b.ucac4_decl, b.jmag, b.hmag, b.kmag, "
            "a.ads_ref, a.arxiv_ref "
            "from interesting_objects a join twomass_ucac4 b "
            "on ((a.hat_field = b.hat_field) and "
            "(a.hat_field_objid = b.hat_field_objid)) "
            "left join hat_lightcurves c on ((a.hat_field = c.hat_field) and "
            "(a.hat_field_objid = c.hat_field_objid)) "
            "where a.candidate_id like %s and a.access_groups "
            "like %s order by a.candidate_id asc"
        ),
        'query_params':('selector','access_group'),
        'query_results':('hatid','ndet',
                         'object_type','candidate_id','final_id',
                         'access_groups',
                         'ra','decl','jmag','hmag','kmag',
                         'ads_ref','arxiv_ref'),
    },
    'field':{
        'dataset_name':'Objects with HAT LCs in HAT field %s',
        'selector_required':True,
        'selector_regex':HATFIELD_REGEX,
        'selector_default':None,
        'search_query':(
            "select a.hat_id, a.ndet, "
            "b.object_type, b.candidate_id, b.final_id, "
            "a.access_groups,"
            "a.ra, a.decl, a.jmag, a.hmag, a.kmag, "
            "a.ads_ref, a.arxiv_ref "
            "from hat_lightcurves a join interesting_objects b on "
            "((a.hat_field = b.hat_field) and "
            "(a.hat_field_objid = b.hat_field_objid)) "
            "where a.hat_field = %s and a.ndet > 0 and a.access_groups like %s "
            "order by a.hat_id fetch first %s rows only"
        ),
        'query_params':('selector','access_group', 'row_limit'),
        'query_results':('hatid','ndet',
                         'object_type','candidate_id','final_id',
                         'access_groups',
                         'ra','decl','jmag','hmag','kmag',
                         'ads_ref','arxiv_ref'),
    },
    'region':{  ## TODO: put in the correct spatial query here
        'dataset_name':'Objects with HAT LCs at %s within %s arcminutes',
        'selector_required':True,
        'selector_regex':COORD_SEARCH_REGEX,
        'selector_default':None,
        'search_query':(
            "select a.hat_id, a.ndet, "
            "b.object_type, b.candidate_id, b.final_id, "
            "a.access_groups,"
            "a.ra, a.decl, a.jmag, a.hmag, a.kmag, "
            "a.ads_ref, a.arxiv_ref "
            "from hat_lightcurves a join interesting_objects b on "
            "((a.hat_field = b.hat_field) and "
            "(a.hat_field_objid = b.hat_field_objid)) "
            "where a.hat_field = %s and a.ndet > 0 and a.access_groups like %s "
            "order by a.hat_id fetch first %s rows only"
        ),
        'query_params':('selector','access_group', 'row_limit', 'radius_limit'),
        'query_results':('hatid','ndet',
                         'object_type','candidate_id','final_id',
                         'access_groups',
                         'ra','decl','jmag','hmag','kmag',
                         'ads_ref','arxiv_ref'),
    },
    'hatids':{
        'dataset_name':'Requested HAT-IDs with HAT LCs',
        'selector_required':True,
        'selector_regex':HATID_REGEX,
        'selector_default':None,
        'search_query':(
            "select a.hat_id, a.ndet, "
            "b.object_type, b.candidate_id, b.final_id, "
            "a.access_groups,"
            "a.ra, a.decl, a.jmag, a.hmag, a.kmag, "
            "a.ads_ref, a.arxiv_ref "
            "from hat_lightcurves a join interesting_objects b on "
            "((a.hat_field = b.hat_field) and "
            "(a.hat_field_objid = b.hat_field_objid)) "
            "where a.hat_id in %s and a.ndet > 0 and a.access_groups like %s "
            "order by a.hat_id fetch first %s rows only"
        ),
        'query_params':('selector','access_group', 'row_limit'),
        'query_results':('hatid','ndet',
                         'object_type','candidate_id','final_id',
                         'access_groups',
                         'ra','decl','jmag','hmag','kmag',
                         'ads_ref','arxiv_ref'),
    },
}




###############################
### 2MASS UTILITY FUNCTIONS ###
###############################

def add_additional_mag_cols(headers, results, additional_mags):
    '''
    This adds the columns listed in additional_mags to results and returns
    the amended results.

    '''

    # find the j, h, kmag cols
    cols = zip(*results)
    jmag = cols[headers.index('jmag')]
    hmag = cols[headers.index('hmag')]
    kmag = cols[headers.index('kmag')]

    outlist = []

    for mag in additional_mags:
        outlist.append(
            [twomass.CONV_FUNCMAP[mag](x,y,z) for (x,y,z) in
             zip(jmag, hmag, kmag)]
            )

    return outlist


def search_results_text(headers,
                        results,
                        orderby,
                        sortorder,
                        querystr,
                        form='ssv'):
    '''
    This returns a text table from the result rows.

    form is either 'ssv' or 'csv' for space-separated values or
    comma-separated values respectively

    '''
    colheaders = headers.split(',')

    outcol_formats = []
    outstr = []
    outheader = ['#']

    # give some information on the results and query used
    metainfo_str = '# %s rows ordered by %s in %s order' % (len(results),
                                                            orderby,
                                                            sortorder)
    outstr.append(metainfo_str)
    query_str = TW.fill(querystr)
    outstr.extend(['#',
                   '# query used:',
                   query_str,
                   '#',
                   '# %s column(s)' % len(colheaders),
                   '#'])


    for colheader in colheaders:

        if colheader in TWOMASS_COLS:

            outcol_formats.append(TWOMASS_COLS[colheader])
            outheader.append(colheader)

        else:

            if LOGGER:
                LOGGER.error('column %s is not present in the twomass table' %
                             colheader)
            else:
                print('searchutils: column %s is not present '
                      'in the twomass table' % colheader)

            return None

    # make it so lines can be modified when we do a sweep
    results = [list(x) for x in results]

    # sweep the results for Nones and convert them to numpy.nans
    for line in results:
        for ind, item in enumerate(line):
            if item is None and 'f' in outcol_formats[ind]:
                line[ind] = nan

    # put together the header and put it in the outstr list
    if form == 'ssv':
        outstr.append(' '.join(outheader))
        rowform = ' '.join(outcol_formats)
    elif form == 'csv':
        headerstr = ','.join(outheader)
        headerstr = headerstr.replace('#,','# ')
        outstr.append(headerstr)
        rowform = ','.join(outcol_formats)

    # put together the rest of the outstr list
    for resultrow in results:

        if form == 'ssv':
            rowstr = rowform % tuple(resultrow)
            outstr.append(rowstr)
        elif form == 'csv':
            rowstr = rowform % tuple(resultrow)
            rowstr = rowstr.replace(' ','')
            outstr.append(rowstr.strip())

    # return the final string
    return '\n'.join(outstr)


def search_results_csv(headers,
                       results,
                       orderby,
                       sortorder,
                       querystr):
    '''
    This just uses search_results_text above in csv mode.

    '''
    return search_results_text(headers,
                               results,
                               orderby,
                               sortorder,
                               querystr,
                               form='csv')


def search_results_json(headers,
                        results,
                        orderby,
                        sortorder,
                        querystr):
    '''
    This returns a JSON object from the result rows.

    '''
    colheaders = headers.split(',')

    nrows = len(results)

    # transpose the results
    resultcols = zip(*results)

    jsondict = {"nresults":nrows,
                "columns":headers,
                "sortedby":orderby,
                "sortorder":sortorder,
                "query":querystr}

    # populate the jsondict
    for header, col in zip(colheaders, resultcols):

        formattedcol = [(TWOMASS_COLS[header] % x).strip() for x in col]
        jsondict[header] = formattedcol

    return json.dumps(jsondict,ensure_ascii=True)


def search_results_dict(headers,
                        results,
                        orderby,
                        sortorder,
                        querystr):
    '''
    This returns a dict from the result rows.

    '''
    colheaders = headers.split(',')

    nrows = len(results)

    # transpose the results
    resultcols = zip(*results)

    returndict = {"nresults":nrows,
                  "columns":headers,
                  "sortedby":orderby,
                  "sortorder":sortorder,
                  "query":querystr}

    # populate the jsondict
    for header, col in zip(colheaders, resultcols):
        returndict[header] = col

    return returndict


def search_results_dict_to_text(results_dict, output_type='text'):
    '''
    This converts a dict returned by search_results_dict to another text format.

    '''
    cols = results_dict['columns']

    results = []

    for col in cols.split(','):
        results.append(results_dict[col])

    results = zip(*results)

    if output_type == 'text':
        form = 'ssv'
    elif output_type == 'csv':
        form = 'csv'

    return search_results_text(cols,
                               results,
                               results_dict['sortedby'],
                               results_dict['sortorder'],
                               results_dict['query'],
                               form=form)



RESULTS_FUNCMAP = {'text':search_results_text,
                   'csv':search_results_csv,
                   'json':search_results_json,
                   'dict':search_results_dict}


##############################
### QUICK-SEARCH FUNCTIONS ###
##############################

def parse_quicksearch_params(search_params,
                             search_radius=300.0,
                             pprint='json',
                             accessuser=None,
                             accessapikey=None):
    '''
    This parses the input to the quicksearch function and returns a set of
    arguments for the twomass_info function to use.

    FIXME: add a quicksearch for arbitrary star tags. This search is run against
    the interesting_objects table. e.g.:

    /lightcurves/search/quick?find=candidate%20planet%20AND%20hat%20planet

    We'll use this to add hyperlinks to the tag labels. The search is carried
    out in the user's session context, so has access to only what they have
    access to. This should then lead to another dataset browser perhaps?

    '''

    # sanitize the search params
    search_params = squeeze(search_params.strip())

    # try all the regexes in series, and use whichever one succeeds first
    twomass_try = re.match(TWOMASSID_REGEX, search_params)
    sdss_try = re.match(SDSSID_REGEX, search_params)
    hatid_try = re.match(HATID_REGEX, search_params)

    coords_sexagesimal_try = re.match(COORD_SEXAGESIMAL_REGEX, search_params)
    coords_decimal_try = re.match(COORD_DECIMAL_REGEX, search_params)

    hncand_try = re.match(HNCAND_REGEX, search_params)
    hscand_try = re.match(HSCAND_REGEX, search_params)
    hnplanet_try = re.match(HNPLANET_REGEX, search_params)
    hsplanet_try = re.match(HSPLANET_REGEX, search_params)

    if twomass_try:

        quicksearch_parsed = True
        quicksearch_function = twomass_info
        quicksearch_func_args = ''.join(twomass_try.groups()[1:])
        quicksearch_func_kwargs = {'identifier_type':'twomass'}
        quicksearch_message = 'params_parsed_twomass'
        quicksearch_func_kwargs['pprint'] = pprint
        quicksearch_type = 'info'

    elif hatid_try:

        quicksearch_parsed = True
        quicksearch_function = twomass_info
        quicksearch_func_args = ''.join(hatid_try.groups()[0])
        quicksearch_func_kwargs = {'identifier_type':'hatid'}
        quicksearch_message = 'params_parsed_hatid'
        quicksearch_func_kwargs['pprint'] = pprint
        quicksearch_type = 'info'

    elif hncand_try:

        quicksearch_parsed = True
        quicksearch_function = object_info
        quicksearch_func_args = ''.join(hncand_try.groups()[0])
        quicksearch_func_kwargs = dict()
        quicksearch_message = 'params_parsed_hncand'
        quicksearch_type = 'info'

    elif hscand_try:

        quicksearch_parsed = True
        quicksearch_function = object_info
        quicksearch_func_args = ''.join(hscand_try.groups()[0])
        quicksearch_func_kwargs = dict()
        quicksearch_message = 'params_parsed_hscand'
        quicksearch_type = 'info'

    elif hnplanet_try:

        quicksearch_parsed = True
        quicksearch_function = object_info
        quicksearch_func_args = ''.join(hnplanet_try.groups()[0])
        quicksearch_func_kwargs = dict()
        quicksearch_message = 'params_parsed_hnplanet'
        quicksearch_type = 'info'

    elif hsplanet_try:

        quicksearch_parsed = True
        quicksearch_function = object_info
        quicksearch_func_args = ''.join(hsplanet_try.groups()[0])
        quicksearch_func_kwargs = dict()
        quicksearch_message = 'params_parsed_hsplanet'
        quicksearch_type = 'info'

    elif coords_sexagesimal_try:

        # first, parse the coordinates

        ra, dec = coords_sexagesimal_try.groups()
        ra_tuple, dec_tuple = hms_str_to_tuple(ra), dms_str_to_tuple(dec)

        ra_hr, ra_min, ra_sec = ra_tuple
        dec_sign, dec_deg, dec_min, dec_sec = dec_tuple

        # make sure the coordinates are all legit
        if ((0 <= ra_hr < 24) and
            (0 <= ra_min < 60) and
            (0 <= ra_sec < 60) and
            (0 <= dec_deg < 90) and
            (0 <= dec_min < 60) and
            (0 <= dec_sec < 60)):

            ra_decimal = hms_to_decimal(ra_hr, ra_min, ra_sec)
            dec_decimal = dms_to_decimal(dec_sign, dec_deg, dec_min, dec_sec)

            quicksearch_parsed = True
            quicksearch_function = twomass_objects_near
            quicksearch_func_args = (ra_decimal, dec_decimal, search_radius)
            quicksearch_func_kwargs = {'twomassids':True}
            quicksearch_message = 'params_parsed_sexagcoords'
            quicksearch_func_kwargs['pprint'] = pprint
            quicksearch_func_kwargs['accessuser'] = accessuser
            quicksearch_func_kwargs['accessapikey'] = accessapikey
            quicksearch_type = 'coords'

        else:

            quicksearch_parsed = False
            quicksearch_function = None
            quicksearch_func_args = None
            quicksearch_func_kwargs = None
            quicksearch_message = 'params_invalid_sexagcoords'
            quicksearch_type = 'coords'


    elif coords_decimal_try:

        ra, dec = coords_decimal_try.groups()

        try:

            ra, dec = float(ra), float(dec)

            if ((abs(ra) < 360.0) and (abs(dec) < 90.0)):

                if ra < 0:
                    ra = 360.0 + ra

                quicksearch_parsed = True
                quicksearch_function = twomass_objects_near
                quicksearch_func_args = (ra, dec, search_radius)
                quicksearch_func_kwargs = {'twomassids':True}
                quicksearch_message = 'params_parsed_decimcoords'
                quicksearch_func_kwargs['pprint'] = pprint
                quicksearch_func_kwargs['accessuser'] = accessuser
                quicksearch_func_kwargs['accessapikey'] = accessapikey
                quicksearch_type = 'coords'

            else:

                quicksearch_parsed = False
                quicksearch_function = None
                quicksearch_func_args = None
                quicksearch_func_kwargs = None
                quicksearch_message = 'params_invalid_decimcoords'
                quicksearch_type = 'coords'


        except Exception as e:

            LOGGER.warning('could not parse quicksearch '
                           'parameters: %s, error was: %s' %
                           (quicksearch_params, e))

            quicksearch_parsed = False
            quicksearch_function = None
            quicksearch_func_args = None
            quicksearch_func_kwargs = None
            quicksearch_message = 'params_invalid_decimcoords'
            quicksearch_type = 'coords'


    elif sdss_try:

        quicksearch_parsed = True
        quicksearch_function = twomass_info

        # get the SDSS ID, but remove the periods and turn it into a 2MASS ID
        quicksearch_func_args = ''.join(sdss_try.groups()[1:])
        quicksearch_func_args = quicksearch_func_args.replace('.','')

        quicksearch_func_kwargs = {'identifier_type':'twomass'}
        quicksearch_message = 'params_parsed_sdss'
        quicksearch_func_kwargs['pprint'] = pprint
        quicksearch_type = 'info'

    else:

        quicksearch_parsed = False
        quicksearch_function = None
        quicksearch_func_args = None
        quicksearch_func_kwargs = None
        quicksearch_message = 'params_invalid'
        quicksearch_type = 'none'

    return (quicksearch_parsed,
            quicksearch_function,
            quicksearch_func_args,
            quicksearch_func_kwargs,
            quicksearch_type,
            quicksearch_message)



def quicksearch(search_params,
                pprint='dict',
                search_radius=300.0,
                accessuser=None,
                accessapikey=None,
                database=None,
                verbose=False):
    '''
    This automatically figures out the type of the input and does a search in
    the twomass database for the corresponding object.

    We have three types of input:

    coordinates
    HAT ID
    2MASS ID

    coordinates:

    - accepted formats:  HH:MM:SS.sss... +/-DD:MM:SS.sss...
                         HH:MM:SS +/-DD:MM:SS
                         HH MM SS.sss... +/-DD MM SS.sss...
                         HH MM SS +/-DD MM SS
                         [D]DD.ddd... +/-DD.ddd...
                         [D]DD +/-DD

    HAT ID:

    - triggers on a HAT-XXX-YYYYYYY identifier

    2MASS ID:

    - triggers on a 2MASS JXXXXXXXX+/-YYYYYYY or JXXXXXXXX+/-YYYYYYY or
      XXXXXXXX+/-YYYYYYY identifier

    SDSS ID:

    - triggers on a SDSS JXXXXXX.XX+/-YYYYYY.Y or JXXXXXX.XX+/-YYYYYY.Y or
      XXXXXX.XX+/-YYYYYY.Y identifier

    This function returns the following (in lists if coordinate search):

    - HAT ID
    - 2MASS ID
    - coordinates
    - mags in B, V, R, I, J, H, K, u, g, r, i, z

    '''

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    # first, parse the parameters
    (parsed_ok, searchfunc, searchfunc_args,
     searchfunc_kwargs, search_type, parsed_msg) = (
        parse_quicksearch_params(search_params,
                                 pprint=pprint,
                                 search_radius=search_radius,
                                 accessuser=accessuser,
                                 accessapikey=accessapikey)
        )


    if parsed_ok:

        # execute the function with generated args and kwargs
        if not isinstance(searchfunc_args, tuple):
            searchfunc_args = (searchfunc_args,)

        searchfunc_kwargs['database'] = database

        if searchfunc is twomass_info:
            searchfunc_kwargs['verbose'] = verbose

        if ('accessuser' not in searchfunc_kwargs and
            searchfunc is not twomass_info):
            searchfunc_kwargs['accessuser'] = accessuser
        if ('accessapikey' not in searchfunc_kwargs and
            searchfunc is not twomass_info):
            searchfunc_kwargs['accessapikey'] = accessapikey

        if searchfunc is object_info and 'identifier_type' in searchfunc_kwargs:
            del searchfunc_kwargs['identifier_type']

        searchresults = searchfunc(*searchfunc_args, **searchfunc_kwargs)

        # if we use the object_info function, we have to reform the dict to make
        # it look like the one for twomass_info
        if searchfunc is object_info:

            resultdict = searchresults[1]

            if resultdict:

                for key in ('candidate_id',
                            'decl_pm',
                            'decl_pm_err',
                            'ra_pm',
                            'ra_pm_err',
                            'final_id',
                            'hat_lastupdated',
                            'hat_lcchecked',
                            'hat_lcfpath',
                            'hat_ndet',
                            'hat_stations',
                            'object_period',
                            'object_epoch',
                            'object_type',
                            'object_accessgroups',
                            'hatlc_accessgroups',
                            'ra_err',
                            'decl_err',
                            'sdssu_err',
                            'sdssg_err',
                            'sdssr_err',
                            'sdssi_err',
                            'sdssz_err',
                            'bmag_err',
                            'vmag_err',
                            'rmag_err',
                            'imag_err',
                            'user_accessgroup',):
                    del resultdict[key]

                for key in ('bmag',
                            'vmag',
                            'rmag',
                            'imag',
                            'sdssu',
                            'sdssg',
                            'sdssr',
                            'sdssi',
                            'sdssz',
                            'jmag',
                            'jmag_err',
                            'hmag',
                            'hmag_err',
                            'kmag',
                            'kmag_err',
                            'qual_flag',
                            'ra',
                            'decl',
                            'hatid',
                            'twomassid'):
                    resultdict[key] = (resultdict[key],)

                resultdict['hat_id'] = resultdict['hatid']
                resultdict['twomass_id'] = resultdict['twomassid']
                resultdict['qual_flag'] = resultdict['qual_flag']
                resultdict['jmag_sig'] = resultdict['jmag_err']
                resultdict['hmag_sig'] = resultdict['hmag_err']
                resultdict['kmag_sig'] = resultdict['kmag_err']
                resultdict['nresults'] = 1

            else:

                resultdict = {'nresults': 0,
                              'query':'find: %s' % search_params}

            resultdict['sortedby'] = 'hatid'
            resultdict['sortorder'] = 'asc'
            resultdict['columns'] = (
                'hat_id,ra,decl,'
                'jmag,jmag_sig,hmag,hmag_sig,kmag,kmag_sig,'
                'qual_flag,twomass_id,bmag,vmag,rmag,imag,'
                'sdssu,sdssg,sdssr,sdssi,sdssz'
                )
            return_tuple = (True, resultdict,
                            search_type, 'quicksearch_ok')


        # otherwise, if the searchfunc is twomass_info, we return its results
        # directly
        else:
            return_tuple = (True, searchresults, search_type, 'quicksearch_ok')

    else:

        return_tuple = (False, None, search_type, parsed_msg)

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple


######################################
### FUNCTION TO RETURN OBJECT INFO ###
######################################

def object_info(identifier,
                accessapikey=None,
                accessuser=None,
                database=None):
    '''
    This function returns all the information needed for the object info pages
    (/lightcurves/object/<identifier>).

    Returns a 4-elem tuple compatible with serviceutils functions.

    task_result is a dictionary with the following keys:

    hatid
    twomass_id
    ucac4_ra + err if in UCAC4 else twomass_ra
    ucac4_decl + err if in UCAC4 else twomass_decl
    ra_propermotion + err
    decl_propermotion + err
    jmag + err
    hmag + err
    kmag + err
    bmag + err from UCAC4 if available, else calculated
    vmag + err from UCAC4 if available, else calculated
    rmag calculated
    imag calculated
    sdssg + err from UCAC4 if available, else calculated
    sdssr + err from UCAC4 if available, else calculated
    sdssi + err from UCAC4 if available, else calculated
    note indicating if we have an LC for this object

    '''

    twomass_try = re.match(TWOMASSID_REGEX, identifier)
    sdss_try = re.match(SDSSID_REGEX, identifier)

    hatid_try = re.match(HATID_REGEX, identifier)
    hncand_try = re.match(HNCAND_REGEX, identifier)
    hscand_try = re.match(HSCAND_REGEX, identifier)
    hnplanet_try = re.match(HNPLANET_REGEX, identifier)
    hsplanet_try = re.match(HSPLANET_REGEX, identifier)

    # make sure the identifier is good
    if not any([twomass_try, sdss_try,
                hatid_try, hncand_try, hscand_try,
                hnplanet_try, hsplanet_try]):
        return (False, None, 'objectinfo', 'objectinfo_invalid_id')

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True


    # if the identifier is a twomassid or hatid
    if twomass_try or hatid_try:

        if twomass_try:
            identifier_type = 'twomass'
        elif hatid_try:
            identifier_type = 'hatid'

        # first, get this object's basic info from the 2MASS table
        obj_twomass_info = twomass_info(
            identifier,
            identifier_type=identifier_type,
            cols=('hat_id,ra,decl,'
                  'jmag,jmag_sig,hmag,hmag_sig,kmag,kmag_sig,'
                  'qual_flag,twomass_id,'
                  'bmag,vmag,rmag,imag,'
                  'sdssu,sdssg,sdssr,sdssi,sdssz'),
            database=database
            )

    # otherwise, it needs to be matched to the final_id in the
    # interesting_objects table
    else:

        if hncand_try and not any([hscand_try,
                                   hnplanet_try,
                                   hsplanet_try]):
            identifier_type = 'hncand'
        elif hscand_try and not any([hncand_try,
                                     hnplanet_try,
                                     hsplanet_try]):
            identifier_type = 'hscand'
        elif hnplanet_try and not any([hncand_try,
                                       hscand_try,
                                       hsplanet_try]):
            identifier_type = 'hnplanet'
        elif hsplanet_try and not any([hncand_try,
                                       hscand_try,
                                       hnplanet_try]):
            identifier_type = 'hsplanet'

        else:
            return (False, None, 'objectinfo', 'invalid_object')


        query = ("select hat_field, hat_field_objid from "
                 "interesting_objects where "
                 "final_id like %s or candidate_id like %s")
        query_params = ('%{identifier}%'.format(identifier=identifier),
                        '%{identifier}%'.format(identifier=identifier),)
        cur.execute(query, query_params)
        results = cur.fetchone()

        # if we get a hatid back, then we're in business, get the twomass info
        # for this object next
        if results and len(results) > 0:
            object_hatid = 'HAT-%s-%s' % results
            obj_twomass_info = twomass_info(
                object_hatid,
                identifier_type='hatid',
                cols=('hat_id,ra,decl,'
                      'jmag,jmag_sig,hmag,hmag_sig,kmag,kmag_sig,'
                      'qual_flag,twomass_id,'
                      'bmag,vmag,rmag,imag,'
                      'sdssu,sdssg,sdssr,sdssi,sdssz'),
                database=database
                )
        # otherwise, we have no idea what this object is
        else:
            return (False, None, 'objectinfo', 'no_such_object')


    # make sure this object actually exists in the database
    if obj_twomass_info['nresults'] > 0:

        hatid = obj_twomass_info['hat_id'][0]

        # prepare the output dictionary
        resultdict = {
            'hatid':obj_twomass_info['hat_id'][0],
            'twomassid':obj_twomass_info['twomass_id'][0],
            'jmag':obj_twomass_info['jmag'][0],
            'jmag_err':obj_twomass_info['jmag_sig'][0],
            'hmag':obj_twomass_info['hmag'][0],
            'hmag_err':obj_twomass_info['hmag_sig'][0],
            'kmag':obj_twomass_info['kmag'][0],
            'kmag_err':obj_twomass_info['kmag_sig'][0],
            'qual_flag':obj_twomass_info['qual_flag'][0],
            'bmag':obj_twomass_info['bmag'][0],
            'bmag_err':None,
            'vmag':obj_twomass_info['vmag'][0],
            'vmag_err':None,
            'rmag':obj_twomass_info['rmag'][0],
            'rmag_err':None,
            'imag':obj_twomass_info['imag'][0],
            'imag_err':None,
            'sdssu':obj_twomass_info['sdssu'][0],
            'sdssu_err':None,
            'sdssg':obj_twomass_info['sdssg'][0],
            'sdssg_err':None,
            'sdssr':obj_twomass_info['sdssr'][0],
            'sdssr_err':None,
            'sdssi':obj_twomass_info['sdssi'][0],
            'sdssi_err':None,
            'sdssz':obj_twomass_info['sdssz'][0],
            'sdssz_err':None,
            'ra':obj_twomass_info['ra'][0],
            'ra_err':None,
            'decl':obj_twomass_info['decl'][0],
            'decl_err':None,
            'ra_pm':None,
            'ra_pm_err':None,
            'decl_pm':None,
            'decl_pm_err':None,
            'hat_ndet':0,
            'hat_stations':None,
            'hat_lcfpath':None,
            'hat_lastupdated':None,
            'candidate_id':None,
            'final_id':None,
            'object_type':None,
            'object_period':None,
            'object_epoch':None,
            'object_accessgroups':None,
            'hatlc_accessgroups':None
            }

        # then, get the extra info from the twomass_ucac4 table
        hat_field, hat_field_objid = [int(x) for x in hatid.split('-')[1:]]

        query = ("select a.twomass_id, "
                 "a.ucac4_ra, a.ucac4_ra_err, a.ucac4_decl, a.ucac4_decl_err, "
                 "a.ra_propermotion, a.ra_propermotion_err, "
                 "a.decl_propermotion, a.decl_propermotion_err, "
                 "a.bmag, a.bmag_err, a.vmag, a.vmag_err, "
                 "a.sdssg, a.sdssg_err, a.sdssr, a.sdssr_err, a.sdssi, a.sdssi_err "
                 "from twomass_ucac4 a "
                 "where a.hat_field = %s and a.hat_field_objid = %s")
        query_params = (hat_field, hat_field_objid)
        cur.execute(query, query_params)

        # add the query to the resultdict
        resultdict['query'] = [cur.query]

        row = cur.fetchone()

        # if this object exists in the UCAC4 table, then grab its info and
        # update the resultdict
        if row and len(row) > 0:

            resultdict['ra'] = row[1] or resultdict['ra']
            resultdict['ra_err'] = row[2] or resultdict['ra_err']
            resultdict['decl'] = row[3] or resultdict['decl']
            resultdict['decl_err'] = row[4] or resultdict['decl_err']

            resultdict['ra_pm'] = row[5] or resultdict['ra_pm']
            resultdict['ra_pm_err'] = row[6] or resultdict['ra_pm_err']
            resultdict['decl_pm'] = row[7] or resultdict['decl_pm']
            resultdict['decl_pm_err'] = row[8] or resultdict['decl_pm_err']

            resultdict['bmag'] = row[9] or resultdict['bmag']
            resultdict['bmag_err'] = row[10] or resultdict['bmag_err']

            resultdict['vmag'] = row[11] or resultdict['vmag']
            resultdict['vmag_err'] = row[12] or resultdict['vmag_err']

            resultdict['sdssg'] = row[13] or resultdict['sdssg']
            resultdict['sdssg_err'] = row[14] or resultdict['sdssg_err']

            resultdict['sdssr'] = row[15] or resultdict['sdssr']
            resultdict['sdssr_err'] = row[16] or resultdict['sdssr_err']

            resultdict['sdssi'] = row[17] or resultdict['sdssi']
            resultdict['sdssi_err'] = row[18] or resultdict['sdssi_err']

        # add this next query to the resultdict and serialize to string
        resultdict['query'].append(cur.query)

        # next, check if this object exists in the interesting_objects table
        query = ("select candidate_id, final_id, object_type, "
                 "object_period, object_epoch, "
                 "access_groups from interesting_objects "
                 "where hat_field = %s and hat_field_objid = %s")
        query_params = (hat_field, hat_field_objid)
        cur.execute(query, query_params)
        row = cur.fetchone()

        # add this next query to the resultdict and serialize to string
        resultdict['query'].append(cur.query)

        # if this object exists in the interesting_objects table, add in its
        # information
        if row and len(row) > 0:

            (resultdict['candidate_id'],
             resultdict['final_id'],
             resultdict['object_type'],
             resultdict['object_period'],
             resultdict['object_epoch'],
             resultdict['object_accessgroups']) = row


        # next, check if the object exists in the hat_lightcurves table
        query = ("select ndet, hat_stations, "
                 "full_lc_fpath, last_updated, access_groups from "
                 "hat_lightcurves where hat_id = %s")
        query_params = (hatid,)
        cur.execute(query, query_params)
        row = cur.fetchone()

        # add this next query to the resultdict and serialize to string
        resultdict['query'].append(cur.query)

        # now, add in the LC information a valid entry for this lightcurve means
        # that the object has been checked to see if it has a lightcurve. it
        # might have zero detections, though
        if row and len(row) > 0:

            resultdict['hat_lcchecked'] = True
            resultdict['hat_ndet'] = row[0] or resultdict['hat_ndet']
            resultdict['hat_stations'] = row[1] or resultdict['hat_stations']
            resultdict['hat_lcfpath'] = row[2] or resultdict['hat_lcfpath']

            # convert LC datetime to UTC and serialize to JSON
            lc_datetime = row[3] - row[3].utcoffset()

            resultdict['hat_lastupdated'] = (
                datetime_to_jsondate(lc_datetime) or
                resultdict['hat_lastupdated']
                )
            resultdict['hatlc_accessgroups'] = row[4]

        # otherwise, we know that we haven't yet attempted to collect an LC for
        # this object.
        else:

            resultdict['hat_lcchecked'] = False


        # concatenate all the queries
        resultdict['query'] = ' / '.join(resultdict['query'])

        # now that we have everything, check what we can return to the user
        hatlc_authinfo = auth_hatlcs([hatid],
                                     username=accessuser,
                                     apikey=accessapikey,
                                     existing_lcs_only=False,
                                     database=database)

        resultdict['user_accessgroup'] = hatlc_authinfo[2]

        # check and restrict access to the whole object
        if ((identifier_type in
             ('hncand','hscand','hnplanet','hsplanet')) and
            (resultdict['user_accessgroup'] not in
             resultdict['object_accessgroups'])):

            return (False, None, 'objectinfo', 'object_not_authorized')

        # check and restrict access to the LCs
        else:

            if hatlc_authinfo[0] is not False:

                if hatid in hatlc_authinfo[1]:
                    returntuple = (True, resultdict,
                                   'objectinfo', 'objectinfo_ok')

                else:
                    resultdict['hat_ndet'] = 0
                    resultdict['hat_stations'] = None
                    resultdict['hat_lcfpath'] = None
                    resultdict['hat_lastupdated'] = None
                    returntuple = (True, resultdict,
                                   'objectinfo', 'objectinfo_ok')

            else:

                if hatlc_authinfo[-1] == 'no_hatlcs_found':
                    returntuple = (True, resultdict, 'objectinfo',
                                   'no_hatlcs_exist')
                elif hatlc_authinfo[-1] == 'no_hatlcs_authorized':
                    returntuple = (True, resultdict, 'objectinfo',
                                   'hatlcs_exist_but_not_authorized')
                else:
                    returntuple = (False, None,
                                   'objectinfo', hatlc_authinfo[2])


    # if the object does not exist, tell the user
    else:
        returntuple = (False, None, 'objectinfo', 'no_such_object')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return returntuple


########################################################
### FUNCTION TO RETURN INFO FOR AN OBJECT COLLECTION ###
########################################################

def collection_info(collection,
                    selector=None,
                    accessapikey=None,
                    accessuser=None,
                    database=None):
    '''This function returns all the information needed for a collection of
    objects.

    colltype is the type of collection, one of the following strings:

    'hnplanets' -> returns all HATNet planets
    'hsplanets' -> returns all HATNet planets
    'hncandidates' -> returns all HATNet planets
    'hscandidates' -> returns all HATNet planets
    'field' -> returns all objects with authorized LCs in field
    'region' -> returns all objects with authorized LCs in coordinate
    'hatid' -> returns all objects with authorized LCs with specified hatids

    selector is used to further refine the search and is usually a string to
    match against. if colltype is 'field' or 'region' or 'hatid', then the
    selector is required.

    Returns a 4-elem tuple compatible with serviceutils functions.

    task_result is a dictionary with the following keys in most cases:

    hatid, ra, dec, jmag, hmag, kmag, ndet
    object_type, candidate_id, final_id, access_groups

    '''

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True


    # get this user's info and query restrictions
    if accessuser and not accessapikey:
        userinfo = check_user(accessuser, getinfo=True, database=database)

    elif not accessuser and accessapikey:
        apikeyinfo = get_apikey_username(accessapikey, database=database)
        if apikeyinfo is not False:
            userinfo = check_user(apikeyinfo[0], getinfo=True, database=database)
        else:
            return (False, None, 'collectioninfo', 'request_unknown_user')

    elif not accessuser and not accessapikey:
        userinfo = check_user('anonuser', getinfo=True, database=database)

    elif accessuser and accessapikey:
        return (False, None, 'collectioninfo', 'request_unknown_user')

    # return failure if no user is selected
    if userinfo is False:
        return (False, None, 'collectioninfo', 'request_unknown_user')


    # some sanity-checking for input params
    if not collection in COLLECTIONS:
        return (False, None, 'collectioninfo', 'no_such_collection')

    if COLLECTIONS[collection]['selector_required'] and not selector:
        return (False, None, 'collectioninfo', 'selector_required')

    # check that the selector matches our regex
    if selector:

        if collection == 'hatids':
            selectortouse = re.findall(HATIDLIST_REGEX, selector)

        else:
            selectortouse = re.match(COLLECTIONS[collection]['selector_regex'],
                                     selector).groups()


        if not selectortouse:
            return (False, None, 'collectioninfo', 'selector_invalid')
        else:
            # use only the first group when query isn't ra/dec
            if 'planets' in collection or 'candidates' in collection:
                selectortouse = selectortouse[0]

            # if the collection is a HAT field, turn it into an int
            elif collection == 'field':
                selectortouse = int(selectortouse[0])

            # if the collection is a list of HATIDs
            elif collection == 'hatids':
                selectortouse = tuple(selectortouse)

            # if the collection is a region
            elif collection == 'region':
                selectortouse = [float(selectortouse[0]),
                                 float(selectortouse[1])]

    else:
        selectortouse = COLLECTIONS[collection]['selector_default']


    # if the collection isn't a region, then we only need one selector and the
    # user's access_group
    if 'planets' in collection or 'candidates' in collection:

        query = COLLECTIONS[collection]['search_query']
        query_params = ('%{selector}%'.format(selector=selectortouse),
                        '%{accessgroup}%'.format(accessgroup=userinfo[3]))
        result_template = COLLECTIONS[collection]['query_results']

    elif collection == 'field':

        query = COLLECTIONS[collection]['search_query']
        query_params = (selectortouse,
                        '%{accessgroup}%'.format(accessgroup=userinfo[3]),
                        userinfo[7])
        result_template = COLLECTIONS[collection]['query_results']

    elif collection == 'hatids':

        query = COLLECTIONS[collection]['search_query']
        query_params = (selectortouse,
                        '%{accessgroup}%'.format(accessgroup=userinfo[3]),
                        userinfo[7])
        result_template = COLLECTIONS[collection]['query_results']


    # elif collection == 'region':

    #     # FIXME: deal with this later (need ra, decl in the right order and
    #     # repeated a couple of times like for a spatial query)
    #     query = COLLECTIONS[collection]['search_query']
    #     query_params = (selectortouse[0], userinfo[3])
    #     result_template = COLLECTIONS[collection]['query_results']

    #     print(query, query_params)

    # deal with unimplemented dataset handling
    else:

        # close the database at the end if we have to
        database.close_cursor(handle)
        if closedb:
            database.close_connection()

        return (False,
                None,
                'collectioninfo',
                'collection_handling_not_implemented')



    # execute the query and get the results
    try:
        cur.execute(query, query_params)
        results = cur.fetchall()

        # if we have results, then process them according to the results
        # template
        if results and len(results) > 0:

            resultdict = {'query': cur.query,
                          'dataset': results,
                          'title': COLLECTIONS[collection]['dataset_name'],
                          'columns': result_template,
                          'accessgroup': userinfo[3]}
            returntuple = (True,
                           resultdict,
                           'collectioninfo',
                           'collectioninfo_ok')

        # if no objects are returned, tell the user
        else:
            returntuple = (False,
                           None,
                           'collectioninfo',
                           'no_collection_objects')


    except Exception as e:

        LOGGER.exception('could not get results from the DB')
        raise


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return returntuple




#####################################
### FUNCTION TO RETURN FIELD INFO ###
#####################################

def field_info(hat_field,
               accessapikey=None,
               accessuser=None,
               database=None):
    '''
    This function returns all the information needed for the field info pages
    (/lightcurves/field/<field>).

    Returns a 4-elem tuple compatible with serviceutils functions.

    task_result is a dictionary with the following keys:

    fieldnum
    nobjects
    nlcs
    racenter
    deccenter

    frontend JS will parse these values from the HTML returned and make a
    request for a field-center plot to /lightcurves/stamps.

    '''



##############################
### 2MASS SEARCH FUNCTIONS ###
##############################

def twomass_info(identifiers,
                 identifier_type='hatid',
                 cols=('hat_id,ra,decl,'
                       'jmag,jmag_sig,hmag,hmag_sig,kmag,kmag_sig,'
                       'qual_flag,twomass_id,'
                       'bmag,vmag,rmag,imag,'
                       'sdssu,sdssg,sdssr,sdssi,sdssz'),
                 pprint='dict',
                 database=None,
                 verbose=False):
    '''
    Returns info for the object with identifiers from the twomass table. If
    identifiers is a list, then will return info on all of the objects in the
    list.

    identifiers is either all HAT-IDs or all 2MASS IDs. Use
    identifier_type='hatid' or identifier_type='twomass' to tell the function
    what they are.

    Extra cols not in the twomass table that can be returned:
    hat_id, bmag, vmag, rmag, imag, sdssu, sdssg, sdssr, sdssi, sdssz

    '''
    # parse and figure out the columns to return for the query
    query_cols = cols.split(',')
    query_cols = [x.strip() for x in query_cols]
    query_col_validate = [(x in TWOMASS_COLS) for x in query_cols]

    # make sure that all of the columns asked for are legit
    if not all(query_col_validate):

        failed_col_inds = [ind for (ind, x) in
                           enumerate(query_col_validate) if x is False]

        failed_cols = [query_cols[ind] for ind in failed_col_inds]
        exception_msg = ('unknown twomass table columns %s' %
                         ', '.join(failed_cols))
        raise ValueError(exception_msg)

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    # take care of the special case of asking for the HAT ID directly
    if 'hat_id' in query_cols:

        hatid_index = query_cols.index('hat_id')
        query_cols[hatid_index] = (
            "('HAT-' ||"
            " to_char(hat_field,'FM000') || '-' ||"
            " to_char(hat_field_objid,'FM0000000')) as hat_id"
            )
    else:
        raise ValueError('need at least one HAT ID to look up')


    # take care of the special case of asking for J-H color directly
    if 'jh_color' in query_cols:

        hatid_index = query_cols.index('jh_color')
        query_cols[hatid_index] = (
            "(jmag - hmag)"
            )

    # take care of the special case of asking for J-K color directly
    if 'jk_color' in query_cols:

        hatid_index = query_cols.index('jk_color')
        query_cols[hatid_index] = (
            "(jmag - kmag)"
            )

    # take care of the special case of asking for H-K color directly
    if 'hk_color' in query_cols:

        hatid_index = query_cols.index('hk_color')
        query_cols[hatid_index] = (
            "(hmag - kmag)"
            )

    # this keeps track of additional mag quantities that need to be calculated
    additional_mags = []

    query_col_set = set(query_cols)
    additional_mag_set = set(('bmag','vmag','rmag','imag',
                              'sdssu','sdssg','sdssr','sdssi','sdssz'))
    query_additional_cols = list(query_col_set.intersection(additional_mag_set))

    # take care of the special case of asking for additional mags directly
    if len(query_additional_cols) > 0:

        for addcol in query_additional_cols:

            col_index  = query_cols.index(addcol)

            if 'kmag' not in query_cols:
                query_cols.insert(col_index,'kmag')
            if 'hmag' not in query_cols:
                query_cols.insert(col_index,'hmag')
            if 'jmag' not in query_cols:
                query_cols.insert(col_index,'jmag')

            query_cols.remove(addcol)
            additional_mags.append(addcol)

    # form the query for the search
    query_col_str = ', '.join(query_cols)

    query_select = 'select %s from twomass where ' % query_col_str

    # convert a single identifier to a list to keep things simple
    if not isinstance(identifiers, list):
        identifiers = [identifiers]


    if identifier_type == 'hatid':
        identifiers = [x.split('-')[1:] for x in identifiers]
        hatids = [[int(x[0]),int(x[1])] for x in identifiers]

        query_condition_list = [
            '(hat_field = %s and hat_field_objid = %s)' % tuple(x)
            for x in hatids
            ]

        query_condition = ' or '.join(query_condition_list)

        query_constraints = []
        for x in hatids:
            query_constraints.extend([x[0],x[1]])

    elif identifier_type == 'twomass':

        twomassids = [x.strip('J').strip() for x in identifiers]

        query_condition = 'twomass_id in %s'
        query_constraints = (tuple(twomassids),)

    # build up and execute the full query
    full_query = (query_select + query_condition)

    if verbose:
        print('searching...')

    qstart = time.time()
    cur.execute(full_query, query_constraints)
    qend = time.time()
    printable_query = cur.query

    if verbose:
        print('query done. time for query: %.3f seconds' % (qend - qstart))
        print('query executed was: %s' % printable_query)

    rows = cur.fetchall()

    if verbose:
        print('fetch complete. rows fetched: %s' % len(rows))


    # now deal with additional magnitude columns
    if len(additional_mags) > 0 and len(rows) > 0:

        addmagcols = add_additional_mag_cols(query_cols, rows, additional_mags)

        resultcols = zip(*rows)
        resultdict = dict((x,y) for (x,y) in zip(query_cols, resultcols))

        # insert the extra columns into the result dictionary
        for addmag, addcol in zip(additional_mags, addmagcols):
            resultdict[addmag] = addcol

        # reform the resultdict to match the 'dist_arcsec' and 'hat_id' keywords
        resultdict_keys = resultdict.keys()
        for key in resultdict_keys:
            if 'dist_arcsec' in key:
                resultdict['dist_arcsec'] = resultdict[key]
                del resultdict[key]
            if 'hat_id' in key:
                resultdict['hat_id'] = resultdict[key]
                del resultdict[key]

        # get the query results to the way they were before adding extra cols
        initial_cols = cols.split(',')

        # get rid of the jhk columns if they're not needed
        if 'jmag' not in initial_cols and 'jmag' in resultdict:
            del resultdict['jmag']
        if 'hmag' not in initial_cols and 'hmag' in resultdict:
            del resultdict['hmag']
        if 'kmag' not in initial_cols and 'kmag' in resultdict:
            del resultdict['kmag']

        # put the resultdict back into a list in the initial query order
        resultcols = []
        for col in initial_cols:
            resultcols.append(resultdict[col])

        # flip the cols back to rows
        rows = zip(*resultcols)

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    # do a pretty-print of the requested data if asked to do so
    if pprint:
        search_results = RESULTS_FUNCMAP[pprint](cols,
                                                 rows,
                                                 'hat_id',
                                                 'asc',
                                                 printable_query)
        return search_results

    # otherwise, just return the rows
    else:
        return rows



def twomass_search(
    ra,
    dec,
    search_rad,
    database=None,
    accessuser=None,
    accessapikey=None,
    searchtype='cone',
    cols='hat_id,ra,decl,jmag,hmag,kmag,qual_flag,twomass_id,dist_arcsec',
    orderby='dist_arcsec',
    sortorder='asc',
    filters=None,
    pprint='dict',
    verbose=False
    ):
    '''
    Does a search in the twomass table around the given coordinates ra, dec,
    using a radius size of search_rad arcseconds. If searchtype = 'box', then
    search_rad is understood to mean half the box width.

    Extra cols not in the twomass table that can be returned:
    hat_id, bmag, vmag, rmag, imag, sdssu, sdssg, sdssr, sdssi, sdssz,

    database = an existing instance of lcdb.LCDB or None

    cols = the columns to return from the database, comma separated string

    searchtype = 'cone' for a cone search
                 'box' for a box search (search_rad = box_width/2)

    orderby = column to sort the result by

    sortorder = either 'asc' or 'desc'

    filters = a sequence of tuples describing the conditions on the
              columns to satisfy. these MUST follow the rules below:

              condition tuple = ('<and|or>','<column_name>','<conditions>')

              if the first item is left out and/or the condition tuple is not
              the first or only one, then all conditions are ANDed
              together. otherwise, the specified condition logical operators are
              used.

              <column_name> MUST be one of those in TWOMASS_COLS above

              <conditions> MUST follow the rules below:

              - the <conditions> string is split using spaces
              - the first token must be an <operator>
              - <operator> is checked to see if it's in COND_OPERATORS
              - <operator> must be followed by a space then <operands> each
                separated by spaces (max two operands for 'between' operator)

              filters that don't meet these rules will be discarded; a warning
              will be logged/printed for each one that ends up this way

              Example filter: (('jmag','< 13.0'),('qual_flag','= AAA'))

    pprint = 'text' to return the results as an ASCII table (space-separated)
             'csv' to return the results as an ASCII text CSV table
             'json' to return the results as a JSON string
             None for no pretty-printing

    verbose = True -> print some info about processing and results to stdout
              False -> work silently
    '''

    # get this user's info and query restrictions
    if accessuser and not accessapikey:
        userinfo = check_user(accessuser, getinfo=True, database=database)

    elif not accessuser and accessapikey:
        apikeyinfo = get_apikey_username(accessapikey, database=database)
        if apikeyinfo is not False:
            userinfo = check_user(apikeyinfo[0], getinfo=True, database=database)
        else:
            return None

    elif not accessuser and not accessapikey:
        userinfo = check_user('anonuser', getinfo=True, database=database)

    elif accessuser and accessapikey:
        return None


    if userinfo is False:
        return None
    else:
        radiuslimit = userinfo[6]*3600.0 # in arcseconds now
        rowlimit = userinfo[7]

    # parse and figure out the columns to return for the query
    query_cols = cols.split(',')
    query_cols = [x.strip() for x in query_cols]
    query_col_validate = [(x in TWOMASS_COLS) for x in query_cols]

    # make sure that all of the columns asked for are legit
    if not all(query_col_validate):

        failed_col_inds = [ind for (ind, x) in
                           enumerate(query_col_validate) if x is False]

        failed_cols = [query_cols[ind] for ind in failed_col_inds]
        exception_msg = ('unknown twomass table columns %s' %
                         ', '.join(failed_cols))
        raise ValueError(exception_msg)

    # make sure the orderby column is legit
    if not orderby in TWOMASS_COLS:
        raise ValueError('%s is not a column in the twomass table' % orderby)

    # make sure sortorder is legit
    if not sortorder in ['asc','desc']:
        raise ValueError('sortorder must be either asc or desc')

    # now proceed to actual processing

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    # take care of the special case of asking for the HAT ID directly
    if 'hat_id' in query_cols:

        hatid_index = query_cols.index('hat_id')
        query_cols[hatid_index] = (
            "('HAT-' || "
            "to_char(hat_field,'FM000') || '-' || "
            "to_char(hat_field_objid,'FM0000000')) "
            "as hat_id"
            )

    # take care of the special case of asking for distance in arcsec from query
    # center coordinates
    if 'dist_arcsec' in query_cols:

        dist_arcsec_index = query_cols.index('dist_arcsec')

        # calculate the great circle distance between the search center
        # coordinates and the ra/deg of each matched object
        gcirc_dist = (
            '(3600.0 * '
            'degrees(2.0*asin(sqrt(pow(sin((radians(decl) - '
            'radians(%s))/2.0),2) +'
            'cos(radians(%s)) * cos(radians(decl)) *'
            'pow(sin((radians(ra) - radians(%s))/2.0),2))))) as dist_arcsec'
            )
        query_cols[dist_arcsec_index] = gcirc_dist % (dec,dec,ra)

    # take care of the special case of asking for J-H color directly
    if 'jh_color' in query_cols:

        hatid_index = query_cols.index('jh_color')
        query_cols[hatid_index] = (
            "(jmag - hmag)"
            )

    # take care of the special case of asking for J-K color directly
    if 'jk_color' in query_cols:

        hatid_index = query_cols.index('jk_color')
        query_cols[hatid_index] = (
            "(jmag - kmag)"
            )

    # take care of the special case of asking for H-K color directly
    if 'hk_color' in query_cols:

        hatid_index = query_cols.index('hk_color')
        query_cols[hatid_index] = (
            "(hmag - kmag)"
            )

    # this keeps track of additional mag quantities that need to be calculated
    additional_mags = []

    query_col_set = set(query_cols)
    additional_mag_set = set(('bmag','vmag','rmag','imag',
                              'sdssu','sdssg','sdssr','sdssi','sdssz'))
    query_additional_cols = list(query_col_set.intersection(additional_mag_set))

    # take care of the special case of asking for additional mags directly
    if len(query_additional_cols) > 0:

        for addcol in query_additional_cols:

            col_index  = query_cols.index(addcol)

            if 'kmag' not in query_cols:
                query_cols.insert(col_index,'kmag')
            if 'hmag' not in query_cols:
                query_cols.insert(col_index,'hmag')
            if 'jmag' not in query_cols:
                query_cols.insert(col_index,'jmag')

            query_cols.remove(addcol)
            additional_mags.append(addcol)


    # form the query for the box search
    query_col_str = ', '.join(query_cols)

    query_select = 'select %s from twomass where ' % query_col_str
    query_sortcond = 'order by %s %s' % (orderby, sortorder)


    # deal with the search_rad limits
    if radiuslimit != -1 and radiuslimit < search_rad:
        search_radius = radiuslimit
    elif ((radiuslimit != -1 and radiuslimit > search_rad)
          or (radiuslimit == -1)):
        search_radius = search_rad

    # if the search type is a box, then just use the simple errbox method
    if searchtype == 'box':

        # figure out the box size to use
        min_ra = ra - search_radius/3600.0
        max_ra = ra + search_radius/3600.0

        min_dec = dec - search_radius/3600.0
        max_dec = dec + search_radius/3600.0

        # these are the initial query data condition strings and constraints,
        # will be added to if there are filters we need to use
        if filters:
            query_datacond = [
                '(box(point(%s,%s),point(%s,%s)) && radec_errbox) and'
                ]
        else:
            query_datacond = [
                '(box(point(%s,%s),point(%s,%s)) && radec_errbox)'
                ]

        query_constraints = [min_ra,min_dec,max_ra,max_dec]

    # if the search type is a cone, then we need to refine the initial errbox
    # selection with constraints on a circle contained inside the errbox
    elif searchtype == 'cone':

        # figure out the box size to use
        min_ra = ra - search_radius/3600.0
        max_ra = ra + search_radius/3600.0

        min_dec = dec - search_radius/3600.0
        max_dec = dec + search_radius/3600.0

        # these are the initial query data condition strings and constraints,
        # will be added to if there are filters we need to use
        query_datacond = ['(box(point(%s,%s),point(%s,%s)) && radec_errbox) and']
        query_constraints = [min_ra,min_dec,max_ra,max_dec]

        # constrain to a circle with radius equal to search_radius
        if filters:

            # from astrogrid (wrong formula?)
            # conesearch_cond = (
            #     '(degrees(2*asin(sqrt(pow(sin(radians(decl - %s)/2),2)) + '
            #     'cos(radians(%s)) * cos(radians(decl)) * '
            #     'pow(sin(radians(ra - %s)/2),2))) < %s) and'
            #     )
            conesearch_cond = (
                '(degrees(2.0*asin(sqrt(pow(sin((radians(decl) - '
                'radians(%s))/2.0),2) +'
                'cos(radians(%s)) * cos(radians(decl)) *'
                'pow(sin((radians(ra) - radians(%s))/2.0),2)))) < %s) and'
                )

        else:

            # from astrogrid (wrong formula?)
            # conesearch_cond = (
            #     '(degrees(2*asin(sqrt(pow(sin(radians(decl - %s)/2),2)) + '
            #     'cos(radians(%s)) * cos(radians(decl)) * '
            #     'pow(sin(radians(ra - %s)/2),2))) < %s)'
            #     )
            conesearch_cond = (
                '(degrees(2.0*asin(sqrt(pow(sin((radians(decl) - '
                'radians(%s))/2.0),2) +'
                'cos(radians(%s)) * cos(radians(decl)) *'
                'pow(sin((radians(ra) - radians(%s))/2.0),2)))) < %s)'
                )

        conesearch_constraints = [dec,dec,ra,search_radius/3600.0]

        query_datacond.append(conesearch_cond)
        query_constraints.extend(conesearch_constraints)


    # parse and deal with the filters if they exist
    if filters:

        for filter_ind, datafilter in enumerate(filters):

            # if the this filter has 3 elems, it needs to be ANDed/ORed with
            # other filters
            if len(datafilter) == 3:

                # the first filter can't have three elems
                if filter_ind == 0:

                    if LOGGER:
                        LOGGER.warning('filter %s is bad, ignoring...' %
                                       datafilter)
                    else:
                        print('searchutils: filter %s is bad, ignoring...' %
                              datafilter)

                    continue

                # deal with three-elem filters as usual
                else:

                    filter_op, filter_col, filter_conds = datafilter

                    if ((filter_op not in ['and','or']) or
                        (filter_col not in TWOMASS_COLS)):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           datafilter)
                        else:
                            print('searchutils: filter %s is bad, ignoring...' %
                                  datafilter)

                        continue

                    else:

                        # parse the filter conditions
                        filter_conds = filter_conds.split()
                        cond_op = filter_conds[0]

                        if cond_op not in COND_OPERATORS:

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, '
                                      'ignoring...' %
                                      datafilter)

                            continue

                        elif ((cond_op == 'between') and
                              (len(filter_conds[1:]) != 2)):

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, '
                                      'ignoring...' %
                                      datafilter)

                            continue

                        elif ((cond_op != 'between') and
                              (len(filter_conds[1:]) != 1)):

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, '
                                      'ignoring...' %
                                      datafilter)

                            continue

                        else:

                            cond_string = COND_OPERATORS[cond_op].format(
                                col=filter_col
                                )
                            cond_string = '{filter_op} {cond_string}'.format(
                                filter_op=filter_op,
                                cond_string=cond_string
                                )
                            query_datacond.append(cond_string)
                            query_constraints.extend(filter_conds[1:])

            elif len(datafilter) == 2:

                # process the first filter
                if filter_ind == 0:

                    filter_col, filter_conds = datafilter

                    if (filter_col not in TWOMASS_COLS):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           datafilter)
                        else:
                            print('searchutils: filter %s is bad, ignoring...' %
                                  datafilter)

                        continue

                    else:

                        # parse the filter conditions
                        filter_conds = filter_conds.split()
                        cond_op = filter_conds[0]

                        if cond_op not in COND_OPERATORS:

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, ignoring...' %
                                      datafilter)

                            continue

                        elif ((cond_op == 'between') and
                              (len(filter_conds[1:]) != 2)):

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, ignoring...' %
                                      datafilter)

                            continue

                        elif ((cond_op != 'between') and
                              (len(filter_conds[1:]) != 1)):

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, ignoring...' %
                                      datafilter)

                            continue

                        else:

                            cond_string = COND_OPERATORS[cond_op].format(
                                col=filter_col
                                )
                            query_datacond.append(cond_string)
                            query_constraints.extend(filter_conds[1:])

                # if this filter isn't the first one and it has only two elems,
                # then it needs to be ANDed with any previous filters
                else:

                    filter_col, filter_conds = datafilter
                    filter_op = 'and'

                    if ((filter_op not in ['and','or']) or
                        (filter_col not in TWOMASS_COLS)):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           datafilter)
                        else:
                            print('searchutils: filter %s is bad, ignoring...' %
                                  datafilter)

                        continue

                    else:

                        # parse the filter conditions
                        filter_conds = filter_conds.split()
                        cond_op = filter_conds[0]

                        if cond_op not in COND_OPERATORS:

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, '
                                      'ignoring...' %
                                      datafilter)

                            continue

                        elif ((cond_op == 'between') and
                              (len(filter_conds[1:]) != 2)):

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, '
                                      'ignoring...' %
                                      datafilter)

                            continue

                        elif ((cond_op != 'between') and
                              (len(filter_conds[1:]) != 1)):

                            if LOGGER:
                                LOGGER.warning('filter %s is bad, ignoring...' %
                                               datafilter)
                            else:
                                print('searchutils: filter %s is bad, '
                                      'ignoring...' %
                                      datafilter)

                            continue

                        else:

                            cond_string = COND_OPERATORS[cond_op].format(
                                col=filter_col
                                )
                            cond_string = '{filter_op} {cond_string}'.format(
                                filter_op=filter_op,
                                cond_string=cond_string
                                )
                            query_datacond.append(cond_string)
                            query_constraints.extend(filter_conds[1:])

            else:

                if LOGGER:
                    LOGGER.warning('filter %s is bad, ignoring...' %
                                   datafilter)
                else:
                    print('searchutils: filter %s is bad, ignoring...' %
                          datafilter)

                continue

    # deal with the query row limit
    if rowlimit == -1:
        query_rowlimit = ''
    else:
        query_rowlimit = ' fetch first %i rows only' % rowlimit

    # build up and execute the full query
    full_query = (query_select + ' '.join(query_datacond) +
                  ' ' + query_sortcond + query_rowlimit)
    printable_query = full_query % tuple(query_constraints)

    if verbose:
        print('query: %s' % printable_query)
        print('searching...')

    qstart = time.time()
    cur.execute(full_query,query_constraints)
    qend = time.time()

    if verbose:
        print('query done. time for query: %.3f seconds' % (qend - qstart))

    rows = cur.fetchall()

    if verbose:
        print('fetch complete. rows fetched: %s' % len(rows))

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    # now deal with additional magnitude columns
    if len(additional_mags) > 0 and len(rows) > 0:

        addmagcols = add_additional_mag_cols(query_cols, rows, additional_mags)

        resultcols = zip(*rows)
        resultdict = dict((x,y) for (x,y) in zip(query_cols, resultcols))

        # insert the extra columns into the result dictionary
        for addmag, addcol in zip(additional_mags, addmagcols):
            resultdict[addmag] = addcol

        # reform the resultdict to match the 'dist_arcsec' and 'hat_id' keywords
        resultdict_keys = resultdict.keys()
        for key in resultdict_keys:
            if 'dist_arcsec' in key:
                resultdict['dist_arcsec'] = resultdict[key]
                del resultdict[key]
            if 'hat_id' in key:
                resultdict['hat_id'] = resultdict[key]
                del resultdict[key]

        # get the query results to the way they were before adding extra cols
        initial_cols = cols.split(',')

        # get rid of the jhk columns if they're not needed
        if 'jmag' not in initial_cols and 'jmag' in resultdict:
            del resultdict['jmag']
        if 'hmag' not in initial_cols and 'hmag' in resultdict:
            del resultdict['hmag']
        if 'kmag' not in initial_cols and 'kmag' in resultdict:
            del resultdict['kmag']

        # put the resultdict back into a list in the initial query order
        resultcols = []
        for col in initial_cols:
            resultcols.append(resultdict[col])

        # flip the cols back to rows
        rows = zip(*resultcols)


    # do a pretty-print of the requested data if asked to do so
    if pprint:
        search_results = RESULTS_FUNCMAP[pprint](cols,
                                                 rows,
                                                 orderby,
                                                 sortorder,
                                                 printable_query)
        return search_results

    # otherwise, just return the rows
    else:
        return rows


def twomass_objects_near(ra, dec, search_radius,
                         accessuser=None,
                         accessapikey=None,
                         database=None,
                         twomassids=False,
                         orderby='dist_arcsec',
                         sortorder='asc',
                         filters=None,
                         pprint='dict',
                         verbose=False,
                         searchtype='cone',
                         returndist=False):
    '''

    This is a simplified function to only return HAT (and optionally, 2MASS IDs)
    within search_radius arcseconds of position (ra, dec). Uses twomass_search
    above with a search type of 'cone'.

    database = an existing instance of lcdb.LCDB or None

    twomassids = set True to return 2MASS IDs as well as HAT-IDs.

    orderby = column to sort the result by

    sortorder = either 'asc' or 'desc'

    filters = a sequence of tuples describing the conditions on the
              columns to satisfy. these MUST follow the rules below:

              condition tuple = ('<and|or>','<column_name>','<conditions>')

              if the first item is left out and/or the condition tuple is not
              the first or only one, then all conditions are ANDed
              together. otherwise, the specified condition logical operators are
              used.

              <column_name> MUST be one of those in TWOMASS_COLS above

              <conditions> MUST follow the rules below:

              - the <conditions> string is split using spaces
              - the first token must be an <operator>
              - <operator> is checked to see if it's in COND_OPERATORS
              - <operator> must be followed by a space then <operands> each
                separated by spaces (max two operands for 'between' operator)

              Example filter: (('jmag','< 13.0'),('qual_flag','= AAA'))

    pprint = 'text' to return the results as an ASCII table
             'json' to return the results as a JSON string
             'csv' to return the results as a CSV table
             None for no pretty-printing

    '''

    if twomassids:
        search_cols = 'hat_id,twomass_id,dist_arcsec'
    else:
        search_cols = 'hat_id,dist_arcsec'

    search_results = twomass_search(
        ra, dec,
        search_radius,
        database=database,
        cols=search_cols,
        orderby=orderby,
        sortorder=sortorder,
        filters=filters,
        pprint=pprint,
        verbose=verbose,
        searchtype=searchtype,
        accessuser=accessuser,
        accessapikey=accessapikey,
        )

    # transpose the results to return columns instead of rows (makes more sense
    # this way)
    if search_results and not pprint:

        # return only the HAT-IDs by default
        if returndist or twomassids:
            search_results = zip(*search_results)
        # otherwise, return all columns
        else:
            search_results = zip(*search_results)[0]

    return search_results



def twomass_search_speedtest(ntrials, mradius_range,
                             accessuser=None):
    '''
    This just runs a bunch of searches at random coordinates with a number of
    matching radii to see how fast the twomass_search function works.

    '''

    trialtimes = dict()

    for radius in mradius_range:

        print('running %s search trials for radius = %s arcsec...' % (ntrials,
                                                                      radius))

        ra_range = numpy.random.rand(ntrials) * 360.0
        dec_range = (numpy.random.randint(-1, 1, ntrials) *
                     numpy.random.rand(ntrials) * 90.0)

        trialtimes[radius] = []

        for ra, dec in zip(ra_range, dec_range):

            starttime = time.time()
            search_results = twomass_search(ra, dec, radius, filters=None,
                                            accessuser=accessuser)
            trialtimes[radius].append(time.time() - starttime)

    return trialtimes


###########################
## HAT OBJECT DB HELPERS ##
###########################

HATCLASSES = {
    'A':'A star',
    'Astar':'A star',
    'BD':'brown dwarf',
    'BP':'bad photometry,false alarm',
    'DTO':'blend with nearby EB',
    'DW':'dwarf star',
    'EB':'eclipsing binary',
    'FM':'F+M binary',
    'G':'giant star',
    'GI':'giant star',
    'GM':'G+M binary',
    'IPC':'inflated planet candidate',
    'KELT':'KELT survey candidate',
    'Mdwarf':'M dwarf',
    'MEB':'M dwarf EB',
    'Neptune':'Neptune-sized candidate',
    'P':'confirmed planet',
    'PC':'strong planet candidate',
    'QES':'QES candidate',
    'QUAD':'quadruple system',
    'ROT':'rapid rotator',
    'SB':'spectroscopic binary',
    'TR':'transit candidate',
    'TRI':'triple system',
    }


def parse_hat_classes(hatclass):
    '''
    This parses the HTRclass or HATSclass column entries in the HATRED.HTR or
    HSCAND.HATS mysql tables. Returns a string suitable for direct insertion
    into the LC server DB's interesting_objects.object_type column.

    '''

    classlist = hatclass.split(',')
    translatedtags = []

    for hatclass in classlist:
        for tag in HATCLASSES:
            if tag in hatclass and '?' in tag:
                translatedtags.append(HATCLASSES[tag] + '?')
            elif tag in hatclass:
                translatedtags.append(HATCLASSES[tag])

    translatedtags = list(set(translatedtags))

    # if HAT-P, or HATS- is in hatclass, then add a hatplanet tag to
    # translatedtags
    if 'HAT-P' in hatclass or 'HATS-' in hatclass:
        translatedtags.append('HAT planet')

    return ','.join(translatedtags)
