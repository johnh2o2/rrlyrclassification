#!/usr/bin/env python

'''
lcutils_hatlc.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Dec 2013

Contains various useful tools for reading, writing, and filtering consolidated
HAT LCs produced by the HAT LC server.

'''

import logging
import os.path
import hashlib
import gzip
import bz2


import numpy as np
import pyfits

###################
## LOCAL IMPORTS ##
###################

#from timeutils import jd_to_bjd, jd_to_hjd # IMPORT TIMEUTILS later; you need CSPICE/PySPICE. Check if Della has those... // They do NOT!
from timeutils_temp import jd_to_bjd, jd_to_hjd


import lcutils_config as conf
import lcutils_formatters as lcform
import lcutils_postprocessing as lcpp


#############
## LOGGING ##
#############

# setup a logger
LOGGER = logging.getLogger('lcutils_hatlc')
LOGGER.addHandler(logging.NullHandler())

# default to debug mode = False
DEBUGMODE = False

def set_debug(debugbool):
    globals()['DEBUGMODE'] = debugbool
    lcform.set_debug(True)


#################################################
## FUNCTIONS TO DEAL WITH THE CONSOLIDATED LCS ##
#################################################

def read_consolidated_hatlc(hatlc,
                            forceformat=None,
                            forcecompression=None):
    '''
    This reads a consolidated HAT LC written by the functions above.

    Returns a dict.

    '''

    lcfname = os.path.basename(hatlc)

    if forceformat and forcecompression:

        if 'gz' in forcecompression:
            lcf = gzip.open(hatlc,'rb')
        elif 'bz2' in forcecompression:
            lcf = bz2.BZ2File(hatlc,'rb')

    else:

        # unzip the files first
        if '.gz' in lcfname:
            lcf = gzip.open(hatlc,'rb')
        elif '.bz2' in lcfname:
            lcf = bz2.BZ2File(hatlc, 'rb')
        else:
            lcf = open(hatlc,'rb')

    if '.fits' in lcfname or (forceformat and forceformat == 'fits'):

        hdulist = pyfits.open(lcf)
        objectinfo = hdulist[0].header
        objectlc = hdulist[1].data
        lccols = objectlc.columns.names
        hdulist.close()
        lcf.close()

        lcdict = {}

        for col in lccols:
            lcdict[col] = np.array(objectlc[col])

        lcdict['hatid'] = objectinfo['hatid']
        lcdict['twomassid'] = objectinfo['2massid']
        lcdict['ra'] = objectinfo['ra']
        lcdict['dec'] = objectinfo['dec']
        lcdict['mags'] = [objectinfo[x] for x in ('vmag','rmag','imag',
                                                  'jmag','hmag','kmag')]
        lcdict['ndet'] = objectinfo['ndet']
        lcdict['hatstations'] = objectinfo['hats']
        lcdict['filters'] = objectinfo['filters']
        lcdict['columns'] = lccols

        return lcdict

    elif '.csv' in lcfname or '.hatlc' in lcfname or (forceformat and
                                                      forceformat == 'csv'):

        lcflines = lcf.readlines()
        lcf.close()

        # now process the read-in LC
        objectdata = [x for x in lcflines if x.startswith('#')]
        objectlc = [x for x in lcflines if not x.startswith('#')]
        objectlc = [x for x in objectlc if len(x) > 1]

        if '.csv' in lcfname:
            objectlc = [x.split(',') for x in objectlc]
        else:
            objectlc = [x.split() for x in objectlc]

        # transpose split rows to get columns
        objectlc = zip(*objectlc)

        # read the header to figure out the object's info and column names
        objectdata = [x.strip('#') for x in objectdata]
        objectdata = [x.strip() for x in objectdata]
        objectdata = [x for x in objectdata if len(x) > 0]

        hatid, twomassid = objectdata[0].split(' - ')
        ra, dec = objectdata[1].split(', ')
        ra = float(ra.split(' = ')[-1].strip(' deg'))
        dec = float(dec.split(' = ')[-1].strip(' deg'))

        vmag, rmag, imag, jmag, hmag, kmag = objectdata[2].split(', ')
        vmag = float(vmag.split(' = ')[-1])
        rmag = float(rmag.split(' = ')[-1])
        imag = float(imag.split(' = ')[-1])
        jmag = float(jmag.split(' = ')[-1])
        hmag = float(hmag.split(' = ')[-1])
        kmag = float(kmag.split(' = ')[-1])

        ndet = int(objectdata[3].split(': ')[-1])
        hatstations = objectdata[4].split(': ')[-1]

        filterhead_ind = objectdata.index('Filters used:')
        columnhead_ind = objectdata.index('Columns:')

        filters = objectdata[filterhead_ind:columnhead_ind]

        columndefs = objectdata[columnhead_ind+1:]

        columns = []
        for line in columndefs:

            colnum, colname, coldesc = line.split(' - ')
            columns.append(colname)

        lcdict = {}

        # now write all the columns to the output dictionary
        for ind, col in enumerate(columns):

            # this formats everything nicely using our existing column
            # definitions
            lcdict[col] = np.array([conf.TEXTLC_OUTPUT_COLUMNS[col][3](x)
                                    for x in objectlc[ind]])

        # write the object metadata to the output dictionary
        lcdict['hatid'] = hatid
        lcdict['twomassid'] = twomassid.strip('2MASS J')
        lcdict['ra'] = ra
        lcdict['dec'] = dec
        lcdict['mags'] = [vmag, rmag, imag, jmag, hmag, kmag]
        lcdict['ndet'] = ndet
        lcdict['hatstations'] = hatstations.split(', ')
        lcdict['filters'] = filters[1:]
        lcdict['cols'] = columns

        return lcdict



def write_consolidated_hatlc(lcdict,
                             outtype='ssv',
                             compress='gz',
                             outdir=None,
                             filterhash=None,
                             normalized=False):
    '''
    This writes a LC dict from a HATLC back to a HATLC. Useful if any
    modifications such as filtering on columns, etc. have been made to the
    lcdict.

    '''

    # generate the objinfo tuple
    objinfo = (lcdict['hatid'],
               lcdict['ra'],
               lcdict['dec'],
               lcdict['mags'][0],
               lcdict['mags'][1],
               lcdict['mags'][2],
               lcdict['mags'][3],
               lcdict['mags'][4],
               lcdict['mags'][5],
               lcdict['twomassid'])


    if outtype == 'json':
        outlc = lcform.OUTPUTLC_FORMATTERS['json'](lcdict,
                                                   objinfo)


    elif outtype is not None:
        outlc = lcform.OUTPUTLC_FORMATTERS[outtype](lcdict,
                                                    objinfo,
                                                    outputdir=outdir,
                                                    compress=compress,
                                                    filterhash=filterhash,
                                                    normalized=normalized)
    else:
        outlc = lcform.OUTPUTLC_FORMATTERS['lcdict'](lcdict,
                                                     objinfo)

    return outlc



def reform_consolidated_hatlc(lcdict,
                              columns='default',
                              outtype='ssv',
                              outdir=None,
                              compress='gz',
                              normalize=False,
                              filterhash=None):
    '''
    This reforms a consolidated HATLC dictionary, returning only the columns
    requested in columns:

    columns = '<default|epdlc|tfalc|redlc>' or a comma separated list of columns
              from TEXTLC_OUTPUT_COLUMNS

    The other parameters are the same as for get_hat_lc.

    '''

    # parse the columns
    if columns in ('default', 'epdlc', 'tfalc', 'redlc', 'fullred', 'full'):

        cols = conf.HATLC_OUTPUT_COLUMNS[columns]

    else:

        cols = columns.split(',')

        if not all([col in conf.TEXTLC_OUTPUT_COLUMNS for col in cols]):
            if LOGGER:
                LOGGER.error('some of the requested columns are invalid')
            if DEBUGMODE:
                print('lcutils: some of the requested columns are invalid')

            return None

    # make a copy of the lcdict and remove the columns we don't want
    modlcdict = lcdict.copy()

    # cols we don't want
    colswantedset = set(cols)
    colspresentset = set(lcdict['cols'])
    colstoremove = colspresentset - colswantedset

    for col in colstoremove:
        del modlcdict[col]

    modlcdict['cols'] = cols

    # now deal with columns to add to the output dict

    # deal with an RSTF column
    if 'RSTF' in cols:

        rstf = ['%s-%s' % (x,y) for (x,y) in zip(lcdict['STF'],
                                                 lcdict['CFN'])]
        modlcdict['RSTF'] = rstf

    # deal with a MJD column
    if ('MJD' in cols and
        'MJD' not in lcdict['cols'] and
        'RJD' in lcdict['cols']):

        rjd = lcdict['RJD']
        rjd = np.array([float(x) for x in rjd])

        # convert to MJD
        mjd = rjd - 0.5

        modlcdict['MJD'] = mjd.tolist()

    # deal with a FJD column
    if ('FJD' in cols and
        'FJD' not in lcdict['cols'] and
        'RJD' in lcdict['cols']):

        rjd = lcdict['RJD']
        rjd = np.array([float(x) for x in rjd])

        fjd = rjd + 2400000.0

        modlcdict['FJD'] = fjd.tolist()

    # get the obslon, obslat, obsalt required for conversion to BJD or HJD, but
    # only if the HJD/BJD columns aren't in the input lcdict
    if 'BJD' in cols or 'HJD' in cols and (not ('BJD' in lcdict['cols'] or
                                                'HJD' in lcdict['cols'])):

        # grab the STF column and get the station ids
        stf = lcdict['STF']
        networks = lcdict['NET']

        # get the station locations to convert to BJD
        obslat, obslon, obsalt = ([conf.HAT_LOCATIONS[x][y][0] for
                                   (x,y) in zip(networks,stf)],
                                  [conf.HAT_LOCATIONS[x][y][1] for
                                   (x,y) in zip(networks,stf)],
                                  [conf.HAT_LOCATIONS[x][y][2] for
                                   (x,y) in zip(networks,stf)])
        objra = lcdict['ra']
        objdec = lcdict['dec']

        # deal with a requested BJD column
        if ('BJD' in cols and
            'BJD' not in lcdict['cols'] and
            'RJD' in lcdict['cols']):

            rjd = lcdict['RJD']
            rjd = np.array([float(x) for x in rjd])

            # convert to full JD
            rjd = rjd + 2400000.0

            # convert to BJD
            bjd = [
                jd_to_bjd(x,objra,objdec,y,z,w)
                for (x,y,z,w) in zip(rjd,obslat,obslon,obsalt)
                ]

            modlcdict['BJD'] = bjd

        # deal with a HJD column
        if ('HJD' in cols and
              'HJD' not in lcdict['cols'] and
              'RJD' in lcdict['cols']):

            rjd = lcdict['RJD']
            rjd = np.array([float(x) for x in rjd])

            # convert to full JD
            rjd = rjd + 2400000.0

            # convert to BJD
            hjd = [jd_to_hjd(x,objra,objdec) for x in rjd]

            modlcdict['HJD'] = hjd


    # if we're supposed to normalize magnitudes, do so
    if normalize:

        # figure out what columns we're normalizing
        if normalize is True:
            normalizecols = 'all'
        elif isinstance(normalize, str) or isinstance(normalize, unicode):
            normalizecols = normalize

        lcpp.normalize_mags_to_group_median(modlcdict,
                                            magcols=normalizecols)



    # now write this modified lcdict as the requested output type
    # note that outtype=None just returns the modified lcdict
    outlc = write_consolidated_hatlc(modlcdict,
                                     outtype=outtype,
                                     outdir=outdir,
                                     compress=compress,
                                     filterhash=filterhash,
                                     normalized=normalize)

    del modlcdict

    return outlc



def parse_hatlc_filters(filters):
    '''
    This parses a filter sequence as defined in filter_consolidated_hatlc below.

    '''

    filter_operators = {'and':'&',
                        'or':'|'}

    # this holds the elements for the np.where statement we're building
    where_eval_strlist = []

    for filter_ind, datafilter in enumerate(filters):

        # if the this filter has 3 elems, it needs to be ANDed/ORed with
        # other filters
        if len(datafilter) == 3:

            # the first filter can't have three elems
            if filter_ind == 0:

                if LOGGER:
                    LOGGER.warning('filter %s is bad, ignoring...' %
                                   list(datafilter))
                if DEBUGMODE:
                    print('lcutils.parse_hatlc_filters: filter %s is bad,'
                          ' ignoring...' %
                          list(datafilter))

                continue

            # deal with three-elem filters as usual
            else:

                filter_op, filter_col, filter_conds = datafilter

                # if the filter operator is not legit or if the filter column
                # doesn't make sense, ignore this filter
                if ((filter_op not in ['and','or']) or
                    (filter_col not in conf.TEXTLC_OUTPUT_COLUMNS)):

                    if LOGGER:
                        LOGGER.warning('filter %s is bad, ignoring...' %
                                       list(datafilter))
                    if DEBUGMODE:
                        print('lcutils.parse_hatlc_filters: '
                              'filter %s is bad, ignoring...' %
                              list(datafilter))

                    continue

                # otherwise, try to process this filter
                else:

                    # parse the filter conditions
                    filter_conds = filter_conds.split()
                    cond_op = filter_conds[0]

                    # if the operator is not legal, ignore this filter
                    if cond_op not in conf.CONDITION_OPERATORS:

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # if the operator is 'between' and there aren't two filter
                    # operands, ignore this filter
                    elif ((cond_op == 'between') and
                          (len(filter_conds[1:]) != 2)):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # if the operator is not 'between' and there's more than one
                    # operand, ignore this filter
                    elif ((cond_op != 'between') and
                          (len(filter_conds[1:]) != 1)):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # if all the checks above pass, then process the filter
                    # definition
                    else:

                        filter_cond_elements = filter_conds
                        filter_cond_operator = filter_cond_elements[0]

                        # deal with the between operator
                        if (len(filter_cond_elements[1:]) > 1 and
                            filter_cond_operator == 'between'):

                            operand1, operand2 = filter_cond_elements[1:]

                            if operand1.isalpha():
                                operand1 = "'%s'" % operand1
                            if operand2.isalpha():
                                operand2 = "'%s'" % operand2

                            # take care of the special case of FLT, FLD having
                            # string format but possibly numerical values
                            if filter_col in ('FLT', 'FLD'):
                                operand1 = "'%s'" % operand1
                                operand2 = "'%s'" % operand2

                            condition_string = conf.CONDITION_OPERATORS[
                                filter_cond_operator
                                ].format(
                                col=filter_col,
                                operand1=operand1,
                                operand2=operand2
                                )

                        # deal with the other operators
                        else:

                            operand1 = filter_cond_elements[1]

                            if operand1.isalpha():
                                operand1 = "'%s'" % operand1

                            # take care of the special case of FLT, FLD having
                            # string format but possibly numerical values
                            if filter_col in ('FLT', 'FLD'):
                                operand1 = "'%s'" % operand1

                            condition_string = conf.CONDITION_OPERATORS[
                                filter_cond_operator
                                ].format(
                                col=filter_col,
                                operand1=operand1
                                )

                        full_filter_string = '%s %s' % (
                            filter_operators[filter_op],
                            condition_string
                            )
                        where_eval_strlist.append(full_filter_string)


        # deal with two element filter definitions
        elif len(list(datafilter)) == 2:

            # process the first filter
            if filter_ind == 0:

                filter_col, filter_conds = datafilter

                # if the filter column is not legit, ignore this filter
                if (filter_col not in conf.TEXTLC_OUTPUT_COLUMNS):

                    if LOGGER:
                        LOGGER.warning('filter %s is bad, ignoring...' %
                                       list(datafilter))
                    if DEBUGMODE:
                        print('lcutils.parse_hatlc_filters:'
                              ' filter %s is bad, ignoring...' %
                              list(datafilter))

                    continue

                # otherwise, try to process this filter
                else:

                    # parse the filter conditions
                    filter_conds = filter_conds.split()
                    cond_op = filter_conds[0]

                    # if the filter condition's operator isn't legit, ignore
                    # this filter
                    if cond_op not in conf.CONDITION_OPERATORS:

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # if the filter condition's operator is between and the
                    # number of filter operands is not 2, then ignore this
                    # filter
                    elif ((cond_op == 'between') and
                          (len(filter_conds[1:]) != 2)):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # otherwise, if we have another filter, but the number of
                    # operands is not 1, ignore this filter
                    elif ((cond_op != 'between') and
                          (len(filter_conds[1:]) != 1)):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # if all of the tests above pass, then process this filter
                    else:

                        filter_cond_elements = filter_conds
                        filter_cond_operator = filter_cond_elements[0]

                        # deal with the between operator
                        if (len(filter_cond_elements[1:]) > 1 and
                            filter_cond_operator == 'between'):

                            operand1, operand2 = filter_cond_elements[1:]

                            if operand1.isalpha():
                                operand1 = "'%s'" % operand1
                            if operand2.isalpha():
                                operand2 = "'%s'" % operand2

                            # take care of the special case of FLT, FLD having
                            # string format but possibly numerical values
                            if filter_col in ('FLT', 'FLD'):
                                operand1 = "'%s'" % operand1
                                operand2 = "'%s'" % operand2

                            condition_string = conf.CONDITION_OPERATORS[
                                filter_cond_operator
                                ].format(
                                col=filter_col,
                                operand1=operand1,
                                operand2=operand2
                                )

                        # deal with the other operators
                        else:

                            operand1 = filter_cond_elements[1]
                            if operand1.isalpha():
                                operand1 = "'%s'" % operand1

                            # take care of the special case of FLT, FLD having
                            # string format but possibly numerical values
                            if filter_col in ('FLT', 'FLD'):
                                operand1 = "'%s'" % operand1

                            condition_string = conf.CONDITION_OPERATORS[
                                filter_cond_operator
                                ].format(
                                col=filter_col,
                                operand1=operand1
                                )

                        full_filter_string = condition_string
                        where_eval_strlist.append(full_filter_string)

            # if this filter isn't the first one and it has only two elems,
            # then it needs to be ANDed with any previous filters
            else:

                filter_col, filter_conds = datafilter
                filter_conds = filter_conds.split()
                filter_op = 'and'

                # if the filter operator or the filter column don't make sense,
                # then ignore this filter
                if ((filter_op not in ['and','or']) or
                    (filter_col not in conf.TEXTLC_OUTPUT_COLUMNS)):

                    if LOGGER:
                        LOGGER.warning('filter %s is bad, ignoring...' %
                                       list(datafilter))
                    if DEBUGMODE:
                        print('lcutils.parse_hatlc_filters:'
                              ' filter %s is bad, ignoring...' %
                              list(datafilter))

                    continue

                # otherwise, try to process this filter
                else:

                    # parse the filter conditions
                    cond_op = filter_conds[0]

                    # if the filter condition's operator isn't legit, then
                    # ignore this filter
                    if cond_op not in conf.CONDITION_OPERATORS:

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # if the filter condition's operator is 'between' but the
                    # number of operands is not 2, then ignore this filter
                    elif ((cond_op == 'between') and
                          (len(filter_conds[1:]) != 2)):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # if the filter condition's operator is not 'between' but
                    # the number of operands is not 1, then ignore this filter
                    elif ((cond_op != 'between') and
                          (len(filter_conds[1:]) != 1)):

                        if LOGGER:
                            LOGGER.warning('filter %s is bad, ignoring...' %
                                           list(datafilter))
                        if DEBUGMODE:
                            print('lcutils.parse_hatlc_filters:'
                                  ' filter %s is bad, ignoring...' %
                                  list(datafilter))

                        continue

                    # if all of the tests above pass, then process this filter
                    else:

                        filter_cond_elements = filter_conds
                        filter_cond_operator = filter_cond_elements[0]

                        # deal with the between operator
                        if (len(filter_cond_elements[1:]) > 1 and
                            filter_cond_operator == 'between'):

                            operand1, operand2 = filter_cond_elements[1:]

                            if operand1.isalpha():
                                operand1 = "'%s'" % operand1
                            if operand2.isalpha():
                                operand2 = "'%s'" % operand2

                            # take care of the special case of FLT, FLD having
                            # string format but possibly numerical values
                            if filter_col in ('FLT', 'FLD'):
                                operand1 = "'%s'" % operand1
                                operand2 = "'%s'" % operand2

                            condition_string = conf.CONDITION_OPERATORS[
                                filter_cond_operator
                                ].format(
                                col=filter_col,
                                operand1=operand1,
                                operand2=operand2
                                )

                        # deal with the other operators
                        else:

                            operand1 = filter_cond_elements[1]

                            if operand1.isalpha():
                                operand1 = "'%s'" % operand1

                            # take care of the special case of FLT, FLD having
                            # string format but possibly numerical values
                            if filter_col in ('FLT', 'FLD'):
                                operand1 = "'%s'" % operand1

                            condition_string = conf.CONDITION_OPERATORS[
                                filter_cond_operator
                                ].format(
                                col=filter_col,
                                operand1=operand1
                                )

                        full_filter_string = '%s %s' % (
                            filter_operators[filter_op],
                            condition_string
                            )
                        where_eval_strlist.append(full_filter_string)

        # if we can't parse this filter at all, then ignore it
        else:

            if LOGGER:
                LOGGER.warning('filter %s is bad, ignoring...' %
                               list(datafilter))
            if DEBUGMODE:
                print('lcutils.parse_hatlc_filters:'
                      ' filter %s is bad, ignoring...' %
                      list(datafilter))

            continue


    # now check for dangling filter operators caused by broken filters in the
    # chain
    where_eval_checklist = where_eval_strlist[:]
    for find, filt in enumerate(where_eval_checklist):
        if find == 0 and filt[0] == '&':
            where_eval_strlist[0] = (where_eval_strlist[0]).lstrip('& ')

    where_eval_string = (
        'np.where({array_expression})'.format(
            array_expression=' '.join(where_eval_strlist)
            )
        )

    return where_eval_string



def filter_consolidated_hatlc(hatlc,
                              filters,
                              hatlc_is_dict=False,
                              outtype='ssv',
                              compress='gz',
                              outdir=None):
    '''
    This filters a consolidated hatlc using the filter list given in filters.

    filters is a sequence of tuples describing the conditions on the columns to
    satisfy. these MUST follow the rules below.

    condition tuple = ('<and|or>','<column_name>','<conditions>')

    if the first item is left out and/or the condition tuple is not the first or
    only one, then all conditions are ANDed together. otherwise, the specified
    condition logical operators are used.

    <column_name> MUST be one of those in lcutils_config.TEXTLC_OUTPUT_COLUMNS

    <conditions> MUST be any one of the following strings:

       '= <operand>',
       '!= <operand>',
       '> <operand>',
       '< <operand>',
       '<= <operand>',
       '>= <operand>',
       'between <operand> <operand>'

    filters that don't meet these rules will be discarded; a warning will be
    logged/printed for each one that ends up this way

    '''

    # if the input is an lcdict instead of a HATLC filename, then just use it
    # directly
    if hatlc_is_dict:
        lc = hatlc
    else:
        lc = read_consolidated_hatlc(hatlc)

    if not isinstance(filters, list):
        filters = [x for x in filters]

    where_eval_string = parse_hatlc_filters(filters)

    # make sure that we have at least one filter to apply
    if where_eval_string != 'np.where()':

        if LOGGER:
            LOGGER.warning('evaling string %s' % where_eval_string)
        if DEBUGMODE:
            print('lcutils: warning: evaling string %s' % where_eval_string)

        # eval the full np.where expression and get the results
        filter_result_indices = eval(where_eval_string)

        # if filter_result_indices has a nonzero length, then use it to filter
        # the lightcurve columns
        if len(filter_result_indices[0]) > 0:

            for col in lc['cols']:
                lc[col] = lc[col][filter_result_indices]

            # amend the ndet key in the LC dict as required
            # after the filtering process
            lc['ndet'] = len(filter_result_indices[0])

            # amend hatstations key in the LC dict as required after the
            # filtering process
            if 'STF' in lc and 'NET' in lc:

                nets_and_stations = ['%s%02i' % (x,y) for
                                     (x,y) in zip(lc['NET'], lc['STF'])]
                remaining_stations = list(set(nets_and_stations))
                lc['hatstations'] = remaining_stations

            else:
                # append a warning saying that the station list is out of date
                # if the STF column is not present
                lc['hatstations'].append('warning: station list is '
                                         'from original LC since '
                                         'no STF/NET columns in original LC')

            # amend filters key in the LC dict as required after the filtering
            # process
            if 'FLT' in lc:

                remaining_filters = np.unique(lc['FLT'])
                lc['filters'] = remaining_filters.tolist()

            else:
                # append a warning saying that the filter list
                # is out of date if the FLT column is not present
                # append a warning saying that the station list is out of date
                # if the STF column is not present
                lc['filters'].append('warning: filter list is '
                                     'from original LC since '
                                     'no FLT column in original LC')


            # generate a filter hash describing what filters were applied
            lcfilters_string = [list(x) for x in filters]
            lcfilters_string = repr(lcfilters_string)
            lcfilterhash = hashlib.md5(lcfilters_string).hexdigest()

            # return a dict if the outtype is none
            if outtype is None:

                return (lc, 'filters_ok', lcfilterhash)

            # otherwise, write the actual LC to disk
            else:

                outfname = write_consolidated_hatlc(lc,
                                                    outtype=outtype,
                                                    compress=compress,
                                                    outdir=outdir,
                                                    filterhash=lcfilterhash)

                return (outfname, 'filters_ok', lcfilterhash)

        else:

            if LOGGER:
                LOGGER.warning('lcutils.filter_consolidated_hatlc: '
                               'applied filters return nothing')
            if DEBUGMODE:
                print('lcutils.filter_consolidated_hatlc: '
                      'applied filters return nothing')

            # generate a filter hash describing what filters were applied
            lcfilters_string = [list(x) for x in filters]
            lcfilters_string = repr(lcfilters_string)
            lcfilterhash = hashlib.md5(lcfilters_string).hexdigest()

            return (hatlc, 'filters_noresults', lcfilterhash)

    else:

        if LOGGER:
            LOGGER.warning('lcutils.filter_consolidated_hatlc: '
                           'no valid filters could be parsed')
        if DEBUGMODE:
            print('lcutils.filter_consolidated_hatlc: '
                  'no valid filters could be parsed')

        return (hatlc, 'filters_invalid', None)
