#!/usr/bin/env python

'''
lcutils.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Apr 2014

Contains various useful tools for post-processing consolidated HAT LCs.

'''

import logging
import os
import os.path
import subprocess
import shlex
import pwd
import glob
import hashlib
import random
import cPickle as pickle

import numpy as np

#############
## LOGGING ##
#############

# setup a logger
LOGGER = logging.getLogger('lcutils_processing')
LOGGER.addHandler(logging.NullHandler())

# default to debug mode = False
DEBUGMODE = False

def set_debug(debugbool):
    globals()['DEBUGMODE'] = debugbool
    lcform.set_debug(True)
    hlc.set_debug(True)



###################
## LOCAL IMPORTS ##
###################

from timeutils_temp import jd_to_bjd, jd_to_hjd

from lcutils_config import *
import lcutils_formatters as lcform
import lcutils_hatlc as hlc

import lcdb


########################
## USEFUL DEFINITIONS ##
########################

LC_MAG_COLUMNS = ('RM1','RM2','RM3',
                  'EP1','EP2','EP3',
                  'TF1','TF2','TF3',
                  'IRM1','IRM2','IRM3',
                  'IEP1','IEP2','IEP3',
                  'ITF1','ITF2','ITF3')


###############
## FUNCTIONS ##
###############

def find_lc_timegroups(lctimes, mingap=4.0):
    '''
    This finds the gaps in the lightcurve, so we can figure out which times are
    for consecutive observations and which represent gaps between
    seasons.

    lctimes is assumed to be in some form of JD.

    min_gap defines how much the difference between consecutive measurements is
    allowed to be to consider them as parts of different timegroups. By default
    it is set to 4.0 days.

    Returns number of groups and Python slice objects for each group like so:

    (ngroups, [slice(start_ind_1, end_ind_1), ...])

    '''

    lc_time_diffs = [(lctimes[x] - lctimes[x-1]) for x in range(1,len(lctimes))]
    lc_time_diffs = np.array(lc_time_diffs)

    group_start_indices = np.where(lc_time_diffs > mingap)[0]

    if len(group_start_indices) > 0:

        group_indices = []

        for i, gindex in enumerate(group_start_indices):

            if i == 0:
                group_indices.append(slice(0,gindex+1))
            else:
                group_indices.append(slice(group_start_indices[i-1]+1,gindex+1))


        # at the end, add the slice for the last group to the end of the times
        # array
        group_indices.append(slice(group_start_indices[-1]+1,len(lctimes)))

    # if there's no large gap in the LC, then there's only one group to worry
    # about
    else:
        group_indices = [slice(0,len(lctimes))]


    return len(group_indices), group_indices


def normalize_mags_to_group_median(lcdict, timecol='RJD', magcols='all', mingap=4.0, returndict=False,normtomedian=True):
    '''
    This normalizes all the magnitudes in a lightcurve to the median magnitude
    of the entire lightcurve. It does this by finding all the 'timegroups' in
    the lightcurve (observation date-series with no gaps larger than mingap),
    normalizing all of those to 0.0, then adding the offset of the full
    lightcurve median magnitude. nans in the mag columns will be ignored.

    timecol = set to the time column to use when finding timegroups

    magcols = set to the columns to normalize mags for. if 'all', normalizes all
              magnitude columns found in the lcdict. otherwise, use a
              comma-separated list of magnitude columns to normalize

    returndict = if True, returns the LC dict for convenience (all operations
                 are done in place, so this isn't necessary)

    Returns an lcdict with the mag columns replaced with their normalized
    versions.

    '''

    # first, get the LC timegroups
    if timecol in lcdict:
        times = lcdict[timecol]
    elif 'BJD' in lcdict:
        times = lcdict['BJD']
    elif 'HJD' in lcdict:
        times = lcdict['HJD']
    elif 'RJD' in lcdict:
        times = lcdict['RJD']

    # if there aren't any time columns in this lcdict, then we can't do any
    # normalization, return it as-is
    else:
        return lcdict

    ngroups, timegroups = find_lc_timegroups(times, mingap=mingap)

    # next, find all the mag columns to normalize
    if magcols == 'all':
        cols_to_normalize = LC_MAG_COLUMNS
    elif magcols == 'redmags':
        cols_to_normalize = ['RM1','RM2','RM3']
    elif magcols == 'epdmags':
        cols_to_normalize = ['EP1','EP2','EP3']
    elif magcols == 'tfamags':
        cols_to_normalize = ['TF1','TF2','TF3']
    elif magcols == 'epdtfa':
        cols_to_normalize = ['EP1','EP2','EP3','TF1','TF2','TF3']
    else:
        cols_to_normalize = magcols.split(',')
        cols_to_normalize = [x.strip() for x in cols_to_normalize]

    # now, normalize each column
    for col in cols_to_normalize:

        if col in LC_MAG_COLUMNS and col in lcdict:

            mags = lcdict[col]

            # find all the non-nan indices
            finite_ind = np.isfinite(mags)

            # find the global median
            global_mag_median = np.median(mags[finite_ind])

            # go through the groups and normalize them to the median for each
            # group
            for tgind, tg in enumerate(timegroups):

                finite_ind = np.isfinite(mags[tg])

                # find this timegroup's median mag and normalize the mags in it
                # to this median
                group_median = np.median((mags[tg])[finite_ind])
                mags[tg] = mags[tg] - group_median

                if DEBUGMODE:
                    print('group %s: elems %s, finite elems %s, median mag %s' %
                          (tgind,
                           len(mags[tg]),
                           len(finite_ind),
                           group_median))


            # now that everything is normalized to 0.0, add the global median
            # offset back to all the mags and write the result back to the dict
            if normtomedian:
                mags = mags + global_mag_median

            lcdict[col] = mags

        else:

            LOGGER.warning('column %s is not present, skipping...' % col)
            if DEBUGMODE:
                print('column %s is not present, skipping...' % col)
            continue

    if returndict:
        return lcdict
