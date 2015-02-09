#!/usr/bin/env python

'''
lcutils.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - May 2013

Contains various useful tools for finding and processing HAT LCs.

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

import MySQLdb as mysql
import psycopg2 as pg

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

from timeutils import jd_to_bjd, jd_to_hjd

from lcutils_config import *
import lcutils_formatters as lcform
import lcutils_hatlc as hlc

import searchutils
import binlc
import lcdb


#########################
## RETRIEVAL FUNCTIONS ##
#########################

def get_hatnet_lc_locations(hatid,
                            photometry_type='AP',
                            project_id=0,
                            duplicate=False,
                            redlc=True,
                            epdlc=True,
                            tfalc=True,
                            database=None):
    '''
    Given a HAT object ID, find all possible locations for a light-curve of the
    object. Searches in the HAT field as well as its neighbors.

    hatid = 'HAT-XXX-YYYYYYY' - the phone number of a HAT object

    photometry_type = 'AP' for aperture photo LCs
                      'ISM' for image-subtraction LCs

    project_id = get LCs in this project ID directory

    duplicate = False -> get only fields designated as not being duplicates
                True -> get LCs from fields regardless of their duplicate status

    redlc = True -> return paths to .rlc files (.ilc for photo_type = 'ISM')
    epdlc = True -> return paths to .epdlc files (.iepd for photo_type = 'ISM')
    tfalc = True -> return paths to .tfalc files (.itfa for photo_type = 'ISM')

    database = LCDB object to reuse an existing and open database connection
               None to set up and use a new database connection

    '''

    EXT_LIST = {'AP':['rlc','epdlc','tfalc'],
                'ISM':['ilc','iepdlc','itfalc']}

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True


    # get the hat field
    hat_field, hat_object = hatid.split('-')[1:]

    # use this query to get paths from the database
    nbr_query = ('select path from field_lc_locations where '
                 '( (photometry_type = %s) and (duplicate = %s) and '
                 '(project_id = %s) and '
                 'field_num = any((select neighbors from hat_fields '
                 'where field_num = %s)::varchar[]) ) or '
                 '( (field_num = %s) and (photometry_type = %s) and '
                 '(project_id = %s) and (duplicate = %s) )')

    cur.execute(nbr_query, (photometry_type,
                            duplicate,
                            project_id,
                            hat_field,
                            hat_field,
                            photometry_type,
                            project_id,
                            duplicate))
    nbr_paths = cur.fetchall()

    # check if the actual field is in the nbr_path list
    if len(nbr_paths) > 0:

        check_paths = zip(*nbr_paths)
        check_field_paths = check_paths[0]
        field_self_present = False

        for fld in check_paths:
            if hat_field in fld:
                field_self_present = True


    if len(nbr_paths) > 0:

        nbr_paths = zip(*nbr_paths)

        if LOGGER:
            LOGGER.debug('probable %s LCs for %s at: %s' % (photometry_type,
                                                           hatid,
                                                           nbr_paths[0]))
        if DEBUGMODE:
            print('probable %s LCs for %s at: %s' % (photometry_type,
                                                     hatid,
                                                     nbr_paths[0]))

        lc_files = []

        for lc_path in nbr_paths[0]:

            search_exts = EXT_LIST[photometry_type]
            for ext in search_exts:
                lc_files.append('%s/%s.%s' % (lc_path,hatid,ext))

    elif len(nbr_paths) == 0 or field_self_present == False:

        lc_files = []

        if LOGGER:
            LOGGER.warning('no LC locations '
                           'of type %s found for %s,'
                           ' trying fallback...' % (photometry_type,
                                                    hatid))
        if DEBUGMODE:
            print('lcutils: no LC locations '
                  'of type %s found for %s, '
                  'trying fallback...' % (photometry_type,
                                          hatid))

        nbr_query = 'select neighbors from hat_fields where field_num = %s'
        cur.execute(nbr_query, (hat_field,))
        rows = cur.fetchone()
        fields_to_check = [hat_field] + rows[0]
        field_paths = ['%s/G%s/LC' % (FALLBACK_HN_LCPATH, x)
                       for x in fields_to_check]

        for lc_path in field_paths:

            search_exts = EXT_LIST[photometry_type]
            for ext in search_exts:
                lc_files.append('%s/%s.%s' % (lc_path,hatid,ext))


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return(lc_files)



def get_hatsouth_lc_locations(hatid,
                              rawlc=True,
                              epdlc=True,
                              tfalc=True,
                              focusframes=False,
                              database=None):
    '''
    This finds the hatsouth LCs based on the field and object name. Similar to
    get_hatnet_lc_locations above.

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

    # get the hat field
    hat_field, hat_object = hatid.split('-')[1:]

    # use this query to get the neighbors of the requested field
    nbr_query = ('select field_num, field_name, neighbors '
                 'from hat_fields '
                 'where field_num = %s')

    cur.execute(nbr_query, (hat_field,))

    # get the target field and its neighbors
    hatfields = cur.fetchone()
    tgtfield, tgtfname, tgtnbrs = hatfields

    if tgtfield not in tgtnbrs:
        tgtnbrs.append(tgtfield)

    # do another query to get all the field numbers and names to form the LC
    # paths
    cur.execute('select field_num, field_name from hat_fields '
                'where field_num in %s',
                (tuple(tgtnbrs),))

    allfieldinfo = cur.fetchall()


    if len(allfieldinfo) > 0:

        fieldnums, fieldnames = zip(*allfieldinfo)

        remote_lc_files = []
        local_lc_files = []

        if rawlc:

            rawlc_paths = ['/S/LC/%s/%s/1/%s.rawlc',
                           '/S/LC/%s/%s/2/%s.rawlc',
                           '/S/LC/%s/%s/3/%s.rawlc',
                           '/S/LC/%s/%s/4/%s.rawlc']

            locallc_paths = ['%s/fetched-hs-%s-%s-1-%s.rawlc',
                             '%s/fetched-hs-%s-%s-2-%s.rawlc',
                             '%s/fetched-hs-%s-%s-3-%s.rawlc',
                             '%s/fetched-hs-%s-%s-4-%s.rawlc']

            rawlcs = [lcp % (x, 'object', hatid)
                      for x in fieldnames
                      for lcp in rawlc_paths]
            remote_lc_files.extend(rawlcs)

            locallcs = [lcp % (LCCACHE, x, 'object', hatid)
                        for x in fieldnums
                        for lcp in locallc_paths]
            local_lc_files.extend(locallcs)

            # decide if we want focus frames too
            if focusframes:

                rawlc_paths = ['/S/LC/%s/%s/1/%s.rawlc',
                               '/S/LC/%s/%s/2/%s.rawlc',
                               '/S/LC/%s/%s/3/%s.rawlc',
                               '/S/LC/%s/%s/4/%s.rawlc']

                locallc_paths = ['%s/fetched-hs-%s-%s-1-%s.rawlc',
                                 '%s/fetched-hs-%s-%s-2-%s.rawlc',
                                 '%s/fetched-hs-%s-%s-3-%s.rawlc',
                                 '%s/fetched-hs-%s-%s-4-%s.rawlc']

                rawlcs = [lcp % (x, 'focus', hatid)
                          for x in fieldnames
                          for lcp in rawlc_paths]
                remote_lc_files.extend(rawlcs)

                locallcs = [lcp % (LCCACHE, x, 'focus', hatid)
                            for x in fieldnums
                            for lcp in locallc_paths]
                local_lc_files.extend(locallcs)


        if epdlc:

            epdlc_paths = ['/S/LC/%s/1/%s.epdlc',
                           '/S/LC/%s/2/%s.epdlc',
                           '/S/LC/%s/3/%s.epdlc',
                           '/S/LC/%s/4/%s.epdlc']

            locallc_paths = ['%s/fetched-hs-%s-%s-1-%s.epdlc',
                             '%s/fetched-hs-%s-%s-2-%s.epdlc',
                             '%s/fetched-hs-%s-%s-3-%s.epdlc',
                             '%s/fetched-hs-%s-%s-4-%s.epdlc']

            epdlcs = [lcp % (x, hatid)
                      for x in fieldnums
                      for lcp in epdlc_paths]
            remote_lc_files.extend(epdlcs)

            locallcs = [lcp % (LCCACHE, x, 'object', hatid)
                        for x in fieldnums
                        for lcp in locallc_paths]
            local_lc_files.extend(locallcs)

            # decide if we want focus frames too
            if focusframes:

                epdlc_paths = ['/S/LC/%s/1.focus/%s.epdlc',
                               '/S/LC/%s/2.focus/%s.epdlc',
                               '/S/LC/%s/3.focus/%s.epdlc',
                               '/S/LC/%s/4.focus/%s.epdlc']

                locallc_paths = ['%s/fetched-hs-%s-%s-1-%s.epdlc',
                                 '%s/fetched-hs-%s-%s-2-%s.epdlc',
                                 '%s/fetched-hs-%s-%s-3-%s.epdlc',
                                 '%s/fetched-hs-%s-%s-4-%s.epdlc']

                epdlcs = [lcp % (x, hatid)
                          for x in fieldnames
                          for lcp in epdlc_paths]
                remote_lc_files.extend(epdlcs)

                locallcs = [lcp % (LCCACHE, x, 'focus', hatid)
                            for x in fieldnums
                            for lcp in locallc_paths]
                local_lc_files.extend(locallcs)

        if tfalc:

            tfalc_paths = ['/S/LC/%s/1/%s.tfalc',
                           '/S/LC/%s/2/%s.tfalc',
                           '/S/LC/%s/3/%s.tfalc',
                           '/S/LC/%s/4/%s.tfalc']

            locallc_paths = ['%s/fetched-hs-%s-%s-1-%s.tfalc',
                             '%s/fetched-hs-%s-%s-2-%s.tfalc',
                             '%s/fetched-hs-%s-%s-3-%s.tfalc',
                             '%s/fetched-hs-%s-%s-4-%s.tfalc']

            tfalcs = [lcp % (x, hatid)
                      for x in fieldnums
                      for lcp in tfalc_paths]
            remote_lc_files.extend(tfalcs)

            locallcs = [lcp % (LCCACHE, x, 'object', hatid)
                        for x in fieldnums
                        for lcp in locallc_paths]
            local_lc_files.extend(locallcs)

            # decide if we want focus frames too
            if focusframes:

                tfalc_paths = ['/S/LC/%s/1.focus/%s.tfalc',
                               '/S/LC/%s/2.focus/%s.tfalc',
                               '/S/LC/%s/3.focus/%s.tfalc',
                               '/S/LC/%s/4.focus/%s.tfalc']

                locallc_paths = ['%s/fetched-hs-%s-%s-1-%s.tfalc',
                                 '%s/fetched-hs-%s-%s-2-%s.tfalc',
                                 '%s/fetched-hs-%s-%s-3-%s.tfalc',
                                 '%s/fetched-hs-%s-%s-4-%s.tfalc']

                tfalcs = [lcp % (x, hatid)
                          for x in fieldnames
                          for lcp in tfalc_paths]
                remote_lc_files.extend(tfalcs)

                locallcs = [lcp % (LCCACHE, x, 'focus', hatid)
                            for x in fieldnums
                            for lcp in locallc_paths]
                local_lc_files.extend(locallcs)

    else:

        local_lc_files = []
        remote_lc_files = []
        fieldnames = []
        fieldnums = []

        if LOGGER:
            LOGGER.warning('no LCs found for %s' % hatid)
        if DEBUGMODE:
            print('lcutils: no LCs found for %s' % hatid)


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return (remote_lc_files, local_lc_files, fieldnames, fieldnums)



def scp_fetch_hatsouth_lcs(remote_files,
                           local_files,
                           field_names,
                           field_nums,
                           hs_lcrepo='/S/LC',
                           user=None,
                           randomizehost=False,
                           useinternalnet=True):
    '''
    This fetches the HATSouth lightcurves given the output of
    get_hatsouth_lc_locations above.

    This uses scp to fetch the probable light-curves from the remote servers.

    NOTE: The user needs to have a password-less login to all the remote hosts
          for this function to work well. Otherwise, be prepared to confront a
          large number of password prompts.

    user = None -> use the currently logged in user as the scp user for the
                   remote servers
           'user' -> use the specified username as the scp user

    useinternalnet -> True, uses the internal network to speed up copying over
                      lightcurves

    '''

    if user is None:
        scpuser = FALLBACK_SSHUSER
    else:
        scpuser = user


    if randomizehost:
        # randomly choose a host to pull the files from.
        # could make this choice in a better way (i.e. choose host with least
        # load or with the wanted files on its actual physical disk), but since
        # most of them are guaranteed to be under heavy load, this is probably a
        # decent compromise.
        hslc_host = random.choice(HSLC_HOST)
    else:
        # this will pick phs7 every time
        hslc_host = HSLC_HOST[0]

    if useinternalnet:
        hslc_host = 'j%s' % hslc_host.split('.')[0]


    # first check if the fieldnames and fieldnums are present in the HATSouth LC
    # repository
    cmdline = 'ssh %s@%s "%s"' % (scpuser, hslc_host, 'ls %s' % hs_lcrepo)

    cmdline = shlex.split(cmdline)

    p = subprocess.Popen(cmdline,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p_stdout, p_stderr = p.communicate()

    repocheck_results = p_stdout.split('\n')

    fnames = [str(x) for x in field_names]
    fnums = [str(x) for x in field_nums]

    wanted_fields_set = set(fnames + fnums)
    present_fields_set = set(repocheck_results)

    wanted_and_present = present_fields_set.intersection(wanted_fields_set)
    wanted_and_present = list(wanted_and_present)

    # if at least some LCs are present for the requested fields, then get them
    if len(wanted_and_present) > 0:

        for remlc, loclc in zip(remote_files, local_files):

            # make sure to only get LCs that are in field directories that
            # actually exist on the server
            check_remote_lc = [(True if x in os.path.dirname(remlc) else False)
                               for x in wanted_and_present]
            if any(check_remote_lc):

                cmdline = 'scp %s@%s:%s %s' % (scpuser, hslc_host,
                                               remlc, loclc)

                cmdline = shlex.split(cmdline)

                p = subprocess.Popen(cmdline,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
                p_stdout, p_stderr = p.communicate()

                if p.returncode == 0:
                    if LOGGER:
                        LOGGER.debug('successfully fetched %s' % remlc)
                    if DEBUGMODE:
                        print('successfully fetched %s' % remlc)

                else:
                    if LOGGER:
                        LOGGER.warning('fetch failed, %s' %
                                       (p_stderr.rstrip('\n'),))
                    if DEBUGMODE:
                        print('fetch failed, %s' % (p_stderr.rstrip('\n'),))

    else:
        if LOGGER:
            LOGGER.debug('not fetching any HS LCs since none were found.')
        if DEBUGMODE:
            print('lcutils: not fetching any HS LCs since none were found.')


    return wanted_and_present



def scp_fetch_hatnet_lcs(lc_locations,
                         user=None,
                         destination=None,
                         useinternalnet=True):
    '''
    This uses scp to fetch the probable light-curves from the remote servers.

    NOTE: The user needs to have a password-less login to all the remote hosts
          for this function to work well. Otherwise, be prepared to confront a
          large number of password prompts.

    lc_locations = output from lcutils.get_lc_locations function

    user = None -> use the currently logged in user as the scp user for the
                   remote servers
           'user' -> use the specified username as the scp user


    destination = None -> put light-curves in the lcserver's LCCACHE directory
                  '/directory/where/lcs/will/go' -> a LC destination path

    useinternalnet = True -> use hostnames prefixed with 'j' to make copying go
                     faster by using the internal network, False -> use usual
                     hostname

    '''

    if user is None:
        scpuser = FALLBACK_SSHUSER
    else:
        scpuser = user

    if destination is None:
        lcdest = LCCACHE
    else:
        lcdest = destination

    # extract the hostname from each lightcurve path
    hostnames = []
    for lc in lc_locations:
        hostname = lc.split('nfs')[1].split('/')[1]

        if useinternalnet:
            hostname = 'j%s' % hostname

        hostnames.append(hostname)

    hostpathmap = {}

    fetched_lcs = []

    # make a hostpathmap to check if the files actually exist on the server
    for h, p in zip(hostnames, lc_locations):
        if h in hostpathmap:
            hostpathmap[h].append(p)
        else:
            hostpathmap[h] = []

    # now for each host, get only LCs that are actually present
    for h in hostpathmap:

        lcpaths = ' '.join(hostpathmap[h])
        ssh_cmdline = 'ssh %s@%s "ls %s"' % (scpuser, h, lcpaths)
        ssh_cmdline = shlex.split(ssh_cmdline)

        p = subprocess.Popen(ssh_cmdline,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        p_stdout, p_stderr = p.communicate()

        repocheck_results = p_stdout.split('\n')

        # only get the files that are actually present
        hnloc_set = set(lc_locations)
        repocheck_set = set(repocheck_results)
        wanted_and_present = list(repocheck_set.intersection(hnloc_set))

        # extract the hostname from each available path
        hostnames = []
        for lc in wanted_and_present:

            hostname = lc.split('nfs')[1].split('/')[1]

            if useinternalnet:
                hostname = 'j%s' % hostname

            hostnames.append(hostname)

        if len(hostnames) > 0 and len(wanted_and_present) > 0:

            # now hit each path and get the light-curve over to the light-curve
            # cache
            for host, lcpath in zip(hostnames, wanted_and_present):

                lcfield = os.path.split(os.path.split(lcpath)[0])[-1]

                # copy example:
                # /nfs/phn15/ar0/H/4KRED/AP/4KAP_LC/150/HAT-194-0001265.rlc will
                # be <LCCACHE>/fetched-hn-phn15-194-HAT-194-0001265.rlc on local
                # server
                cmdline = 'scp %s@%s:%s %s/fetched-hn-%s-%s-%s' % (
                    scpuser,
                    host,
                    lcpath,
                    lcdest,
                    host.replace('j',''),
                    lcfield,
                    os.path.basename(lcpath)
                )
                cmdline = shlex.split(cmdline)

                p = subprocess.Popen(cmdline,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
                p_stdout, p_stderr = p.communicate()

                if p.returncode == 0:
                    if LOGGER:
                        LOGGER.debug('successfully fetched %s' % lcpath)
                    if DEBUGMODE:
                        print('successfully fetched %s' % lcpath)
                    fetched_lcs.append(lcpath)

                else:
                    if LOGGER:
                        LOGGER.warning(
                            'fetch failed, %s' % (p_stderr.rstrip('\n'),)
                            )
                    if DEBUGMODE:
                        print('fetch failed, %s' % (p_stderr.rstrip('\n'),))

        else:

            if LOGGER:
                LOGGER.debug('not fetching any HN LCs since none were found.')
            if DEBUGMODE:
                print('lcutils: not fetching any HN LCs since none were found.')

    return fetched_lcs


def collect_fetched_lightcurves(hatid,
                                lcdir=None,
                                outdir=None,
                                removefetched=False,
                                saveinfo=False,
                                ignorecollected=True):
    '''
    This collects the lightcurves copied over from a remote server and placed in
    the LCCACHE directory for a single HAT object. The collection is done in the
    following way:

    * all lightcurves with the same extension are copied into one file of that
      same extension

    * they are then sorted by the time/stid column if possible

    * collected light-curves are returned as a dict of filepaths by lc type


    hatid = collect light-curves for this HAT ID

    lcdir = None -> use LCCACHE as directory to look for object LCs
            '/other/lc/path' -> use this directory to look for object LCs

    outdir = None -> use LCCACHE as directory to hold output LCs
             '/other/lc/path' -> use this directory to hold output LCs

    removefetched = True -> remove all fetched LCs for this object
                    False -> keep the fetched LCs for this object

    saveinfo = True -> saves the collected file dictionary to a pickled file
               False -> does not save anything

    The output LC is of the form:

    collected-<HAT ID>.<rlc|epdlc|tfalc, ..., etc.>

    The output collected file dictionary picked file is of the form:

    <HAT-ID>-collected-lcs.pkl

    '''

    # first, check if this LC has already been collected
    if not ignorecollected:

        if DEBUGMODE:
            print('checking for existing collected LC info file for %s' % hatid)

        checkfname = '%s-collected-lcs.pkl' % hatid
        if lcdir:
            checkpath = os.path.join(lcdir, checkfname)
        else:
            checkpath = os.path.join(LCCACHE, checkfname)

        if os.path.exists(checkpath):

            collectedinfo = open(checkpath,'rb')
            collecteddict = pickle.load(collectedinfo)
            collectedinfo.close()

            if LOGGER:
                LOGGER.warning('%s LCs already collected, returning collected '
                               'LC dict from file %s' % (hatid, checkpath))
            if DEBUGMODE:
                print('%s LCs already collected, returning collected '
                      'LC dict from file %s' % (hatid, checkpath))

            return collecteddict



    # networks to group by
    networks = {'hn':['rlc','epdlc','tfalc',
                      'ilc','iepdlc','itfalc'],
                'hs':['rawlc','epdlc','tfalc',
                      'ilc','iepd','itfa']}

    outdict = {'hn':{},
               'hs':{}}

    # collect by network then by extension
    for net in networks:

        extensions = networks[net]

        for ext in extensions:

            search_path = 'fetched-%s-*%s.%s' % (net, hatid, ext)
            outlc_path = 'collected-%s-%s.%s' % (net, hatid, ext)

            if lcdir is None:
                search_path = os.path.join(LCCACHE,search_path)
            else:
                search_path = os.path.join(lcdir,search_path)

            if outdir is None:
                outlc_path = os.path.join(LCCACHE,outlc_path)
            else:
                outlc_path = os.path.join(outdir,outlc_path)

            # find the LCs for this extension
            extlcs = glob.glob(search_path)

            # exclude already collected LCs
            extlcs = [x for x in extlcs if 'collected' not in x]

            # there are three types of LCs to worry about:
            # 1. HATSouth binary rawlcs
            # 2. HATSouth text epdlcs and tfalcs
            # 3. HATNet text rlcs, epdlcs, and tfalcs

            # 1. collect the HATSouth binary rawlc lightcurves
            if len(extlcs) > 0 and ext == 'rawlc' and net == 'hs':

                if LOGGER:
                    LOGGER.debug('found %s %s LC(s) for %s' % (len(extlcs),
                                                              ext,
                                                              hatid))
                if DEBUGMODE:
                    print('lcutils: found %s %s LC(s) for %s' % (len(extlcs),
                                                        ext,
                                                        hatid))

                outlc_contents = {}

                # get all the columns from all the files together
                for lcind, lc in enumerate(extlcs):

                    rawlc = binlc.read_hat_binlc(lc)

                    if lcind == 0:
                        for col in rawlc['cols']:
                            outlc_contents[col] = rawlc[col]

                    else:
                        for col in rawlc['cols']:
                            outlc_contents[col].extend(rawlc[col])

                    del rawlc


                # now sort the outlc columns by BJD
                sort_col = outlc_contents['BJD']
                unique_sort_col, unique_inds = np.unique(sort_col,
                                                         return_index=True)

                outlc_cols = []
                outlc_lineformat = []

                # prepare the output file
                for col in HATLC_COL_DEFS['hs']['rawlc']:
                    outcol = (np.array(outlc_contents[col]))[unique_inds]
                    outlc_cols.append(outcol)
                    outlc_lineformat.append(TEXTLC_OUTPUT_COLUMNS[col][1])

                outlc_lineformat = ' '.join(outlc_lineformat)
                outlc_cols = zip(*outlc_cols)

                del outlc_contents

                # once the sort is done, write these columns to a collected LC
                # file
                outf = open(outlc_path,'wb')
                for line in outlc_cols:
                    outline = outlc_lineformat % line
                    outline = outline + '\n'
                    outf.write(outline)
                outf.close()

                del outlc_cols

                if LOGGER:
                    LOGGER.debug('wrote collected %s LC to %s' % (ext,
                                                                 outlc_path))
                if DEBUGMODE:
                    print('lcutils: wrote collected %s LC to %s' % (ext,
                                                                    outlc_path))


                # if we're supposed to remove the component LCs, do so
                if removefetched:
                    for lc in extlcs:
                        os.remove(lc)

                outdict[net][ext] = outlc_path


            # 2. collect the HATSouth text epdlc and tfalc lightcurves
            elif len(extlcs) > 0 and ext != 'rawlc' and net == 'hs':

                if LOGGER:
                    LOGGER.debug('found %s %s LC(s) for %s' % (len(extlcs),
                                                              ext,
                                                              hatid))
                if DEBUGMODE:
                    print('lcutils: found %s %s LC(s) for %s' % (len(extlcs),
                                                        ext,
                                                        hatid))

                outlc_contents = []

                # collect all LCs
                for lc in extlcs:

                    lcf = open(lc,'rb')
                    lc_contents = lcf.readlines()

                    # filter out the comment lines
                    lc_contents = [x for x in lc_contents
                                   if not(x.startswith('#'))]

                    # filter out lines longer than the defined lengths for each
                    # LC type
                    if '-hn-' in lc:
                        maxlinelen = TEXTLC_LINE_LENGTHS[ext]
                    elif '-hs-' in lc and ext == 'epdlc':
                        maxlinelen = TEXTLC_LINE_LENGTHS['hsepd']
                    elif '-hs-' in lc and ext == 'tfalc':
                        maxlinelen = TEXTLC_LINE_LENGTHS['hstfa']

                    lc_contents = [x for x in lc_contents
                                   if len(x) < maxlinelen]

                    lcf.close()
                    outlc_contents.extend(lc_contents)

                    del lc_contents

                # now that we've collected all the LC lines, we need to sort
                # them, the sort column for HS text LCs is column 4 (the BJD)
                sort_col_ind = 3

                if len(outlc_contents) > 0:

                    split_outlc_contents = [x.split() for x in outlc_contents]
                    transposed_outlc_contents = zip(*split_outlc_contents)

                    # get the unique sorted indices of the sort column
                    sort_col = np.array(transposed_outlc_contents[sort_col_ind])
                    unique_sort_col, unique_inds = np.unique(sort_col,
                                                             return_index=True)

                    # sort the lines in the same order using the indices of the
                    # unique elements only
                    outlc_contents = (np.array(outlc_contents))[unique_inds]
                    outlc_contents = [x for x in outlc_contents]

                    # now write the collected LC to the output file
                    outf = open(outlc_path,'wb')
                    for line in outlc_contents:
                        outf.write(line)
                    outf.close()

                    del outlc_contents
                    del split_outlc_contents
                    del transposed_outlc_contents

                    if LOGGER:
                        LOGGER.debug('wrote collected %s LC to %s' % (ext,
                                                                     outlc_path))
                    if DEBUGMODE:
                        print('lcutils: wrote collected %s LC to %s' % (ext,
                                                                        outlc_path))

                    outdict[net][ext] = outlc_path

                else:

                    if LOGGER:
                        LOGGER.warning('no HS LC found for %s' % hatid)
                    if DEBUGMODE:
                        print('lcutils: no HS LC found for %s' % hatid)


                # if we're supposed to remove the component LCs, do so
                if removefetched:
                    for lc in extlcs:
                        os.remove(lc)



            # 3. collect the HATNet text rlc, epdlc, and tfalc lightcurves
            elif len(extlcs) > 0 and ext != 'rawlc' and net == 'hn':

                if LOGGER:
                    LOGGER.debug('found %s %s LC(s) for %s' % (len(extlcs),
                                                              ext,
                                                              hatid))
                if DEBUGMODE:
                    print('lcutils: found %s %s LC(s) for %s' % (len(extlcs),
                                                                 ext,
                                                                 hatid))

                outlc_contents = []

                # collect all LCs
                for lc in extlcs:

                    lcf = open(lc,'rb')
                    lc_contents = lcf.readlines()
                    lcf.close()

                    # filter out the comment lines
                    lc_contents = [x for x in lc_contents
                                   if not(x.startswith('#'))]

                    # set the longphoto to false initially
                    longphoto = False

                    if len(lc_contents) == 0:
                        if DEBUGMODE:
                            print('%s appears to be empty, skipping...' %
                                  lc)
                        continue

                    # check if this is the new version of the lightcurves using
                    # the HS text LC format.

                    # first, pick a random line from lc_contents
                    lc_randomline = random.choice(lc_contents)

                    # then, check if '_' is in the first column (stid-framenum
                    # key)
                    lc_randomline = lc_randomline.split()

                    if ext == 'rlc':
                        hs_keycolind = 0
                    elif ext == 'epdlc':
                        hs_keycolind = 1
                    elif ext == 'tfalc':
                        hs_keycolind = 1
                    else:
                        hs_keycolind = 0

                    if '_' in lc_randomline[hs_keycolind]:
                        hsformat = True
                        if DEBUGMODE:
                            print('%s is in hsformat' % lc)
                    else:
                        hsformat = False

                    # check if any LC line contains long version photometry
                    # flags
                    for line in lc_contents:
                        if 'G---' in line:
                            longphoto = True
                            break

                    # decide how many characters a line is supposed to have
                    if ext == 'rlc':
                        if not longphoto and not hsformat:
                            maxlinelen = TEXTLC_LINE_LENGTHS[ext]
                        elif longphoto and not hsformat:
                            maxlinelen = 260
                        elif not longphoto and hsformat:
                            maxlinelen = TEXTLC_LINE_LENGTHS['hsrlc']

                    elif ext == 'epdlc':
                        if hsformat:
                            maxlinelen = TEXTLC_LINE_LENGTHS['hsepd']
                        else:
                            maxlinelen = TEXTLC_LINE_LENGTHS[ext]

                    elif ext == 'tfalc':
                        if hsformat:
                            maxlinelen = TEXTLC_LINE_LENGTHS['hstfa']
                        else:
                            maxlinelen = TEXTLC_LINE_LENGTHS[ext]

                    # handle ISM LCs
                    elif ext in ('ilc','iepdlc','itfalc'):
                        maxlinelen = TEXTLC_LINE_LENGTHS[ext]


                    # now filter out lines longer than allowed
                    lc_contents = [x for x in lc_contents if
                                   len(x) < maxlinelen]

                    # if this is an ISM LC, replace all NULLs with NaNs
                    if ext in ('ilc','iepdlc','itfalc'):
                        lc_contents = [x.replace('NULL','nan') for
                                       x in lc_contents]

                    split_lc_contents = [x.split() for x in lc_contents]


                    # here, check if we're using a tfalc in HN format and if the
                    # number of columns is ok
                    if ext == 'tfalc' and not hsformat:

                        lc_line_ncols = [len(x) for x in split_lc_contents]

                        # this means we're using a 7-col tfalc, instead of the
                        # usual 9-col tfalc, and we need to pad with NaNs
                        if all([(x == 18) for x in lc_line_ncols]):

                            split_lc_contents = [
                                (x[:17] + ['NaN','NaN'] + [x[-1]]) for x in
                                split_lc_contents
                                ]
                            lc_contents = ['%s\n' % (' '.join(x))
                                           for x in split_lc_contents]


                    # if this is an hsformat HN LC, we need to move some columns
                    # around and do some formatting to get them to match all the
                    # others
                    if hsformat:

                        # get the column definitions
                        lc_col_definitions = HATLC_COL_DEFS['hs'][ext]

                        hsformat_lc_cols = [x.split() for x in lc_contents]
                        hsformat_lc_cols = zip(*hsformat_lc_cols)
                        hnformat_lc_cols = []

                        if ext == 'rlc':

                            # get the columns except for the first one
                            cols_to_get = HATLC_COL_DEFS['hn']['rlc'][2:]

                            for col in cols_to_get:
                                hnformat_lc_cols.append(
                                    hsformat_lc_cols[
                                        lc_col_definitions.index(col)
                                        ]
                                    )

                            # fix the RSTF column to remove the _1 part and add
                            # that column in
                            rstf_col = (
                                hsformat_lc_cols[
                                    lc_col_definitions.index('RSTFC')
                                    ]
                                )
                            rstf_col = [x.rstrip('_1') for x in rstf_col]
                            hnformat_lc_cols.insert(0,rstf_col)

                            # add the HAT-ID column in
                            hatid_col = [
                                hatid for x in range(len(hnformat_lc_cols[0]))
                                ]
                            hnformat_lc_cols.insert(0,hatid_col)

                            # now format back to a list of lines
                            hnformat_lc_rows = zip(*hnformat_lc_cols)
                            lc_contents = ['%s\n' % (' '.join(x))
                                           for x in hnformat_lc_rows]

                        elif ext == 'epdlc':

                            # get the columns except for the first one
                            cols_to_get = HATLC_COL_DEFS['hn']['epdlc'][1:]
                            for col in cols_to_get:
                                hnformat_lc_cols.append(
                                    hsformat_lc_cols[
                                        lc_col_definitions.index(col)
                                        ]
                                    )
                            # fix the RSTF column to remove the _1 part and add
                            # that column in
                            rstf_col = (
                                hsformat_lc_cols[
                                    lc_col_definitions.index('ESTFC')
                                    ]
                                )
                            rstf_col = [x.rstrip('_1') for x in rstf_col]

                            hnformat_lc_cols.insert(0,rstf_col)

                            # now format back to a list of lines
                            hnformat_lc_rows = zip(*hnformat_lc_cols)
                            lc_contents = ['%s\n' % (' '.join(x))
                                           for x in hnformat_lc_rows]

                        elif ext == 'tfalc':

                            # get the columns except for the first one
                            cols_to_get = HATLC_COL_DEFS['hn']['tfalc'][1:]
                            for col in cols_to_get:
                                hnformat_lc_cols.append(
                                    hsformat_lc_cols[
                                        lc_col_definitions.index(col)
                                        ]
                                    )
                            # fix the RSTF column to remove the _1 part and add
                            # that column in
                            rstf_col = (
                                hsformat_lc_cols[
                                    lc_col_definitions.index('TSTFC')
                                    ]
                                )
                            rstf_col = [x.rstrip('_1') for x in rstf_col]
                            hnformat_lc_cols.insert(0,rstf_col)

                            # now format back to a list of lines
                            hnformat_lc_rows = zip(*hnformat_lc_cols)
                            lc_contents = ['%s\n' % (' '.join(x))
                                           for x in hnformat_lc_rows]


                    # now extend the final lightcurve line contents by the just
                    # read lines
                    outlc_contents.extend(lc_contents)


                ## END OF PER LC PROCESSING
                # now that we've collected all the LC lines, we need to sort and
                # uniquify them
                if ext == 'rlc':
                    sort_col_ind = 1
                elif ext in ('epdlc','tfalc','ilc','iepdlc','itfalc'):
                    sort_col_ind = 0

                # make sure that outlc_contents is not empty
                if len(outlc_contents) > 0:

                    split_outlc_contents = [x.split() for x in outlc_contents]

                    # we need to guard against different numbers of columns
                    # potentially truncating things (e.g. most tfalcs have 9
                    # cols after the instrumental mags, which is the default,
                    # but some have only 7 columns after the instrumental
                    # mags). we should check the number of cols against the
                    # extension, and remove the lines which have ncols different
                    # from what is expected.

                    # filter out lines that don't match expected column numbers
                    split_outlc_contents = [
                        x for x in split_outlc_contents
                        if len(x) == len(HATLC_COL_DEFS[net][ext])
                    ]

                    # reform outlc_contents after this filtering
                    outlc_contents = ['%s\n' % (' '.join(x))
                                      for x in split_outlc_contents]

                    # make sure outlc_contents is not empty after ncol filtering
                    if len(outlc_contents) > 0:

                        # transpose to check for sorting by date and to check for
                        # long format photometry flags
                        transposed_outlc_contents = zip(*split_outlc_contents)

                        # also check for long photometry flags and strip to first
                        # character if found
                        if longphoto and ext == 'rlc':

                            # get the indices for the photometry flags
                            iq1flag_index = HATLC_COL_DEFS[net][ext].index('IQ1')
                            iq2flag_index = HATLC_COL_DEFS[net][ext].index('IQ2')
                            iq3flag_index = HATLC_COL_DEFS[net][ext].index('IQ3')

                            iq1_flags = transposed_outlc_contents[iq1flag_index]
                            iq2_flags = transposed_outlc_contents[iq2flag_index]
                            iq3_flags = transposed_outlc_contents[iq3flag_index]

                            transposed_outlc_contents[iq1flag_index] = (
                                [x[0] for x in iq1_flags]
                                )
                            transposed_outlc_contents[iq2flag_index] = (
                                [x[0] for x in iq2_flags]
                                )
                            transposed_outlc_contents[iq3flag_index] = (
                                [x[0] for x in iq3_flags]
                                )

                            outlc_contents = zip(*transposed_outlc_contents)
                            outlc_contents = ['%s\n' % (' '.join(x))
                                              for x in outlc_contents]

                        # get the unique sorted indices of the sort column
                        sort_col = np.array(transposed_outlc_contents[sort_col_ind])
                        unique_sort_col, unique_inds = np.unique(sort_col,
                                                                 return_index=True)

                        # sort the lines in the same order using the indices of the
                        # unique elements only
                        outlc_contents = (np.array(outlc_contents))[unique_inds]
                        outlc_contents = outlc_contents.tolist()

                        # now write the collected LC to the output file
                        outf = open(outlc_path,'wb')
                        for line in outlc_contents:
                            outf.write(line)
                        outf.close()

                        del outlc_contents
                        del split_outlc_contents
                        del transposed_outlc_contents

                        if LOGGER:
                            LOGGER.debug(
                                'wrote collected %s LC to %s' % (
                                    ext,
                                    outlc_path
                                    )
                                )
                        if DEBUGMODE:
                            print('lcutils: wrote collected %s LC to %s' %
                                  (ext,outlc_path))


                        # if we're supposed to remove the component LCs, do so
                        if removefetched:
                            for lc in extlcs:
                                os.remove(lc)

                        outdict[net][ext] = outlc_path

                    else:
                        if LOGGER:
                            LOGGER.debug(
                                'collecting %s LC failed '
                                'because ncols is weird! skipping...' % ext
                                )
                        if DEBUGMODE:
                            print(
                                'collecting %s LC failed '
                                'because ncols is weird! skipping...' % ext
                                )

    outdict['hatid'] = hatid

    # if we're supposed to save collected info for this HAT-ID
    if (outdict['hn'] or outdict['hs']) and saveinfo:

        # dump this dictionary to a pickle format file to be read and used
        # by get_hatlcs
        pklfname = '%s-collected-lcs.pkl' % hatid
        pklfpath = os.path.join(LCCACHE, pklfname)
        pklf = open(pklfpath,'wb')
        pickle.dump(outdict, pklf)
        pklf.close()


    return outdict



def consolidate_hatnet_lightcurves(collected_lc_dict,
                                   removecollected=False):
    '''
    This pulls together HN LCs into a datatable.

    '''

    datatable = []

    # get all of the information from all of the lightcurves
    # leave the date in FJD form since we can easily convert it to any other
    # format later
    lc_cols = ['RSTF','ESTF','TSTF',
               'XCC','YCC','BGV','BGE','FSV','FDV','FKV',
               'IM1','IE1','IQ1','IM2','IE2','IQ2','IM3','IE3','IQ3',
               'RM1','RM2','RM3',
               'EP1','EP2','EP3',
               'TF1','TF2','TF3']

    # connect to the HATMASTER database
    db = mysql.connect(user=MYSQL_USER_HN,
                       passwd=MYSQL_PASS_HN,
                       db=MYSQL_DATA_HN,
                       host=MYSQL_HOST_HN)
    cur = db.cursor()

    hnsources = collected_lc_dict['hn']
    hn_stfs = []
    lcdata = {}

    # get all data into a dict
    for source in hnsources:

        lcf = open(hnsources[source],'rb')
        databuf = lcf.readlines()
        databuf = [x.split() for x in databuf]
        databuf = zip(*databuf)
        lcf.close()

        lcdata[source] = {}

        # get all the data from this source, taking care to ignore columns
        # that are not needed
        for colind, col in enumerate(HATLC_COL_DEFS['hn'][source]):
            if col is not None:
                lcdata[source][col] = np.array(databuf[colind])

        # get the unique stid-framenum keys from each type of LC
        if source in ('rlc','ilc'):
            hn_stfs.extend(lcdata[source]['RSTF'].tolist())
        elif source in ('epdlc','iepdlc'):
            hn_stfs.extend(lcdata[source]['ESTF'].tolist())
        elif source in ('tfalc','itfalc'):
            hn_stfs.extend(lcdata[source]['TSTF'].tolist())

    # these are the unique stf keys
    hn_stfs = list(set(hn_stfs))

    if DEBUGMODE:
        print('%s total unique detections across all LCs of type: %s' %
              (len(hn_stfs), ', '.join(hnsources.keys())))

    # now lookup all the JDs corresponding to all the station ids and frame
    # numbers from the HATMASTER database

    # split the stid-framenum key, while taking care to remove any weird
    # keys
    stidframenum = [x.split('-') for x in hn_stfs if '-' in x]
    stids, framenums = (np.array([int(x[0]) for x in stidframenum]),
                      np.array([int(x[1]) for x in stidframenum]))

    unique_stids = np.unique(stids)

    query_constraints = []
    query_conditions = []

    for st in unique_stids:

        st_frn_ind = np.where(stids == st)

        # we need to do this because MySQLdb doesn't map tuples directly to
        # SQL arrays as far as I know
        framenum_condition = '%s,' * len(st_frn_ind[0])
        framenum_condition = framenum_condition.strip(',')

        constraint = [st] + [int(x) for x in framenums[st_frn_ind]]
        query_constraints.extend(constraint)

        framenum_condition = 'IMfnum in (%s)' % framenum_condition
        full_condition = '(IMstid = %s and ' + framenum_condition + ')'
        query_conditions.append(full_condition)

    query_constraints = tuple(query_constraints)

    # for HATNet we need to get the frame JD, filter, field name, CCD, HA
    # and Z from the HATMASTER database
    stf_query = ('select a.IMstid, a.IMfnum, a.IMjd, a.IMfilid, '
                 'a.IMobject, a.IMccd, b.IAnha, b.IAnz, a.IMtexp '
                 'from Images a join ImAstrom b on '
                 '((a.IMstid = b.IAstid) and '
                 '(a.IMfnum = b.IAfnum) and (a.IMccd = b.IAccd)) where %s' %
                 (' or '.join(query_conditions)))

    if LOGGER:
        LOGGER.debug('querying HATMASTER DB for metadata...')
    if DEBUGMODE:
        print('lcutils: querying HATMASTER DB for metadata...')

    cur.execute(stf_query, query_constraints)
    jd_rows = cur.fetchall()

    if DEBUGMODE:
        print('%s rows in HATMASTER DB '
              'corresponding to these detections.' % len(jd_rows))

    cur.close()
    db.close()

    stfs = np.array(['-'.join([str(x[0]),str(x[1])]) for x in jd_rows])
    rjds = np.array([x[2] for x in jd_rows])
    filters = np.array([x[3] for x in jd_rows])
    fields = np.array([x[4] for x in jd_rows])
    ccds = np.array([x[5] for x in jd_rows])
    hourangles = np.array([x[6] for x in jd_rows])
    zenithdists = np.array([x[7] for x in jd_rows])
    exptimes = np.array([x[8] for x in jd_rows])

    # sort by rjd
    rjd_sort_ind = np.argsort(rjds)
    stfs = stfs[rjd_sort_ind]
    rjds = rjds[rjd_sort_ind]
    filters = filters[rjd_sort_ind]
    fields = fields[rjd_sort_ind]
    ccds = ccds[rjd_sort_ind]
    hourangles = hourangles[rjd_sort_ind]
    zenithdists = zenithdists[rjd_sort_ind]
    exptimes = exptimes[rjd_sort_ind]

    if LOGGER:
        LOGGER.info('consolidating HN lightcurves...')
    if DEBUGMODE:
        print('lcutils: consolidating HN lightcurves...')

    # get the columns from the rlc LC matching all the STF keys
    if 'rlc' in lcdata and 'RSTF' in lcdata['rlc']:
        # make a mask and find indices where rlc stfs == all stfs
        rlc_stf_mask = np.in1d(lcdata['rlc']['RSTF'],
                               stfs,
                               assume_unique=True)
    else:
        rlc_stf_mask = None

    # get the columns from the epdlc LC matching all the STF keys
    if 'epdlc' in lcdata and 'RSTF' in lcdata['epdlc']:
        # make a mask and find indices where epdlc stfs == all stfs
        epdlc_stf_mask = np.in1d(lcdata['epdlc']['RSTF'],
                               stfs,
                               assume_unique=True)
    else:
        epdlc_stf_mask = None

    # get the columns from the tfalc LC matching all the STF keys
    if 'tfalc' in lcdata and 'RSTF' in lcdata['tfalc']:
        # make a mask and find indices where tfalc stfs == all stfs
        tfalc_stf_mask = np.in1d(lcdata['tfalc']['RSTF'],
                               stfs,
                               assume_unique=True)
    else:
        tfalc_stf_mask = None

    # get the columns from the ilc LC matching all the STF keys
    if 'ilc' in lcdata and 'RSTF' in lcdata['ilc']:
        # make a mask and find indices where ilc stfs == all stfs
        ilc_stf_mask = np.in1d(lcdata['ilc']['RSTF'],
                               stfs,
                               assume_unique=True)
    else:
        ilc_stf_mask = None

    # get the columns from the iepdlc LC matching all the STF keys
    if 'iepdlc' in lcdata and 'RSTF' in lcdata['iepdlc']:
        # make a mask and find indices where iepdlc stfs == all stfs
        iepdlc_stf_mask = np.in1d(lcdata['iepdlc']['RSTF'],
                               stfs,
                               assume_unique=True)
    else:
        iepdlc_stf_mask = None

    # get the columns from the itfalc LC matching all the STF keys
    if 'itfalc' in lcdata and 'RSTF' in lcdata['itfalc']:

        # make a mask and find indices where itfalc stfs == all stfs
        itfalc_stf_mask = np.in1d(lcdata['itfalc']['RSTF'],
                               stfs,
                               assume_unique=True)
    else:
        itfalc_stf_mask = None


    # global metadata
    obsfield = [(x.split('_')[-1] if '_' in x else x) for
                x in fields]
    ha = [(TEXTLC_OUTPUT_COLUMNS['IHA'][1] % x.strip() if x else 'NaN')
          for x in hourangles]
    zd = [(TEXTLC_OUTPUT_COLUMNS['IZD'][1] % x.strip() if x else 'NaN')
          for x in zenithdists]
    exp = [(TEXTLC_OUTPUT_COLUMNS['EXP'][1] % x.strip() if x else 'NaN')
          for x in zenithdists]

    # now get the data
    # CCD METADATA ONLY FROM RLC


    return lcdata, stfs



def consolidate_all_lightcurves(collected_lc_dict,
                                removecollected=False):
    '''
    This pulls together all of the columns from all of the LCs in
    collected_lc_dict and writes out a data table with all of the info in
    appropriate order.

    if removecollected = True, will remove the LCs in collected_lc_dict, and the
    collected LC info pickle file as well

    '''

    hatid = collected_lc_dict['hatid']

    # we put together the final data-table containing all of the data (even if
    # it's missing)
    datatable = []

    ################################################
    # FIRST: collect all the available HATNet data #
    ################################################

    if collected_lc_dict['hn']:

        # get all of the information from all of the lightcurves
        # leave the date in FJD form since we can easily convert it to any other
        # format later
        lc_cols = ['RSTF','ESTF','TSTF',
                   'XCC','YCC','BGV','BGE','FSV','FDV','FKV',
                   'IM1','IE1','IQ1','IM2','IE2','IQ2','IM3','IE3','IQ3',
                   'RM1','RM2','RM3',
                   'EP1','EP2','EP3',
                   'TF1','TF2','TF3']

        # connect to the HATMASTER database
        db = mysql.connect(user=MYSQL_USER_HN,
                           passwd=MYSQL_PASS_HN,
                           db=MYSQL_DATA_HN,
                           host=MYSQL_HOST_HN)
        cur = db.cursor()

        hnsources = collected_lc_dict['hn']
        hn_stfs = []
        lcdata = {}

        # get all data into a dict
        for source in hnsources:

            lcf = open(hnsources[source],'rb')
            databuf = lcf.readlines()
            databuf = [x.split() for x in databuf]
            databuf = zip(*databuf)
            lcf.close()

            lcdata[source] = {}

            # get all the data from this source, taking care to ignore columns
            # that are not needed
            for colind, col in enumerate(HATLC_COL_DEFS['hn'][source]):

                if col is not None:
                    lcdata[source][col] = np.array(databuf[colind])

            # get the unique stid-framenum keys from each type of LC
            if source in ('rlc','ilc'):
                hn_stfs.extend(lcdata[source]['RSTF'].tolist())
            elif source in ('epdlc','iepdlc'):
                hn_stfs.extend(lcdata[source]['ESTF'].tolist())
            elif source in ('tfalc','itfalc'):
                hn_stfs.extend(lcdata[source]['TSTF'].tolist())

        # these are the unique stf keys
        hn_stfs = [x.strip() for x in hn_stfs]
        hn_stfs = list(set(hn_stfs))

        if DEBUGMODE:
            print('%s total unique detections across all LCs of type: %s' %
                  (len(hn_stfs), ', '.join(hnsources.keys())))

        # now lookup all the JDs corresponding to all the station ids and frame
        # numbers from the HATMASTER database

        # split the stid-framenum key, while taking care to remove any weird
        # keys
        stidframenum = [x.split('-') for x in hn_stfs if '-' in x]
        stids, framenums = (np.array([int(x[0]) for x in stidframenum]),
                          np.array([int(x[1]) for x in stidframenum]))

        unique_stids = np.unique(stids)

        query_constraints = []
        query_conditions = []

        for st in unique_stids:

            st_frn_ind = np.where(stids == st)

            # we need to do this because MySQLdb doesn't map tuples directly to
            # SQL arrays as far as I know
            framenum_condition = '%s,' * len(st_frn_ind[0])
            framenum_condition = framenum_condition.strip(',')

            constraint = [st] + [int(x) for x in framenums[st_frn_ind]]
            query_constraints.extend(constraint)

            framenum_condition = 'IMfnum in (%s)' % framenum_condition
            full_condition = '(IMstid = %s and ' + framenum_condition + ')'
            query_conditions.append(full_condition)

        query_constraints = tuple(query_constraints)

        # for HATNet we need to get the frame JD, filter, field name, CCD, HA
        # and Z from the HATMASTER database
        stf_query = ('select a.IMstid, a.IMfnum, a.IMjd, a.IMfilid, '
                     'a.IMobject, a.IMccd, b.IAnha, b.IAnz, a.IMtexp '
                     'from Images a join ImAstrom b on '
                     '((a.IMstid = b.IAstid) and '
                     '(a.IMfnum = b.IAfnum) and (a.IMccd = b.IAccd)) where %s' %
                     (' or '.join(query_conditions)))

        if LOGGER:
            LOGGER.debug('querying HATMASTER DB for metadata...')
        if DEBUGMODE:
            print('lcutils: querying HATMASTER DB for metadata...')

        cur.execute(stf_query, query_constraints)
        jd_rows = cur.fetchall()

        if DEBUGMODE:
            print('%s rows in HATMASTER DB '
                  'corresponding to these detections.' % len(jd_rows))

        cur.close()
        db.close()

        stfs = np.array(['-'.join([str(x[0]),str(x[1])]) for x in jd_rows])
        rjds = np.array([x[2] for x in jd_rows])
        filters = np.array([x[3] for x in jd_rows])
        fields = np.array([x[4] for x in jd_rows])
        ccds = np.array([x[5] for x in jd_rows])
        hourangles = np.array([x[6] for x in jd_rows])
        zenithdists = np.array([x[7] for x in jd_rows])
        exptimes = np.array([x[8] for x in jd_rows])

        # sort by rjd
        rjd_sort_ind = np.argsort(rjds)
        stfs = stfs[rjd_sort_ind]
        rjds = rjds[rjd_sort_ind]
        filters = filters[rjd_sort_ind]
        fields = fields[rjd_sort_ind]
        ccds = ccds[rjd_sort_ind]
        hourangles = hourangles[rjd_sort_ind]
        zenithdists = zenithdists[rjd_sort_ind]
        exptimes = exptimes[rjd_sort_ind]

        if LOGGER:
            LOGGER.info('consolidating HN lightcurves...')
        if DEBUGMODE:
            print('lcutils: consolidating HN lightcurves...')

        for data_ind, data_stf in enumerate(stfs):

            obsfield = fields[data_ind]
            if '_' in obsfield:
                obsfield = obsfield.split('_')[-1]

            ha = ((TEXTLC_OUTPUT_COLUMNS['IHA'][1] %
                   hourangles[data_ind]).strip()
                  if hourangles[data_ind] else 'NaN')
            zd = ((TEXTLC_OUTPUT_COLUMNS['IZD'][1] %
                   zenithdists[data_ind]).strip()
                  if zenithdists[data_ind] else 'NaN')
            exp = ((TEXTLC_OUTPUT_COLUMNS['EXP'][1] %
                    exptimes[data_ind]).strip()
                   if exptimes[data_ind] else 'NaN')

            stid, framenum = stfs[data_ind].split('-')

            # get the data from any available source

            # AP LCs
            if 'rlc' in lcdata:
                rstf_ind = np.where(lcdata['rlc']['RSTF'] == data_stf)
            else:
                rstf_ind = ([],)

            if 'epdlc' in lcdata:
                estf_ind = np.where(lcdata['epdlc']['ESTF'] == data_stf)
            else:
                estf_ind = ([],)

            if 'tfalc' in lcdata:
                tstf_ind = np.where(lcdata['tfalc']['TSTF'] == data_stf)
            else:
                tstf_ind = ([],)

            # ISM LCs
            if 'ilc' in lcdata:
                irstf_ind = np.where(lcdata['ilc']['RSTF'] == data_stf)
            else:
                irstf_ind = ([],)

            if 'iepdlc' in lcdata:
                iestf_ind = np.where(lcdata['iepdlc']['ESTF'] == data_stf)
            else:
                iestf_ind = ([],)

            if 'itfalc' in lcdata:
                itstf_ind = np.where(lcdata['itfalc']['TSTF'] == data_stf)
            else:
                itstf_ind = ([],)


            # CASCADE LC COLUMNS ACCORDING TO AVAILABILITY

            # CCD METADATA FROM RLC ONLY
            if len(rstf_ind[0]) > 0:
                xcc = np.asscalar(lcdata['rlc']['XCC'][rstf_ind])
            else:
                xcc = 'NaN'
            if len(rstf_ind[0]) > 0:
                ycc = np.asscalar(lcdata['rlc']['YCC'][rstf_ind])
            else:
                ycc = 'NaN'
            if len(rstf_ind[0]) > 0:
                bgv = np.asscalar(lcdata['rlc']['BGV'][rstf_ind])
            else:
                bgv = 'NaN'
            if len(rstf_ind[0]) > 0:
                bge = np.asscalar(lcdata['rlc']['BGE'][rstf_ind])
            else:
                bge = 'NaN'
            if len(rstf_ind[0]) > 0:
                fsv = np.asscalar(lcdata['rlc']['FSV'][rstf_ind])
            else:
                fsv = 'NaN'
            if len(rstf_ind[0]) > 0:
                fdv = np.asscalar(lcdata['rlc']['FDV'][rstf_ind])
            else:
                fdv = 'NaN'
            if len(rstf_ind[0]) > 0:
                fkv = np.asscalar(lcdata['rlc']['FKV'][rstf_ind])
            else:
                fkv = 'NaN'

            # INSTRUMENTAL MAGNITUDES FROM RLC > EPDLC > TFALC
            if len(rstf_ind[0]) > 0:
                im1 = np.asscalar(lcdata['rlc']['IM1'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                im1 = np.asscalar(lcdata['epdlc']['IM1'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                im1 = np.asscalar(lcdata['tfalc']['IM1'][tstf_ind])
            else:
                im1 = 'NaN'
            if len(rstf_ind[0]) > 0:
                ie1 = np.asscalar(lcdata['rlc']['IE1'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                ie1 = np.asscalar(lcdata['epdlc']['IE1'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                ie1 = np.asscalar(lcdata['tfalc']['IE1'][tstf_ind])
            else:
                ie1 = 'NaN'
            if len(rstf_ind[0]) > 0:
                iq1 = np.asscalar(lcdata['rlc']['IQ1'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                iq1 = np.asscalar(lcdata['epdlc']['IQ1'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                iq1 = np.asscalar(lcdata['tfalc']['IQ1'][tstf_ind])
            else:
                iq1 = 'NaN'

            if len(rstf_ind[0]) > 0:
                im2 = np.asscalar(lcdata['rlc']['IM2'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                im2 = np.asscalar(lcdata['epdlc']['IM2'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                im2 = np.asscalar(lcdata['tfalc']['IM2'][tstf_ind])
            else:
                im2 = 'NaN'
            if len(rstf_ind[0]) > 0:
                ie2 = np.asscalar(lcdata['rlc']['IE2'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                ie2 = np.asscalar(lcdata['epdlc']['IE2'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                ie2 = np.asscalar(lcdata['tfalc']['IE2'][tstf_ind])
            else:
                ie2 = 'NaN'
            if len(rstf_ind[0]) > 0:
                iq2 = np.asscalar(lcdata['rlc']['IQ2'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                iq2 = np.asscalar(lcdata['epdlc']['IQ2'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                iq2 = np.asscalar(lcdata['tfalc']['IQ2'][tstf_ind])
            else:
                iq2 = 'NaN'

            if len(rstf_ind[0]) > 0:
                im3 = np.asscalar(lcdata['rlc']['IM3'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                im3 = np.asscalar(lcdata['epdlc']['IM3'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                im3 = np.asscalar(lcdata['tfalc']['IM3'][tstf_ind])
            else:
                im3 = 'NaN'
            if len(rstf_ind[0]) > 0:
                ie3 = np.asscalar(lcdata['rlc']['IE3'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                ie3 = np.asscalar(lcdata['epdlc']['IE3'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                ie3 = np.asscalar(lcdata['tfalc']['IE3'][tstf_ind])
            else:
                ie3 = 'NaN'
            if len(rstf_ind[0]) > 0:
                iq3 = np.asscalar(lcdata['rlc']['IQ3'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                iq3 = np.asscalar(lcdata['epdlc']['IQ3'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                iq3 = np.asscalar(lcdata['tfalc']['IQ3'][tstf_ind])
            else:
                iq3 = 'NaN'

            # REDUCED MAGNITUDES FROM RLC > EPDLC > TFALC
            if len(rstf_ind[0]) > 0:
                rm1 = np.asscalar(lcdata['rlc']['RM1'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                rm1 = np.asscalar(lcdata['epdlc']['RM1'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                rm1 = np.asscalar(lcdata['tfalc']['RM1'][tstf_ind])
            else:
                rm1 = 'NaN'

            if len(rstf_ind[0]) > 0:
                rm2 = np.asscalar(lcdata['rlc']['RM2'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                rm2 = np.asscalar(lcdata['epdlc']['RM2'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                rm2 = np.asscalar(lcdata['tfalc']['RM2'][tstf_ind])
            else:
                rm2 = 'NaN'

            if len(rstf_ind[0]) > 0:
                rm3 = np.asscalar(lcdata['rlc']['RM3'][rstf_ind])
            elif len(estf_ind[0]) > 0:
                rm3 = np.asscalar(lcdata['epdlc']['RM3'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                rm3 = np.asscalar(lcdata['tfalc']['RM3'][tstf_ind])
            else:
                rm3 = 'NaN'

            # EPD MAGNITUDES FROM EPDLC > TFALC
            if len(estf_ind[0]) > 0:
                ep1 = np.asscalar(lcdata['epdlc']['EP1'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                ep1 = np.asscalar(lcdata['tfalc']['EP1'][tstf_ind])
            else:
                ep1 = 'NaN'

            if len(estf_ind[0]) > 0:
                ep2 = np.asscalar(lcdata['epdlc']['EP2'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                ep2 = np.asscalar(lcdata['tfalc']['EP2'][tstf_ind])
            else:
                ep2 = 'NaN'

            if len(estf_ind[0]) > 0:
                ep3 = np.asscalar(lcdata['epdlc']['EP3'][estf_ind])
            elif len(tstf_ind[0]) > 0:
                ep3 = np.asscalar(lcdata['tfalc']['EP3'][tstf_ind])
            else:
                ep3 = 'NaN'

            # TFA MAGNITUDES FROM TFALC ONLY
            if len(tstf_ind[0]) > 0:
                tf1 = np.asscalar(lcdata['tfalc']['TF1'][tstf_ind])
            else:
                tf1 = 'NaN'

            if len(tstf_ind[0]) > 0:
                tf2 = np.asscalar(lcdata['tfalc']['TF2'][tstf_ind])
            else:
                tf2 = 'NaN'

            if len(tstf_ind[0]) > 0:
                tf3 = np.asscalar(lcdata['tfalc']['TF3'][tstf_ind])
            else:
                tf3 = 'NaN'


            # ISM REDUCED MAGS FROM ILC > IEPDLC > ITFALC
            # ISM REDUCED MAG ERRS FROM ILC > IEPDLC > ITFALC
            # ISM REDUCED MAG FLAGS FROM ILC > IEPDLC > ITFALC

            if len(irstf_ind[0]) > 0:
                irm1 = np.asscalar(lcdata['ilc']['IRM1'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                irm1 = np.asscalar(lcdata['iepdlc']['IRM1'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                irm1 = np.asscalar(lcdata['itfalc']['IRM1'][itstf_ind])
            else:
                irm1 = 'NaN'
            if len(irstf_ind[0]) > 0:
                ire1 = np.asscalar(lcdata['ilc']['IRE1'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                ire1 = np.asscalar(lcdata['iepdlc']['IRE1'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                ire1 = np.asscalar(lcdata['itfalc']['IRE1'][itstf_ind])
            else:
                ire1 = 'NaN'
            if len(irstf_ind[0]) > 0:
                irq1 = np.asscalar(lcdata['ilc']['IRQ1'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                irq1 = np.asscalar(lcdata['iepdlc']['IRQ1'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                irq1 = np.asscalar(lcdata['itfalc']['IRQ1'][itstf_ind])
            else:
                irq1 = 'NaN'

            if len(irstf_ind[0]) > 0:
                irm2 = np.asscalar(lcdata['ilc']['IRM2'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                irm2 = np.asscalar(lcdata['iepdlc']['IRM2'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                irm2 = np.asscalar(lcdata['itfalc']['IRM2'][itstf_ind])
            else:
                irm2 = 'NaN'
            if len(irstf_ind[0]) > 0:
                ire2 = np.asscalar(lcdata['ilc']['IRE2'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                ire2 = np.asscalar(lcdata['iepdlc']['IRE2'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                ire2 = np.asscalar(lcdata['itfalc']['IRE2'][itstf_ind])
            else:
                ire2 = 'NaN'
            if len(irstf_ind[0]) > 0:
                irq2 = np.asscalar(lcdata['ilc']['IRQ2'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                irq2 = np.asscalar(lcdata['iepdlc']['IRQ2'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                irq2 = np.asscalar(lcdata['itfalc']['IRQ2'][itstf_ind])
            else:
                irq2 = 'NaN'

            if len(irstf_ind[0]) > 0:
                irm3 = np.asscalar(lcdata['ilc']['IRM3'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                irm3 = np.asscalar(lcdata['iepdlc']['IRM3'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                irm3 = np.asscalar(lcdata['itfalc']['IRM3'][itstf_ind])
            else:
                irm3 = 'NaN'
            if len(irstf_ind[0]) > 0:
                ire3 = np.asscalar(lcdata['ilc']['IRE3'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                ire3 = np.asscalar(lcdata['iepdlc']['IRE3'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                ire3 = np.asscalar(lcdata['itfalc']['IRE3'][itstf_ind])
            else:
                ire3 = 'NaN'
            if len(irstf_ind[0]) > 0:
                irq3 = np.asscalar(lcdata['ilc']['IRQ3'][irstf_ind])
            elif len(iestf_ind[0]) > 0:
                irq3 = np.asscalar(lcdata['iepdlc']['IRQ3'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                irq3 = np.asscalar(lcdata['itfalc']['IRQ3'][itstf_ind])
            else:
                irq3 = 'NaN'

            # ISM EPD MAGS FROM IEPDLC > ITFALC
            if len(iestf_ind[0]) > 0:
                iep1 = np.asscalar(lcdata['iepdlc']['IEP1'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                iep1 = np.asscalar(lcdata['itfalc']['IEP1'][itstf_ind])
            else:
                iep1 = 'NaN'

            if len(iestf_ind[0]) > 0:
                iep2 = np.asscalar(lcdata['iepdlc']['IEP2'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                iep2 = np.asscalar(lcdata['itfalc']['IEP2'][itstf_ind])
            else:
                iep2 = 'NaN'

            if len(iestf_ind[0]) > 0:
                iep3 = np.asscalar(lcdata['iepdlc']['IEP3'][iestf_ind])
            elif len(itstf_ind[0]) > 0:
                iep3 = np.asscalar(lcdata['itfalc']['IEP3'][itstf_ind])
            else:
                iep3 = 'NaN'

            # ISM TFA MAGS FROM ITFALC ONLY
            if len(itstf_ind[0]) > 0:
                itf1 = np.asscalar(lcdata['itfalc']['ITF1'][itstf_ind])
            else:
                itf1 = 'NaN'

            if len(itstf_ind[0]) > 0:
                itf2 = np.asscalar(lcdata['itfalc']['ITF2'][itstf_ind])
            else:
                itf2 = 'NaN'

            if len(itstf_ind[0]) > 0:
                itf3 = np.asscalar(lcdata['itfalc']['ITF3'][itstf_ind])
            else:
                itf3 = 'NaN'



            # the dataline containing ALL of the data for this object
            dataline = [rjds[data_ind],
                        'HN', # the HAT network collecting this data point
                        stid,
                        framenum,
                        ccds[data_ind],
                        filters[data_ind],
                        obsfield,
                        exp,
                        xcc,
                        ycc,
                        bgv,
                        bge,
                        fsv,
                        fdv,
                        fkv,
                        ha,
                        zd,
                        im1,ie1,iq1,
                        im2,ie2,iq2,
                        im3,ie3,iq3,
                        rm1,rm2,rm3,
                        ep1,ep2,ep3,
                        tf1,tf2,tf3,
                        irm1,ire1,irq1,
                        irm2,ire2,irq2,
                        irm3,ire3,irq3,
                        iep1,iep2,iep3,
                        itf1,itf2,itf3]

            datatable.append(dataline)

        # get rid of useless things to try to save memory
        del lcdata
        del stfs
        del rjds
        del filters
        del fields
        del ccds
        del hourangles
        del zenithdists

        # remove the collected LCs once we have no use for them
        if removecollected:
            for lckey in collected_lc_dict['hn']:
                os.remove(collected_lc_dict['hn'][lckey])


    ##########################################################
    # SECOND: collect all the available HATSouth data        #
    ##########################################################
    # !!!!WARNING!!!!! THIS IS NOW OUT OF DATE, NEEDS FIXING #
    ##########################################################
    # FIXME: add ISM phot here as well

    if collected_lc_dict['hs']:

        # get all of the information to be put into the datatable from all of
        # the lightcurves
        lc_cols = ['RJD','STF','CFN','CCD','FLD','ESTFC','TSTFC',
                   'XCC','YCC','BGV','BGE','FSV','FDV','FKV','IHA','IZD',
                   'IM1','IE1','IQ1','IM2','IE2','IQ2','IM3','IE3','IQ3',
                   'RM1','RM2','RM3',
                   'EP1','EP2','EP3',
                   'TF1','TF2','TF3']

        # connect to the HATSOUTH database
        db = mysql.connect(user=MYSQL_USER_HS,
                           passwd=MYSQL_PASS_HS,
                           db=MYSQL_DATA_HS,
                           host=MYSQL_HOST_HS)
        cur = db.cursor()

        column_sources = {}

        # this tortured expression consolidates the columns and sources
        for col in lc_cols:
            if (col in TEXTLC_CONSOLIDATION_SOURCES['hs']):
                for source in TEXTLC_CONSOLIDATION_SOURCES['hs'][col]:
                    if source in collected_lc_dict['hs']:
                        column_sources[col] = collected_lc_dict['hs'][source]
                        break
            else:
                if LOGGER:
                    LOGGER.warning('column %s is not a valid column,'
                                   ' skipping...' %
                                   col)
                if DEBUGMODE:
                    print('lcutils: column %s is not a valid column, '
                          'skipping...' %
                          col)

        lcdata = {}

        # load up all the LCs in collected_lc_dict
        for key in collected_lc_dict['hs']:

            lcf = open(collected_lc_dict['hs'][key],'r')
            lcdata[key] = lcf.readlines()
            lcf.close()

            # split, strip, transpose to get data in columns for fast retrieval
            lcdata[key] = [x.rstrip('\n') for x in lcdata[key]]
            lcdata[key] = [x.split() for x in lcdata[key]]
            lcdata[key] = zip(*lcdata[key])

        # this dict holds all of the requested data that gets turned into the
        # datatable later on
        allcols = {}

        # do a sanity check to make sure we have all the columns we asked for
        for col in lc_cols:

            if col not in column_sources:

                if LOGGER:
                    LOGGER.warning('column %s could not be retrieved, '
                                   'skipping...' %
                                   col)
                if DEBUGMODE:
                    print('lcutils: column %s could not be retrieved, '
                          'skipping...' %
                          col)

                continue

            # if we have this column, then prepare it for output
            else:

                # grab the column index from the file
                lcfpath = column_sources[col]
                lcftype = lcfpath.split('.')[-1]

                lcfcol_index = HATLC_COL_DEFS['hs'][lcftype].index(col)

                # grab the column data
                allcols[col] = np.array(lcdata[lcftype][lcfcol_index])


        # now we do the rest of the processing for HATSouth:
        # 1. form the RSTFC column from the rawlc's STF, CFN, and CCD columns
        # 2. use this column to match across the detections in the epd/tfalcs
        # 3. form the datalines with all of the data
        # 4. append our datalines to the end of the datatable

        if 'STF' in allcols and 'CFN' in allcols and 'CCD' in allcols:

            if DEBUGMODE:
                print('%s total unique detections across all LCs of type: %s' %
                      (len(rstfcs), ', '.join(collected_lc_dict['hs'].keys())))

            # get the RSFTC column for everything in the rawlc
            rstfcs = np.array(
                ['%s-%s_%s' % (x,y,z) for (x,y,z) in zip(allcols['STF'],
                                                         allcols['CFN'],
                                                         allcols['CCD'])]
                )

            rstfs = np.array(
                ['%s-%s' % (x,y) for (x,y) in zip(allcols['STF'],
                                                  allcols['CFN'])]
                )

            if DEBUGMODE:
                print('%s total unique detections across all LCs of type: %s' %
                      (len(rstfcs), ', '.join(collected_lc_dict['hs'].keys())))

            unique_stids = list(set(allcols['STF']))
            stids = np.array(allcols['STF'])
            framenums = np.array(allcols['CFN'])

            # query the images table for the values of:
            # filters
            # exptimes

            query_constraints = []
            query_conditions = []

            for st in unique_stids:

                st_frn_ind = np.where(stids == st)

                # we need to do this because MySQLdb doesn't map tuples directly
                # to SQL arrays as far as I know
                framenum_condition = '%s,' * len(st_frn_ind[0])
                framenum_condition = framenum_condition.strip(',')

                constraint = [st] + [int(x) for x in framenums[st_frn_ind]]
                query_constraints.extend(constraint)

                framenum_condition = 'IMfnum in (%s)' % framenum_condition
                full_condition = '(IMstid = %s and ' + framenum_condition + ')'
                query_conditions.append(full_condition)

            query_constraints = tuple(query_constraints)

            hs_query = ('select IMstid, IMfnum, IMccd, IMfilid, '
                        'IMtexp from Images where %s' %
                        ' or '.join(query_conditions))

            if LOGGER:
                LOGGER.debug('querying HATSOUTH DB for metadata...')
            if DEBUGMODE:
                print('lcutils: querying HATSOUTH DB for metadata...')

            cur.execute(hs_query, query_constraints)
            rows = cur.fetchall()

            cur.close()
            db.close()

            if DEBUGMODE:
                print('%s rows in HSCALIB DB '
                      'corresponding to these detections.' % len(rows))

            db_rstfcs = np.array(['%s-%s_%s' % (x[0],x[1],x[2]) for x in rows])
            filters = np.array([x[3] for x in rows])
            exptimes = np.array([x[4] for x in rows])

            # get the RJDs for everything in the rawlc
            rjds = allcols['RJD']

            # sort by rjd
            rjd_sort_ind = np.argsort(rjds)
            rstfcs = rstfcs[rjd_sort_ind]
            rstfs = rstfs[rjd_sort_ind]
            rjds = rjds[rjd_sort_ind]

            if LOGGER:
                LOGGER.info('consolidating HS lightcurves...')
            if DEBUGMODE:
                print('lcutils: consolidating HS lightcurves...')

            for data_ind, data_rstfc in enumerate(rstfcs):

                # get the indices for reduced, EPD, and TFA mags while handling
                # missing values

                rstfc_ind = np.where(rstfcs == data_rstfc)

                if 'ESTFC' in allcols:
                    estfc_ind = np.where(allcols['ESTFC'] == data_rstfc)
                else:
                    estfc_ind = ([],)

                if 'TSTFC' in allcols:
                    tstfc_ind = np.where(allcols['TSTFC'] == data_rstfc)
                else:
                    tstfc_ind = ([],)

                hsdb_ind = np.where(db_rstfcs == data_rstfc)
                filt = (np.asscalar(filters[hsdb_ind])
                        if len(hsdb_ind[0]) > 0 else 'NaN')

                # need to turn the float into a str to check for busted floats
                expt = ((TEXTLC_OUTPUT_COLUMNS['EXP'][1] %
                         np.asscalar(exptimes[hsdb_ind]))
                        if len(hsdb_ind[0]) > 0 else 'NaN')

                # the dataline containing ALL of the data for this object
                dataline = [rjds[data_ind],
                            'HS', # the HAT network collecting this data point
                            (np.asscalar(allcols['STF'][data_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['CFN'][data_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['CCD'][data_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            filt,
                            (np.asscalar(allcols['FLD'][data_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            expt,
                            (np.asscalar(allcols['XCC'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['YCC'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['BGV'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['BGE'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['FSV'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['FDV'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['FKV'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IHA'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IZD'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IM1'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IE1'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IQ1'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IM2'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IE2'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IQ2'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IM3'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IE3'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['IQ3'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['RM1'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['RM2'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['RM3'][rstfc_ind])
                             if len(rstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['EP1'][estfc_ind])
                             if len(estfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['EP2'][estfc_ind])
                             if len(estfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['EP3'][estfc_ind])
                             if len(estfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['TF1'][tstfc_ind])
                             if len(tstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['TF2'][tstfc_ind])
                             if len(tstfc_ind[0]) > 0 else 'NaN'),
                            (np.asscalar(allcols['TF3'][tstfc_ind])
                             if len(tstfc_ind[0]) > 0 else 'NaN'),
                            'NaN','NaN','NaN', # ISM reduced mags aper 1
                            'NaN','NaN','NaN', # ISM reduced mags aper 2
                            'NaN','NaN','NaN', # ISM reduced mags aper 3
                            'NaN','NaN','NaN', # ISM EPD mags
                            'NaN','NaN','NaN'] # ISM TFA mags

                datatable.append(dataline)

        else:

            if LOGGER:
                LOGGER.error('STF, CFN, or CCD column missing in collected '
                             'LCs for %s' % hatid)
            if DEBUGMODE:
                print('STF, CFN, or CCD column missing in collected '
                      'LCs for %s' % hatid)



        del lcdata
        del allcols
        del rstfcs
        del rstfs
        del rjds

        # remove the collected LCs once we have no use for them
        if removecollected:
            for lckey in collected_lc_dict['hs']:
                os.remove(collected_lc_dict['hs'][lckey])

    ##############################################################
    # THIRD: resort the full datatable if we have HN AND HS data #
    ##############################################################
    # FIXME: fix this once HS ISM phot is added

    if collected_lc_dict['hn'] and collected_lc_dict['hs']:

        datatable_sortcol = np.array([x[0] for x in datatable])
        datatable = np.array(datatable)

        datatable_sortind = np.argsort(datatable_sortcol)

        # NOTE: we leave this as an ndarray in the case of HN+HS data
        datatable = datatable[datatable_sortind]

        del datatable_sortcol
        del datatable_sortind

    ################################################################
    # FINALLY: return the full datatable and clean up if necessary #
    ################################################################

    if removecollected:

        collectedlc_info_fname = '%s-collected-lcs.pkl' % hatid
        collectedlc_info_fpath = os.path.join(LCCACHE,
                                              collectedlc_info_fname)
        if os.path.exists(collectedlc_info_fpath):
            os.remove(collectedlc_info_fpath)


    return datatable



def process_consolidated_lightcurve(
    hatid,
    datatable,
    requestedcols,
    database=None,
    ):
    '''
    This does the final processing of the data table returned by
    consolidate_all_lightcurves above.

    requestedcols = a list of col names to pull out of the datatable and write
                    to the output, must all be in TEXTLC_OUTPUT_COLUMNS

    '''

    # make sure all of the requested columns are legit
    reqcols = requestedcols

    if not all([x in TEXTLC_OUTPUT_COLUMNS for x in reqcols]):

        if LOGGER:
            LOGGER.error('some requested columns cannot be retrieved')
        if DEBUGMODE:
            print('lcutils: some requested columns cannot be retrieved')

        return None


    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    # grab the STF column and get the station ids
    stf = [int(x[2]) for x in datatable]
    networks = [x[1] for x in datatable]

    all_stations = ['%s%02i' % (x,y) for (x,y) in zip(networks, stf)]
    unique_stations = list(set(all_stations))


    if 'HJD' in reqcols or 'BJD' in reqcols:

        # get the station locations to convert to BJD
        obslat, obslon, obsalt = ([HAT_LOCATIONS[x][y][0] for
                                   (x,y) in zip(networks,stf)],
                                  [HAT_LOCATIONS[x][y][1] for
                                   (x,y) in zip(networks,stf)],
                                  [HAT_LOCATIONS[x][y][2] for
                                   (x,y) in zip(networks,stf)])


    # grab info for this object from the 2MASS table
    objinfo = searchutils.twomass_info(
        hatid,
        cols='hat_id,ra,decl,vmag,rmag,imag,jmag,hmag,kmag,twomass_id',
        database=database,
        pprint=None
        )

    # if this object doesn't exist in the HAT 2MASS table (horribly unlikely),
    # then return immediately with None because we're not sure what went wrong
    if len(objinfo) == 0:
        return None
    else:
        objra, objdec, objv, objr, obji, objj, objh, objk, tmid = objinfo[0][1:]

    # put together a dictionary of the final output values
    outlc = {}

    # put the list of HAT stations into the dict
    outlc['hatstations'] = unique_stations

    # put the list of unique filters into the dict
    allfilters = [('%s-%s' % (x[1], x[5])) for x in datatable]
    uniquefilters = list(set(allfilters))

    filterdescs = []

    for x in uniquefilters:

        filtnet, filtid = x.split('-')
        filtid = int(filtid) if filtid.isdigit() else 0
        if filtnet == 'HN':
            filter_description = ('%s - %s - %s' %
                                  (filtid,
                                   HATNET_FILTER_DESCRIPTIONS[filtid][0],
                                   HATNET_FILTER_DESCRIPTIONS[filtid][1]))
            filterdescs.append(filter_description)
        elif filtnet == 'HS':
            filter_description = ('%s - %s - %s' %
                                  (filtid,
                                   HATSOUTH_FILTER_DESCRIPTIONS[filtid][0],
                                   HATSOUTH_FILTER_DESCRIPTIONS[filtid][1]))
            filterdescs.append(filter_description)

    outlc['filters'] = filterdescs


    for col in reqcols:

        if DEBUGMODE:
            print(col)

        if col == 'RSTF':

            outlc['RSTF'] = ['%s-%s' % (x[2],x[3]) for x in datatable]

        elif col == 'STF':

            outlc['STF'] = stf

        elif col == 'CFN':

            outlc['CFN'] = [int(x[3]) for x in datatable]

        # deal with a BJD column
        elif col == 'BJD':

            rjd = [x[0] for x in datatable]
            rjd = np.array([float(x) for x in rjd])

            # convert to full JD
            rjd = rjd + 2400000.0

            # convert to BJD
            bjd = [
                jd_to_bjd(x,objra,objdec,y,z,w)
                for (x,y,z,w) in zip(rjd,obslat,obslon,obsalt)
                ]

            outlc['BJD'] = bjd

        # deal with a HJD column
        elif col == 'HJD':

            rjd = [x[0] for x in datatable]
            rjd = np.array([float(x) for x in rjd])

            # convert to full JD
            rjd = rjd + 2400000.0

            # convert to BJD
            hjd = [jd_to_hjd(x,objra,objdec) for x in rjd]

            outlc['HJD'] = hjd

        # deal with a MJD column
        elif col == 'MJD':

            rjd = [x[0] for x in datatable]
            rjd = np.array([float(x) for x in rjd])

            # convert to MJD
            mjd = rjd - 0.5
            outlc['MJD'] = mjd.tolist()

        # deal with a RJD column
        elif col == 'RJD':

            outlc['RJD'] = [float(x[0]) for x in datatable]

        # deal with a FJD column
        elif col == 'FJD':

            rjd = [x[0] for x in datatable]
            rjd = np.array([float(x) for x in rjd])

            fjd = rjd + 2400000.0

            outlc['FJD'] = fjd.tolist()

        elif col in ['IQ1','IQ2','IQ3',
                     'FLD','CCD','NET',
                     'IRQ1','IRQ2','IRQ3']:

            coldataind = TEXTLC_DATATABLE_COLUMNS.index(col)
            outlc[col] = [x[coldataind] for x in datatable]

        elif col == 'FLT':

            coldataind = TEXTLC_DATATABLE_COLUMNS.index(col)
            outlc[col] = [(int(x[coldataind]) if x[coldataind] is not None else 0)
                          for x in datatable]

        # now process the rest of the columns (these don't need special
        # treatment)
        else:

            coldataind = TEXTLC_DATATABLE_COLUMNS.index(col)

            # this is to make sure all elements in this column can be turned
            # into floats, if any can't, we'll turn them into NaNs
            coldata = [x[coldataind] for x in datatable]

            floatcheck = [(len(x.split('.')) == 2) if x else False
                          for x in coldata]
            coldatachecked = [float(x) if y else float('NaN')
                              for (x,y) in zip(coldata,floatcheck)]
            outlc[col] = coldatachecked


    # add the column header list to the dict so we know what order to put the
    # columns in for the output LC
    outlc['cols'] = reqcols

    # this is the master output LC dictionary
    outlcdict = lcform.OUTPUTLC_FORMATTERS['lcdict'](outlc, objinfo[0])

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    # delete useless things
    del allfilters
    del coldata

    return outlcdict



def check_hatlc_database(hatids,
                         database=None):
    '''
    This checks to see if the objects in list hatids are present in the
    hat_lightcurves table of the LC DB. Returns the full_lc_fpath and
    last_updated columns if present, otherwise returns None.

    database is an instance of the LCDB object from utils/lcdb.py. If not passed
    as input, will create its own database connection using the credentials in
    conf/lcserver.conf.

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

    sql_statement = ('select hat_id, full_lc_fpath, last_updated, '
                     'access_groups '
                     'from hat_lightcurves '
                     'where hat_id in %s')
    cur.execute(sql_statement, (tuple(hatids),))
    rows = cur.fetchall()

    if len(rows) > 0:
        results = rows
    else:
        results = None

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return results



def update_hatlc_database(lcdict,
                          full_lc_fpath,
                          lc_coord_errbox=0.5,
                          lc_accessgroup='superuser,hatgroup',
                          nolc_hatid=None,
                          database=None):
    '''
    This updates the hat_lightcurves table in the LC DB. Returns True if update
    was successful, otherwise returns False.

    lcdict -> from lcutils_hatlc functions

    full_lc_fpath -> path to the full LC on disk

    lc_coord_errbox -> used to construct the coordinates errbox (in arcsec)

    lc_is_public -> sets the made_public column for this hatid

    database is an instance of the LCDB object from utils/lcdb.py. If not passed
    as input, will create its own database connection using the credentials in
    conf/lcserver.conf.

    '''

    closedb = False
    success = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    # if we have a valid lcdict and a valid full_lc_fpath, go ahead and grab
    # this object's info from the lcdict
    if lcdict is not None and full_lc_fpath is not None and full_lc_fpath != '':

        hatid, ra, dec, vmag, rmag, imag, jmag, hmag, kmag, tmid = (
            lcdict['hatid'],
            lcdict['ra'],
            lcdict['dec'],
            lcdict['mags'][0],
            lcdict['mags'][1],
            lcdict['mags'][2],
            lcdict['mags'][3],
            lcdict['mags'][4],
            lcdict['mags'][5],
            lcdict['twomassid']
            )

        hat_ndet = lcdict['ndet']
        hat_stations = ','.join(lcdict['hatstations'])

        hatfield, hatfieldobjid = hatid.split('-')[1:]

        errbox = (ra + lc_coord_errbox/(3600.0*np.cos(np.radians(dec))),
                  dec + lc_coord_errbox/3600.0,
                  ra - lc_coord_errbox/(3600.0*np.cos(np.radians(dec))),
                  dec - lc_coord_errbox/3600.0)

        sql_params = (hatid, int(hatfield), int(hatfieldobjid),
                      tmid, ra, dec,
                      errbox[0], errbox[1], errbox[2], errbox[3],
                      vmag, rmag, imag, jmag, hmag, kmag,
                      hat_ndet, hat_stations, lc_accessgroup, full_lc_fpath)

        sql_statement = ('insert into hat_lightcurves '
                         '(hat_id, hat_field, hat_field_objid, '
                         'twomass_id, ra, decl, '
                         'errbox, '
                         'vmag, rmag, imag, jmag, hmag, kmag, '
                         'ndet, hat_stations, access_groups, full_lc_fpath) '
                         'values '
                         '(%s, %s, %s, '
                         ' %s, %s, %s, '
                         'box(point(%s,%s),point(%s,%s)), '
                         '%s, %s, %s, %s, %s, %s, '
                         '%s, %s, %s, %s)')

    # if the lcdict is None and full_lc_fpath is either None or '', then the LC
    # for this object does not exist in our archives yet. update the database
    # with this fact
    else:

        # get this object's information using the searchutils.twomass_info
        # function, since we don't have it from the lcdict
        objinfo = searchutils.twomass_info(
            nolc_hatid,
            cols='hat_id,ra,decl,vmag,rmag,imag,jmag,hmag,kmag,twomass_id',
            database=database,
            pprint=None
            )

        if objinfo and len(objinfo) > 0:
            hatid, ra, dec, vmag, rmag, imag, jmag, hmag, kmag, tmid = objinfo[0]
        else:
            return False

        hat_ndet = 0
        hat_stations = None

        hatfield, hatfieldobjid = hatid.split('-')[1:]

        errbox = (ra + lc_coord_errbox/(3600.0*np.cos(np.radians(dec))),
                  dec + lc_coord_errbox/3600.0,
                  ra - lc_coord_errbox/(3600.0*np.cos(np.radians(dec))),
                  dec - lc_coord_errbox/3600.0)


        sql_params = (hatid, int(hatfield), int(hatfieldobjid),
                      tmid, ra, dec,
                      errbox[0], errbox[1], errbox[2], errbox[3],
                      vmag, rmag, imag, jmag, hmag, kmag,
                      hat_ndet, hat_stations, lc_accessgroup, None)

        sql_statement = ('insert into hat_lightcurves '
                         '(hat_id, hat_field, hat_field_objid, '
                         'twomass_id, ra, decl, '
                         'errbox, '
                         'vmag, rmag, imag, jmag, hmag, kmag, '
                         'ndet, hat_stations, access_groups, full_lc_fpath) '
                         'values '
                         '(%s, %s, %s, '
                         '%s, %s, %s, '
                         'box(point(%s,%s),point(%s,%s)), '
                         '%s, %s, %s, %s, %s, %s, '
                         '%s, %s, %s, %s)')


    ##
    # now do the database work below
    ##

    try:
        cur.execute(sql_statement, sql_params)
        database.commit()
        success = True

    # if this hatid already exists in the database, then update its records
    except pg.IntegrityError as e:

        if LOGGER:
            LOGGER.warning('%s already exists, updating...' % (hatid,))
        if DEBUGMODE:
            print('%s already exists, updating...' % (hatid,))

        database.rollback()
        sql_statement = (
            'update hat_lightcurves set '
            'ndet = %s, hat_stations = %s, '
            'access_groups = %s, full_lc_fpath = %s, '
            'last_updated = current_timestamp '
            'where hat_id = %s'
            )
        sql_params = (hat_ndet,
                      hat_stations,
                      lc_accessgroup,
                      full_lc_fpath,
                      hatid)
        cur.execute(sql_statement, sql_params)
        database.commit()
        success = True

    except pg.InterfaceError as e:

        if LOGGER:
            LOGGER.error('inserting %s into LC DB failed. error was: %s'
                         % (hatid, e))
        if DEBUGMODE:
            print('inserting %s into LC DB failed. error was: %s'
                  % (hatid, e))
        database.rollback()

    except pg.DatabaseError as e:
        if LOGGER:
            LOGGER.error('inserting %s into LC DB failed. error was: %s' %
                         (hatid,e))
        if DEBUGMODE:
            print('inserting %s into LC DB failed. error was: %s' % (hatid,e))
        database.rollback()

    finally:
        if not success:
            database.rollback()

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return success


def update_database_lcstatus(hatid,
                             set_fetched=None,
                             set_collected=None,
                             set_processed=None,
                             database=None):
    '''
    This updates the status flags for a HATID in the database.

    '''

    closedb = False
    success = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    flag_cols = []
    flag_vals = []

    if set_fetched is not None:

        flag_cols.append('fetched = %s')
        flag_vals.append(set_fetched)

    if set_collected is not None:

        flag_cols.append('collected = %s')
        flag_vals.append(set_collected)

    if set_processed is not None:

        flag_cols.append('processed = %s')
        flag_vals.append(set_processed)


    if (set_fetched is not None or
        set_collected is not None or
        set_processed is not None):

        query = ('update hat_lightcurves set {flag_cols} where hat_id = %s')
        query = query.format(flag_cols=', '.join(flag_cols))
        query_params = tuple(flag_vals + [hatid])

        try:

            cur.execute(query, query_params)
            database.commit()
            success = True

        except Exception as e:
            LOGGER.error('could not update LC status for %s, error was: %s'
                         % (e,hatid))
            if DEBUGMODE:
                print('could not update LC status for %s, error was: %s'
                       % (e,hatid))
            success = False

    else:
        success = False


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return success


def check_hatlc_filesystem(hatid,
                           columns,
                           filters,
                           compress,
                           outtype,
                           normalized=False,
                           lcdir=None):
    '''
    This just checks the filesystem to see if there's an existing file with the
    requested columns, filters, compression type, and output type for hatid. If
    so, returns the path to it. Also checks if an LC file with all columns in
    the canonical format (csv.gz) also exists, and returns its path. If neither
    file exists, returns (False, False)

    '''

    if filters:
        lcfilters_string = [list(x) for x in filters]
        lcfilters_string = repr(lcfilters_string)
        lcfilterhash = hashlib.md5(lcfilters_string).hexdigest()
    else:
        lcfilterhash = None

    requested_lcfname_hash = lcform.lcfname_hash(hatid,
                                                 columns,
                                                 outtype,
                                                 filterhash=lcfilterhash,
                                                 normalized=normalized,
                                                 compressed=compress)

    full_lcfname_hash = lcform.lcfname_hash(hatid,
                                            HATLC_OUTPUT_COLUMNS['full'],
                                            'csv',
                                            filterhash=None,
                                            normalized=False,
                                            compressed=compress)

    if lcdir:
        lcsearchdir = lcdir
    else:
        hatfield = hatid.split('-')[1]
        lcsearchdir = os.path.join(LCCACHE, hatfield)

    lcsearch_paths = [os.path.join(lcsearchdir, requested_lcfname_hash),
                      os.path.join(lcsearchdir, full_lcfname_hash)]

    lcf_exists = [os.path.exists(x) for x in lcsearch_paths]

    lcfile_paths = [(os.path.abspath(x) if y else None)
                    for (x,y) in zip(lcsearch_paths, lcf_exists)]

    return lcfile_paths
