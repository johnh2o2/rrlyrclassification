#!/usr/bin/env python

'''
lcutils_formatters.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Dec 2013

Contains various useful tools for formatting LC files.

'''

import logging
import os.path
import subprocess
import shlex
import hashlib
from json import dumps

from math import isnan

import numpy as np
import pyfits

###################
## LOCAL IMPORTS ##
###################

import lcutils_config as conf

#############
## LOGGING ##
#############

# setup a logger
LOGGER = logging.getLogger('lcutils_formatters')
LOGGER.addHandler(logging.NullHandler())

# default to debug mode = False
DEBUGMODE = False

def set_debug(debugbool):
    globals()['DEBUGMODE'] = debugbool


##################################
## OUTPUT LIGHTCURVE FORMATTERS ##
##################################

def lcfname_hash(hatid,
                 cols,
                 form,
                 filterhash=None,
                 normalized=False,
                 compressed='gz'):
    '''
    This generates a filename encoding the LC's HATID, columns, and format.

    '''

    colstr = ','.join(cols)

    colhash = hashlib.md5(colstr).hexdigest()

    basename = '{hatid}-{colhash}{filterhash}{normflag}{form}{compression}'

    if filterhash:
        lcfilterhash = '-filtered-%s' % filterhash
    else:
        lcfilterhash = ''

    if normalized:
        lcnormflag = '-normalized'
    else:
        lcnormflag = ''

    if compressed:
        lccompression = '.%s' % compressed
    else:
        lccompression = ''

    if form == 'ssv':
        lcform = '.hatlc'
    elif form == 'csv':
        lcform = '-hatlc.csv'
    elif form == 'fits':
        lcform = '-hatlc.fits'
    elif form == 'json':
        lcform = '-hatlc.json'
    else:
        lcform = '-hatlc.csv'

    # now combine everything into the filename
    outfname = basename.format(hatid=hatid,
                               colhash=colhash,
                               filterhash=lcfilterhash,
                               normflag=lcnormflag,
                               form=lcform,
                               compression=lccompression)

    return outfname


def compress_lcfile(lcfilepath, method='gz'):
    '''
    This compresses the lightcurve using the method specified.

    method = 'gz', 'bz2'

    '''

    zip_cmds = {'gz':'gzip -f %s' % lcfilepath,
                'bz2':'bzip2 -f %s' % lcfilepath}

    outfile = lcfilepath + '.' + method

    proccmd = shlex.split(zip_cmds[method])

    try:
        retcode = subprocess.check_call(proccmd)
        return outfile
    except subprocess.CalledProcessError:
        return None


def decompress_lcfile(lcfilepath):
    '''
    This decompresses the lightcurve.

    '''

    unzip_cmds = {'gz':'gunzip %s' % lcfilepath,
                  'bz2':'bunzip2 %s' % lcfilepath}

    if '.bz2' in lcfilepath:
        outfile = lcfilepath.rstrip('.bz2')
        method = 'bz2'
    elif '.gz' in lcfilepath:
        outfile = lcfilepath.rstrip('.gz')
        method = 'gz'

    proccmd = shlex.split(unzip_cmds[method])

    try:
        retcode = subprocess.check_call(proccmd)
        return outfile
    except subprocess.CalledProcessError:
        return None


def lcdict_to_text(lcdict,
                   objinfo,
                   form='csv',
                   outputdir=None,
                   compress='gz',
                   filterhash=None,
                   normalized=False):
    '''
    This returns a text table file from the light-curve dictionary.

    form = 'ssv' -> space-separated values
           'csv' -> comma-separated values

    output .hatlc LCs have the following filename schema:

    <HAT ID>.hatlc

    The header of the LC looks like the following:

    # HAT-XXX-YYYYYYY
    # 2MASS JXXXXXX.XX+/-YYYYYY.Y
    # RA = XXX.XXX deg, DEC = XX.XXX deg
    # V = XX.YY, R = XX.YY, J = XX.YY, H = XX.YY, K = XX.YY
    #
    # total LC points: XXXXX
    # HAT stations: X, Y, Z, etc.
    #
    # col XX - <col header> - <col description>
    # col YY - <col header> - <col description>
    # col ZZ - <col header> - <col description>
    # ...,
    # etc.



    '''

    colheaders = lcdict['cols']

    hatid, ra, dec, vmag, rmag, imag, jmag, hmag, kmag, tmid = objinfo
    ndet = len(lcdict[colheaders[0]])
    hatstations = sorted(lcdict['hatstations'])
    hatstations = ', '.join([x for x in hatstations])

    # format the column list
    col_list = []

    for i, col in enumerate(colheaders):
        col_line = '# col %02i - %s - %s' % (i+1,
                                             col,
                                             conf.TEXTLC_OUTPUT_COLUMNS[col][0])
        col_list.append(col_line)

    col_list = '\n'.join(col_list)

    filterlist = ['# %s' % x for x in lcdict['filters']]
    filterlist = '\n'.join(filterlist)

    # format the LC header
    lc_header = conf.TEXTLC_HEADER_TEMPLATE.format(
        hatid=hatid,
        twomassid=tmid.strip(),
        ra=ra or -999.0,
        dec=dec or -999.0,
        vmag=vmag or -999.0,
        rmag=rmag or -999.0,
        imag=imag or -999.0,
        jmag=jmag or -999.0,
        hmag=hmag or -999.0,
        kmag=kmag or -999.0,
        tmid=tmid,
        ndet=ndet,
        hatstations=hatstations,
        filterlist=filterlist,
        columnlist=col_list
        )

    # collect the lightcurve data for writing to the output file
    lc_lines = [lcdict[x] for x in colheaders]
    lc_lines = zip(*lc_lines)

    # figure out the format for each line
    lc_formatted_line = [conf.TEXTLC_OUTPUT_COLUMNS[x][1] for x in colheaders]

    if form == 'ssv':
        lc_formatted_line = ' '.join(lc_formatted_line)
    elif form == 'csv':
        lc_formatted_line = ','.join(lc_formatted_line)

    outfname = lcfname_hash(hatid,
                            lcdict['cols'],
                            form,
                            filterhash=filterhash,
                            normalized=normalized,
                            compressed=None) # don't use compress here since we
                                             # do it below

    if outputdir:
        outfname = os.path.join(outputdir, outfname)
    else:
        outfname = os.path.join(conf.LCCACHE, outfname)

    # open and write the header to the output file
    outf = open(outfname,'w')
    outf.write(lc_header)

    # now write the rest of the output LC
    for lineind, line in enumerate(lc_lines):

        try:

            outline = lc_formatted_line % line
            outline = outline.strip()

            if form == 'csv':
                outline = outline.replace(' ','')

            outf.write('%s\n' % outline)

        except:
            print('HATID %s: formatting failed for '
                  'line %s: %s, line skipped...' %
                  (hatid, lineind, line))
            continue

    outf.close()

    if compress:
        outfname = compress_lcfile(outfname,
                                   method=compress)

    return outfname



def lcdict_to_csv(lcdict,
                  objinfo,
                  outputdir=None,
                  compress='gz',
                  filterhash=None,
                  normalized=False):
    '''
    This returns a CSV text table from the light-curve dictionary.

    '''

    return lcdict_to_text(lcdict,
                          objinfo,
                          form='csv',
                          outputdir=outputdir,
                          compress=compress,
                          filterhash=filterhash,
                          normalized=normalized)



def lcdict_to_fits(lcdict,
                   objinfo,
                   outputdir=None,
                   compress='gz',
                   filterhash=None,
                   normalized=False):
    '''
    This returns a FITS file from the light-curve dictionary. The FITS file is
    written to the conf.LCCACHE directory by default.

    '''

    colheaders = lcdict['cols']

    hatid, ra, dec, vmag, rmag, imag, jmag, hmag, kmag, tmid = objinfo
    ndet = len(lcdict[colheaders[0]])
    hatstations = sorted(lcdict['hatstations'])
    hatstations = ', '.join([x for x in hatstations])
    filterlist = ', '.join(lcdict['filters'])

    # this is stuff for the primary HDU of the FITS header
    fitsheader_list = [('hatid', hatid, 'HAT ID of the object'),
                       ('2massid', tmid, '2MASS ID of the object'),
                       ('ra', ra, 'Right ascension (J2000) [deg]'),
                       ('dec', dec, 'Declination (J2000) [deg]'),
                       ('vmag', vmag, 'V-band magnitude'),
                       ('rmag', rmag, 'R-band magnitude'),
                       ('imag', imag, 'I-band magnitude'),
                       ('jmag', jmag, 'J-band magnitude'),
                       ('hmag', hmag, 'H-band magnitude'),
                       ('kmag', jmag, 'K-band magnitude'),
                       ('ndet', ndet, 'Number of detections in LC'),
                       ('hats', hatstations, ('HAT stations contributing '
                                            'data for LC')),
                       ('filters', filterlist, 'Filters used for LC data'),]

    # the output name of the FITS file, a .gz will be added at the end if we
    # decide to compress it
    outfname = lcfname_hash(hatid,
                            lcdict['cols'],
                            'fits',
                            filterhash=filterhash,
                            normalized=normalized,
                            compressed=None) # don't use compress here since we
                                             # do it below

    if outputdir:
        outfname = os.path.join(outputdir, outfname)
    else:
        outfname = os.path.join(conf.LCCACHE, outfname)

    # define the table columns
    fitstable_columns = []

    for col in colheaders:

        data = lcdict[col]
        dataformat = conf.TEXTLC_OUTPUT_COLUMNS[col][2]
        dataname = col

        fitscol = pyfits.Column(name=dataname,
                                format=dataformat,
                                array=np.array(data))
        fitstable_columns.append(fitscol)

    # generate the table HDU
    fitstable_hdu = pyfits.new_table(fitstable_columns)

    # generate the primary HDU and add our object's info to its header
    primaryhdu = pyfits.PrimaryHDU()
    for card in fitsheader_list:
        primaryhdu.header[card[0]] = (card[1], card[2])

    hdulist = pyfits.HDUList([primaryhdu, fitstable_hdu])
    hdulist.writeto(outfname)

    if compress:
        outfname = compress_lcfile(outfname,
                                   method=compress)

    return outfname


def lcdict_to_hatlcdict(lcdict,
                        objinfo):
    '''
    This just adds the objinfo keys and values to the dictionary lcdict and
    returns it.

    '''
    colheaders = lcdict['cols']

    hatid, ra, dec, vmag, rmag, imag, jmag, hmag, kmag, tmid = objinfo
    ndet = len(lcdict[colheaders[0]])
    hatstations = sorted(lcdict['hatstations'])
    hatstations = ', '.join([x for x in hatstations])
    filterlist = ', '.join(lcdict['filters'])

    lcdict['hatid'] = hatid
    lcdict['ra'] = ra
    lcdict['dec'] = dec
    lcdict['mags'] = [vmag, rmag, imag, jmag, hmag, kmag]
    lcdict['twomassid'] = tmid
    lcdict['ndet'] = ndet

    # convert all the lists in the lcdict to ndarrays to keep the lcdict the
    # same format as the output of read_consolidated_hatlc
    for col in colheaders:
        lcdict[col] = np.array(lcdict[col])

    return lcdict


def lcdict_to_json(lcdict,
                   objinfo):
    '''
    This converts an lcdict to JSON format (while also converting to
    consolidated hatlc format).

    '''

    colheaders = lcdict['cols']

    hatid, ra, dec, vmag, rmag, imag, jmag, hmag, kmag, tmid = objinfo
    ndet = len(lcdict[colheaders[0]])
    hatstations = sorted(lcdict['hatstations'])
    hatstations = ', '.join([x for x in hatstations])
    filterlist = ', '.join(lcdict['filters'])

    lcdict['hatid'] = hatid
    lcdict['ra'] = ra
    lcdict['dec'] = dec
    lcdict['mags'] = [vmag, rmag, imag, jmag, hmag, kmag]
    lcdict['twomassid'] = tmid
    lcdict['ndet'] = ndet

    # dump the columns to lists
    for col in colheaders:
        lcdict[col] = lcdict[col].tolist()
        # get rid of the NaNs
        lcdict[col] = [x if not isnan(x) else None for x in lcdict[col]]

    return dumps(lcdict, ensure_ascii=True)


def hatlc_dict_to_json(hatlcdict):
    '''
    This converts an consolidated hatlc dict to JSON format.

    '''

    colheaders = hatlcdict['cols']

    for col in colheaders:
        hatlcdict[col] = hatlcdict[col].tolist()

    return dumps(lcdict, ensure_ascii=True)


OUTPUTLC_FORMATTERS = {'ssv':lcdict_to_text,
                       'csv':lcdict_to_csv,
                       'fits':lcdict_to_fits,
                       'lcdict':lcdict_to_hatlcdict,
                       'json':lcdict_to_json}
