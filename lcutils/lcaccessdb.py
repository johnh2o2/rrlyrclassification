#!/usr/bin/env python

'''
lcaccessdb.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Jan 2014

This provides a DB interface for checking and administrating access to HAT
lightcurves and HAT fields for specific users and apikeys.

'''

import logging
import lcdb

from lcauthdb import check_user, get_apikey_username

#############
## LOGGING ##
#############

# setup a logger
LOGGER = logging.getLogger('lcaccessdb')
LOGGER.addHandler(logging.NullHandler())

# default to debug mode = False
DEBUGMODE = False

def set_debug(debugbool):
    globals()['DEBUGMODE'] = debugbool


##############################
## AUTHENTICATION FUNCTIONS ##
##############################

def auth_hatlcs(hat_list,
                username=None,
                apikey=None,
                existing_lcs_only=True,
                database=None):
    '''
    This checks if either the username or the apikey provided is authorized to
    access the lightcurves for the objects in hat_list. hat_list is a list of
    HATIDs.

    Returns a list of the following form:

    ['hatid1', 'hatid2', 'hatid3', ...]

    if these hatids from the provided hat_list are authorized for this user. If
    none of the provided hatids in hat_list are authorized, returns a tuple
    (False, None, 'no_hatlcs_authorized')

    '''

    # set up a database connection
    closedb = False

    # if both username and apikey are provided, then return a failure
    if username and apikey:

        return (False, None, None, 'apikey_and_username_provided')

    # we need either a username or an apikey, but not both. if neither is
    # provided, we use the anonuser to try access
    elif ((username and not apikey) or
          (apikey and not username) or
          (not apikey and not username)):

        if database:
            cur, handle = database.newcursor()
        else:
            database = lcdb.LCDB()
            database.open_default()
            cur, handle = database.newcursor()
            closedb = True

        # figure out the user/apikey info to later compare the accessgroup with
        # the hatlc accessgroups. if neither apikey nor username are provided,
        # use the anonuser
        if not apikey and not username:
            user_info = check_user('anonuser',
                                   getinfo=True,
                                   database=database)
            access_user = 'anonuser'

        # otherwise, if one of them is provided, figure out the corresponding
        # user info
        else:

            if apikey and not username:
                access_user = get_apikey_username(apikey,
                                                  database=database)
                if access_user is not False:
                    user_info = check_user(access_user,
                                           getinfo=True,
                                           database=database)
                else:
                    user_info = False

            elif username and not apikey:
                access_user = username
                user_info = check_user(access_user,
                                       getinfo=True,
                                       database=database)


        # make sure this user exists
        if user_info is not False:

            user_accessgroup = user_info[3]

            # generate the query
            if existing_lcs_only:
                query = ('select hat_id, access_groups '
                         'from hat_lightcurves where '
                         'hat_id in %s and full_lc_fpath is not null')
            else:
                query = ('select hat_id, access_groups '
                         'from hat_lightcurves where '
                         'hat_id in %s')

            # enforce uniqueness of hat_ids in input hat_list
            query_params = tuple(set(hat_list))

            try:

                cur.execute(query, (query_params,))
                rows = cur.fetchall()

                if rows and len(rows) > 0:

                    authorized_hatids = [x[0] for x in rows
                                         if user_accessgroup in x[1]]
                    if len(authorized_hatids) == 0:
                        return_tuple = (False,
                                        None,
                                        user_accessgroup,
                                        'no_hatlcs_authorized')
                    else:
                        return_tuple = (True,
                                        authorized_hatids,
                                        user_accessgroup,
                                        'hatlc_auth_ok')
                else:
                    return_tuple = (False,
                                    None,
                                    user_accessgroup,
                                    'no_hatlcs_found')

            except Exception as e:

                if LOGGER:
                    LOGGER.error('authorizing hatlcs for '
                                 'user %s and apikey %s failed, '
                                 'error was: %s' % (access_user, apikey, e))
                if DEBUGMODE:
                    print('authorizing hatlcs for '
                          'user %s and apikey %s failed, '
                          'error was: %s' % (access_user, apikey, e))
                database.rollback()

                return_tuple = (False, None, None, 'database_error')

        # if we can't figure out this user's accessgroups, then we can't
        # authorize access
        else:
            return_tuple = (False, None, None, 'user_not_found')

        # close the database cursor
        database.close_cursor(handle)

    # close the database at the end if we have to
    if closedb:
        database.close_connection()

    return return_tuple



def auth_hatfields(field_list,
                   username=None,
                   apikey=None,
                   database=None):
    '''
    This checks if either the username or the apikey provided is authorized to
    access the HAT fields in hat_list.

    Returns a list of the following form:

    ['hatfield1', 'hatfield2', 'hatfield3', ...]

    hatfieldX MUST be a string! (because we have fields like 088, etc.)

    if these hatfields from the provided field_list are authorized for this
    user. If none of the provided hatfields in field_list are authorized,
    returns a tuple (False, None, 'no_hatfields_authorized')

    '''

    # set up a database connection
    closedb = False

    # if both username and apikey are provided, then return a failure
    if username and apikey:

        return (False, None, 'apikey_and_username_provided')

    # we need either a username or an apikey, but not both. if neither is
    # provided, we use the anonuser to try access
    elif ((username and not apikey) or
          (apikey and not username) or
          (not apikey and not username)):

        if database:
            cur, handle = database.newcursor()
        else:
            database = lcdb.LCDB()
            database.open_default()
            cur, handle = database.newcursor()
            closedb = True

        # figure out the user/apikey info to later compare the accessgroup with
        # the hatfield accessgroups. if neither apikey nor username are
        # provided, use the anonuser
        if not apikey and not username:
            user_info = check_user('anonuser',
                                   getinfo=True,
                                   database=database)
            access_user = 'anonuser'

        # otherwise, if one of them is provided, figure out the corresponding
        # user info
        else:

            if apikey and not username:
                access_user = get_apikey_username(apikey,
                                                  database=database)
                if access_user is not False:
                    user_info = check_user(access_user,
                                           getinfo=True,
                                           database=database)
                else:
                    user_info = False

            elif username and not apikey:
                access_user = username
                user_info = check_user(access_user,
                                       getinfo=True,
                                       database=database)


        # make sure this user exists
        if user_info is not False:

            user_accessgroup = user_info[3]

            # generate the query
            query = ('select field_num, access_groups from hat_fields where '
                     'field_num in %s')
            # enforce uniqueness of fields in input field_list
            query_params = tuple(set(field_list))

            try:

                cur.execute(query, (query_params,))
                rows = cur.fetchall()

                if rows and len(rows) > 0:

                    authorized_hatids = [x[0] for x in rows
                                         if user_accessgroup in x[1]]
                    if len(authorized_hatids) == 0:
                        return_tuple = (False, None, 'no_hatfields_authorized')
                    else:
                        return_tuple = (True, authorized_hatids,
                                        'hatfield_auth_ok')
                else:
                    return_tuple = (False, None, 'no_hatfields_found')

            except Exception as e:

                if LOGGER:
                    LOGGER.error('authorizing hatfields for '
                                 'user %s and apikey %s failed, '
                                 'error was: %s' % (access_user, apikey, e))
                if DEBUGMODE:
                    print('authorizing hatfields for '
                          'user %s and apikey %s failed, '
                          'error was: %s' % (access_user, apikey, e))
                database.rollback()

                return_tuple = (False, None, 'database_error')

        # if we can't figure out this user's accessgroups, then we can't
        # authorize access
        else:
            return_tuple = (False, None, 'user_not_found')

        # close the database cursor
        database.close_cursor(handle)

    # close the database at the end if we have to
    if closedb:
        database.close_connection()

    return return_tuple



def auth_hatlcs_for_hatfields(field_list,
                              username=None,
                              apikey=None,
                              existing_lcs_only=True,
                              database=None):
    '''
    This checks if either the username or the apikey provided is authorized to
    access the HAT LCs in field and returns the authorized HAT IDs.

    Returns a dict of the following form:

    {'hatfield1':[list of hatids],
     'hatfield2':[list of hatids],
     'hatfield3':[list of hatids], ...}

    hatfieldX MUST be a string! (because we have fields like 088, etc.)

    if the hatlcs from the hatfields from the provided field_list are authorized
    for this user. If none of the hatlcs from a field are authorized for this
    user/apikey, it will not be returned as a key in the dict.

    If none of the provided hatfields in field_list are authorized, returns a
    tuple (False, None, 'no_hatfields_authorized')

    '''

    # set up a database connection
    closedb = False

    # if both username and apikey are provided, then return a failure
    if username and apikey:

        return (False, None, 'apikey_and_username_provided')

    # we need either a username or an apikey, but not both. if neither is
    # provided, we use the anonuser to try access
    elif ((username and not apikey) or
          (apikey and not username) or
          (not apikey and not username)):

        if database:
            cur, handle = database.newcursor()
        else:
            database = lcdb.LCDB()
            database.open_default()
            cur, handle = database.newcursor()
            closedb = True

        # figure out the user/apikey info to later compare the accessgroup with
        # the hatfield accessgroups. if neither apikey nor username are
        # provided, use the anonuser
        if not apikey and not username:
            user_info = check_user('anonuser',
                                   getinfo=True,
                                   database=database)
            access_user = 'anonuser'

        # otherwise, if one of them is provided, figure out the corresponding
        # user info
        else:

            if apikey and not username:
                access_user = get_apikey_username(apikey,
                                                  database=database)
                if access_user is not False:
                    user_info = check_user(access_user,
                                           getinfo=True,
                                           database=database)
                else:
                    user_info = False

            elif username and not apikey:
                access_user = username
                user_info = check_user(access_user,
                                       getinfo=True,
                                       database=database)


        # make sure this user exists
        if user_info is not False:

            user_accessgroup = user_info[3]

            # enforce uniqueness of fields in input field_list
            query_params = list(set(field_list))

            # check if we can access all of these fields first
            fieldaccess_info = auth_hatfields(query_params,
                                              username=username,
                                              apikey=apikey,
                                              database=database)

            if fieldaccess_info[0] == True:

                allowed_fields = fieldaccess_info[1]

                # generate the query
                if existing_lcs_only:
                    query = ('select a.field_num, b.hat_id,  '
                             'b.access_groups as lc_access_groups from '
                             'hat_fields a join hat_lightcurves b '
                             'on (a.field_num::integer = b.hat_field) where '
                             '(a.field_num in %s) and '
                             '(b.full_lc_fpath is not null) '
                             'order by a.field_num')
                else:
                    query = ('select a.field_num, b.hat_id,  '
                             'b.access_groups as lc_access_groups from '
                             'hat_fields a join hat_lightcurves b '
                             'on (a.field_num::integer = b.hat_field) where '
                             '(a.field_num in %s) order by a.field_num')


                try:

                    cur.execute(query, (tuple(allowed_fields),))
                    rows = cur.fetchall()

                    if rows and len(rows) > 0:

                        authorized_hatlcs = {}

                        print(allowed_fields)
                        print(user_accessgroup)

                        for field in allowed_fields:

                            hatlcs_for_field = [x[1] for x in rows if
                                                (user_accessgroup in x[2] and
                                                 field == x[0])]
                            if len(hatlcs_for_field) > 0:
                                authorized_hatlcs[field] = hatlcs_for_field
                            else:
                                authorized_hatlcs[field] = None

                        return_tuple = (True, authorized_hatlcs,
                                        'hatfield_and_hatlc_auth_ok')
                    else:
                        return_tuple = (False, None,
                                        'no_authorized_hatlcs_found')

                except Exception as e:

                    if LOGGER:
                        LOGGER.error('authorizing hatfields for '
                                     'user %s and apikey %s failed, '
                                     'error was: %s' % (access_user, apikey, e))
                    if DEBUGMODE:
                        print('authorizing hatfields for '
                              'user %s and apikey %s failed, '
                              'error was: %s' % (access_user, apikey, e))
                    database.rollback()

                    return_tuple = (False, None, 'database_error')

            # if we can't access any of the fields, then return nothing
            else:
                return_tuple = (False, None, 'no_hatfields_authorized')

        # if we can't figure out this user's accessgroups, then we can't
        # authorize access
        else:
            return_tuple = (False, None, 'user_not_found')

        # close the database cursor
        database.close_cursor(handle)

    # close the database at the end if we have to
    if closedb:
        database.close_connection()

    return return_tuple



##############################
## ADMINISTRATION FUNCTIONS ##
##############################

def change_hatlc_accessgroups(hatlc_list,
                              accessgroups='superuser,hatgroup',
                              database=None):
    '''
    This is used to add/remove access restrictions on HAT LCs in hatlc_list.

    accessgroups = string comma-delimited list of access groups allowed to
                   access this LC

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

    query = ('update hat_lightcurves set access_groups = %s '
             'where hat_id in %s returning hat_id, access_groups')

    # enforce uniqueness in input hatlc_list
    query_params = tuple(set(hatlc_list))

    try:

        cur.execute(query, (accessgroups, query_params))
        rows = cur.fetchall()
        database.commit()

        if rows and len(rows) > 0:

            change_success = all(x[1] == accessgroups for x in rows)
            if change_success == True:
                return_tuple = (True, 'all_hatlcs_accessgroups_updated')
            else:
                return_tuple = (False, 'not_all_accessgroups_updated')
        else:
            return_tuple = (False, 'no_matching_hatlcs_found')

    except Exception as e:

        if LOGGER:
            LOGGER.error('changing hatlc accessgroups failed, '
                         'error was: %s' % e)
        if DEBUGMODE:
            print('changing hatlc accessgroups failed, '
                  'error was: %s' % e)
        database.rollback()

        return_tuple = (False, 'database_error')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def change_hatfield_accessgroups(field_list,
                                 accessgroups='superuser,hatgroup',
                                 propagate_to_hatlcs=False,
                                 database=None):
    '''
    This is used to add/remove access restrictions on HAT fields in field_list.

    fields in field_list MUST be strings.

    accessgroups = string comma-delimited list of access groups allowed to
                   access this field

    propagate_to_hatlcs = if True, also sets all HAT LCs for this field to the
                          access_groups as listed in the accessgroups parameter

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

    field_query = ('update hat_fields set access_groups = %s '
                   'where field_num in %s returning field_num, access_groups')
    hatlc_query = ('update hat_lightcurves set access_groups = %s '
                   'where hat_field in %s returning hat_id, access_groups')

    # enforce uniqueness in input field_list
    field_query_params = (accessgroups,
                          tuple(set(field_list)))
    hatlc_query_params = (accessgroups,
                          tuple(set([int(x) for x in field_list])))
    try:

        # first do the field update query

        cur.execute(field_query, field_query_params)
        rows = cur.fetchall()
        database.commit()

        if rows and len(rows) > 0:

            change_success = all(x[1] == accessgroups for x in rows)
            if change_success == True:
                return_tuple = (True, 'all_hatfields_accessgroups_updated')
            else:
                return_tuple = (False, 'not_all_accessgroups_updated')
        else:
            return_tuple = (False, 'no_matching_hatfields_found')

        # then do the hatlc update query if required. if this operation fails,
        # consider the whole thing a failure
        if propagate_to_hatlcs:

            cur.execute(hatlc_query, hatlc_query_params)
            rows = cur.fetchall()
            database.commit()

            if rows and len(rows) > 0:

                change_success = all(x[1] == accessgroups for x in rows)
                if change_success == True:
                    return_tuple = (
                        True,
                        'all_hatfield_and_hatlc_accessgroups_updated'
                        )
                else:
                    return_tuple = (False, 'not_all_accessgroups_updated')
            else:
                return_tuple = (False, 'no_matching_hatlcs_in_hatfields_found')

    except Exception as e:

        if LOGGER:
            LOGGER.error('changing hatlc accessgroups failed, '
                         'error was: %s' % e)
        if DEBUGMODE:
            print('changing hatlc accessgroups failed, '
                  'error was: %s' % e)
        database.rollback()

        return_tuple = (False, 'database_error')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple
