#!/usr/bin/env python

'''
lcauthdb.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Jan 2014

This provides a DB interface for authentication for various lcserver services,
including the lcfinder-service, search-service, and the web frontends.

'''

import logging
import os
import os.path
import subprocess, shlex
import time
import hashlib
import base64

from json import dumps

from passlib.hash import bcrypt_sha256

import lcdb

#############
## LOGGING ##
#############

# setup a logger
LOGGER = logging.getLogger('lcauthdb')
LOGGER.addHandler(logging.NullHandler())

# default to debug mode = False
DEBUGMODE = False

def set_debug(debugbool):
    globals()['DEBUGMODE'] = debugbool


########################
## USEFUL CONFIG DATA ##
########################

def load_access_config(access_conf_file):
    '''
    This loads the access groups configuration info from the file
    access_conf_file.

    Returns a dict.

    '''

# FIXME: move these to a conf file (conf/lcserver-accessgroups.conf) so they can
# be dynamically loaded by the function load_access_config above. This variable
# will then just become the results of a call to that function (the same dict
# format will be returned).
USER_ACCESSGROUPS = {
    'superuser':{'allowed_services':('lcfrontend,'
                                     'searchfrontend,'
                                     'lcdirect,'
                                     'searchdirect,'
                                     'stampsdirect,'
                                     'userpref,'
                                     'lcdash,'
                                     'admindash'),
                 'requests_per_day':-1,
                 'stamp_size':600.0,
                 'search_radius':-1,
                 'search_rows':-1},
    'hatgroup':{'allowed_services':('lcfrontend,'
                                    'searchfrontend,'
                                    'lcdirect,'
                                    'stampsdirect,'
                                    'searchdirect,'
                                    'userpref,'
                                    'lcdash'),
                'requests_per_day':-1,
                'stamp_size':600.0,
                'search_radius':-1,
                'search_rows':-1},
    'hatnet':{'allowed_services':('lcfrontend,'
                                    'searchfrontend,'
                                    'lcdirect,'
                                    'stampsdirect,'
                                    'searchdirect,'
                                    'userpref,'
                                    'lcdash'),
                'requests_per_day':-1,
                'stamp_size':600.0,
                'search_radius':-1,
                'search_rows':-1},
    'hatsouth':{'allowed_services':('lcfrontend,'
                                    'searchfrontend,'
                                    'lcdirect,'
                                    'stampsdirect,'
                                    'searchdirect,'
                                    'userpref,'
                                    'lcdash'),
                'requests_per_day':-1,
                'stamp_size':600.0,
                'search_radius':-1,
                'search_rows':-1},
    'princeton':{'allowed_services':('lcfrontend,'
                                     'searchfrontend,'
                                     'lcdirect,'
                                     'stampsdirect,'
                                     'searchdirect,'
                                     'userpref'),
                 'requests_per_day':1000000,
                 'stamp_size':300.0,
                 'search_radius':30.0,
                 'search_rows':-1},
    'collaborator':{'allowed_services':('lcfrontend,'
                                        'searchfrontend,'
                                        'lcdirect,'
                                        'stampsdirect,'
                                        'searchdirect,'
                                        'userpref'),
                    'requests_per_day':500000,
                    'stamp_size':300.0,
                    'search_radius':10.0,
                    'search_rows':30000000},
    'hnstudent':{'allowed_services':('lcfrontend,'
                                     'searchfrontend,'
                                     'lcdirect,'
                                     'stampsdirect,'
                                     'searchdirect,'
                                     'userpref'),
               'requests_per_day':250000,
               'stamp_size':300.0,
               'search_radius':10.0,
               'search_rows':30000000},
    'hsstudent':{'allowed_services':('lcfrontend,'
                                     'searchfrontend,'
                                     'lcdirect,'
                                     'stampsdirect,'
                                     'searchdirect,'
                                     'userpref'),
               'requests_per_day':250000,
               'stamp_size':300.0,
               'search_radius':10.0,
               'search_rows':30000000},
    'public':{'allowed_services':('lcfrontend,'
                                  'searchfrontend,'
                                  'lcdirect,'
                                  'stampsdirect,'
                                  'searchdirect,'
                                  'userpref'),
              'requests_per_day':50000,
              'stamp_size':60.0,
              'search_radius':5.0,
              'search_rows':500000},
    'anonymous':{'allowed_services':('lcfrontend,'
                                     'searchfrontend,'
                                     'lcdirect,'
                                     'stampsdirect,'
                                     'searchdirect'),
                 'requests_per_day':10000,
                 'stamp_size':60.0,
                 'search_radius':1.0,
                 'search_rows':100000},
    'nouser':{'allowed_services':'none',
              'requests_per_day':0,
              'stamp_size':0.0,
              'search_radius':0.0,
              'search_rows':0},
    }



# defined in order from lowest to highest
ACCESSLEVELS = {'anonymous':1,
                'public':2,
                'princeton':3,
                'hnstudent':4,
                'hsstudent':4,
                'collaborator':5,
                'hatnet':6,
                'hatsouth':6,
                'hatgroup':7,
                'superuser':8}


def access_hierarchy(usergroup):
    '''
    This returns a string representing the levels including and above the
    usergroup.

    '''

    if usergroup in ACCESSLEVELS:
        thislevel = ACCESSLEVELS[usergroup]
        hierarchy = (x for x in ACCESSLEVELS if ACCESSLEVELS[x] >= thislevel)
        return ','.join(hierarchy)
    else:
        return None



###################################################
## AUTHENTICATION AGAINST PASSWORDS AND SERVICES ##
###################################################

def auth_password(username,
                  password,
                  database=None):
    '''
    This returns a tuple of the form:

    (<user_check>, <credential_check>)

    where:

    user_check = True if both user and password exist in the database, False
    otherwise

    credential_check = True if (1) user and password match an existing
    combination in the database and (2) user's email has been verified and (3)
    the user's account is not locked

    A successful authentication requires a returned tuple of (True, True)

    database is an LCDB object containing an already activated Postgres DB
    connection. if None, a new database connection will be created and used.

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

    # handle missing credentials. if any are missing, then substitute fake
    # values. these are compared to the credentials of the 'nouser' user in the
    # database; this test will always fail. this is a bit convoluted, since we
    # want to avoid timing based attacks on the authentication services.
    if not username or len(username) == 0:
        authuser = 'nouser'
    else:
        authuser = username

    if not password or len(username) == 0:
        authpass = 'nopassword'
    else:
        authpass = password

    # get the username and password from the database. we also get the password
    # for the 'nouser' user for a fake comparison later if the specified
    # username isn't in the database.
    sql_query = ("select username, password, email_verified, account_locked "
                 "from lcsweb_users "
                 "where username = %s union all "
                 "select username, password, email_verified, account_locked "
                 "from lcsweb_users "
                 "where username = 'nouser'")

    cur.execute(sql_query, (authuser,))
    rows = cur.fetchall()

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    if rows and len(rows) > 1:

        dbuser, dbpass, emailverified, accountlocked = rows[0]
        passverified = bcrypt_sha256.verify(authpass, dbpass)

        if passverified and emailverified and not accountlocked:
            return (True, True)
        else:
            return (True, False)

    else:

        # if the database returns no rows, then the user doesn't exist in the
        # database. we do a fake password comparison anyway to avoid the leakage
        # of this information.
        nouser, nopass, noemail, nolock = rows[0]
        passverify = bcrypt_sha256.verify('nopassword', nopass)
        return (False, False)



def auth_service_user(session_token,
                      ipaddress,
                      clientheader,
                      service,
                      logrequest=True,
                      database=None):
    '''
    This checks if the session_token associated with a user is allowed to use
    the specified service.

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

    # check if the provided session token exists and is currently active
    session_info = check_user_session(session_token,
                                      getinfo=True,
                                      onlyactive=True,
                                      database=database)

    if session_info is not False:

        # get the username and check if it's authorized to use the requested
        # service
        username = session_info[1]
        user_info = check_user(username, getinfo=True, database=database)
        allowed_services = user_info[4].split(',')
        user_locked = user_info[8]

        if service in allowed_services and not user_locked:
            return_tuple = (True, 'service_auth_ok')

            if logrequest:
                log_auth_request(username,
                                 'session',
                                 ipaddress,
                                 clientheader,
                                 service,
                                 True,
                                 'service_auth_ok',
                                 database=database)

        else:
            return_tuple = (False, 'service_not_allowed')

            if logrequest:
                log_auth_request(username,
                                 'session',
                                 ipaddress,
                                 clientheader,
                                 service,
                                 False,
                                 'service_not_allowed',
                                 database=database)

    else:
        return_tuple = (False, 'session_not_active')

        if logrequest:
            log_auth_request('unknown',
                             'session',
                             ipaddress,
                             clientheader,
                             service,
                             False,
                             'session_not_active',
                             database=database)

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple


def auth_service_apikey(username,
                        apikey,
                        service,
                        database=None):
    '''
    This checks if the apikey associated with a user is allowed to use the
    specified service.

    '''

    if service in ('userprefs','admindash'):
        return (False, 'apikey_denied_for_service')

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    # get info for the apikey
    apikey_info = check_apikey(username,
                               service,
                               getinfo=True,
                               database=database)

    if apikey_info is not False:

        # check if the apikey isn't locked and isn't expired
        apikey_locked = apikey_info[5]
        apikey_expirationdate = apikey_info[2]

        # get the username and check if it's authorized to use the requested
        # service
        user_info = check_user(username, getinfo=True, database=database)
        allowed_services = user_info[4].split(',')
        user_locked = user_info[8]

        # make sure the apikey in the DB matches the one provided, is authorized
        # for the service, its user isn't locked, isn't locked itself, and has
        # not expired
        if (apikey == apikey_info[1] and
            service in allowed_services and
            not user_locked and
            not apikey_locked and
            apikey_expirationdate > time.time()):
            return_tuple = (True, 'service_auth_ok')
        else:
            return_tuple = (False, 'service_not_allowed')

    else:
        return_tuple = (False, 'apikey_not_found')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



###########################################
## USER CREATION, MODIFICATION, DELETION ##
###########################################

def generate_password(maxchar=15):
    '''
    This uses pwgen to come up with a password.

    '''

    pwgen_cmd = shlex.split('pwgen -1 -c -n %s 1' % maxchar)
    pwgen_proc = subprocess.Popen(pwgen_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    pwgen_results, pwgen_stderr = pwgen_proc.communicate()
    password_to_use = pwgen_results.strip('\n')

    return password_to_use



def add_user(username,
             email,
             password=None,
             access_group='public',
             database=None
             ):
    '''This creates a user for LC server webservices.

    username = string containing username (max = 256 characters)

    email = string containing email (max = 1024 characters)

    password = string containing password (max = 2048 characters)
               if None, a 15 character password will be generated

    access_group = a string indicating which access group this user belongs
                    to. access groups are from the following:

                    LOGGED IN AND VALIDATED MANUALLY AND EMAIL VERIFIED

                    hatgroup: unlim reqs/day, unlim searchrad, unlim maxrows
                    princeton: 1M reqs/day, 30deg searchrad, unlim maxrows
                    collaborator: 500K reqs/day, 10deg searchrad, 30M maxrows

                    LOGGED IN and EMAIL VERIFIED

                    public: 100K reqs/day, 5deg searchrad, 5M maxrows

                    NOT SIGNED IN

                    anonymous: 10K reqs/day, 1deg searchrad, 100K maxrows

                    access_groups are also used to enforce restrictions on
                    lightcurve and field access. the user's access_group must
                    match one of the access_groups of a field/lightcurve to be
                    able to see it at all.

    allowed_services = a string containing a comma separated list of LC server
                       web-services this user is allowed to use. these services
                       are:

                       lcfrontend -> frontend for LC retrieval
                       lcdirect -> direct retrieval for LCs
                       searchfrontend -> frontend for search services
                       searchdirect -> direct search request handler
                       userpref -> handler for user preferences
                       admindash -> administrative dashboard
                       lcdash -> lightcurve dashboard for censoring LCs

    Returns (True, <username>, email, <gen_password> if password was None else
             None, <access_group>, <allowed_services>, <profile_id>,
             <auth_message>)

             or

             (False, None, None, None, None, None, None, <auth_message>)
             if user creation fails

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

    # make sure we have a username and an email to process
    if username and email:

        # first, check if this user already exists
        user_exists = check_user(username,
                                 database=database)

        if user_exists:
            return_tuple = (False, None, None, None, None, None, None,
                            'user_already_exists')

        # if the user doesn't exist, then make one
        else:

            # first check if we need to generate a password
            if not password:

                password_to_use = generate_password()
                password_generated = True

            else:

                password_to_use = password
                password_generated = False

            # get this account's service level
            allowed_services = (
                USER_ACCESSGROUPS[access_group]['allowed_services']
                )
            requests_per_day = (
                USER_ACCESSGROUPS[access_group]['requests_per_day']
                )
            search_radius = (
                USER_ACCESSGROUPS[access_group]['search_radius']
                )
            search_rows = (
                USER_ACCESSGROUPS[access_group]['search_rows']
                )

            # hash the password
            pwhash = bcrypt_sha256.encrypt(password_to_use)

            # generate a user_profile_id
            profileid = hashlib.sha512(
                '%s-%s-%s-%s' % (username, email,time.time(),
                                 allowed_services)
                ).hexdigest()

            profileid = str(int(profileid, base=16))
            profileid = int(profileid[:7])

            query = ('insert into lcsweb_users '
                     '(username, email, access_group, allowed_services, '
                     'max_requests_per_day, max_search_radius, '
                     'max_search_rows, password, profile_id) values '
                     '(%s, %s, %s, %s, %s, %s, %s, %s, %s)')

            try:
                cur.execute(query,
                            (username, email, access_group, allowed_services,
                             requests_per_day, search_radius, search_rows,
                             pwhash, profileid))
                database.commit()

                return_tuple = (True,
                                username,
                                email,
                                password_to_use if password_generated else None,
                                access_group,
                                allowed_services,
                                profileid,
                                'user_created_ok')

            except Exception as e:

                if LOGGER:
                    LOGGER.error('creating user %s failed, error was: %s'
                                 % (username, e))
                if DEBUGMODE:
                    print('creating user %s failed, error was: %s'
                                 % (username, e))
                database.rollback()

                return_tuple = (False, None, None, None, None, None, None,
                                'database_error')

    else:
        return_tuple = (False, None, None, None, None, None, None,
                        'no_user_or_password')


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def check_user(username,
               getinfo=False,
               getpass=False,
               database=None):
    '''
    This checks if a username exists in the database already. If getinfo is
    True, returns the user's information.

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

    if getinfo:
        query = ('select username, email, email_verified, access_group, '
                 'allowed_services, max_requests_per_day, max_search_radius, '
                 'max_search_rows, account_locked, failed_logins, '
                 'account_createdate, account_lockdate, '
                 'profile_id, password_reset_token from lcsweb_users '
                 'where username = %s')
    elif getpass:
        query = ('select username, password, email, email_verified, '
                 'account_locked, failed_logins, account_createdate, '
                 'account_lockdate, password_reset_token from lcsweb_users '
                 'where username = %s')
    else:
        query = 'select username from lcsweb_users where username = %s'

    try:

        cur.execute(query, (username,))
        rows = cur.fetchall()

        if rows and len(rows) > 0:
            if getinfo or getpass:
                returnval = rows[0]
            else:
                returnval = True
        else:
            returnval = False

    except Exception as e:

        if LOGGER:
            LOGGER.error('checking for user %s failed, error was %s'
                         % (username, e))
        if DEBUGMODE:
            print('checking for user %s failed, error was %s'
                         % (username, e))
        database.rollback()

        returnval = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return returnval



def remove_user(username,
                database=None):
    '''
    This gets rid of a user entirely.

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

    # delete ALL traces of this user except for session/auth histories
    queries = ['delete from lcsweb_users where username = %s',
               'delete from lcsweb_apikeys where username = %s',
               'delete from lcsweb_apikey_usage_history where username = %s',
               'delete from lcsweb_user_task_history where username = %s']


    try:

        for query in queries:
            cur.execute(query, (username,))
        database.commit()
        user_exists = check_user(username,
                                 database=database)
        if user_exists:
            removal_ok = False
        else:
            removal_ok = True

    except Exception as e:

        if LOGGER:
            LOGGER.error('deleting user %s failed, error was %s'
                         % (username, e))
        if DEBUGMODE:
            print('deleting user %s failed, error was %s'
                         % (username, e))
        database.rollback()

        removal_ok = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return removal_ok


def locktoggle_user(username,
                    lock=True,
                    database=None):
    '''
    This locks/unlocks a user out from access to the LC server.

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

    if lock == True:
        query = ('update lcsweb_users set account_locked = %s, '
                 'account_lockdate = current_timestamp '
                 'where username = %s returning account_locked')
    elif lock == False:
        query = ('update lcsweb_users set account_locked = %s, '
                 'account_lockdate = NULL '
                 'where username = %s returning account_locked')


    try:

        cur.execute(query, (lock, username))
        rows = cur.fetchall()
        database.commit()

        account_locked = rows[0][0]

        if account_locked == lock:
            locktoggle_ok = True
        else:
            locktoggle_ok = False

    except Exception as e:

        if LOGGER:
            LOGGER.error('lock toggle for user %s failed, error was: %s'
                         % (username, e))
        if DEBUGMODE:
            print('lock toggle for user %s failed, error was %s'
                         % (username, e))
        database.rollback()

        locktoggle_ok = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return locktoggle_ok



def update_user_password(username,
                         old_password,
                         new_password=None,
                         database=None):
    '''
    Updates the password for a user.

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

    # first, we check if the user's old password matches the user's entered
    # password
    user_ok, pass_ok = auth_password(username, old_password)

    if user_ok and pass_ok:

        if new_password:
            new_password_to_use = new_password
            password_generated = False
        else:
            new_password_to_use = generate_password()
            password_generated = True

        hashed_password = bcrypt_sha256.encrypt(new_password_to_use)

        query = ('update lcsweb_users set password = %s '
                 'where username = %s')

        try:

            cur.execute(query, (hashed_password, username))
            database.commit()

            return_tuple = (True,
                            new_password_to_use if password_generated else None)

        except Exception as e:

            if LOGGER:
                LOGGER.error('lock toggle for user %s failed, error was: %s'
                             % (username, e))
            if DEBUGMODE:
                print('lock toggle for user %s failed, error was %s'
                             % (username, e))
            database.rollback()

            return_tuple = (False, None)

    else:

            return_tuple = (False, None)

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple


def reset_user_password(username,
                        resettoken,
                        newpassword,
                        database=None):
    '''
    This is used to reset a user's password.

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

    # first, we check if this user exists
    user_info = check_user(username, getpass=True, database=database)
    fake_info = check_user('nouser', getpass=True, database=database)

    exec_fake_hashes = False

    # if this user exists, then check the provided reset token against the
    # stored reset token
    if user_info is not False:

        # check the password reset token against what we have on file
        passreset_ok = bcrypt_sha256.verify(resettoken, user_info[8])

        # check if the provided new password isn't the same as the old password
        newpass_is_same = bcrypt_sha256.verify(newpassword, user_info[1])

        # if the reset token checks out, then set the password to the new
        # password entered by the user
        if passreset_ok and not newpass_is_same:

            # hash the new password
            newpasshash = bcrypt_sha256.encrypt(newpassword)

            # generate a new password reset token for the user
            newpassresettoken = generate_password(maxchar=20)
            resettoken_hash = bcrypt_sha256.encrypt(newpassresettoken)

            # prepare the query
            passchange_query = ('update lcsweb_users set password = %s, '
                                'password_reset_token = %s where '
                                'username = %s')
            passchange_params = (newpasshash, resettoken_hash, username)

            # execute the query
            try:
                cur.execute(passchange_query, passchange_params)
                database.commit()

                return_tuple = (True, newpassresettoken, 'password_reset_ok')

            except Exception as e:

                if LOGGER:
                    LOGGER.error('changing user %s failed, error was: %s'
                                 % (username, e))
                if DEBUGMODE:
                    print('creating user %s failed, error was: %s'
                                 % (username, e))
                database.rollback()

                return_tuple = (False, None, 'database_error')
                exec_fake_hashes = True

        else:

            return_tuple = (False, None, 'reset_token_or_newpass_invalid')
            exec_fake_hashes = True

    else:

       return_tuple = (False, None, 'user_does_not_exist')
       exec_fake_hashes = True


    if exec_fake_hashes:

        # do fake hashes to let no information leak to the outside

        # verifying current reset token
        fake_hash = bcrypt_sha256.verify(resettoken, fake_info[1])

        # verifying new password != old password
        fake_hash = bcrypt_sha256.verify(resettoken, fake_info[1])

        # hashing the new password
        fake_hash = bcrypt_sha256.verify(resettoken, fake_info[1])

        # hashing the new reset token
        fake_hash = bcrypt_sha256.verify(resettoken, fake_info[1])


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def update_user_info(username,
                     email=None,
                     access_group=None,
                     allowed_services=None,
                     email_verified=None,
                     requests_per_day=None,
                     search_radius=None,
                     search_rows=None,
                     profile_id=None,
                     password_reset_token=None,
                     database=None):
    '''
    This updates a user's email, access_groups, or allowed_services.

    Returns a tuple where the first element is the result (True or False) of the
    update, and the other elements are the newly updated columns for the user in
    the users table.

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

    # first, check if this user exists
    user_info = check_user(username,
                             getinfo=True,
                             database=database)

    if user_info is not False:

        columns_to_update = []
        query_params = []
        new_column_values = []

        if email is not None:
            query_params.extend(['email = %s', 'email_verified = %s'])
            columns_to_update.extend(['email','email_verified'])
            new_column_values.extend([email, False])

        if access_group is not None:
            query_params.append('access_group = %s')
            columns_to_update.append('access_group')
            new_column_values.append(access_group)

        if allowed_services is not None:
            query_params.append('allowed_services = %s')
            columns_to_update.append('allowed_services')
            new_column_values.append(allowed_services)

        if email_verified is not None:
            query_params.append('email_verified = %s')
            columns_to_update.append('email_verified')
            new_column_values.append(email_verified)

        if requests_per_day is not None:
            query_params.append('max_requests_per_day = %s')
            columns_to_update.append('max_requests_per_day')
            new_column_values.append(requests_per_day)

        if search_radius is not None:
            query_params.append('max_search_radius = %s')
            columns_to_update.append('max_search_radius')
            new_column_values.append(search_radius)

        if search_rows is not None:
            query_params.append('max_search_rows = %s')
            columns_to_update.append('max_search_rows')
            new_column_values.append(search_rows)

        # update or generate a new profile_id
        if profile_id is not None:
            query_params.append('profile_id = %s')
            columns_to_update.append('profile_id')
            if profile_id == 'generatenew':
                new_profile_id = hashlib.sha512(
                    '%s-%s-%s-%s' % (user_info[0],
                                     user_info[1],
                                     user_info[10],
                                     user_info[4])
                    ).hexdigest()
                new_profile_id = str(int(new_profile_id, base=16))
                new_profile_id = int(new_profile_id[:7])
            else:
                new_profile_id = profile_id
            new_column_values.append(new_profile_id)

        # update or generate a new password_reset_token
        if password_reset_token is not None:
            query_params.append('password_reset_token = %s')
            columns_to_update.append('password_reset_token')
            if password_reset_token == 'generatenew':
                new_password_reset_token = generate_password(maxchar=20)
            else:
                new_password_reset_token = password_reset_token

            # hash this in the same way as a password since it's effectively a
            # password
            new_password_reset_token_hash = bcrypt_sha256.encrypt(
                new_password_reset_token
                )

            new_column_values.append(new_password_reset_token_hash)

        if (len(columns_to_update) > 0 and
            len(query_params) > 0 and
            len(new_column_values) > 0):

            update_placeholders = ', '.join(query_params)
            column_placeholders = ', '.join(columns_to_update)

            query = ('update lcsweb_users set {update_placeholders} '
                     'where username = %s '
                     'returning {column_placeholders}')

            query = query.format(update_placeholders=update_placeholders,
                                 column_placeholders=column_placeholders)

            # execute the query
            try:

                cur.execute(query, tuple(new_column_values + [username]))
                rows = cur.fetchall()
                database.commit()

                return_tuple = list(rows[0])

                # we need to handle the password_reset_token specially
                if 'password_reset_token' in columns_to_update:
                    reset_token_index = columns_to_update.index(
                        'password_reset_token'
                        )
                    return_tuple[reset_token_index] = new_password_reset_token

                return_tuple = tuple([True] + return_tuple)

            except Exception as e:

                if LOGGER:
                    LOGGER.error('updating info for user %s failed, '
                                 'error was: %s'
                                 % (username, e))
                    LOGGER.error('query was: %s' % cur.query)
                if DEBUGMODE:
                    print('updating info for user %s failed, error was %s'
                                 % (username, e))
                    print('query was: %s' % cur.query)
                database.rollback()

                return_tuple = (False,)

        # nothing was changed, successfully return
        else:
            return_tuple = (True,)

    # the user does not exist, return a failure
    else:
        return_tuple = (False,)

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



#############################################
## APIKEY CREATION, MODIFICATION, DELETION ##
#############################################

def make_apikey():
    '''
    This just generates an apikey when called.

    '''

    return (
        base64.b64encode(hashlib.sha256(os.urandom(12)).hexdigest()).strip('=')
        )



def generate_apikey(username,
                    service,
                    expires=365,
                    database=None):
    '''
    This generates an apikey for a user and a specific service. By default, the
    apikey will expire after 365 days, otherwise: expires is number of days
    indicating when the apikey should expire. The following services cannot be
    accessed with an apikey:

    userpref
    admindash
    lcfrontend
    searchfrontend

    The following services allow apikey usage:

    lcdirect
    searchdirect


    '''

    # make sure no one can generate an apikey for user-interactive services, we
    # may change this later, so admins can change lots of people's info at once
    # using an API call, but for now, this is more suitable
    if service in ('userpref', 'admindash','lcfrontend','searchfrontend'):
        return (False, None, None, 'apikey_denied_for_service')

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    # first check if this user is legit
    user_info = check_user(username, getinfo=True, database=database)

    if user_info is not False:

        account_locked = user_info[8]
        allowed_services = user_info[4].split(',')

        # if the account's locked or the requested service isn't allowed, we
        # can't make an apikey
        if account_locked or service not in allowed_services:

            return_tuple = (False, None, None,
                            'account_locked_or_service_denied')

        else:

            apikey_info = check_apikey(username, service,
                                       getinfo=True,
                                       database=database)

            # if an apikey for this user and service already exists, is not
            # locked, and has not expired, return it since we don't need to do
            # anything else
            if (apikey_info is not False and
                apikey_info[2] > time.time() and
                apikey_info[5] == False):

                return_tuple = (True, apikey_info[1], service,
                                'apikey_exists_and_ok')

            # if an apikey for this user and service already exists, may be
            # expired, but is locked, return a failure since it's probably been
            # locked for abuse
            elif (apikey_info is not False and apikey_info[5] == True):

                return_tuple = (False, None, None, 'apikey_exists_but_locked')

            # otherwise, the apikey doesn't exist; order up a new apikey for
            # this user and service
            elif apikey_info is False:

                reqs_per_day = user_info[5]
                expiration_date = time.time() + expires*86400.0

                # generate an apikey
                user_apikey = make_apikey()

                # prepare the database query
                query = ('insert into lcsweb_apikeys '
                         '(username, apikey, expires_unixtime, '
                         'for_service, requests_per_day) values '
                         '(%s, %s, %s, %s, %s)')

                query_params = (username,
                                user_apikey,
                                expiration_date,
                                service,
                                reqs_per_day)

                try:

                    cur.execute(query, query_params)
                    database.commit()

                    return_tuple = (True, user_apikey, service,
                                    'apikey_created_ok')

                except Exception as e:

                    if LOGGER:
                        LOGGER.error('creating apikey for user %s failed, '
                                     'error was: %s'
                                     % (username, e))
                    if DEBUGMODE:
                        print('creating apikey for user %s failed, '
                              'error was: %s'
                              % (username, e))
                    database.rollback()

                    return_tuple = (False, None, None, 'database_error')

            # finally, if the apikey does exist but is expired, return a
            # failure, since we need to remove the expired key first
            elif (apikey_info is not False and
                  apikey_info[2] < time.time()):

                return_tuple = (False, None, None, 'apikey_exists_but_expired')

            # we should never get here
            else:

                return_tuple = (False, None, None, 'wtf_error')

    # this user doesn't exist, return a failure
    else:

        return_tuple = (False, None, None, 'nonexistent_user')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def check_apikey(username, service,
                 getinfo=False,
                 database=None):
    '''
    This checks if an apikey exists for the user username and requested
    service. Returns the apikey if it exists, otherwise returns False.

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

    if getinfo:
        query = ('select username, apikey, expires_unixtime, '
                 'for_service, requests_per_day, apikey_locked '
                 'from lcsweb_apikeys '
                 'where username = %s and for_service = %s')
    else:
        query = ('select apikey from lcsweb_apikeys where '
                 'username = %s and for_service = %s')

    try:

        cur.execute(query, (username, service))
        rows = cur.fetchall()

        if rows and len(rows) > 0:
            if getinfo:
                returnval = rows[0]
            else:
                returnval = True
        else:
            returnval = False

    except Exception as e:

        if LOGGER:
            LOGGER.error('checking for apikey for user %s failed, '
                         'error was %s'
                         % (username, e))
        if DEBUGMODE:
            print('checking for apikey for user %s failed, '
                  'error was %s'
                  % (username, e))
        database.rollback()

        returnval = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return returnval



def get_apikey_username(apikey, database=None):
    '''
    This just gets the username associated with the given apikey.

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

    query = ('select username from lcsweb_apikeys '
             'where apikey = %s')

    try:

        cur.execute(query, (apikey,))
        rows = cur.fetchall()

        if rows and len(rows) > 0:
            returnval = rows[0]
        else:
            returnval = False

    except Exception as e:

        if LOGGER:
            LOGGER.error('checking for username for apikey %s failed, '
                         'error was %s'
                         % (apikey, e))
        if DEBUGMODE:
            print('checking for username for apikey %s failed, '
                  'error was %s'
                  % (apikey, e))
        database.rollback()

        returnval = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return returnval



def locktoggle_apikey(username,
                      service,
                      lock=True,
                      database=None):
    '''
    This locks/unlocks an apikey for a user.

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

    if lock == True:
        query = ('update lcsweb_apikeys set apikey_locked = %s '
                 'where username = %s and for_service = %s '
                 'returning apikey_locked')
    elif lock == False:
        query = ('update lcsweb_apikeys set apikey_locked = %s '
                 'where username = %s and for_service = %s '
                 'returning apikey_locked')

    try:

        cur.execute(query, (lock, username, service))
        rows = cur.fetchall()
        database.commit()

        apikey_locked = rows[0][0]

        if apikey_locked == lock:
            locktoggle_ok = True
        else:
            locktoggle_ok = False

    except Exception as e:

        if LOGGER:
            LOGGER.error('lock toggle for apikey for user %s, '
                         'service %s failed, '
                         'error was: %s'
                         % (username, service, e))
        if DEBUGMODE:
            print('lock toggle for apikey for user %s, '
                  'service %s failed, '
                  'error was: %s'
                  % (username, service, e))
        database.rollback()

        locktoggle_ok = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return locktoggle_ok



def remove_apikey(username,
                  service,
                  database=None):
    '''
    This remove an apikey for a user. If this is manually done, then the user
    must be logged in. The only other way to remove an apikey is by a superuser,
    or to let it expire on its own.

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

    query = ('delete from lcsweb_apikeys where '
             'username = %s and for_service = %s')

    try:

        cur.execute(query, (username, service))
        database.commit()
        apikey_exists = check_apikey(username, service, database=database)
        if apikey_exists:
            removal_ok = False
        else:
            removal_ok = True

    except Exception as e:

        if LOGGER:
            LOGGER.error('deleting apikey for user %s, '
                         'service %s failed, error was %s'
                         % (username, service, e))
        if DEBUGMODE:
            print('deleting apikey for user %s, '
                  'service %s failed, error was %s'
                  % (username, service, e))
        database.rollback()

        removal_ok = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return removal_ok



###########################
## USER SESSION HANDLING ##
###########################

def make_session_token(username,
                       ipaddress,
                       clientheader):
    '''
    This generates a new session token.

    '''

    tokeninfo = '%s|%s|%s|%s' % (username,
                                 ipaddress,
                                 clientheader,
                                 time.time())

    return hashlib.sha512(tokeninfo).hexdigest()



def initiate_anonymous_session(ipaddress,
                               clientheader,
                               database=None):
    '''This function is used to initiate an anonymous session for the HAT LC
    server.

    The anon user workflow is like so:

    1. user hits any page of the lcserver. we check for the hatlcserver_session
       cookie.

    2. if hatlcserver_session is None, we run this function to insert a new
       session into the database, and return the session token to the frontend,
       which sets hatlcserver_session cookie = session_token, with expiry set
       for session end (upon which we call logout_user_session) and clear the
       cookie.

       if the hatlcserver_session is not None, we check it against the
       lcsweb_user_sessions table. if nothing is found, the frontend clears the
       hatlcserver_session cookie, calls this function to get a new anonymous
       session cookie, and sets it. if something is found, we do nothing.

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

    session_user = 'anonuser'

    # get a session token
    user_sessiontoken = make_session_token(session_user,
                                           ipaddress,
                                           clientheader)
    user_logintime = time.time()

    # prepare the query to insert a new session into the
    # lcsweb_user_sessions table
    sessions_query = ('insert into lcsweb_user_sessions '
                      '(username, login_unixtime, session_token, '
                      'ip_address, client_header, session_active) '
                      'values '
                      '(%s, %s, %s, %s, %s, %s)')
    sessions_query_params = (session_user,
                             user_logintime,
                             user_sessiontoken,
                             ipaddress,
                             clientheader,
                             True)

    # prepare the query to insert a new session into the
    # lcsweb_user_session_history table
    sessionhistory_query = ('insert into lcsweb_user_session_history '
                            '(session_token, username, login_unixtime, '
                            'ip_address, client_header) values '
                            '(%s, %s, %s, %s, %s)')
    sessionhistory_query_params = (user_sessiontoken,
                                   session_user,
                                   user_logintime,
                                   ipaddress,
                                   clientheader)

    try:

        # execute the queries
        cur.execute(sessions_query, sessions_query_params)
        cur.execute(sessionhistory_query, sessionhistory_query_params)
        database.commit()

        return_tuple = (True, user_sessiontoken)

    except Exception as e:

        if LOGGER:
            LOGGER.error('initiating session for user %s failed,'
                         ' error was: %s' % (session_user, e))
            LOGGER.error('query was: %s' % cur.query)
        if DEBUGMODE:
            print('initiating session for user %s failed, '
                  'error was %s' % (session_user, e))
            print('query was: %s' % cur.query)

        database.rollback()
        return_tuple = (False, None)

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def authenticate_user_session(ipaddress,
                              clientheader,
                              username=None,
                              password=None,
                              database=None):
    '''
    This authenticates a user against the auth database using either the
    password. If username is None, then the anonuser user is used (which only
    has access to publicly and anonymously available data and services). This
    function is to process logins or anonymous use and initiates a user session
    only. Each service will process its access separately using the
    auth_service_user function above.

    Returns a tuple:

    (<auth_result>, <session_token>, <profile_id>)

    <auth_result> is True or False
    <session_token> is the session token generated for this user
    <profile_id> is the user's ID to get them to their personal profile page

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

    # if we have both a username and password
    if username and password:

        session_user = username
        session_pass = password

    # if we have neither, then use the anonuser
    elif not username and not password:

        session_user = 'anonuser'
        session_pass = 'lohrizohph0Quei'

    # if one of the values is missing, then sub dummy values in
    else:

        session_user = 'nouser'
        session_pass = 'nopassword'

    # check the user's credentials
    user_ok, pass_ok = auth_password(session_user,
                                     session_pass,
                                     database=database)

    # make sure we have both the user and password correct
    if user_ok and pass_ok:

        # get a session token
        user_sessiontoken = make_session_token(session_user,
                                               ipaddress,
                                               clientheader)
        user_logintime = time.time()

        # prepare the query to insert a new session into the
        # lcsweb_user_sessions table
        sessions_query = ('insert into lcsweb_user_sessions '
                          '(username, login_unixtime, session_token, '
                          'ip_address, client_header, session_active) '
                          'values '
                          '(%s, %s, %s, %s, %s, %s)')
        sessions_query_params = (session_user,
                                 user_logintime,
                                 user_sessiontoken,
                                 ipaddress,
                                 clientheader,
                                 True)

        # prepare the query to insert a new session into the
        # lcsweb_user_session_history table
        sessionhistory_query = ('insert into lcsweb_user_session_history '
                                '(session_token, username, login_unixtime, '
                                'ip_address, client_header) values '
                                '(%s, %s, %s, %s, %s)')
        sessionhistory_query_params = (user_sessiontoken,
                                       session_user,
                                       user_logintime,
                                       ipaddress,
                                       clientheader)

        # prepare the query to get the user's profile id
        profileid_query = ('select profile_id from lcsweb_users where '
                           'username = %s')
        profileid_query_params = (session_user,)

        try:

            # execute the queries
            cur.execute(sessions_query, sessions_query_params)
            cur.execute(sessionhistory_query, sessionhistory_query_params)
            database.commit()

            cur.execute(profileid_query, profileid_query_params)
            rows = cur.fetchone()
            user_profileid = rows[0]

            return_tuple = (True, user_sessiontoken, user_profileid)

        except Exception as e:

            if LOGGER:
                LOGGER.error('initiating session for user %s failed,'
                             ' error was: %s' % (session_user, e))
                LOGGER.error('query was: %s' % cur.query)
            if DEBUGMODE:
                print('initiating session for user %s failed, '
                      'error was %s' % (session_user, e))
                print('query was: %s' % cur.query)

            database.rollback()
            return_tuple = (False, None, None)

    # otherwise, we return a denial
    else:
        return_tuple = (False, None, None)


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def logout_user_session(session_token,
                        database=None):
    '''
    This logs out the user with session_token from all services. If the user is
    anonymous, then the user is automatically logged out when they leave the
    page (this function will be called automatically when on-page JS fires an
    AJAX request indicating that the session was closed).

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

    # we need a query to remove the session_token from the lcsweb_user_sessions
    # table, and another to update the lcsweb_user_session_history table
    sessions_query = ('delete from lcsweb_user_sessions '
                      'where session_token = %s')
    sessionhistory_query = ('update lcsweb_user_session_history set '
                            'logout_unixtime = %s where session_token = %s')

    try:

        # execute the queries
        cur.execute(sessions_query, (session_token,))
        cur.execute(sessionhistory_query, (time.time(),
                                           session_token))
        database.commit()

        session_exists = check_user_session(session_token)

        if session_exists:
            logout_ok = False
        else:
            logout_ok = True

    except Exception as e:

        if LOGGER:
            LOGGER.error('removing session %s failed, error was %s'
                         % (session_token, e))
        if DEBUGMODE:
            print('removing session %s failed, error was %s'
                         % (session_token, e))
        database.rollback()

        logout_ok = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return logout_ok



def check_user_session(session_token,
                       getinfo=False,
                       onlyactive=False,
                       database=None):
    '''
    This checks if the session_token provided is present in the
    lcsweb_user_sessions table. If onlyactive = True, then the session's
    session_active column must be True for this function to return the session
    as valid. If getinfo = True, returns all session info for this
    session_token.

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

    if getinfo:
        query = ('select session_token, username, login_unixtime, '
                 'ip_address, client_header, session_active '
                 'from lcsweb_user_sessions where session_token = %s '
                 '{onlyactive_param}')
    else:
        query = ('select session_token from lcsweb_user_sessions '
                 'where session_token = %s {onlyactive_param}')

    if onlyactive:
        query = query.format(onlyactive_param='and session_active = TRUE')
    else:
        query = query.format(onlyactive_param='')

    try:

        cur.execute(query, (session_token,))
        rows = cur.fetchall()

        if rows and len(rows) > 0:
            if getinfo:
                returnval = rows[0]
            else:
                returnval = True

        else:
            returnval = False

    except Exception as e:

        if LOGGER:
            LOGGER.error('checking for session %s failed, error was %s'
                         % (session_token, e))
        if DEBUGMODE:
            print('checking for session %s failed, error was %s'
                         % (session_token, e))
        database.rollback()

        returnval = False

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return returnval



def update_user_session(session_token,
                        ipaddress=None,
                        clientheader=None,
                        sessionactive=None,
                        database=None):
    '''
    This updates the session information for a session_token. It can be used to
    update the client's ip_address, client_header, or toggle between the active
    and nonactive states of the session. 'active' means the user is currently on
    the page; this will be set to False when the user leaves the page via an
    AJAX notification from the frontend's JS.

    Returns a tuple where the first element is the result (True or False) of the
    update, and the other elements are the newly updated columns for the session
    in the user_sessions table.

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


    # first, check if this session exists
    session_exists = check_user_session(session_token,
                                        database=database)

    if session_exists:

        columns_to_update = []
        query_params = []
        new_column_values = []

        if ipaddress is not None:
            query_params.append('ip_address = %s')
            columns_to_update.append('ip_address')
            new_column_values.append(ipaddress)

        if clientheader is not None:
            query_params.append('client_header = %s')
            columns_to_update.append('client_header')
            new_column_values.append(clientheader)

        if sessionactive is not None:
            query_params.append('session_active = %s')
            columns_to_update.append('session_active')
            new_column_values.append(sessionactive)

        if (len(columns_to_update) > 0 and
            len(query_params) > 0 and
            len(new_column_values) > 0):

            update_placeholders = ', '.join(query_params)
            column_placeholders = ', '.join(columns_to_update)

            sessions_query = ('update lcsweb_user_sessions set '
                              '{update_placeholders} where '
                              'session_token = %s '
                              'returning {column_placeholders}')

            sessions_query = sessions_query.format(
                update_placeholders=update_placeholders,
                column_placeholders=column_placeholders
                )

            # execute the query
            try:

                cur.execute(sessions_query,
                            tuple(new_column_values + [session_token]))
                rows = cur.fetchall()
                database.commit()

                return_tuple = [True] + list(rows[0])
                return_tuple = tuple(return_tuple)

            except Exception as e:

                if LOGGER:
                    LOGGER.error('updating info for session %s failed, '
                                 'error was: %s'
                                 % (session_token, e))
                    LOGGER.error('query was: %s' % cur.query)
                if DEBUGMODE:
                    print('updating info for session %s failed, error was %s'
                                 % (session_token, e))
                    print('query was: %s' % cur.query)
                database.rollback()

                return_tuple = (False,)

        # nothing to update, return a success
        else:
            return_tuple = (True,)

    # the session doesn't exist, return False
    else:
        return_tuple = (False,)

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def sweep_old_sessions(minage_hours=24,
                       keep_active=True,
                       database=None):

    '''
    This removes sessions older than minage_hours from the lcsweb_user_sessions
    table. Used by the periodic session sweep functionality.

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

    query = ('delete from lcsweb_user_sessions where '
             '((%s - login_unixtime) < %s){keepactive_option}')
    if keep_active:
        keepactive_option = ' and (session_active = False)'
    else:
        keepactive_option = ''

    query = query.format(keepactive_option=keepactive_option)
    query_params = (time.time(),
                    minage_hours * 3600.0)

    # delete old sessions and return the count of existing sessions
    try:
        cur.execute(query, query_params)
        database.commit()

        query = 'select count(*) from lcsweb_user_sessions'
        cur.execute(query)
        rows = cur.fetchone()

        if rows and len(rows) > 0:
            return_tuple = (True, rows[0], 'session_sweep_ok')
        else:
            return_tuple = (False, None, 'session_sweep_failed_db_error')

    except Exception as e:

        if LOGGER:
            LOGGER.error('cannot sweep sessions, error was: %s' % e)
        if DEBUGMODE:
            print('cannot sweep sessions, error was: %s' % e)
        database.rollback()

        return_tuple = (False, None, 'database_error')


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple


#####################
## APIKEY HANDLING ##
#####################

def authenticate_apikey(service,
                        ipaddress,
                        clientheader,
                        username=None,
                        apikey=None,
                        database=None,
                        logrequest=True):

    '''
    This authenticates a user against the auth database using the apikey. If the
    apikey is None, then the apikey associated with the anonymous user is used
    (which only has access to public and anonymously available data and
    services). This uses the auth_service_apikey function above and updates the
    lcsweb_apikey_usage_history table in the database. Also updates the
    lcsweb_authentication_request_history table if logrequest = True.

    Returns a tuple:

    (<auth_result>, <apikey>, <service>, <status_message>)

    <auth_result> is True or False
    <apikey> is the API key of the user
    <service> is the service this apikey was authenticated for
    <status_message> is a message indicating success or reason for failure

    '''

    # we need both username and apikey or neither of them, can't proceed if only
    # one of them is missing and the other is provided
    if (not username and apikey) or (apikey and not username):

        if logrequest:
            log_auth_request(username,
                             'apikey',
                             ipaddress,
                             clientheader,
                             service,
                             False,
                             'user_or_apikey_missing',
                             database=database)

        return (False, None, None, 'user_or_apikey_missing')

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    # if no credentials provided, use the anonuser's apikey
    if not username and not apikey:
        apikey_info = check_apikey('anonuser',service,
                                   getinfo=True,
                                   database=database)

        if apikey_info is not False:
            user_to_check = 'anonuser'
            apikey_to_check = apikey_info[1]
        else:
            user_to_check = 'nouser',
            apikey_to_check = 'noapikey'

    else:

        user_to_check = username
        apikey_to_check = apikey

    # check the user, apikey, service combo
    apikey_auth_result = auth_service_apikey(user_to_check,
                                             apikey_to_check,
                                             service,
                                             database=database)

    # if the authentication for this apikey succeeds, then update the
    # appropriate tables in the DB
    if apikey_auth_result[0] == True:

        query = ('insert into lcsweb_apikey_usage_history '
                 '(apikey, username, service_used, used_unixtime, '
                 'ip_address, client_header) values '
                 '(%s, %s, %s, %s, %s, %s)')
        query_params = (apikey_to_check,
                        user_to_check,
                        service,
                        time.time(),
                        ipaddress,
                        clientheader)

        try:

            cur.execute(query, query_params)
            database.commit()



            return_tuple = (True, apikey_to_check,
                            service, apikey_auth_result[1])

            if logrequest:
                log_auth_request(username,
                                 'apikey',
                                 ipaddress,
                                 clientheader,
                                 service,
                                 True,
                                 apikey_auth_result[1],
                                 database=database)

        except Exception as e:

            if LOGGER:
                LOGGER.error('apikey usage for user %s, service %s failed,'
                             ' error was: %s' % (apikey_to_check,
                                                 service,
                                                 e))
                LOGGER.error('query was: %s' % cur.query)
            if DEBUGMODE:
                print('apikey usage for user %s, service %s failed,'
                      ' error was: %s' % (apikey_to_check,
                                          service,
                                          e))
                print('query was: %s' % cur.query)

            database.rollback()
            return_tuple = (False, None, None, 'database_error')

            if logrequest:
                log_auth_request(username,
                                 'apikey',
                                 ipaddress,
                                 clientheader,
                                 service,
                                 False,
                                 'database_error',
                                 database=database)


    # apikey auth failed, tell the caller
    else:
        return_tuple = (False, None, None, apikey_auth_result[1])

        if logrequest:
            log_auth_request(username,
                             'apikey',
                             ipaddress,
                             clientheader,
                             service,
                             False,
                             apikey_auth_result[1],
                             database=database)

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple


##################################################
## VERIFICATION TOKENS FOR SIGNING UP NEW USERS ##
##################################################

def make_verification_token():
    '''
    This just generates a verification token in URL safe text.

    '''

    return hashlib.sha256(os.urandom(16)).hexdigest()



def new_verification_token(email,
                           ipaddress,
                           clientheader,
                           database=None):
    '''
    This generates a verification token for email and inserts it into the
    database. The token expires in 15 minutes by default (expires is minutes for
    expiry). Uses the lcsweb_user_verifications table in the DB.

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

    # first check if there's an existing user with this email address
    # (remember we're using email addresses as usernames)
    usercheck = check_user(email,
                           getinfo=True,
                           database=database)

    # if the user exists, then we need to check if the account has been verified
    if usercheck is not False:

        # if the user's email has been verified, but the account is locked, the
        # user has been active in the past, but now is locked out; must contact
        # us to unlock
        if usercheck[2] == True and usercheck[8] == True:
            return_tuple = (False, None, None,
                            'user_exists_verified_locked')

        # if the user's email not verified, and account is locked, then this is
        # an newly created account but verification token might have expired
        elif usercheck[2] == False and usercheck[8] == True:
            return_tuple = (False, None, None,
                            'user_exists_notverified_locked')

        # if the user's email is verified and account is not locked, this is an
        # active account
        elif usercheck[2] == True and usercheck[8] == False:
            return_tuple = (False, None, None,
                            'user_exists_verified_notlocked')

        # if the user's email is not verified and the account is not locked,
        # something has gone wrong. we'll deny access
        elif usercheck[2] == True and usercheck[8] == False:
            return_tuple = (False, None, None,
                            'user_exists_notverified_notlocked')

        # if the user's email is not verified and the account is not locked,
        # something has gone wrong. we'll deny access
        else:
            return_tuple = (False, None, None, 'verification_failed')

    # if there is no existing user with this username, then we can generate a
    # verification token
    elif usercheck is False:

        verification_token = make_verification_token()
        verification_time = time.time()

        query = ('insert into lcsweb_user_verifications '
                 '(email, client_header, ip_address, '
                 'verification_token, sent_unixtime) '
                 ' values (%s, %s, %s,%s,%s)')

        query_params = (email, clientheader, ipaddress,
                        verification_token, verification_time)

        try:

            cur.execute(query, query_params)
            database.commit()

            return_tuple = (True, email,
                            verification_token, 'verification_logged')

        except Exception as e:

            if LOGGER:
                LOGGER.error('could not insert verification token into DB, '
                             'error was: %s' % e)
            if DEBUGMODE:
                print('could not insert verification token into DB, '
                      'error was: %s' % e)
            database.rollback()

            return_tuple = (False, None, None, 'database_error')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def check_verification_token(email,
                             verification_token,
                             max_expiry_minutes=15,
                             database=None):
    '''
    This checks if a verification token for the given email is valid and not
    expired. If expired or not matching, returns a failure. If successfully
    verified, returns a success and updates the verified column in the
    table. Uses the lcsweb_user_verifications table in the DB.

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


    # here we only check tokens that have not yet been verified because there's
    # no point in checking already verified tokens
    query = ('select email, verification_token, sent_unixtime from '
             'lcsweb_user_verifications where '
             '(email = %s) and (verification_token = %s) '
             'and (verified = FALSE)')

    query_params = (email,
                    verification_token)

    try:

        cur.execute(query, query_params)
        rows = cur.fetchall()

        if rows and len(rows) > 0:

            # compare the present time to the sent_time to see if the
            # verification token is still valid
            db_email, db_token, db_senttime = rows[0]

            if (time.time() - db_senttime) < 60.0*max_expiry_minutes:


                query = ('update lcsweb_user_verifications '
                         'set verified = TRUE '
                         'where email = %s and '
                         'verification_token = %s '
                         'returning verified')
                query_params = (email, verification_token)

                cur.execute(query, query_params)
                rows = cur.fetchall()
                database.commit()

                if rows and len(rows) > 0 and rows[0][0] == True:
                    return_tuple = (True, 'verification_successful')
                else:
                    return_tuple = (False, 'verification_failed_db_error')

            else:

                # if the verification token expired, remove it from the database
                query = ('delete from lcsweb_user_verifications '
                         'where email = %s and verification_token = %s')
                query_params = (email, verification_token)
                cur.execute(query, query_params)
                database.commit()

                return_tuple = (False, 'verification_token_expired')

        else:

            return_tuple = (False, 'verification_token_not_found')

    except Exception as e:

        if LOGGER:
            LOGGER.error('could not check verification token from DB, '
                         'error was: %s' % e)
        if DEBUGMODE:
            print('could not check verification token from DB, '
                  'error was: %s' % e)
        database.rollback()

        return_tuple = (False, 'database_error')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



####################
## LIMIT HANDLING ##
####################

def check_request_limits(username,
                         ipaddress=None,
                         apikey=False,
                         serviceperiod=86400.0,
                         database=None):
    '''
    This checks if the user has run over the request limits of their account.

    if apikey = False, checks for the user's session requests. If True, checks
    for apikey requests. Returns the total number of requests if available.

    serviceperiod defines the time interval in which to check for user
    requests. By default it's set at 24 hours (86400.0 seconds).

    Uses the lcsweb_authentication_request_history table in the DB.

    total requests for a user =
              (auth requests via logged-in user session for all services +
               auth requests via all apikeys for all services)

    '''

    if apikey == True:
        request_method = 'apikey'
    else:
        request_method = 'session'

    # the query method depends on if we're checking on the anonuser or not. the
    # anonuser can log in from any ip address, so we need to limit the
    # combination of the username 'anonuser' and the ip address they're
    # reporting. this is easily defeated, however, if the client just changes
    # their IP address.
    # FIXME: think about how to deal with this vulnerability!
    if username == 'anonuser' and ipaddress is not None:
        query = ("select count(*) from lcsweb_authentication_request_history "
                 "where username = %s and (%s - auth_unixtime) < %s and "
                 "ip_address = %s")
        query_params = (username, time.time(), serviceperiod, ipaddress)
    elif username != 'anonuser':
        query = ("select count(*) from lcsweb_authentication_request_history "
                 "where username = %s and (%s - auth_unixtime) < %s")
        query_params = (username, time.time(), serviceperiod)
    else:
        return (False, 'not_anonuser_and_no_ipaddress', None)

    # set up a database connection
    closedb = False

    if database:
        cur, handle = database.newcursor()
    else:
        database = lcdb.LCDB()
        database.open_default()
        cur, handle = database.newcursor()
        closedb = True

    try:
        cur.execute(query, query_params)
        rows = cur.fetchone()

        # if we have a result
        if rows and len(rows) > 0:

            nrequests = rows[0]

            # check the allowed requests for this user
            user_info = check_user(username,
                                   getinfo=True,
                                   database=database)

            # make sure this user exists
            if user_info is not False:

                user_allowed_requests = user_info[5]

                if user_allowed_requests == -1:
                    return_tuple = (True, 'auth_requests_below_limit',
                                    nrequests)
                elif nrequests <= user_allowed_requests:
                    return_tuple = (True, 'auth_requests_below_limit',
                                    nrequests)
                else:
                    return_tuple = (False, 'auth_requests_exceed_limit',
                                    nrequests)

            else:
                return_tuple = (False, 'user_not_found')
        else:
            # otherwise, we can't figure this out, 'fail deadly', i.e. deny this
            # request if there's an error
            return_tuple = (False,
                            'cannot_find_request_limits_database_error',
                            None)

    except Exception as e:

        if LOGGER:
            LOGGER.error('cannot check request limits '
                         'for user %s, error was: %s' % (username, e))
        if DEBUGMODE:
            print('cannot check request limits '
                  'for user %s, error was: %s' % (username, e))
        database.rollback()

        return_tuple = (False, 'database_error', None)


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple


###########################
## LOGGING AUTH REQUESTS ##
###########################

def log_auth_request(username,
                     requestmethod,
                     ipaddress,
                     clientheader,
                     service,
                     outcome,
                     message,
                     database=None):
    '''
    This logs the authentication request and its outcome to the
    lcsweb_authentication_request_history table.

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


    if username is None:
        username = 'anonuser'

    db_outcome = 'successful' if outcome == True else 'denied'

    query = ('insert into lcsweb_authentication_request_history '
             '(username, auth_unixtime, ip_address, '
             'client_header, request_method, service_requested, '
             'auth_outcome, auth_message) values '
             '(%s, %s, %s, %s, %s, %s, %s, %s)')

    query_params = (username, time.time(), ipaddress,
                    clientheader, requestmethod, service,
                    db_outcome, message)

    try:

        cur.execute(query, query_params)
        database.commit()

        returnval = (True, 'authrequest_logged_ok')

    except Exception as e:
        if LOGGER:
            LOGGER.error('failed to log authrequest, error was: %s' % e)
        if DEBUGMODE:
            print('failed to log authrequest, error was: %s' % e)
        database.rollback()

        returnval = (False, 'database_error')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return returnval



#####################################################
## ALL-IN-ONE AUTHENTICATION OF LC SERVER REQUESTS ##
#####################################################

def authenticate_and_log_request(ipaddress,
                                 clientheader,
                                 service,
                                 session_token=None,
                                 apikey=None,
                                 database=None):
    '''This does the following:

    1. checks if the user/apikey is allowed to use service
    2. check if the user/apikey is under its request limits
    3. logs the request auth outcome
    4. returns the authresult and authmessage

    Requires either session_token to be not None OR apikey to be not False.

    (session_token = <token>, apikey = False) -> interactive logged in user
    (session_token = None, apikey = False) -> interactive logged in user

    (session_token = None, apikey = None)     -> anonuser's apikey
    (session_token = None, apikey = <apikey>) -> noninteractive apikey

    anything else is not recognized

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

    # (session_token is None, apikey is not False) -> anonuser's / actual apikey
    # use the anonuser's apikey in the this case
    if session_token is None and apikey is not False:

        # if the apikey is none, use the anonuser's apikey
        if apikey is None:
            apikey_info = check_apikey('anonuser',
                                       service,
                                       getinfo=True,
                                       database=database)

            if apikey_info is not False:
                apikey_user = 'anonuser'
                apikey_to_use = apikey_info[1]
                apikey_check = True
            else:
                apikey_check = False

        # if the apikey is not None, get the user associated with it
        else:

            apikey_user = get_apikey_username(apikey, database=database)

            if apikey_user is not False:
                apikey_to_use = apikey
                apikey_check = True

            else:
                apikey_check = False

        if apikey_check is True:

            # first, check if the apikey can use this service
            auth_result = authenticate_apikey(service,
                                              ipaddress,
                                              clientheader,
                                              username=apikey_user,
                                              apikey=apikey_to_use,
                                              database=database,
                                              logrequest=True)

            # then, if the auth_result is True, check the apikey's limits
            if auth_result[0] is True:

                limits_result = check_request_limits(apikey_user,
                                                     ipaddress=ipaddress,
                                                     apikey=True,
                                                     database=database)

                if limits_result[0] is True:
                    return_tuple = (True,
                                    apikey_user,
                                    'apikey',
                                    apikey_to_use,
                                    'request_approved')
                else:
                    return_tuple = (False,
                                    apikey_user,
                                    'apikey',
                                    apikey_to_use,
                                    'request_over_limit')

            else:

                return_tuple = (False,
                                apikey_user,
                                'apikey',
                                apikey_to_use,
                                'request_auth_failed')

        else:

            return_tuple = (False,
                            None,
                            'apikey',
                            None,
                            'request_unknown_user')


    # (session_token is None/<token>, apikey is False) -> interactive user
    # use the session_token in this case
    elif apikey is False:

        # if there is actually a session token available, use it
        if session_token is not None:

            # check the user associated with this session_token
            session_info = check_user_session(session_token, getinfo=True,
                                              onlyactive=True, database=database)

        # if there isn't one, then we're in anonymous mode, generate session
        # token by ourselves and make the session_info tuple. the function that
        # calls this function is then responsible for closing the anonymous
        # session later
        else:

            # insert the session into the lcsweb_user_sessions table
            anon_session_ok, anon_session_token = initiate_anonymous_session(
                ipaddress,
                clientheader,
                database=database
                )

            if anon_session_ok:
                session_info = (anon_session_token, 'anonuser')
            else:
                session_info = False


        # now proceed with the rest of the processing by checking if
        # session_info is not False (i.e. no session exists)
        if session_info is not False:

            session_token_to_use = session_info[0]
            session_user = session_info[1]
            session_check = True

        else:

            session_check = False


        if session_check is True:

            # first, check if this user can use this service
            auth_result = auth_service_user(session_token_to_use,
                                            ipaddress,
                                            clientheader,
                                            service,
                                            logrequest=True,
                                            database=database)

            # if the user is allowed to use the service, check their limits
            if auth_result[0] is True:

                limits_result = check_request_limits(session_user,
                                                     ipaddress=ipaddress,
                                                     apikey=False,
                                                     database=database)

                if limits_result[0] is True:
                    return_tuple = (True,
                                    session_user,
                                    'session',
                                    session_token_to_use,
                                    'request_approved')
                else:
                    return_tuple = (False,
                                    session_user,
                                    'session',
                                    session_token_to_use,
                                    'request_over_limit')
            else:

                return_tuple = (False,
                                session_user,
                                'session',
                                session_token_to_use,
                                'request_auth_failed')

        else:

            return_tuple = (False,
                            None,
                            'session',
                            None,
                            'request_unknown_user')


    # all other combinations of session_token and apikey are invalid
    else:

        return_tuple = (False, None, None, None, 'request_type_unknown')


    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple


######################
## LOGGING OF TASKS ##
######################

def log_task_start(username,
                   ipaddress,
                   clientheader,
                   task_type,
                   request_method,
                   request_unixtime,
                   request_args,
                   request_kwargs,
                   database=None):
    '''
    This logs starting stats about the task to the lcsweb_user_task_history
    table. All items from this table are shown at /lightcurves/users/home/tasks
    and the latest 5 items are shown in the tasks panel at
    /lightcurves/users/home.

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

    query = ('insert into lcsweb_user_task_history '
             '(username, ipaddress, clientheader,'
             'tasktype, taskmethod, unixtime, taskargs, taskkwargs) values '
             '(%s, %s, %s, %s, %s, %s, %s, %s) returning taskid')

    query_params = (username, ipaddress, clientheader,
                    task_type, request_method, request_unixtime,
                    dumps(request_args,ensure_ascii=True),
                    dumps(request_kwargs,ensure_ascii=True))

    try:

        cur.execute(query, query_params)
        rows = cur.fetchone()
        database.commit()

        if len(rows) > 0:
            return_tuple = (True, rows[0], 'task_start_logged')
        else:
            return_tuple = (False, -1, 'task_start_logging_failed')

    except Exception as e:

        LOGGER.exception('failed to log task request start, '
                         'error was: %s' % e)
        database.rollback()
        return_tuple = (False, -1, 'task_start_logging_failed')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple



def log_task_end(
    task_id, # use this as primary key to find task
    task_starttime, # if that fails, use this + next 3 to find task
    task_username,
    task_ipaddress,
    task_clientheader,
    task_ok,
    task_message,
    task_elapsedtime,
    task_query,
    database=None
    ):
    '''
    This logs final stats about the task to the lcsweb_user_task_history
    table. All items from this table are shown at /lightcurves/users/home/tasks
    and the latest 5 items are shown in the tasks panel at
    /lightcurves/users/home.

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

    if task_id > -1:
        query = ('update lcsweb_user_task_history set '
                 'returnok = %s, returnmsg = %s, elapsedtime = %s, '
                 'sqlquery = %s where taskid = %s returning taskid')
        query_params = (task_ok, task_message, task_elapsedtime,
                        task_query, task_id)
    else:
        query = ('update lcsweb_user_task_history set '
                 'returnok = %s, returnmsg = %s, elapsedtime = %s, '
                 'sqlquery = %s where unixtime = %s and '
                 'username = %s and ipaddress = %s and clientheader = %s '
                 'returning taskid')
        query_params = (task_ok, task_message, task_elapsedtime,
                        task_query, task_starttime, task_username,
                        task_ipaddress, task_clientheader)

    try:

        cur.execute(query, query_params)
        rows = cur.fetchone()
        database.commit()

        if len(rows) > 0:
            return_tuple = (True, rows[0], 'task_end_logged')
        else:
            return_tuple = (False, -1, 'task_end_logging_failed')

    except Exception as e:

        LOGGER.exception('failed to log task request, error was: %s' % e)
        database.rollback()
        return_tuple = (False, -1, 'task_end_loggging_failed')

    # close the database at the end if we have to
    database.close_cursor(handle)
    if closedb:
        database.close_connection()

    return return_tuple
