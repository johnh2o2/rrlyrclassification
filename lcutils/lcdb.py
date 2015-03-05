#!/usr/bin/env python

'''
lcdb.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - May 2013

Serves as a lightweight PostgreSQL DB interface for other modules in this
project.

'''

import logging
import ConfigParser
import os
import hashlib

import psycopg2 as pg
import psycopg2.extras

# setup a logger
LOGGER = logging.getLogger('lcdb')
LOGGER.addHandler(logging.NullHandler())

# parse the configuration file to get the default database credentials
'''
CONF_FILE = 'lcserver.conf'

CONF = ConfigParser.ConfigParser()
CONF.read(CONF_FILE)

# database config
DBUSER = CONF.get('database','user')
DBPASS = CONF.get('database','password')
DBDATA = CONF.get('database','database')
DBHOST = CONF.get('database','host')
'''

class LCDB(object):
    '''
    This is an object serving as an interface to the postgres DB. It implements
    the following methods:

    LCDB.open(database,
              user,
              password) -> open a new connection to postgres using the provided
                           credentials

    LCDB.cursor(handle,
                dictcursor=<True/False>) -> return a postgres DB cursor with the
                                            handle specified. if the handle does
                                            not exist, a new one will be created
                                            for the cursor. if dictcursor is
                                            True, will return a cursor allowing
                                            addressing columns as a dictionary

    LCDB.commit() -> commits any pending transactions to the DB

    LCDB.close_cursor(handle) -> closes cursor specified in handle

    LCDB.close_connection() -> close all cursors currently active, and then
                               close the connection to the database

    LCDB's main purpose is to avoid creating new postgres connections throughout
    the course of lcserver's work; these are relatively expensive. Instead, we
    get new cursors when needed, and then pass these around as needed. The
    connection also remains open for the whole lifetime of the lcserver's
    runtime, keeping things simple.

    '''

    def __init__(self,
                 database=None,
                 user=None,
                 password=None,
                 host=None):

        self.connection = None
        self.user = None
        self.database = None
        self.host = None
        self.cursors = {}

        if database and user and password and host:
            self.open(database, user, password, host)



    def open(self, database, user, password, host):
        '''
        This just creates the database connection and stores it in
        self.connection.

        '''
        try:

            self.connection = pg.connect(user=user,
                                         password=password,
                                         database=database,
                                         host=host)
            if LOGGER:
                LOGGER.debug('postgres connection successfully '
                            'created, using DB %s, user %s' % (database,
                                                               user))
            else:
                print('postgres connection successfully '
                      'created, using DB %s, user %s' % (database,
                                                         user))
            self.database = database
            self.user = user

        except:

            if LOGGER:
                LOGGER.debug('postgres connection failed, '
                            'using DB %s, user %s' % (database,
                                                      user))
            else:
                print('postgres connection failed, '
                      'using DB %s, user %s' % (database, user))
                self.database = None
                self.user = None



    def open_default(self):
        '''
        This opens the database connection using the default database parameters
        given in the lcserver.conf file.

        '''

        self.open(DBDATA, DBUSER, DBPASS, DBHOST)


    def autocommit(self):
        '''
        This sets the database connection to autocommit. Must be called before
        any cursors have been instantiated.

        '''

        if len(self.cursors.keys()) == 0:
            self.connection.autocommit = True
        else:
            raise AttributeError('database cursors are already active, '
                                 'cannot switch to autocommit now')


    def cursor(self, handle, dictcursor=False):
        '''
        This gets or creates a DB cursor for the current DB connection.

        dictcursor = True -> use a cursor where each returned row can be
                             addressed as a dictionary by column name

        '''

        if handle in self.cursors:

            return self.cursors[handle]

        else:
            if dictcursor:
                self.cursors[handle] = self.connection.cursor(
                    cursor_factory=psycopg2.extras.DictCursor
                    )
            else:
                self.cursors[handle] = self.connection.cursor()

            return self.cursors[handle]


    def newcursor(self, dictcursor=False):
        '''
        This creates a DB cursor for the current DB connection using a
        randomly generated handle. Returns a tuple with cursor and handle.

        dictcursor = True -> use a cursor where each returned row can be
                             addressed as a dictionary by column name

        '''

        handle = hashlib.sha256(os.urandom(12)).hexdigest()

        if dictcursor:
            self.cursors[handle] = self.connection.cursor(
                cursor_factory=psycopg2.extras.DictCursor
                )
        else:
            self.cursors[handle] = self.connection.cursor()

            return (self.cursors[handle], handle)



    def commit(self):
        '''
        This just calls the connection's commit method.

        '''
        if not self.connection.closed:
            self.connection.commit()
        else:
            raise AttributeError('postgres connection to %s is closed' %
                                 self.database)


    def rollback(self):
        '''
        This just calls the connection's commit method.

        '''
        if not self.connection.closed:
            self.connection.rollback()
        else:
            raise AttributeError('postgres connection to %s is closed' %
                                 self.database)



    def close_cursor(self, handle):
        '''
        Closes the cursor specified and removes it from the self.cursors
        dictionary.

        '''

        if handle in self.cursors:
            self.cursors[handle].close()
        else:
            raise KeyError('cursor with handle %s was not found' % handle)



    def close_connection(self):
        '''
        This closes all cursors currently in use, and then closes the DB
        connection.

        '''

        self.connection.close()
        if LOGGER:
            LOGGER.debug('postgres connection closed for DB %s' % self.database)
        else:
            print('postgres connection closed for DB %s' % self.database)
