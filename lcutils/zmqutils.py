#!/usr/bin/env python

'''
zmqutils.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - May 2013

This file contains some utility functions useful for message passing over the
zeroMQ message bus.

'''

from datetime import datetime
from json import dumps, loads
from time import time
from socket import gethostname
from hashlib import sha256

CURRENT_HOSTNAME = gethostname()

def to_zmq_msg(msg_to,
               msg_from,
               msg_subject,
               msg_dict,
               resptoken=False):
    '''
    This takes strings for the message sender, recipient, and subject and a dict
    for the content of a message to be sent over the ZMQ wire, and turns them
    into a string format that allows filtering by the ZMQ sockets.

    WARNING: All numbers are turned into strings, so look out for precision
    loss.

    The message format is multipart (suitable for
    zmq.context.socket.send_multipart):

    ['<msg_to>|<msg_subject>','<msg_from>','{<msg_dict JSON>}']

    NOTE: if the msg_dict does not contain a unix_time field, one will be added
    automatically and populated with the current time

    if resptoken = True, a reponse token will be generated based on the contents
    of the msg_dict pre-formatting by this function, and will be put into the
    msg_dict.

    '''

    # if we have to add a response token
    # add the unix_time key to the msg_dict
    if 'unix_time' not in msg_dict:
        msg_dict['unix_time'] = '%.3f' % time()

    # add the hostname key to the msg_dict
    if 'sender_hostname' not in msg_dict:
        msg_dict['sender_hostname'] = CURRENT_HOSTNAME

    # generate a response token if requested. we do this after inserting the
    # time so each response token generated is guaranteed to be unique
    if resptoken:
        json_msgdict = dumps(msg_dict, ensure_ascii=True)
        tokenhash = sha256(json_msgdict).hexdigest()
        msg_dict['response_token'] = tokenhash

    msg_key = '%s|%s' % (msg_to, msg_subject)

    msg_contents = dumps(msg_dict,ensure_ascii=True)
    msg = [str(msg_key), str(msg_from), str(msg_contents)]

    if resptoken:
        return (msg, tokenhash)
    else:
        return msg



def from_zmq_msg(msg):
    '''
    This takes a multipart msg sent over the ZMQ wire, and turns it into a dict
    containing the fields, and the address metadata.

    '''

    msg_key, msg_sender, msg_contents = msg

    msg_recipient, msg_subject = msg_key.split('|')

    msg_sender = msg_sender.strip()
    msg_recipient = msg_recipient.strip()
    msg_subject = msg_subject.strip()
    msg_dict = loads(msg_contents)

    return msg_recipient, msg_sender, msg_subject, msg_dict



def datetime_to_jsondate(date_time, dateonly=False, timeonly=False):
    '''
    This turns a datetime object into a string so that it can be used in a JSON
    message. Note: datetime is expected in UTC.

    '''

    if dateonly and not timeonly:
        return date_time.strftime('%Y-%m-%d')
    elif timeonly and not dateonly:
        return date_time.strftime('%H:%M:%S.%f')
    elif not timeonly and not dateonly:
        return date_time.strftime('%Y-%m-%d %H:%M:%S.%f')
    else:
        return date_time.strftime('%Y-%m-%d %H:%M:%S.%f')


def jsondate_to_datetime(json_date):
    '''
    This turns a JSON date string into a datetime object. Note: json_date is
    expected in UTC.

    '''
    return datetime.strptime(json_date,'%Y-%m-%d %H:%M:%S.%f')
