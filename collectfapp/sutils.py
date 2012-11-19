# Session Utils
# Python module to exploit Django's session utility

def sput(session, key, val):
    # put (key,val) to session data
    session[key] = val
    session.modified = True

def sget(session, key):
    # read from session data
    return session[key]

def sdel(session, key):
    # delete from session data
    del session[key]

def sin(session, key):
    return key in session


def clear(session):
    for key,val in session.items():
        if not key.startswith('_'):
            del session[key]

