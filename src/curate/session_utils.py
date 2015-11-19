"""
Session utils
Python module to use Django's session utility in curation forms.
"""

PREFIX = 'curation_'


def put(session, key, val):
    """Puts (key,val) to session data."""
    session[PREFIX + key] = val
    session.modified = True


def get(session, key):
    """Reads from session data."""
    return session.get(PREFIX + key)


def remove(session, key):
    """Deletes from session data."""
    if PREFIX + key in session:
        del session[PREFIX + key]


def has(session, key):
    """Checks if the key is in session."""
    return PREFIX + key in session


def clear(session):
    """Clears session object."""
    for key, val in session.items():
        if key.startswith(PREFIX):
            del session[key]


def session_print(session):
    """Prints the session."""
    sep = '*' * 80
    print sep
    print 'session data:'
    for k, v in session.items():
        print k, v
    print sep
