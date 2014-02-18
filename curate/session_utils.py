"""
Session utils

Python module to use Django's session utility.

"""

def put(session, key, val):
    """Put (key,val) to session data."""
    session[key] = val
    session.modified = True

def get(session, key):
    """Read from session data."""
    return session.get(key)

def remove(session, key):
    """Delete from session data."""
    del session[key]

def has(session, key):
    """Check if the key is in session."""
    return key in session

def clear(session):
    """Clear session object."""
    for key,val in session.items():
        if not key.startswith('_'):
            del session[key]
