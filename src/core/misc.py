"""These are some miscellaneous functions that are used here and there."""


def is_float(s):
    """Given a string s, return if it is a floating point number or not."""
    try:
        float(s)
        return True
    except ValueError:
        return False
