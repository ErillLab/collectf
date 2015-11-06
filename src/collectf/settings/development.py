from .base import *             # NOQA

# Database
# https://docs.djangoproject.com/en/dev/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'collectfdb_devel',
        'USER': 'sefa',
        'PASSWORD': '46544654',
        'HOST': 'localhost',
        'PORT': ''
    }
}

# INSTALLED_APPS += ['debug_toolbar']
DEBUG_TOOLBAR_PANELS = [
    'debug_toolbar.panels.sql.SQLPanel',
    'debug_toolbar.panels.profiling.ProfilingPanel'
]
