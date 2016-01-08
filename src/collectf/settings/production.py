from .base import *             # NOQA

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'collectfdb',
        'USER': 'erilllab',
        'PASSWORD': '',
        'HOST': 'localhost',
        'PORT': ''
    }
}

DEBUG = False

ALLOWED_HOSTS = ['.collectf-v2.umbc.edu',
                 '.collectf.org',
                 '.collectf.umbc.edu']

# Gather static files here.
STATIC_ROOT = os.path.join(BASE_DIR, 'collectstatic')

SITE_ID = 1

# Silence "Invalid HTTP_HOST header" messages
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'null': {
            'class': 'logging.NullHandler',
        },
    },
    'loggers': {
        'django.security.DisallowedHost': {
            'handlers': ['null'],
            'propagate': False,
        },
    },
}

WEBLOGO_BIN = '/usr/bin/weblogo'

