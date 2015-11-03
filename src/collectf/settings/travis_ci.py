from .base import *             # NOQA

# Database
# https://docs.djangoproject.com/en/dev/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'collectfdb_devel',
        'USER': 'travis',
        'PASSWORD': '',
        'HOST': 'localhost',
        'PORT': ''
    }
}
