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
