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
