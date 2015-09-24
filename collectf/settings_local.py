"""
Django settings for collectf project.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.6/ref/settings/
"""
import os

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.6/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'egl#16o5n_=$j6()aqg47dahi2y%*li&fgrb63krv5_xnxy!k%'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

TEMPLATE_DIRS = (
    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'base', 'templates'),
    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'homepage', 'templates'),
    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'browse', 'templates'),
    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'curate', 'templates'),
    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'ncbi', 'templates'),
)

ALLOWED_HOSTS = []

# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    # CollecTF apps
    'base',
    'homepage',
    'browse',
    'curate',
    # dependencies
    'django_extensions',
    'south',
    'bootstrapform',
    'registration',
#    'debug_toolbar',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'collectf.urls'

WSGI_APPLICATION = 'collectf.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'STORAGE_ENGINE': 'MyISAM',
        'NAME': 'collectfdb_test',
        'USER': 'sefa',
        'PASSWORD': '46544654',
    }
}

# Change database username/password if build on Travis CI.
if os.getenv('BUILD_ON_TRAVIS', None):
    DATABASES['default']['USER'] = 'travis'
    DATABASES['default']['PASSWORD'] = ''

# Internationalization
# https://docs.djangoproject.com/en/1.6/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.6/howto/static-files/

STATIC_URL = '/static/'

STATICFILES_DIRS = (
    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'browse', 'static'),
)

STATIC_ROOT = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sitestatic')

# Redirect to here when logged in
LOGIN_REDIRECT_URL = '/'

# Make session-serializer use pickle, instead of JSON
SESSION_SERIALIZER = 'django.contrib.sessions.serializers.PickleSerializer'

# Absolute path to the directory of pickle files that are used.
# Don't put anything in this directory yourself.
PICKLE_ROOT = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pickles")


# Curation submission settings
NUMBER_OF_GENOME_ACCESSION_FIELDS = 8
NUMBER_OF_TF_ACCESSION_FIELDS = 5
NUMBER_OF_EXTERNAL_DATABASE_FIELDS = 5

# required for django-registration
ACCOUNT_ACTIVATION_DAYS = 7

# Profiling settings
DEBUG_TOOLBAR_PANELS = [
    #'debug_toolbar.panels.versions.VersionsPanel',
    'debug_toolbar.panels.timer.TimerPanel',
    #'debug_toolbar.panels.settings.SettingsPanel',
    #'debug_toolbar.panels.headers.HeadersPanel',
    #'debug_toolbar.panels.request.RequestPanel',
    'debug_toolbar.panels.sql.SQLPanel',
   # 'debug_toolbar.panels.staticfiles.StaticFilesPanel',
   # 'debug_toolbar.panels.templates.TemplatesPanel',
   # 'debug_toolbar.panels.cache.CachePanel',
   # 'debug_toolbar.panels.signals.SignalsPanel',
   # 'debug_toolbar.panels.logging.LoggingPanel',
   # 'debug_toolbar.panels.redirects.RedirectsPanel',
    'debug_toolbar.panels.profiling.ProfilingPanel',
]

