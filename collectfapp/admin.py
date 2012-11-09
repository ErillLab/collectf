from collectfapp.models import *
from django.contrib import admin
from django.db.models import get_models
from django.db.models import get_app

# curation admin view
class CurationAdmin(admin.ModelAdmin):
    filter_horizontal = ("experimental_techniques",)
admin.site.register(Curation, CurationAdmin) # register Curation

# register rest of models
app = get_app("collectfapp")
for model in get_models(app):
    try:
        admin.site.register(model)
    except admin.sites.AlreadyRegistered:
        pass


# unregister Regulation
# don't show any field. Gene field in the Regulation model
# ruins everything. As default behavior, Django tries to
# generate a dropdown box with __ALL__ genes in the database,
# which kills the server.
admin.site.unregister(Regulation)
