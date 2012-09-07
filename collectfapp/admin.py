from collectfapp.models import *
from django.contrib import admin
from django.db.models import get_models
from django.db.models import get_app

# curator admin view
class CuratorAdmin(admin.ModelAdmin):
    filter_horizontal = ("assigned_papers",)

admin.site.register(Curator, CuratorAdmin) # register Curator
# register rest of models
app = get_app("collectfapp")
for model in get_models(app):
    try:
        admin.site.register(model)
    except admin.sites.AlreadyRegistered:
        pass

