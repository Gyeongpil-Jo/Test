from django.contrib import admin
from django.urls import path, include


urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('apps.home.urls')),
    path('graphene/', include('apps.graphene.urls')),
    path('dna/', include('apps.dna.urls')),
]
