from django.contrib import admin
from django.urls import path, include


urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('apps.home.urls')),
    path('graphene/', include('apps.graphene.urls')),
    path('SA/', include('apps.SA.urls')),
    path('dna/', include('apps.dna.urls')),
]
