from django.urls import path
from . import views


app_name = 'graphene'


urlpatterns = [
    path('', views.input_graphene, name='input'),
    path('<int:job_id>/', views.output_graphene, name='output'),
    path('file_download/<int:job_id>/', views.file_download, name='file_download'),
]
