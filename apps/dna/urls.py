from django.urls import path
from . import views


app_name = 'dna'


urlpatterns = [
    path('', views.input_dna, name='input'),
    path('<int:job_id>/', views.output_dna, name='output'),
    path('file_download/<int:job_id>/', views.file_download, name='file_download'),
]
