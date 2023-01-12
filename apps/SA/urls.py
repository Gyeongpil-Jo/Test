from django.urls import path
from . import views


app_name = 'SA'


urlpatterns = [
    path('', views.input_sa, name='input'),
    path('<int:job_id>/', views.output_sa, name='output'),
    path('<int:job_id>/multi/', views.output_multi_sa, name='output_multi'),
]
