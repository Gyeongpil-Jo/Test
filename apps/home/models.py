from django.db import models
from django.contrib.auth.models import User


class Generating(models.Model):
    file = models.FileField(upload_to='files/', blank=False)
    type = models.TextField()
    create_date = models.DateTimeField()