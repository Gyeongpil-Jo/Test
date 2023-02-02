
from django.db import models


class Dna(models.Model):
    seq = models.TextField(default='')
    seq_num = models.TextField(default='', blank=True)
