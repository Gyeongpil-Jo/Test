
from django.db import models


class Dna(models.Model):
    seq = models.TextField(default='')
    seq_num = models.IntegerField(default='', blank=True)
