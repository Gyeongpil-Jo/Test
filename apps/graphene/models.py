from django.db import models


class Graphene(models.Model):
    size_x = models.TextField(default='')
    size_y = models.TextField(default='')
    pbc = models.TextField(default='None', blank=True)
    cnt = models.TextField(default='None', blank=True)
    hole = models.TextField(default='None')
    radius = models.TextField()
