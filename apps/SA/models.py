from django.db import models


class Sequence(models.Model):
    S_1 = models.TextField()
    S_2 = models.TextField()


class MultipleSequence(models.Model):
    fasta_file = models.FileField(upload_to='fasta/', blank=False)