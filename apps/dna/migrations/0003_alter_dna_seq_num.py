# Generated by Django 4.1 on 2023-02-02 12:20

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dna', '0002_dna_seq_num'),
    ]

    operations = [
        migrations.AlterField(
            model_name='dna',
            name='seq_num',
            field=models.TextField(default=''),
        ),
    ]
