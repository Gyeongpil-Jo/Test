# Generated by Django 4.1 on 2022-11-07 11:28

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('graphene', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='graphene',
            name='size',
        ),
        migrations.AddField(
            model_name='graphene',
            name='size_x',
            field=models.TextField(default=''),
        ),
        migrations.AddField(
            model_name='graphene',
            name='size_y',
            field=models.TextField(default=''),
        ),
    ]
