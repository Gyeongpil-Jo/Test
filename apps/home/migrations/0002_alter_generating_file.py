# Generated by Django 4.1 on 2022-11-14 07:35

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('home', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='generating',
            name='file',
            field=models.FileField(upload_to='files/'),
        ),
    ]