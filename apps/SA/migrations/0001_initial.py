# Generated by Django 4.1 on 2022-12-22 08:37

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('S_1', models.TextField()),
                ('S_2', models.TextField()),
            ],
        ),
    ]
