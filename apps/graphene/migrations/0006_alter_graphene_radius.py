# Generated by Django 4.1 on 2022-12-22 08:32

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('graphene', '0005_alter_graphene_radius'),
    ]

    operations = [
        migrations.AlterField(
            model_name='graphene',
            name='radius',
            field=models.TextField(),
        ),
    ]
