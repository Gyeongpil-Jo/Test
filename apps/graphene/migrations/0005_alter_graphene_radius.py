# Generated by Django 4.1 on 2022-12-21 02:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('graphene', '0004_graphene_radius_alter_graphene_cnt_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='graphene',
            name='radius',
            field=models.TextField(blank=True, default=''),
        ),
    ]
