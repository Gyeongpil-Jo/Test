# Generated by Django 4.1 on 2022-12-23 08:29

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('SA', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='MultipleSequence',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fasta_file', models.FileField(upload_to='fasta/')),
            ],
        ),
    ]