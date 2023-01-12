from django.shortcuts import render, get_object_or_404, redirect
from django.http import HttpResponseNotAllowed, HttpResponse
from django.utils import timezone
from django.core.paginator import Paginator
from django.contrib.auth.decorators import login_required
from .forms import GeneratingForm
from django.contrib import messages
from django.core import serializers
from .models import Generating
import subprocess


def generate(request, generating):
    with open(generating.file.path, 'r') as f:
        text = f.readlines()
    context = {'result': generating.type, 'text': text}
    return render(request, 'home/generate.html', context)


def main(request):
    if request.method == 'POST':
        form = GeneratingForm(request.POST, request.FILES)
        if form.is_valid():
            generating = form.save(commit=False)
            generating.create_date = timezone.now()
            generating.save()
            return generate(request, generating)

    else:
        form = GeneratingForm()
    context = {'form': form}
    return render(request, 'home/main.html', context)




