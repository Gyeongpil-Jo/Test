from django.shortcuts import render, get_object_or_404, redirect
from django.http import FileResponse
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponseNotAllowed, HttpResponse
from django.core.paginator import Paginator
from django.contrib.auth.decorators import login_required
from django.contrib import messages

from .forms import GrapheneForm
from .models import Graphene
from Test.settings import BASE_DIR, MEDIA_ROOT

import subprocess
import tempfile
import tarfile
import os


def input_graphene(request):
    if request.method == 'POST':
        form = GrapheneForm(request.POST)

        if form.is_valid():
            job = form.save(commit=False)
            job.save()
            build_graphene(request, job_id=job.id)

            return redirect('graphene:output', job_id=job.id)

    else:
        form = GrapheneForm()
    context = {'form': form}
    return render(request, 'graphene/input.html', context)


def output_graphene(request, job_id):
    job = get_object_or_404(Graphene, pk=job_id)
    context = {'job': job}
    return render(request, 'graphene/output.html', context)


def build_graphene(request, job_id):
    job = get_object_or_404(Graphene, pk=job_id)

    args = ['/usr/bin/perl', BASE_DIR / 'build_graphene.pl',
            f"--size={job.size_x},{job.size_y}"]

    if job.pbc != 'None':
        args.append(f"--pbc={job.pbc}")

    if job.cnt != 'None':
        args.append(f"--cnt={job.cnt}")

    if job.hole != 'None':
        args.append(f"--hole={job.radius}")

    with tempfile.TemporaryDirectory() as temp_dir:
        proc_result = subprocess.run(
            args=args,
            cwd=temp_dir,
            stdout=subprocess.PIPE,
            text=True,
        )
        if proc_result.returncode == 0:
            with tarfile.open(BASE_DIR / 'media' / 'jobs' / f'{job.id:06d}.tar.gz', 'w:gz') as tar:
                cur_dir = os.getcwd()
                os.chdir(temp_dir)
                for f in ('graphene.pdb', 'graphene.itp', 'graphene.posres.itp'):
                    tar.add(f, recursive=False)
                os.chdir(cur_dir)
        else:
            pass


def file_download(request, job_id):
    file_path = os.path.join(MEDIA_ROOT, 'jobs')
    fs = FileSystemStorage(file_path)
    response = FileResponse(fs.open(f'{job_id:06d}.tar.gz', 'rb'))
    response['Content-Disposition'] = 'attachment; filename={}'.format("graphene.tar.gz")

    return response
