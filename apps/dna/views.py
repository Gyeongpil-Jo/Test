from django.shortcuts import render, get_object_or_404, redirect
from django.http import FileResponse, JsonResponse
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponseNotAllowed, HttpResponse
from django.core.paginator import Paginator
from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.core import serializers

from .forms import DnaForm
from .models import Dna
from Test.settings.base import BASE_DIR, MEDIA_ROOT

import subprocess
import tempfile
import tarfile
import json
import os


def input_dna(request):
    if request.method == 'POST':
        form = DnaForm(request.POST)

        if form.is_valid():
            job = form.save(commit=False)
            job.seq_num = len(job.seq)
            job.save()

            return redirect('dna:output', job_id=job.id)

    else:
        form = DnaForm()
    context = {'form': form}
    return render(request, 'dna/input.html', context)


def output_dna(request, job_id):
    job = get_object_or_404(Dna, pk=job_id)
    context = {'job': job}
    return render(request, 'dna/output.html', context)


def build_dna(request, job_id):
    job = get_object_or_404(Dna, pk=job_id)
    program_dir = BASE_DIR + "/dna"

    with tempfile.TemporaryDirectory() as temp_dir:
        os.system(f"cp -r {program_dir}/* {temp_dir}")
        proc_result = subprocess.run(
            args=[f"{temp_dir}/loop.sh {job.seq} {job.seq_num}"],
            stdout=subprocess.PIPE,
            cwd=temp_dir,
            text=True,
            shell=True,
        )

        if proc_result.returncode == 0:
            with tarfile.open(BASE_DIR / 'media' / 'jobs_dna' / f'{job_id:06d}.tar.gz', 'w:gz') as tar:
                cur_dir = os.getcwd()
                os.chdir(temp_dir)
                for f in ('dna.pdb', 'dna.itp', 'program_1.nab', 'template.nab', 'loop.sh'):
                    tar.add(f, recursive=False)
                os.chdir(cur_dir)
        else:
            pass

def file_download(request, job_id):
    file_path = os.path.join(MEDIA_ROOT, 'jobs_dna')
    fs = FileSystemStorage(file_path)
    response = FileResponse(fs.open(f'{job_id:06d}.tar.gz', 'rb'))
    response['Content-Disposition'] = 'attachment; filename={}'.format("dna.tar.gz")

    return response
