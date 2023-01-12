from django.shortcuts import render, redirect, get_object_or_404
from .models import Sequence, MultipleSequence
from .forms import SequenceForm, MultipleSequenceForm

from Bio import pairwise2, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.pairwise2 import format_alignment


def input_sa(request):
    form = SequenceForm()
    multi_form = MultipleSequenceForm()
    if request.method == 'POST':
        if request.POST.get('form_type') == 'double':
            double_form = SequenceForm(request.POST)

            if double_form.is_valid():
                job = double_form.save(commit=False)
                job.save()

                return redirect('SA:output', job_id=job.id)

        if request.POST.get('form_type') == 'multi':
            multi_form = MultipleSequenceForm(request.POST, request.FILES)

            if multi_form.is_valid():
                job = multi_form.save(commit=False)
                job.save()

                return redirect('SA:output_multi', job_id=job.id)

    else:
        form = SequenceForm()
        multi_form = MultipleSequenceForm()

    context = {'form': form, 'multi_form': multi_form}

    return render(request, 'SA/input.html', context)


def output_sa(request, job_id):
    job = get_object_or_404(Sequence, pk=job_id)
    result = sequence_alignment(job.S_1, job.S_2)
    context = {'result': result}

    return render(request, 'SA/output.html', context)


def sequence_alignment(s1, s2):
    seq1 = Seq(s1)
    seq2 = Seq(s2)
    alignment = pairwise2.align.globalxx(seq1, seq2)

    result = []
    for alignments in alignment:
        result.append(format_alignment(*alignments))

    return result


def output_multi_sa(request, job_id):
    result = multiple_sequence_alignment(job_id)
    context = {'result': result}

    return render(request, 'SA/output_multi.html', context)


def multiple_sequence_alignment(job_id):
    job = get_object_or_404(MultipleSequence, pk=job_id)
    with open(job.fasta_file.path, 'r') as fasta:
        alignment = AlignIO.read(fasta, 'fasta')
        print(alignment)
    return alignment


