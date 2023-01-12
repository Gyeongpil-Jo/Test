from django import forms
from django.core.exceptions import NON_FIELD_ERRORS
from apps.SA.models import Sequence, MultipleSequence


class SequenceForm(forms.ModelForm):
    class Meta:
        model = Sequence
        fields = ['S_1', 'S_2']
        labels = {
            'S_1': 'Sequence 1',
            'S_2': 'Sequence 2',
        }


class MultipleSequenceForm(forms.ModelForm):
    class Meta:
        model = MultipleSequence
        fields = ['fasta_file']

