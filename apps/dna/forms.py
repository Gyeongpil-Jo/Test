
from django import forms
from django.core.exceptions import NON_FIELD_ERRORS
from apps.dna.models import Dna


class DnaForm(forms.ModelForm):
    class Meta:
        model = Dna
        fields = ['seq', 'seq_num']
