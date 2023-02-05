
from django import forms
from django.core.exceptions import NON_FIELD_ERRORS
from apps.dna.models import Dna


class DnaForm(forms.ModelForm):
    class Meta:
        model = Dna
        fields = ['seq', 'seq_num']
        error_messages = {
            'seq': {
                'blank': 'Your inputs must be A, T, G, C DNA sequence'
            }
        }

    def clean_seq(self):
        seq = self.cleaned_data.get('seq')
        seq_tmp = seq.upper()

        for s in ('A', 'T', 'G', 'C'):
            seq_tmp = seq_tmp.replace(s, '')

        if len(seq_tmp) != 0:
            seq = ''

        return seq
