from django import forms
from apps.home.models import Generating


class GeneratingForm(forms.ModelForm):
    class Meta:
        model = Generating
        fields = ['file', 'type']
        labels = {
            'file': 'file',
            'type': 'type',
        }


class DateFormatString(forms.Form):
    format_string = forms.CharField()
