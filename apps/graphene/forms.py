from django import forms
from django.core.exceptions import NON_FIELD_ERRORS
from apps.graphene.models import Graphene


class GrapheneForm(forms.ModelForm):
    class Meta:
        model = Graphene
        fields = ['size_x', 'size_y', 'pbc', 'cnt', 'hole', 'radius']
        labels = {
            'size_x': 'x size',
            'size_y': 'y size',
            'pbc': 'pbc',
            'cnt': 'cnt',
            'hole': 'hole',
            'radius': 'radius',
        }
        widgets = {
            'pbc': forms.RadioSelect(choices=())
        }
        error_messages = {
            'radius': {
                'blank': 'The hole must be smaller than x and y size'
            }
        }

    def clean_radius(self):
        r = self.cleaned_data.get('radius')
        x = self.cleaned_data.get('size_x')
        y = self.cleaned_data.get('size_y')

        if float(r) >= float(x)/2.0 or float(r) >= float(y)/2.0:
            r = ''
        return r
