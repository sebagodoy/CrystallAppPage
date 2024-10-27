from django import forms
from . import models

class NP_C5form(forms.ModelForm):
    class Meta:
        model = models.NP_C5
        fields = ['name','NPtype','Composition','InterAtomDistance','NumberRings']

class NP_C5typeform(forms.ModelForm):
    class Meta:
        mode = models.NP_C5
        fields = ['NPtype']

class NP_deca(forms.ModelForm):
    class Meta:
        model = models.NP_deca
        fields = ['name','Composition','InterAtomDistance','NumberRings']

class NP_ino(forms.ModelForm):
    class Meta:
        model = models.NP_ino
        fields = ['name','Composition','InterAtomDistance','NumberRings','CilinderExtend']

class NP_marks(forms.ModelForm):
    class Meta:
        model = models.NP_marks
        fields = ['name','Composition','InterAtomDistance','NumberRings','CilinderExtend', 'Throughs']



#### ---- Surface models
class Surf_ManualCell(forms.ModelForm):
    #AtomCoords = forms.CharField(widget=forms.Textarea)
    class Meta:
        model = models.Surf_ManualCell
        fields = ['name',
                  'BaseVectors', 'BaseAngles', 'AtomCoords',
                  'NumberLayers',
                  'ExploreUnits','ExploreAngles',
                  'ExploreIsotMax',
                  'UnitsInSurface', 'AddedVacuum', 'FixZbool', 'GenSurfA', 'GenSurfB']
        #widgets = {
        #    'BaseVectors': forms.CharField(attrs={'style': 'width:500px'})
        #}
        #def __init__(self, *args, **kwargs):
        #    super(Surf_ManualCell, self).__init__(*args, **kwargs)  # Call to ModelForm constructor
        #    self.fields['BaseVectors'].widget.attrs['cols'] = 5
        #    self.fields['BaseVectors'].widget.attrs['rows'] = 5


