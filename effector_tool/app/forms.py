from django import forms

class ContactCourse(forms.Form):
    sequencia = forms.CharField(label='Sequência',
                                widget=forms.Textarea
                                )