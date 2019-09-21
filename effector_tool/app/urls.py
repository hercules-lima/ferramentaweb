from django.contrib import admin
from django.urls import path
from .views import *

urlpatterns = [
    path('classificar_proteina/', classificar_proteina, name='classificar_proteina' ),
]
