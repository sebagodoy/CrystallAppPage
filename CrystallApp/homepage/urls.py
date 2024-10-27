from django.urls import path
from . import views


urlpatterns = [
    path('', views.home),
    path('NPsHome', views.NPsHome),
    path('NPsDeca', views.NPsDeca),
    path('NPsIno', views.NPsIno),
    path('NPsMarks', views.NPsMarks),
    path('SurfsManualCell', views.SurfsManualCell)
]