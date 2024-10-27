from django.db import models
from django import forms

# Create your models here.
class NP_C5(models.Model):
    name = models.CharField(max_length=50, default="My NP")
    NPtype = models.CharField(max_length=20,
                              choices=[('Deca','Decahedron'),
                                       ('Ino','Ino-Decahedron'),
                                       ('Marks','Marks-decahedron')]
                              )
    # cilinder specification
    Composition = models.CharField(max_length=50, default="Xx")
    InterAtomDistance = models.DecimalField(max_digits=2, decimal_places=2, default=2.58)
    NumberRings = models.PositiveSmallIntegerField(default=4)

    def __str__(self):
        return self.name

class NP_deca(models.Model):
    name = models.CharField(max_length=50, default="My NP")
    Composition = models.CharField(max_length=50, default="Yy")
    InterAtomDistance = models.DecimalField(max_digits=2, decimal_places=2, default=2.38)
    NumberRings = models.PositiveSmallIntegerField(default=4)

class NP_ino(models.Model):
    name = models.CharField(max_length=50, default="My NP")
    Composition = models.CharField(max_length=50, default="Yy")
    InterAtomDistance = models.DecimalField(max_digits=2, decimal_places=2, default=2.38)
    NumberRings = models.PositiveSmallIntegerField(default=4)
    CilinderExtend = models.PositiveSmallIntegerField(default=4)

class NP_marks(models.Model):
    name = models.CharField(max_length=50, default="My NP")
    Composition = models.CharField(max_length=50, default="Yy")
    InterAtomDistance = models.DecimalField(max_digits=2, decimal_places=2, default=2.38)
    NumberRings = models.PositiveSmallIntegerField(default=4)
    Throughs = models.PositiveSmallIntegerField(default=1)
    CilinderExtend = models.PositiveSmallIntegerField(default=4)




#### ---- Surface forms
class Surf_ManualCell(models.Model):
    name = models.CharField(max_length=50, default="My Surf")
    # Basis definition
    BaseVectors = models.CharField(max_length=50, default="2.4620 2.4620 3.4818")
    BaseAngles = models.CharField(max_length=50, default="80.0 80.0 90.0")
    # Atomic coordinate box
    AtomCoords = models.TextField(default="Ni 0.00 0.00 0.00 \nCo 0.50 0.50 0.50 ")
    #### ---- Construction of surfaces
    # Surface definition
    NumberLayers = models.PositiveSmallIntegerField(default=4)
    # exploration ranges
    ExploreUnits = models.CharField(max_length=50, default="4 10")
    ExploreAngles = models.CharField(max_length=50, default="20.0 160.0")
    # explore Isotropy
    ExploreIsotMax = models.DecimalField(default=3.5, decimal_places=1, max_digits=2)
    # Resulting surface
    UnitsInSurface = models.IntegerField(default=4)
    AddedVacuum = models.DecimalField(default=10.0, decimal_places=2, max_digits=4)
    # Generate Surface
    GenSurfA = models.CharField(max_length=50, default="0 2")
    GenSurfB = models.CharField(max_length=50, default="2 2")
    # FixZ
    FixZbool = models.BooleanField(default=True)



