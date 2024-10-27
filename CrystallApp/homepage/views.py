from django.shortcuts import render
from datetime  import datetime
import numpy as np
from . import CellFunctions as cfcn
from. import NPsC5Functions as c5fcn
from django.http import HttpResponse
from . import OutputGeom as ogeom

xx =  datetime.now()
CurrentAtomCoords = 'Default text'

from . import forms

# Create your views here.

def home(request):
    return render(request, 'homepage/home.html', {'date':xx})

def NPsHome(request):
    NPform = forms.NP_C5form()
    return render(request, 'homepage/NPsHome.html', {'date':xx, 'form':NPform})

def NPsDeca(request):
    # form information
    NPdecaform = forms.NP_deca(initial={"name":'My_NP'})
    context = {'date':xx, 'form':NPdecaform, 'ShowResult':False}
    # post information
    if request.method == 'POST':
        # update default values
        NPdecaform = forms.NP_deca(initial={"name": request.POST["name"],
                                            "Composition": request.POST["Composition"],
                                            'InterAtomDistance': request.POST["InterAtomDistance"],
                                            "NumberRings": request.POST["NumberRings"]
                                            })
        context = {'date': xx, 'form': NPdecaform, 'ShowResult': False}

        print("Submited request to NPs:deca")
        # Activate showing results
        context['ShowResult'] = request.POST["name"]
        # Pass on data
        context['name'] = request.POST["name"]
        context['Composition'] = request.POST["Composition"]
        context['InterAtomDistance'] = request.POST["InterAtomDistance"]
        context["NumberRings"] = request.POST["NumberRings"]

        # Create the NP
        NPDeca = c5fcn.CreateNP_C5(dRef=float(request.POST["InterAtomDistance"]),
                                   rings=int(request.POST["NumberRings"]),
                                   heigh=0,
                                   marks=0.,
                                   element=request.POST["Composition"],
                                   facet_extension=1.,
                                   column_extension=1.)
        CurrentAtomCoords = NPDeca
        NPContainer = c5fcn.NPcontainer(NPDeca)
        # print("NPDeca:")
        # Deca NP special parameters
        context["NAtoms"] = len(NPDeca)
        context["boxX"] = "{:.2f}".format(NPContainer[0])
        context["boxY"] = "{:.2f}".format(NPContainer[1])
        context["boxZ"] = "{:.2f}".format(NPContainer[2])
        context["XYZfile"] = ogeom.write_XYZ(NPDeca,
                                             context['name'],
                                             [context['Composition'] for i in NPDeca])
        context["POSCARfile"] = ogeom.write_POSCAR([[i+j/2+2 for i,j in zip(iAt, [float(context["boxX"]),
                                                                              float(context["boxY"]),
                                                                              0]
                                                                        )
                                                     ]
                                                    for iAt in NPDeca],
                                                   [[float(context["boxX"])+4, 0., 0.],
                                                    [0., float(context["boxY"])+4,0.],
                                                    [0., 0., float(context["boxZ"])+4]],
                                             [context['Composition'] for i in NPDeca])
        if 'submit' in request.POST:
            # render computed page
            return render(request, 'homepage/NPsDeca.html', context)

        elif "POSCARfile" in request.POST:
            # Download XYZ file response
            [print("**"*50) for i in range(2)]
            print("POSCAR has been requested")
            # Generate and download file
            filename = "_".join([i for i in context['name'].split()])+"_POSCAR"
            content = context["POSCARfile"]
            response = HttpResponse(content, content_type='text/plain')
            response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
            return response

        elif "XYZfile" in request.POST:
            # Download XYZ file response
            [print("**"*50) for i in range(2)]
            print("XYZ has been requested")
            # Generate and download file
            filename = context['name']+".xyz"
            content = context["XYZfile"]
            response = HttpResponse(content, content_type='text/plain')
            response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
            return response
    else:
        print("Request method is not POST")

    # this is the default, clean page response
    return render(request, 'homepage/NPsDeca.html', context)



def NPsIno(request):
    # form information
    NPinoform = forms.NP_ino(initial={"name":'My_NP'})
    context = {'date':xx, 'form':NPinoform, 'ShowResult':False}
    # post information
    if request.method == 'POST':
        # update default values
        NPinoform = forms.NP_ino(initial={"name": request.POST["name"],
                                          "Composition": request.POST["Composition"],
                                          'InterAtomDistance': request.POST["InterAtomDistance"],
                                          "NumberRings": request.POST["NumberRings"],
                                          "CilinderExtend": request.POST["CilinderExtend"]
                                          })
        context = {'date': xx, 'form': NPinoform, 'ShowResult': False}

        print("Submited request to NPs:ino")
        # Activate showing results
        context['ShowResult'] = request.POST["name"]
        # Pass on data
        context['name'] = request.POST["name"]
        context['Composition'] = request.POST["Composition"]
        context['InterAtomDistance'] = request.POST["InterAtomDistance"]
        context["NumberRings"] = request.POST["NumberRings"]
        context["CilinderExtend"] = request.POST["CilinderExtend"]

        # Create the NP
        NPIno = c5fcn.CreateNP_C5(dRef=float(request.POST["InterAtomDistance"]),
                                   rings=int(request.POST["NumberRings"]),
                                   heigh=int(request.POST["CilinderExtend"]),
                                   marks=0.,
                                   element=request.POST["Composition"],
                                   facet_extension=1.,
                                   column_extension=1.)
        CurrentAtomCoords = NPIno
        NPContainer = c5fcn.NPcontainer(NPIno)
        # print("NPDeca:")
        # Deca NP special parameters
        context["NAtoms"] = len(NPIno)
        context["boxX"] = "{:.2f}".format(NPContainer[0])
        context["boxY"] = "{:.2f}".format(NPContainer[1])
        context["boxZ"] = "{:.2f}".format(NPContainer[2])
        context["XYZfile"] = ogeom.write_XYZ(NPIno,
                                             context['name'],
                                             [context['Composition'] for i in NPIno])
        context["POSCARfile"] = ogeom.write_POSCAR([[i+j/2+2 for i,j in zip(iAt, [float(context["boxX"]),
                                                                              float(context["boxY"]),
                                                                              0]
                                                                        )
                                                     ]
                                                    for iAt in NPIno],
                                                   [[float(context["boxX"])+4, 0., 0.],
                                                    [0., float(context["boxY"])+4,0.],
                                                    [0., 0., float(context["boxZ"])+4]],
                                             [context['Composition'] for i in NPIno])
        if 'submit' in request.POST:
            # render computed page
            return render(request, 'homepage/NPsIno.html', context)

        elif "POSCARfile" in request.POST:
            # Download XYZ file response
            [print("**"*50) for i in range(2)]
            print("POSCAR has been requested")
            # Generate and download file
            filename = "_".join([i for i in context['name'].split()])+"_POSCAR"
            content = context["POSCARfile"]
            response = HttpResponse(content, content_type='text/plain')
            response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
            return response

        elif "XYZfile" in request.POST:
            # Download XYZ file response
            [print("**"*50) for i in range(2)]
            print("XYZ has been requested")
            # Generate and download file
            filename = context['name']+".xyz"
            content = context["XYZfile"]
            response = HttpResponse(content, content_type='text/plain')
            response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
            return response
    else:
        print("Request method is not POST")

    # this is the default, clean page response
    return render(request, 'homepage/NPsIno.html', context)

def NPsMarks(request):
    # form information
    NPmarksform = forms.NP_marks(initial={"name":'My_NP'})
    context = {'date':xx, 'form':NPmarksform, 'ShowResult':False}
    # post information
    if request.method == 'POST':
        # update default values
        NPmarksform = forms.NP_marks(initial={"name": request.POST["name"],
                                          "Composition": request.POST["Composition"],
                                          'InterAtomDistance': request.POST["InterAtomDistance"],
                                          "NumberRings": request.POST["NumberRings"],
                                          "CilinderExtend": request.POST["CilinderExtend"],
                                          "Throughs": request.POST["Throughs"]
                                          })
        context = {'date': xx, 'form': NPmarksform, 'ShowResult': False}

        print("Submited request to NPs:ino")
        # Activate showing results
        context['ShowResult'] = request.POST["name"]
        # Pass on data
        context['name'] = request.POST["name"]
        context['Composition'] = request.POST["Composition"]
        context['InterAtomDistance'] = request.POST["InterAtomDistance"]
        context["NumberRings"] = request.POST["NumberRings"]
        context["CilinderExtend"] = request.POST["CilinderExtend"]
        context["Throughs"] = request.POST["Throughs"]

        # Create the NP
        NPIno = c5fcn.CreateNP_C5(dRef=float(request.POST["InterAtomDistance"]),
                                   rings=int(request.POST["NumberRings"]),
                                   heigh=int(request.POST["CilinderExtend"]),
                                   marks=int(request.POST["Throughs"]),
                                   element=request.POST["Composition"],
                                   facet_extension=1.,
                                   column_extension=1.)
        CurrentAtomCoords = NPIno
        NPContainer = c5fcn.NPcontainer(NPIno)
        # print("NPDeca:")
        # Deca NP special parameters
        context["NAtoms"] = len(NPIno)
        context["boxX"] = "{:.2f}".format(NPContainer[0])
        context["boxY"] = "{:.2f}".format(NPContainer[1])
        context["boxZ"] = "{:.2f}".format(NPContainer[2])
        context["XYZfile"] = ogeom.write_XYZ(NPIno,
                                             context['name'],
                                             [context['Composition'] for i in NPIno])
        context["POSCARfile"] = ogeom.write_POSCAR([[i+j/2+2 for i,j in zip(iAt, [float(context["boxX"]),
                                                                              float(context["boxY"]),
                                                                              0]
                                                                        )
                                                     ]
                                                    for iAt in NPIno],
                                                   [[float(context["boxX"])+4, 0., 0.],
                                                    [0., float(context["boxY"])+4,0.],
                                                    [0., 0., float(context["boxZ"])+4]],
                                             [context['Composition'] for i in NPIno])
        if 'submit' in request.POST:
            # render computed page
            return render(request, 'homepage/NPsMarks.html', context)

        elif "POSCARfile" in request.POST:
            # Download XYZ file response
            [print("**"*50) for i in range(2)]
            print("POSCAR has been requested")
            # Generate and download file
            filename = "_".join([i for i in context['name'].split()])+"_POSCAR"
            content = context["POSCARfile"]
            response = HttpResponse(content, content_type='text/plain')
            response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
            return response

        elif "XYZfile" in request.POST:
            # Download XYZ file response
            [print("**"*50) for i in range(2)]
            print("XYZ has been requested")
            # Generate and download file
            filename = context['name']+".xyz"
            content = context["XYZfile"]
            response = HttpResponse(content, content_type='text/plain')
            response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
            return response
    else:
        print("Request method is not POST")

    # this is the default, clean page response
    return render(request, 'homepage/NPsMarks.html', context)



def SurfsHome(request):
    #create object
    response = HttpResponse("str(context['name'])", content_type="text/plain")
    filename = "my-file.txt"
    content = 'any string generated by django'
    response = HttpResponse(content, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)

    return response
    #return render(request, 'homepage/SurfsHome.html', {'date':xx})




def SurfsManualCell(request):
    # form information
    print("*"*90)
    print("*"*5 + " Requested render of SurfsManualCell page")

    # post information
    if request.method == 'POST':
        print("request method is POST")
        # update default values
        SurfManualform = forms.Surf_ManualCell(
            initial={"name": request.POST["name"],
                     "BaseVectors": request.POST["BaseVectors"],
                     "BaseAngles": request.POST["BaseAngles"],
                     "AtomCoords": request.POST["AtomCoords"],
                     "NumberLayers": request.POST["NumberLayers"],
                     "ExploreUnits": request.POST["ExploreUnits"],
                     "ExploreAngles": request.POST["ExploreAngles"],
                     "ExploreIsotMax": request.POST["ExploreIsotMax"],
                     "UnitsInSurface": request.POST["UnitsInSurface"],
                     "FixZbool": bool(request.POST.get("FixZbool", False)),
                     "AddedVacuum": request.POST["AddedVacuum"],
                     "GenSurfA": request.POST["GenSurfA"],
                     "GenSurfB": request.POST["GenSurfB"]
                     }
        )
        print("Surfs Manual form is created")
        context = {'date': xx, 'form': SurfManualform,
                   'ShowResultSingle': False, 'ShowResultList': False , 'ShowResultVectors':False}

        print("Submited request to Surface creation")
        # Pass on data
        context['name'] = request.POST["name"]
        context['BaseVectors'] = request.POST["BaseVectors"]
        context['BaseAngles'] = request.POST["BaseAngles"]
        context['AtomCoords'] = request.POST["AtomCoords"]
        context['NumberLayers'] = request.POST["NumberLayers"]
        context['ExploreUnits'] = request.POST["ExploreUnits"]
        context['ExploreAngles'] = request.POST["ExploreAngles"]
        context['ExploreIsotMax'] = request.POST["ExploreIsotMax"]
        context['UnitsInSurface'] = request.POST["UnitsInSurface"]
        context['FixZbool'] = bool(request.POST.get("FixZbool", False))
        context['AddedVacuum'] = request.POST['AddedVacuum']
        context['GenSurfA'] = request.POST['GenSurfA']
        context['GenSurfB'] = request.POST['GenSurfB']

        #### Work on cell
        # treat data base cell
        #latticebase = cfcn.abcCell2vectCell([float(i) for i in context['BaseVectors'].split()],
        #                                    [float(i) for i in context['BaseAngles'].split()],
        #                                   angleformat='deg')

        try:

            if 'SuperCellList' in request.POST:
                print(" Requested to create cell list")
                # Activate showing results
                context['ShowResultList'] = request.POST["name"]
                context['ShowResultSingle'] = False
                context['ShowResultVectors'] = False

                # compute cell options
                SupercellsListOut = cfcn.Supercells([float(i) for i in context['BaseVectors'].split()],
                                                   float(context['BaseAngles'].split()[2]) * np.pi / 180,
                                                   A_min=int(context['ExploreUnits'].split()[0]),
                                                   A_max=int(context['ExploreUnits'].split()[1]),
                                                   theta_min=float(context['ExploreAngles'].split()[0]),
                                                   theta_max=float(context['ExploreAngles'].split()[1]),
                                                   AnsyMax=float(context['ExploreIsotMax']))

                print("*" * 90)
                print(SupercellsListOut[0])
                context['ListElementsNumber'] = len(SupercellsListOut[0])


            elif 'SuperCellListDownload' in request.POST: # Generate and download file

                # compute cell options
                SupercellsListOut = cfcn.Supercells([float(i) for i in context['BaseVectors'].split()],
                                                   float(context['BaseAngles'].split()[2]) * np.pi / 180,
                                                   A_min=int(context['ExploreUnits'].split()[0]),
                                                   A_max=int(context['ExploreUnits'].split()[1]),
                                                   theta_min=float(context['ExploreAngles'].split()[0]),
                                                   theta_max=float(context['ExploreAngles'].split()[1]),
                                                   AnsyMax=float(context['ExploreIsotMax']))

                print("*" * 90)
                print(SupercellsListOut[0])
                print('-'*80)
                print(SupercellsListOut[1])

                # Parse list for output
                OUTtitles = ['Base Units',
                             'Vector A', 'Vector B', 'angle (grad)',
                             'Area (A2)', 'Isotropic']
                OUT = '#'*4 + ' Recomended supercells\n'
                OUT += OUTtitles[0].rjust(11).ljust(12) # units
                OUT += OUTtitles[1].rjust(11).ljust(12)
                OUT += OUTtitles[2].rjust(11).ljust(12)
                OUT += OUTtitles[3].rjust(13).ljust(14) # angle
                OUT += OUTtitles[4].rjust(10).ljust(11) # area
                OUT += OUTtitles[5].rjust(10).ljust(11) # Isotropy
                OUT += '\n'
                OUT += '-'*90+'\n'
                for iLine in SupercellsListOut[0]:
                    OUT += str(int(iLine[0])).rjust(11).ljust(12)
                    OUT += str(iLine[3][0]).rjust(11).ljust(12)
                    OUT += str(iLine[3][1]).rjust(11).ljust(12)
                    OUT += str('{:.2f}'.format(iLine[2])).rjust(13).ljust(14)
                    OUT += str('{:.4f}'.format(iLine[5])).rjust(10).ljust(11)
                    OUT += str(bool(iLine[4])).rjust(10).ljust(11)
                    #OUT += str(iLine)
                    OUT += '\n'

                OUT += '\n'+'='*90+'\n'

                OUT += '#'*4 + ' All supercells combinations\n'
                OUT += OUTtitles[0].rjust(11).ljust(12) # units
                OUT += OUTtitles[1].rjust(11).ljust(12)
                OUT += OUTtitles[2].rjust(11).ljust(12)
                OUT += OUTtitles[3].rjust(13).ljust(14) # angle
                OUT += OUTtitles[4].rjust(10).ljust(11) # area
                OUT += OUTtitles[5].rjust(10).ljust(11) # Isotropy
                OUT += '\n'
                for iLine in SupercellsListOut[1].keys():
                    print('key:')
                    print(iLine)
                    print(SupercellsListOut[1][iLine])
                    OUT += '-' * 90 + '\n'
                    for j in SupercellsListOut[1][iLine]:
                        # OUT += str(j)
                        OUT += str(int(j['N'])).rjust(11).ljust(12)
                        OUT += str(j['Vects'][0]).rjust(11).ljust(12)
                        OUT += str(j['Vects'][1]).rjust(11).ljust(12)
                        OUT += str('{:.2f}'.format(j['angle'])).rjust(13).ljust(14)
                        OUT += str('{:.4f}'.format(j['area'])).rjust(10).ljust(11)
                        OUT += str(bool(j['Iso'])).rjust(10).ljust(11)
                        OUT += '\n'
                    OUT += '-'*30



                filename = "_".join([i for i in context['name'].split()])+".txt"
                content = OUT
                response = HttpResponse(content, content_type='text/plain')
                response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
                return response

            elif 'SingleSuperCell' in request.POST:
                print(" Requested single cell creation")
                # Activate showing results
                context['ShowResultSingle'] = request.POST["name"]
                context['ShowResultList'] = False
                context['ShowResultVectors'] = False

                # compute cell options
                SupercellsListOut = cfcn.Supercells([float(i) for i in context['BaseVectors'].split()],
                                                   float(context['BaseAngles'].split()[2]) * np.pi / 180,
                                                   A_min=int(context['UnitsInSurface']),
                                                   A_max=int(context['UnitsInSurface']),
                                                   theta_min=10,
                                                   theta_max=190,
                                                   AnsyMax=float(context['ExploreIsotMax']))
                context['SingleCell'] = 0.

            elif 'POSCARfileSingle' in request.POST:
                print(" Requested single cell creation")
                # Activate showing results
                context['ShowResultSingle'] = request.POST["name"]
                context['ShowResultList'] = False
                context['ShowResultVectors'] = False

                # compute cell options
                SupercellsListOut = cfcn.Supercells([float(i) for i in context['BaseVectors'].split()],
                                                   float(context['BaseAngles'].split()[2]) * np.pi / 180,
                                                   A_min=int(context['UnitsInSurface']),
                                                   A_max=int(context['UnitsInSurface']),
                                                   theta_min=10,
                                                   theta_max=190,
                                                   AnsyMax=float(context['ExploreIsotMax']))

                Vect1 = str(SupercellsListOut[0][0][3][0])[1:-1]
                Vect2 = str(SupercellsListOut[0][0][3][1])[1:-1]

                # Define lattice vectors
                latticebase = cfcn.abcCell2vectCell([float(i) for i in context['BaseVectors'].split()],
                                                    [float(i) for i in context['BaseAngles'].split()],
                                                    angleformat='deg')
                # Parse atons in base
                print("Parse lattice base")
                AtomSitesBase = [[i.split()[0] for i in context['AtomCoords'].split('\n')],
                                 [[float(j) for j in i.split()[1:]] for i in context['AtomCoords'].split('\n')]
                                 ]
                print(f'results: {AtomSitesBase}\n\n')

                # Generate Supercell
                ThisSuperCell =cfcn.GenerateSupercell(
                    [int(i) for i in Vect1.split()],
                    [int(i) for i in Vect2.split()],
                    AtomsSites=AtomSitesBase,
                    CellParams= latticebase,
                    nLayers=int(request.POST['NumberLayers']),
                    vacuum=float(request.POST['AddedVacuum']),
                    MakeZortho=context['FixZbool']
                )
                # Generate POSCAR
                filename = "_".join([i for i in context['name'].split()])+"_POSCAR"
                content = cfcn.WritePoscarMini(ThisSuperCell['lattice'],
                                               ThisSuperCell['Atoms'],
                                               ThisSuperCell['Names'],
                                               Coordinates='Direct')
                response = HttpResponse(content, content_type='text/plain')
                response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
                return response

            # Genearte supercell from base vector definition
            elif 'GenerateSuperCellVectors' in request.POST:
                print(" Requested single cell creation")
                # Activate showing results
                context['ShowResultVectors'] = request.POST["name"]
                context['ShowResultList'] = False
                context['ShowResultSingle'] = False

                # Define lattice vectors
                latticebase = cfcn.abcCell2vectCell([float(i) for i in context['BaseVectors'].split()],
                                                   [float(i) for i in context['BaseAngles'].split()],
                                                   angleformat='deg')
                # Parse atons in base
                print("Parse lattice base")
                AtomSitesBase = [[i.split()[0] for i in context['AtomCoords'].split('\n')],
                                 [[float(j) for j in i.split()[1:]] for i in context['AtomCoords'].split('\n')]
                                 ]
                print(f'results: {AtomSitesBase}\n\n')

                # Generate Supercell
                ThisSuperCell =cfcn.GenerateSupercell(
                    [int(i) for i in request.POST['GenSurfA'].split()],
                    [int(i) for i in request.POST['GenSurfB'].split()],
                    AtomsSites=AtomSitesBase,
                    CellParams= latticebase,
                    nLayers=int(request.POST['NumberLayers']),
                    vacuum=float(request.POST['AddedVacuum']),
                    MakeZortho=context['FixZbool']
                )
                context['GenSurfaceNumAtoms'] = len(ThisSuperCell['Names'])
                context['GenSurfaceArea'] = '{:.2f}'.format(cfcn.LatticeABarea(ThisSuperCell['lattice']))

            elif 'POSCARfileVectors' in request.POST:
                print(" Requested single cell creation to POSCAR download")
                # Activate showing results
                context['ShowResultVectors'] = request.POST["name"]
                context['ShowResultList'] = False
                context['ShowResultSingle'] = False

                # Define lattice vectors
                latticebase = cfcn.abcCell2vectCell([float(i) for i in context['BaseVectors'].split()],
                                                   [float(i) for i in context['BaseAngles'].split()],
                                                   angleformat='deg')
                # Parse atons in base
                print("Parse lattice base")
                AtomSitesBase = [[i.split()[0] for i in context['AtomCoords'].split('\n')],
                                 [[float(j) for j in i.split()[1:]] for i in context['AtomCoords'].split('\n')]
                                 ]
                print(f'results: {AtomSitesBase}\n\n')

                # Generate Supercell
                ThisSuperCell =cfcn.GenerateSupercell(
                    [int(i) for i in request.POST['GenSurfA'].split()],
                    [int(i) for i in request.POST['GenSurfB'].split()],
                    AtomsSites=AtomSitesBase,
                    CellParams= latticebase,
                    nLayers=int(request.POST['NumberLayers']),
                    vacuum=float(request.POST['AddedVacuum']),
                    MakeZortho=context['FixZbool']
                )
                # Generate POSCAR
                filename = "_".join([i for i in context['name'].split()])+"_POSCAR"
                content = cfcn.WritePoscarMini(ThisSuperCell['lattice'],
                                               ThisSuperCell['Atoms'],
                                               ThisSuperCell['Names'],
                                               Coordinates='Direct')
                response = HttpResponse(content, content_type='text/plain')
                response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
                return response

            elif 'XYZfileVectors' in request.POST:
                print(" Requested single cell creation to POSCAR download")
                # Activate showing results
                context['ShowResultVectors'] = request.POST["name"]
                context['ShowResultList'] = False
                context['ShowResultSingle'] = False

                # Define lattice vectors
                latticebase = cfcn.abcCell2vectCell([float(i) for i in context['BaseVectors'].split()],
                                                   [float(i) for i in context['BaseAngles'].split()],
                                                   angleformat='deg')
                # Parse atons in base
                print("Parse lattice base")
                AtomSitesBase = [[i.split()[0] for i in context['AtomCoords'].split('\n')],
                                 [[float(j) for j in i.split()[1:]] for i in context['AtomCoords'].split('\n')]
                                 ]
                print(f'results: {AtomSitesBase}\n\n')

                # Generate Supercell
                ThisSuperCell =cfcn.GenerateSupercell(
                    [int(i) for i in request.POST['GenSurfA'].split()],
                    [int(i) for i in request.POST['GenSurfB'].split()],
                    AtomsSites=AtomSitesBase,
                    CellParams= latticebase,
                    nLayers=int(request.POST['NumberLayers']),
                    vacuum=float(request.POST['AddedVacuum']),
                    MakeZortho=context['FixZbool']
                )
                # Generate POSCAR
                filename = "_".join([i for i in context['name'].split()])+".xyz"
                content = cfcn.WriteXYZMini(ThisSuperCell['lattice'],
                                               ThisSuperCell['Atoms'],
                                               ThisSuperCell['Names'])
                response = HttpResponse(content, content_type='text/plain')
                response['Content-Disposition'] = 'attachment; filename={0}'.format(filename)
                return response

        except:
            return render(request, 'homepage/SomethingWrong.html', context)




        print("proceding to render with context")
        #print(context)
        return render(request, 'homepage/SurfsManualCell.html', context)


    else:
        print("Request method is not POST")
        print("Returning default response with clear page")

    # this is the default, clean page response
    SurfManualform = forms.Surf_ManualCell(initial={"name":'My_surf'})
    context = {'date':xx, 'form':SurfManualform,
               'ShowResultSingle':False, 'ShowResultList':False, 'ShowResultVectors':False}
    return render(request, 'homepage/SurfsManualCell.html', context)

