import numpy as np


debugbol = True


def debug(intext, init=False):
    if init:
        print("=" * 90 + "\n")
    if debugbol:
        print(intext)
    return None


### --- parse hkl
def parsehkl(in_hkl):
    '''Check hkl text'''
    try:
        hkl = [int(i) for i in in_hkl.split()]
        return True if len(hkl) == 2 else False
    except:
        return False




### --- cell format transformations
def abcCell2vectCell(abc, AlphaBetaGamma, angleformat="rad"):
    '''Helper function: takes (a, b, c) , (alpha, beta, gamma in rad) defined cell
    and returns a [vector1, vector2, vector3] defined cell.'''
    debug(f"Called abcCell2vectCell with abc:{abc} and angles {AlphaBetaGamma} in format {angleformat}", init=True)
    if angleformat == "deg":
        AlphaBetaGamma = [i * np.pi / 180 for i in AlphaBetaGamma]
    avect = [abc[0], 0., 0.]
    bvect = [abc[1] * np.cos(AlphaBetaGamma[2]),
             abc[1] * np.sin(AlphaBetaGamma[2]),
             0.]
    _1 = abc[2] * np.cos(AlphaBetaGamma[0])
    _2 = abc[2] * np.cos(AlphaBetaGamma[1])
    _3base = np.sqrt(_1 ** 2 + _2 ** 2)
    cvect = [_1, _2, np.sqrt(abc[2] ** 2 - _3base ** 2)]
    print(f"returns: {[avect, bvect, cvect]}\n\n")
    return np.array([avect, bvect, cvect])


#### ---- angle
def MyAngle(_u, _v, format='deg'):
    '''Helper function: input vectors (2D or 3D) and returns
    angle (format='deg' or 'rad').'''
    debug(f"Called angle between v1:{_u}, v2:{_v}, out format:{format}")
    _ang = np.arccos(np.dot(_u, _v) / (np.linalg.norm(_u) * np.linalg.norm(_v)))
    if format == 'deg':
        return _ang * (180 / np.pi)
    if format == 'rad':
        return _ang
    else:
        raise TypeError("Requested weird angle format")


def abcGammaArea(abc, gamma):
    '''Helper function: computes area from abc iterable vector and gamma (rad) angle'''
    debug(f"Computed area for abc:{abc}, angle:{gamma} (rad)")
    return abc[0] * abc[1] * np.sin(gamma)

def LatticeABarea(lattice):
    '''Helper function: computes a,b area from lattice vectors'''
    return np.linalg.norm(np.cross(lattice[0], lattice[1]))

########################################################################################################################

#### ---- Find supercells
def Supercells(ab,  # ab, tuple of baseCell vector lengths
               angle,  # angle(a,b), base vectors, float radians
               n=10,  # grid search size
               f=1.5,
               A_min=.5, A_max=5,  # Range of acceptable units
               theta_min=59, theta_max=121,  # Range of aceptable angles
               AnsyMax=3.
               ):
    '''Comp. Function: searches possible supercells, returns list os possibilities'''
    # La fonction supercells prend comme argument les vecteurs u et v qui définissent le plan atomique
    # sur lequel on cherche des supercellules les plus isotropes possibles. n définit le domaine du plan
    # avec y>0 à explorer, f correspond à l'anisotropie maximum tolérée, A_max le nombre maximum d'atomes
    # souhaités par supercellule, theta_min et theta_max correspondent à l'intervalle d'angle toléré pour
    # la supercellule.
    debug(f"Called Supercells using ab {ab} : angle {angle}", init=True)

    # Base vectors
    slab_a = ab[0]
    slab_b = ab[1]
    slab_parameter_ratio = slab_b / slab_a
    # Cell basis base
    u = [1., 0.]
    v = [np.cos(angle), np.sin(angle)]
    # Orthonormal base
    u_p = np.array([1, 0])
    v_p = np.array([0, 1])

    # 1) Définition d'un maillage sur la partie y>0 avec des vecteurs w non nuls exprimés dans la base M
    mesh = []
    for i in range(-n, n + 1, 1):
        for j in range(0, n + 1, 1):
            w = [i * u_p[0] + j * v_p[0], i * u_p[1] + j * v_p[1]]
            if w[0] != 0 or w[1] != 0:
                mesh.append(np.array(w))
    # mesh: maps pairs x(-n -> n) / y(0 -> n), orthonormal space to work in
    #      [ array([-n, 0]) , array([-n,1]) , ... , array([n, n])]

    # Définition de la matrice de changement de base entre le Maillage M (coordonnées entières) et
    # le plan orthonormé P
    P = np.array([u, v]).T

    # 2) Création de la liste de supercellules qui respectent les contraintes de dimension posées
    results = []  # container
    AllSurfaceOptions = []  # Container
    isotropic = False
    for a in mesh:
        for b in mesh:
            # a et b définissent une supercellule, on les convertit dans la base orthonormée
            # et on calcule le déterminant dans la base orthonormée (= aire normalisée)
            a_p = P.dot(a).T
            b_p = P.dot(b).T
            A_p = abs(np.linalg.det(np.array([a_p, b_p])))

            # calcul du déterminant dans la base M (= nb d'atomes compris dans la cellule)
            A = abs(np.linalg.det(np.array([a, b]).T))
            bA = True if A >= A_min and A <= A_max else False

            # Calcul de l'angle entre les deux vecteurs
            #print(f"a: {a_p} , b: {b_p}, a.b: {np.dot(a_p, b_p)}")
            theta = abs(np.arccos(np.dot(a_p, b_p) / (np.linalg.norm(a_p) * np.linalg.norm(b_p))) * 180 / np.pi)
            bTheta = True if theta > theta_min and theta < theta_max else False

            # TODO: Check this definition
            # Is anisotropy acceptable?
            AnIsotropy = max(np.linalg.norm(a_p), np.linalg.norm(b_p)) / min(np.linalg.norm(a_p), np.linalg.norm(b_p))
            # bAnsy = True if round(AnIsotropy, 5) <=round(f*slab_parameter_ratio,5) else False
            bAnsy = True if round(AnIsotropy, 5) < AnsyMax else False

            # if bA:
            #    print(f"{bAnsy and bA and bTheta}: Option Ansi:{AnIsotropy}({bAnsy}), Angle:{theta}({bTheta}), containing {A} (in range:{bA})")

            # Condition de sélection de la supercellule :
            #       1) nombre d'atomes compris dans la supercellule (bA: A >=1, < A_max),
            #       2) intervalle d'angle autorisé (bTheta: theta_min < theta < theta_max)
            #       3) anisotropie des normes tolérée
            if bA and bTheta and bAnsy:
                # print(f"{bAnsy and bA and bTheta}: Option Ansi:{AnIsotropy}({bAnsy}), Angle:{theta}({bTheta}), containing {A} (in range:{bA})")
                # La cellule trouvée est-elle isotrope ?
                isotropic = True if round(np.linalg.norm(a_p) * slab_parameter_ratio, 5) == round(np.linalg.norm(b_p),
                                                                                                  5) else False

                # On ajoute à la liste préliminaire des résultats les paramètres calculés :
                # calcul de l'aire dans la base orthonormée P en tenant compte des valeurs des paramètres de mailles
                aire = round(100 / (A_p * slab_parameter_ratio * slab_a ** 2 * abs(
                    np.sin(abs(np.arccos(np.dot(u, v)) * 180 / np.pi)))), 4)
                results.append([round(A, 1),  # 0 number of atoms: float?
                                (round(np.linalg.norm(a_p), 5),  # 1
                                 round(np.linalg.norm(b_p), 5)
                                 ),
                                theta,  # 2 angle: float
                                (a, b),  # 3 coordinates: int, int
                                isotropic,  # 4 Is isotropic: bool
                                aire  # 5 1/Area: float
                                ])
                AllSurfaceOptions.append({"N": round(A, 1),
                                          "vLengths": (round(np.linalg.norm(a_p), 5), round(np.linalg.norm(b_p), 5)),
                                          "angle": round(theta),
                                          "area": 1 / aire,
                                          "Vects": (a, b),
                                          "Iso": isotropic})
    print(f" Got how many results {len(results)}")
    if len(results) == 0:
        raise NotImplementedError("Supercells function did not produced any useful results")

    # 3) Tri et sélection des résultats
    results_f = []
    results = sorted(results, key=lambda tup: tup[0])  # Tri selon A
    p = results[0][0]
    choice = results[0]  # Choix du premier candidat
    i = 0
    while i < len(results):  # loop results
        # Pour un même A, on en choisit un unique préféré
        while i < len(results) and results[i][0] == A:
            if results[i][4] == True:
                if choice[4] == False:  # On a un candidat isotrope qui devient meilleur candidat
                    choice = results[i]
                elif abs(choice[2] - 90) > abs(results[i][2] - 90):
                    # On a un candidat isotrope dont l'angle est plus proche de 90° qui devient meilleur candidat.
                    choice = results[i]

            elif choice[4] == False and max(choice[1][0], choice[1][1]) / min(choice[1][0], choice[1][1]) > max(
                    results[i][1][0], results[i][1][1]) / min(results[i][1][0], results[i][1][1]) and abs(
                    choice[2] - 90) > abs(results[i][2] - 90):
                # On a un candidat anisotrope dont l'anisotropie est plus faible et dont l'angle est plus proche de 90° qui devient meilleur candidat.
                choice = results[i]
            i += 1
        if A == choice[0]:  # On ajoute le meilleur candidat à la liste finale de résultats
            results_f.append(choice)
        A += 1
        if i < len(results):
            choice = results[i]  # Choix du premier candidat une fois que A est itéré.

    # 4) group all posible options
    SurfaceOptions = {i[0]: [] for i in results}
    for iOpt in AllSurfaceOptions:
        if theta_min < iOpt["angle"] < theta_max:
            SurfaceOptions[iOpt["N"]].append(iOpt)

    return results_f, SurfaceOptions

def SurfaceListToText(input):
    txt = ""
    for iLine in input:
        txt += str(iLine)+"\n"
    return txt

########################################################################################################################


def Cartesian2Direct(Lattice, Sites):
    '''Transforms sites from cartesian to direct coordinates in a lattice
    Input:  Lattice [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
            Sites [ [Names, Names, ...], [ [at1,xyz], [at2, xyz], ... ](cartesian) ]
    Output: Sites [ [Names, Names, ...], [ [at1,xyz], [at2, xyz], ... ](direct in lattice) ]
    '''
    MC2D = np.linalg.inv(np.array(Lattice).T)
    nSites = [MC2D.dot(iAt) for iAt in Sites[1]]
    return [Sites[0], nSites]


def MakeCellZortho(Lattice, Sites):
    ''' Makes lattice orthogonal in z, so the POSCAR is properly written
    Input:  Lattice1 [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]](non-z-orthogonal)
            Sites [ [Names, Names, ...], [ [at1,xyz], [at2, xyz], ... ](direct-in-Lattice) ]
    Output: Lattice2 [[x1, y1, z1], [x2, y2, z2], [0., 0., z3]](z-orthogonal)
            Sites [ [Names, Names, ...], [ [at1,xyz], [at2, xyz], ... ](direct-in-Lattice2) ]'''
    # Fix lattice
    debug(f"Called makeLatticeOrthoZ with lattice:{Lattice} and sites:{Sites}", init=True)
    debug(f"    Used par of Lattice:{Lattice[:2]}")
    nLattice = np.append(Lattice[:2], [[0., 0., Lattice[2][2]]], axis=0)
    debug(f"    nLattice: {nLattice}")
    # fix sites
    mC2R = np.array(Lattice).T  # matrix transf., direct coord. in original cell to real space
    mR2C2 = np.linalg.inv(nLattice.T)  # matrix transf., real space to new cell space
    RealAtoms = [mC2R.dot(iAt) for iAt in Sites[1]]
    DirectAtoms = [mR2C2.dot(iAt) for iAt in RealAtoms]
    # fix in range [0, 1]
    for i, iAt in enumerate(DirectAtoms):
        for j, jCoord in enumerate(iAt):
            while jCoord >= 1.:
                jCoord -= 1.
            while jCoord < 0.:
                jCoord += 1.
            DirectAtoms[i][j] = jCoord

    debug(f"Return Lattice:{nLattice} sites:{[Sites[0], DirectAtoms]}")
    return nLattice, [Sites[0], DirectAtoms]


def WritePoscarMini(Lattice, Sites, Names, Coordinates="Cartesian"):
    '''Writes POSCAR type file from cell and sites information
    Input:  Lattice=[[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
            Sites= [[x1,y1,z1], [...], ... ]
            Names= [Xx, Y, Z, ... ]
    Output: POSCAR string'''
    # start POSCAR string
    POSCAR = "Supercell \n" + "{:.16f}".format(1.0).rjust(20) + "\n"
    # write cell
    for iCoord in Lattice:
        POSCAR += " ".join(["{:.15f}".format(round(i, 16)).rjust(18) for i in iCoord]) + "\n"
    # Sort by atom type
    print(Names)
    Types = {i: [] for i in Names}
    for idx, iType in enumerate(Names):
        Types[iType].append(Sites[idx])
    # writting atomic content
    for i in list(Types.keys()):
        POSCAR += "   " + i
    POSCAR += "\n"
    for i in list(Types.keys()):
        POSCAR += "   " + str(len(Types[i]))
    POSCAR += "\n"
    POSCAR += "Selective Dynamics\n"
    POSCAR += Coordinates + "\n"
    # Write atoms
    for iType in list(Types.keys()):
        # order atoms in vertical
        iAtomTypeList = Types[iType]
        iAtomTypeList.sort(key=lambda tup: tup[2])
        for iAt in iAtomTypeList:
            POSCAR += " ".join(["{:.16f}".format(k).rjust(20)
                                for k in iAt]) + " F F F \n"
    return POSCAR

def WriteXYZMini(Lattice, Sites, Names):
    '''Writes POSCAR type file from cell and sites information
    Input:  Lattice=[[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
            Sites= [[x1,y1,z1], [...], ... ]
            Names= [Xx, Y, Z, ... ]
    Output: POSCAR string'''
    LatticeArray = np.array(Lattice)
    OUT = str(len(Sites))+'\n'
    OUT += 'Surface\n'
    for iSite, iName in zip(Sites, Names):
        OUT += iName.rjust(3) + ' '
        OUT += ' '.join(['{:.12f}'.format(i).rjust(18) for i in LatticeArray.T.dot(iSite)]) + '\n'

    return OUT

def GenerateSupercell(
        aVect,  # SlabVector-1 in the space of the base cell
        bVect,  # SlabVector-2 in the space of the base cell
        AtomsSites,  # list of atoms
        CellParams,  # cell vectors [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
        nLayers,  # vertical repetitions of the basis cell
        vacuum,  # Vacuum to add, float
        MakeZortho=False
):
    '''Output function: Generates POSCAR file'''

    debug(f"Called GenerateSupercell, aVect:{aVect}, bVect:{bVect}, sites:{AtomsSites}, Lattice:{CellParams}",
          init=True)

    #### ---- Decide order of base vectors, aVect closer to angle zero in plane
    if MyAngle(aVect, [1, 0]) > MyAngle(bVect, [1, 0]):
        aVect, bVect = (bVect, aVect)

    #### ---- Construct minimal mesh to scan atoms
    minXmesh = min(0, aVect[0], bVect[0], aVect[0] + bVect[0]) - 1
    maxXmesh = max(0, aVect[0], bVect[0], aVect[0] + bVect[0]) + 1
    maxYmesh = max(0, aVect[1], bVect[1], aVect[1] + bVect[1]) + 1
    mesh = [np.array([i, j, 0]) for i in range(minXmesh, maxXmesh + 1) for j in range(0, maxYmesh + 1)]
    debug(f"Contructed  mesh 3D X:{minXmesh} to {maxXmesh}, {maxYmesh} from origin, Z=0")

    #### ---- Map onto real space mesh
    aVect3D = np.append(aVect, 0.)
    bVect3D = np.append(bVect, 0.)
    # mC2R maps 3D direct coordinates (cell space) to real space
    mC2R = CellParams[:2, :2].T
    mC2R3D = CellParams.T
    AtomsInMesh = []  # Container, sites - real coordinates
    NamesInMesh = []  # Container, sites - names
    for mi in mesh:
        for iName, iAt in enumerate(AtomsSites[1]):
            # print(f"iName : {iName}")
            # print(f"iAt index: {iAt}")
            AtomsInMesh.append(mC2R3D.dot(mi) +  # cell translation
                               mC2R3D.dot(iAt)  # site within cell
                               )
            NamesInMesh.append(AtomsSites[0][iName])

    # map supercell vectors onto real space
    aVect3Dreal = mC2R3D.dot(aVect3D)
    bVect3Dreal = mC2R3D.dot(bVect3D)
    cVect3Dreal = CellParams[2]

    # rotate c vector from [1, 0] plane to fit with aVectReal
    # cVectAngle = MyAngle(aVectReal, [1., 0.], format='rad')
    # Mrotc = np.array([[np.cos(2 * np.pi - cVectAngle), -np.sin(2 * np.pi - cVectAngle)],
    #                 [np.sin(2 * np.pi - cVectAngle), np.cos(2 * np.pi - cVectAngle)]])
    # cVectReal3D = np.append(Mrotc.dot(CellParams[2][:2]), CellParams[2][2])
    # cVectReal3D = CellParams[2]

    #### ---- Map atoms onto the Supercell mesh and select
    # mapping matrix 2D: map(aVect, bVect) to real space with mC2R,use those as base, inverse
    # to get mapping matrix from real space to the (aVect, bVect) space.
    mR2SC = np.linalg.inv(np.array([aVect3Dreal,
                                    bVect3Dreal,
                                    cVect3Dreal]).T)
    # Select if 2D relative coordinates [0,0[
    AtomsInSuperCell = []  # Container, accepted atom coordinates (cartesian)
    NamesInSuperCell = []  # Container, accepted atom names
    for j, iAt in enumerate(AtomsInMesh):
        _iAt = mR2SC.dot(iAt)  # Map
        if min([round(k, 15) for k in _iAt]) < 0 or max([round(k, 15) for k in _iAt]) >= 1:
            continue
        else:
            debug(f"Included: {iAt} with {_iAt}")
            AtomsInSuperCell.append(iAt)  # almacena coord. cartesian
            NamesInSuperCell.append(NamesInMesh[j])  # almacena names

    # Multiply a few layers
    AtomsInSuperCellLayers = []  # Container, accepted atoms in multiple layers
    NamesInSuperCellLayers = []  # Container, accepted Names in multiple layers

    for iLy in range(nLayers):
        for j, iAt in enumerate(AtomsInSuperCell):
            AtomsInSuperCellLayers.append(iAt + iLy * cVect3Dreal)
            NamesInSuperCellLayers.append(NamesInSuperCell[j])
    # Fix cell height by nLayers and vacuum
    MaxHeigh = max([i[2] for i in AtomsInSuperCellLayers])
    VerticalAngle = MyAngle(CellParams[2], [0., 0., 1.], format="rad")
    UnitVertical = CellParams[2] / np.linalg.norm(CellParams[2])

    # In case vacuum zero, need to use CellParams[2]*nLayers to preserve periodicity (though makeZortho would break it)
    _cAddedVacc = UnitVertical * (MaxHeigh + vacuum) / np.cos(VerticalAngle)
    _cnLayers = CellParams[2] * nLayers
    if np.linalg.norm(_cAddedVacc) > np.linalg.norm(_cnLayers):
        cVect3DrealLayers = _cAddedVacc
    else:
        cVect3DrealLayers = _cnLayers

    #### Move to direct coordinates
    _, AtomsInSuperCellLayers_Direct = Cartesian2Direct([aVect3Dreal, bVect3Dreal, cVect3DrealLayers],
                                                        [[], AtomsInSuperCellLayers])

    #### Rotate base cells (sites in direct coordinates are ok)
    #### Rotate cell
    cellangle = MyAngle(aVect3Dreal[:2], [1., 0.], format='rad')
    Mrot = np.array([[np.cos(2 * np.pi - cellangle), -np.sin(2 * np.pi - cellangle), 0.],
                     [np.sin(2 * np.pi - cellangle), np.cos(2 * np.pi - cellangle), 0.],
                     [0., 0., 1.]])
    aVect3DrealRot = Mrot.dot(aVect3Dreal)
    bVect3DrealRot = Mrot.dot(bVect3Dreal)
    cVect3DrealLayersRot = Mrot.dot(cVect3DrealLayers)

    fLattice = [aVect3DrealRot, bVect3DrealRot, cVect3DrealLayersRot]

    #### Make orthonormal (if required)
    if MakeZortho:
        debug(
            f"Forcing lattice to be z-orthogonal, initial lattice: {fLattice}, sites: {AtomsInSuperCellLayers_Direct}")
        fLattice, _Sites = MakeCellZortho(fLattice, [[], AtomsInSuperCellLayers_Direct])
        AtomsInSuperCellLayers_Direct = _Sites[1]
    debug(f"Writing lattice:{fLattice}, sites:{AtomsInSuperCellLayers_Direct}")

    return {'lattice':fLattice,
            'Atoms':AtomsInSuperCellLayers_Direct,
            'Names':NamesInSuperCellLayers}


    #return POSCAR
