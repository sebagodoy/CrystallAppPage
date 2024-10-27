#### ---- export XYZ

def write_XYZ_file(AtomCollection, OutFile = 'Out.xyz', msg = "Geometry output", Element = 'X'):
    # AtomCollection: list atommic coordinates, cartessian [ [], [], ... ]
    with open(OutFile, 'w') as f:
        f.write(f"{len(AtomCollection)}\n")
        f.write(msg + "\n")
        [ f.write(Element + ' '.join(["{:.6f}".format(k).rjust(11) for k in i]) + "\n") for i in AtomCollection]

def write_XYZ(AtomCollection, msg = "Geometry output", Elements = 'X'):
    # AtomCollection: list atommic coordinates, cartessian [ [], [], ... ]
    print("WRITTING XYZ")
    #print(msg)
    #print(Elements)
    OutFile = ""
    OutFile += f"{len(AtomCollection)}\n"
    OutFile += msg + "\n"
    for iElement, iPosition in zip(Elements, AtomCollection):
        #print(f'STEP: {iElement}, {iPosition}')
        #OutFile += iElement + ' '.join(["{:.6f}".format(k).rjust(11) for k in AtomCollection]) + "\n"
        OutFile += iElement + ' '.join(["{:.6f}".format(k).rjust(11) for k in iPosition]) + "\n"
    return OutFile

def write_POSCAR(AtomCollection ,Lattice, Names, Coordinates="Cartesian"):
    '''Writes POSCAR type file from cell and sites information
    Input:  Lattice=[[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
            AtomCollection= [[x1,y1,z1], [...], ... ]
            Names= [Xx, Y, Z, ... ]
    Output: POSCAR string'''
    print(AtomCollection)
    print("-"*80)
    print(Lattice)
    print("-"*80)
    print(Names)
    # start POSCAR string
    POSCAR = "Supercell \n" + "{:.16f}".format(1.0).rjust(20) + "\n"
    # write cell
    print(f"LATTICE : {Lattice}")
    for iCoord in Lattice:
        POSCAR += " ".join(["{:.15f}".format(round(i, 16)).rjust(18) for i in iCoord]) + "\n"
    # Sort by atom type
    Types = {i: [] for i in Names}
    for idx, iType in enumerate(Names):
        Types[iType].append(AtomCollection[idx])
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


