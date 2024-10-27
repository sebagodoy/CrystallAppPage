import numpy as np

def NPcontainer(siteCoords):
    allX = [i[0] for i in siteCoords]
    allY = [i[1] for i in siteCoords]
    allZ = [i[2] for i in siteCoords]
    #print(allX)
    #print(allY)
    #print(allZ)
    return [abs(max(i))+abs(min(i)) for i in [allX, allY, allZ]]

def CreateNP_C5(dRef, rings, heigh, marks,
                element='X',
                facet_extension = 1., column_extension = 1.):
    print("Creating NP-C5")
    #[print(i) for i in [dRef, rings, heigh, marks,element,
    #            facet_extension,column_extension]]
    # helper functions
    def CilToCart(vectCil):
        return np.array([vectCil[0] * np.cos(vectCil[1]), vectCil[0] * np.sin(vectCil[1]), vectCil[2]])
    def CartToCil(vectCil):
        return np.array([np.linalg.norm([vectCil[0], vectCil[1], 0]), np.arctan(vectCil[1] / vectCil[0]), vectCil[2]])

    def GenerateColumn(ir, ith, iz0, n, id):
        # Cilindrical system
        # ir, ith column coodr, iz0 cilindrical sistem
        # id: atom distance
        return [np.array([ir, ith, iz0 + k * id]) for k in range(n)]

    #### ---- Creating NP sites

    #Axis = 10
    #dAA = 2.39774  # 2.462

    #column_length = 5  # extend the column to Ino's shappe

    #truncated_corner = 3
    #zurco_width = 0

    AtomCollection = []

    # -- Central column
    [AtomCollection.append(i) for i in GenerateColumn(0, 0, 0, rings + int(heigh), dRef)]

    # -- Rings
    AtomCollectionRadial = []
    for iRing in range(1, rings):
        #### ---- Direct distribution
        # Two point in ring, separate
        oTheta = (2 * np.pi / 5)
        aPoint = CilToCart(np.array([iRing * dRef * facet_extension * np.sqrt(3) / 2, 0, iRing * dRef / 2]))
        bPoint = CilToCart(np.array([iRing * dRef * facet_extension * np.sqrt(3) / 2, oTheta, iRing * dRef / 2]))
        abPoint = bPoint - aPoint
        abPointMag = np.linalg.norm(abPoint)
        abPointNorm = abPoint * (1 / abPointMag)

        for jAt in range(iRing):
            Position = aPoint + abPointNorm * jAt * abPointMag / iRing
            PositionCil = CartToCil(Position)

            if iRing + marks >= rings:
                if jAt + rings - iRing - 1 < marks:
                    continue
                if rings - jAt - 1 < marks:
                    continue

            [AtomCollectionRadial.append(i)
             for i in GenerateColumn(PositionCil[0],  # corrected distance from origin
                                     PositionCil[1],  # angle from origin
                                     PositionCil[2],  # z displacement
                                     rings - iRing + int(heigh),  # How many atoms in this column
                                     dRef * column_extension  # Atomic distance (column)
                                     )
             ]
        pass
    # First slice
    [AtomCollection.append(i) for i in AtomCollectionRadial]
    # Other slices, C5 symmetry
    SymmAng = 2 * np.pi / 5
    for iSym in range(1, 5):
        [AtomCollection.append(i + np.array([0, iSym * SymmAng, 0])) for i in AtomCollectionRadial]

    AtomCollectionCart = [np.array(CilToCart(i)) for i in AtomCollection]

    return AtomCollectionCart