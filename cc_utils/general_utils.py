__author__ = 'clyde'

import math
import numpy as np
import pandas


def norm(vector):
    """ Returns the norm (length) of the vector."""
    # note: this is a very hot function, hence the odd optimization
    # Unoptimized it is: return np.sqrt(np.sum(np.square(vector)))
    return np.sqrt(np.dot(vector, vector))

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if math.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

def gen_energy_table(products, reactants, delta=None):
    from ase import atoms
    rxn_Es = []
    names = []

    if delta:
        delt_Es = []

    for i in range(len(products)):
        if isinstance(products[i], atoms.Atoms):
            main_product = products[i]
            try:
                product_E = products[i].calc.energy_zero
            except AttributeError:
                product_E = float('nan')
        else:
            main_product = products[i][0]
            try:
                product_E = sum([p.calc.energy_zero for p in products[i]])
            except AttributeError:
                product_E = float('nan')

        if isinstance(reactants[i], atoms.Atoms):
            try:
                reactant_E = reactants[i].calc.energy_zero
            except AttributeError:
                reactant_E = float('nan')
        else:
            try:
                reactant_E = sum([r.calc.energy_zero for r in reactants[i]])
            except AttributeError:
                reactant_E = float('nan')

        rxn_E = 23.060542301388647 * (product_E - reactant_E)

        if delta:
            delt_E = rxn_E - delta
            delt_Es.append(delt_E)

        rxn_Es.append(rxn_E)
        names.append(main_product.calc.label)

    o_data = pandas.Series(rxn_Es, names)
    d = {'Rxn Energy': o_data}

    if delta:
        o_dft_delt_data = pandas.Series(delt_Es, names)
        d.update({'Rxn Energy Delta': o_dft_delt_data})
    return pandas.DataFrame(d)

