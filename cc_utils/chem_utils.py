__author__ = 'clyde'

def mol2_parse(fn,flag='converted'):
    """extracts atoms and bonds from a mol2 file, flag='converted' recovers the atom translations,
    flag='original' recovers the original atom types"""
    if flag=='converted':
        atom_ind = 5
    elif flag == 'original':
        atom_ind = 1
    else:
        raise RuntimeError('Invalid selection')

    with open(fn) as f:
        data = f.read()
        section_data = data.split('@')[1:]
        keys = [e.replace('@<TRIPOS>', '') for e in data.split() if '@' in e]
        master_dict = {k: section_data[i].split('\n')[1:-1] for i, k in enumerate(keys)}

        symbols = [l.split()[atom_ind] for l in master_dict.get('ATOM')]
        bonds = [l.split()[1:3] for l in master_dict.get('BOND')]

        return symbols, bonds

def mol2_replace(fn, new_fn, atom_types, flag='converted'):
    """replaces the atom types in a mol2 file with those provided to the function, flag='converted' replaces the atom translations
    flag='original' replaces the original atom types"""
    if flag=='converted':
        replace_ind = 5
    elif flag == 'original':
        replace_ind = 1
    else:
        raise RuntimeError('Invalid selection')

    with open(fn) as f:
        data = f.read()
        section_data = data.split('@')[1:]
        atom_section_index, atom_section = next((i,s.split('\n')) for i,s in enumerate(section_data) if '<TRIPOS>ATOM' in s.split('\n'))
        atom_header, atom_lines = atom_section[0], atom_section[1:-1]

        new_atom_lines = []
        assert len(atom_lines) == len(atom_types)

        for i,line in enumerate(atom_lines):
            line_contents = line.split()
            line_contents[replace_ind] = atom_types[i]
#            new_line = ' '.join(line_contents)
            new_line = '{:>7} {:<3} {:>16} {:>9} {:>9} {:<3} {:>8} {:<2} {:>12}'.format(*line_contents)
            new_atom_lines.append(new_line)

        new_atom_section = "\n".join([atom_header] + new_atom_lines) + '\n'
        new_section_data = section_data[0:atom_section_index] + [new_atom_section] + section_data[atom_section_index+1:]
        new_data = '@' + '@'.join(new_section_data)

    with open(new_fn, 'w') as f:
        f.write(new_data)

#todo
def atom_type_map(data, amber_mol_f, gaff_mol_f):
    with open(amber_mol_f) as amber_f:
        amber_data = amber_f.readlines()
    with open(gaff_mol_f) as gaff_f:
        gaff_data = gaff_f.readlines()

    cleaned_amber = [e.strip() for e in data[2:-1]]
    cleaned_gaff = [e.strip() for e in data[2:-1]]
    amber_atoms = [e.split()[5] for e in cleaned_amber if len(e.split()) == 9]
    gaff_atoms = [e.split()[5] for e in cleaned_gaff if len(e.split()) == 9]

    atom_dict = {}
    alt_dict = {}

    for i in range(len(amber_atoms)):
        if amber_atoms[i] not in atom_dict:
            atom_dict.update({amber_atoms[i]:gaff_atoms[i]})

        elif atom_dict[amber_atoms[i]] != gaff_atoms[i] and amber_atoms[i] in alt_dict and alt_dict[amber_atoms[i]] != gaff_atoms[i]:
            alt_dict[amber_atoms[i]].append(gaff_atoms[i])

        elif atom_dict[amber_atoms[i]] != gaff_atoms[i] and not amber_atoms[i] in alt_dict:
            alt_dict.update({amber_atoms[i]:[gaff_atoms[i]]})

    return atom_dict

def set_frag_atom_types(fragment_mol2, new_fragment_mol2, full_pdb, atom_types='amber'):
    import os
    import Biskit as B
    fragment_orig_atoms, bonds = mol2_parse(fragment_mol2, flag='original')

    m = B.PDBModel(full_pdb)
    orig_atoms = [a['name'] for a in m]

    fragment_indices = None
    for i,n in enumerate(orig_atoms):
        try:
            if orig_atoms[i:i+len(fragment_orig_atoms)] == fragment_orig_atoms:
                fragment_indices = range(i, i+len(fragment_orig_atoms))
                continue
        except IndexError:
            pass

    if not fragment_indices:
        raise RuntimeError('Failed to find fragment')

    if atom_types == 'amber':
        os.system('antechamber -fi pdb -i {p} -fo mol2 -o {m} -at amber'.format(p=full_pdb, m=full_pdb.replace('.pdb','_amber.mol2')))
        amber_atoms, bonds = mol2_parse(full_pdb.replace('.pdb', '_amber.mol2'))
        fragment_atoms = [amber_atoms[i] for i in fragment_indices]
    elif atom_types == 'gaff':
        os.system('antechamber -fi pdb -i {p} -fo mol2 -o {m}'.format(p=full_pdb, m=full_pdb.replace('.pdb','_gaff.mol2')))
        gaff_atoms, bonds = mol2_parse(full_pdb.replace('.pdb', '_gaff.mol2'))
        fragment_atoms = [gaff_atoms[i] for i in fragment_indices]
    else:
        raise RuntimeError('Invalid atom_type selected')

    mol2_replace(fragment_mol2, new_fragment_mol2, fragment_atoms)

def set_nonstd_name(fragment_mol2, res_name):
    """Sets the same of the residue from LIG to the string specified by res_name"""
    with open(fragment_mol2) as mol2_f:
        content = mol2_f.read()

    with open(fragment_mol2, 'w') as mol2_f:
        mol2_f.write(content.replace('LIG', res_name))


def get_bad(mol, gaff_atoms, amber_atoms, non_std_atom_indexes):
    """Computes the atoms involved in bonds, angles and dihedrals spanning both
    the non standard residue specified and and the rest of the protein.
    Returns mixed amber-gaff, gaff-gaff and amber-amber combinations"""

    from cc_notebook_utils import unwind
    std_nn_ind = unwind([[n for n in mol.graph.neighbors[i] if n not in non_std_atom_indexes] for i in non_std_atom_indexes])
    nstd_nn_ind = [[n for n in mol.graph.neighbors[i] if n in non_std_atom_indexes] for i in std_nn_ind]
    std_nnn_ind = [[n for n in mol.graph.neighbors[i] if n not in non_std_atom_indexes] for i in std_nn_ind]
    nstd_nnn_ind = [[n for n in mol.graph.neighbors[i] if n in non_std_atom_indexes] for i in unwind(nstd_nn_ind)]

    amb_std_nn = [amber_atoms[i] for i in std_nn_ind]
    amb_std_nnn = [[amber_atoms[i] for i in nnn_ind] for nnn_ind in std_nnn_ind]
    gaff_std_nn = [gaff_atoms[i] for i in std_nn_ind]
    gaff_std_nnn = [[gaff_atoms[i] for i in nnn_ind] for nnn_ind in std_nnn_ind]

    amb_nstd_nn = [[amber_atoms[i] for i in nn_ind] for nn_ind in nstd_nn_ind]
    amb_nstd_nnn = [[amber_atoms[i] for i in nnn_ind] for nnn_ind in nstd_nnn_ind]
    gaff_nstd_nn = [[gaff_atoms[i] for i in nn_ind] for nn_ind in nstd_nn_ind]
    gaff_nstd_nnn = [[gaff_atoms[i] for i in nnn_ind] for nnn_ind in nstd_nnn_ind]

    bonds, gbonds, abonds = [], [], []
    for i in range(len(amb_std_nn)):
        for k in range(len(gaff_nstd_nn[i])):
            bonds.append([amb_std_nn[i], gaff_nstd_nn[i][k]])
            gbonds.append([gaff_std_nn[i], gaff_nstd_nn[i][k]])
            abonds.append([amb_std_nn[i], amb_nstd_nn[i][k]])

    dihedrals, gdihedrals, adihedrals = [], [], []
    for i in range(len(amb_std_nn)):
        for j in range(len(amb_std_nnn[i])):
            for k in range(len(gaff_nstd_nn[i])):
                for l in range(len(gaff_nstd_nnn[i])):
                    dihedrals.append([amb_std_nnn[i][j], amb_std_nn[i], gaff_nstd_nn[i][k], gaff_nstd_nnn[i][l]])
                    gdihedrals.append([gaff_std_nnn[i][j], gaff_std_nn[i], gaff_nstd_nn[i][k], gaff_nstd_nnn[i][l]])
                    adihedrals.append([amb_std_nnn[i][j], amb_std_nn[i], amb_nstd_nn[i][k], amb_nstd_nnn[i][l]])

    angles, gangles, aangles = [], [], []
    for i in range(len(dihedrals)):
        ang1, ang2 = dihedrals[i][0:-1], dihedrals[i][1:]
        if ang1 not in angles:
            angles.append(ang1)
            gangles.append(gdihedrals[i][0:-1])
            aangles.append(adihedrals[i][0:-1])
        if ang2 not in angles:
            angles.append(ang2)
            gangles.append(gdihedrals[i][1:])
            aangles.append(adihedrals[i][1:])

    return {'mixed':[bonds, angles, dihedrals], 'gaff': [gbonds, gangles, gdihedrals] , 'amber': [abonds, aangles, adihedrals]}

def get_ht_ind(mol, non_std_atom_indexes):
    """returns the head and the tail index of the nonstandard residue specified"""
    from cc_notebook_utils import unwind
    std_nn_ind = unwind([[n for n in mol.graph.neighbors[i] if n not in non_std_atom_indexes] for i in non_std_atom_indexes])
    nstd_nn_ind = [[n for n in mol.graph.neighbors[i] if n in non_std_atom_indexes] for i in std_nn_ind]
    return unwind(nstd_nn_ind) - non_std_atom_indexes[0]

def get_amber_params_from_xml(vers='parm99'):
    """returns a dictionary of parameters corresponding to the amber forcefield specified.
    Parameters are taken from an xml file within the ambertools mtkpp module.
    WARNING: These parameters have errors in them! Use get_amber_params_from_dat instead"""
    from lxml import etree
    import os

    amb_path = os.environ['AMBERHOME']
    xml_param_path = amb_path + '/dat/mtkpp'

    try:
        params = etree.parse(xml_param_path + '/{v}.xml'.format(v=vers))
    except IOError:
        raise IOError('Parameter xml file {f} missing from {p}'.format(f=vers+'.xml', p=xml_param_path))

    root = params.getroot()
    content = root.getchildren()

    final_dict = {}
    for info in content:
        key = info.tag
        if key != 'equivalentAtoms':
            values = [dict(e.items()) for e in info.getchildren() if e.items()]
            final_dict.update({key:values})
        else:
            equiv_atoms = info.getchildren()
            atom_keys = [orig.items()[0][1] for orig in equiv_atoms]
            atom_values = [[equiv.items()[0][1] for equiv in orig.getchildren()] for orig in equiv_atoms]
            equiv_dict = dict([(e, atom_keys[i]) for i in range(len(atom_keys)) for e in atom_values[i]])

    #sort out the vdw entries corresponding to equivalent atom types by pasting in the suitable parameter values
    vdw_keys = sorted(list(set([k for e in final_dict.get('types', []) for k in e.keys()])))
    missing_vdw_entries = [entry for entry in final_dict.get('types', []) if entry and sorted(entry.keys()) != vdw_keys]

    for missing_entry in missing_vdw_entries:
        missing_keys = [k for k in vdw_keys if k not in missing_entry.keys()]
        #if this next(...) statement breaks it means we have missing entries but no equivalent atoms defined
        equiv_entry = next(entry for entry in final_dict.get('types', []) if entry['name'] == equiv_dict[missing_entry['name']])
        missing_entry.update({k: equiv_entry[k] for k in missing_keys})

    return {k:final_dict[k] for k in final_dict if final_dict[k]}


#taken from acpype: https://ccpn.svn.sourceforge.net/svnroot/ccpn/branches/stable/ccpn/python/acpype/acpype.py
def get_amber_params_from_frcmod(vers='ff99SB'):
    """returns a dictionary of parameters corresponding to the amber parameters specified.
     Parameters are taken from a frcmod file generated by ambertools tleap program."""
    from collections import OrderedDict
    import os

    if not os.path.isabs(vers):
        amb_path = os.environ['AMBERHOME']
        dat_param_path = amb_path + '/dat/leap/parm'
        frcmod_fn = 'frcmod.' + vers
    else:
        dat_param_path, frcmod_fn = os.path.split(vers)

    try:
        with open(dat_param_path + '/' + frcmod_fn) as f:
            lista = f.readlines()
    except IOError:
        raise IOError('Parameter file {f} missing from {p}'.format(f=vers+'.frcmod', p=dat_param_path))

    heads = ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'HBON', 'NONB']
    d = OrderedDict()
    for line in lista[1:]:
        line = line.strip()
        if line[:4] in heads:
            head = line[:4]
            d[head] = []
            dd = OrderedDict()
            continue
        elif line:
            key = line.replace(' -', '-').replace('- ', '-').split()[0]
            if key in dd:
                if not dd[key].count(line):
                    dd[key].append(line)
            else:
                dd[key] = [line]
            d[head] = dd

    for k in d.keys():
        if not d[k]: d.pop(k)

    final_dict = {'types': [],'bondLengths': [], 'bondAngles': [], 'bondTorsions': [], 'bondImpropers': []}

    for types_k in d.get('NONB',[]):
        for param_str in d['NONB'][types_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            element, vdwRadius, potentialWellDepth = param_str.split()[0:3]
            description = " ".join(param_str.split()[3:])
            mass_str = d['MASS'][types_k][0].replace(' -', '-').replace('- ', '-')
            mass = mass_str.split()[1]
            final_dict['types'].append({'description': description, 'element': element, 'hybridization': 'unknown', 'mass': mass,
                                        'name': element, 'potentialWellDepth': potentialWellDepth, 'vdwRadius': vdwRadius})

    for length_k in d.get('BOND',[]):
        for param_str in d['BOND'][length_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            t1,t2 = length_k.split('-')
            keq,req = param_str.split()[1:3]
            final_dict['bondLengths'].append({'t1': t1, 't2': t2, 'keq': keq, 'req': req})

    for angle_k in d.get('ANGL',[]):
        for param_str in d['ANGL'][angle_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            t1,t2,t3 = angle_k.split('-')
            keq, req = param_str.split()[1:3]
            final_dict['bondAngles'].append({'t1': t1, 't2': t2, 't3': t3, 'keq': keq, 'req': req})

    for dihedral_k in d.get('DIHE',[]):
        for param_str in d['DIHE'][dihedral_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            t1,t2,t3,t4 = dihedral_k.split('-')
            npth, Vn, gamma, Nt = param_str.split()[1:5]
            final_dict['bondTorsions'].append({'t1': t1, 't2': t2, 't3': t3, 't4': t4, 'npth': npth, 'Vn': Vn, 'gamma': gamma, 'Nt': Nt})

    for improper_k in d.get('IMPR',[]):
        for param_str in d['IMPR'][improper_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            t1,t2,t3,t4 = improper_k.split('-')
            Vn, gamma, Nt = param_str.split()[1:4]
            final_dict['bondImpropers'].append({'t1': t1, 't2': t2, 't3': t3, 't4': t4, 'Vn': Vn, 'gamma': gamma, 'Nt': Nt})

    return final_dict


def get_amber_params_from_dat(vers='parm99'):
    """returns a dictionary of parameters corresponding to the amber parameters specified.
     Parameters are taken from the main amber dat file."""
    import os
    from collections import OrderedDict
    amb_path = os.environ['AMBERHOME']
    dat_param_path = amb_path + '/dat/leap/parm'

    try:
        with open(dat_param_path + '/' + vers + '.dat') as f:
            dat_data = f.read()
    except IOError:
        raise IOError('Parameter file {f} missing from {p}'.format(f=vers+'.dat', p=dat_param_path))


    mass, blengths, bangs, btorsions, bimpropers = [s.split('\n') for s in dat_data.split('\nEND')[0].split('\n\n')[0:5]]
    mass = mass[1:]
    blengths = blengths[1:]

    nonbnds = dat_data.split('\nEND')[0].split('MOD4')[1].split('\n\n')[0].split('\n')[1:-1]

    init_data = {'MASS': mass, 'BOND': blengths, 'ANGL': bangs, 'DIHE': btorsions, 'IMPR': bimpropers, 'NONB': nonbnds}
    d= OrderedDict()

    for k in init_data:
        head = k
        d[head] = []
        dd = OrderedDict()

        for line in init_data[head]:
            key = line.replace(' -', '-').replace('- ', '-').split()[0]

            if key in dd:
                if not dd[key].count(line):
                    dd[key].append(line)
            else:
                dd[key] = [line]
            d[head] = dd

    for k in d.keys():
        if not d[k]: d.pop(k)

    final_dict = {'types': [], 'bondLengths': [], 'bondAngles': [], 'bondTorsions': [], 'bondImpropers': []}
    for types_k in d['NONB']:
        for param_str in d['NONB'][types_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            element, vdwRadius, potentialWellDepth = param_str.split()[0:3]
            description = " ".join(param_str.split()[3:])
            mass_str = d['MASS'][types_k][0].replace(' -', '-').replace('- ', '-')
            mass = mass_str.split()[1]
            final_dict['types'].append({'description': description, 'element': element, 'hybridization': 'unknown', 'mass': mass,
                                        'name': element, 'potentialWellDepth': potentialWellDepth, 'vdwRadius': vdwRadius})

    for length_k in d['BOND']:
        for param_str in d['BOND'][length_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            t1,t2 = length_k.split('-')
            keq,req = param_str.split()[1:3]
            final_dict['bondLengths'].append({'t1': t1, 't2': t2, 'keq': keq, 'req': req})

    for angle_k in d['ANGL']:
        for param_str in d['ANGL'][angle_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            t1,t2,t3 = angle_k.split('-')
            keq, req = param_str.split()[1:3]
            final_dict['bondAngles'].append({'t1': t1, 't2': t2, 't3': t3, 'keq': keq, 'req': req})

    for dihedral_k in d['DIHE']:
        for param_str in d['DIHE'][dihedral_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            t1,t2,t3,t4 = dihedral_k.split('-')
            try:
                npth, Vn, gamma, Nt = param_str.split()[1:5]
            except ValueError:
                print(param_str)
            final_dict['bondTorsions'].append({'t1': t1, 't2': t2, 't3': t3, 't4': t4, 'npth': npth, 'Vn': Vn, 'gamma': gamma, 'Nt': Nt})

    for improper_k in d['IMPR']:
        for param_str in d['IMPR'][improper_k]:
            param_str = param_str.split('  ')[0].replace(' -', '-').replace('- ', '-') + '  ' + '  '.join(param_str.split('  ')[1:])
            t1,t2,t3,t4 = improper_k.split('-')
            Vn, gamma, Nt = param_str.split()[1:4]
            final_dict['bondImpropers'].append({'t1': t1, 't2': t2, 't3': t3, 't4': t4, 'Vn': Vn, 'gamma': gamma, 'Nt': Nt})

    return final_dict

#see p61 of the Ambertools13 manual
def sum_dihedrals(unsummed_dihedrals):
    """combines several torsion parameter values for corresponding to a particular atom combination together
    For torsions with multiple harmonic potentials with different periodicities the gaussian input format and
    the amber format differ so here we convert from amber_amber to gaussian_amber"""

    import copy
    new_dihedrals = copy.deepcopy(unsummed_dihedrals)
    # Nt gives the periodicity, in the amber files it can be negative or positive,
    # negative indicates there are other harmonic potentials associated with the same atoms
    for torsion in new_dihedrals:
        comp=abs(int(torsion.pop('Nt').split('.')[0]))
        torsion.update({'gamma_{c}'.format(c=comp):torsion.pop('gamma')})
        torsion.update({'Vn_{c}'.format(c=comp):torsion.pop('Vn')})

    unique_dihedrals = set([(t['t1'],t['t2'],t['t3'],t['t4']) for t in new_dihedrals])
    for unique_dihedral in unique_dihedrals:
        torsion = next(t for t in new_dihedrals if t and (t['t1'], t['t2'], t['t3'], t['t4']) == unique_dihedral)
        torsion_inds = [i for i,t in enumerate(new_dihedrals) if t and (t['t1'], t['t2'], t['t3'], t['t4']) == unique_dihedral]

        for torsion_ind in torsion_inds:
            torsion_2 = new_dihedrals[torsion_ind]
            if (torsion['t1'], torsion['t2'], torsion['t3'], torsion['t4']) == (torsion_2['t1'], torsion_2['t2'], torsion_2['t3'], torsion_2['t4']) and torsion is not torsion_2:
                torsion.update(torsion_2)
                new_dihedrals[torsion_ind] = None
    return [t for t in new_dihedrals if t]


def gen_dihedral_str(torsion):
    """generates a gaussian readable string from a torsion parameter dictionary"""
    if torsion['t1'] == 'X':
        torsion['t1'] = '*'
    if torsion['t2'] == 'X':
        torsion['t2'] = '*'
    if torsion['t3'] == 'X':
        torsion['t3'] = '*'
    if torsion['t4'] == 'X':
        torsion['t4'] = '*'

    torsion['Vn_1'] = torsion.get('Vn_1','0.00')
    torsion['Vn_2'] = torsion.get('Vn_2','0.00')
    torsion['Vn_3'] = torsion.get('Vn_3','0.00')
    torsion['Vn_4'] = torsion.get('Vn_4','0.00')
    torsion['gamma_1'] = str(int(float(torsion.get('gamma_1','0'))))
    torsion['gamma_2'] = str(int(float(torsion.get('gamma_2','0'))))
    torsion['gamma_3'] = str(int(float(torsion.get('gamma_3','0'))))
    torsion['gamma_4'] = str(int(float(torsion.get('gamma_4','0'))))


    return 'AmbTrs {a1:<2} {a2:<3} {a3:<3} {a4:<3} {g1:<5} {g2:<5} {g3:<5} {g4:<5} {v1:<6} {v2:<6} {v3:<6} {v4:<6} {p}.0'.format(a1=torsion['t1'], a2=torsion['t2'], a3=torsion['t3'], a4=torsion['t4'],
                                                                                                                                 g1=torsion['gamma_1'], g2=torsion['gamma_2'], g3=torsion['gamma_3'], g4=torsion['gamma_4'],
                                                                                                                                 v1=torsion['Vn_1'], v2=torsion['Vn_2'], v3=torsion['Vn_3'], v4=torsion['Vn_4'],
                                                                                                                                 p=torsion['npth'])


def gen_improper_str(improper):
    """generates a gaussian readable string from an improper torsion parameter dictionary"""
    if improper['t1'] == 'X':
        improper['t1'] = '*'
    if improper['t2'] == 'X':
        improper['t2'] = '*'
    if improper['t3'] == 'X':
        improper['t3'] = '*'
    if improper['t4'] == 'X':
        improper['t4'] = '*'

    return 'ImpTrs {a1:<2} {a2:<2} {a3:<3} {a4:<3} {v:<6} {g:<5} {n}'.format(a1=improper['t1'], a2=improper['t2'], a3=improper['t3'], a4=improper['t4'],
                                                                             g=improper['gamma'], v=improper['Vn'], n=improper['Nt'])


def write_gauss_amber_params(params, para_file=None):
    """genererates gaussian readable forcefield file from an amber parameter dictionary"""
    nbnds = ['NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.000']

    stretches = []
    for stretch in params.get('bondLengths', []):
        stretches.append('HrmStr1 {a1:<2} {a2:<3} {k:<6} {bl}'.format(a1=stretch['t1'], a2=stretch['t2'], k=stretch['keq'], bl=stretch['req']))

    angles = []
    for angle in params.get('bondAngles', []):
        angles.append('HrmBnd1 {a1:<2} {a2:<3} {a3:<3} {k:<6} {ang}'.format(a1= angle['t1'], a2 = angle['t2'], a3 = angle['t3'], k=angle['keq'], ang=angle['req']))

    torsions = []
    combined_torsions = sum_dihedrals(params.get('bondTorsions', []))
    for torsion in combined_torsions:
        torsions.append(gen_dihedral_str(torsion))

    impropers = []
    for improper in params.get('bondImpropers', []):
        impropers.append(gen_improper_str(improper))

    vdws = []
    for vdw in params.get('types', []):
        vdws.append('VDW {a:<2} {r:<6} {d}'.format(a=vdw['name'], r=vdw['vdwRadius'], d=vdw['potentialWellDepth']))

    gauss_amb = '\n\n'.join(['\n'.join(nbnds), '\n'.join(stretches), '\n'.join(angles), '\n'.join(torsions), '\n'.join(impropers), '\n'.join(vdws)])

    if para_file:
        with open(para_file, 'w') as f:
            f.write(gauss_amb)
    else:
        return gauss_amb


#todo worry about reversed order parameters that are infact equivalent e.g. t1,t2,t3,t4 == t4,t3,t2,t1
def update_amber_params(orig_params, mod_params, overwrite=True):
    """updates a dictionary of amber parameters with the values present in a second dictionary of amber parameters.
     This is useful for modifications of a parameter set, e.g. for ff99 -> ff99SB we load parm99 then update it with frcmod.ff99SB"""
    import copy
    final_params = copy.deepcopy(orig_params)

    for k in mod_params.keys():
        for mod_p in mod_params[k]:
            # for each mod param find equivalent orig param and if we are overwriting update it,
            # if no equivalent orig param found just add mod param to the new param stack
            try:
                replace_p = next(org_p for org_p in final_params[k] if (abs(int(float(org_p.get('Nt','0')))), org_p.get('element'), org_p.get('t1'), org_p.get('t2'), org_p.get('t3'), org_p.get('t4'))
                                                                    == (abs(int(float(mod_p.get('Nt','0')))), mod_p.get('element'), mod_p.get('t1'), mod_p.get('t2'), mod_p.get('t3'), mod_p.get('t4')))
                if overwrite:
                    replace_p.update(mod_p)
            except StopIteration:
                final_params[k].append(mod_p)

    return final_params
