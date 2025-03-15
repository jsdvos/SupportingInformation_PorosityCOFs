import os

import numpy as np
np.random.seed(1912180624)

from yaff import System
from molmod.units import angstrom, deg

def get_indices(n):
    high = 96
    for i in np.random.choice(high, size = n, replace = False):
        yield i

sys = System.from_file('PA_TP.chk') # Unphysical interlayer distance of 10A
rvecs = sys.cell.rvecs.copy()
rvecs[2][2] = 3.0*angstrom
sys.cell.update_rvecs(rvecs)
assert sys.natom == 72
sys = sys.supercell(2,2,1)

# Create system with six different offsets
sys6 = sys.supercell(1,1,6)
d = np.array([0.0, 0.0, 0.0])
for i in range(6):
    sys6.pos[i*288:(i+1)*288:] += d
    angle = (30+i*60)*deg
    d += np.array([np.cos(angle), np.sin(angle), 0.0])*2.5*angstrom

# System 1: ABCD.. inclination
rvecs = sys.cell.rvecs.copy()
rvecs[2] = np.array([0.0, 2.5*angstrom, 3.0*angstrom])
sys.cell.update_rvecs(rvecs)
abcdef = sys6.supercell(1,1,2) # ABCDEF inclination

# Convert TP to TFB
def iter_todo():
    def indices(n):
        result = []
        high = 96
        for i in np.random.choice(high, size = n, replace = False):
            result.extend([x for x in indices_tp(i)])
        return result
    def indices_tp(n):
        yield n
    
    # Mixed COFs
    for i in [10, 19, 29, 38, 48, 58, 67, 77, 86]:
        yield 'full{}'.format(i), indices(i)
    
    # TFB-PA COF
    yield 'all', [i for i in range(96)]
    # TP-PA COF
    yield 'none', []

def convert_system(sys, indices):
    def iter_tp(index):
        if index % 2 == 0:
            for i in range(15):
                yield 72*(index//2) + i
        elif index % 2 == 1:
            for i in range(57, 72):
                yield 72*(index//2) + i

    def iter_pa(index):
        if index % 3 == 0:
            for i in range(15, 29):
                yield 72*(index//3) + i
        elif index % 3 == 1:
            for i in range(29, 43):
                yield 72*(index//3) + i
        elif index % 3 == 2:
            for i in range(43, 57):
                yield 72*(index//3) + i

    def get_neightype(i):
        for j in sys.neighs1[i]:
            if not ffatypes[j].endswith('PA11'):
                return {'TP': 1, 'TPB': 2}[ffatypes[j].split('_')[-1]]

    bonds = sys.bonds.copy()
    pos = sys.pos.copy()
    rvecs = sys.cell.rvecs.copy()
    ffatypes = [sys.get_ffatype(i) for i in range(sys.natom)]
    numbers = sys.numbers.copy()
    # Convert TP to TFB
    for index in indices:
        for count, i in enumerate(iter_tp(index)):
            ffatype = ffatypes[i]
            assert ffatype.endswith('_TP')
            if ffatype == 'C_C3_TP':
                ffatypes[i] = 'C_C3_H3C2N_TPB'
            elif ffatype == 'C_C2O_TP':
                ffatypes[i] = 'C_HC2_C4_TPB'
            elif ffatype == 'O_TP':
                numbers[i] = 1
                ffatypes[i] = 'H_C_C2_TPB'
            elif ffatype == 'C_HCN_TP':
                ffatypes[i] = 'C_HCN_C3_TPB'
            elif ffatype == 'H_C_CN_TP':
                ffatypes[i] = 'H_C_CN_TPB'
            else:
                raise RuntimeError('Did not expect ffatype ' + ffatype)
            assert not ffatypes[i] == ffatype

    # Update PA11
    for index in range(144):
        new_indices = {}
        for count, i in enumerate(iter_pa(index)):
            new_indices[count] = i
        target = 'PA{}{}'.format(get_neightype(new_indices[10]), get_neightype(new_indices[11]))
        for count, i in new_indices.items():
            ffatypes[i] = {
                    0: {
                        'PA11': 'C_C2N_0_PA11',
                        'PA12': 'C_C2N_H3C3_0_PA12',
                        'PA21': 'C_C2N_H2C3_PA12',
                        'PA22': 'C_C2N_PA22'
                        },
                    1: {
                        'PA11': 'C_C2N_0_PA11',
                        'PA12': 'C_C2N_H2C3_PA12',
                        'PA21': 'C_C2N_H3C3_0_PA12',
                        'PA22': 'C_C2N_PA22'
                        },
                    2: {
                        'PA11': 'C_HC2_HC2N_0_PA11',
                        'PA12': 'C_HC2_HC2N_1_PA12',
                        'PA21': 'C_HC2_HC2N_0_PA12',
                        'PA22': 'C_HC2_PA22'
                        },
                    3: {
                        'PA11': 'C_HC2_HC2N_0_PA11',
                        'PA12': 'C_HC2_HC2N_0_PA12',
                        'PA21': 'C_HC2_HC2N_1_PA12',
                        'PA22': 'C_HC2_PA22'
                        },
                    4: {
                        'PA11': 'C_HC2_HC2N_0_PA11',
                        'PA12': 'C_HC2_HC2N_0_PA12',
                        'PA21': 'C_HC2_HC2N_1_PA12',
                        'PA22': 'C_HC2_PA22'
                        },
                    5: {
                        'PA11': 'C_HC2_HC2N_0_PA11',
                        'PA12': 'C_HC2_HC2N_1_PA12',
                        'PA21': 'C_HC2_HC2N_0_PA12',
                        'PA22': 'C_HC2_PA22'
                        },
                    6: {
                        'PA11': 'H_C_C2_0_PA11',
                        'PA12': 'H_C_C2_1_PA12',
                        'PA21': 'H_C_C2_0_PA12',
                        'PA22': 'H_C_PA22'
                        },
                    7: {
                        'PA11': 'H_C_C2_0_PA11',
                        'PA12': 'H_C_C2_0_PA12',
                        'PA21': 'H_C_C2_1_PA12',
                        'PA22': 'H_C_PA22'
                        },
                    8: {
                        'PA11': 'H_C_C2_0_PA11',
                        'PA12': 'H_C_C2_0_PA12',
                        'PA21': 'H_C_C2_1_PA12',
                        'PA22': 'H_C_PA22'
                        },
                    9: {
                        'PA11': 'H_C_C2_0_PA11',
                        'PA12': 'H_C_C2_1_PA12',
                        'PA21': 'H_C_C2_0_PA12',
                        'PA22': 'H_C_PA22'
                        },
                    10: {
                        'PA11': 'N_0_PA11',
                        'PA12': 'N_HC2_0_PA12',
                        'PA21': 'N_C2_PA12',
                        'PA22': 'N_C2_PA22'
                        },
                    11: {
                        'PA11': 'N_0_PA11',
                        'PA12': 'N_C2_PA12',
                        'PA21': 'N_HC2_0_PA12',
                        'PA22': 'N_C2_PA22'
                        },
                    12: {
                        'PA11': 'H_N_0_PA11',
                        'PA12': 'DEL',
                        'PA21': 'H_N_0_PA12',
                        'PA22': 'DEL'
                        },
                    13: {
                        'PA11': 'H_N_0_PA11',
                        'PA12': 'H_N_0_PA12',
                        'PA21': 'DEL',
                        'PA22': 'DEL'
                        }
                    }[count][target]
            if target == 'PA21':
                assert ffatypes[i].split('_')[-1] in ['PA12', 'DEL']
            else:
                assert ffatypes[i].split('_')[-1] in [target, 'DEL']

        
    # Remove imine hydrogens assigned ffatype DEL
    system_indices = [i for i in range(sys.natom) if not ffatypes[i] == 'DEL']
    new_system = System(numbers = numbers, pos = pos, bonds = bonds, ffatypes = ffatypes, rvecs = rvecs)
    return new_system.subsystem(system_indices)
    
for name, system in [('ABCDEF', abcdef)]:
    for label, indices in iter_todo():
        print('{}_{} ({}): {}'.format(name, label, len(indices), indices))
        new_system = convert_system(system, indices)
        folder = 'output'
        if not os.path.exists(folder):
            os.makedirs(folder)
        new_system.to_file('output/{}_{}.chk'.format(name, label))

        

