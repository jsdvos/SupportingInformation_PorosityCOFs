import os
import h5py

from yaff.system import System, log
log.set_level(0)
from molmod.periodic import periodic as pt
from molmod.units import angstrom

def write_frame(system, counter, path):
    folder = os.path.join(path, 'frame{}'.format(counter))
    if not os.path.exists(folder):
        os.mkdir(folder)
    system.to_file(os.path.join(folder, 'frame{}.chk'.format(counter)))

for temp in ['77K', '293K', '400K', '500K']:
    for struct in os.listdir('../../MDrun/{}'.format(temp)):
        print(struct)
        fn_h5 = '../../MDrun/{}/{}/traj.h5'.format(temp, struct)
        if not os.path.exists(fn_h5): continue
        system = System.from_file(fn_h5)
        write_frame(system, 0, '../calculations/0K/' + struct)
        elements = [pt[i].symbol for i in system.numbers]
        f = h5py.File(fn_h5, 'r')
        assert len(f['trajectory']['counter']) == 15001
        for i in range(7500, 15001, 15):
            pos = f['trajectory']['pos'][i]
            rvecs = f['trajectory']['cell'][i]
            system.pos[:] = pos
            system.cell.update_rvecs(rvecs)
            write_frame(system, i, '../calculations/{}/{}'.format(temp, struct))
