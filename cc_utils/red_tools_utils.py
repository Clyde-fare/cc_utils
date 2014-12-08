__author__ = 'clyde'

import gaussian_job_manager
from gaussian_job_manager import Job
from pbs_util import pbs
import os

def set_red(**kwargs):
    """sets the parameters for red tools

    xred boolean determines if the xred gui is used default false
    nodes int determines number of processors used in qm calculation
    qm str determines which package is used to perform the calculation
    opt determines whether geometry optimization will be performed
    mep determines whether charge fitting will be carried out
    refit determines whether we are rebuilding from a previous R.E.D. job
    type determines the charge derivation model
    cor determines the numerical accuracy of the charge values"""

    xred = kwargs.get('xred', False)
    xred = 'On' if xred else 'OFF'

    nodes = kwargs.get('nodes', 8)
    nodes = str(nodes)

    qm = kwargs.get('qm', 'GAUSSIAN')
    data_dir = kwargs.get('data_dir', 'Data-RED')

    opt = kwargs.get('opt', True)
    opt = 'On' if opt else 'OFF'

    mep = kwargs.get('mep', True)
    mep = 'On' if mep else 'False'

    refit = kwargs.get('refit', False)
    refit = 'On' if refit else 'Off'

    calc_type = kwargs.get('calc_type', 'RESP-A1')

    cor = kwargs.get('cor', 4)
    cor = str(cor)

    red_dict = {'XRED': xred, 'NP': nodes, 'QMSOFT': qm, 'DIR': data_dir, 'OPT_Calc': opt, 'MEPCHR_Calc': mep, 'Re_Fit': refit, 'CHR_TYP': calc_type, 'COR_CHR': cor}

    red_home = os.environ['REDHOME']
    red_std_fn = 'RED-vIII.5.pl'
    red_tools_fn = red_home + '/' + red_std_fn
    with open(red_tools_fn) as red_f:
        red_contents = red_f.readlines()

    param_start_index = next(i for i, l in enumerate(red_contents) if 'MAIN PROGRAM' in l)
    final_red_contents = red_contents[param_start_index:]

    for param in red_dict.keys():
        param_line_index = next(i for i, l in enumerate(final_red_contents) if '$' + param + '=' in l.replace(' ', ''))
        param_line = red_contents[param_start_index + param_line_index]
        new_param_line = param_line.split('=')[0] + ' =  "{p}";\n'.format(p=red_dict[param])
        red_contents[param_start_index + param_line_index] = new_param_line

    mod_red_tools_fn = 'mod_' + red_std_fn
    with open(mod_red_tools_fn, 'w') as mod_red_f:
        mod_red_f.writelines(red_contents)

    os.system('chmod 755 {mod_red}'.format(mod_red=mod_red_tools_fn))
    return mod_red_tools_fn


def set_red_job(red_fn, job):
    job_script = job.gen_header() + """module load ambertools
module load gaussian
export GAUSS_SCRDIR=$TMPDIR
current_path=`pwd`
echo $GAUSS_SCRDIR
cd $(echo $PBS_O_WORKDIR)
$current_path/{rd_fn} > $current_path/RED_out.log
""".format(rd_fn=red_fn)

    with open('red_job.sh', 'w') as red_job_f:
        red_job_f.write(job_script)

    return 'red_job.sh'


def red_on_server(p2n_f_nm, **kwargs):
    nodes = kwargs.get('nodes', 8)
    memory = kwargs.get('memory', 8*1400)
    time = kwargs.get('time', 5)
    queue = kwargs.get('queue')

    job = Job(procs=nodes, memory=memory, walltime=time, queue=queue)

    xred = kwargs.get('xred', False)
    opt = kwargs.get('opt', True)
    mep = kwargs.get('mep', True)
    refit = kwargs.get('refit', False)
    qm = kwargs.get('qm', 'GAUSSIAN')
    calc_type = kwargs.get('calc_type', 'RESP-A1')
    cor = kwargs.get('cor', 4)
    data_dir = kwargs.get('data_dir', 'Data-RED')

    red_fn = set_red(xred=xred, nodes=nodes, opt=opt, mep=mep, refit=refit, qm=qm, data_dir=data_dir, calc_type=calc_type, cor=cor)
    red_job_fn = set_red_job(red_fn, job)

    os.system('mv {fn}  Mol_red1.p2n'.format(fn=p2n_f_nm))
    gaussian_job_manager.send_to_cx1('Mol_red1.p2n')

    gaussian_job_manager.send_to_cx1(red_job_fn)
    gaussian_job_manager.send_to_cx1(red_fn)

    qsub_command = job.exec_command(red_job_fn)

    ssh = pbs.connect_server(ssh=True)
    stdin, stdout, stderr = ssh.exec_command(qsub_command)
    qid = stdout.read().split('.')[0]
    ssh.close()
    return qid


def run_red(p2n_f_nm, **kwargs):
    nodes = kwargs.get('nodes', 8)
    xred = kwargs.get('xred', False)
    opt = kwargs.get('opt', True)
    mep = kwargs.get('mep', True)
    refit = kwargs.get('refit', False)
    qm = kwargs.get('qm', 'GAUSSIAN')
    calc_type = kwargs.get('calc_type', 'RESP-A1')
    cor = kwargs.get('cor', 4)

    red_fn = set_red(xred=xred, nodes=nodes, opt=opt, mep=mep, refit=refit, qm=qm, calc_type=calc_type, cor=cor)
    os.system('cp {p2n} Mol_red1.p2n'.format(p2n=p2n_f_nm))
    os.system('./{red} &'.format(red=red_fn))


from pbs_util.pbs import PBSUtilQStatError
def get_red_data_from_server(id):
    local_home = os.environ['ASE_HOME']
    try:
        active_dir = os.getcwd().split(local_home)[1]
    except IndexError:
        raise RuntimeError('Not running from within ASE_HOME')

    status=True
    try:
        status = pbs.qstat(id)
    except PBSUtilQStatError:
        os.system('rsync -av {l} {s}'.format(l=os.getcwd(), s=os.environ['GAUSS_HOST'] + active_dir))

    return status