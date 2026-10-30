"""
This Eophis script manages coupling between toy_earth.py and models.py. OASIS interface is already hard-coded in toy_earth.py.
Purposes are to write coupling namelist during preproduction mode, configure OASIS interface for models.py, and orchestrate connexions between exchanged data and models.py.
Checkout docstrings in toy_earth.py and documentation to learn more about what coupling is supposed to achieve.

"""
# eophis API
import eophis
from eophis import Freqs, Domains
# other modules
import argparse
import os
import torch

def earth_info():

    # namelist info
    with open("run.def") as f:
        d = dict(line.strip().split("=", 1) for line in f if "=" in line and not line.lstrip().startswith("#"))
        nday = int(d["nday"])

    step = 21600
    ite_0 = 1
    ite_end = int(24*3600*nday / step)

    # tunnel configuration
    tunnel_config = list()
    tunnel_config.append( { 'label' : 'TO_EARTH', \
                            'grids' : { 'lmdz' :  {'npts' : (96,96)} }, \
                            'exchs' : [ {'freq' : step, 'grd' : 'lmdz', 'lvl' : 79, 'in' : ['u','v','tpot','pres','cos_day','sin_day'], 'out' : ['du','dv'] } , \
                                        {'freq' : step, 'grd' : 'lmdz', 'lvl' : 1, 'in' : ['psol','topo'], 'out' : [] } ] }
                        )

    return tunnel_config, step, ite_end, ite_0


def preproduction():
    """ Preproduction features: register Tunnel, write the coupling namelist from Tunnel registration and parameters in 'earth_namelist'. """
    eophis.info('========= EOPHIS: Pre-Production =========')
    eophis.info('  Aim: write coupling namelist\n')

    # earth info
    tunnel_config, step, it_end, it_0 = earth_info()
    total_time = (it_end - it_0 + 1) * step

    # tunnel registration (lazy) compulsory to update namelist
    eophis.register_tunnels( tunnel_config )
    
    # write updated namelist
    eophis.write_coupling_namelist( simulation_time=total_time )


def production():
    """
    Production mode features: register and create Tunnel object (will compare coherence between Tunnel and namcouple content). Tunnel opening starts coupling.
    Import add_100() model from models.py
    Perform a static reception (msk)
    Assemble a Loop and a Router
    
    """
    eophis.info('========= EOPHIS: Production =========')
    eophis.info('  Aim: execute coupled simulation\n')

    #  Earth Coupling
    # ++++++++++++++++
    tunnel_config, step, it_end, it_0 = earth_info()
    niter = it_end - it_0 + 1
    total_time = niter * step

    # tunnel registration (lazy)
    earth, = eophis.register_tunnels( tunnel_config )

    # link all tunnels (beware, dormant errors will likely appear here)
    eophis.open_tunnels()

    #  Models
    # ++++++++
    from models import inference_pytorch_lightning

    #  Assemble
    # ++++++++++
    @eophis.all_in_all_out(geo_model=earth, step=step, niter=niter)
    def loop_core(**inputs):
        """
        Loop is defined with the decorator and time step informations from 'earth_namelist'.
        Content of loop_core is the Router. Connexions between exchanged variables and models are defined here.
        inputs dictionnary contains the variables that Eophis automatically received from toy_earth.py through OASIS.
        sst_var is the result of add_100() with sst passed as argument, same with svt and svt_var.
        Fill 'outputs' dictionnary  will let Eophis knows what to send back to toy_earth.py through OASIS under corresponding key
        
        """
        outputs = {}
        stacked_res = inference_pytorch_lightning( inputs['u'], inputs['v'], inputs['tpot'], inputs['pres'],
                inputs['psol'], inputs['topo'], inputs['cos_day'], inputs['sin_day'] )

        outputs['du'], outputs['dv']  = stacked_res
        return outputs

    #  Run
    # +++++
    eophis.starter(loop_core)
    
# ============================ #
if __name__=='__main__':
    """ Switch between preproduction and production mode. Call different instructions depending on the selected mode. """

    parser = argparse.ArgumentParser()
    parser.add_argument('--exec', dest='exec', type=str, default='prod', help='Execution type: preprod or prod')
    args = parser.parse_args()

    eophis.set_mode(args.exec)
    
    if args.exec == 'preprod':
        preproduction()
    elif args.exec == 'prod':
        production()
    else:
        eophis.abort(f'Unknown execution mode {args.exec}, use "preprod" or "prod"')
