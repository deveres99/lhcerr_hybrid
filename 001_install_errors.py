import argparse
import yaml
import numpy as np

import xtrack as xt
import xpart as xp
from cpymad.madx import Madx


from aperture_tools import patch_aperture_with


parser = argparse.ArgumentParser()
parser.add_argument("--config_path", 
                    nargs='?', 
                    help="path to yaml config file", 
                    type=str, 
                    default='./example_injection.yaml')


def install_errors(beam, optics_scenario, errors_scenario, path_errors, errors, seed, energy, save_path):
    if errors_scenario not in ["injection", "collision"]:
        raise ValueError(f"Undefined error scenario: {errors_scenario}!")
    if optics_scenario == "injection" and errors_scenario != "injection":
        raise ValueError(f"Mismatched optics and error scenario! ({optics_scenario} vs {errors_scenario})")
    
    mad = Madx()

    save_name = f"lhcb{beam}_errors{seed}_{optics_scenario}"
    load_name = f"lhcb{beam}_clean_{optics_scenario}"

    mad.call(save_path + load_name + ".seq")

    mad.input(f"Beam,particle=proton,sequence=lhcb{beam},energy={energy*1e-9};")

    mybeam = 1 if beam == 1 else 4 
    mad.input(f'''
System,"ln -fns /afs/cern.ch/eng/lhc/optics/runIII run3opt";
call,file="run3opt/toolkit/macro.madx";
              
mybeam={mybeam};

use, sequence=lhcb{beam};

on_disp = 0;

myseed               = {seed};
ver_lhc_run          = 3;
on_errors_LHC        = 1;
call, file="{path_errors}LHC/Msubroutines.madx";
call, file="{path_errors}LHC/Msubroutines_MS_MSS_MO.madx";
call, file="{path_errors}/Orbit_Routines.madx";
call, file="{path_errors}/SelectLHCMonCor_v1.madx";
call, file="{path_errors}HL-LHC/macro_error_v1.madx";   ! macros for error generation in the new IT/D1's
exec, ON_ALL_MULT;
ON_A1s =  0 ; ON_A1r =  0 ; ON_B1s =  0 ; ON_B1r =  0 ;
ON_A2s =  0 ; ON_A2r =  0 ; ON_B2s =  0 ; ON_B2r =  0 ;
ON_A3s =  1 ; ON_A3r =  1 ; ON_B3s =  1 ; ON_B3r =  1 ;
ON_A4s =  1 ; ON_A4r =  1 ; ON_B4s =  1 ; ON_B4r =  1 ;
ON_A5s =  1 ; ON_A5r =  1 ; ON_B5s =  1 ; ON_B5r =  1 ;
ON_A6s =  1 ; ON_A6r =  1 ; ON_B6s =  1 ; ON_B6r =  1 ;
ON_A7s =  1 ; ON_A7r =  1 ; ON_B7s =  1 ; ON_B7r =  1 ;
ON_A8s =  1 ; ON_A8r =  1 ; ON_B8s =  1 ; ON_B8r =  1 ;
ON_A9s =  1 ; ON_A9r =  1 ; ON_B9s =  1 ; ON_B9r =  1 ;
ON_A10s = 1 ; ON_A10r = 1 ; ON_B10s = 1 ; ON_B10r = 1 ;
ON_A11s = 1 ; ON_A11r = 1 ; ON_B11s = 1 ; ON_B11r = 1 ;
ON_A12s = 1 ; ON_A12r = 1 ; ON_B12s = 1 ; ON_B12r = 1 ;
ON_A13s = 1 ; ON_A13r = 1 ; ON_B13s = 1 ; ON_B13r = 1 ;
ON_A14s = 1 ; ON_A14r = 1 ; ON_B14s = 1 ; ON_B14r = 1 ;
ON_A15s = 1 ; ON_A15r = 1 ; ON_B15s = 1 ; ON_B15r = 1 ;

twiss, table=nominal;   // used by orbit correction
beta.ip1=table(twiss,IP1,betx);value,beta.ip1;
''')

    all_wise_types = [
        "MB", "MBRB", "MBRC", "MBRS", "MBX", "MBXW", "MBW", 
        "MQW", "MQTL", "MQMC", "MQX", "MQY", "MQM", "MQML", "MQ"
    ]
    all_fidel_types = ["MS", "MSS", "MO"]
    requested_bend_types = [err for err in errors if err.startswith("MB")]
    requested_quad_types = [err for err in errors if err.startswith("MQ")]
    requested_fidel_types = [err for err in errors if err.startswith("MS") or err.startswith("MO")]
    if np.any(np.isin(errors, all_wise_types)):
        mad.input(f'''
! disable crossing bumps
exec, crossing_save;
exec, crossing_disable;

readtable, file="{path_errors}LHC/rotations_Q2_integral.tab";
readtable, file="{path_errors}LHC/wise/{errors_scenario}_errors-emfqcs-{seed}.tfs" ;
''')
        for err in requested_bend_types:
            mad.input(f'''
call, file="{path_errors}LHC/Efcomp_{err}.madx"  ;
''')
        mad.input("ON_B2Saux=on_B2S;on_B2S=0;")

        for err in requested_quad_types:
            mad.input(f'''
call, file="{path_errors}LHC/Efcomp_{err}.madx"  ;
''')
        mad.input("on_B2S=ON_B2Saux;")

    if np.any(np.isin(errors, all_fidel_types)):
        mad.input(f'''
readtable, file="{path_errors}LHC/fidel/injection_errors-emfqcs-{seed}.tfs" ;
''')
        for err in requested_fidel_types:
            mad.input(f'''
call, file="{path_errors}LHC/Efcomp_{err}.madx"  ;
''')
    
    mad.input(f"exec, crossing_restore;")
    mad.input(f"save, sequence=lhcb{beam}, file={save_path + save_name + '.seq'};")

    if( beam == 1):
        mad_sequence = mad.sequence.lhcb1
    else:
        mad_sequence = mad.sequence.lhcb2

    # With apertures
    line = xt.Line.from_madx_sequence(mad_sequence, apply_madx_errors=True, enable_field_errors=True, install_apertures=True, deferred_expressions=True)
    line.particle_ref = xp.Particles(p0c=energy, q0=1, mass0=xp.PROTON_MASS_EV)

    collimators = [name for name in line.element_names
                    if (name.startswith('tc') or name.startswith('td'))
                    and not '_aper' in name and not name[-4:-2]=='mk' and not name[:4] == 'tcds'
                    and not name[:4] == 'tcdd' and not name[:5] == 'tclim' and not name[:3] == 'tca'
                    and not (name[-5]=='.' and name[-3]=='.') and not name[:5] == 'tcdqm'
                ]
    df = line.check_aperture(needs_aperture=collimators)
    missing_apertures = df.loc[df.has_aperture_problem, 'name'].values
    if beam == 1:
        patch_aperture_with(['mo.28r3.b1', 'mo.32r3.b1'], 'mo.22r1.b1_mken_aper', missing_apertures, line)
        patch_aperture_with(['mqwa.f5l7.b1..1', 'mqwa.f5l7.b1..2', 'mqwa.f5l7.b1..3',
                            'mqwa.f5l7.b1..4', 'mqwa.f5r7.b1..1', 'mqwa.f5r7.b1..2',
                            'mqwa.f5r7.b1..3', 'mqwa.f5r7.b1..4'
                            ], 'mqwa.e5l3.b1_mken_aper', missing_apertures, line)
        patch_aperture_with(['tdisa.a4l2.b1', 'tdisb.a4l2.b1', 'tdisc.a4l2.b1'
                            ], xt.LimitRect(min_x=-0.043, max_x=0.043, min_y=-0.055, max_y=0.055), missing_apertures, line)
        patch_aperture_with('tcld.a11r2.b1', xt.LimitEllipse(a=4e-2, b=4e-2), missing_apertures, line)
        patch_aperture_with(['tcspm.b4l7.b1', 'tcspm.e5r7.b1', 'tcspm.6r7.b1'
                            ], xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-06), 
                            missing_apertures, line)
        patch_aperture_with(['tcpch.a4l7.b1', 'tcpcv.a6l7.b1'
                            ], xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-06), 
                            missing_apertures, line)
    else:
        patch_aperture_with(['mo.32r3.b2', 'mo.28r3.b2'], 'mo.22l1.b2_mken_aper', missing_apertures, line)
        patch_aperture_with(['mqwa.f5r7.b2..1', 'mqwa.f5r7.b2..2', 'mqwa.f5r7.b2..3',
                            'mqwa.f5r7.b2..4', 'mqwa.f5l7.b2..1', 'mqwa.f5l7.b2..2',
                            'mqwa.f5l7.b2..3', 'mqwa.f5l7.b2..4'
                            ], 'mqwa.e5r3.b2_mken_aper', missing_apertures, line)
        patch_aperture_with(['tdisa.a4r8.b2', 'tdisb.a4r8.b2', 'tdisc.a4r8.b2'
                            ], xt.LimitRect(min_x=-0.043, max_x=0.043, min_y=-0.055, max_y=0.055), missing_apertures, line)
        patch_aperture_with('tcld.a11l2.b2', xt.LimitEllipse(a=4e-2, b=4e-2), missing_apertures, line)
        patch_aperture_with(['tcspm.d4r7.b2', 'tcspm.b4r7.b2', 'tcspm.e5l7.b2', 'tcspm.6l7.b2'
                            ], xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-06), 
                            missing_apertures, line)
        patch_aperture_with(['tcpch.a5r7.b2', 'tcpcv.a6r7.b2'
                            ], xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-06), 
                            missing_apertures, line)

    line.to_json(save_path + save_name + ".json")
    
    # Without apertures
    line = xt.Line.from_madx_sequence(mad_sequence, apply_madx_errors=True, enable_field_errors=True, install_apertures=False, deferred_expressions=True)
    line.particle_ref = xp.Particles(p0c=energy, q0=1, mass0=xp.PROTON_MASS_EV)

    line.to_json(save_path + save_name + "_noaper.json")


def main():
    args = parser.parse_args()
    config_path = args.config_path
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    if config["beam"] == 0:
        install_errors(
            beam=1, 
            optics_scenario=config["optics_scenario"], 
            errors_scenario=config["errors_scenario"], 
            path_errors=config["path_errors"], 
            errors=config["errors"], 
            seed=config["error_seed"], 
            energy=float(config["energy"]), 
            save_path=config["save_path"]
        )
        install_errors(
            beam=2, 
            optics_scenario=config["optics_scenario"], 
            errors_scenario=config["errors_scenario"], 
            path_errors=config["path_errors"], 
            errors=config["errors"], 
            seed=config["error_seed"], 
            energy=float(config["energy"]), 
            save_path=config["save_path"]
        )
    elif config["beam"] == 1:
        install_errors(
            beam=1, 
            optics_scenario=config["optics_scenario"], 
            errors_scenario=config["errors_scenario"], 
            path_errors=config["path_errors"], 
            errors=config["errors"], 
            seed=config["error_seed"], 
            energy=float(config["energy"]), 
            save_path=config["save_path"]
        )
    elif config["beam"] == 2:
        install_errors(
            beam=2, 
            optics_scenario=config["optics_scenario"], 
            errors_scenario=config["errors_scenario"], 
            path_errors=config["path_errors"], 
            errors=config["errors"], 
            seed=config["error_seed"], 
            energy=float(config["energy"]), 
            save_path=config["save_path"]
        )
    else:
        raise ValueError(f"Beam must be 0, 1 or 2, {config['beam']} is not accepted!")


if __name__ == "__main__":
    main()
