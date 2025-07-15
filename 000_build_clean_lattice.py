import argparse
import yaml

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


def build_clean_lattice(beam, path_machine_config, scenario, energy, save_path):
    mad = Madx()

    mad.call(path_machine_config + scenario + f"/track_{scenario}_b{beam}.madx")
    if( beam == 1):
        mad_sequence = mad.sequence.lhcb1
        save_name = f"lhcb1_clean_{scenario}"
    else:
        mad_sequence = mad.sequence.lhcb2
        save_name = f"lhcb2_clean_{scenario}"
    mad.input(f"save, sequence=lhcb{beam}, file={save_path + save_name + '.seq'};")

    # With apertures
    line = xt.Line.from_madx_sequence(mad_sequence, apply_madx_errors=False, install_apertures=True, deferred_expressions=True)
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
    line = xt.Line.from_madx_sequence(mad_sequence, apply_madx_errors=False, install_apertures=False, deferred_expressions=True)
    line.particle_ref = xp.Particles(p0c=energy, q0=1, mass0=xp.PROTON_MASS_EV)

    line.to_json(save_path + save_name + "noaper.json")

def main():
    args = parser.parse_args()
    config_path = args.config_path
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    if config["beam"] == 0:
        build_clean_lattice(
            beam=1, 
            path_machine_config=config["path_machine_config"], 
            scenario=config["optics_scenario"], 
            energy=config["energy"], 
            save_path=config["save_path"]
        )
        build_clean_lattice(
            beam=2, 
            path_machine_config=config["path_machine_config"], 
            scenario=config["optics_scenario"], 
            energy=config["energy"], 
            save_path=config["save_path"]
        )
    elif config["beam"] == 1:
        build_clean_lattice(
            beam=1, 
            path_machine_config=config["path_machine_config"], 
            scenario=config["optics_scenario"], 
            energy=config["energy"], 
            save_path=config["save_path"]
        )
    elif config["beam"] == 2:
        build_clean_lattice(
            beam=2, 
            path_machine_config=config["path_machine_config"], 
            scenario=config["optics_scenario"], 
            energy=config["energy"], 
            save_path=config["save_path"]
        )
    else:
        raise ValueError(f"Beam must be 0, 1 or 2, {config['beam']} is not accepted!")


if __name__ == "__main__":
    main()