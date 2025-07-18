import argparse
import yaml
import numpy as np

import xtrack as xt
import xpart as xp

from correction_tools import *


parser = argparse.ArgumentParser()
parser.add_argument("--config_path", 
                    nargs='?', 
                    help="path to yaml config file", 
                    type=str, 
                    default='./example_injection.yaml')


def correct_and_tune(beam, optics_scenario, errors_scenario, seed, save_path):
    line_ref = xt.Line.from_json(save_path + f"lhcb{beam}_clean_{optics_scenario}.json")
    
    line = xt.Line.from_json(save_path + f"lhcb{beam}_errors{seed}_{optics_scenario}.json")
    
    line_noaper = xt.Line.from_json(save_path + f"lhcb{beam}_errors{seed}_{optics_scenario}_noaper.json")

    line_ref.twiss_default["method"] = '4d'
    
    line.twiss_default["method"] = '4d'
    
    line_noaper.twiss_default["method"] = '4d'

    assign_spool_pieces(line, beam)
    assign_spool_pieces(line_noaper, beam)

    # save_crossing(line_ref)
    # disable_crossing(line_ref)
    # tw_ref = line_ref.twiss()
    # enable_crossing(line_ref)
    
    # set_correctors(line)
    # save_crossing(line)
    # disable_crossing(line)
    # line.correct_trajectory(twiss_table=tw_ref)
    # enable_crossing(line)
    
    # set_correctors(line_noaper)
    # save_crossing(line_noaper)
    # disable_crossing(line_noaper)
    # line_noaper.correct_trajectory(twiss_table=tw_ref)
    # enable_crossing(line_noaper)

    tw_ref = line_ref.twiss()

    match_tune_chrom(line, tw_ref.qx, tw_ref.qy, tw_ref.dqx, tw_ref.dqy, tol=[1e-4, 5e-5, 1e-5, 5e-6])
    match_coupling(line, beam, c_minus=0.001, tol=1e-5)
    match_tune_chrom(line, tw_ref.qx, tw_ref.qy, tw_ref.dqx, tw_ref.dqy, tol=[1e-4, 5e-5, 1e-5, 5e-6])

    match_tune_chrom(line_noaper, tw_ref.qx, tw_ref.qy, tw_ref.dqx, tw_ref.dqy, tol=[1e-4, 5e-5, 1e-5, 5e-6])
    match_coupling(line_noaper, beam, c_minus=0.001, tol=1e-5)
    match_tune_chrom(line_noaper, tw_ref.qx, tw_ref.qy, tw_ref.dqx, tw_ref.dqy, tol=[1e-4, 5e-5, 1e-5, 5e-6])

    line.to_json(save_path + f"lhcb{beam}_errors{seed}corrected_{optics_scenario}.json")

    line_noaper.to_json(save_path + f"lhcb{beam}_errors{seed}corrected_{optics_scenario}_noaper.json")


def main():
    args = parser.parse_args()
    config_path = args.config_path
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    if config["beam"] == 0:
        correct_and_tune(
            beam=1, 
            optics_scenario=config["optics_scenario"], 
            errors_scenario=config["errors_scenario"], 
            seed=config["error_seed"], 
            save_path=config["save_path"]
        )
        correct_and_tune(
            beam=2, 
            optics_scenario=config["optics_scenario"], 
            errors_scenario=config["errors_scenario"], 
            seed=config["error_seed"], 
            save_path=config["save_path"]
        )
    elif config["beam"] == 1:
        correct_and_tune(
            beam=1, 
            optics_scenario=config["optics_scenario"], 
            errors_scenario=config["errors_scenario"], 
            seed=config["error_seed"], 
            save_path=config["save_path"]
        )
    elif config["beam"] == 2:
        correct_and_tune(
            beam=2, 
            optics_scenario=config["optics_scenario"], 
            errors_scenario=config["errors_scenario"], 
            seed=config["error_seed"], 
            save_path=config["save_path"]
        )
    else:
        raise ValueError(f"Beam must be 0, 1 or 2, {config['beam']} is not accepted!")


if __name__ == "__main__":
    main()
