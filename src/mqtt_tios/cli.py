""" cli.py: command line tools: tios_ls and tios_write """

from .tiosutils import TiosXTCWriter, TiosNCWriter
from ._version import __version__
from .config import config
from tqdm import tqdm
from argparse import ArgumentParser
from pathlib import Path
from .control import killer


# sys.tracebacklimit = 0  # Suppress traceback unless in debug mode

def tios_collect_cli():
    parser = ArgumentParser(description="Save simulation data to"
                            " an XTC or NetCDF file using MQTT.")
    parser.add_argument("--broker", type=str,
                        default=config.broker,
                        help="The address of the MQTT broker.")
    parser.add_argument("sim_id", type=str,
                        help="The unique identifier for the simulation.")
    parser.add_argument("output_file", type=str,
                        help="The output XTC/NC file path.")
    parser.add_argument("--timeout", type=int, default=60,
                        help="Timeout in seconds for waiting for new frames.")
    parser.add_argument("--port", type=int, default=config.port,
                        help="The port number of the MQTT broker.")
    parser.add_argument("--max_frames", type=int, default=None,
                        help="Maximum number of frames to write.")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    mqtt_broker = args.broker
    sim_id = args.sim_id
    output_file = args.output_file
    timeout = args.timeout
    port = args.port
    max_frames = args.max_frames

    ext = Path(output_file).suffix.lower()
    if ext == '.xtc':
        writer = TiosXTCWriter
    elif ext == '.nc':
        writer = TiosNCWriter
    else:
        raise ValueError(f"Unsupported file extension: {ext}. Use .xtc or .nc")

    with writer(sim_id, output_file,
                mqtt_broker=mqtt_broker,
                timeout=timeout, port=port) as f:
        t = tqdm(unit=" frames", total=max_frames)
        n_frames = 0
        while not (f.timedout or killer.kill_now):
            f.write_frame()
            n_frames += 1
            if max_frames and n_frames >= max_frames:
                killer.kill_now = True
            t.update()
        t.close()


def tios_ls_cli():
    parser = ArgumentParser(description="List available simulations"
                            " from an MQTT broker.")
    parser.add_argument("--broker", type=str,
                        default=config.broker,
                        help="The address of the MQTT broker.")
    parser.add_argument("--port", type=int, default=config.port,
                        help="The port number of the MQTT broker.")
    parser.add_argument("--timeout", type=int, default=5,
                        help="Timeout in seconds for waiting for simulations.")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    mqtt_broker = args.broker
    port = args.port

    from .tiosutils import get_simulations

    simulations = get_simulations(mqtt_broker, port=port,
                                  timeout=args.timeout)
    if len(simulations) > 0:
        print("Available simulations:")
        sorted_sims = {k: v for k, v in sorted(
            simulations.items(),
            key=lambda item: item[1]['status'],
            reverse=True)}

        for simId, simData in sorted_sims.items():
            print(f" - {simId}: {simData['summary']}: {simData['status']}")
    else:
        print("No simulations found.")
