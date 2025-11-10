#!/usr/bin/env python3

import datetime
from .tiosutils import TiosXTCWriter, TiosNCWriter
from ._version import __version__
from tqdm import tqdm
from argparse import ArgumentParser
from pathlib import Path
import sys
import signal
from datetime import datetime, timedelta


#. sys.tracebacklimit = 0  # Suppress traceback unless in debug mode


class GracefulKiller:
    kill_now = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, signum, frame):
        self.kill_now = True


def tios_write_cli():
    parser = ArgumentParser(description="Write simulation data to"
                            " an XTC or NetCDF file using MQTT.")
    parser.add_argument("mqtt_broker", type=str,
                        help="The address of the MQTT broker.")
    parser.add_argument("sim_id", type=str,
                        help="The unique identifier for the simulation.")
    parser.add_argument("output_file", type=str,
                        help="The output XTC/NC file path.")
    parser.add_argument("--timeout", type=int, default=60,
                        help="Timeout in seconds for waiting for new frames.")
    parser.add_argument("--port", type=int, default=1883,
                        help="The port number of the MQTT broker.")
    parser.add_argument("--max_frames", type=int, default=None,
                        help="Maximum number of frames to write.")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    mqtt_broker = args.mqtt_broker
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

    killer = GracefulKiller()

    with writer(mqtt_broker, sim_id, output_file,
                timeout=timeout, port=port) as f:
        t = tqdm(unit=" frames", total=max_frames)
        while not (f.timedout() or killer.kill_now):
            f.write_frame()
            t.update()
        t.close()


def tios_ls_cli():
    parser = ArgumentParser(description="List available simulations"
                            " from an MQTT broker.")
    parser.add_argument("mqtt_broker", type=str,
                        help="The address of the MQTT broker.")
    parser.add_argument("--port", type=int, default=1883,
                        help="The port number of the MQTT broker.")
    parser.add_argument("--timeout", type=int, default=5,
                        help="Timeout in seconds for waiting for simulations.")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    mqtt_broker = args.mqtt_broker
    port = args.port

    from .tiosutils import get_simulations
    now = datetime.now()
    simulations = get_simulations(mqtt_broker, port=port,
                                   timeout=args.timeout,
                                   patient=False)
    if len(simulations ) > 0:
        print("Available simulations:")
        running = {}
        stopped = {}
        for s in simulations:
            if simulations[s]['is_running']:
                running[s] = simulations[s]
            else:
                stopped[s] = simulations[s]
        if len(running) > 0:
            print("Running simulations:")
            running = {k: v for k, v in sorted(running.items(), key=lambda item: item[1]['last_update'], reverse=True)}
            for simId, simData in running.items():
                summary = simData['summary']
                last_update = datetime.fromtimestamp(simData['last_update'])
                print(f" - {simId}: {summary} (last update: {last_update})")
        if len(stopped) > 0:
            print("Stopped simulations:")
            stopped = {k: v for k, v in sorted(stopped.items(), key=lambda item: item[1]['last_update'], reverse=True)}
            for simId, simData in stopped.items():
                summary = simData['summary']
                last_update = datetime.fromtimestamp(simData['last_update'])
                interval = (now - last_update)
                if interval.days > 0:
                    interval = f"{interval.days} days"
                elif interval.seconds >= 3600:
                    hours = interval.seconds // 3600
                    interval = f"{hours} hours"
                elif interval.seconds >= 60:
                    minutes = interval.seconds // 60
                    interval = f"{minutes} minutes"
                else:
                    interval = f"{interval.seconds} seconds"
                print(f" - {simId}: {summary} (last update: {interval} ago)")
    else:
        print("No simulations found.")
