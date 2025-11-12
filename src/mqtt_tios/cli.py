#!/usr/bin/env python3


from .tiosutils import TiosXTCWriter, TiosNCWriter
from ._version import __version__
from tqdm import tqdm
from argparse import ArgumentParser
from pathlib import Path
import signal
from datetime import datetime


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
                                   timeout=args.timeout)
    if len(simulations ) > 0:
        print("Available simulations:")
        sorted_sims = {k: v for k, v in sorted(simulations.items(), key=lambda item: item[1]['last_update'], reverse=True)}
        
        
        for simId, simData in sorted_sims.items():
            summary = simData['summary']
            last_update = datetime.fromtimestamp(simData['last_update'])
            delta = (now - last_update)
            if delta.days > 0:
                interval = f"{delta.days} day"
                if delta.days > 1:
                    interval += "s"
            elif delta.seconds >= 3600:
                hours = delta.seconds // 3600
                interval = f"{hours} hour"
                if hours > 1:
                    interval += "s"
            elif delta.seconds >= 60:
                minutes = delta.seconds // 60
                interval = f"{minutes} minute"
                if minutes > 1:
                    interval += "s"
            else:
                interval = f"{delta.seconds} second"
                if delta.seconds > 1:
                    interval += "s"
            status = []
            if simData['is_running']:
                status.append("running")
            else:
                status.append("stopped")
            if simData['has_checkpoint']:
                status.append("checkpointed")
            status_str = ", ".join(status)
            print(f" - {simId}: {summary}: last updated {interval} ago, {status_str}")
    else:
        print("No simulations found.")
