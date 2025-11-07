#!/usr/bin/env python3

from tiosutils import TiosXTCWriter, TiosNCWriter
from tqdm import tqdm
from argparse import ArgumentParser
from pathlib import Path
import sys
import signal

sys.tracebacklimit = 0
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


if __name__ == "__main__":
    tios_write_cli()
