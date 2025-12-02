""" tiosutils.py: lower-level routines for command-line tools """
from time import sleep
import zlib
import numpy as np
from mdtraj.utils import box_vectors_to_lengths_and_angles
from mdtraj.formats import NetCDFTrajectoryFile, XTCTrajectoryFile
from .clients import TiosSubscriber, TiosMonitor
from .config import config
from .control import killer

EOT = 'EOT'.encode('utf-8')


def interruptable_get(client, timeout=None):
    """ A get that can be interrupted"""
    time_left = timeout or 1
    while time_left > 0 and client.states.qsize() == 0 and not killer.kill_now:
        sleep(1)
        client.poke()
        if timeout:
            time_left -= 1
    if client.states.empty():
        return None
    else:
        return client.states.get()


class TiosXTCWriter():
    def __init__(self, sim_id, xtcfilename, mqtt_broker=None, port=None,
                 timeout=60):
        # Set the broker address and port
        self.broker_address = mqtt_broker or config.broker
        self.port = port or config.port
        self.xtcfilename = xtcfilename
        self.timeout = timeout

        self._client = TiosSubscriber(sim_id,
                                      broker_address=self.broker_address,
                                      port=port)
        self.timedout = False
        self.xtcfile = None
        self.framebuffer = None
        self.saved_frames = 0

    def write_frame(self):

        zdata = interruptable_get(self._client, timeout=self.timeout)
        if zdata is None:
            print('Timeout waiting for frame')
            self.timedout = True
            return

        if zdata == EOT:
            print('End of transmission')
            self.timedout = True
            return

        try:
            data = zlib.decompress(zdata)
        except:
            print(f'Error decompressing {zdata}')
            raise
        self.framebuffer = np.frombuffer(
            data,
            dtype=np.float32).reshape((-1, 3))
        if self.xtcfile is None:
            self.xtcfile = XTCTrajectoryFile(self.xtcfilename, 'w')
        xyz = self.framebuffer[4:]
        box = self.framebuffer[1:4]
        t = self.framebuffer[0, 0]
        self.saved_frames += 1
        self.xtcfile.write(xyz, time=t, step=self.saved_frames, box=box)

    def close(self):
        self.xtcfile.close()
        self._client.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class TiosNCWriter():
    def __init__(self, sim_id, ncfilename, mqtt_broker=None, port=None,
                 timeout=60):
        # Set the broker address and port
        self.broker_address = mqtt_broker or config.broker
        self.port = port or config.port
        self.ncfilename = ncfilename
        self.timeout = timeout

        self._client = TiosSubscriber(sim_id,
                                      broker_address=self.broker_address,
                                      port=port)
        self.timedout = False
        self.ncfile = None
        self.framebuffer = None
        self.saved_frames = 0

    def write_frame(self):
        
        zdata = interruptable_get(self._client, timeout=self.timeout)
        if zdata is None:
            print('Timeout waiting for frame')
            self.timedout = True
            return

        if zdata == EOT:
            print('End of transmission')
            self.timedout = True
            return
        try:
            data = zlib.decompress(zdata)
        except:
            print('Error decompressing {zdata}')
            raise
        self.framebuffer = np.frombuffer(
            data,
            dtype=np.float32).reshape((-1, 3))
        if self.ncfile is None:
            self.ncfile = NetCDFTrajectoryFile(self.ncfilename, 'w')
        xyz = self.framebuffer[4:]
        box = self.framebuffer[1:4]
        a, b, c, alpha, beta, gamma = box_vectors_to_lengths_and_angles(*box)
        t = self.framebuffer[0, 0]
        self.saved_frames += 1
        self.ncfile.write(xyz, time=t,
                          cell_lengths=(a, b, c),
                          cell_angles=(alpha, beta, gamma))

    def close(self):
        self.ncfile.close()
        self._client.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def get_simulations(broker_address=None, port=None, timeout=10):

    broker_address = broker_address or config.broker
    port = port or config.port

    client = TiosMonitor(broker_address=broker_address, port=port)
    sleep(timeout)
    simulations = {}
    for k in client.status:
        simulations[k] = {'status': client.status[k].decode(), 'summary': '(not available)'}
    for k in client.summary:
        if not k in simulations:
            simulations[k] = {'status': '(unknown)', 'summary': client.summary[k].decode()}
        else:
            simulations[k]['summary'] = client.summary[k].decode()
    
    return simulations
