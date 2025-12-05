""" tiosutils.py: lower-level routines for command-line tools """

from time import sleep
import zlib
import json
from tempfile import NamedTemporaryFile
import numpy as np
import mdtraj as mdt
from mdtraj.utils import box_vectors_to_lengths_and_angles
from mdtraj.formats import NetCDFTrajectoryFile, XTCTrajectoryFile
from .clients import TiosSubscriber, TiosMonitor
from .config import config
from .control import killer

EOT = 'EOT'.encode('utf-8')


def get_topology(client, timeout=10):
    """ Get the (MDTraj) topology from the publisher via MQTT

    Parameters
    ----------
    client : TiosSubscriber
        The MQTT client subscribed to the simulation.
    timeout : int
        Timeout in seconds to wait for the topology checkpoint.

    Returns
    -------
    topology : mdtraj.Topology
        The MDTraj topology object obtained from the publisher.
    """
    time_left = timeout
    while time_left > 0 and client.checkpoint is None and not killer.kill_now:
        sleep(1)
        time_left -= 1
    if client.checkpoint is None:
        raise TimeoutError('Timeout waiting for topology')
    metadata = json.loads(client.checkpoint)
    with NamedTemporaryFile('w+', suffix='.pdb', delete=True) as f:
        f.write(metadata['pdbdata'])
        f.flush()
        f.seek(0)
        topology = mdt.load_pdb(f.name).topology
    return topology


def _interruptable_get(client, timeout=None):
    """ A get that can be interrupted"""
    time_left = timeout or 1
    while time_left > 0 and client.states.qsize() == 0 and not killer.kill_now:
        sleep(1)
        if timeout:
            time_left -= 1
    if client.states.empty():
        return None
    else:
        return client.states.get()


class TiosPDBWriter():
    """ Write frames from a TIOS simulation to PDB files.

    Parameters
    ----------
    sim_id : str
        The simulation ID to subscribe to.
    pdbfilename : str
        The PDB filename pattern to write frames to. Can include
        a '{frame}' placeholder for the frame number. Otherwise the same
        filename is used for all frames.
    mqtt_broker : str, optional
        The address of the MQTT broker. If None, uses the default from config.
    port : int, optional
        The port number of the MQTT broker. If None, uses the default from
        config.
    timeout : int, optional
        Timeout in seconds for waiting for new frames. Default is 60.

    Methods:
        write_frame(): write the next frame to a PDB file.
        close(): close the MQTT client.

    """
    def __init__(self, sim_id, pdbfilename, mqtt_broker=None, port=None,
                 timeout=60):
        # Set the broker address and port
        self.broker_address = mqtt_broker or config.broker
        self.port = port or config.port
        self.pdbfilename = pdbfilename
        self.timeout = timeout

        self._client = TiosSubscriber(sim_id,
                                      broker_address=self.broker_address,
                                      port=port)
        self.timedout = False
        self.framebuffer = None
        self.saved_frames = 0
        self.topology = get_topology(self._client)

    def write_frame(self):

        zdata = _interruptable_get(self._client, timeout=self.timeout)
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
        except Exception as e:
            print(f'Error decompressing {zdata}')
            raise e
        self.framebuffer = np.frombuffer(
            data,
            dtype=np.float32).reshape((-1, 3))
        xyz = self.framebuffer[4:]
        box = self.framebuffer[1:4]
        a, b, c, alpha, beta, gamma = box_vectors_to_lengths_and_angles(*box)

        traj = mdt.Trajectory(xyz, topology=self.topology,
                              unitcell_lengths=[a, b, c],
                              unitcell_angles=[alpha, beta, gamma])
        traj.save_pdb(self.pdbfilename.format(frame=self.saved_frames))
        self.saved_frames += 1

    def close(self):
        self._client.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


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

        zdata = _interruptable_get(self._client, timeout=self.timeout)
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
        except Exception as e:
            print(f'Error decompressing {zdata}')
            raise e
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
        if self.xtcfile is not None:
            self.xtcfile.close()
        self._client.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class TiosNCWriter():
    """ Write frames from a TIOS simulation to NetCDF files.

    Parameters
    ----------
    sim_id : str
        The simulation ID to subscribe to.
    ncfilename : str
        The NetCDF filename to write frames to.
    mqtt_broker : str, optional
        The address of the MQTT broker. If None, uses the default from config.
    port : int, optional
        The port number of the MQTT broker. If None, uses the default from
        config.
    timeout : int, optional
        Timeout in seconds for waiting for new frames. Default is 60.
    Methods:
        write_frame(): write the next frame to the NetCDF file.
        close(): close the MQTT client and the NetCDF file."""

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
        zdata = _interruptable_get(self._client, timeout=self.timeout)
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
        except Exception as e:
            print(f'Error decompressing {zdata}')
            raise e
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
        if self.ncfile is not None:
            self.ncfile.close()
        self._client.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def get_simulations(broker_address=None, port=None, timeout=10):
    """ Get the list of available simulations from an MQTT broker.
    Parameters
    ----------
    broker_address : str, optional
        The address of the MQTT broker. If None, uses the default from config.
    port : int, optional
        The port number of the MQTT broker. If None, uses the default from
        config.
    timeout : int, optional
        Timeout in seconds to wait for simulations. Default is 10.
    Returns
    -------
    simulations : dict
        A dictionary with simulation IDs as keys and dictionaries with
        'status', 'summary', and 'subscribed' as values.
    """

    broker_address = broker_address or config.broker
    port = port or config.port

    client = TiosMonitor(broker_address=broker_address, port=port)
    sleep(timeout)
    simulations = {}
    for k in client.status:
        simulations[k] = {'status': client.status[k].decode(),
                          'summary': '(not available)',
                          'subscribed': False}
    for k in client.summary:
        if k not in simulations:
            simulations[k] = {'status': '(unknown)',
                              'summary': client.summary[k].decode()}
        else:
            simulations[k]['summary'] = client.summary[k].decode()

    for k in client.subscribed:
        if client.subscribed[k]:
            state = 'running'
        else:
            state = 'stopped'
        if k not in simulations:
            simulations[k] = {'status': state,
                              'summary': '(not available)'}
        else:
            if simulations[k]['status'] == 'online':
                simulations[k]['status'] = state
    client.close()
    return simulations
