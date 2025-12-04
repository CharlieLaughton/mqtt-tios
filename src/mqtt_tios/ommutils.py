""" ommutils.py: tios integrations for OpenMM simulations """
from time import sleep
import zlib
import json
from io import StringIO

import numpy as np
from .clients import TiosPublisher
import sys

from openmm import XmlSerializer
from openmm.app import PDBFile, Simulation
from openmm.unit import picosecond

sys.tracebacklimit = None  # suppress traceback for ConnectionError


def serialize_simulation(simulation):
    """Serialize an OpenMM simulation.

    Used for reconstructing simulations.
    Args:
        simulation: An OpenMM simulation object.
    Returns:
        bytes: serialized simulation with current state.

    """
    st = simulation.context.getState(positions=True)
    f = StringIO('')
    PDBFile.writeFile(simulation.topology, st.getPositions(), f)
    f.seek(0)
    pdbdata = f.read()

    f = StringIO('')
    simulation.saveState(f)
    f.seek(0)
    statedata = f.read()

    metadata = {'integrator': XmlSerializer.serialize(simulation.integrator),
                'system': XmlSerializer.serialize(simulation.system),
                'pdbdata': pdbdata,
                'statedata': statedata
                }

    return json.dumps(metadata).encode('utf-8')


def deserialize_simulation(simulation_data):
    """Get the uncompressed simulation from metadata.
    Args:
        simulation_data (bytes): serialized simulation data.
    Returns:
        simulation: An OpenMM simulation object.

    """
    metadata = json.loads(simulation_data)
    integrator = XmlSerializer.deserialize(metadata['integrator'])
    system = XmlSerializer.deserialize(metadata['system'])
    f = StringIO(metadata['pdbdata'])
    f.seek(0)
    pdb = PDBFile(f)
    simulation = Simulation(pdb.topology, system, integrator)
    f = StringIO(metadata['statedata'])
    f.seek(0)
    simulation.loadState(f)
    return simulation


def retrieve_checkpoint(client):
    """Retrieve a checkpoint of a simulation from the broker"""
    if client.checkpoint is None:
        raise ValueError(f'Error: Simulation {client.simId} has no checkpoint')
    return deserialize_simulation(client.checkpoint)


class TiosMqttReporter():
    """A reporter that sends OpenMM simulation data via MQTT."""
    def __init__(self, client, reportInterval,
                 checkpointInterval=None,
                 summary=None,
                 enforcePeriodicBox=None,
                 exists_ok=False,
                 on_demand=True,
                 verbose=False):
        """Initialize the MQTT reporter.
        Parameters
        ----------
        client : TiosPublisher
            The MQTT client
        reportInterval : int
            The interval (in steps) at which to report simulation data.
        checkpointInterval : int, optional
            The interval (in steps) at which to checkpoint the simulation.
        summary : str, optional
            An optional summary description of the simulation.
        enforcePeriodicBox : bool, optional
            If True, wrap coordinates to be within the periodic box.
        exists_ok : bool, optional
            If True, any existing data for this simulation will be overwritten
        on_demand : bool, optional
            If True, the simulation will pause if no subscribers are connected.
        """
        if not isinstance(client, TiosPublisher):
            raise ValueError('Error: client must be a TiosPublisher')

        if checkpointInterval is not None:
            if checkpointInterval % reportInterval != 0:
                raise ValueError('Error: checkpointInterval must be a '
                                 'multiple of reportInterval')
        self._client = client
        self._reportInterval = reportInterval
        self._checkpointInterval = checkpointInterval
        self._summary = summary
        self._enforcePeriodicBox = enforcePeriodicBox
        self._on_demand = on_demand
        self.verbose = verbose

        if self._client.summary is not None:
            self._new = False
            if not exists_ok:
                raise ValueError(
                    f'Simulation ID {self._client.simId} already exists.')
        else:
            self._new = True

        self._framebuffer = None
        self._first_report = True

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        dict
            A dictionary describing the required information for the
            next report
        """
        ri = self._reportInterval
        steps = ri - simulation.currentStep % ri
        return {'steps': steps,
                'periodic': self._enforcePeriodicBox,
                'include': ['positions']}

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation

        """
        if self._new:
            self.register_simulation(simulation, summary=self._summary)
            self._new = False

        if self._first_report:
            self.checkpoint_simulation(simulation)
            self._last_checkpoint_step = simulation.currentStep
            self._first_report = False
        elif (self._checkpointInterval is not None and
              (simulation.currentStep - self._last_checkpoint_step
               >= self._checkpointInterval)):
            self.checkpoint_simulation(simulation)
            self._last_checkpoint_step = simulation.currentStep
        if self._on_demand:
            if self.verbose and not self._client.has_subscribers:
                print('Waiting for subscribers to connect...')
                while not self._client.has_subscribers:
                    sleep(1)
                print('Subscriber connected, resuming simulation.')
        positions = state.getPositions(asNumpy=True)
        box = state.getPeriodicBoxVectors(asNumpy=True)
        t = state.getTime().value_in_unit(picosecond)

        if self._framebuffer is None:
            n_atoms, _ = positions.shape
            self._framebuffer = np.zeros((n_atoms + 4, 3), dtype=np.float32)
        self._framebuffer[0, 0] = t
        self._framebuffer[1:4] = box
        self._framebuffer[4:] = positions
        data = zlib.compress(self._framebuffer.tobytes())
        self._client.state = data

    def close(self):
        """ Close the MQTT connection

        Note: this method is not called automatically by OpenMM.

        """
        self._client.close()

    def register_simulation(self, simulation, summary=None):
        """Register a new simulation.

        Parameters
        ----------
        simulation : OpenMM Simulation
            The OpenMM Simulation object to register.
        summary : str, optional
            An optional summary description of the simulation.

        """
        if summary is None:
            n_atoms = simulation.context.getSystem().getNumParticles()
            summary = f'OpenMM simulation with {n_atoms} atoms.'
        self._client.summary = summary.encode('utf-8')
        self.simId = self._client.simId
        self.checkpoint_simulation(simulation)

    def checkpoint_simulation(self, simulation):
        """Checkpoint an existing simulation.
        Parameters
        ----------
        simulation : OpenMM Simulation
            The OpenMM Simulation object to checkpoint.
        """

        data = serialize_simulation(simulation)
        self._client.checkpoint = data
