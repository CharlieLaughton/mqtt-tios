import zlib
import json
from io import StringIO

import numpy as np
from .mqttutils import MqttReader, MqttWriter
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


def retrieve_simulation(brokerAddress, simId, port=1883, username=None, password=None):
    """Retrieve an existing simulation.
    Parameters
    ----------
    brokerAddress : str
        The address of the MQTT broker.
    simId : str
        A unique identifier for the simulation.
    port : int, optional
        The port number of the MQTT broker.
    username : str, optional
        The username for MQTT authentication.
    password : str, optional
        The password for MQTT authentication.
    Returns
    ------
    simulation : OpenMM Simulation
        The OpenMM Simulation object for the specified simulation ID.
    """
    with MqttReader(brokerAddress,
                    f'tios/{simId}/checkpoint',
                    port=port,
                    username=username,
                    password=password,
                    client_id='retriever') as f:
        try:
            msg = f.readmessage()
            if msg is None:
                raise ConnectionError('No checkpoint data received.')
            elif msg.timestamp < 0:
                raise ConnectionError('No checkpoint data available for this simulation.')
            data = zlib.decompress(msg.payload)
            simulation = deserialize_simulation(data)
            return simulation
        except ConnectionError as e:
            print('Error retrieving simulation:', e, flush=True,
                  file=sys.stderr)


class TiosMqttReporter():
    """A reporter that sends OpenMM simulation data via MQTT."""
    def __init__(self, brokerAddress, simId, reportInterval,
                 checkpointInterval=None,
                 summary=None,
                 enforcePeriodicBox=None, port=1883,
                 username=None, password=None,
                 exists_ok=False):
        """Initialize the MQTT reporter.
        Parameters
        ----------
        brokerAddress : str
            The address of the MQTT broker.
        simId : str
            A unique identifier for the simulation.
        reportInterval : int
            The interval (in steps) at which to report simulation data.
        checkpointInterval : int, optional
            The interval (in steps) at which to checkpoint the simulation.
        summary : str, optional
            An optional summary description of the simulation.
        enforcePeriodicBox : bool, optional
            If True, wrap coordinates to be within the periodic box.
        port : int, optional
            The port number of the MQTT broker.
        username : str, optional
            The username for MQTT authentication.
        password : str, optional
            The password for MQTT authentication.
        """

        self._reportInterval = reportInterval
        self._checkpointInterval = checkpointInterval
        self._summary = summary
        self._enforcePeriodicBox = enforcePeriodicBox
        self._brokerAddress = brokerAddress
        self._port = port
        self._simId = simId
        self._username = username
        self._password = password

        if self.check_exists(simId):
            self._new = False
            if not exists_ok:
                raise ValueError(f'Simulation ID {simId} already exists.')
        else:
            self._new = True

        self._framebuffer = None
        self._first_report = True
        self._report_writer = MqttWriter(brokerAddress,
                                        f'tios/{simId}/state',
                                        port=port,
                                        username=username,
                                        password=password,
                                        status_topic=f'tios/{simId}/status',
                                        client_id='report_writer')
        if checkpointInterval is not None:
            self._checkpoint_writer = MqttWriter(brokerAddress,
                                                 f'tios/{simId}/checkpoint',
                                                 port=port,
                                                 username=username,
                                                 password=password,
                                                 client_id='checkpoint_writer')

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
            self.register_simulation(self._simId, simulation, summary=self._summary)
            self._new = False

        if self._first_report:
            self.checkpoint_simulation(simulation)
            self._last_checkpoint_step = simulation.currentStep
            self._first_report = False
        elif (self._checkpointInterval is not None and
              simulation.currentStep - self._last_checkpoint_step >= self._checkpointInterval):
            self.checkpoint_simulation(simulation)
            self._last_checkpoint_step = simulation.currentStep

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
        try:
            self._report_writer.writemessage(data)
        except ConnectionError as e:
            print('Error publishing simulation snapshot:', e, flush=True,
                  file=sys.stderr)

    def close(self):
        """ Close the MQTT connection

        Note: this method is not called automatically by OpenMM.

        """
        self._report_writer.close()

    def check_exists(self, simId):
        """Check if a simulation with the given ID exists.

        Parameters
        ----------
        simId : str
            A unique identifier for the simulation.

        Returns
        -------
        bool
            True if the simulation exists, False otherwise.

        """
        with MqttReader(self._brokerAddress,
                        f'tios/{simId}/#',
                        port=self._port,
                        username=self._username,
                        password=self._password,
                        patient=False, timeout=2,
                        client_id='checker') as f:
            try:
                msg = f.readmessage()
                if msg is None:
                    return False
                return True
            except ConnectionError:
                return False

    def register_simulation(self, simId, simulation, summary=None):
        """Register a new simulation.

        Parameters
        ----------
        simId : str
            A unique identifier for the simulation.
        simulation : OpenMM Simulation
            The OpenMM Simulation object to register.
        summary : str, optional
            An optional summary description of the simulation.

        """
        if self.check_exists(simId):
            raise ValueError(f'Simulation ID {simId} is already in use.')
        if summary is None:
            summary = f'OpenMM simulation with {simulation.context.getSystem().getNumParticles()} atoms.'
        with MqttWriter(self._brokerAddress,
                        f'tios/{simId}/summary',
                        port=self._port,
                        username=self._username,
                        password=self._password,
                        client_id='registrar') as f:
            f.writemessage(summary.encode('utf-8'), retain=True)
        self.simId = simId
        self.checkpoint_simulation(simulation)

    def checkpoint_simulation(self, simulation):
        """Checkpoint an existing simulation.
        Parameters
        ----------
        simulation : OpenMM Simulation
            The OpenMM Simulation object to checkpoint.
        """

        data = serialize_simulation(simulation)
        try:
            self._checkpoint_writer.writemessage(zlib.compress(data), retain=True)
        except ConnectionError as e:
            print('Error checkpointing simulation:', e, flush=True,
                  file=sys.stderr)
