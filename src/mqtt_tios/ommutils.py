import time
import zlib
import json
from io import StringIO

import numpy as np
from .mqttutils import MqttWriter
import sys

from openmm import XmlSerializer, LangevinMiddleIntegrator
from openmm.app import PDBFile, Simulation
from openmm.unit import kelvin, picosecond

sys.tracebacklimit = None  # suppress traceback for ConnectionError


def get_compressed_state(simulation):
    """Get the compressed state of the simulation.

    Used for checkpointing and restarting simulations.

    Args:
        simulation: An OpenMM simulation object.
    Returns:
        bytes: Compressed state data.

    """
    f = StringIO()
    simulation.saveState(f)
    f.seek(0)
    state = {
        'statedata': f.read(),
        'timestamp': time.time()
    }
    statez = zlib.compress(json.dumps(state).encode('utf-8'))
    return statez


def get_compressed_simulation(simulation):
    """Get the compressed metadata of the simulation.

    Used for reconstructing simulations.
    Args:
        simulation: An OpenMM simulation object.
    Returns:
        bytes: Compressed metadata.

    """
    st = simulation.context.getState(positions=True)
    f = StringIO('')
    PDBFile.writeFile(simulation.topology, st.getPositions(), f)
    f.seek(0)
    pdbdata = f.read()

    if isinstance(simulation.integrator, LangevinMiddleIntegrator):
        integrator_code = 'LMI'
    else:
        integrator_code = 'unknown'

    metadata = {'integrator': integrator_code,
                'temperature': simulation.integrator.getTemperature() / kelvin,
                'friction': simulation.integrator.getFriction() * picosecond,
                'timestep': simulation.integrator.getStepSize() / picosecond,
                'system': XmlSerializer.serialize(simulation.system),
                'pdb': pdbdata
                }

    metadataz = zlib.compress(json.dumps(metadata).encode('utf-8'))
    return metadataz


def get_uncompressed_simulation(simulationz):
    """Get the uncompressed simulation from metadata.
    Args:
        simulationz (bytes): Compressed simulation metadata.
    Returns:
        simulation: An OpenMM simulation object.

    """
    metadata = json.loads(zlib.decompress(simulationz))
    assert metadata['integrator'] == 'LMI'
    integrator = LangevinMiddleIntegrator(metadata['temperature']*kelvin,
                                          metadata['friction']/picosecond,
                                          metadata['timestep']*picosecond)
    system = XmlSerializer.deserialize(metadata['system'])
    f = StringIO(metadata['pdb'])
    f.seek(0)
    pdb = PDBFile(f)
    simulation = Simulation(pdb.topology, system, integrator)
    return simulation


def get_uncompressed_state(statez):
    """Get the uncompressed state from a compressed state.

    Args:
        statez (bytes): Compressed state data.

    Returns:
        statedata: The uncompressed state data.
        timestamp: The timestamp of the state.

    """
    state = json.loads(zlib.decompress(statez))
    statedata = state['statedata']
    ts = state['timestamp']
    return statedata, ts


class MqttReporter():
    """A reporter that sends OpenMM simulation data via MQTT."""
    def __init__(self, brokerAddress, simId, reportInterval,
                 summary=None, enforcePeriodicBox=None, port=1883,
                 username=None, password=None):
        """Initialize the MQTT reporter.
        Parameters
        ----------
        brokerAddress : str
            The address of the MQTT broker.
        simId : str
            A unique identifier for the simulation.
        reportInterval : int
            The interval (in steps) at which to report simulation data.
        summary : str, optional
            A summary description of the simulation.
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
        self._enforcePeriodicBox = enforcePeriodicBox
        self._brokerAddress = brokerAddress
        self._port = port
        self._simId = simId
        self._summary = summary

        self._framebuffer = None
        self.sim_topic = f'tios/{self._simId}/simulation'
        self.state_topic = f'tios/{self._simId}/state'
        self.summ_topic = f'tios/{self._simId}/summary'
        try:
            self._writer = MqttWriter(self._brokerAddress, self.state_topic,
                                      port=self._port,
                                      username=username, password=password,
                                      client_id=self._simId)
        except ConnectionError as e:
            raise e.with_traceback(None) from None

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
        positions = state.getPositions(asNumpy=True)
        box = state.getPeriodicBoxVectors(asNumpy=True)
        time = state.getTime().value_in_unit(picosecond)

        if self._framebuffer is None:
            simz = get_compressed_simulation(simulation)
            try:
                self._writer.writemessage(simz, topic=self.sim_topic,
                                          retain=True)
            except ConnectionError as e:
                print('Error saving simulation:', e, flush=True,
                      file=sys.stderr)
            print('saved simulation')
            n_atoms, _ = positions.shape
            if self._summary is None:
                self._summary = f'OpenMM simulation of {self._simId},'
                self._summary += f' n_atoms={n_atoms}'
            try:
                self._writer.writemessage(self._summary, topic=self.summ_topic,
                                          retain=True)
            except ConnectionError as e:
                print('Error saving summary:', e, flush=True, file=sys.stderr)
            print('saved summary')
            self._framebuffer = np.zeros((n_atoms + 4, 3), dtype=np.float32)
        self._framebuffer[0, 0] = time
        self._framebuffer[1:4] = box
        self._framebuffer[4:] = positions
        data = zlib.compress(self._framebuffer.tobytes())
        try:
            self._writer.writemessage(data)
        except ConnectionError as e:
            print('Error saving simulation snapshot:', e, flush=True,
                  file=sys.stderr)

    def close(self):
        """ Close the MQTT connection

        Removes the retained simulation and summary data too.
        Note: this method is not called automatically by OpenMM.

        """
        self._writer.writemessage('', topic=self.sim_topic, retain=True)
        self._writer.writemessage('', topic=self.summ_topic, retain=True)
        self._writer.close()
