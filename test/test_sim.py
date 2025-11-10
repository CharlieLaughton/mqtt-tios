#!/usr/bin/env python3

"""
Test script to run an OpenMM simulation and publish snapshots via MQTT.
"""

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, Simulation
from openmm.app import HBonds, PME
from openmm.unit import nanometer, picosecond, kelvin, bar

inpcrdfile = '5fdr_A.inpcrd'
prmtopfile = '5fdr_A.prmtop'


inpcrd = AmberInpcrdFile(inpcrdfile)
prmtop = AmberPrmtopFile(prmtopfile, periodicBoxVectors=inpcrd.boxVectors)

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)
system.addForce(MonteCarloBarostat(1*bar, 310*kelvin))

integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond,
                                      0.002*picosecond)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
print('System built, now minimizing energy...')

simulation.minimizeEnergy()
print('Done.')

""" This is where the mqtt_tios integration starts """
from mqtt_tios import ommutils
sim_id = '5fdr_A'
# mqtt_broker = 'localhost'
mqtt_broker = 'broker.hivemq.com'
mqtt_port = 1883
save_int = 500
print(f'Mqtt-reporter will publish snapshots with id "{sim_id}" to broker "{mqtt_broker}" on port {mqtt_port}.')
client = ommutils.TiosMqttClient(mqtt_broker, port=mqtt_port)
client.register_simulation(
    sim_id,
    simulation,
    summary='MD simulation of 5fdr_A protein in water')
print('Simulation registered with mqtt_tios broker.')
mqtt_reporter = client.create_reporter(save_int)
simulation.reporters.append(mqtt_reporter)
print('Mqtt-reporter added to simulation.')
print('Snapshots will be published every', save_int, 'steps.')
""" End of mqtt_tios integration """

n_steps = int(input('Enter number of MD steps to run (0 = terminate): '))
while n_steps > 0:
    print(f'Running {n_steps} MD steps, saving data every {save_int} steps')
    simulation.step(n_steps)
    n_steps = int(input('Enter no. of further steps to run (0 = terminate): '))

""" A final mqtt_tios cleanup, removing the simulation from the registry """
mqtt_reporter.close()
