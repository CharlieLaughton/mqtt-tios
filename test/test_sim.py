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
report_interval = 500
checkpoint_interval = 5000
mqtt_reporter = ommutils.TiosMqttReporter(
    sim_id,
    report_interval,
    checkpointInterval=checkpoint_interval,
    exists_ok=True)

simulation.reporters.append(mqtt_reporter)
""" End of mqtt_tios integration """

n_steps = int(input('Enter number of MD steps to run (0 = terminate): '))
while n_steps > 0:
    print(f'Running {n_steps} MD steps, saving data every {report_interval} steps')
    simulation.step(n_steps)
    n_steps = int(input('Enter no. of further steps to run (0 = terminate): '))

""" A final mqtt_tios cleanup, removing the simulation from the registry """
mqtt_reporter.close()
