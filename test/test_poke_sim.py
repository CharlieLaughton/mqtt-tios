#!/usr/bin/env python3

"""
Test script to run an OpenMM simulation and publish snapshots via MQTT.

This version demonstrates using a client to 'poke' the simulation into action.

Before you run this, make sure you have access to an MQTT broker and
have set the environment variables TIOSBROKER and TIOSPORT, e.g.:

export TIOSBROKER=localhost
export TIOSPORT=1883

If your broker needs username/password authentication, also set:

export TIOSUSERNAME=<username>
export TIOSPASWORD=<password>

"""
# The start is just standard OpenMM:

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

# This is where the mqtt_tios integration starts

from mqtt_tios import ommutils, clients
import time

sim_id = '5fdr_A'
report_interval = 500
checkpoint_interval = 5000

# Create a client and wait a couple of seconds
# for it to be populated with data from the broker
client = clients.TiosPublisher(sim_id)
time.sleep(2)

# Create an OpenMM Reporter that is "subscriber sensitive":
tios_reporter = ommutils.TiosMqttReporter(
    client,
    report_interval,
    checkpointInterval=checkpoint_interval,
    exists_ok=True,
    on_demand=True)

simulation.reporters.append(tios_reporter)

n_steps = 20000
while True:
        print(f'Running {n_steps} MD steps, saving data every {report_interval} steps')
        simulation.step(n_steps)

