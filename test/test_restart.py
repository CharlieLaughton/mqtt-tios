#!/usr/bin/env python3

"""
Test script to run an OpenMM simulation and publish snapshots via MQTT.

Demonstrates restarting the simulation from a saved state.

Assumes that 'test_sim.py' has been run previously to register the simulation
with the MQTT (Tios) broker.
"""

from mqtt_tios import ommutils, clients
import time

sim_id = '5fdr_A'
report_interval = 500
checkpoint_interval = 5000

# Create a client, and allow a couple of seconds for data
# to be retrieved from the broker:
client = clients.TiosPublisher(sim_id)
time.sleep(2)

# Create a reporter:
tios_reporter = ommutils.TiosMqttReporter(
    client,
    report_interval,
    checkpointInterval=checkpoint_interval,
    exists_ok=True)

# Retrieve the latest checkpoint for this simulation:
simulation = ommutils.retrieve_checkpoint(client)
print(f'Simulation {sim_id} retrieved.')

# Add the reporter to the simulation and then start MD cycles:
simulation.reporters.append(tios_reporter)
print('tios-reporter added to simulation.')
print('Snapshots will be published every', report_interval, 'steps.')
print('A checkpoint will be saved every', checkpoint_interval, 'steps.')

n_steps = int(input('Enter number of MD steps to run (0 = terminate): '))
while n_steps > 0:
    print(f'Running {n_steps} MD steps, saving data every {report_interval} steps')
    simulation.step(n_steps)
    n_steps = int(input('Enter no. of further steps to run (0 = terminate): '))

# A final clean-up (optional, but polite):
tios_reporter.close()
