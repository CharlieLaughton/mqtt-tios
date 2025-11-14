#!/usr/bin/env python3

"""
Test script to run an OpenMM simulation and publish snapshots via MQTT.

Demonstrates restarting the simulation from a saved state.

Assumes that 'test_sim.py' has been run previously to register the simulation
with the MQTT broker.
"""

from mqtt_tios import ommutils
sim_id = '5fdr_A'
report_interval = 500
checkpoint_interval = 5000
mqtt_reporter = ommutils.TiosMqttReporter(
    sim_id,
    report_interval,
    checkpointInterval=checkpoint_interval,
    exists_ok=True)
simulation = ommutils.retrieve_simulation(sim_id)
print(f'Simulation {sim_id} retrieved.')

simulation.reporters.append(mqtt_reporter)
print('Mqtt-reporter added to simulation.')
print('Snapshots will be published every', report_interval, 'steps.')

n_steps = int(input('Enter number of MD steps to run (0 = terminate): '))
while n_steps > 0:
    print(f'Running {n_steps} MD steps, saving data every {report_interval} steps')
    simulation.step(n_steps)
    n_steps = int(input('Enter no. of further steps to run (0 = terminate): '))

""" A final mqtt_tios cleanup, removing the simulation from the registry """
mqtt_reporter.close()
