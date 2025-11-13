#!/usr/bin/env python3

"""
Test script to run an OpenMM simulation and publish snapshots via MQTT.

Demonstrates restarting the simulation from a saved state.

Assumes that 'test_sim.py' has been run previously to register the simulation
with the MQTT broker.
"""

from mqtt_tios import ommutils
sim_id = '5fdr_A'
# mqtt_broker = 'localhost'
mqtt_broker = 'charlielaughton.com'
# mqtt_port = 1883
mqtt_port = 8080 
save_int = 500
client = ommutils.TiosMqttClient(mqtt_broker, port=mqtt_port, username='tios_publisher', password='publisher_tios')
simulation = client.retrieve_simulation(sim_id)
print(f'Simulation {sim_id} retrieved from mqtt_tios broker {mqtt_broker}.')
print(f'client simId = {client.simId}')

mqtt_reporter = client.create_reporter(save_int)
simulation.reporters.append(mqtt_reporter)
print('Mqtt-reporter added to simulation.')
print('Snapshots will be published every', save_int, 'steps.')

n_steps = int(input('Enter number of MD steps to run (0 = terminate): '))
while n_steps > 0:
    print(f'Running {n_steps} MD steps, saving data every {save_int} steps')
    simulation.step(n_steps)
    n_steps = int(input('Enter no. of further steps to run (0 = terminate): '))

""" A final mqtt_tios cleanup, removing the simulation from the registry """
mqtt_reporter.close()
