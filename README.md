# MQTT-Tios

MQTT-Tios is a version of [Tios](https://bitbucket.org/claughton/tios/wiki/Overview) that communicates using the industry-standard
[MQTT protocol](https://mqtt.org) for IoT devices.

This is a complete rewrite of the original Tios code, with currently a very different user interface. The focus is on simulations using [OpenMM](https://openmm.org).


MQTT-Tios provides:

1. A `Reporter` class for OpenMM simulations that, rather than saving snapshots of the simulation to a file, publishes them to the internet.
2. Two command-line utilities: `tios_ls` to list streaming simulations currently running anywhere; and `tios_collect` to tap into a stream and collect snapshots as a conventional trajectory file.

## How to use

1. Identify an mqtt broker you can make use of: there are free ones out there that you can use for testing; installing a local [mosquitto broker](https://mosquitto.org/download/) is easy too.
2. Clone/download this repo, `cd` into it and then install the `mqtt-tios` Python package:
```
  git clone git@github.com:CharlieLaughton/mqtt-tios.git
  cd mqtt-tios
  pip install .
  ```
4. Set environment variables for your broker:
```
  export TIOSBROKER=localhost # or wherever yur broker is
  export TIOSPORT=1883 # or whatever port it is using
```
5. Add a tios-reporter to your OpenMM simulation, e.g.:
   ```
   ...
   from mqtt-tios import ommutils, clients
   simId = 'test_simulation' # your choice of name
   tios_client = clients.TiosPublisher(simId)
   reportInterval = 500
   tios_reporter = ommutils.TiosMqttReporter(
        tios_client,
        reportInterval,
        checkpointInterval=5000,
        exists_ok=True)
   simulation.reporters.append(tios_reporter)
   ...
   ```
6. Once your simulation is running, it will publish snaphots to the broker.
7. It also saves checkpoints of the whole simulation state at intervals, this permits simulations to be restarted seamlessly, even on different resources from the ones they were originally running on.
8. From some other terminal window or computer with mqtt-tios installed and the environment variables set, use `tios_ls` to check the simulations available, then connect to the running simulation and save published snapshots to a trajectory file:
  ```
  % tios_ls
  Available simulations:
   - test_simulation: OpenMM simulation with 20962 atoms.: online
   - other_simulation: OpenMM simulation with 16545 atoms: offline

  % tios_collect test_simulation test_sim.nc
  49 frames [00:52,  1.05s/ frames]
  ```
8. Hit `Ctrl-C` to quit collecting when you have enough. Frames can be saved in GROMACS .xtc or AMBER .nc format.
9. If the simulation finishes, the collector will stall but not quit. If the simulation is restarted, collection will restart.

## Examples
See `test_sim.py` and `restart_sim.py` in the `/tests` folder for inspiration.

## Considerations
Anyone who can access your mqtt broker can see and collect data from your simulations. Mqtt-tios supports mqtt's username/password authentication and access control mechanisms so if you have control over your broker you can limit write access to chosen users, but it does not currently support SSL/TLS security. 
   

