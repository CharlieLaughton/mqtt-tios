# MQTT-Tios

**MQTT-Tios** is a version of [Tios](https://bitbucket.org/claughton/tios/wiki/Overview) that communicates using the industry-standard
[MQTT protocol](https://mqtt.org) for IoT devices.

This is a complete rewrite of the original Tios code, with currently a very different user interface. The focus is on simulations using [OpenMM](https://openmm.org).


**MQTT-Tios** provides:

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
  export TIOSBROKER=localhost # or wherever your broker is
  export TIOSPORT=1883 # or whatever port it is using
```
5. Add a tios-reporter to your OpenMM simulation code, e.g.:
   ```
   ...
   from mqtt-tios import ommutils, clients
   # Create a client:
   simId = 'test_simulation' # your choice of name
   tios_client = clients.TiosPublisher(simId)
   # Create a reporter:
   reportInterval = 500
   tios_reporter = ommutils.TiosMqttReporter(
        tios_client,
        reportInterval,
        exists_ok=True)
   # Attach the reporter to your simulation:
   simulation.reporters.append(tios_reporter)
   ...
   ```
6. Start your simulation in the normal way.
7. From some other terminal window or computer with `mqtt-tios` installed and the environment variables set, connect to the running simulation and save published snapshots to a trajectory file:
  ```
  % tios_collect test_simulation test_sim.nc
  49 frames [00:52,  1.05s/ frames]
  ```
8. Hit `Ctrl-C` to quit collecting when you have enough. Frames can be saved in GROMACS .xtc or AMBER .nc format. Single snapshots can be saved in .pdb format.

## "On demand" simulations
In the standard pub/sub model, publishers publish data on topics irrespective of whether there are any subscribers to that topic. However **MQTT-Tios** uses an "on demand" model - 
simulations are held in a sleeping state until at least one subscriber is registered, so expensive MD data is not wasted. The `tios_ls` command provides the state of all 
registered simulations:

```
    % tios_ls
Available simulations:
   - test_simulation: OpenMM simulation with 20962 atoms.: stopped
   - second_simulation: OpenMM simulation with 31415 atoms.: running
   - other_simulation: OpenMM simulation with 16545 atoms: offline
```
"Offline" simulations are those not currently attached to an OpenMM simulator; "stopped" ones are waiting for their first subscriber to register, while "running" simulations
are already delivering snapshots to a collector somewhere. New collectors can be attached to "running" simulations as well as "stopped" ones.

## Checkpointing and restarting
The **MQTT-Tios** reporter can also publish checkpoints of the entire simulation state at intervals. This allows simulations to be restarted seamlessly. If a collector is subscribed to 
a simulation when it ends, it can wait until the restart occurs and then pick up where it left off, appending further snapshots to the same trajectory file.

## Examples
See `test_sim.py` and `restart_sim.py` in the `/tests` folder for inspiration.

## Considerations
Anyone who can access your mqtt broker can see and collect data from your simulations. **MQTT-Tios** supports mqtt's username/password authentication and access control mechanisms so if you have control over your broker you can limit write access to chosen users, but it does not currently support SSL/TLS security. 


