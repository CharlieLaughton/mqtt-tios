# MQTT-Tios

MQTT-Tios is a version of [Tios](https://bitbucket.org/claughton/tios/wiki/Overview) that communicates using the indusstry-standard
[MQTT protocol](https://mqtt.org) for IoT devices.

This is a complete rewrite of the original Tios code, with currently a very different user interface. The focus is on simulations using [OpenMM](https://openmm.org).


## Thumbnail sketch

1. Identify an mqtt broker you can make use of: there are free ones out there that you can use for testing; installing a local [mosquitto broker](https://mosquitto.org/download/) is easy too.
2. Clone/download this repo, `cd` into it and then install the `mqtt-tios` Python package:

  git clone
  cd mqtt-tios
  pip install .
  
4. Set environment variables for your broker:
  export TIOSBROKER=localhost # or wherever yur broker is
  export TIOSPORT=1883 # or whatever port it is using
5. Add a tios-reporter to your OpenMM simulation, e.g.:
6. 
   ...
   from mqtt-tios import ommutils, clients
   simId = 'test_simulation' # your choice of name
   tios_client = clients.TiosPublisher(simId)
   reportInterval = 500
   tios_reporter = ommutils.TiosMqttReporter(
        client,
        reportInterval,
        checkpointInterval=5000,
        exists_ok=True)
   simulation.reporters.append(tios_reporter)
   ...
   
8. Once your simulation is running, it will publish snaphots to the broker.
9. From some other terminal window or computer with mqtt-tios installed and the environment variables set, connect to the running simulation and save published snapshots to a trajectory file:
    
  tios_collect test_simulation test_sim.nc
  49 frames [00:52,  1.05s/ frames]
  
11. Hit `Ctrl-C` to quit collecting when you have enough. Frames can be saved in GROMACS .xtc or AMBER .nc format.
12. If the simulation finishes, the collector will stall but not quit. If the simulation is restarted, collection will restart.

## Examples
See `test_sim.py` and `restart_sim.py` in the `/tests` folder for inspiration.
   

