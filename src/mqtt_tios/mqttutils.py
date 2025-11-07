import paho.mqtt.client as mqtt
from queue import SimpleQueue, Empty


class MqttWriter():
    """Class for writing messages to an MQTT broker."""

    def __init__(self, broker_address, default_topic,
                 port=1883, verbose=False,
                 username=None, password=None):
        """Initialize the MQTT writer.

        Args:
            broker_address (str): Address of the MQTT broker.
            default_topic (str): Default topic to publish messages to.
            port (int): Port number of the MQTT broker.
            verbose (bool): If True, print debug information.
            username (str): Username for MQTT authentication.
            password (str): Password for MQTT authentication.
        """
        self._broker_address = broker_address
        self._default_topic = default_topic
        self._port = port
        self.verbose = verbose

        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2)
        self._client.on_connect = self._on_connect
        self._client.username_pw_set(username=username, password=password)

        # Connect to the broker
        self._client.connect(self._broker_address, self._port, 60)
        self._client.loop_start()

    def _on_connect(self, client, userdata, mid, reasoncode, properties):
        if reasoncode != 0:
            raise ConnectionError(f'Error - connection failed, reason code={reasoncode}')

    def writemessage(self, payload, topic=None, retain=False):
        """Write a message to the MQTT broker.
        Args:
            payload (bytes): The message payload to send.
            topic (str): The topic to publish the message to. If None, uses default topic.
            retain (bool): If True, the message will be retained by the broker.
        """
        if topic is None:
            topic = self._default_topic
        result = self._client.publish(topic, payload, retain=retain)
        if result[0] != mqtt.MQTT_ERR_SUCCESS:
            raise ConnectionError(f'Error - publish failed, result={result[0]}')
        if self.verbose:
            print(f'wrote message to {topic}: {len(payload)} bytes')

    def close(self):
        """Close the MQTT connection."""
        self._client.disconnect()
        self._client.loop_stop()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class MqttReader():
    """Class for reading messages from an MQTT broker."""

    def __init__(self, broker_address, topic,
                 port=1883, verbose=False,
                 username=None, password=None,
                 patient=True, timeout=60):
        """Initialize the MQTT reader.

        Args:
            broker_address (str): Address of the MQTT broker.
            topic (str): Topic to subscribe to.
            port (int): Port number of the MQTT broker.
            verbose (bool): If True, print debug information.
            username (str): Username for MQTT authentication.
            password (str): Password for MQTT authentication.
            patient (bool): If True, wait indefinitely for the first message.
            timeout (int): Timeout in seconds for reading messages.
        """
        self._broker_address = broker_address
        self._topic = topic
        self._port = port
        self.verbose = verbose
        self.patient = patient
        self.timeout = timeout

        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2)
        self._client.on_connect = self._on_connect
        self._client.on_message = self._on_message
        self._client.username_pw_set(username=username, password=password),

        # Connect to the broker
        self._client.connect(self._broker_address, self._port, 60)
        self._client.loop_start()
        result = self._client.subscribe(self._topic)
        if result[0] != mqtt.MQTT_ERR_SUCCESS:
            raise ConnectionError(f'Error - subscription failed, result={result[0]}')
        self._queue = SimpleQueue()
        self.timedout = False
        self.firstmessage = True

    def _on_connect(self, client, userdata, mid, reasoncode, properties):
        if self.verbose:
            print(f"Connected with reason code {reasoncode}")
        if reasoncode != 0:
            raise ConnectionError(f'Error - connection failed, reason code={reasoncode}')

    def _on_message(self, client, userdata, msg):
        if self.verbose:
            print(f"Received message on topic '{msg.topic}': {len(msg.payload)} bytes")
        self._queue.put(msg)

    def readmessage(self, timeout=None):
        """Read a message from the MQTT broker.
        Args:
            timeout (int): Timeout in seconds for reading the message.
        Returns:
            msg (MQTTMessage or None): The received message or None if timed out.
        """
        if timeout is None:
            timeout = self.timeout
        if self.firstmessage and self.patient:
            timeout = None
            self.firstmessage = False
        self.timedout = False
        try:
            msg = self._queue.get(timeout=timeout)
        except Empty:
            self.timedout = True
            msg = None
        return msg

    def close(self):
        """Close the MQTT connection."""
        self._client.disconnect()
        self._client.loop_stop()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
