import paho.mqtt.client as mqtt
from queue import SimpleQueue, Empty
import string
import random
import time
import struct

# Generate a random client ID
alphabet = string.ascii_lowercase + string.digits


def random_choice():
    return ''.join(random.choices(alphabet, k=6))


def add_timestamp(payload):
    """Add a timestamp to the payload."""
    if not isinstance(payload, bytes):
        raise ValueError('Payload must be bytes')
    return payload + struct.pack('f', time.time())


def add_eot(payload):
    """Add an end-of-transmission (EOT) marker to the payload.

    This is a negative-values timestamp."""
    if not isinstance(payload, bytes):
        raise ValueError('Payload must be bytes')
    return payload + struct.pack('f', -time.time())


def remove_timestamp(payload):
    """Remove the timestamp from the payload."""
    if not isinstance(payload, bytes) or len(payload) < 4:
        raise ValueError('Payload must be bytes with at least 4 bytes')
    data = payload[:-4]
    timestamp_bytes = payload[-4:]
    timestamp = struct.unpack('f', timestamp_bytes)[0]
    return data, timestamp


class MqttMessage():
    """Class representing an MQTT message with topic and payload."""
    def __init__(self, topic, payload, timestamp):
        self.topic = topic
        self.payload = payload
        self.timestamp = timestamp

    def __repr__(self):
        return f'MqttMessage(topic={self.topic}, payload_length={len(self.payload)}, timestamp={self.timestamp})'

    def __str__(self):
        return self.__repr__()


class MqttWriter():
    """Class for writing messages to an MQTT broker."""

    def __init__(self, broker_address, default_topic,
                 port=1883, verbose=False,
                 username=None, password=None,
                 client_id=None):
        """Initialize the MQTT writer.

        Args:
            broker_address (str): Address of the MQTT broker.
            default_topic (str): Default topic to publish messages to.
            port (int): Port number of the MQTT broker.
            verbose (bool): If True, print debug information.
            username (str): Username for MQTT authentication.
            password (str): Password for MQTT authentication.
            client_id (str): Client ID for MQTT connection, a random ID is appended.
        """
        self._broker_address = broker_address
        self._default_topic = default_topic
        self._port = port
        self.verbose = verbose
        if client_id is None:
            client_id = 'mqtt_writer-' + random_choice()
        else:
            client_id += '-' + random_choice()

        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2,
                                   client_id=client_id)
        self._client.on_connect = self._on_connect
        self._client.username_pw_set(username=username, password=password)
        self._client.will_set(self._default_topic, payload=add_eot(b''), retain=True)

        # Connect to the broker
        self._client.connect(self._broker_address, self._port, 60)
        self._client.loop_start()

    def _on_connect(self, client, userdata, mid, reasoncode, properties):
        if reasoncode != 0:
            raise ConnectionError(f'Error - connection failed, reason code={reasoncode}')

    def writemessage(self, payload, topic=None, retain=False):
        """Write a message to the MQTT broker.
        Args:
            payload (bytes): message payload to send.
            topic (str): The topic to publish the message to. If None, uses default topic.
            retain (bool): If True, the message will be retained by the broker.
        """
        if topic is None:
            topic = self._default_topic
        pt = add_timestamp(payload)
        result = self._client.publish(topic, pt, retain=retain)
        if result[0] != mqtt.MQTT_ERR_SUCCESS:
            raise ConnectionError(f'Error - publish failed, result={result[0]}')
        if self.verbose:
            print(f'wrote message to {topic}: {len(payload)} bytes')

    def writeeot(self, payload, topic=None):
        """Write an end-of-transmission (EOT) message to the MQTT broker.
        Args:
            topic (str): The topic to publish the EOT message to. If None, uses default topic.
        """
        if topic is None:
            topic = self._default_topic
        pt = add_eot(payload)
        result = self._client.publish(topic, pt, retain=True)
        if result[0] != mqtt.MQTT_ERR_SUCCESS:
            raise ConnectionError(f'Error - publish failed, result={result[0]}')
        if self.verbose:
            print(f'wrote EOT message to {topic}')

    def close(self):
        """Close the MQTT connection."""
        self.writeeot(b'')
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
                 patient=True, timeout=60,
                 client_id=None):
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
            client_id (str): Client ID for MQTT connection, a random ID is appended.
        """
        self._broker_address = broker_address
        self._topic = topic
        self._port = port
        self.verbose = verbose
        self.patient = patient
        self.timeout = timeout
        if client_id is None:
            client_id = 'mqtt_reader-' + random_choice()
        else:
            client_id += '-' + random_choice()

        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2,
                                   client_id=client_id)
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
        self.eot = False
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
            msg
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
            raise TimeoutError('Timeout while waiting for MQTT message.')
        payload, timestamp = remove_timestamp(msg.payload)
        if timestamp < 0:
            # End-of-transmission marker
            raise EOFError('End of transmission received.')
        return MqttMessage(msg.topic, payload, timestamp)

    def close(self):
        """Close the MQTT connection."""
        self._client.disconnect()
        self._client.loop_stop()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
