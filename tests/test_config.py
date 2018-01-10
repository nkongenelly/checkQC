
import unittest

from checkQC.config import Config, ConfigFactory


class TestConfig(unittest.TestCase):

    def setUp(self):
        self.first_handler = {"name": "first_handler"}
        self.second_handler = {"name": "second_handler"}
        self.default_handler = {"name": "default_handler", "config": "some_value"}
        config_dict = {'instrument_type_mappings': {'SN': 'hiseq2000', 'M': 'miseq', 'D': 'hiseq2500', 'ST': 'hiseqx'},
                       'miseq_v3': {
                           300: {'handlers': [
                               self.first_handler]},
                           '150-299': {'handlers': [
                               self.second_handler]}
                       },
                       "default_handlers": [
                           self.default_handler,
                           self.first_handler
                       ]}
        self.config = Config(config_dict)

    def test_exact_match(self):
        handlers = self.config.get_handler_configs('miseq_v3', 300)
        self.assertListEqual(handlers, [self.first_handler, self.default_handler])

    def test_interval_match(self):
        handlers = self.config.get_handler_configs('miseq_v3', 175)
        self.assertListEqual(handlers, [self.second_handler, self.default_handler, self.first_handler])

    def test_no_match(self):
        with self.assertRaises(KeyError):
            self.config.get_handler_configs('miseq_v3', 999)

    def test_call_with_str(self):
        handlers = self.config.get_handler_configs('miseq_v3', "300")
        self.assertListEqual(handlers, [self.first_handler, self.default_handler])


class TestConfigFactory(unittest.TestCase):

    def test_get_logging_config_file_default(self):
        result = ConfigFactory.get_logging_config_dict(None)
        default_config = {'version': 1, 'disable_existing_loggers': False,
                          'formatters':
                              {'simple': {'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s'}},
                          'handlers':
                              {'console':
                                   {'class': 'logging.StreamHandler', 'level': 'DEBUG', 'formatter': 'simple',
                                    'stream': 'ext://sys.stdout'},
                               'file_handler': {'class': 'logging.handlers.RotatingFileHandler', 'level': 'DEBUG',
                                                'formatter': 'simple', 'filename': 'checkqc-ws.log',
                                                'maxBytes': 10485760, 'backupCount': 20,
                                                'encoding': 'utf8'}},
                          'root':
                              {'level': 'DEBUG', 'handlers': ['console', 'file_handler']}}
        self.assertEqual(result, default_config)


if __name__ == '__main__':
    unittest.main()
