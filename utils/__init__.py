import logging

format_string = '%(asctime)s %(name)-15s [%(funcName)-45s]: %(levelname)-8s %(message)s'

console = logging.StreamHandler()
formatter = logging.Formatter(format_string)
console.setFormatter(formatter)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(console)