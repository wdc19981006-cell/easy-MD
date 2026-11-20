"""Run upstream CLI with the Python 3.12 ConfigParser API rename shim.

readfp was the old name of read_file. This changes configuration I/O only;
the LightDock GSO, DFIRE and candidate-generation implementations are untouched.
"""
import configparser
from pathlib import Path
import runpy
import sys

if not hasattr(configparser.ConfigParser, 'readfp'):
    configparser.ConfigParser.readfp = configparser.ConfigParser.read_file

if __name__ == '__main__':
    script = str(Path(sys.argv[1]).resolve())
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(script, run_name='__main__')
