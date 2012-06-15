#!/user/bin/env python
"""This is a testing script used to run semiconductor on multiple sized problems.
    Usage: size_test.py <id_of_test>"""

import ConfigParser
import subprocess
import time
import sys

program_location = "Debug/semiconductor"
config_location = "semiconductor.conf"

node = sys.argv[1]
config = ConfigParser.ConfigParser()
config.read("semiconductor.conf")
with open("ppn" + node + "-acml.txt", "w") as f:
    for x in [8]:
        for y in range(1, 12):
            num = str(y) + "00e-" + str(x)
            config.set("Problem", "length", num + " " + num)
            with open("semiconductor.conf", "wb") as c:
                config.write(c)
            start_time = time.time()
            out = subprocess.check_output([program_location, config_location])
            f.write(out)
            msg = "Run " + num + " took " + str(time.time() - start_time) + "seconds\n\n"
            f.write(msg)
            print msg
