#!/usr/bin/env python

import sys
import shlex
import subprocess

def try_get_int_or_float(value):
    try:
        return int(value)
    except:
        return float(value)

class StatResult(object):
    def __init__(self, tag, num_points, num_polygons_limit, num_polygons_actual, score):
        self.tag = tag
        self.num_points = num_points
        self.num_polygons_limit = num_polygons_limit
        self.num_polygons_actual = num_polygons_actual
        self.score = score

    def __str__(self):
        return "tag={0} np={1} npl={2} npa={3} score={4}".format(
                    self.tag,
                    self.num_points,
                    self.num_polygons_limit,
                    self.num_polygons_actual,
                    self.score)

    @staticmethod
    def parse(tag, fd):
        values = []
        for line in fd.readlines():
            words = line.split("=")
            if len(words) == 1:
                values.append(try_get_int_or_float(words[0]))
            else:
                values.append(try_get_int_or_float(words[1]))
        values[0] = values[0] / 2
        return StatResult(tag, *values)


class Tester(object):
    CMD = "java -cp /home/stanislav/Documents/Code/Contests/TopCoder/Marathon/tco15/SmallPolygons/out/production/SmallPolygons/ SmallPolygons"
    def run(self, seed):
        cmd_line = "java -jar tester.jar -exec \"{0}\" -seed {1} -vis".format(self.CMD, seed)
        args = shlex.split(cmd_line)
        res = subprocess.Popen(args, stdout=subprocess.PIPE)
        res.wait()
        return StatResult.parse("simple", res.stdout)        


def main(seed):
    tester = Tester()
    if seed is not None:
        stat_res = tester.run(seed)
        print stat_res


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
