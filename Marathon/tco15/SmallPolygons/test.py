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
            print line
            words = line.split("=")
            if len(words) == 1:
                values.append(try_get_int_or_float(words[0]))
            else:
                values.append(try_get_int_or_float(words[1]))
        values[0] = values[0] / 2
        return StatResult(tag, *values)


class Tester(object):
    CMD_TEMPLATE = "java -cp /home/stanislav/Documents/Code/Contests/TopCoder/Marathon/tco15/SmallPolygons/out/production/SmallPolygons/ SmallPolygons {0}"
    def __init__(self, vis=False, tag="baseline", debug=False):
        self.vis = vis
        self.tag = tag
        self.debug = debug
        self.cmd = self.CMD_TEMPLATE.format(tag)

    def run(self, seed):
        cmd_line = "java -jar tester.jar -exec \"{0}\" -seed {1} {2} {3}".format(
                    self.cmd,
                    seed,
                    "-vis" if self.vis else "",
                    "-manual -debug" if self.debug else "")
        print cmd_line
        args = shlex.split(cmd_line)
        res = subprocess.Popen(args, stdout=subprocess.PIPE)
        res.wait()
        print res.stdout
        return StatResult.parse(self.tag, res.stdout)        


def main(tag, debug, seed):
    if seed is not None:
        tester = Tester(vis=True, tag=tag, debug=debug)
        stat_res = tester.run(seed)
        print stat_res
    else:
        num_tests = 100
        score = 0.0
        tester = Tester(tag=tag)
        for seed in range(1, num_tests + 1):
            stat_res = tester.run(seed)
            print "Test {0}: {1}".format(seed, stat_res)
            score += stat_res.score
        print score / num_tests

if __name__ == "__main__":
    main(
        sys.argv[1] if len(sys.argv) > 1 else "baseline",
        "debug" == sys.argv[2] if len(sys.argv) > 2 else False,
        int(sys.argv[3]) if len(sys.argv) > 3 else None)
