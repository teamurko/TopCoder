#!/bin/sh
java -jar tester.jar -f test_data -test test_data/example_tests.txt -gtf test_data/example_gtf.txt -exec build/detector -seed 3 -vis output.png
