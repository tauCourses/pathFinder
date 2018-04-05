from subprocess import Popen, PIPE
from os import listdir
from os.path import isfile, join
from filecmp import cmp
import time

dir = 'tests'
print("run tests")
robots = [join(dir, f) for f in listdir(dir) if isfile(join(dir, f)) and join(dir, f).startswith(dir + '/robot')]
obstacles = [join(dir, f) for f in listdir(dir) if isfile(join(dir, f)) and join(dir, f).startswith(dir + '/obstacles')]

robots.sort()
obstacles.sort()

for test in zip(robots, obstacles):
    p = Popen(['./PathFinder', test[0], test[1], 'output1'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    time.sleep(0.1)

    print('test %s' % str(test))
    if p.returncode != 0:
        print('FAILED %s' % str(err.decode()))
        print('error code %d' % p.returncode)
        continue

    print('PASSED')
    p = Popen(['python3.6', 'PreviewPy.py', test[0], test[1], 'output1'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    time.sleep(1)
