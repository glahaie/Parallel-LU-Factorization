import numpy

A = 100* numpy.random.randn(30, 30)

for line in A:
    print reduce(lambda x, y: x + " "+y,
             map(lambda x: "{"+str(x)+":.3f}",
                 range(0, len(line)))).format(*line)
