import pprint
import scipy
import scipy.linalg   # SciPy Linear Algebra Library


A = []
with open("../test/matrice3.txt", "r") as matrice:

    for line in matrice:
        A.append(map(float, line.strip().split(" ")))


P, L, U = scipy.linalg.lu(A)

for line in U:
    print reduce(lambda x, y: x + " "+y,
             map(lambda x: "{"+str(x)+":.3f}",
                 range(0, len(line)))).format(*line)

print "\n"
for line in L:
    print reduce(lambda x, y: x + " "+y,
             map(lambda x: "{"+str(x)+":.3f}",
                 range(0, len(line)))).format(*line)


