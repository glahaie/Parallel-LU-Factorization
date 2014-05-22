#! /usr/bin/python
# -*- coding:utf8 -*-
import pprint


def gauss(A, n):

    for i in range(0, n):
        for j in range(i+1, n):
            A[j][i] = A[j][i] / A[i][i]
            for k in range(i+1, n):
                A[j][k] = A[j][k] - A[j][i]*A[i][k]



def main():

    A=[[5.00, -9.00, 3.00, 14.00, 22.00, 23.00],
       [15.00, 17.00, -29.00, 23.00, 2.00, 1.00],
       [44.00, 39.00, 51.00, 13.00, 12.00, 11.00],
       [59.00, -8.00, 14.00, -17.00, 3.00, 41.00],
       [21.00, -27.00, 64.00, 44.00, 13.00, 19.00],
       [47.00, 52.00, -32.00, -5.00, 6.00, 12.00]]
#    A = []
    #with open("matrice4.txt", "r") as matrice:
        #for line in matrice:
            #A.append(map(float, line.strip().split(" ")))

    gauss(A, 6)

    for line in A:
        print reduce(lambda x, y: x + " "+y,
                map(lambda x: "{"+str(x)+":.3f}",
                    range(0, len(line)))).format(*line)

    #for i in b:
        #print i

    #gauss(b, 4)

    #for i in b:
        #print i

if __name__ == "__main__":
    main()
