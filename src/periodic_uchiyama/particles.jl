# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:09:32 2016

@author: nayi
"""

from numpy import zeros, array, dot, sqrt, matrix, floor
from numpy.random import rand, seed
from numpy.linalg import norm
import matplotlib.pyplot as plt


def init_particles(n, epsilon,iseed=1):
    """ Tirages des positions initiales des particules sans overlap """
    q = zeros((n, 2))
    l1 = 0.5
    l2 = 0.5
    q[0]=array([0.5,0.5])

    J4 = matrix([[1, 0], [0, 1]])

    seed(iseed)
    # tirage carre avec départ au centre pour la première particule
    # for klm = 2:N
    for klm in range(1, n):
        frein = 0
        overlap = 1
        while (overlap == 1) and frein < 10000:

            p0 = array([[2 * epsilon + (1 - 4 * epsilon) * rand(1), 2 * epsilon + (1 - 4 * epsilon) * rand(1)]])
            p = dot(J4, p0)
            overlap2 = 0
            # for subh in 1:(klm-1)
            for subh in range(0, klm):
                if norm(p - q[subh], 1) < 2 * epsilon:
                    overlap2 = 1

            overlap = overlap2
            frein += 1
        if frein == 10000:
            print("Echec tirage initial")
        q[klm] = p

    """ Tirage uniforme sur les vitesses """
    vitesses = array([[1, 0], [0, 1], [-1, 0], [0, -1]])
    Ia = floor(4 * rand(n, 1)).astype(int)
    v = zeros((n, 2))
    for j in range(n):
        v[j] = vitesses[Ia[j]]

#    """ Tirage non uniforme sur les vitesses """
#    vitesses = array([[1, 0], [0, 1], [-1, 0], [0, -1]])
#    v = zeros((n, 2))
#    for j in range(n):
#        r = rand();
#        if r < 0.5:
#           v[j] = vitesses[0]
#        elif  r < 0.75 and r > 0.5:
#           v[j] = vitesses[2]
#        elif r < 0.90 and r > 0.75:
#           v[j] = vitesses[1]
#        elif r > 0.9:
#           v[j] = vitesses[3]

    return q, v


class Lozenge:
    def __init__(self, center, h=1, fc='b'):
        self.center = center
        self.h = h
        x, y = center
        xy = [[x + h, y], [x, y + h], [x - h, y], [x, y - h]]
        self.patch = plt.Polygon(xy, closed=True, edgecolor=None, fc=fc)

    def set_center(self, center):
        x, y = center
        self.center = center
        xy = [[x + self.h, y], [x, y + self.h], [x - self.h, y], [x, y - self.h]]
        self.patch.set_xy(xy)
