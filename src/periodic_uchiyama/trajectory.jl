import time

class Trajectory:
    """ Classe contenant les donnees """
    def __init__(self, npart, epsilon,iseed=1):

        """ Initialize particles and build collision matrices """
        self.npart = npart
        self.epsilon = epsilon
        self.q, self.v = init_particles(npart, epsilon,iseed)
        self.collisions = PCollisionMatrix(npart, self.q, self.v, epsilon)
        self.marker = 0
        self.part_a = []
        self.part_b = []
        self.absc=[0.5]
        self.ord=[0.5]
        self.absc_ind=[0.5]
        self.ord_ind=[0.5]
        self.v_x=[self.v[0][0]]
        self.v_y=[self.v[0][1]]
        self.v_x_ind=[self.v[0][0]]
        self.v_y_ind=[self.v[0][1]]
        self.tcoll  = 0.0
        self.ctime  = []
        self.etime  = 0.0



function step(self)

    if self.tcoll == 0:
        self.etime = 0.0

    start_time = time.time()

    """ Compute small time step between two collisions """

    i1, i2 = self.collisions.dt_min_position

    tempsp = self.collisions.dt_min
    num_fant = self.collisions.fant-1

    self.q = (self.q + tempsp * self.v)%1

    """ Algorithme de collision inter particules"""

    if self.v[i1] * self.v[i2] == 0

        self.v[[i1, i2]] = self.v[[i2, i1]]  # swap velocities
        if i1 == self.marker:
            self.marker = i2
        elseif i2 == self.marker:
            self.marker = i1

    elseif dot(self.v[i1], self.v[i2]) == -1:

        rot = array([[0, -1], [1, 0]])

        if (rot * self.v[i1]) * (self.q[i2] + offset[num_fant] - self.q[i1]) < 0
            self.v[i1] = rot * self.v[i1]
            self.v[i2] = rot * self.v[i2]
        elseif (rot * self.v[i1]) * (self.q[i2] + offset[num_fant] - self.q[i1]) > 0
            self.v[i1] = - rot * self.v[i1]
            self.v[i2] = - rot * self.v[i2]
        end

    self.collisions.dt = self.collisions.dt - tempsp
    self.collisions.reset(i1)
    self.collisions.reset(i2)


    compute_dt(i1, self.q, self.v, self.epsilon, self.collisions.dt, self.collisions.fantome)

    compute_dt(i2, self.q, self.v, self.epsilon, self.collisions.dt, self.collisions.fantome)


    self.part_a.append(i1+1)
    self.part_b.append(i2+1)
    self.absc.append(self.q[0][0])
    self.ord.append(self.q[0][1])
    self.absc_ind.append(self.q[self.marker][0])
    self.ord_ind.append(self.q[self.marker][1])
    self.v_x.append(self.v[0][0])
    self.v_y.append(self.v[0][1])
    self.v_x_ind.append(self.v[self.marker][0])
    self.v_y_ind.append(self.v[self.marker][1])


    self.tcoll += tempsp

    self.ctime.append(self.tcoll)

    end_time = time.time()

    self.etime += end_time - start_time

end

function run(self, nstep=15000, ctime=2)

    istep = 0
    while istep < nstep && self.tcoll < ctime

        istep += 1
        step(self)

        return istep, self.tcoll, self.etime

    end

end
