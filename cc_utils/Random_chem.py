
#not entirely sure what this does - splits a polymer into monomers? And then looks at the orbitals associated with the monomoers?
#this is from http://jarvistmoorefrost.wordpress.com/2012/01/30/donor-acceptor-character/

import re
from math import sqrt
from copy import deepcopy
import sys

#from http://physical-biochemistry.blogspot.co.uk/2013/09/non-interacting-particles.html
def Non_interacting_particles():
    #####  Non-Interacting Particles  #####
    #####        VERSION 1.3          #####
    #####       VIDEO VERSION         #####

    # Import of random, numpy and pylab library
    import pylab
    import numpy
    import copy
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    def save_video(X, Y, box_width, filename):
      """
        Takes 2D coordinates X[i] and Y[i], saves it into
        filename.mp4 with len(X) steps.
      """

      # Variables
      #box_width = 10.0
      video_frames = len(X)

      # Setup figure
      fig = plt.figure()
      fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
      ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                          xlim=(-(box_width+1.2), (box_width+1.2)), ylim=(-(box_width+0.4), box_width+0.4))

      # Particles holds the locations of the particles
      particles, = ax.plot([], [], 'ko', ms=6)

      # Set an Retangle around the box
      rect = plt.Rectangle((-box_width, -box_width), 2*box_width, 2*box_width, fill=False, linewidth=2.0, edgecolor='k')
      ax.add_patch(rect)

      # Remove default axis for cinematic effect
      plt.axis('off')

      def init():
          """initialize animation"""
          particles.set_data([], [])
          return particles

      def animate(i):
          """perform animation step"""
          particles.set_data(X[i], Y[i])
          return particles

      # Create the animation object
      ani = animation.FuncAnimation(fig, animate, frames=video_frames,
                                    interval=10, blit=True, init_func=init)

      # Save the animation to disk
      # Change fps for another framerate
      ani.save(filename+'.mp4', fps=25)

    ### initializing program ###

    # Definition of variables
    n_particles = 100
    n_steps = 10000
    time_difference = 0.001

    # initialize empty frame lists
    x_frames = []
    y_frames = []

    # Definition of position lists where each particle start in (0,0)
    position_x = [0 for i in range (n_particles)]
    position_y = [0 for i in range (n_particles)]

    # Definition of velocities (the change in movement to the time diffrence)
    velocity_x = [0 for i in range (n_particles)]
    velocity_y = [0 for i in range (n_particles)]
    speed =[0 for i in range (n_particles)]
    fixme_scale = 0.3 #FIXME scale.

    # The speed are given from the Maxwell-Boltzmann distribution
    for n in range (n_particles):
        speed_x = numpy.random.randn()
        speed_y = numpy.random.randn()
        speed_z = numpy.random.randn()
        speed[n]= fixme_scale * pylab.sqrt(speed_x*speed_x + speed_y*speed_y + speed_z*speed_z)
        # the speed should be Maxwell-Boltzmann distributed
        x = (2*numpy.random.rand() - 1) * speed[n]
        y = (2*numpy.random.rand() - 1) * speed[n]
        velocity_x[n] = x
        velocity_y[n] = y

    # Doing a histrogram for the distribution
    pylab.figure()
    n, bins, patches = pylab.hist(speed, 20, histtype='bar')
    pylab.xlabel("Speed")
    pylab.ylabel("Number of species")
    pylab.show()

    # Looping for n_steps
    for n in range (n_steps):
        # save a frame every 10 steps
        if n % 10 == 0:
            x_frames.append(copy.copy(position_x))
            y_frames.append(copy.copy(position_y))
        # Looping for n_particles
        for i in range (n_particles):
            # Doing the calculations and updating the position lists
            position_x [i] = position_x [i] + time_difference * velocity_x [i]
            position_y [i] = position_y [i] + time_difference * velocity_y [i]
            # Ensureing particles does not leave the box for the x and y position
            if abs(position_x [i])>1:
                velocity_x [i] = - velocity_x [i]
            if abs(position_y [i])>1:
                velocity_y [i] = - velocity_y [i]

    save_video(x_frames, y_frames, 1.0, "my_video3")

def GetMolCoords(filename):
    """	Read in the coordinates from an xyz file

    at the moment this is quite a lazy code that will fail if the lines start with an empty char.

    """
    filein = open(filename, 'r')
    res = []
    for line in filein:
        words = re.split(' ', line)
        coord = [float(words[1]), float(words[2]), float(words[3])]
        res.append(coord)
    return res

def GetDist(i,j):
    """ compute the distance between two 3d vectors ij represented as lists """
    res = 0.
    for k in range(3):
        res += (i[k]-j[k])*(i[k]-j[k])
    return sqrt(res)

def GetAllBonds(coords, dn =1.7):
    """ works out all the neighbours within a distance dn

    a pretty primitive n^2 scaling nearest neighbour search

    """
    edges = list()
    vertex = dict()
    n = len(coords)
    for i in range(n):
        vertex[i] = set()
        for j in range(i+1,n):
            if GetDist(coords[i],coords[j]) < dn:
                edges.append([i,j])
                vertex[i].add(j)
    return [edges,vertex]

def Partition(vertex):
    """ computes the partition class defined by the list of neighbour vertex[i] of points i

    I think I copied this from NR. The output is a list of indeces naming the partition that each vertex i belongs to.
    The input is a list of equivalence relation for each vertex i. Example:

    0_1_2  3__4
    vertex = [[1],[2],[],[4],[] ]
    class = [1,1,1,2,2]

    """
    #stick each atom in it's own class:
    classes = range(len(vertex))

    for v in vertex:
        ns = vertex[v] # the list of neighbours of v
        for n in ns:
            j = classes[v] # find the ancestor of v
            while classes[j] != j :
                j = classes[j]
            k = classes[n]   #find the ancestor of n
            while classes[k] != k :
                k = classes[k]
            if j!=k: #if they are not related, make them so
                classes[j]=k
    for i in range(len(vertex)) :
        while classes[i] != classes[classes[i]] :
            classes[i]=classes[classes[i]]
    return classes

def ConvertPartitionIntoDictionary(parts):
    """ takes the output from the Partition function and turns it into a dictionary
     grouping together  vertices from the same partition

    ex: part = [1,1,1,2,2]
    output={1:[0,1,2],2:[3,4]}
    """
    #convert parts into list of parts
    classes = dict()
    for i in range(len(parts)):
        key = parts[i]
        if  not key in classes:
            classes[key]= []
        classes[key].append(i)
    return classes

#eliminate a particular bond
#61 - 64
def EliminateABond(ind, vs,n, nmin=1):
    """ determines if the bond 'ind' is between different monomer units. vs is the current list
    of bonds and n is the current number of subparts of the molecule
            ___        ___
           /   \      /   \
        __/     \_XX_/     \___
          \     /    \     /
           \___/      \___/


    If bond XX is cut, then the molecule is split into two parts. Therefore we can
    define the monomer units as those places separated by a bond that - if cut-
    will result in a molecule with an increased number of molecules. Molecules
    are defined as partition classes of the equivalence relation "being bonded to".
    Of course, it is necessary to ignore those cases where the new partiion class
    created is of size nmin or less, for example C-H bonds if cut will creae a lone
    H atom.
    """
    [mini, maxi]= [min(ind), max(ind)]
    if not maxi in vs[mini]:
        print "bond dont exist"
        return
    atoms2 = deepcopy(vs)
    atoms2[mini].remove(maxi)
    partscheck = Partition(atoms2)
    # convert the partition information into a dictionary of which vertex belongs to which class
    classes = ConvertPartitionIntoDictionary(partscheck)

    # check if the cut was a small cut, i.e. if the size of the new partition is smaller <= nmin
    smallcut = False
    for pi, pn in classes.iteritems():
        if len(pn) <= nmin :
            smallcut = True
    #check if the # of partitions has changed
    lenafter = len(classes)
    if smallcut or lenafter==n :
        return [vs,n] # if not, return the same as before
    return [atoms2,n+1] # if yes, increment n and return the "cut" list of atoms

def GetHL(namein):
    """ read the natural transition orbitals from a log file

    """
    filein = open(namein, 'r')
    ntopercol = 5

    HOMO = []
    LUMO = []
    while filein:
        line = filein.readline()
        if re.match(' NBasis=', line):
            words = re.split('\W+', line)
            nbasis = int(words[2]) # confusing, the first word is an empty string cos the line starts with a separator

        if re.match(' Alpha spin Natural Transition Orbitals for state', line):
            #read the HOMO
            filein.readline()
            filein.readline()
            filein.readline()
            for i in range(nbasis):
                line = filein.readline()
                words = re.split('[ \t]+', line)
                if len(words) == ntopercol+5:
                    HOMO.append(0)
                HOMO[-1] += float(words[-1])*float(words[-1])

            #read the LUMO
            filein.readline()
            filein.readline()
            filein.readline()
            for i in range(nbasis):
                line = filein.readline()
                words = re.split('[ \t]+', line)
                if len(words) == ntopercol+5:
                    LUMO.append(0)
                LUMO[-1] += float(words[-ntopercol])*float(words[-ntopercol])
            break
    """
    print "#****ATOM LABEL | OCC NTO^2 | VIRT NTO^2****"
    for i in range(len(HOMO)):
        print "#", i+1, HOMO[i], LUMO[i]
    """
    return [HOMO, LUMO]

DEBUG = False
if __name__== '__main__'  and DEBUG == True:
    namegeom = sys.argv[1]
    namelog  = sys.argv[2]

    coords = GetMolCoords(namegeom)
    [NTOO,NTOV] = GetHL(namelog)
    totnto =0
    totntv =0
    for i in NTOO:
        totnto += i
    for i in NTOV:
        totntv += i

    [bonds, atoms] = GetAllBonds(coords)
    n=1
    [newat , newn ] =EliminateABond([62,65], atoms,n)
    print newn


elif __name__== '__main__' :
    namegeom = sys.argv[1]
    namelog  = sys.argv[2]

    coords = GetMolCoords(namegeom)
    [NTOO,NTOV] = GetHL(namelog)
    totnto =0
    totntv =0
    for i in NTOO:
        totnto += i
    for i in NTOV:
        totntv += i

    [bonds, atoms] = GetAllBonds(coords)
    n=1
    for b in bonds:
        [atoms,n] = EliminateABond(b, atoms,n)
    parts = Partition(atoms)
    """
    #this section prints out which atom is in which monomer unit
    print "#*** ATOM LABEL | PARTITION LABEL***"
    print "#n: ", n
    for i in range(len(parts)):
        print "#", i+1, parts[i]
    """
    #convert parts into list of parts
    classes = ConvertPartitionIntoDictionary(parts)

    print "#*************************************************"
    print "#*****************tot ntos per unit:**************"
    print "# PARTITION LABEL, OCC NTO%, VIRT NTO%, POS*******"
    for k, listat in classes.iteritems():
        toth =0.
        totl =0.
        totx =0.
        toty =0.
        totz =0.
        for at in listat:
            totx+=coords[at][0]
            toty+=coords[at][1]
            totz+=coords[at][2]
            toth+=NTOO[at]
            totl+=NTOV[at]
        print k, toth*100/totnto, totl*100/totntv, '<', totx/len(listat), ',' , toty/len(listat), ',' , totz/len(listat), '>'