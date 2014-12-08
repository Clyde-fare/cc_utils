import numpy as np
import math

def unwind(lst_of_lsts):
    return [e for l in lst_of_lsts for e in l]

def norm(vector):
    """ Returns the norm (length) of the vector."""
    # note: this is a very hot function, hence the odd optimization
    # Unoptimized it is: return np.sqrt(np.sum(np.square(vector)))
    return np.sqrt(np.dot(vector, vector))

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if math.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle


##analysis of Lee's mode character
# etd = np.array([-0.5082, 1.6651, 2.1033])
# with open('DDs.txt') as data_f:
#     data = data_f.readlines()
# proc_data = [l.split(':')[1] for l in data]
# proc_data = [l.strip() for l in proc_data]
# proc_data_2 = [l.replace('D', 'E') for l in proc_data]
# proc_data_3 = [l.split() for l in proc_data_2]
#
# for l in proc_data_3:
#     for i in range(len(l)):
#         l[i] = float(l[i])
#
# proc_data_4 = [np.array(l) for l in proc_data_3]
# dds = proc_data_4
#
# angles = [angle_between(etd, d) for d in dds]
# angles = [abs(a-(math.pi/2)) for a in angles]
# freqs=np.loadtxt('freqs.txt')
# freqs = [f[1] for f in freqs]
#
#
# def custavg(angles, freqs, cutoff):
#     """average of angles based on frequency"""
#     masks = []
#     for i in range(len(freqs)):
#         mask = []
#         backward_range = range(len(freqs))[:i]
#         backward_range.reverse()
#
#         for j in backward_range:
#             if abs(freqs[j]-freqs[i]) < cutoff:
#                 mask.append(j)
#             else:
#                 break
#         mask.reverse()
#
#         forward_range = range(len(freqs))[i:]
#         for j in forward_range:
#             if abs(freqs[j]-freqs[i]) <= cutoff:
#                 mask.append(j)
#             else:
#                 break
#
#         masks.append(mask)
#     avged_angles = [mean([angles[i] for i in mask]) for mask in masks]
#     return avged_angles

# %pylab
# n=10
#
# #plot(movavg(freqs,n), movavg(angles,n), 'ro')
# plot(freqs,angles,'ro')
# xlabel('frequency')
# ylabel('mode angle')
# show()
#
# import matplotlib.animation as animation
# fig = plt.figure()
# #ax = plt.axes(xlim=(0, 4000), ylim=(0, 3))
# ax = plt.axes(xlim=(1250, 1800), ylim=(0, math.pi/2))
# plt.xlabel('frequency')
# plt.ylabel('mode angle')
# line, = ax.plot([], [], lw=2)
#
# #range 1200-1850 = [7552:10556]
# scale_factor = 0.96
# freqs = orig_freqs[7552:10556]
# freqs = [f*scale_factor for f in freqs]
# angles = orig_angles[7552:10556]
#
# max_freq_window = 2
# frames = 100
# freq_windows = [i * float(max_freq_window)/frames for i in range(frames)]
# angle_data = [custavg(angles, freqs, w) for w in freq_windows]

# # initialization function: plot the background of each frame
# def init():
#     line.set_data([], [])
#     return line,
#
# #update frame
# def animate(i):
# #    x = movavg(freqs,i+1)
# #    y = movavg(angles,i+1)
#     x = freqs
#     y = angle_data[i]
#     line.set_data(x, y)
#     return line,
#
#
#
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=frames, interval=10, blit=True)
# anim.save('movavg_Data_v3.mp4')