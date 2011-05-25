
import numpy as np
import pyMPCD
from matplotlib.pyplot import subplot, plot, hist, xlabel, ylabel, legend, show

# take a testcase from the library
box = pyMPCD.test_cases.eight_PBC()

# run for 10 times 10 steps
kin_data = []
for i in range(10):
    for t in range(10):
        box.one_full_step()
        # replacing the above line by box.one_full_step()
        # makes use of the Fortran version instead.
        kin_data.append(np.sum(box.so_v**2))
    print "Total kinetic energy   : ", np.sum(box.so_v**2)
    print "Total momentum         : ", np.sum(box.so_v, axis=0)

kin_data = np.array(kin_data)
kin_0 = kin_data[0]

subplot(211)
# The kinetic energy in the course time, measured in MPCD steps
plot(kin_data-kin_0)
xlabel('MPCD steps')
ylabel('kinetic energy deviation')
subplot(212)
# The velocity distribution along the three axis
hist(box.so_v, normed=True, bins=20, label=['x-axis','y-axis','z-axis'])
xlabel('velocity along the three axis')
ylabel('velocity')
legend()
show()
