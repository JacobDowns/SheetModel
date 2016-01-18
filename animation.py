import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Spacial coordinate
xs = np.linspace(0, 60.0, 3)
# Seconds per year
spy = 60.0 * 60.0 * 24.0 * 365.0 
# End time
T = 1.5 * spy
# Time step
dt = 60.0 * 60.0 * 8.0
#ts = np.linspace(0, 1.5 * spy, 200)
ts = np.arange(0, 1.5 * spy, dt)

print len(ts)

def update_melt(num, xs, ts, l):
  #Get current time
  t = ts[num]
  # Compute melt pattern 
  ms = 0.5 * (3.0 * (60.0 - xs) / 60.0) * (np.cos( ((2.0 * np.pi) / spy) * t) + 1.0)
  # Update plot
  l.set_data(xs, ms)
  return l,

# Set up plot and animation
fig1 = plt.figure(num = None, figsize=(10, 5))


l, = plt.plot([], [], 'k-')
plt.xlim(0, 60.0)
plt.ylim(0.0, 3.1)
plt.xlabel('X (km)')
plt.ylabel('Melt (m/a)')
plt.title('Melt Rate')
line_ani = animation.FuncAnimation(fig1, update_melt, len(ts), fargs=(xs, ts, l),
    interval = 66.6, blit = False)
line_ani.save('melt.mp4')

plt.show()

#animation.FuncAnimation()



