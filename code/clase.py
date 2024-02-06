import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk
import matplotlib.animation as animation
from IPython import display
from matplotlib.widgets import Slider, Button


class pendulo:
    global anim
    def __init__(self,fig):
        #Slider de la masa de la esfera 1
        axfreq = plt.axes([0.1, 0.18, 0.0225, 0.63])
        freq_slider = Slider(
        ax=axfreq,
        label='m1',
        valmin=0.0,
        valmax=10.0,
        valinit=5.0,
        orientation="vertical"
        )
        #al movernos en el slider actualizamos con la funcion update
        freq_slider.on_changed(self.updatem1)

        #Slider de la esfera de la masa 2
        axfreq2 = plt.axes([0.05, 0.18, 0.0225, 0.63])
        freq_slider2 = Slider(
        ax=axfreq2,
        label='m2',
        valmin=0.0,
        valmax=10.0,
        valinit=5.0,
        orientation="vertical"
        )
        #al movernos en el slider actualizamos con la funcion update
        freq_slider2.on_changed(self.updatem2)

        #Slider del largo de la primera cuerda
        axfreq3 = plt.axes([0.85, 0.18, 0.0225, 0.63])
        freq_slider3 = Slider(
        ax=axfreq3,
        label='l1',
        valmin=0.0,
        valmax=20.0,
        valinit=10.0,
        orientation="vertical"
        )
        #al movernos en el slider actualizamos con la funcion update
        freq_slider3.on_changed(self.updatel1)

        #Slider del largo de la segunda cuerda
        axfreq4 = plt.axes([0.9, 0.18, 0.0225, 0.63])
        freq_slider4 = Slider(
        ax=axfreq4,
        label='l2',
        valmin=0.0,
        valmax=20.0,
        valinit=10.0,
        orientation="vertical"
        )
        #al movernos en el slider actualizamos con la funcion update
        freq_slider4.on_changed(self.updatel2)

        #Slider del valor de theta1
        axfreq5 = plt.axes([0.2, 0.95, 0.65, 0.03])
        freq_slider5 = Slider(
        ax=axfreq5,
        label='theta1',
        valmin=-90.0,
        valmax=90.0,
        valinit=0.0,
        )
        #al movernos en el slider actualizamos con la funcion update
        freq_slider5.on_changed(self.updatetheta1)

        #Slider del valor de theta2
        axfreq6 = plt.axes([0.2, 0.9, 0.65, 0.03])
        freq_slider6 = Slider(
        ax=axfreq6,
        label='theta2',
        valmin=-90.0,
        valmax=90.0,
        valinit=0.0,
        )
        #al movernos en el slider actualizamos con la funcion update
        freq_slider6.on_changed(self.updatetheta2)

        self.x1_data = list()
        self.y1_data = list()
        self.x2_data = list()
        self.y2_data = list()
        ax1 = fig.add_axes([0.835, 0.02, 0.1, 0.075])
        startButton = Button(ax1, 'Start')
        startButton.on_clicked(self.begin)
        self.x0 = 0
        self.y0 = 0
        self.x1, self.y1 = 0, 0
        self.x2, self.y2 = 0, 0
        self.g = 9.81
        self.m1 = 5
        self.m2 = 10
        self.l1 = 15
        self.l2 = 10
        self.theta1 = 80*np.pi/180
        self.theta2 = 20*np.pi/180     
        self.trayectory1, = ax.plot([0,0], [0,0], lw=1, c='b')
        self.trayectory2, = ax.plot([0,0], [0,0], lw=1, c='r')   
        self.m = self.m1/self.m2
        self.r1 = self.radius(self.m1)
        self.r2 = self.radius(self.m2)
        self.aspect = 1000
        self.line1, = ax.plot([self.x0, self.x1], [self.y0, self.y1], lw=3, c='k')
        self.line2, = ax.plot([self.x1, self.x2], [self.y1, self.y2], lw=3, c='k')
        self.circle1 = ax.add_patch(plt.Circle((self.x1, self.y1), 0, fc='b', zorder=3))
        self.circle2 = ax.add_patch(plt.Circle((self.x2, self.y2), 0, fc='r', zorder=3))
        plt.show()    

    def updatem1(self, t):
        self.m1 = t
        self.m = self.m1/self.m2
        self.r1 = self.radius(self.m1)

    def updatem2(self, t):
        self.m2 = t
        self.m = self.m1/self.m2
        self.r2 = self.radius(self.m2)   

    def updatel1(self, t):
        self.l1 = t

    def updatel2(self, t):
        self.l2 = t

    def updatetheta1(self, t):
        self.theta1 = t*np.pi/180

    def updatetheta2(self, t):
        self.theta2 = t*np.pi/180 


    def changeAngle(self, theta):
        return theta-0.5*np.pi

    def convertAngle(self, theta):
        return theta-0.5*np.pi

    def firstPendulumPos(self,t1):  
        return self.x0 + np.cos(self.changeAngle(t1))*self.l1, self.y0 + np.sin(self.changeAngle(t1))*self.l1

    def secondPendulumPos(self, a1, b1, t2):
        return a1 + np.cos(self.changeAngle(t2))*self.l2, b1 + np.sin(self.changeAngle(t2))*self.l2

    def lerp(x0, y0, x1, y1):
        return np.linspace(x0, x1, 1000), np.linspace(y0, y1, 1000)

    def updateDoublePendulumPos(self, a1, b1, a2, b2, t1, t2):
        a1,b1 = self.firstPendulumPos(t1)
        a2,b2 = self.secondPendulumPos(a1, b1, t2)
        return a1, b1, a2, b2

    def saveDoublePendulumPos(self, a1, b1, a2, b2):
        self.x1_data.append(a1)
        self.y1_data.append(b1)
        self.x2_data.append(a2)
        self.y2_data.append(b2)
    def radius(self, m):
        return 0.3*m    
    def animate(self, i):
        minimumi = 0
        if i > 200:
            minimumi = i-200
        x1, y1 = self.x1_data[i], self.y1_data[i]
        x2, y2 = self.x2_data[i], self.y2_data[i]
        self.line1.set_data([self.x0, x1], [self.y0, y1])
        self.line2.set_data([x1, x2], [y1, y2])
        self.circle1.set_center((x1, y1))
        self.circle2.set_center((x2, y2))
        self.trayectory1.set_xdata(self.trayectory1_x[minimumi:i])
        self.trayectory1.set_ydata(self.trayectory1_y[minimumi:i])
        self.trayectory2.set_xdata(self.trayectory2_x[minimumi:i])
        self.trayectory2.set_ydata(self.trayectory2_y[minimumi:i])

        if i == self.nsteps-1:
            self.trayectory1.set_xdata(self.trayectory1_x[i-1:i])
            self.trayectory1.set_ydata(self.trayectory1_y[i-1:i])
            self.trayectory2.set_xdata(self.trayectory2_x[i-1:i])
            self.trayectory2.set_ydata(self.trayectory2_y[i-1:i])
            self.line1.set_data([self.x0,self.x1_data[0]], [self.y0, self.y1_data[0]])
            self.line2.set_data([self.x1_data[0], self.x2_data[0]], [self.y1_data[0], self.y2_data[0]]) 
            self.circle1.set_center((self.x1_data[0], self.y1_data[0]))
            self.circle2.set_center((self.x2_data[0], self.y2_data[0]))  
            self.x1, self.y1 = self.x1_data[0], self.y1_data[0]
            self.x2, self.y2 = self.x2_data[0], self.y2_data[0]

    def begin(self, event):
        self.x1, self.y1 = 0, 0
        self.x2, self.y2 = 0, 0
        self.circle1.remove()
        self.circle1 = ax.add_patch(plt.Circle((self.x1, self.y1), self.r1, fc='b', zorder=3))
        self.circle2.remove()
        self.circle2 = ax.add_patch(plt.Circle((self.x2, self.y2), self.r2, fc='r', zorder=3))
        self.x1_data.clear()
        self.y1_data.clear()
        self.x2_data.clear()
        self.y2_data.clear()

        theta1_dot = 0
        theta2_dot = 0

        theta1 = self.theta1
        theta2 = self.theta2
        t = 0
        dt = 1/60

        while t<20:
            a1 = (self.m+1)*self.l1/self.l2
            a2 = self.l2/self.l1
            b = np.cos(theta1-theta2)
            c1 = (theta2_dot**2)*np.sin(theta1-theta2)+(self.m+1)*self.g*np.sin(theta1)/self.l2
            c2 = -(theta1_dot**2)*np.sin(theta1-theta2)+self.g*np.sin(theta2)/self.l1
            
            theta1_ddot = (a2*c1 - b*c2)/(b**2 - a1*a2)
            theta2_ddot = (a1*c2 - b*c1)/(b**2 - a1*a2)
            
            theta1_dot += theta1_ddot*dt
            theta2_dot += theta2_ddot*dt

            theta1 += theta1_dot*dt
            theta2 += theta2_dot*dt
            
            self.x1,self.y1,self.x2,self.y2 = self.updateDoublePendulumPos(self.x1, 
            self.y1, self.x2, self.y2, theta1, theta2)
            self.saveDoublePendulumPos(self.x1,self.y1,self.x2,self.y2)

            t += dt
        self.nsteps = len(self.x1_data)
        self.trayectory1_x = [[self.x1_data[i],self.x1_data[i+1]] for i in range(self.nsteps-1)]
        self.trayectory1_y = [[self.y1_data[i],self.y1_data[i+1]] for i in range(self.nsteps-1)]
        self.trayectory2_x = [[self.x2_data[i],self.x2_data[i+1]] for i in range(self.nsteps-1)]
        self.trayectory2_y = [[self.y2_data[i],self.y2_data[i+1]] for i in range(self.nsteps-1)]

        self.circle1.set_center((self.x1_data[0], self.y1_data[0]))
        self.circle2.set_center((self.x2_data[0], self.y2_data[0]))  

        nframes = self.nsteps
        interval = 1
        global anim
        anim = animation.FuncAnimation(fig, self.animate, frames=nframes, 
        repeat=False, interval=interval, blit = False)          



# Initialize the animation plot. Make the aspect ratio equal so it looks right.
fig = plt.figure()
ax = fig.add_subplot(aspect='equal')  
# Set the plot limits so that the pendulum has room to swing!
ax.set_xlim( -1.25*(40) , 1.25*(40) )
ax.set_ylim( -1.25*(40) , 1.25*(40) )

pendulo(fig)