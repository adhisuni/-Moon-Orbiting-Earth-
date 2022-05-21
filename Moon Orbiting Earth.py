#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import vpython as vp
M_E = 5.97*10**24

M_M = 7.35*10**22
Radius_E = 6.4*10**6
Radius_M = 1.35*10**6

R_apo = 4.054 * 10**8
R_peri = 3.64*10**8


def angVel( time ):
    return 2 * vp.pi / time
from vpython import *
from scipy.constants import G, pi
scene = canvas()

### Condition given by question:
t = 0
dt = 0.5
T = 27.3               ## average period of moon to orbit the earth

X_E = vp.vec(0,0,0)    #### this is the initial displacement of the earth
X_M = vp.vec(0,0,R_apo)

theta = angVel( 86400 ) * dt                           
phi = angVel( 2360534 ) * dt  
### We need to think about the vector and dimensions while creating this simulation. So, for acceleration or velocity we need to focus on vectors
v_e = vp.vec(0,0,0) # this is the earth's initial velocity when the earth is not moving at all
v_m = vp.vec(962,0,0) # the moon is travelling towards the positive x-axis initially

F = - G * M_E * M_M / mag2(X_M)  #### this is force of gravitation formula which will further assis us to find 'a'


a_E = vp.vec(0,0,0) ### no acceleration for earth
a_M = vp.vec(F/M_M,0,0)

L = vp.vec( 0, 0, M_M *(X_M.z * v_m.x - X_M.x * v_m.z) ) #### this is the linear momentum of the Moon

### Now for the energies:
KE = 0.5 * M_M * mag2(v_m) ####since moon is the only one in motion, we have its KE but not that of earth
PE = - G * M_E *M_M / R_apo  ### using formula kind of similar to F we had earlier
TE = KE + PE

### starting woking on vpython:

Earth = sphere(pos=X_E,
              vel = v_e,
              a = a_E,
              texture = textures.earth,
              radius = Radius_E)
#### doing exactly same for the moon

Moon = sphere(pos=X_M,
             vel = v_m,
             a = a_M,
             radius=Radius_M,
             texture=textures.metal,
             make_trail=True,
             trail_radius=0.55* Radius_M,
             retain=10000,
             interval=100)

Apogee_Mark = arrow( pos=X_M + vp.vec( 0, 6 * Radius_E, 0 ),
                     axis=vp.vec( 0, - 4 * Radius_E, 0 ),
                     shaftwidth=1.5 * Radius_M,
                     texture=textures.metal )


Perigee_Mark = arrow( pos=vp.vec( 0, 10 * Radius_E, -360254141.7010045 ), 
                      axis=vp.vec( 0, - 8 * Radius_E, 0 ),
                      shaftwidth=3.5 * Radius_M,
                      texture=textures.metal )

#### IN order to see the trail, movement of the moon, apogee and perigee point, I have to make the arial view.
scene.camera.pos = Earth.pos - ( ( Moon.pos - Earth.pos ) / 15 ) + vec( 0, 1.4 * Radius_E, 0 )
scene.camera.axis = Moon.pos - scene.camera.pos

scene.camera.pos = Moon.pos + ( Moon.pos / 50 ) + vec( 0, 1.4 * Radius_M, 0 )
scene.camera.axis = Earth.pos - scene.camera.pos

scene.camera.pos = vp.vec( 1.644e8, 1.37829e8, 5.94024e8 )
scene.camera.axis = vp.vec( -1.91006e8, -1.60134e8, -6.90157e8 )

####Now to plot with the given condition:
plot_now = vp.graph( x=0, y=0, width=500, height = 500,
                xmin=-4.5*10**8, xmam=4.5*10**8,
                ymin=-4.5*10**8, ymax=4.5*10**8,
                foreground=color.black,
                background=color.white,
                title='X vs. Y Position',
                xtitle='x(t) [m]',
                ytitile='y(t) [m]')

f1 = vp.gcurve(color=color.red)
Energy_plot = graph( x=0, y =0, width =500, height=500,
                   xmin=0, xmax=27.3,
                   ymin=-1*10**29, ymax = 1*10**29,
                   foreground=color.green,
                   background=color.white,
                   xtitle='time [days]',
                   ytitile='energy [J]')

KEC = vp.gcurve(color=color.black)
PEC = vp.gcurve(color=color.blue)
TEC = vp.gcurve(color=color.cyan)

T_plot = graph(x=0, y=0, width=500, height = 500,
              xmin=0, xmax=27.3,
              ymin=1*10**34, ymax=1*10**35,
              foreground=color.blue,
              background=color.white,
              title='Angular Momentum Vs. Time',
              xtitle='time[days]',
              ytitle='L[J S]')
T1= gcurve(color=color.red)


while True:
    vp.rate(900000)
    
    ### gravitational force and angular momentum updates first (especially force)
    FG = - G * M_E * M_M / mag2( Moon.pos )
    L = M_M * vec( 0, 0, ( Moon.pos.z * Moon.vel.x ) - ( Moon.pos.x * Moon.vel.z ) )

   
   
    Moon.trail_color = vp.vec( abs( Moon.pos.x ), abs( Moon.pos.z ), abs( Moon.pos.x - Moon.pos.z ) ) / R_apo
    
    ### Euler-Cromer updates for the Moon ###
    Moon.acc = ( FG / M_M / mag( Moon.pos ) ) * vp.vec( Moon.pos.x, 0, Moon.pos.z )
    Moon.vel = Moon.vel + Moon.acc * dt
    Moon_past_pos = Moon.pos
    Moon.pos = Moon.pos + Moon.vel * dt
    
    ### energies of the system ###
    KE = 0.5 * M_M * mag( Moon.vel )**2
    PE = - G * M_E * M_M / mag( Moon.pos )
    TE = KE + PE
   
            
    ### plots ###
    f1.plot( Moon.pos.x, Moon.pos.z )
    KEC.plot( t, KE )
    PEC.plot( t, PE )
    TEC.plot( t, TE )
    T1.plot( t, L.z )
    
    t =  t+ dt
  


# In[ ]:




