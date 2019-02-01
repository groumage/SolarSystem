from matplotlib import pyplot as plt
from numpy import *
from mpl_toolkits.mplot3d import Axes3D

masse = array([1.00000597682, 0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 1/(1.3e8)])

pos_ini = array([[0, 0, 0],[2.2911898857058, 4.4742096428394, -0.0698526085811],[-8.4077274422046, -4.9386528867475, 0.4206194745212],[19.9506411858858, -2.1189594536568, -0.2505524004432],[26.3841113101207, -14.2637141476417, -0.3141826415529],[4.7240342661351, -31.8776619468451, 2.0454833905119]])

vit_ini = array([[0, 0, 0],[-0.0068140485228, 0.0038022184356, 0.0001366932116],[0.0025181864149, -0.0048203311629, -0.000016441567],[-0.0004495229028, 0.0037304364206, 0.0000197241717],[0.0014665671695, 0.0027825194375, -0.0000910521187],[0.0031588952992, -0.0001620298531, -0.0008976519916]])

g = 2.95912208286e-4

planet = ["Soleil", "Jupiter", "Saturne", "Uranus", "Neptune", "Pluton"]
colors = ["yellow","blue","green","red","orange","brown"]

nb_planets = len(planet)

def distance_soleil_ini():
    for i in range(1, nb_planets):
        print("La distance au 1er septembre 2012 a 00h00 entre " + planet[i] + " et le Soleil est : " + str(linalg.norm(pos_ini[i])) + " ua.")

def calc_energy(pos, vit):
    a = 0
    b = 0
    for i in range(nb_planets):
        a = a + masse[i] * linalg.norm(vit[i])**2
    for i in range(5):
        for j in range(i+1, 6):
            b = b + g * (masse[i] * masse[j])/(linalg.norm(pos[i] - pos[j]))
    return 1/2 * a - b

def calc_vit(i, pos):
    res = array([0.,0.,0.])
    for j in range (0, nb_planets):
        if i != j:
            res = res - g * masse[j] * (pos[i] - pos[j])/(linalg.norm(pos[i] - pos[j])**3)
    return res
    
def solaire_euler(dt, N):
    length = (int)(math.floor(N/dt))
    pos = zeros((length, nb_planets, 3))
    energy = zeros(length)
    pos_old = pos_ini
    vit_old = vit_ini
    pos[0] = pos_ini
    energy[0] = calc_energy(pos_ini, vit_ini)
    for i in range(1, length):
        pos_next = pos_old + dt * vit_old
        vit_next = vit_old + dt * array([calc_vit(0,pos_old),calc_vit(1,pos_old),calc_vit(2,pos_old),calc_vit(3,pos_old),calc_vit(4,pos_old),calc_vit(5,pos_old)])
        energy[i] = calc_energy(pos_next, vit_next)
        pos[i] = pos_next
        pos_old = pos_next
        vit_old = vit_next

    fig = plt.figure('Simulation_euler')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d([-40.,40.])
    ax.set_ylim3d([-40.,40.])
    ax.set_zlim3d([-20.,20.])
    for i in range(0, nb_planets):
        x_list = zeros(length)
        y_list = zeros(length)
        z_list = zeros(length)
        for j in range(0, length):
            x_list[j] = pos[j][i][0]
            y_list[j] = pos[j][i][1]
            z_list[j] = pos[j][i][2]
        ax.plot(x_list, y_list, z_list, '-', color=colors[i], label=planet[i])
        ax.legend()
    plt.show()
    
    fig = plt.figure('Hamiltonien_euler')
    t = linspace(0, N, length)
    plt.plot(t, energy, '-', label="Energie")
    plt.legend()
    plt.show()
    
def solaire_stormer_verlet(dt, N):
    length = (int)(math.floor(N/dt))
    pos = zeros((length, nb_planets, 3))
    energy = zeros(length)
    pos_old = pos_ini
    vit_old = vit_ini
    pos[0] = pos_ini
    energy[0] = calc_energy(pos_ini, vit_ini)
    for i in range(1, length):
        pos_inter = pos_old + dt/2 * vit_old
        vit_next = vit_old + dt * array([calc_vit(0,pos_old),calc_vit(1,pos_old),calc_vit(2,pos_old),calc_vit(3,pos_old),calc_vit(4,pos_old),calc_vit(5,pos_old)])
        pos_next = pos_inter + dt/2 * vit_next
        energy[i] = calc_energy(pos_next, vit_next)
        pos[i] = pos_next
        pos_old = pos_next
        vit_old = vit_next

    fig = plt.figure('Simulation_stormer_verlet')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d([-40.,40.])
    ax.set_ylim3d([-40.,40.])
    ax.set_zlim3d([-20.,20.])
    for i in range(0, 6):
        x_list = zeros(length)
        y_list = zeros(length)
        z_list = zeros(length)
        for j in range(0, length):
            x_list[j] = pos[j][i][0]
            y_list[j] = pos[j][i][1]
            z_list[j] = pos[j][i][2]
        ax.plot(x_list, y_list, z_list, '-', color=colors[i], label=planet[i])
        ax.legend()
    
    fig = plt.figure('Hamiltonien_stormer_verlet')
    t = linspace(0, N, length)
    plt.plot(t, energy, '-', label="Energie")
    plt.legend()
    plt.show()
    
def indice_min(tab):
    min = tab[100]
    ind_min = 100
    for i in range(101, len(tab)):
        if tab[i] < min:
            ind_min = i
            min = tab[i]
    return ind_min
    
def revolution(dt, N):  
    length = (int)(math.floor(N/dt))
    pos = zeros((length, nb_planets, 3))
    pos_old = pos_ini
    vit_old = vit_ini
    pos[0] = pos_ini
    distance_nep = zeros(length)
    distance_plu = zeros(length)
    
    for i in range(1, length):
        pos_inter = pos_old + dt/2 * vit_old
        vit_next = vit_old + dt * array([calc_vit(0,pos_old),calc_vit(1,pos_old),calc_vit(2,pos_old),calc_vit(3,pos_old),calc_vit(4,pos_old),calc_vit(5,pos_old)])
        pos_next = pos_inter + dt/2 * vit_next
        pos[i] = pos_next
        pos_old = pos_next
        vit_old = vit_next
        distance_nep[i] = linalg.norm(pos_next[4] - pos_ini[4])
        distance_plu[i] = linalg.norm(pos_next[5] - pos_ini[5])
    print("La periode de revolution de Neptune est de " + str(indice_min(distance_nep) * dt) + " jours.")
    print("La periode de revolution de Pluton est de " + str(indice_min(distance_plu) * dt) + " jours.")

def solaire_anime(dt, N):
    length = (int)(math.floor(N/dt))
    pos_old = pos_ini
    vit_old = vit_ini
    x = zeros((length, nb_planets))
    y = zeros((length, nb_planets))
    z = zeros((length, nb_planets))
    p = zeros((3, nb_planets))

    x[0] = array([pos_ini[0][0],pos_ini[1][0],pos_ini[2][0],pos_ini[3][0],pos_ini[4][0],pos_ini[5][0]])
    y[0] = array([pos_ini[0][1],pos_ini[1][1],pos_ini[2][1],pos_ini[3][1],pos_ini[4][1],pos_ini[5][1]])
    z[0] = array([pos_ini[0][2],pos_ini[1][2],pos_ini[2][2],pos_ini[3][2],pos_ini[4][2],pos_ini[5][2]])
    
    fig = plt.figure('Simulation_anime')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d([-40.,40.])
    ax.set_ylim3d([-40.,40.])
    ax.set_zlim3d([-20.,20.])
    
    lines = []
    pts = []
    
    for i in range(nb_planets):
        lines += ax.plot([], [], [],color=colors[i])
        pts += ax.plot([], [], [], marker="o",color=colors[i],label=planet[i])
        
    for n in range(1, length):
        pos_inter = pos_old + dt/2 * vit_old
        vit_next = vit_old + dt * array([calc_vit(0,pos_old),calc_vit(1,pos_old),calc_vit(2,pos_old),calc_vit(3,pos_old),calc_vit(4,pos_old),calc_vit(5,pos_old)])
        pos_next = pos_inter + dt/2 * vit_next
        
        x[n] = array([pos_next[0][0],pos_next[1][0],pos_next[2][0],pos_next[3][0],pos_next[4][0],pos_next[5][0]])
        y[n] = array([pos_next[0][1],pos_next[1][1],pos_next[2][1],pos_next[3][1],pos_next[4][1],pos_next[5][1]])
        z[n] = array([pos_next[0][2],pos_next[1][2],pos_next[2][2],pos_next[3][2],pos_next[4][2],pos_next[5][2]])
        
        pos_old = pos_next
        vit_old = vit_next
        
        for pt,line,i in zip(pts,lines,range(6)):
            line.set_data(x[:n,i],y[:n,i])
            line.set_3d_properties(z[:n,i])
            pt.set_data(x[n,i],y[n,i])
            pt.set_3d_properties(z[n,i])
            
        plt.pause(0.00001)
        ax.legend()
    
#distance_soleil_ini()
#solaire_euler(10,200000)
#solaire_euler_2(1000,20000)
solaire_stormer_verlet(10,200000)
#solaire_anime(10, 20000)
#revolution(10,200000)