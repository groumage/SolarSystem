from matplotlib import pyplot as plt
from numpy import *
import sys

# sun and 4 first planet are mix into a single element
masse = array(
    [1.00000597682, 0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 1 / (1.3e8)])

pos_ini = array([[0, 0, 0], [2.2911898857058, 4.4742096428394, -0.0698526085811],
                 [-8.4077274422046, -4.9386528867475, 0.4206194745212],
                 [19.9506411858858, -2.1189594536568, -0.2505524004432],
                 [26.3841113101207, -14.2637141476417, -0.3141826415529],
                 [4.7240342661351, -31.8776619468451, 2.0454833905119]])

vit_ini = array([[0, 0, 0], [-0.0068140485228, 0.0038022184356, 0.0001366932116],
                 [0.0025181864149, -0.0048203311629, -0.000016441567],
                 [-0.0004495229028, 0.0037304364206, 0.0000197241717],
                 [0.0014665671695, 0.0027825194375, -0.0000910521187],
                 [0.0031588952992, -0.0001620298531, -0.0008976519916]])

g = 2.95912208286e-4
x_y_limit = 40.0
z_limit = 20.0

planets = ["Sun", "Jupiter", "Saturne", "Uranus", "Neptune", "Pluton"]
colors = ["yellow", "blue", "green", "red", "orange", "brown"]


def distance_sun():
    for i in range(1, len(planets)):
        print("The distance on September 1st 2012 at 00:00 between the sun and " + planets[i] + " is : " + str(
            linalg.norm(pos_ini[i])) + " au (astronomic unit).")


def energy(pos, vit):
    term1 = 0
    term2 = 0
    for i in range(len(planets)):
        term1 = term1 + masse[i] * linalg.norm(vit[i]) ** 2
    for i in range(5):
        for j in range(i + 1, 6):
            term2 = term2 + g * (masse[i] * masse[j]) / (linalg.norm(pos[i] - pos[j]))
    return 1 / 2 * term1 - term2


def vitesse(i, pos):
    res = array([0., 0., 0.])
    for j in range(len(planets)):
        if i != j:
            res = res - g * masse[j] * (pos[i] - pos[j]) / (linalg.norm(pos[i] - pos[j]) ** 3)
    return res


# euler schema
def solaire_euler(dt, N, save_figure):
    length = int(math.floor(N / dt))
    pos = zeros((length, len(planets), 3))
    storage_energy = zeros(length)
    pos_old = pos_ini
    vit_old = vit_ini
    pos[0] = pos_ini
    storage_energy[0] = energy(pos_ini, vit_ini)

    for i in range(length):
        pos_next = pos_old + dt * vit_old
        vit_next = vit_old + dt * array([vitesse(i, pos_old) for i in range(len(planets))])
        storage_energy[i] = energy(pos_next, vit_next)
        pos[i] = pos_next
        pos_old = pos_next
        vit_old = vit_next

    fig = plt.figure('Simulation euler explicite')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d([-40., 40.])
    ax.set_ylim3d([-40., 40.])
    ax.set_zlim3d([-20., 20.])

    for i in range(len(planets)):
        x_list = zeros(length)
        y_list = zeros(length)
        z_list = zeros(length)
        for j in range(length):
            x_list[j] = pos[j][i][0]
            y_list[j] = pos[j][i][1]
            z_list[j] = pos[j][i][2]
        ax.plot(x_list, y_list, z_list, '-', color=colors[i], label=planets[i])
        ax.legend()

    if save_figure:
        plt.savefig('../ressources/euler-explicite_' + str(dt) + '_' + str(N) + '_simu.eps', format='eps')
    plt.figure('Hamiltonien euler explicite')
    t = linspace(0, N, length)
    plt.plot(t, storage_energy, '-', label="Energie")
    plt.xlabel("Jours")
    plt.legend()
    if save_figure:
        plt.savefig('../ressources/euler-explicite_' + str(dt) + '_' + str(N) + '_energie.eps', format='eps')
    plt.show()


# stormer-verlet schema
def solaire_stormer_verlet(dt, N, save_figure):
    length = int(math.floor(N / dt))
    pos = zeros((length, len(planets), 3))
    storage_energy = zeros(length)
    pos_old = pos_ini
    vit_old = vit_ini
    pos[0] = pos_ini
    storage_energy[0] = energy(pos_ini, vit_ini)
    for i in range(length):
        pos_inter = pos_old + dt / 2 * vit_old
        vit_next = vit_old + dt * array([vitesse(i, pos_old) for i in range(len(planets))])
        pos_next = pos_inter + dt / 2 * vit_next
        storage_energy[i] = energy(pos_next, vit_next)
        pos[i] = pos_next
        pos_old = pos_next
        vit_old = vit_next

    fig = plt.figure('Simulation stormer verlet')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d([-x_y_limit, x_y_limit])
    ax.set_ylim3d([-x_y_limit, x_y_limit])
    ax.set_zlim3d([-z_limit, z_limit])

    for i in range(len(planets)):
        x_list = zeros(length)
        y_list = zeros(length)
        z_list = zeros(length)
        for j in range(length):
            x_list[j] = pos[j][i][0]
            y_list[j] = pos[j][i][1]
            z_list[j] = pos[j][i][2]
        ax.plot(x_list, y_list, z_list, '-', color=colors[i], label=planets[i])
        ax.legend()

    if save_figure:
        plt.savefig('../ressources/stormer-verlet_' + str(dt) + '_' + str(N) + '_simu.eps', format='eps')
    fig = plt.figure('Hamiltonien stormer verlet')
    t = linspace(0, N, length)
    plt.plot(t, storage_energy, '-', label="Energie")
    plt.xlabel("Jours")
    plt.legend()
    if save_figure:
        plt.savefig('../ressources/stormer-verlet_' + str(dt) + '_' + str(N) + '_energie.eps', format='eps')
    plt.show()


def revolution(dt, N):
    length = int(math.floor(N / dt))
    pos = zeros((length, len(planets), 3))
    pos_old = pos_ini
    vit_old = vit_ini
    pos[0] = pos_ini
    distance_from_init = zeros((len(planets), length))
    for i in range(length):
        pos_inter = pos_old + dt / 2 * vit_old
        vit_next = vit_old + dt * array([vitesse(i, pos_old) for i in range(len(planets))])
        pos_next = pos_inter + dt / 2 * vit_next
        pos[i] = pos_next
        pos_old = pos_next
        vit_old = vit_next
        for j in range(len(planets)):
            distance_from_init[j][i] = linalg.norm(pos_next[j] - pos_ini[j])
    for i in range(1, len(planets)):
        idx = where(distance_from_init[i] == distance_from_init[i][int(length * 0.01):].min())[0][0]
        print("Revolution of " + str(planets[i]) + " is " + str(idx * dt) + " days.")


# use stormer-verlet schema
def solaire_anime(dt, N, dt_anim):
    length = (int)(math.floor(N / dt))
    pos_old = pos_ini
    vit_old = vit_ini
    x = zeros((length, len(planets)))
    y = zeros((length, len(planets)))
    z = zeros((length, len(planets)))
    p = zeros((3, len(planets)))

    x[0] = array([pos_ini[0][0], pos_ini[1][0], pos_ini[2][0], pos_ini[3][0], pos_ini[4][0], pos_ini[5][0]])
    y[0] = array([pos_ini[0][1], pos_ini[1][1], pos_ini[2][1], pos_ini[3][1], pos_ini[4][1], pos_ini[5][1]])
    z[0] = array([pos_ini[0][2], pos_ini[1][2], pos_ini[2][2], pos_ini[3][2], pos_ini[4][2], pos_ini[5][2]])

    fig = plt.figure('Simulation_anime')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d([-40., 40.])
    ax.set_ylim3d([-40., 40.])
    ax.set_zlim3d([-20., 20.])

    lines = []
    pts = []

    for i in range(len(planets)):
        lines += ax.plot([], [], [], color=colors[i])
        pts += ax.plot([], [], [], marker="o", color=colors[i], label=planets[i])

    for n in range(1, length):
        pos_inter = pos_old + dt / 2 * vit_old
        vit_next = vit_old + dt * array([vitesse(i, pos_old) for i in range(len(planets))])
        pos_next = pos_inter + dt / 2 * vit_next

        x[n] = array([pos_next[0][0], pos_next[1][0], pos_next[2][0], pos_next[3][0], pos_next[4][0], pos_next[5][0]])
        y[n] = array([pos_next[0][1], pos_next[1][1], pos_next[2][1], pos_next[3][1], pos_next[4][1], pos_next[5][1]])
        z[n] = array([pos_next[0][2], pos_next[1][2], pos_next[2][2], pos_next[3][2], pos_next[4][2], pos_next[5][2]])

        pos_old = pos_next
        vit_old = vit_next

        if n % dt_anim == 0:
            ax.legend()
            for pt, line, i in zip(pts, lines, range(6)):
                line.set_data(x[:n, i], y[:n, i])
                line.set_3d_properties(z[:n, i])
                pt.set_data(x[n, i], y[n, i])
                pt.set_3d_properties(z[n, i])
            plt.pause(0.00001)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("usage:", sys.argv[0], "soleil | euler step_time time_interval [--save] | stormer-verlet step_time time_interval [--save] | anim step_time time_interval interval_print | revolution step_time time_interval")
        print("\tstep_time: the dt value (for instance: 10)")
        print("\ttime_interval: the N value (for instance: 200 000)")
        print("\tinterval_anim: time between frame of anim (for instance: 200)")
        sys.exit(1)
    if sys.argv[1] == "--help":
        print("usage:", sys.argv[0], "sun | euler step_time time_interval [--save] | stormer-verlet step_time time_interval [--save] | anim step_time time_interval interval_anim | revolution step_time time_interval")
        print("\tstep_time: the dt value (for instance: 10)")
        print("\ttime_interval: the N value (for instance: 200 000)")
        print("\tinterval_anim: time between frame of anim (for instance: 200)")
        sys.exit(1)
    if sys.argv[1] == "sun":
        distance_sun()
    elif sys.argv[1] == "euler":
        if len(sys.argv) != 4 and len(sys.argv) != 5:
            print("usage: euler step_time time_interval [--save]")
            sys.exit(1)
        if len(sys.argv) == 5 and sys.argv[4] == "--save":
            solaire_euler(int(sys.argv[2]), int(sys.argv[3]), True)
        elif len(sys.argv) == 5 and sys.argv[4] != "--save":
            print("usage: euler step_time time_interval [--save]")
            sys.exit(1)
        else:
            solaire_euler(int(sys.argv[2]), int(sys.argv[3]), False)
    elif sys.argv[1] == "stormer-verlet":
        if len(sys.argv) != 4 and len(sys.argv) != 5:
            print("usage: stormer-verlet step_time time_interval [--save]")
            sys.exit(1)
        if len(sys.argv) == 5 and sys.argv[4] == "--save":
            solaire_stormer_verlet(int(sys.argv[2]), int(sys.argv[3]), True)
        elif len(sys.argv) == 5 and sys.argv[4] != "--save":
            print("usage: stormer-verlet step_time time_interval [--save]")
            sys.exit(1)
        else:
            solaire_stormer_verlet(int(sys.argv[2]), int(sys.argv[3]), False)
    elif sys.argv[1] == "anim":
        if len(sys.argv) != 5:
            print("usage: anim step_time time_interval interval_print")
            sys.exit(1)
        solaire_anime(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    elif sys.argv[1] == "revolution":
        if len(sys.argv) != 4:
            print("usage: revolution step_time time_interval")
            sys.exit(1)
        revolution(int(sys.argv[2]), int(sys.argv[3]))
