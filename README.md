# Simulation of the Solar System

---
This is a student project in numerical analysis. The aim is to simulate the Solar System. We simulate the 5 main body
(Sun, Jupiter, Saturn and Uranus) of the solar system and Pluto. So this is a 6-body problem. The only forces taken into account are
gravitational interactions.

Through this project we manage to have a first experience with Python and write the implementation of numerical schemas.
We used two schemas :
* forward Euler method
* Stormer-Verlet method.

![Image of solar system](https://github.com/groumage/SolarSystem/blob/master/ressources/solar_system_simulation.png?raw=true)
The above screenshot is generate using ``python3 solar-system.py stormer-verlet 1 45000``.

This project is implemented in Python, using mainly ``numpy`` and ``matplotlib``. The
[subject](https://github.com/groumage/SolarSystem/blob/master/ressources/sujet.pdf) of this project is available, as well as
the [report](https://github.com/groumage/SolarSystem/blob/master/ressources/rapport.pdf) (written in French).

The report is also available at [Overleaf](https://www.overleaf.com/read/kybhbxbcxtsd).

Quick execution (assuming all necessary module installed):  
``python3 src/solar-system.py stormer-verlet 100 200000``  
``python3 src/solar-system.py anim 1 100000 200``