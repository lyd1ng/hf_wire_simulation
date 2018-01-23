# What is the hf_wire_simulation?

The hf_wire_simulation is a simulation programm written in python2
to simulate the voltage, current and roh of a hf signal.
The result can be exported to animated gifs using gnuplot.
This project will be converted to python3 soon.

```
usage: main.py [-h] [-y GAMMA] [-R R] [-L L] [-C C] [-G G] [-l L] [-r R]
               [--Z_L Z_L] [--Zl ZL] [-f F] [-w W] [-U U] [-T T] [--dt DT]
               [-s NUMBER_OF_SAMPLES] [--plot_u PLOT_U] [--plot_i] [--plot_r]
               [-o OUTPUT_PATH] [--gif] [--gif_delay GIF_DELAY]
               [--gif_view_x GIF_VIEW_X] [--gif_view_z GIF_VIEW_Z]
               [--gif_scale GIF_SCALE] [--gif_scale_z GIF_SCALE_Z]

optional arguments:
  -h, --help            show this help message and exit
  -y GAMMA              Complex wire parameter
  -R R                  Ohm per cm
  -L L                  Inductivity per cm
  -C C                  Conductivty per cm
  -G G                  Moh per cm
  -l L                  Wire length
  -r R                  Reflection factor
  --Z_L Z_L             Characteristic impedance
  --Zl ZL               Terminal resistance
  -f F                  Signal frequency
  -w W                  Circular signal frequency
  -U U                  Voltage amplitude
  -T T                  Maximal simulation time
  --dt DT               Simulation time step
  -s NUMBER_OF_SAMPLES  Number of samples
  --plot_u PLOT_U       Plot voltage
  --plot_i              Plot current
  --plot_r              Plot Roh
  -o OUTPUT_PATH        Where to store the output files
  --gif                 Should the files be postprocessed using gnuplot
  --gif_delay GIF_DELAY
                        The gif animate delay
  --gif_view_x GIF_VIEW_X
                        The x view angle of the gif
  --gif_view_z GIF_VIEW_Z
                        The z view angle of the gif
  --gif_scale GIF_SCALE
                        The view scale of the gif
  --gif_scale_z GIF_SCALE_Z
                        The z view scale of the gif
```
