#######################
#  ######################
#    #   Leitungs-    #####
#    #   Simulation   #####
#  ######################
#######################
import os
import cmath
from argparse import ArgumentParser
import gplot_vector3 as gp

# Replace j with a more concret version of the imaginary entity
j = 0 + 1j

# The lists of gplot 3d vectors
u_vectors = list()
i_vectors = list()
r_vectors = list()

# Cause of a bug in the argparse module negative complexe numbers
# have to be seperated by a whitespace (in the parameter list)
# which cant be parsed thats  why their are stripped bevore
# using complex(..) to parse them


def parse_complex(cmplx_str):
    if cmplx_str is None:
        return None
    s = "".join(cmplx_str.split())
    return complex(s)

# This Function is needed to scale the gnuplot axes the right way
# if the postprocess flag is enabled


def get_max_abs(vectors):
    maximum = 0
    for v in vectors:
        if abs(v) > maximum:
            maximum = abs(v)
    return maximum

# Some parameters are redundant but at least one version has to
# be specified so they can be converted into the recommend one
# Unused parameters needs a value inequal to None cause the
# programme will terminate if one parameter is still None


def post_process_parameters(p):
    p.gamma = parse_complex(p.gamma)
    p.r = parse_complex(p.r)
    p.Z_L = parse_complex(p.Z_L)
    p.Zl = parse_complex(p.Zl)

    if p.T < 0:
        if p.dt is None:
            return False
        else:
            p.dt = p.dt * -1

    if p.w is None:
        if p.f is None:
            return False
        else:
            p.w = float(2 * cmath.pi * float(p.f))
    else:
        p.f = "unused"

    if p.gamma is None:
        if p.R is None or p.L is None or p.C is None or p.G is None:
            return False
        else:
            p.gamma = cmath.sqrt((p.R + p.w * p.L * j) * (p.G + p.w * p.C * j))
    else:
        p.R = p.L = p.C = p.G = "unused"

    if p.r is None:
        if p.Z_L is None or p.Zl is None:
            return False
        else:
            p.r = (p.Zl - p.Z_L) / (p.Zl + p.Z_L)
    else:
        p.Zl = "unused"

    if p.gif is True:
        if p.gif_delay < 0:
            return False

    for entry in vars(p).values():
        if entry is None:
            return False
    return True


def calculate_u(z, t):
    # First calculate the complex returning voltage
    Ur = (args.U * cmath.e**(-args.gamma * args.l) * args.r)
    U = (args.U * cmath.e**(-args.gamma * z) + Ur * cmath.e**(args.gamma * z)) * cmath.e**(j * args.w * t)
    return U.real, U.imag


def calculate_i(z, t):
    # First calculate the complex returning voltage
    Ur = (args.U * cmath.e**(-args.gamma * args.l) * args.r)
    Current = (1.0 / args.Z_L) * (args.U * cmath.e**(-args.gamma * z) - Ur * cmath.e**(args.gamma * z)) * cmath.e**(j * args.w * t)
    return Current.real, Current.imag


def calculate_r(z):
    r = args.r * cmath.e ** (-2 * args.gamma * (args.l - z))
    return r.real, r.imag


# Add all possible parameters
parser = ArgumentParser()
parser.add_argument("-y", dest="gamma",
                    help="Complex wire parameter", type=str)
parser.add_argument("-R", dest="R", help="Ohm per cm", type=float)
parser.add_argument("-L", dest="L", help="Inductivity per cm", type=float)
parser.add_argument("-C", dest="C", help="Conductivty per cm", type=float)
parser.add_argument("-G", dest="G", help="Moh per cm", type=float)
parser.add_argument("-l", dest="l", help="Wire length", type=float)
parser.add_argument("-r", dest="r", help="Reflection factor", type=str)
parser.add_argument("--Z_L", dest="Z_L",
                    help="Characteristic impedance", type=str)
parser.add_argument("--Zl", dest="Zl", help="Terminal resistance", type=str)
parser.add_argument("-f", dest="f", help="Signal frequency", type=float)
parser.add_argument(
    "-w", dest="w", help="Circular signal frequency", type=float)
parser.add_argument("-U", dest="U", default="0",
                    help="Voltage amplitude", type=float)
parser.add_argument("-T", dest="T", help="Maximal simulation time", type=float)
parser.add_argument("--dt", dest="dt", help="Simulation time step", type=float)
parser.add_argument("-s", dest="number_of_samples",
                    help="Number of samples", type=int)
parser.add_argument("--plot_u", dest="plot_u", default=True,
                    type=bool, help="Plot voltage")
parser.add_argument("--plot_i", dest="plot_i", default=False,
                    action="store_true", help="Plot current")
parser.add_argument("--plot_r", dest="plot_r", default=False,
                    action="store_true", help="Plot Roh")
parser.add_argument("-o", dest="output_path",
                    help="Where to store the output files")
parser.add_argument("--gif", dest="gif", default=False, action="store_true",
                    help="Should the files be postprocessed using gnuplot")
parser.add_argument("--gif_delay", dest="gif_delay",
                    type=float, default=-1, help="The gif animate delay")
parser.add_argument("--gif_view_x", dest="gif_view_x",
                    default="45.0", help="The x view angle of the gif")
parser.add_argument("--gif_view_z", dest="gif_view_z",
                    default="45.0", help="The z view angle of the gif")
parser.add_argument("--gif_scale", dest="gif_scale",
                    default="1.0", help="The view scale of the gif")
parser.add_argument("--gif_scale_z", dest="gif_scale_z",
                    default="1.0", help="The z view scale of the gif")

# Parse parameters
args = parser.parse_args()
if post_process_parameters(args) is False:
    print("Parameters missing!")
    print("use -h/--help for further help")
    exit(-1)

# Calculate all vectors and store them in lists of gnuplot_vector3 instances
t = 0
distance_of_vectors = args.l / args.number_of_samples
while t < args.T:
    for i in range(0, args.number_of_samples):
        if args.plot_u:
            u_vectors.append(gp.gplot_vector3(
                i * distance_of_vectors, 0, 0, 0, 0, 0))
            u_vectors[-1].dz, u_vectors[-1].dy = calculate_u(
                i * distance_of_vectors, t)
        if args.plot_i:
            i_vectors.append(gp.gplot_vector3(
                i * distance_of_vectors, 0, 0, 0, 0, 0))
            i_vectors[-1].dz, i_vectors[-1].dy = calculate_i(
                i * distance_of_vectors, t)
        if args.plot_r:
            r_vectors.append(gp.gplot_vector3(
                i * distance_of_vectors, 0, 0, 0, 0, 0))
            r_vectors[-1].dz, r_vectors[-1].dy = calculate_r(
                i * distance_of_vectors)
    t = t + args.dt

# Write voltage vectors to file
if args.plot_u:
    path = args.output_path + ".voltage"
    fd = open(path, "w")
    for i in range(1, len(u_vectors)):
        fd.write(str(u_vectors[i - 1]) + "\n")
        if i % args.number_of_samples == 0:
            fd.write("\n\n")
    fd.close()

# Write current vectors to file
if args.plot_i:
    path = args.output_path + ".current"
    fd = open(path, "w")
    for i in range(1, len(i_vectors)):
        fd.write(str(i_vectors[i - 1]) + "\n")
        if i % args.number_of_samples == 0:
            fd.write("\n\n")
    fd.close()

# Write roh vectors to file
if args.plot_r:
    path = args.output_path + ".roh"
    fd = open(path, "w")
    for i in range(1, len(r_vectors)):
        fd.write(str(r_vectors[i - 1]) + "\n")
        if i % args.number_of_samples == 0:
            fd.write("\n\n")
    fd.close()

# Invoke the postprocess.plg script if desired
if args.gif:
    if args.plot_u:
        output_path = args.output_path + ".voltage"
        gif_path = args.output_path + ".voltage.gif"
        maximum = str(get_max_abs(u_vectors))
        os.system("gnuplot -e \"\
            in1='" + output_path + "';\
            in2='0';\
            in3='0';\
            out='" + gif_path + "';\
            range='" + maximum + "';\
            color1='blue';\
            color2='blue';\
            color3='blue';\
            delay='" + str(args.gif_delay) + "';\
            x_angle='" + str(args.gif_view_x) + "';\
            z_angle='" + str(args.gif_view_z) + "';\
            scale='" + str(args.gif_scale) + "';\
            scale_z='" + str(args.gif_scale_z) + "';\
            omega='" + str(args.w) + "';\
            gamma='" + str(args.gamma) + "';\
            r='" + str(args.r) + "'\" postprocess.plg")

    if args.plot_i:
        output_path = args.output_path + ".current"
        gif_path = args.output_path + ".current.gif"
        maximum = str(get_max_abs(i_vectors))
        os.system("gnuplot -e \"\
            in1='" + output_path + "';\
            in2='0';\
            in3='0';\
            out='" + gif_path + "';\
            range='" + maximum + "';\
            color1='red';\
            color2='red';\
            color3='red';\
            delay='" + str(args.gif_delay) + "';\
            x_angle='" + str(args.gif_view_x) + "';\
            z_angle='" + str(args.gif_view_z) + "';\
            scale='" + str(args.gif_scale) + "';\
            scale_z='" + str(args.gif_scale_z) + "';\
            omega='" + str(args.w) + "';\
            gamma='" + str(args.gamma) + "';\
            r='" + str(args.r) + "'\" postprocess.plg")

    if args.plot_r:
        output_path = args.output_path + ".roh"
        gif_path = args.output_path + ".roh.gif"
        maximum = str(get_max_abs(r_vectors))
        os.system("gnuplot -e \"\
            in1='" + output_path + "';\
            in2='0';\
            in3='0';\
            out='" + gif_path + "';\
            range='" + maximum + "';\
            color1='black';\
            color2='black';\
            color3='black';\
            delay='" + str(args.gif_delay) + "';\
            x_angle='" + str(args.gif_view_x) + "';\
            z_angle='" + str(args.gif_view_z) + "';\
            scale='" + str(args.gif_scale) + "';\
            scale_z='" + str(args.gif_scale_z) + "';\
            omega='" + str(args.w) + "';\
            gamma='" + str(args.gamma) + "';\
            r='" + str(args.r) + "'\" postprocess.plg")
