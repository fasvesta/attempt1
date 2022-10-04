import numpy as np
from cpymad.madx import Madx
import xtrack as xt
import xpart as xp
from statisticalEmittance import *
# import json
import xobjects as xo
import xfields as xf

####################
# Choose a context #
####################

context = xo.ContextCpu()
# context = xo.ContextCupy()
#context = xo.ContextPyopencl('0.0')

print(context)

mad = Madx()
mad.call('PSB/madx/psb_injection_example.madx')

line= xt.Line.from_madx_sequence(mad.sequence['psb'])
line.particle_ref=xp.Particles(mass0=xp.PROTON_MASS_EV,
                               gamma0=mad.sequence.psb.beam.gamma)

nemitt_x=5.4e-6
nemitt_y=3.5e-6
bunch_intensity=5e11
sigma_z=10.67
n_part = int(3e3)

# from space charge example
num_turns=600 # is this the number of turns to track?

num_spacecharge_interactions = 202 # is this interactions per turn?
tol_spacecharge_position = 1e-2 # is this the minimum/maximum space between sc elements?

# Available modes: frozen/quasi-frozen/pic
mode = 'frozen'

#############################################
# Install spacecharge interactions (frozen) #
#############################################

lprofile = xf.LongitudinalProfileQGaussian(
        number_of_particles=bunch_intensity,
        sigma_z=sigma_z,
        z0=0.,
        q_parameter=1.)

xf.install_spacecharge_frozen(line=line,
                   particle_ref=line.particle_ref,
                   longitudinal_profile=lprofile,
                   nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                   sigma_z=sigma_z,
                   num_spacecharge_interactions=num_spacecharge_interactions,
                   tol_spacecharge_position=tol_spacecharge_position)


#################
# Build Tracker #
#################

tracker = xt.Tracker(_context=context,
                    line=line)
tracker_sc_off = tracker.filter_elements(exclude_types_starting_with='SpaceCh')

######################
# Generate particles #
######################

# particles = xp.generate_matched_gaussian_bunch(_context=context, num_particles=n_part,
#                             total_intensity_particles=bunch_intensity,
#                             nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
#                             particle_ref=line.particle_ref,
#                             tracker=tracker_sc_off)

x_norm, y_norm, _, _ = xp.generate_2D_polar_grid(
    theta_range=(0.01, np.pi/2-0.01),
    ntheta = 20,
    r_range = (0.1, 3),
    nr = 30)

particles = xp.build_particles(tracker=tracker, particle_ref=line.particle_ref,
                               x_norm=x_norm, y_norm=y_norm, delta=0,
                               scale_with_transverse_norm_emitt=(nemitt_x, nemitt_y))

tracker.track(particles, num_turns=num_turns, turn_by_turn_monitor=True )

np.save('x',tracker.record_last_track.x)
np.save('px',tracker.record_last_track.px)
np.save('y',tracker.record_last_track.y)
np.save('py',tracker.record_last_track.py)
np.save('z',tracker.record_last_track.zeta)
np.save('d',tracker.record_last_track.delta)

