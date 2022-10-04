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

nemitt_x=3e-6
nemitt_y=1e-6
bunch_intensity=5e11
sigma_z=15.
n_part = 3e3

# from space charge example
num_turns=5 # is this the number of turns to track?

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
                   particle_ref=particle_ref,
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

p_gaussian = xp.generate_matched_gaussian_bunch(_context=context, num_particles=n_part,
                            total_intensity_particles=bunch_intensity,
                            nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
                            particle_ref=line.particle_ref,
                            tracker=tracker_sc_off)

#########
# Track #
#########

# r=StatisticalEmittance()
# epsn_x = []
# epsn_y = []
tracker.track(p_gaussian, num_turns=num_turns, turn_by_turn_monitor=True )


for ii in range(num_turns):
    bunch_moments=r.measure_bunch_moments(p_gaussian)
    epsn_x.append(bunch_moments['nemitt_x'])
    epsn_y.append(bunch_moments['nemitt_y'])

print('epsn_x = ',epsn_x)
print('epsn_y = ',epsn_y)
