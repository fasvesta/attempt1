# - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * 
"""
#
#  StatisticalEmittance
#   A module to calculate the transverse emittances from the beam
#   based on PyORBIT & desy-thesis-05-014
#
#   Version : 1.0.1
#   Author  : F. Asvesta
#   Contact : fasvesta .at. cern .dot. ch
"""
# - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

import numpy as np

class StatisticalEmittance(object):
    """
    class for the statistical emittance estimation
    Returns:StatisticalEmittance instance
    """

    def __init__(self, particles=None):
        """
        Initialization function
        Input: particles:  [xpart distribution object]
        Returns: void
        """ 
        self.coordinate_matrix=None  

        if particles:
            self.set_particles(particles)
        else:
            self.beam_matrix=None
            self.beta0=None
            self.gamma0=None
            print("#StatisticalEmittance : Provide distribution in [set_particles]")
        self.dx=None
        self.coordinate_matrix_betatronic = None
        self.emitt_x = None
        self.betx = None
        self.emitt_4d = None
        self.coupling = None

    def set_particles(self,particles):
        """
        Provide distribution to calculate emittances and optics
        Input:  particles:  [xpart distribution object]
        Returns: void
        """
        if self.coordinate_matrix is not None:
            self.dx=None
            self.coordinate_matrix_betatronic = None
            self.emitt_x = None
            self.betx = None
            self.emitt_4d = None
            self.coupling = None

        mask_alive = particles.state>=1           
        self.coordinate_matrix=np.array([particles.x[mask_alive],particles.px[mask_alive],particles.y[mask_alive],particles.py[mask_alive],particles.zeta[mask_alive],particles.delta[mask_alive]])
        self.beam_matrix=np.matmul(self.coordinate_matrix,self.coordinate_matrix.T)/len(particles.x[mask_alive])
        self.beta0=particles.beta0[0]
        self.gamma0=particles.gamma0[0]

    def correlation(self,par1,par2, betatronic=True):
        """
        Calculation of the correlation for the beam matrices
        Inputs: par1 : [0|1|2|3|4|5]
                par2 : [0|1|2|3|4|5]
                integers corresponding to coordinates (0->x), (1->px), (2->y), (3->py), (4->z), (5->dp)
                betatronic : [bool] if True the betatronic matrices are considered (default=True)
        Returns: <(a-<a>)*(b-<b>)> = <a*b> - <a>*<b>
        """
        if par1 in range(0,6) and par2 in range(0,6):
            if betatronic:
                if par1>3 or par2>3:
                    raise IOError('#StatisticalEmittance::correlation: if betatronic par1 and par2 need to be [0|1|2|3]')
                elif self.coordinate_matrix_betatronic is None:
                    self.betatronic_matrices()
                return self.beam_matrix_betatronic[par1,par2]-np.nanmean(self.coordinate_matrix_betatronic[par1])*np.nanmean(self.coordinate_matrix_betatronic[par2])
            else:
                return self.beam_matrix[par1,par2]-np.nanmean(self.coordinate_matrix[par1])*np.nanmean(self.coordinate_matrix[par2])
        else:
            raise IOError('#StatisticalEmittance::correlation: par1 and par2 need to be [0|1|2|3|4|5]')
    
    def betatronic_matrices(self):
        """
        Evaluation of the coordinates and beam matrix excluding dispersive components
        Returns: void
        """
        if self.dx is None:
            self.calculate_dispersion()

        x_betatronic=self.coordinate_matrix[0]-self.dx*self.coordinate_matrix[5]
        px_betatronic=self.coordinate_matrix[1]-self.dpx*self.coordinate_matrix[5]
        y_betatronic=self.coordinate_matrix[2]-self.dy*self.coordinate_matrix[5]
        py_betatronic=self.coordinate_matrix[3]-self.dpy*self.coordinate_matrix[5]

        self.coordinate_matrix_betatronic=np.array([x_betatronic,px_betatronic,y_betatronic,py_betatronic])
        self.beam_matrix_betatronic=self.beam_matrix-self.dispersion_table*self.corr5

    def calculate_dispersion(self):
        """
        Statistical dispersion evaluation
        Returns: void
        """
        self.corr5=self.correlation(5,5, betatronic=False)
        self.dx=self.correlation(0,5, betatronic=False)/self.corr5
        self.dpx=self.correlation(1,5, betatronic=False)/self.corr5
        self.dy=self.correlation(2,5, betatronic=False)/self.corr5
        self.dpy=self.correlation(3,5, betatronic=False)/self.corr5
        disp_table=np.array([[self.dx],[self.dpx],[self.dy],[self.dpy],[0],[0]])
        self.dispersion_table=np.matmul(disp_table,disp_table.T)
    
    def calculate_emittance(self, fourD=False):
        """
        Transverse emittance evaluation 
        Returns: void
        """
        if self.emitt_x is None:
            self.x_matrix=np.array([[self.correlation(0,0, betatronic=True),self.correlation(0,1, betatronic=True)],[self.correlation(1,0, betatronic=True),self.correlation(1,1, betatronic=True)]])
            self.emitt_x=np.sqrt(np.abs(np.linalg.det(self.x_matrix)))
            self.y_matrix=np.array([[self.correlation(2,2, betatronic=True),self.correlation(2,3, betatronic=True)],[self.correlation(3,2, betatronic=True),self.correlation(3,3, betatronic=True)]])
            self.emitt_y=np.sqrt(np.abs(np.linalg.det(self.y_matrix)))
        if fourD:
            x_y_matrix=np.array([[self.correlation(0,2, betatronic=True),self.correlation(0,3, betatronic=True)],[self.correlation(1,2, betatronic=True),self.correlation(1,3, betatronic=True)]])
            full_matrix=np.append(np.append(self.x_matrix,x_y_matrix,axis=1),np.append(x_y_matrix.T,self.y_matrix,axis=1),axis=0)
            self.emitt_4d=np.sqrt(np.linalg.det(full_matrix))

    def calculate_coupling_factor(self):
        """
        Coupling evaluation as described in doi: 10.18429/JACoW-LINAC2018-THPO118
        Returns: void
        """
        if self.emitt_4d is None:
            self.calculate_emittance(fourD=True)
        self.coupling=self.emitt_x*self.emitt_y/self.emitt_4d - 1.
    
    def calculate_twiss_functions(self):
        """
        Twiss functions evaluation 
        Returns: void
        """
        if self.emitt_x is None:
            self.calculate_emittance()
        self.betx = self.x_matrix[0,0]/self.emitt_x
        self.alfx = - self.x_matrix[0,1]/self.emitt_x
        self.gamx = self.x_matrix[1,1]/self.emitt_x
        self.bety = self.y_matrix[0,0]/self.emitt_y
        self.alfy = - self.y_matrix[0,1]/self.emitt_y
        self.gamy = self.y_matrix[1,1]/self.emitt_y
    
    def measure_bunch_moments(self, particles, coupling = False):
        self.set_particles(particles)
        if coupling:
            self.calculate_emittance(fourD=True)
            self.calculate_coupling_factor()
        else:
            self.calculate_emittance()
        self.calculate_twiss_functions()
        self.bunch_moments={'nemitt_x': self.emitt_x*self.beta0*self.gamma0, 'nemitt_y': self.emitt_y*self.beta0*self.gamma0, 
                            'betx': self.betx, 'bety': self.bety, 
                            'alfx': self.alfx , 'alfy': self.alfy,
                            'gamx': self.gamx , 'gamy': self.gamy,
                            'dx': self.dx, 'dy': self.dy,
                            'dpx': self.dpx , 'dpy': self.dpy}
        if coupling:
            self.bunch_moments['coupling']=self.coupling
            self.bunch_moments['emitt_4d']=self.emitt_4d
        return self.bunch_moments
