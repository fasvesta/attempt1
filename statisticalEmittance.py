# - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * 
"""
#
#   statisticalEmittance
#   A module to calculate the transverse emittances from the beam
#   based on pyORBIT & desy-thesis-05-014
#
#   Version : 1.0
#   Author  : F. Asvesta
#   Contact : fasvesta .at. cern .dot. ch
#
"""
# - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - *

import numpy as np

class statisticalEmittance(object):
    """
    class for the statistical emittance estimation
    Returns: statisticalEmittance instance
    """

    def __init__(self, x=None, px=None, y=None, py=None, z=None, dp=None):
        """
        Initialization function
        Input : particle coordinates as 1D numpy arrays
        Returns: void
        """   
        self.coordinateMatrix=np.array([x,px,y,py,z,dp])
        self.beamMatrix=np.matmul(self.coordinateMatrix,self.coordinateMatrix.T)/len(self.x)
        self.dispersionX=None
        self.coordinateMatrixBetatronic = None

        def correlation(self,par1,par2, betatronic=True):
            """
            Calculation of the correlation for the beam matrices
            Inputs:   par1 : [0|1|2|3|4|5]
                      par2 : [0|1|2|3|4|5]
                             integers corresponding to coordinates (0->x),(1->px), (2->y), (3->py),(4->z), (5->dp)
                      betatronic : [bool] if True the betatronic matrices are considered (default=True)
            Returns: <(a-<a>)*(b-<b>)> = <a*b> - <a>*<b>
            """
            if par1 in range(0,6) and par2 in range(0,6):
                if betatronic:
                    if self.coordinateMatrixBetatronic is None:
                        self.betatronicMatrices()
                    return self.beamMatrixBetatronic[par1,par2]-np.nanmean(self.coordinateMatrixBetatronic[par1])*np.nanmean(self.coordinateMatrixBetatronic[par2])
                else:
                    return self.beamMatrix[par1,par2]-np.nanmean(self.coordinateMatrix[par1])*np.nanmean(self.coordinateMatrix[par2])
            else:
                raise IOError('# statisticalEmittance::correlation: par1 and par2 need to be [0|1|2|3|4|5]')

        #in pyorbit the dispersion is returned *p*b and then for the betatronic they do /E/b**2 / for me it's not needed we can have the same without the momentum similar to http://www-library.desy.de/preparch/desy/thesis/desy-thesis-05-014.pdf (13) 

        def calculateDispersion(self):
            """
            Statistical dispersion evaluation
            Returns: void
            """
            self.corr5=self.correlation(5,5, betatronic=False)
            self.dispersionX=self.correlation(0,5, betatronic=False)/self.corr5
            self.dispersionPx=self.correlation(1,5, betatronic=False)/self.corr5
            self.dispersionY=self.correlation(2,5, betatronic=False)/self.corr5
            self.dispersionPy=self.correlation(3,5, betatronic=False)/self.corr5

        def betatronicMatrices(self):
            """
            Evaluation of the coordinates and beam matrix excluding dispersive components
            Returns: void
            """
            if self.dispersionX is None:
                self.calculateDispersion()

            xBetatronic=self.coordinateMatrix[0]-self.dispersionX*self.coordinateMatrix[5]
            pxBetatronic=self.coordinateMatrix[1]-self.dispersionPx*self.coordinateMatrix[5]
            yBetatronic=self.coordinateMatrix[2]-self.dispersionY*self.coordinateMatrix[5]
            pyBetatronic=self.coordinateMatrix[3]-self.dispersionPy*self.coordinateMatrix[5]

            self.coordinateMatrixBetatronic=np.array([xBetatronic,pxBetatronic,yBetatronic,pyBetatronic,self.coordinateMatrix[4],self.coordinateMatrix[5]])
            self.beamMatrixBetatronic=np.matmul(self.coordinateMatrixBetatronic,self.coordinateMatrixBetatronic.T)/len(xBetatronic)

        # self.twissBetax=self.coordinateMatrixBetatronic[0,0]

    

    # def estimateBetaFunction(particles):
