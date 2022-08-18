# statistical emittance estimation
# python implementation of the emittance calculation from pyOrbit

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
        self.x=x
        self.px=px
        self.y=y
        self.py=py
        self.z=z
        self.dp=dp        
        self.coordinateMatrix=np.array([self.x,self.px,self.y,self.py,self.z,self.dp])
        self.beamMatrix=np.matmul(self.coordinateMatrix,self.coordinateMatrix.T)/len(self.x)
        self.dispersionX=None
        self.xBetatronic=None

        #in pyorbit the dispersion is returned *p*b and then for the betatronic they do /E/b**2 / for me it's not needed we can have the same without the momentum similar to http://www-library.desy.de/preparch/desy/thesis/desy-thesis-05-014.pdf (13) 

        def correlation(self,par1,par2, betatronic=True):
            """
            Calculation of the correlation for the beam matrices
            Inputs:   par1 : [0|1|2|3|4|5]
                      par2 : [0|1|2|3|4|5]
                             integers corresponding to coordinates (0->x),(1->px), (2->y), (3->py),(4->z), (5->dp)
                      betatronic : [bool] if True the betatronic matrices are considered (default=True)
            Returns: <(x-<x>)*(y-<y>)> = <x*y> - <x>*<y>
            """
            if par1 in range(0,6):
                if par2 in range(0,6):
                    if betatronic:
                        if self.xBetatronic is None:
                            self.betatronicMatrices()
                        return self.beamMatrixBetatronic[par1,par2]-np.nanmean(self.coordinateMatrixBetatronic[par1])*np.nanmean(self.coordinateMatrixBetatronic[par2])
                    else:
                        return self.beamMatrix[par1,par2]-np.nanmean(self.coordinateMatrix[par1])*np.nanmean(self.coordinateMatrix[par2])
                else:
                    raise IOError('# statisticalEmittance::correlation: par1 and par2 need to be [0|1|2|3|4|5]')
                    
        def calculateDispersion(self):
            """
            Statistical dispersion evaluation
            Returns: void
            """
            self.dispersionX=self.beamMatrix[0,5]/self.beamMatrix[5,5]
            self.dispersionPx=self.beamMatrix[1,5]/self.beamMatrix[5,5]
            self.dispersionY=self.beamMatrix[2,5]/self.beamMatrix[5,5]
            self.dispersionPy=self.beamMatrix[3,5]/self.beamMatrix[5,5]

        def betatronicMatrices(self):
            """
            Evaluation of the coordinates and beam matrix excluding dispersive components
            Returns: void
            """
            if self.dispersionX is None:
                self.calculateDispersion()

            self.xBetatronic=self.x-self.dispersionX*self.dp
            self.pxBetatronic=self.px-self.dispersionPx*self.dp 
            self.yBetatronic=self.y-self.dispersionY*self.dp 
            self.pyBetatronic=self.py-self.dispersionPy*self.dp

            self.coordinateMatrixBetatronic=np.array([self.xBetatronic,self.pxBetatronic,self.yBetatronic,self.pyBetatronic,self.z,self.dp])
            self.beamMatrixBetatronic=np.matmul(self.coordinateMatrixBetatronic,self.coordinateMatrixBetatronic.T)/len(self.xBetatronic)

        #na checkarw an xreiazetai to full i mono to transverse kommati
        

        self.twissBetax=self.coordinateMatrixBetatronic[0,0]
        
    def centeredCorrelation(self, co1,co2):
        #<(x-<x>)*(y-<y>)> = <x*y> - <x>*<y>
        
    

    def estimateBetaFunction(particles):
        
    def estimateDispersionFunction(particles):
        if(ic < 0 || ic > 1 ) return 0.;
	double x_dE_avg = this->getCorrelation(2*ic, 5);
	double dE2_avg = fabs(this->getCorrelation(5, 5));
	double dispersion = x_dE_avg/dE2_avg * bunch_momentum * bunch_beta;
	return dispersion;
            
    def estimateEmittance(particles):
    
