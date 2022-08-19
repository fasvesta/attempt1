# - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * - - * 
"""
#
#   statisticalEmittance
#   A module to calculate the transverse emittances from the beam
#   based on PyORBIT & desy-thesis-05-014
#
#   Version : 1.0
#   Author  : F. Asvesta
#   Contact : fasvesta .at. cern .dot. ch
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
        Input:  x:  [1D numpy array]
                px: [1D numpy array]
                y:  [1D numpy array]
                py: [1D numpy array]
                z:  [1D numpy array]
                dp: [1D numpy array]
                particles coordinates 
        Returns: void
        """   
        self.coordinateMatrix=np.array([x,px,y,py,z,dp])
        self.beamMatrix=np.matmul(self.coordinateMatrix,self.coordinateMatrix.T)/len(self.x)
        self.dispersionX=None
        self.coordinateMatrixBetatronic = None
        self.emittanceX = None

    def correlation(self,par1,par2, betatronic=True):
        """
        Calculation of the correlation for the beam matrices
        Inputs: par1 : [0|1|2|3|4|5]
                par2 : [0|1|2|3|4|5]
                integers corresponding to coordinates (0->x),(1->px), (2->y), (3->py),(4->z), (5->dp)
                betatronic : [bool] if True the betatronic matrices are considered (default=True)
        Returns: <(a-<a>)*(b-<b>)> = <a*b> - <a>*<b>
        """
        if par1 in range(0,6) and par2 in range(0,6):
            if betatronic:
                if par1>3 or par2>3:
                    raise IOError('# statisticalEmittance::correlation: if betatronic par1 and par2 need to be [0|1|2|3]')
                elif self.coordinateMatrixBetatronic is None:
                    self.betatronicMatrices()
                return np.abs(self.beamMatrixBetatronic[par1,par2]-np.nanmean(self.coordinateMatrixBetatronic[par1])*np.nanmean(self.coordinateMatrixBetatronic[par2]))
            else:
                return np.abs(self.beamMatrix[par1,par2]-np.nanmean(self.coordinateMatrix[par1])*np.nanmean(self.coordinateMatrix[par2]))
        else:
            raise IOError('# statisticalEmittance::correlation: par1 and par2 need to be [0|1|2|3|4|5]')

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

        self.coordinateMatrixBetatronic=np.array([xBetatronic,pxBetatronic,yBetatronic,pyBetatronic])
        self.beamMatrixBetatronic=np.matmul(self.coordinateMatrixBetatronic,self.coordinateMatrixBetatronic.T)/len(xBetatronic)
    
    def calculateEmittance(self):
        """
        Transverse emittance evaluation 
        Returns: void
        """
        self.xMatrix=np.array([[self.correlation(0,0, betatronic=True),self.correlation(0,1, betatronic=True)],[self.correlation(1,0, betatronic=True),self.correlation(1,1, betatronic=True)]])
        self.emittanceX=np.sqrt(np.abs(np.linalg.det(self.xMatrix)))
        self.yMatrix=np.array([[self.correlation(2,2, betatronic=True),self.correlation(2,3, betatronic=True)],[self.correlation(3,2, betatronic=True),self.correlation(3,3, betatronic=True)]])
        self.emittanceY=np.sqrt(np.abs(np.linalg.det(self.yMatrix)))
        xYMatrix=np.array([[self.correlation(0,2, betatronic=True),self.correlation(0,3, betatronic=True)],[self.correlation(1,2, betatronic=True),self.correlation(1,3, betatronic=True)]])
        fullMatrix=np.append(np.append(self.xMatrix,xYMatrix,axis=1),np.append(xYMatrix.T,self.yMatrix,axis=1),axis=0)
        self.fourDEmittance=np.sprt(np.linalg.det(fullMatrix))

    def calculateCouplingFactor(self):
        """
        Coupling evaluation as described in doi: 10.18429/JACoW-LINAC2018-THPO118
        Returns: void
        """
        if self.emittanceX is None:
            self.calculateEmittance()
        self.coupling=self.emittanceX*self.emittanceY/self.fourDEmittance - 1.

    
    def calculateTwissFunctions(self):
        """
        Twiss functions evaluation 
        Returns: void
        """
        if self.emittanceX is None:
            self.calculateEmittance()
        self.betaX = self.xMatrix[0,0]/self.emittanceX
        self.alfaX = - self.xMatrix[0,1]/self.emittanceX
        self.gammaX = self.xMatrix[1,1]/self.emittanceX
        self.betaY = self.yMatrix[0,0]/self.emittanceY
        self.alfaY = - self.yMatrix[0,1]/self.emittanceY
        self.gammaY = self.yMatrix[1,1]/self.emittanceY

    def getEmittanceX(self):
        """
        Returns horizontal emittance
        Returns: [float]
        """
        if self.emittanceX is None:
           self.calculateEmittance()
        return self.emittanceX

    def getEmittanceY(self):
        """
        Returns vertical emittance
        Returns: [float]
        """
        if self.emittanceX is None:
           self.calculateEmittance()
        return self.emittanceY

    def getFourDEmittance(self):
        """
        Returns vertical emittance
        Returns: [float]
        """
        if self.emittanceX is None:
           self.calculateEmittance()
        return self.fourDEmittance

    def getNormalizedEmittanceX(self, beta, gamma):
        """
        Returns normalized horizontal emittance
        Inputs: beta:   [float] bunch beta
                gamma:  [float] bunch gamma
        Returns: [float]
        """
        if self.emittanceX is None:
           self.calculateEmittance()
        return self.emittanceX*beta*gamma

    def getNormalizedEmittanceY(self, beta, gamma):
        """
        Returns normalized vertical emittance
        Inputs: beta:   [float] bunch beta
                gamma:  [float] bunch gamma
        Returns: [float]
        """
        if self.emittanceX is None:
           self.calculateEmittance()
        return self.emittanceY*beta*gamma

    def getBetaX(self):
        """
        Returns beta function x
        Returns: [float]
        """
        if self.betaX is None:
           self.calculateTwissFunctions()
        return self.betaX

    def getBetaY(self):
        """
        Returns beta function y
        Returns: [float]
        """
        if self.betaY is None:
           self.calculateTwissFunctions()
        return self.betaY

    def getAlfX(self):
        """
        Returns alfa function x
        Returns: [float]
        """
        if self.betaX is None:
           self.calculateTwissFunctions()
        return self.alfaX

    def getAlfY(self):
        """
        Returns alfa function y
        Returns: [float]
        """
        if self.betaY is None:
           self.calculateTwissFunctions()
        return self.alfaY

    def getGammaX(self):
        """
        Returns alfa function x
        Returns: [float]
        """
        if self.betaX is None:
           self.calculateTwissFunctions()
        return self.gammaX

    def getGammaY(self):
        """
        Returns alfa function y
        Returns: [float]
        """
        if self.betaY is None:
           self.calculateTwissFunctions()
        return self.gammaY

    def getCouplingFactor(self):
        """
        Returns the coupling factor as in doi: 10.18429/JACoW-LINAC2018-THPO118
        coupling factor is 0 for fully uncoupled beams.
        Returns: [float]
        """
        if self.betaY is None:
           self.calculateTwissFunctions()
        return self.gammaY