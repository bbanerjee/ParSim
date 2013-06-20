import numpy as np


#===============================================================================
class MaterialModel:
    # Defines material models - accessed using getStress 
    #  - actual computation done in static methods   
    # Returns stress tensor and jacobian of deformation
    def __init__(self, modelName, props):
        self.modelName = modelName               # Selects Material Model
        self.props = props
        
    def getStress( self, F ):
        model = getattr( self, self.modelName )
        S,Ja = model(self.props, F);    
        return (S,Ja)
    
    def changeProps( self, props ):
        self.props = props


    @staticmethod
    def planeStrainNeoHookean( props, F ):
        # Props - poisson, E
        I2 = F*0.
        I2[0,0] = I2[1,1] = 1.
        v = props['poisson']
        E = props['modulus']
        l = E * v / ((1.+v)*(1.-2.*v))
        m = 0.5 * E / (1.+v)
        Ja = F[0,0]*F[1,1] - F[1,0]*F[0,1]
        S = I2*l*np.log(Ja)/Ja + m/Ja * (np.dot(F, F.T) - I2)
        
        return (S,Ja)
    
    @staticmethod
    def planeStrainNeoHookeanMaxStress( props, F ):
        # Props - poisson, E, maxStress
        I2 = np.eye(2)
        v = props['poisson']
        E = props['modulus']
        sMax = props['maxStress']
        l = E * v / ((1.+v)*(1.-2.*v))
        m = 0.5 * E / (1.+v)
        Ja = np.linalg.det(F)
        S = I2*l*np.log(Ja)/Ja + m/Ja * (np.dot(F, F.T) - I2)
        if vonMises(S) > sMax: 
            S = I2*0.
            Ja = 1.
        
        return (S,Ja)
    
    @staticmethod
    def vonMises( S ):
        return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
            3.*S[1,0]*S[0,1] )