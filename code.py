import numpy as np
import matplotlib.pyplot as plt
import random as rand
from collections import OrderedDict as Odict
import math as mt

class Modele_vocal_cords:
    
    def __init__(self, lambda_min, lambda_max):
        self.plage = np.linspace(lambda_min, lambda_max, 100)
    
    
    def _compute_coeff(self):
        lambd=10
        
        '''Melanosome'''
        mua_mel = 519* (lambd/500)**(-3) #cm^[-1] valide de [200 - 1200nm]
        
        '''Eau'''
    
    def _init_couches(self):
        self.layers = Odict({})
        
        self.layers['epithelium'] = {'epaisseur': 20*1e-6,
                                     'absorption': 10000,
                                     'diffusion': 1000000,
                                     'g': 0.9}
        self.layers['lamina'] = {'epaisseur': 50*1e-6,
                                     'absorption': 10000,
                                     'diffusion': 1000000,
                                     'g': 0.9}   
                                 
        return self.layers
                                 
class simulation:
    
    def __init__(self, layers):
        self.nb_photons = 10000
        self.thickness = 1*1e-3
        self.delta_z = 1*1e-5
        self.nb_pts = int(self.thickness/self.delta_z)
        self.couche = layers
    

    def transmission_photons(self):

        for i in range(len(self.couche.keys())):
            name = list(self.couche.keys())[i]
            if i == 0:
                self.propagation(self.couche[name])
            else:
                self.propagation(self.couche[name], self.dir)
                
        return self.data
        
    def propagation(self, layer, mu = [0,0,1]):
        absorption = []
        mu_x = mu[0]
        mu_y = mu[1]
        mu_z = mu[2]
        self.absorption_position = {}
        z_photons = 0
        x_photons = 0
        y_photons = 0
        pos_x = [0]
        pos_z = [0]
        
        iteration = 0            
        while iteration < 200:
            iteration +=1
            des_step = rand.random()
            step = -mt.log(des_step)/(layer['absorption']+layer['diffusion'])
            x_photons += mu_x*step
            y_photons += mu_y*step
            z_photons -= mu_z*step
            
            if 0>= pos_z[-1] > -layer['epaisseur']:
            
                absorp = self.nb_photons * (1- mt.exp(-self.delta_z*layer['absorption']))
                self.nb_photons -= absorp
     
                des_theta = rand.random()
                des_phi = rand.random()
                phi = 2* mt.pi*des_phi
                g = layer['g']
            
                anis = (1/(2*g))*(1+g**2-(1-g**2)**2/(1-g+2*g*des_theta)**2)
                theta = mt.acos(anis)
            
            
                if mu_z == 1:
                    mu_x = mt.sin(theta) * mt.cos(phi)
                    mu_y = mt.sin(theta) * mt.cos(phi)
                    mu_z = mt.cos(theta)
                elif mu_z == -1:
                    mu_x = mt.sin(theta) * mt.cos(phi)
                    mu_y = -mt.sin(theta) * mt.cos(phi)
                    mu_z = -mt.cos(theta)
                else:
                    mu_x = mt.sin(theta)*(mu_x*mu_z*mt.cos(theta) - mu_y*mt.sin(phi))/mt.sqrt(1-mu_z**2) + mu_x*mt.cos(theta)
                    mu_y = mt.sin(theta)*(mu_y*mu_z*mt.cos(theta) + mu_x*mt.sin(phi))/mt.sqrt(1-mu_z**2) + mu_y*mt.cos(theta)
                    mu_z = mt.sqrt(1 - mu_z**2)*mt.sin(theta)*mt.cos(phi) + mu_z * mt.cos(theta)
            else:
                pass
                
            pos_z.append(z_photons)
            pos_x.append(x_photons)
            absorption.append(absorp)
            self.dir = [mu_x, mu_y, mu_z]
            self.data = tuple([pos_x, pos_z, absorption])
            

            


class plots:
    
    def __init__(self, data, layers):
        self.thickness = []
        self.data = data
        self.layers = layers
        self.thickness.append(-self.layers['epithelium']['epaisseur'])
                       
    plt.figure()
    def plot_layers(self):
        dat1 = self.data[0]
        dat2 = self.data[1]
        thickness = 0
        for i in range(len(self.layers.keys())):
            name = list(self.layers)[i]
            layer = self.layers[name]
            thickness -= layer['epaisseur']
            plt.plot([min(dat1),max(dat1)],[thickness, thickness])
            
        plt.plot(dat1, dat2)
        plt.plot([min(dat1),max(dat1)],[0, 0])
        plt.xlabel('Position X')
        plt.ylabel('Position Z')
        plt.show()
        
        
mod = Modele_vocal_cords(0, 1)
layers = mod._init_couches()
sim = simulation(layers)
data = sim.transmission_photons()

pl = plots(list(data), layers)
pl.plot_layers()






def read_txt(self,fileID):
    data = []
    separate = [0]
    self.data_x = []
    self.data_y = []
    n=0;
    f = open(fileID,'r')
    f_read = f.read()
    for j,i in enumerate(f_read):
        if i in ["\t", "\n"]:
            n+=1
            separate.append(j)
            data.append(f_read[separate[n-1]:separate[n]])
    
    for i,j in enumerate(data):
        if i%2 == 0:
            j = j[1:]
            data_x.append(j)
        else:
            j = j[1:]
            data_y.append(j)

    data_x = data_x[2:]
    data_y = data_y[2:]
    i=0
    for i in range(10):
        data_x.append(float(data_x[i]))
        data_y.append(float(data_y[i]))




















    
                                     