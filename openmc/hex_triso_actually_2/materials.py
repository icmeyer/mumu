from math import pi

import openmc

avogadro = 6.0221409e+23

def build_fuel_material(T, pack_frac):
        """ 
        Return openmc material that represents
        a homogenous fuel element

        TRISO specification (density in g/cc, r in microns)
        UCO density = 10.5 r = 175
        C  density = 1.0 r = 275
        C  density = 1.9 r = 315
        SiC density = 3.2 r = 350
        C  density = 1.9 r = 390
        """
        # Working from outer to inner
        # Region 5
        triso_c_n_density   = num_density(350, 390,  1.9,  12,   1)
        # Region 4
        triso_c_n_density  += num_density(315, 350,  3.2,  40,   1)
        triso_si_n_density  = num_density(315, 350,  3.2,  40,   1)
        # Region 3
        triso_c_n_density  += num_density(275, 315,  1.9,  12,   1)
        # Region 2
        triso_c_n_density  += num_density(175, 275,  1.0,  12,   1)
        # Region 1
        triso_c_n_density  += num_density(  0, 175, 10.5, 268, 0.5)
        triso_u_n_density   = num_density(  0, 175, 10.5, 268,   1)
        triso_o_n_density   = num_density(  0, 175, 10.5, 268, 1.5)

        # Scale to a 1 cc volume
        vol_ratio = 1/sphere_vol(0, 390)
        triso_c_n_density   = triso_c_n_density  * vol_ratio
        triso_si_n_density  = triso_si_n_density * vol_ratio
        triso_u_n_density   = triso_u_n_density  * vol_ratio
        triso_o_n_density   = triso_o_n_density  * vol_ratio


        # pack_frac = 0.40
        graphite_n = 1.1995/12 * avogadro
        carbon_density = pack_frac * triso_c_n_density \
                         + (1 - pack_frac)*graphite_n

        homogeneous_fuel =  openmc.Material()
        homogeneous_fuel.set_density('g/cm3', 1.1995)
        homogeneous_fuel.add_element('U', pack_frac*triso_u_n_density, 
                                     enrichment=20)
        homogeneous_fuel.add_element('C', carbon_density)
        homogeneous_fuel.add_element('O', pack_frac*triso_o_n_density)
        homogeneous_fuel.add_element('Si', pack_frac*triso_si_n_density)
        homogeneous_fuel.temperature = T

        return homogeneous_fuel



def sphere_vol(r1, r2):
    return (4/3)*pi*(r2**3 - r1**3)


def num_density(r1, r2, rho, molecule_mass, mols):
    n_density = rho * sphere_vol(r1, r2) * (mols/molecule_mass) * avogadro
    return n_density

    

    



