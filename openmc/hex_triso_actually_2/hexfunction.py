import openmc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from math import pi
import copy

from materials import build_fuel_material

def hexfunction(pitch, pack_frac): 

    settings = openmc.Settings()
    # Set high tolerance to allow use of lower temperature xs
    settings.temperature['tolerance']=1000
    settings.temperature['method'] = 'nearest'
    settings.temperature['multipole'] = True
    settings.cutoff={'energy': 1e-8} #energy cutoff in eV
    
    #############################
    ###       MATERIALS       ###
    #############################
    mat_list = []
    enrichment = 20.0

    uo2 = openmc.Material(1, "uo2")
    uo2.add_element('U', 1.0, enrichment=enrichment)
    uo2.add_element('O', 2.0)
    uo2.set_density('g/cm3', 10.97)
    uo2.temperature = 900 #kelvin
    mat_list.append(uo2)
    
    graphite = openmc.Material(2, "graphite")
    graphite.set_density('g/cm3', 1.1995)
    graphite.add_element('C', 1.0)
    graphite.add_s_alpha_beta('c_Graphite')
    graphite.temperature = 900 #kelvin
    mat_list.append(graphite)

    # sodium = openmc.Material(3, "sodium")
    sodium = openmc.Material()
    sodium.set_density('g/cm3', 0.8017) # 900 K
    sodium.add_element('Na', 1.0)
    # sodium.add_s_alpha_beta('c_Graphite')
    sodium.temperature = 900 #kelvin
    mat_list.append(sodium)

    # naoh = openmc.Material(6, "naoh")
    naoh = openmc.Material()
    naoh.set_density('g/cm3', 1.5) # 900 K
    naoh.add_element('Na', 1.0)
    naoh.add_element('O', 1.0)
    naoh.add_element('H', 1.0)
    # sodium.add_s_alpha_beta('c_Graphite')
    naoh.temperature = 900 #kelvin
    mat_list.append(naoh)

    # TRISO Materials
    fuel = openmc.Material(name='Fuel')
    fuel.set_density('g/cm3', 10.5)
    # fuel.add_nuclide('U235', 4.6716e-02)
    fuel.add_nuclide('U235', 0.0667372)
    # fuel.add_nuclide('U238', 2.8697e-01)
    fuel.add_nuclide('U238', 0.2669488)
    fuel.add_nuclide('O16',  5.0000e-01)
    fuel.add_element('C', 1.6667e-01)
    mat_list.append(fuel)

    buff = openmc.Material(name='Buffer')
    buff.set_density('g/cm3', 1.0)
    buff.add_element('C', 1.0)
    buff.add_s_alpha_beta('c_Graphite')
    mat_list.append(buff)

    PyC1 = openmc.Material(name='PyC1')
    PyC1.set_density('g/cm3', 1.9)
    PyC1.add_element('C', 1.0)
    PyC1.add_s_alpha_beta('c_Graphite')
    mat_list.append(PyC1)

    PyC2 = openmc.Material(name='PyC2')
    PyC2.set_density('g/cm3', 1.87)
    PyC2.add_element('C', 1.0)
    PyC2.add_s_alpha_beta('c_Graphite')
    mat_list.append(PyC2)

    SiC = openmc.Material(name='SiC')
    SiC.set_density('g/cm3', 3.2)
    SiC.add_element('C', 0.5)
    SiC.add_element('Si', 0.5)
    mat_list.append(SiC)

    fuel_temp = 900
    homogeneous_fuel = build_fuel_material(fuel_temp, pack_frac)
    mat_list.append(homogeneous_fuel)
    
    mats = openmc.Materials(mat_list)
    mats.export_to_xml()
    
    #############################
    ###       GEOMETRY        ###
    #############################
    pitch = 17.4
    fuel_bottom = -pitch/2
    fuel_top = pitch/2
    coolant_r = 4
    # fuel_r = (pitch/2 - coolant_r)/2
    fuel_r = (pitch/(3**(1/2)) - coolant_r)/2

    hex_universe = openmc.Universe()

    top = openmc.ZPlane(z0=fuel_top, boundary_type='reflective')
    bottom = openmc.ZPlane(z0=fuel_bottom, boundary_type='reflective')
    surf_fuel = openmc.ZCylinder(r=fuel_r)

    # Make TRISOS to be filled in fuel cylinders by chopping up
    # fuel cylinder into segments
    n_cyls = 40
    fuel_segment_heights = np.linspace(fuel_bottom, fuel_top, n_cyls)
    segment_height = fuel_segment_heights[1] - fuel_segment_heights[0]
    fuel_planes = [bottom]
    fuel_cells = []
    for i, height in enumerate(fuel_segment_heights[1:-1]):
        this_plane = openmc.ZPlane(z0=height)
        fuel_planes.append(this_plane)
        this_cell = openmc.Cell()
        this_cell.region = +fuel_planes[i] & -fuel_planes[i+1] & -surf_fuel
        fuel_cells.append(copy.deepcopy(this_cell))
    # last cell
    fuel_planes.append(top)
    this_cell = openmc.Cell()
    this_cell.region = +fuel_planes[-2] & -fuel_planes[-1] & -surf_fuel
    fuel_cells.append(copy.deepcopy(this_cell))

    # Make fuel cylinder
    fuel_cyl_top = openmc.ZPlane(z0=segment_height/2)
    fuel_cyl_bottom = openmc.ZPlane(z0=-segment_height/2)
    fuel_triso_region = -surf_fuel & +fuel_cyl_bottom & -fuel_cyl_top
    outer_radius = 425.*1e-4
    # openmc.model.triso._Cylinder.from_region(fuel_region, outer_radius)

    spheres = [openmc.Sphere(r=r*1e-4) for r in [215., 315., 350., 385.]]
    cells = [openmc.Cell(fill=fuel, region=-spheres[0]),
             openmc.Cell(fill=buff, region=+spheres[0] & -spheres[1]),
             openmc.Cell(fill=PyC1, region=+spheres[1] & -spheres[2]),
             openmc.Cell(fill=SiC, region=+spheres[2] & -spheres[3]),
             openmc.Cell(fill=PyC2, region=+spheres[3])]
    triso_univ = openmc.Universe(cells=cells)
    outer_radius = 425.*1e-4
    centers = openmc.model.pack_spheres(radius=outer_radius,
                                        region=fuel_triso_region,
                                        pf=pack_frac)
    trisos = [openmc.model.TRISO(outer_radius, triso_univ, c) for c in centers]
    outside_trisos = openmc.Intersection(~t.region for t in trisos)
    # background_region = outside_trisos & +fuel_cyl_bottom & \
                        # -fuel_cyl_top & -surf_fuel
    background_region = outside_trisos
    background_cell = openmc.Cell(fill=graphite, region=background_region)

    fuel_triso_univ = openmc.Universe()
    fuel_triso_univ.add_cell(background_cell)
    for idx, triso in enumerate(trisos):
        fuel_triso_univ.add_cell(triso)

    # Fill in fuel cells with triso cells and translate to location
    for i, cell in enumerate(fuel_cells):
        cell_height = segment_height*(i+1/2) + fuel_bottom
        cell.translation = [0, 0, cell_height]
        cell.fill = fuel_triso_univ
    fuel_cell_univ = openmc.Universe(cells=fuel_cells)

    # For testing solid fuel
    # test_region = +bottom & -top & -surf_fuel
    # fuel_cell = openmc.Cell(region=test_region, fill=fuel)
    # fuel_cell_univ = openmc.Universe(cells=[fuel_cell])

    coolant_cyl = openmc.ZCylinder(r=coolant_r)
    coolant_region = -coolant_cyl
    coolant_cell = openmc.Cell()
    coolant_cell.fill = naoh
    coolant_cell.region = coolant_region
    hex_universe.add_cell(coolant_cell)

    hex_prism = openmc.get_hexagonal_prism(edge_length=pitch/(3**1/2),
                                             boundary_type='reflective')

    graphite_region = hex_prism & +coolant_cyl & -top & +bottom
    graphite_cell = openmc.Cell()
    graphite_cell.fill = graphite
    graphite_cell.region = graphite_region
    hex_universe.add_cell(graphite_cell)

    fuel_cells = []
    root3 = 3**(1/2)
    half_to_vertex =  pitch/root3/2
    half_to_edge = pitch/4

    # fuel_id = 100
    offset_angle = 30
    n_pins = 6
    for i in range(n_pins):
        theta = (offset_angle + i/n_pins*360)*pi/180
        r = coolant_r + fuel_r + 0.01
        x = r*np.cos(theta)
        y = r*np.sin(theta)

        fuel_cyl_bound = openmc.ZCylinder(x0=x, y0=y, r=fuel_r)
        graphite_cell.region &= +fuel_cyl_bound

        fuel_cell = openmc.Cell()
        fuel_cell.fill = copy.deepcopy(fuel_cell_univ)
        fuel_cell.translation = [x, y, 0]
        fuel_cell.region = -fuel_cyl_bound & -top & +bottom
        # fuel_cell.id = fuel_id
        # fuel_id += 1

        fuel_cells.append(fuel_cell)
        hex_universe.add_cell(fuel_cell)
    
    geom = openmc.Geometry(hex_universe)
    # geom = openmc.Geometry(fuel_cell_univ)
    geom.export_to_xml()
    
    #####################################
    ###        SOURCE/BATCHES         ###
    #####################################
    point = openmc.stats.Point((0, 0, 0))
    src = openmc.Source(space=point)
    
    settings.source = src
    settings.batches = 50
    settings.inactive = 10
    settings.particles = 200
    
    settings.export_to_xml()
    
    
    #############################
    ###       TALLIES         ###
    #############################
    # Instantiate an empty Tallies object
    tallies_file = openmc.Tallies()
    
    # K-Eigenvalue (infinity) tallies
    fiss_rate = openmc.Tally(name='fiss. rate')
    fiss_rate.scores = ['nu-fission']
    tallies_file.append(fiss_rate)
    
    abs_rate = openmc.Tally(name='abs. rate')
    abs_rate.scores = ['absorption']
    tallies_file.append(abs_rate)
    
    # Resonance Escape Probability tallies
    therm_abs_rate = openmc.Tally(name='therm. abs. rate')
    therm_abs_rate.scores = ['absorption']
    therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625])]
    tallies_file.append(therm_abs_rate)
    
    # Thermal Flux Utilization tallies
    # fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
    # fuel_therm_abs_rate.scores = ['absorption']
    # fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625]),
    #                                        openmc.CellFilter([fuel_cell])]
    # tallies_file.append(fuel_therm_abs_rate)
    
    # Fast Fission Factor tallies
    therm_fiss_rate = openmc.Tally(name='therm. fiss. rate')
    therm_fiss_rate.scores = ['nu-fission']
    therm_fiss_rate.filters = [openmc.EnergyFilter([0., 0.625])]
    tallies_file.append(therm_fiss_rate)
    
    tallies_file.export_to_xml()
    
    #############################
    ###       PLOTTING        ###
    #############################
    zs = np.linspace(0,1,2)
    plots = []
    for z in zs:
        p = openmc.Plot()
        p.filename = 'pinplot'+str(z)
        p.width = (1.4*pitch, 1.4*pitch)
        p.pixels = (2000, 2000)
        p.color_by = 'material'
        p.origin = [0, 0, z]
        # p.color_by = 'cell'
        # p.colors = {homogeneous_fuel: 'yellow', naoh: 'grey', graphite: 'black'}
        p.colors = {fuel: 'yellow', naoh: 'grey', graphite: 'black'}
        plots.append(copy.deepcopy(p))
    
    plots = openmc.Plots(plots)
    plots.export_to_xml()
    
    # openmc.plot_geometry(output = False)
    openmc.plot_geometry()
    # pngstring = 'pinplot{}.png'.format(str(pitch))
    # subprocess.call(['convert','pinplot.ppm',pngstring])
    # subprocess.call(['mv',pngstring,'figures/'+pngstring])
    
    
    #############################
    ###       EXECUTION       ###
    #############################
    # openmc.run(output=False)
    openmc.run()
    sp = openmc.StatePoint('statepoint.{}.h5'.format(settings.batches))
    # Collect all the tallies 
    fiss_rate = sp.get_tally(name='fiss. rate')
    fiss_rate_df = fiss_rate.get_pandas_dataframe()
    abs_rate = sp.get_tally(name='abs. rate')
    abs_rate_df = abs_rate.get_pandas_dataframe()
    therm_abs_rate = sp.get_tally(name='therm. abs. rate')
    therm_abs_rate_df = therm_abs_rate.get_pandas_dataframe()
    fuel_therm_abs_rate = sp.get_tally(name='fuel therm. abs. rate')
    fuel_therm_abs_rate_df = fuel_therm_abs_rate.get_pandas_dataframe()
    therm_fiss_rate = sp.get_tally(name='therm. fiss. rate')
    therm_fiss_rate_df = therm_fiss_rate.get_pandas_dataframe()
    
    
    # Compute k-infinity
    kinf = fiss_rate / abs_rate
    kinf_df = kinf.get_pandas_dataframe()
    
    # Compute resonance escape probability
    res_esc = (therm_abs_rate)/(abs_rate)
    res_esc_df = res_esc.get_pandas_dataframe()
    
    # Compute fast fission factor
    fast_fiss = fiss_rate/therm_fiss_rate
    fast_fiss_df = fast_fiss.get_pandas_dataframe()
    
    # Compute thermal flux utilization
    therm_util = fuel_therm_abs_rate / therm_abs_rate
    therm_util_df = therm_util.get_pandas_dataframe()
    
    # Compute neutrons produced per absorption
    eta = therm_fiss_rate / fuel_therm_abs_rate
    eta_df = eta.get_pandas_dataframe()
    

    columns = ['pitch','enrichment',
             'kinf mean','kinf sd','res_esc mean','res_esc sd',
             'fast_fiss mean','fast_fiss sd','therm_util mean',
             'therm_util sd','eta mean','eta sd']
    data = [[pitch,enrichment,
            kinf_df['mean'][0],kinf_df['std. dev.'][0],
            res_esc_df['mean'][0],res_esc_df['std. dev.'][0],
            fast_fiss_df['mean'][0],fast_fiss_df['std. dev.'][0],
            therm_util_df['mean'][0],therm_util_df['std. dev.'][0],
            eta_df['mean'][0],eta_df['std. dev.'][0]]]
    all_tallies = pd.DataFrame(data, columns = columns)
    
    return all_tallies

    
    # How to extract values from dataframes
    # print(kinf_df['mean'])
    # print(kinf_df['std. dev.'])


if __name__ == '__main__':
    pitch = 17.4
    pack_frac = 0.10
    all_tallies = hexfunction(pitch, pack_frac)
