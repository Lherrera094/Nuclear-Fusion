#!/usr/bin/env python
# coding: utf-8
# %%


# Developer: Luis Carlos Herrera Quesada
# Date: 23/05/2023
# Universidad Carlos III de Madrid


# ## Surface Representations from VMEC file ("wout")

# Necessary libraries

# %%


import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
import plotly.graph_objects as go
from matplotlib import cm
import os


# %%


#Creates file to save matrix values
def File_Data_Save(path):
    file = "magnetic_coordinates.txt"
    
    if file not in path:
        new_file = open(file, "w")
        
    return new_file


# %%


#Creates folder to save images and data
def Folder_Save(device):
    
    if not os.path.exists(f"{device}/"):
        os.makedirs(f"{device}")


# Extract main data: 
# <br>
# * rmnc = Fourier coefficients for the radial values
# * zmnc = Fourier coefficients for the z values
# * r1 = Major radius
# * r2 = Minor radius
# * itor = toroidal modes
# * ipol = poloidal modes
# * n_surf = number of surfaces
# * Nf = number of field periods for device

# %%


def get_Values(vmec_file):
    
    # Extract data
    rmnc = vmec_file.variables['rmnc'][:]
    zmns = vmec_file.variables['zmns'][:]
    lmns = vmec_file.variables['lmns'][:] #lambda
    bmnc = vmec_file.variables['bmnc'][:]

    r1 = vmec_file.variables['Rmajor_p'].getValue() 
    r2 = vmec_file.variables['Aminor_p'].getValue() 

    itor = vmec_file.variables['xn'][:]  
    ipol = vmec_file.variables['xm'][:]  
    n_surf = vmec_file.variables['ns'].getValue()  
    Nf = vmec_file.variables['nfp'].getValue() 
    n_imode = len(itor)
    
    return rmnc, zmns, bmnc, itor, ipol, n_surf, n_imode


# %%


def mesh_Toroidal_Coordinates(res):
    # Construct flux coordinates grid
    phi = np.linspace(0, 2*np.pi, res)
    theta = np.linspace(0, 2*np.pi, res)
    Theta, Phi = np.meshgrid(theta, phi)
    
    return Theta, Phi


# # Magnetic Surface

# The next function creates a grid of ($\varphi$,$\vartheta$) coordinates, then creates empty matrices for coordinates R and Z and finally computes the magnetic surface data

# %%


def Magnetic_Surface_Matrix(rmnc, zmns, bmnc, ipol, itor, n_imode, res, isurf, n_surf):
    
    Theta, Phi = mesh_Toroidal_Coordinates(res)
    #defines the matrix for (R,Z)
    R,Z,B = np.empty((n_surf,res,res)), np.empty((n_surf,res,res)), np.empty((n_surf,res,res))
    
    for i_phi in range(res): 
        for i_theta in range(res):
            R[isurf,i_phi,i_theta] = 0
            Z[isurf,i_phi,i_theta] = 0
            B[isurf,i_phi,i_theta] = 0
            for imode in range(n_imode):
                R[isurf,i_phi,i_theta] += rmnc[isurf,imode] * np.cos( (ipol[imode] * Theta[i_phi,i_theta]) - (itor[imode] * Phi[i_phi,i_theta]) )
                Z[isurf,i_phi,i_theta] += zmns[isurf,imode] * np.sin( (ipol[imode] * Theta[i_phi,i_theta]) - (itor[imode] * Phi[i_phi,i_theta]) )
                B[isurf,i_phi,i_theta] += bmnc[isurf,imode] * np.cos( (ipol[imode] * Theta[i_phi,i_theta]) - (itor[imode] * Phi[i_phi,i_theta]) )
            
            
    return R,Z,B


# Next function obtains the cartesian coordinates. Creates empty matrices for X and Y and transforms R and Z to cartesian

# %%


def Cartesian_Coordinates_Change(r,z,res,n_surf):
    Theta, Phi = mesh_Toroidal_Coordinates(res)
    X,Y = np.empty((n_surf,res,res)), np.empty((n_surf,res,res))
    X = r * np.cos(Phi)
    Y = r * np.sin(Phi)
    
    return X,Y


# #When paralelization is available.
# <br>
# for isurf in range(n_surf):
#     for i_phi in range(len(phi)): 
#         for i_theta in range(len(theta)):
#             R[isurf,i_phi,i_theta] = 0
#             Z[isurf,i_phi,i_theta] = 0
#             for imode in range(n_imode):
#                 R[isurf,i_phi,i_theta] += rmnc[isurf,imode] * np.cos( (ipol[imode] * Theta[i_phi,i_theta]) - (itor[imode] * Phi[i_phi,i_theta]) )
#                 Z[isurf,i_phi,i_theta] += zmns[isurf,imode] * np.sin( (ipol[imode] * Theta[i_phi,i_theta]) - (itor[imode] * Phi[i_phi,i_theta]) )

# # Poloidal proyection in given toroidal angle

# %%


def Poloidal_Proyection(rmnc, zmns, ipol, itor, n_imode, res, i_phi, n_surf):
    
    Theta, Phi = mesh_Toroidal_Coordinates(res)
    
    #defines the matris for (r,z)
    R,Z = np.empty((n_surf,res,res)), np.empty(((n_surf,res,res)))
    
    for isurf in range(n_surf): 
        for i_theta in range(res):
            R[isurf,i_phi,i_theta] = 0
            Z[isurf,i_phi,i_theta] = 0
            for imode in range(n_imode):
                R[isurf,i_phi,i_theta] += rmnc[isurf,imode] * np.cos( (ipol[imode] * Theta[i_phi,i_theta]) - (itor[imode] * Phi[i_phi,i_theta]) )
                Z[isurf,i_phi,i_theta] += zmns[isurf,imode] * np.sin( (ipol[imode] * Theta[i_phi,i_theta]) - (itor[imode] * Phi[i_phi,i_theta]) )
                
    return R,Z


# %%




