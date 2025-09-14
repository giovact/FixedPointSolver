
using PyPlot, LaTeXStrings, PyCall, JLD


@pyimport matplotlib.gridspec as gspec
@pyimport matplotlib.cm as cms
@pyimport matplotlib.colors as pltcolors

#include("plot_free_energies.jl")
include("plot_res_scanning.jl")
include("plot_res_exploration.jl")