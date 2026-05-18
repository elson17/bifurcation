cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(PlacentArtAnalysis)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/PlacentArtAnalysis")
set(EXTENSION_CATEGORY "PlacentArtAnalysis")
set(EXTENSION_CONTRIBUTORS "Antoine KNEIB (IADI. CIC-IT)")
set(EXTENSION_DESCRIPTION "Placenta Vasculature Analysis Tool and more ...")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/PlacentArtAnalysis.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/PlacentArtAnalysis/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(BifurcationsAnalysis)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
