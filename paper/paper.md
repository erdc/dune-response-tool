---
title: 'DRT: A rapid coastal foredune evolution prediction tool'
tags:
  - matlab
  - coastal
  - morphology
  - dunes
  - storms
authors:
  - name: Nicholas Cohn
    orcid: 0000-0003-4287-039X
    affiliation: 1
  - name: Katherine Brodie
    orcid: 0000-0003-1945-5112
    affiliation: 1
affiliations:
 - name: US Army Engineer Research and Development Center, Coastal and Hydraulics Laboratory - Field Research Facility, Duck, NC
   index: 1
date: 22 June 2021
bibliography: paper.bib
---

# Summary

Coastal foredunes are the first line of defense against ocean-side flooding along many low-lying sandy coastlines. For this reason,
understanding when dunes will be eroded during major storm events and under which conditions these morphological
features will naturally grow from wind-driven processes is important for quantifying long-term coastal risk. Numerous numerical models of varying complexity exist to simulate the accretion and erosion of coastal dunes. All such models require inputs such as wave, wind, and tidal forcing conditions, as well as information regarding the local topographic profile. However, acquiring these necessary site-specific characteristics, setting up numerical models with this information in the appropriate formats,
compiling and running models, and visualizing model output can take significant time and computational effort. These
barriers often limit the ability for non-experts to access or utilize dune modeling tools for planning or prediction purposes.

![DRT code workflow showing connections between the graphical user interface, functions to download and process environmental data, gridding functions, model simulations, and visualization.\label{fig:fig1}](DRT_Workflow.jpg){ width=80% }

The Dune Response Tool ('DRT') was designed to be used by coastal practitioners for planning purposes and
for research applications, with minimal external input required. 'DRT' includes a series of Matlab codes (\autoref{fig:fig1}) to download
relevant oceanographic and meteorological variables from public data sources for any ocean-fronting coastline in the continental US, which are used to force simple models for simulating dune erosion from
waves and for quantifying volumetric dune accretion from winds. These tools can be run in both hindcasting and forecasting modes, with
model simulations typically taking seconds to minutes to run for typical applications. The back-end codes are wrapped within
a graphical user interface (GUI, \autoref{fig:fig2}) created within Matlab 2020b to allow for straightforward field site selection and for visualizing model outputs.

![Screengrab of the DRT graphical user interface.\label{fig:fig2}](DRT_GUI.jpg){ width=60% }

# Statement of Need
Rapid predictive tools of coastal sediment transport and morphology change, even if they lack a full depiction of the relevant
hydrodynamic or morphdynamics processes, are important for understand coastal change hazards. In 'DRT' all
necessary pre-processing steps, model simulations, and visualizations for simple erosional and accretional dune volume change models have been included within a single package.
The codes pull publicly accessible environmental data sources, with these datasources as outlined in \autoref{fig:fig1}, differing for hindcast
and forecast modes based on (1) national scale data availability, (2) data accessibility (e.g., optimizing for fast download speeds and consistent server access),
and (3) data quality. End to end tools that are able to generate appropriate model boundary conditions through to generate output graphics are either cost prohibitive or out of
reach for the non-expert. This tool aims to provide an easy to use means to complete planning level assessments of dune erosion and growth, with datasets pre-loaded or automatically downloaded for the majority of the US coastline.


# 'DRT' Workflow
The code includes pre-processing scripts which automatically download relevant oceanographic and meteorological variables from public data sources for
any ocean-fronting coastline in the continental US in both hindcasting (any date from 1980 to 2017) and
forecasting (now to 4 days from now) mode as shown in \autoref{fig:fig1}.

Simplified beach-dune profiles can be generated from a database of morphometric data hat covers the majority of the U.S. oceanfront coastline that was derived from datasets from @USGS:2017 for the Gulf and Atlantic US coastlines and @Mull:2014 for the US Pacific Northwest coastline. Data is stored in an excel file included with the code that can be modified to add new morphometric data in case of updated data availability in the future. Where pre-compiled beach and dune profile data does not exist or where the user wishes to input their own morphology information, this information can be input into the tool. An additional radio button is added to call the ArcGIS server for the Joint Airborne Lidar Bathymetry Technical Center of Expertise (JALBTCX) which provides dune crest and dune toe elevation data for portions of the US coastline as well. A gridded profile based on the user or database provided morphometrics is then generated  for input into the volume change models.

The downloaded environmental forcings are subsequently used to calculate volumetric dune erosion for a user-defined time interval at hourly time steps utilizing a 1D dune erosion model based on @Palmsten:2012. A set of equations for estimating wind-blown transport into dunes following the general approach of @Delgado:2011, implementing the transport equation of @Kawamura:1951, is also included.

The tool purposely has relatively few input parameters, although the user may modify relevant wind- and wave- related transport coefficients, as well as local grain size, wave runup modifiers, and morphologic parameters in the "Model Attributes" or "Advanced" tab within the GUI. A simple ensemble approach is added where a range of input parameters can be included to bound the solution and give a better sense of uncertainty of predictions for the end-user.

Although these various tools do not include all relevant morphodynamic processes controlled dune evolution during storms, the combination of approaches is purposely meant to provide a rapid answer for planning purposes. For more complicated cases, the code can be modified to ingest different data inputs or allow for new capabilities that are not yet incorporated.


# 'DRT 'Example
An example of a model setup initialized on 29 October 2012 at 39.42°N, 74.32°W near Brigantine, NJ represents the model predicted erosion and accretion from Hurricane Sandy, as shown in \autoref{fig:fig3}. Consistent with observations from the storm [@Sopkin:2014], the tool shows destruction of the dune system over the time period. On a standard PC the tool took less than a minute to download relevant environmental data from the web, run the model, and display the output.

![Example model hindcast model output for Hurricane Sandy for a site in Brigantine, New Jersey, USA. Time series of wind, wave, and water levels are provided in the upper panels. Simulated erosional profile change is shown in the bottom left panel. Volumetric wave-driven dune erosion (blue line) and wind-driven volumetric dune growth (red line) are shown in the bottom right panel. \label{fig:fig3}](DRT_Example.jpg)

# Installation
The code consists of a matlab-based .mlapp interface which calls a series of individual .m scripts, as described in \autoref{fig:fig1}. These scripts can also be run independently outside of the GUI.

# Future Work
This initial tool development utilizes available datasets for a subset of the continental US. Future iterations could include morphologic datasets that also cover California, Hawaii, Alaska, and the Great Lakes region, as well as overseas regions. This initial domain was based on ease of access to data availability.

This tool relies of web access to various environmental datasets as well. Most of the these datasets reside on either

 THREDDS or ERDDAP servers, with the later platform generally being more stable. Based on model testing there is not 100% uptime in all of these data products which provides complications for running DRT if servers are down. Future iterations of the tool could add in new web services and/or allow for easier ingestion of non-web-based environmental datasets.

# Acknowledgements
Development of 'DRT' was funded by the U.S. Army Engineer Research and Development Center's (ERDC) Coastal Inlets Research Program (CIRP).

# References
