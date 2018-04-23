# script to get an example of everything provided in the postMINX program
# can run with cropped versions of plumes as well # note to self to set option for cropping the plume area

# import things that are necessary
# import postMINX as pm

# get list of all MINX outputs wanted for a given orbit
# probably something like
# txts = glob.glob('/fs/site2/dev/eccc/aq/r1/nod001/MINX/Working/AB_SK_*2017/Digitize_Output/*/Plumes_O*-B*-SPWB*.txt')
minxPlumePath = ''
minxPlumes = [] #list of testfiles stored in a location
modelPath = '' #path to all model files, postMINX will deal with picking which file to use hopefully

# TODO: get all the relevant datetime firepixels all in one location so can process giant chunks at once

# Get all plume data
# Get all plume and associated model data (maybe change to opening/closing fst only once?), note that model info can be added in at the same time as initializating the plume
plumes = []
# this should get all key features for the plume to do things such as stats and plots
for ind, plumeText in enumerate(minxPlumes):
    p = pm.Plume.from_filename(plumeText, model_dir = modelPath)
    # should by now have lon/lat/distance/direction/terrain/nowindheight/windheight/filteredheight/wind vars/modelheight/val/points
    # get remaining content desired (products)
    # example
    p.Model.get_product(['GZ', 'AF', 'TO3'])
    plumes.append(p)

# For each plume/model plume
# Plot and save the following comparisons
# NOTE: postMINX has changed so that user has to ask for show and whatnot, just save all of them and close all as needed

#  want a function to call plume function and then save
methods = ['contour', 'colour', 'scatter_3D', 'hist']
for plm in plumes:


# Can have a few plots on each image for comparison (say by plumes on the same day if there's more than one plume a day)

# Save all statistics to a table then export to a csv
methods = ['max', 'min', 'median', 'mean', 'std']
# additional functions for a cropped plume. Will probably have to be more selective about the plumes picked for this
