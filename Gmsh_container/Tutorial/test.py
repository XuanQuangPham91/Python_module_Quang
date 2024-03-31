# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 1
#
#  Geometry basics, elementary entities, physical groups
#
# ------------------------------------------------------------------------------

# The Python API is entirely defined in the `gmsh.py' module (which contains the
# full documentation of all the functions in the API):
import gmsh
import sys

# factory = gmsh.model.geo
factory = gmsh.model.occ

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()

# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("t1")

# The Python API provides direct access to each supported geometry (CAD)
# kernel. The built-in kernel is used in this first tutorial: the corresponding
# API functions have the `gmsh.model.geo' prefix.

# The first type of `elementary entity' in Gmsh is a `Point'. To create a point
# with the built-in CAD kernel, the Python API function is
# gmsh.model.geo.addPoint():
# - the first 3 arguments are the point coordinates (x, y, z)
# - the next (optional) argument is the target mesh size close to the point
# - the last (optional) argument is the point tag (a stricly positive integer
#   that uniquely identifies the point)
lc = 1e-2
r_ref = 0.34
r_out = 0.48
# center
# Points =======================================================================
factory.addPoint(0, 0, 0, lc, 1)
## Rectangular points
factory.addPoint(-0.5, -0.5, 0, lc, 2)
factory.addPoint(-0.5, +0.5, 0, lc, 3)
factory.addPoint(+0.5, +0.5, 0, lc, 4)
factory.addPoint(+0.5, -0.5, 0, lc, 5)

## Circle points
# factory.addPoint(-r_ref, 0.0, 0, lc, 6)
# factory.addPoint(-r_out, 0.0, 0, lc, 7)

# Lines ========================================================================
factory.addLine(2, 3, 1)
factory.addLine(3, 4, 2)
factory.addLine(4, 5, 3)
factory.addLine(5, 2, 4)

#===============================================================================
# gmsh.model.geo.addCircleArc(startTag=, centerTag=, endTag=, tag=-1, nx=0., ny=0., nz=0.)
# Volumes can be constructed from (closed) curve loops thanks to the
# `addThruSections()' function
factory.addCircle(0, 0, 0, r_out, 6)
factory.addCurveLoop([6], 6)
factory.addCircle(0, 0, 0, r_ref, 7)
factory.addCurveLoop([7], 7)
factory.synchronize()


# The third elementary entity is the surface. In order to define a simple
# rectangular surface from the four curves defined above, a curve loop has first
# to be defined. A curve loop is defined by an ordered list of connected curves,
# a sign being associated with each curve (depending on the orientation of the
# curve to form a loop). The API function to create curve loops takes a list
# of integers as first argument, and the curve loop tag (which must be unique
# amongst curve loops) as the second (optional) argument:
factory.addCurveLoop([1, 2, 3, 4], 1)


# We can then define the surface as a list of curve loops (only one here,
# representing the external contour, since there are no holes--see `t4.py' for
# an example of a surface with a hole):
factory.addPlaneSurface([1], 1)
factory.addPlaneSurface([6], 6)
factory.addPlaneSurface([7], 7)

# Before they can be meshed (and, more generally, before they can be used by API
# functions outside of the built-in CAD kernel functions), the CAD entities must
# be synchronized with the Gmsh model, which will create the relevant Gmsh data
# structures. This is achieved by the gmsh.model.geo.synchronize() API call for
# the built-in CAD kernel. Synchronizations can be called at any time, but they
# involve a non trivial amount of processing; so while you could synchronize the
# internal CAD data after every CAD command, it is usually better to minimize
# the number of synchronization points.
factory.synchronize()

# At this level, Gmsh knows everything to display the rectangular surface 1 and
# to mesh it. An optional step is needed if we want to group elementary
# geometrical entities into more meaningful groups, e.g. to define some
# mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
# material ("steel", "carbon") properties.
#
# Such groups are called "Physical Groups" in Gmsh. By default, if physical
# groups are defined, Gmsh will export in output files only mesh elements that
# belong to at least one physical group. (To force Gmsh to save all elements,
# whether they belong to physical groups or not, set the `Mesh.SaveAll' option
# to 1.) Physical groups are also identified by tags, i.e. stricly positive
# integers, that should be unique per dimension (0D, 1D, 2D or 3D). Physical
# groups can also be given names.
#
# Here we define a physical curve that groups the left, bottom and right curves
# in a single group (with prescribed tag 5); and a physical surface with name
# "My surface" (with an automatic tag) containing the geometrical surface 1:
gmsh.model.addPhysicalGroup(1, [1, 2, 4], 5)
gmsh.model.addPhysicalGroup(2, [1], name="My surface")

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

# ... and save it to disk
gmsh.write("t1.msh")

# Remember that by default, if physical groups are defined, Gmsh will export in
# the output mesh file only those elements that belong to at least one physical
# group. To force Gmsh to save all elements, you can use
#
# gmsh.option.setNumber("Mesh.SaveAll", 1)

# By default, Gmsh saves meshes in the latest version of the Gmsh mesh file
# format (the `MSH' format). You can save meshes in other mesh formats by
# specifying a filename with a different extension. For example
#
#   gmsh.write("t1.unv")
#
# will save the mesh in the UNV format. You can also save the mesh in older
# versions of the MSH format: simply set
#
#   gmsh.option.setNumber("Mesh.MshFileVersion", x)
#
# for any version number `x'. As an alternative, you can also not specify the
# format explicitly, and just choose a filename with the `.msh2' or `.msh4'
# extension.

# To visualize the model we can run the graphical user interface with
# `gmsh.fltk.run()'. Here we run it only if "-nopopup" is not provided in the
# command line arguments:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Note that starting with Gmsh 3.0, models can be built using other geometry
# kernels than the default "built-in" kernel. To use the OpenCASCADE CAD kernel
# instead of the built-in kernel, you should use the functions with the
# `gmsh.model.occ' prefix.
#
# Different CAD kernels have different features. With OpenCASCADE, instead of
# defining the surface by successively defining 4 points, 4 curves and 1 curve
# loop, one can define the rectangular surface directly with
#
# gmsh.model.occ.addRectangle(.2, 0, 0, .1, .3)
#
# After synchronization with the Gmsh model with
#
# gmsh.model.occ.synchronize()
#
# the underlying curves and points could be accessed with
# gmsh.model.getBoundary().
#
# See e.g. `t16.py', `t18.py', `t19.py' or `t20.py' for complete examples based
# on OpenCASCADE, and `examples/api' for more.

# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()
