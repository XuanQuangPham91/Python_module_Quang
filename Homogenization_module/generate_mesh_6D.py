from dolfin import *
from mshr import *
import matplotlib.pyplot as plt

# Define domain
d_out = 0.5
r_ref = 0.375
outer_rectangle = Rectangle(Point(-d_out, -d_out), Point(d_out, d_out))
inner_rectangle = Rectangle(Point(-r_ref, -r_ref), Point(r_ref, r_ref))
sub_3 = Rectangle(Point(-d_out, -d_out), Point(-r_ref, r_ref))
sub_4 = Rectangle(Point(-d_out, r_ref), Point(r_ref, d_out))
sub_5 = Rectangle(Point(r_ref, -r_ref), Point(d_out, d_out))
sub_6 = Rectangle(Point(-r_ref, -d_out), Point(d_out, -r_ref))
circle = Circle(Point(0., 0.), r_ref, segments=30)
domain = outer_rectangle

# Create mesh
domain.set_subdomain(1, circle)
domain.set_subdomain(2, inner_rectangle-circle)
domain.set_subdomain(3, sub_3)
domain.set_subdomain(4, sub_4)
domain.set_subdomain(5, sub_5)
domain.set_subdomain(6, sub_6)

# interver_num = 78
interver_num = 30
mesh = generate_mesh(domain, interver_num)
print(mesh.num_vertices())

# Create subdomains
subdomains = MeshFunction("size_t", mesh, 2, mesh.domains())
plt.jet()
plt.colorbar(plot(subdomains))
plt.show()

# Create boundaries
class Left(SubDomain):
    def __init__(self, y_min, y_max):
       SubDomain.__init__(self)
       self.y_min = y_min
       self.y_max = y_max

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], -d_out) and x[1] >= self.y_min and x[1] <= self.y_max


class Right(SubDomain):
    def __init__(self, y_min, y_max):
        SubDomain.__init__(self)
        self.y_min = y_min
        self.y_max = y_max

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], d_out) and x[1] >= self.y_min and x[1] <= self.y_max


class Bottom(SubDomain):
    def __init__(self, x_min, x_max):
        SubDomain.__init__(self)
        self.x_min = x_min
        self.x_max = x_max

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], -d_out) and x[0] >= self.x_min and x[0] <= self.x_max
        # return on_boundary and near(x[1], -1.0)


class Top(SubDomain):
    def __init__(self, x_min, x_max):
        SubDomain.__init__(self)
        self.x_min = x_min
        self.x_max = x_max

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], d_out) and x[0] >= self.x_min and x[0] <= self.x_max


boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)
bottomRight = Bottom(r_ref, d_out)
bottomRight.mark(boundaries, 1)
bottomMid = Bottom(-r_ref,r_ref)
bottomMid.mark(boundaries, 2)
bottomLeft = Bottom(-d_out, -r_ref)
bottomLeft.mark(boundaries, 3)

leftBot = Left(-d_out, -r_ref)
leftBot.mark(boundaries, 4)
leftMid = Left(-r_ref,r_ref)
leftMid.mark(boundaries, 5)
leftTop = Left(r_ref, d_out)
leftTop.mark(boundaries, 6)

topLeft = Top(-d_out, -r_ref)
topLeft.mark(boundaries, 7)
topMid = Top(-r_ref,r_ref)
topMid.mark(boundaries, 8)
topRight = Top(r_ref, d_out)
topRight.mark(boundaries, 9)

rightTop = Right(r_ref, d_out)
rightTop.mark(boundaries, 10)
rightMid = Right(-r_ref,r_ref)
rightMid.mark(boundaries, 11)
rightBot = Right(-d_out, -r_ref)
rightBot.mark(boundaries, 12)

# Save
File("mesh_6B/thermal_block_R%s.xml" % str(r_ref)) << mesh
File("mesh_6B/thermal_block_physical_region_R%s.xml" % str(r_ref)) << subdomains
File("mesh_6B/thermal_block_facet_region_R%s.xml" % str(r_ref)) << boundaries
# XDMFFile("thermal_block_R%s.xdmf" % str(r_ref)).write(mesh)
# XDMFFile("thermal_block_physical_region_R%s.xdmf" % str(r_ref)).write(subdomains)
# XDMFFile("thermal_block_facet_region_R%s.xdmf" % str(r_ref)).write(boundaries)

# Visualize
num_node = mesh.num_vertices()
print(num_node)
plot(subdomains)
plot(mesh)
plt.show()
print(mesh.num_vertices())
print(mesh.num_cells())
