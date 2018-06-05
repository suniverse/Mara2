#include <cassert>
#include <cmath>
#include "CubedSphereMeshGeometry.hpp"




// ============================================================================
CubedSphereMeshGeometry::CubedSphereMeshGeometry (Location loc, MappingMode map) : location (loc), mappingMode (map)
{
    shape = {{ 128, 1, 1, 1, 1 }};
    lower = {{ -1.0, -1.0, 1.0 }};
    upper = {{ 1.0, 1.0, 10. }};
    location = loc;
    mappingMode = map;
}

CubedSphereMeshGeometry::CubedSphereMeshGeometry(Cow::Shape S, Location loc, MappingMode map) : location (loc), mappingMode (map)
{
    shape = {{ S[0], S[1], S[2], 1, 1 }};
    lower = {{ -1.0, -1.0, 1.0 }};
    upper = {{ 1.0, 1.0, 10. }};
    location = loc;
    mappingMode = map;

}

CubedSphereMeshGeometry::CubedSphereMeshGeometry (int ni, int nj, int nk, Location loc, MappingMode map) : location (loc), mappingMode (map)
{
    shape = {{ ni, nj, nk, 1, 1 }};
    lower = {{ -1.0, -1.0, 1.0 }};
    upper = {{ 1.0, 1.0, 10. }};
    location = loc;
    mappingMode = map;

}

void CubedSphereMeshGeometry::setCellsShape (Cow::Shape S)
{
    shape[0] = S[0];
    shape[1] = S[1];
    shape[2] = S[2];
}

void CubedSphereMeshGeometry::setLowerUpper (Coordinate L, Coordinate U)
{
    lower = L;
    upper = U;

    std::cout << "upper: " << " " << U[0] << " " << U[1] << " " << U[2] << std::endl;
}

Cow::Shape CubedSphereMeshGeometry::cellsShape() const
{
    return shape;
}

Cow::Index CubedSphereMeshGeometry::indexAtCoordinate (Coordinate x) const
{
    throw std::logic_error ("CubedSphereMeshGeometry::indexAtCoordinate not implemented");
    return Cow::Index();
}

Coordinate CubedSphereMeshGeometry::coordinateAtIndex (double i, double j, double k) const
{

    const double r = getEdge (i + 0.5, 0);
    const double q = getEdge (j + 0.5, 1);
    const double p = getEdge (k + 0.5, 2);
    return Coordinate ({{ r, q, p }});
}

unsigned long CubedSphereMeshGeometry::totalCellsInMesh() const
{
    return shape[0] * shape[1] * shape[2];
}

double CubedSphereMeshGeometry::cellLength (int n1, int n2, int nr, int axis) const // modified
{
    const double x0 = getEdge (n1, 0);
    const double y0 = getEdge (n2, 1);
    const double r = getEdge (nr, 2);

    const double x1 = getEdge (n1 + 1, 0);
    const double y1 = getEdge (n2 + 1, 1);

    switch (mappingMode)
    {
        case MappingMode::equidistance :
        {
            switch (axis)
            {
                case 0: return r * (x1 - x0) * (std::sqrt(y0 * y0 + 1) * rpow(x0, y0, -1.0) + std::sqrt(y0 * y0 + 1) * rpow(x1, y0, -1.0))/2.;
                case 1: return r * (y1 - y0) * (std::sqrt(x0 * x0 + 1) * rpow(x0, y0, -1.0) + std::sqrt(x0 * x0 + 1) * rpow(x0, y1, -1.0))/2.;
                case 2: return (getEdge (nr + 1, 2) - getEdge (nr, 2));
                default: throw std::logic_error ("Axis argument not 0, 1, or 2");
            }
        }
        case MappingMode::equiangle :
        {
            double xs = tan(x0);
            double xe = tan(x1);
            double ys = tan(y0);
            double ye = tan(y1);
            switch (axis)
            {
                case 0: return r * (xs - xe) * (std::sqrt(ys * ys + 1) * rpow(xs, ys, -1.0) / cossq(x0)
                 + std::sqrt(ys * ys + 1) * rpow(xe, ys, -1.0) / cossq(x1) )/2.;
                case 1: return r * (ys - ye) * (std::sqrt(xs * xs + 1) * rpow(xs, ys, -1.0) / cossq(y0)
                 + std::sqrt(xs * xs + 1) * rpow(xs, ye, -1.0) / cossq(y1) )/2.;
                case 2: return (getEdge (nr + 1, 2) - getEdge (nr, 2));
                default: throw std::logic_error ("Axis argument not 0, 1, or 2");
            }
        }
        default: throw std::logic_error ("Mapping Mode not gnomonic");
    }
}

double CubedSphereMeshGeometry::cellVolume (int n1, int n2, int nr) const // modified
{
    const double s0 = faceArea(n1, n2, nr, 2);
    const double s1 = faceArea(n1, n2, nr+1, 2);
    const double dr = edgeLength(n1, n2, nr, 2);

    return (s0 + s1) * dr /2.;
}

double CubedSphereMeshGeometry::meshVolume() const // modified
{
    const double r0 = lower[2];
    const double r1 = upper[2];
    const double PI = 3.14159265358979;

    return 2 * PI / 9 * (r1 * r1 * r1 - r0 * r0 * r0);
}

double CubedSphereMeshGeometry::faceArea (int n1, int n2, int nr, int axis) const // modified
{
    const double x0 = getEdge (n1, 0);
    const double y0 = getEdge (n2, 1);
    const double r = getEdge (nr, 2);

    const double x1 = getEdge (n1 + 1, 0);
    const double y1 = getEdge (n2 + 1, 1);

    switch (mappingMode)
    {
        case MappingMode::equidistance :
        {
            switch (axis)
            {
                case 0: return (edgeLength(n1, n2, nr, 1) + edgeLength(n1, n2, nr+1, 1)) / 2. * edgeLength(n1, n2, nr, 2);
                case 1: return (edgeLength(n1, n2, nr, 0) + edgeLength(n1, n2, nr+1, 0)) / 2. * edgeLength(n1, n2, nr, 2);
                case 2: return r * r * (x1 - x0) * (y1 - y0) * (rpow(x0, y0,-3.0) + rpow(x0, y1, -3.0) + rpow(x1, y0, -3.0) + rpow(x1,y1,-3.0)) / 4.;
                default: throw std::logic_error ("Axis argument not 0, 1, or 2");
            }
        }
        case MappingMode::equiangle :
        {
            double xs = tan(x0);
            double xe = tan(x1);
            double ys = tan(y0);
            double ye = tan(y1);
            switch (axis)
            {
                case 0: return (edgeLength(n1, n2, nr, 1) + edgeLength(n1, n2, nr+1, 1)) / 2. * edgeLength(n1, n2, nr, 2);
                case 1: return (edgeLength(n1, n2, nr, 0) + edgeLength(n1, n2, nr+1, 0)) / 2. * edgeLength(n1, n2, nr, 2);
                case 2: return r * r * (x1 - x0) * (y1 - y0) * (rpow(xs, ys,-3.0)/(cossq(x0) * cossq(y0)) + rpow(xs, ye, -3.0)/(cossq(x0) * cossq(y1)) + rpow(xe, ys, -3.0)/(cossq(x1) * cossq(y0)) + rpow(xe, ye, -3.0)/(cossq(x1) * cossq(y1))) / 4.;
                default: throw std::logic_error ("Axis argument not 0, 1, or 2");
            }
        }
        default: throw std::logic_error ("Mapping Mode not gnomonic");
    }
}

UnitVector CubedSphereMeshGeometry::faceNormal (int i, int j, int k, int axis) const
{
    switch (axis)
    {
        case 0: return UnitVector::xhat;
        case 1: return UnitVector::yhat;
        case 2: return UnitVector::zhat;
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

double CubedSphereMeshGeometry::edgeLength (int n1, int n2, int nr, int axis) const // modified
{
    const double x0 = getEdge (n1, 0);
    const double y0 = getEdge (n2, 1);
    const double r = getEdge (nr, 2);

    const double x1 = getEdge (n1 + 1, 0);
    const double y1 = getEdge (n2 + 1, 1);

    switch (mappingMode)
    {
        case MappingMode::equidistance :
        {
            switch (axis)
            {
                case 0: return r * (x1 - x0) * (std::sqrt(y0 * y0 + 1) * rpow(x0, y0, -1.0) + std::sqrt(y0 * y0 + 1) * rpow(x1, y0, -1.0))/2.;
                case 1: return r * (y1 - y0) * (std::sqrt(x0 * x0 + 1) * rpow(x0, y0, -1.0) + std::sqrt(x0 * x0 + 1) * rpow(x0, y1, -1.0))/2.;
                case 2: return (getEdge (nr + 1, 2) - getEdge (nr, 2));
                default: throw std::logic_error ("Axis argument not 0, 1, or 2");
            }
        }
        case MappingMode::equiangle :
        {
            double xs = tan(x0);
            double xe = tan(x1);
            double ys = tan(y0);
            double ye = tan(y1);
            switch (axis)
            {
                case 0: return r * (xs - xe) * (std::sqrt(ys * ys + 1) * rpow(xs, ys, -1.0) / cossq(x0)
                 + std::sqrt(ys * ys + 1) * rpow(xe, ys, -1.0) / cossq(x1) )/2.;
                case 1: return r * (ys - ye) * (std::sqrt(xs * xs + 1) * rpow(xs, ys, -1.0) / cossq(y0)
                 + std::sqrt(xs * xs + 1) * rpow(xs, ye, -1.0) / cossq(y1) )/2.;
                case 2: return (getEdge (nr + 1, 2) - getEdge (nr, 2));
                default: throw std::logic_error ("Axis argument not 0, 1, or 2");
            }
        }
        default: throw std::logic_error ("Mapping Mode not gnomonic");
    }
}

UnitVector CubedSphereMeshGeometry::edgeVector (int i, int j, int k, int axis) const
{
    switch (axis)
    {
        case 0: return UnitVector::xhat;
        case 1: return UnitVector::yhat;
        case 2: return UnitVector::zhat;
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

Cow::Array CubedSphereMeshGeometry::getPointCoordinates (int axis) const
{
    if (axis < 0 || 2 < axis) throw std::logic_error ("Invalid axis");
    auto coords = Cow::Array (shape[axis] + 1);

    for (int n = 0; n < coords.size(); ++n)
    {
        coords[n] = coordinateAtIndex (n - 0.5, n - 0.5, n - 0.5)[axis];
    }
    return coords;
}

std::shared_ptr<MeshGeometry> CubedSphereMeshGeometry::duplicate() const
{
    auto mg = new CubedSphereMeshGeometry(location, mappingMode);
    *mg = *this;
    return std::shared_ptr<MeshGeometry> (mg);
}

std::string CubedSphereMeshGeometry::getType() const
{
    return "CubedSphere";
}

double CubedSphereMeshGeometry::getEdge (double n, int axis) const
{
    switch (axis)
    {
        case 0: return lower[0] + (upper[0] - lower[0]) * n / shape[0];
        case 1: return lower[1] + (upper[1] - lower[1]) * n / shape[1];
        case 2: return lower[2] * std::pow (upper[2] / lower[2], n / shape[2]);

        default: assert (false);
    }
}

CubedSphereMeshGeometry::Location CubedSphereMeshGeometry::getLocation() const
{ 
    return location; 
}

Coordinate CubedSphereMeshGeometry::getVertexCoordinate(int n1, int n2, int nr) const // modified
{
    //bound check
    if (n1 < 0 || n1 > shape[0]) throw std::logic_error ("Vertex out of mesh in the i direction");
    if (n2 < 0 || n2 > shape[1]) throw std::logic_error ("Vertex out of mesh in the j direction");
    if (nr < 0 || nr > shape[2]) throw std::logic_error ("Vertex out of mesh in the radial direction");

    const double PI = 3.14159265358979;
    Coordinate coord = {{0,0,0}};
    double e1, e2, R, Rnorm;
    R = lower[2] * std::pow (upper[2] / lower[2], nr / shape[2]);
    switch(mappingMode)
    {
        case(MappingMode::equiangle):     
            e1 = tan(-PI/4. + n1 * PI/2. / shape[1]);
            e2 = tan(-PI/4. + n2 * PI/2. / shape[2]);
        case(MappingMode::equidistance):  
            e1 = -1. + n1 * 2. / shape[1];
            e2 = -1. + n2 * 2. / shape[2];
        default: throw std::logic_error ("Mapping Mode not gnomonic");

    }
    Rnorm = sqrt(e1*e1 + e2*e2 + 1);
    switch(location)
    {
        case(Location::top): // local coord: ( x,  y,  z)
            coord[0] = e1 * R / Rnorm;
            coord[1] = e2 * R / Rnorm;
            coord[2] = R / Rnorm;     
        case(Location::bottom):  // local coord: (-x,  y, -z)
            coord[0] = - e1 * R / Rnorm;
            coord[1] = e2 * R / Rnorm;
            coord[2] = - R / Rnorm;
        case(Location::left): // local coord: (-z,  x, -y)
            coord[0] = e2 * R / Rnorm;
            coord[1] = - R / Rnorm;     
            coord[2] = - e1 * R / Rnorm;
        case(Location::right):  // local coord: (-z, -x,  y)
            coord[0] = - e2 * R / Rnorm;
            coord[1] = R / Rnorm;  
            coord[2] = - e1 * R / Rnorm;
        case(Location::front): // local coord: (-z,  y,  x)
            coord[0] = R / Rnorm; 
            coord[1] = e2 * R / Rnorm;  
            coord[2] = - e1 * R / Rnorm;  
        case(Location::back): // local coord: (-z, -y, -x)
            coord[0] = - R / Rnorm; 
            coord[1] = - e2 * R / Rnorm;
            coord[2] = - e1 * R / Rnorm;
        default: throw std::logic_error ("Location not on one of six faces");
    }

    return coord;
}
