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
            break;
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
            break;
        }
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
            break;
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
            break;
        }
    }
}

UnitVector CubedSphereMeshGeometry::faceNormal (int n1, int n2, int nr, int axis) const // modified
{
    Coordinate start = getVertexCoordinate(n1, n2, nr);
    Coordinate end1 = getVertexCoordinate(n1+1, n2, nr);
    Coordinate end2 = getVertexCoordinate(n1, n2+1, nr);
    Coordinate end3 = getVertexCoordinate(n1, n2, nr+1);

    Coordinate ihat = {{0,0,0}};
    Coordinate jhat = {{0,0,0}};
    Coordinate rhat = {{0,0,0}};
    Coordinate rad = getVertexCoordinate(n1+0.5, n2+0.5, nr);

    for (int d=0; d<3; ++d)
    {
        ihat[d] = end1[d] - start[d];
        jhat[d] = end2[d] - start[d];
        rhat[d] = end3[d] - start[d];
    }

    switch (axis)
    {
        case 0: return UnitVector::normalizeFrom(jhat[1] * rhat[2] - jhat[2] * rhat[1], jhat[2] * rhat[0] - jhat[0] * rhat[2], jhat[0] * rhat[1] - jhat[1] * rhat[0]);
        case 1: return UnitVector::normalizeFrom(rhat[1] * ihat[2] - rhat[2] * ihat[1], rhat[2] * ihat[0] - rhat[0] * ihat[2], rhat[0] * ihat[1] - rhat[1] * ihat[0]);
        case 2: return UnitVector::normalizeFrom(rad[0], rad[1], rad[2]);
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

double CubedSphereMeshGeometry::edgeLength (int n1, int n2, int nr, int axis) const // modified
{
    Coordinate start = getVertexCoordinate(n1, n2, nr);
    Coordinate end = {{0, 0, 0}};
    double dis = 0;

    switch (axis)
    {
        case 0:     
            end = getVertexCoordinate(n1+1, n2, nr);
            break;
        case 1:             
            end = getVertexCoordinate(n1, n2+1, nr);
            break;
        case 2:             
            end = getVertexCoordinate(n1, n2, nr+1);
            break;
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }

    for (int d=0; d<3; ++d)
    {
        dis += (start[d] - end[d]) * (start[d] - end[d]);
    }
    return std::sqrt(dis);
}

UnitVector CubedSphereMeshGeometry::edgeVector (int n1, int n2, int nr, int axis) const //modified
{
    Coordinate start = getVertexCoordinate(n1, n2, nr);
    Coordinate end = {{0, 0, 0}};

    switch (axis)
    {
        case 0: 
            end = getVertexCoordinate(n1+1, n2, nr);
            break;
        case 1: 
            end = getVertexCoordinate(n1, n2+1, nr);
            break;
        case 2: 
            end = getVertexCoordinate(n1, n2, nr+1);
            break;
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
    return UnitVector::normalizeFrom(end[0]-start[0], end[1]-start[1], end[2]-start[2]);

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

double CubedSphereMeshGeometry::getEdge (double n, int axis) const // modified
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
            break;
        case(MappingMode::equidistance):  
            e1 = -1. + n1 * 2. / shape[1];
            e2 = -1. + n2 * 2. / shape[2];
            break;
    }
    Rnorm = sqrt(e1*e1 + e2*e2 + 1);
    switch(location)
    {
        case(Location::top): // local coord: ( x,  y,  z)
            coord[0] = e1 * R / Rnorm;
            coord[1] = e2 * R / Rnorm;
            coord[2] = R / Rnorm;
            break;     
        case(Location::bottom):  // local coord: (-x,  y, -z)
            coord[0] = - e1 * R / Rnorm;
            coord[1] = e2 * R / Rnorm;
            coord[2] = - R / Rnorm;
            break;
        case(Location::left): // local coord: (-z,  x, -y)
            coord[0] = e2 * R / Rnorm;
            coord[1] = - R / Rnorm;     
            coord[2] = - e1 * R / Rnorm;
            break;
        case(Location::right):  // local coord: (-z, -x,  y)
            coord[0] = - e2 * R / Rnorm;
            coord[1] = R / Rnorm;  
            coord[2] = - e1 * R / Rnorm;
            break;
        case(Location::front): // local coord: (-z,  y,  x)
            coord[0] = R / Rnorm; 
            coord[1] = e2 * R / Rnorm;  
            coord[2] = - e1 * R / Rnorm;
            break;  
        case(Location::back): // local coord: (-z, -y, -x)
            coord[0] = - R / Rnorm; 
            coord[1] = - e2 * R / Rnorm;
            coord[2] = - e1 * R / Rnorm;
            break;
    }

    return coord;
}
