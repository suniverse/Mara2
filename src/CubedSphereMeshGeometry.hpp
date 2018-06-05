#pragma once
#include "Mara.hpp"

/* Cubed sphere mesh Geometry using gnomonic projection, for each of the six patches
 Use both equidistant and equiangular discretizations, set using discretization mode
 The shape is (# of cells on the cube faces, # in radial direction)
 The local coordinate system is set to be right-handed on each faces, the conventions are
 Top:       ( x,  y,  z)
 Bottom:    (-x,  y, -z)
 Left:      (-z,  x, -y)
 Right:     (-z, -x,  y)
 Front:     (-z,  y,  x)
 Back:      (-z, -y, -x)
*/

class CubedSphereMeshGeometry : public MeshGeometry
{
public:
    enum class Location { top, bottom, left, right, front, back, };
    enum class MappingMode { equiangle, equidistance, };

    using Shape = Cow::Shape;
    using Index = Cow::Index;

    // ========================================================================
    CubedSphereMeshGeometry (Location loc, MappingMode map = MappingMode::equiangle);
    CubedSphereMeshGeometry (Cow::Shape shape, Location loc, MappingMode map = MappingMode::equiangle);
    CubedSphereMeshGeometry (int ni, int nj, int nk, Location loc, MappingMode map = MappingMode::equiangle);
    Location getLocation() const;
    Coordinate getVertexCoordinate(int nr, int n1, int n2) const;

    // ========================================================================
    void setCellsShape (Cow::Shape S) override;
    void setLowerUpper (Coordinate L, Coordinate U) override;
    Cow::Shape cellsShape() const override;
    Cow::Index indexAtCoordinate (Coordinate x) const override;
    Coordinate coordinateAtIndex (double i, double j, double k) const override;
    unsigned long totalCellsInMesh() const override;
    double cellLength (int i, int j, int k, int axis) const override;
    double cellVolume (int i, int j, int k) const override;
    double meshVolume() const override;
    double faceArea (int i, int j, int k, int axis) const override;
    UnitVector faceNormal (int i, int j, int k, int axis) const override;
    double edgeLength (int i, int j, int k, int axis) const override;
    UnitVector edgeVector (int i, int j, int k, int axis) const override;
    Cow::Array getPointCoordinates (int axis) const override;
    std::shared_ptr<MeshGeometry> duplicate() const override;
    std::string getType() const override;

private:
    double getEdge (double i, int axis) const;
    double rpow(double x, double y, double n) const
    {
        return pow(std::sqrt(x * x + y * y + 1), n);
    };
    double cossq(double x) const
    {
        return cos(x) * cos(x);
    };

    Location location;
    MappingMode mappingMode;

    Coordinate lower;
    Coordinate upper;
    Shape shape;
};
