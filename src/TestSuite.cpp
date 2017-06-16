#define CATCH_CONFIG_FAST_COMPILE
#include "Catch.hpp"
#include "Array.hpp"

using namespace Cow;




// NOTE: TestSuite class implementation is in Mara.cpp to speed compilation
// time.




// ============================================================================
#include "Stencil.hpp"

SCENARIO ("Stencil operations should behave correctly")
{
    GIVEN ("A stencil with shape [2, 3, 4]")
    {
        auto stencilShape = Shape();
        auto stencil = Stencil();
        auto operation = [&] (const Array& sd, Array& od) { stencilShape = sd.shape(); };
        stencil.setFootprintLower (-2,-3,-4);
        stencil.setFootprintUpper ( 2, 3, 4);
        stencil.setCodomainRank (4, 5);

        WHEN ("The input array is [12, 12, 12] and the codomain rank is [4, 5]")
        {
            auto source = Array (12, 12, 12);
            auto result = stencil.evaluate (operation, source);

            THEN ("The output array shape is [8, 6, 4, 4, 5]")
            {
                CHECK (result.size(0) == 8);
                CHECK (result.size(1) == 6);
                CHECK (result.size(2) == 4);
                CHECK (result.size(3) == 4);
                CHECK (result.size(4) == 5);
            }

            THEN ("The array of stencil data has shape [5, 7, 9, 1, 1]")
            {
                CHECK (stencilShape[0] == 5);
                CHECK (stencilShape[1] == 7);
                CHECK (stencilShape[2] == 9);
                CHECK (stencilShape[3] == 1);
                CHECK (stencilShape[4] == 1);
            }
        }
    }

    GIVEN ("A stencil with shape [2, 2, 0]")
    {
        auto stencilShape = Shape();
        auto stencil = Stencil();
        auto operation = [&] (const Array& sd, Array& od) { stencilShape = sd.shape(); };
        auto source = Array (12, 12, 1, 8);
        stencil.setFootprintLower (-2,-2, 0);
        stencil.setFootprintUpper ( 2, 2, 0);
        stencil.setCodomainRank (3, 1);

        WHEN ("The input array is [12, 12, 1, 8] and the codomain rank is [3, 1]")
        {
            auto result = stencil.evaluate (operation, source);

            THEN ("The output array shape is [8, 8, 1, 3, 1]")
            {
                CHECK (result.size(0) == 8);
                CHECK (result.size(1) == 8);
                CHECK (result.size(2) == 1);
                CHECK (result.size(3) == 3);
                CHECK (result.size(4) == 1);
            }

            THEN ("The array of stencil data has shape [5, 5, 1, 8, 1]")
            {
                CHECK (stencilShape[0] == 5);
                CHECK (stencilShape[1] == 5);
                CHECK (stencilShape[2] == 1);
                CHECK (stencilShape[3] == 8);
                CHECK (stencilShape[4] == 1);
            }
        }

        WHEN ("The stencil footprint is changed to [2, 2, 1]")
        {
            stencil.setFootprintLower (-2,-2,-1);
            stencil.setFootprintUpper ( 2, 2, 1);

            THEN ("An exception is raised because the source data shape is too small")
            {
                CHECK_THROWS_AS (stencil.evaluate (operation, source), std::logic_error);
            }
        }
    }

    GIVEN ("A second order finite difference stencil")
    {
        auto stencil = Stencil();
        auto x = Array (256);
        auto y = Array (256);

        for (int i = 0; i < 256; ++i)
        {
            x[i] = (i + 0.5) / 256;
            y[i] = std::sin (x[i]);
        }
        stencil.setCodomainRank (2, 1);

        WHEN ("The stencil is centered, of the type [-1, 0, +1]")
        {
            stencil.setFootprintLower (-1, 0, 0);
            stencil.setFootprintUpper ( 1, 0, 0);

            auto operation = [] (const Array& x0, const Array& y0, Array& deriv)
            {
                deriv[0] = (y0[2] - y0[0]) / (x0[2] - x0[0]);
                deriv[1] = std::cos (x0[1]);
            };
            auto deriv = stencil.evaluate (operation, x, y);

            THEN ("The resulting array is smaller by two data points")
            {
                CHECK (deriv.size(0) == 254);
            }

            THEN ("The derivative approximation is correct to one part in 10^5")
            {
                CHECK (deriv (  0, 0, 0, 0) == Approx (deriv (  0, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (128, 0, 0, 0) == Approx (deriv (128, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (253, 0, 0, 0) == Approx (deriv (253, 0, 0, 1)).epsilon (1e-5));
            }
        }

        WHEN ("The stencil is centered, of the type [0, +1]")
        {
            stencil.setFootprintLower (0, 0, 0);
            stencil.setFootprintUpper (1, 0, 0);

            auto operation = [] (const Array& x0, const Array& y0, Array& deriv)
            {
                deriv[0] = (y0[1] - y0[0]) / (x0[1] - x0[0]);
                deriv[1] = std::cos (0.5 * (x0[0] + x0[1]));
            };
            auto deriv = stencil.evaluate (operation, x, y);

            THEN ("The resulting array is smaller by one data point")
            {
                CHECK (deriv.size(0) == 255);
            }

            THEN ("The derivative approximation is correct to one part in 10^5")
            {
                CHECK (deriv (  0, 0, 0, 0) == Approx (deriv (  0, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (128, 0, 0, 0) == Approx (deriv (128, 0, 0, 1)).epsilon (1e-5));
                CHECK (deriv (254, 0, 0, 0) == Approx (deriv (254, 0, 0, 1)).epsilon (1e-5));
            }
        }
    }
}




// ============================================================================
#include "Mara.hpp"
#include "CartesianMeshGeometry.hpp"
#include "BoundaryConditions.hpp"
#include "ConservationLaws.hpp"
#include "ConstrainedTransport.hpp"
#include "IntercellFluxSchemes.hpp"
#include "RiemannSolvers.hpp"

SCENARIO ("Mara session should launch if given minimal setup")
{
    auto session = MaraSession();
    auto setup = SimulationSetup();

    setup.finalTime               = 0.001;
    setup.checkpointInterval      = 0.0;
    setup.vtkOutputInterval       = 0.0;
    setup.boundaryCondition       = std::make_shared<PeriodicBoundaryCondition>();
    setup.conservationLaw         = std::make_shared<NewtonianHydro>();
    setup.riemannSolver           = std::make_shared<HlleRiemannSolver>();
    setup.meshGeometry            = std::make_shared<CartesianMeshGeometry>();
    setup.intercellFluxScheme     = std::make_shared<MethodOfLinesPlm>();
    setup.constrainedTransport    = std::make_shared<UniformCartesianCT>();
    setup.initialDataFunction     = [] (double x, double, double)
    {
        return std::vector<double> { 1.0 + 0.1 * std::sin (2 * M_PI * x), 0., 0., 0., 1. };
    };

    session.getLogger()->setLogToNull();
    CHECK (session.launch (setup).simulationIter == 1);
}




// ============================================================================
#include "HDF5.hpp"
#include "TimeSeriesManager.hpp"

SCENARIO ("Time series manager should behave reasonably")
{
    GIVEN ("A time series manager populated with data series, with unequal lengths")
    {
        TimeSeriesManager timeSeriesManager;
        timeSeriesManager.getLogger()->setLogToNull();
        timeSeriesManager.append ("mean_energy", 1.1);
        timeSeriesManager.append ("num_unhealthy_zones", 12);
        timeSeriesManager.append ("mean_pressure", 2.2);
        timeSeriesManager.append ("num_unhealthy_zones", 13);

        WHEN ("Columns are queried")
        {
            THEN ("They have the expected size")
            {
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 2);
                CHECK (timeSeriesManager.getSeriesNamesInt().size() == 1);
                CHECK (timeSeriesManager.getSeriesDouble ("mean_energy").size() == 1);
                CHECK (timeSeriesManager.getSeriesDouble ("mean_pressure").size() == 1);
                CHECK (timeSeriesManager.getSeriesInt ("num_unhealthy_zones").size() == 2);
                CHECK_THROWS_AS (timeSeriesManager.getSeriesInt ("not_there"), std::out_of_range);
            }
        }

        WHEN ("Time series data is written to HDF5")
        {
            auto hdf5File = Cow::H5::File ("time_series_test.h5", "w");
            timeSeriesManager.write (hdf5File);

            THEN ("There are the expected number of columns in the file")
            {
                int numDatasets = 0;
                hdf5File.iterate ([&] (std::string) { numDatasets++; });
                CHECK (numDatasets == 3);
            }

            THEN ("Time series data can be cleared and restored from the file")
            {
                timeSeriesManager.clear();
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 0);

                timeSeriesManager.load (hdf5File);
                CHECK (timeSeriesManager.getSeriesNamesDouble().size() == 2);
                CHECK (timeSeriesManager.getSeriesInt ("num_unhealthy_zones").size() == 2);
            }
        }
    }
}




// ============================================================================
#include "MeshOperator.hpp"
#include "BoundaryConditions.hpp"

SCENARIO ("Mesh operator should generate data, take div and curl", "[MeshOperator]")
{
    GIVEN ("A 1D mesh operator and cartesian geometry")
    {
        auto fn = [] (double x, double, double)
        {
            return std::vector<double> { std::sin (2 * M_PI * x), std::cos (2 * M_PI * x) };
        };
        auto mo = std::make_shared<MeshOperator>();
        auto mg = std::make_shared<CartesianMeshGeometry>();
        auto bc = std::make_shared<PeriodicBoundaryCondition>();

        mo->setMeshGeometry (mg);

        WHEN ("Data is generated in volumes")
        {
            auto P = mo->generate (fn, MeshLocation::cell);
    
            THEN ("The returned array has the expected size")
            {
                CHECK (P.size(0) == mg->cellsShape()[0]);
                CHECK (P.size(3) == 2);
            }

            THEN ("The returned array has the expected data")
            {
                CHECK (P (0, 0, 0, 0) == fn (mg->coordinateAtIndex(0, 0, 0)[0], 0., 0.)[0]);
                CHECK (P (0, 0, 0, 1) == fn (mg->coordinateAtIndex(0, 0, 0)[0], 0., 0.)[1]);
            }
        }

        WHEN ("Two guard zones are left untouched by the initial data function")
        {
            int numGuard = 2;
            auto fullShape = mg->cellsShape();
            fullShape[0] += 2 * numGuard;
            fullShape[3] = 2;

            auto interiorRegion = Region().withLower (0, numGuard).withUpper (0, -numGuard);
            auto P = Array (fullShape);

            P[interiorRegion] = mo->generate (fn, MeshLocation::cell);
            bc->apply (P, MeshLocation::cell, MeshBoundary::left,  0, numGuard);
            bc->apply (P, MeshLocation::cell, MeshBoundary::right, 0, numGuard);

            THEN ("Periodic boundary conditions assign the expected values to the guard zones")
            {
                CHECK (P (0, 0, 0, 0) == P (mg->cellsShape()[0], 0, 0, 0));
                CHECK (P (0, 0, 0, 1) == P (mg->cellsShape()[0], 0, 0, 1));
                CHECK (P (mg->cellsShape()[0] + numGuard, 0, 0, 0) == P (numGuard, 0, 0, 0));
                CHECK (P (mg->cellsShape()[0] + numGuard, 0, 0, 1) == P (numGuard, 0, 0, 1));
            }
        }

        WHEN ("Data is generated on faces")
        {
            auto B = mo->generate (fn, MeshLocation::face);
            auto M = mo->divergence (B);

            THEN ("The returned array has the expected shape")
            {
                CHECK (B.size(0) == mg->cellsShape()[0] + 1);
                CHECK (B.size(3) == 2);
                CHECK (B.size(4) == 3);
            }

            THEN ("The array of corresponding divergences has the correct shape")
            {
                CHECK (M.size(0) == mg->cellsShape()[0]);
            }
        }
    }

    GIVEN ("A 3D mesh operator and a [16, 32, 8] cartesian geometry")
    {
        auto abcField = [] (double x, double y, double z)
        {
            const double k = 2 * M_PI;
            const double A = 1.0;
            const double B = 1.0;
            const double C = 1.0;
            const double b1 = C * std::cos (k * z) - B * std::sin (k * y);
            const double b2 = A * std::cos (k * x) - C * std::sin (k * z);
            const double b3 = B * std::cos (k * y) - A * std::sin (k * x);
            return std::vector<double> { b1, b2, b3 };
        };
        auto divLikeField = [] (double x, double y, double z)
        {
            const double k = 2 * M_PI;
            const double b1 = std::cos (k * x);
            const double b2 = std::cos (k * y);
            const double b3 = std::cos (k * z);
            return std::vector<double> { b1, b2, b3 };
        };
        auto mo = std::make_shared<MeshOperator>();
        auto mg = std::make_shared<CartesianMeshGeometry>();

        mg->setCellsShape ({{16, 32, 8}});
        mo->setMeshGeometry (mg);

        WHEN ("An ABC field is generated on faces in flux mode")
        {
            auto B = mo->generate (abcField, MeshLocation::face, MeshOperator::VectorMode::fluxish);
            auto M = mo->divergence (B);

            THEN ("The B field has the expected shape [17, 33, 9] with one rank (1, 3)")
            {
                CHECK (B.size(0) == 17);
                CHECK (B.size(1) == 33);
                CHECK (B.size(2) == 9);
                CHECK (B.size(3) == 1);
                CHECK (B.size(4) == 3);
            }

            THEN ("The monopole field has the expected shape [16, 32, 8] with rank (1, 1)")
            {
                CHECK (M.size(0) == 16);
                CHECK (M.size(1) == 32);
                CHECK (M.size(2) == 8);
                CHECK (M.size(3) == 1);
                CHECK (M.size(4) == 1);
            }

            THEN ("The monopole field is zero at some randomly checked locations")
            {
                CHECK (M (2, 6, 6) == Approx (0.0));
                CHECK (M (6, 0, 6) == Approx (0.0));
                CHECK (M (6, 6, 2) == Approx (0.0));
                CHECK (M (0, 1, 0) == Approx (0.0));
            }
        }

        WHEN ("A non-solenoidal field is generated on faces in flux mode")
        {
            auto V = mo->generate (divLikeField, MeshLocation::face, MeshOperator::VectorMode::fluxish);
            auto D = mo->divergence (V);

            THEN ("The divergence field is non-zero at some randomly checked locations")
            {
                CHECK (D (2, 6, 6) != Approx (0.0));
                CHECK (D (6, 0, 6) != Approx (0.0));
                CHECK (D (6, 6, 2) != Approx (0.0));
                CHECK (D (0, 1, 0) != Approx (0.0));
            }
        }

        WHEN ("An ABC field is generated on edges in scalars mode")
        {
            auto A = mo->generate (abcField, MeshLocation::edge, MeshOperator::VectorMode::scalars);

            THEN ("The A field has the expected shape [17, 33, 9] with rank (3, 3)")
            {
                CHECK (A.size(0) == 17);
                CHECK (A.size(1) == 33);
                CHECK (A.size(2) == 9);
                CHECK (A.size(3) == 3);
                CHECK (A.size(4) == 3);
            }
        }

        WHEN ("An ABC field is generated on edges in EMF mode")
        {
            auto A = mo->generate (abcField, MeshLocation::edge, MeshOperator::VectorMode::emflike);
            auto B = mo->curl (A);
            auto M = mo->divergence (B);

            THEN ("The A field has the expected shape [17, 33, 9] with rank (1, 3)")
            {
                CHECK (A.size(0) == 17);
                CHECK (A.size(1) == 33);
                CHECK (A.size(2) == 9);
                CHECK (A.size(3) == 1);
                CHECK (A.size(4) == 3);
            }

            THEN ("The curl of A has the expected shape [17, 33, 9] with rank (1, 3)")
            {
                CHECK (B.size(0) == 17);
                CHECK (B.size(1) == 33);
                CHECK (B.size(2) == 9);
                CHECK (B.size(3) == 1);
                CHECK (B.size(4) == 3);
            }

            THEN ("The div of B has the expected shape [16, 32, 8] with rank (1, 1)")
            {
                CHECK (M.size(0) == 16);
                CHECK (M.size(1) == 32);
                CHECK (M.size(2) == 8);
                CHECK (M.size(3) == 1);
                CHECK (M.size(4) == 1);
            }

            THEN ("The div of B is zero at some randomly checked locations")
            {
                CHECK (M (2, 6, 6) == Approx (0.0));
                CHECK (M (6, 0, 6) == Approx (0.0));
                CHECK (M (6, 6, 2) == Approx (0.0));
                CHECK (M (0, 1, 0) == Approx (0.0));
            }
        }
    }
}




// ============================================================================
#include "MeshOperator.hpp"

SCENARIO ("Mesh operator should run flux sweeps", "[MeshOperator]")
{
    GIVEN ("A 1D mesh operator, cartesian geometry, and a trivial Godunov flux function")
    {
        auto mg = std::make_shared<CartesianMeshGeometry>(8, 1, 1);
        auto mo = std::make_shared<MeshOperator>(mg);

        auto cd = [] (double x, double, double) { return std::vector<double> { std::sin (M_PI * x) }; };
        auto fd = [] (double x, double, double) { return std::vector<double> { 0., 1., 2. }; };
        auto gd = [] (GodunovStencil& stencil)
        {
            if (stencil.faceNormal == UnitVector::xhat)
            {
                CHECK (stencil.cellData.size(0) == 2);
                CHECK (stencil.cellData.size(1) == 1);
                CHECK (stencil.faceData.size(0) == 3); // is a 1D array of data on a single face
                CHECK (stencil.faceData[0] == 0.);
                CHECK (stencil.faceData[1] == 1.);
                CHECK (stencil.faceData[2] == 2.);
            }
            else
            {
                CHECK (stencil.cellData.size(0) == 0);
                CHECK (stencil.cellData.size(1) == 1);
                CHECK (stencil.faceData.size(0) == 3);
            }
        };

        auto cellData = mo->generate (cd, MeshLocation::cell);
        auto faceData = mo->generate (fd, MeshLocation::face);

        THEN ("Initial data (face and cell) has the expected shape")
        {
            CHECK (cellData.size(0) == 8);
            CHECK (cellData.size(3) == 1);
            CHECK (faceData.size(0) == 9);
            CHECK (faceData.size(3) == 3);
        }

        auto footprint = Shape {{ 2, 0, 0 }};
        auto faceFlux = mo->godunov (gd, cellData, faceData, footprint);

        THEN ("Godunov fluxes have the expected shape")
        {
            CHECK (faceFlux.size(0) == 9);
            CHECK (faceFlux.size(1) == 2);
            CHECK (faceFlux.size(2) == 2);
            CHECK (faceFlux.size(3) == 1);
        }
    }

    GIVEN ("A 3D mesh operator, cartesian geometry, and a trivial Godunov flux function")
    {
        auto mg = std::make_shared<CartesianMeshGeometry>(8, 8, 8);
        auto mo = std::make_shared<MeshOperator>(mg);

        auto id = [] (double x, double y, double z) { return std::vector<double> { 0., 1., 2. }; };
        auto cd = [] (double x, double y, double z) { return std::vector<double> { }; };
        auto gd = [] (GodunovStencil& stencil)
        {
            CHECK (stencil.cellData.size(0) == 2);
            CHECK (stencil.cellData.size(1) == 3);
            CHECK (stencil.cellCoords.size() == 2);
            CHECK (stencil.godunovFlux.size(0) == 3);
            stencil.godunovFlux[0] = stencil.cellCoords[0][0]; // Fhat = x of cell to left of face
            stencil.godunovFlux[1] = stencil.cellCoords[0][1]; // Fhat = y of cell to left of face
            stencil.godunovFlux[2] = stencil.cellCoords[0][2]; // Fhat = z of cell to left of face
        };

        auto cellData = mo->generate (id, MeshLocation::cell);
        auto faceData = mo->generate (cd, MeshLocation::face);

        THEN ("Initial data (face and cell) has the expected shape")
        {
            CHECK (cellData.size(0) == 8);
            CHECK (cellData.size(1) == 8);
            CHECK (cellData.size(2) == 8);
            CHECK (cellData.size(3) == 3);
            CHECK (faceData.size(0) == 9);
            CHECK (faceData.size(1) == 9);
            CHECK (faceData.size(2) == 9);
            CHECK (faceData.size(3) == 0);
        }

        auto footprint = Shape {{ 2, 2, 2 }};
        auto start = Shape {{ -1, 0, +1 }}; // choose an arbitrary start location
        auto godunovFlux = mo->godunov (gd, cellData, faceData, footprint, start);

        THEN ("Godunov fluxes have the expected shape")
        {
            CHECK (godunovFlux.size(0) == 9);
            CHECK (godunovFlux.size(1) == 9);
            CHECK (godunovFlux.size(2) == 9);
            CHECK (godunovFlux.size(3) == 3);
        }

        THEN ("Godunov fluxes have values consistent with start index")
        {
            CHECK (godunovFlux (1, 0, 0, 0, 0) == mg->coordinateAtIndex (start[0], 0, 0)[0]);
            CHECK (godunovFlux (0, 1, 0, 1, 1) == mg->coordinateAtIndex (0, start[1], 0)[1]);
            CHECK (godunovFlux (0, 0, 1, 2, 2) == mg->coordinateAtIndex (0, 0, start[2])[2]);
        }
    }
}




// ============================================================================
#include "FieldOperator.hpp"

SCENARIO ("Field operator should convert prim -> cons", "[FieldOperator]")
{
    GIVEN ("A field operator and a 1D mesh")
    {
        auto mg = std::make_shared<CartesianMeshGeometry>();
        auto mo = std::make_shared<MeshOperator>(mg);
        auto cl = std::make_shared<NewtonianHydro>();
        auto fo = std::make_shared<FieldOperator>(cl);
        auto id = [] (double x, double, double)
        {
            return std::vector<double> { 1. + 0.5 * std::cos (M_PI * x), 1., 2., 3., 1. };
        };

        WHEN ("The mesh operator is used to generate primitives and measures")
        {
            auto P = mo->generate (id, MeshLocation::cell);
            auto V = mo->measure (MeshLocation::cell);
            auto L = mo->linearCellDimension();

            THEN ("The cell linear dimension is correct")
            {
                CHECK (L.shape() == mg->cellsShape());
                CHECK (L (0, 0, 0) == mg->cellLength (0, 0, 0, 0));
            }

            THEN ("The cell volume is correct")
            {
                CHECK (V.shape() == mg->cellsShape());
                CHECK (V (0, 0, 0) == mg->cellVolume (0, 0, 0));
            }
        }

        WHEN ("The field operator is used to generate conserved variables")
        {
            auto P = mo->generate (id, MeshLocation::cell);
            auto X = mo->cellCentroidCoordinates();
            auto U = fo->generateConserved (P);
            auto Q = fo->recoverPrimitive (U);

            THEN ("The conserved variables are as expected")
            {
                for (int i = 0; i < U.size(0); ++i)
                {
                    auto R = ConservationLaw::Request();
                    auto Pi = id (X (i), 0, 0);
                    auto Si = cl->fromPrimitive (R, &Pi[0]);

                    CHECK (U (i, 0, 0, 0) == Si.U[0]);
                    CHECK (U (i, 0, 0, 1) == Si.U[1]);
                    CHECK (U (i, 0, 0, 2) == Si.U[2]);
                    CHECK (U (i, 0, 0, 3) == Si.U[3]);
                    CHECK (U (i, 0, 0, 4) == Si.U[4]);
                }
            }

            THEN ("The re-generated primatives are the same as the originals")
            {
                CHECK (P.size() == Q.size());
                
                for (int n = 0; n < P.size(); ++n)
                {
                    CHECK (P[n] == Approx (Q[n]));
                }
            }
        }

        WHEN ("The and mesh field operators are used to generate the Courant time and diagnostics")
        {
            auto id = [] (double, double, double) { return std::vector<double> { 1., 1., 0., 0., 0. }; };
            auto P = mo->generate (id, MeshLocation::cell);
            auto V = mo->measure (MeshLocation::cell);
            auto L = mo->linearCellDimension();
            auto D = fo->volumeIntegratedDiagnostics (P, V);
            auto N = cl->getDiagnosticNames();

            THEN ("The Courant time is as expected")
            {
                CHECK (fo->getCourantTimestep (P, L) == mg->cellLength (0, 0, 0, 0));
            }

            THEN ("The diagnostics are reasonable")
            {
                CHECK (N[0] == "mass");
                CHECK (D[0] == mg->meshVolume());
            }
        }
    }
}




// ============================================================================
#include "UnitVector.hpp"
#include "RotationMatrix.hpp"

SCENARIO ("Rotation matrices should behave appropriately with unit vectors", "[UnitVector]")
{
    GIVEN ("The unit vector zhat")
    {
        UnitVector xhat = UnitVector::fromCartesian (1, 0, 0);
        UnitVector yhat = UnitVector::fromCartesian (0, 1, 0);
        UnitVector zhat = UnitVector::fromCartesian (0, 0, 1);

        THEN ("Y (pi / 2) zhat = xhat")
        {
            REQUIRE ((RotationMatrix::aboutY (M_PI / 2) * zhat).pitchAngleMu == Approx (xhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutY (M_PI / 2) * zhat).azimuthalAnglePhi == Approx (xhat.azimuthalAnglePhi));
        }

        THEN ("X (-pi / 2) zhat = yhat")
        {
            REQUIRE ((RotationMatrix::aboutX (-M_PI / 2) * zhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutX (-M_PI / 2) * zhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("Z (pi / 2) xhat = yhat")
        {
            REQUIRE ((RotationMatrix::aboutZ (M_PI / 2) * xhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutZ (M_PI / 2) * xhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (xhat) = xhat")
        {
            REQUIRE (zhat.withPolarAxis (xhat).pitchAngleMu == Approx (xhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (xhat).azimuthalAnglePhi == Approx (xhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (yhat) = yhat")
        {
            REQUIRE (zhat.withPolarAxis (yhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (yhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (zhat) = zhat")
        {
            REQUIRE (zhat.withPolarAxis (zhat).pitchAngleMu == Approx (zhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (zhat).azimuthalAnglePhi == Approx (zhat.azimuthalAnglePhi));
        }
    }

    REQUIRE (UnitVector::xhat != UnitVector::yhat);
    REQUIRE (UnitVector::xhat == UnitVector::xhat);
    REQUIRE (UnitVector::yhat == UnitVector::yhat);
    REQUIRE (UnitVector::zhat == UnitVector::zhat);
}