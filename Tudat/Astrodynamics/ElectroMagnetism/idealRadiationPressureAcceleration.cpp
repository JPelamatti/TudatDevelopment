/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      121123    D. Dirkx          File created.
 *      130124    K. Kumar          Added missing file header; updated layout; migrated force
 *                                  free function to separate file; added acceleration free
 *                                  function.
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureForce.h"

#include <cmath>

#include <vector>
#include <Eigen/Core>



namespace tudat
{
namespace electro_magnetism
{

Eigen::Vector3d normalToSail;
Eigen::Vector3d lightOrientedSolution;
Eigen::MatrixXd R_12(3,3);

const double gamma;
const double beta;


//! Compute radiation pressure acceleration using an ideal sail model.

Eigen::Vector3d computeIdealRadiationPressureAcceleration(
        const double radiationPressure,
        const Eigen::Vector3d& vectorFromSource,
        const Eigen::Vector2d& sailAngles,
        const double area,
        const double radiationPressureCoefficient,
        const double mass )
{

    normalToSail.push_back(cos(sailAngles(0)));
    normalToSail.push_back(sin(sailAngles(0))*sin(sailAngles(1)));
    normalToSail.push_back(sin(sailAngles(0))*cos(sailAngles(1)));

    gamma = asin(vectorFromSource(2)/vectorFromSource.norm());
    beta = atan(vectorFromSource(1)/vectorFromSource(0));

    R_12(0,0) = cos(beta);
    R_12(0,1) = -sin(beta);
    R_12(0,2) = -cos(beta)*sin(gamma);
    R_12(1,0) = sin(beta)*cos(gamma);
    R_12(1,1) = cos(beta);
    R_12(1,2) = -sin(beta)*sin(gamma);
    R_12(2,0) = sin(gamma);
    R_12(2,1) = 0;
    R_12(2,2) = cos(gamma);

    lightOrientedSolution = computeIdealRadiationPressureForce(
                radiationPressure, normalToSail, area, radiationPressureCoefficient, sailAngles) / mass;

    return R_12*lightOrientedSolution;

}

} // namespace electro_magnetism
} // namespace tudat
