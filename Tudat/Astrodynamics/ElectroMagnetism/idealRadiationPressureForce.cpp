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
 *      130124    K. Kumar          Migrated force free function to new file.
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/Astrodynamics/ElectroMagnetism/idealRadiationPressureForce.h"

#include <iostream>

#include <cmath>


namespace tudat
{
namespace electro_magnetism
{

//! Compute radiation pressure force using an ideal solar sail model.
Eigen::Vector3d computeIdealRadiationPressureForce(
        const double radiationPressure,
        const Eigen::Vector3d& vectorFromSource,
        const Eigen::Vector3d& normalToSail,
        const double area,
        const double emissivity )
{
    double dotProduct = vectorFromSource.dot( normalToSail );

    if( dotProduct > 1.0 )
    {
        dotProduct = 1.0;
    }
    else if(dotProduct < -1.0)
    {
        dotProduct = -1.0;
    }

    double cone_Angle = std::acos(dotProduct);
    //double cone_Angle = 0.0;

    std::cout << "\n cone angle  \n " << cone_Angle;
    std::cout << "\n vector from source \n" << vectorFromSource;
    std::cout << "\n normal to sail \n" << normalToSail;
    std::cout << "\n pressure \n" << radiationPressure * area * std::cos( cone_Angle ) * ( ( 1.0 - emissivity ) * vectorFromSource +
             2.0 * emissivity  * std::cos( cone_Angle ) * normalToSail );

    return radiationPressure * area * std::cos( cone_Angle ) * ( ( 1.0 - emissivity ) * vectorFromSource +
            2.0 * emissivity  * std::cos( cone_Angle ) * normalToSail );

}

} // namespace electro_magnetism
} // namespace tudat
