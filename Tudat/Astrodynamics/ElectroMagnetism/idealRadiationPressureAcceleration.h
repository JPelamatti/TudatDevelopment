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
 *                                  function; added missing Doxygen comments.
 *
 *    References
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *
 *    Notes
 *
 */

#ifndef TUDAT_IDEAL_RADIATION_PRESSURE_ACCELERATION_H
#define TUDAT_IDEAL_RADIATION_PRESSURE_ACCELERATION_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/idealRadiationPressureForce.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute radiation pressure acceleration using a cannon-ball model.
/*!
 * Computes radiation pressure acceleration using a cannon-ball model, i.e. assuming force to be in
 * opposite direction of the vector to the source. This function is essentially a wrapper for the
 * function that computes the force.
 * opposite direction of the vector to the source.
 * \param radiationPressure Radiation pressure at target. N.B: the usual way of computing the
 *          radiation pressure at the target, in case the source is the Sun, is to take the
 *          radiation presssure at 1 AU and scale it based on the distance from the Sun     [N/m^2]
 * \param vectorToSource Vector pointing from target to source. N.B: this must be a unit
 *          vector! To compute the unit vector based on a given position vector, you can
 *          use the .normalize() or .normalized() member functions of an Eigen::Vector3d
 *          object.                                                                             [-]
 * \param area Area on which radiation pressure is assumed to act.                            [m^2]
 * \param radiationPressureCoefficient Coefficient to scale effective force. Equal to 1 +
 *          emmisivitty, assuming no diffuse reflection.                                        [-]
 * \param mass Mass of accelerated body.                                                       [kg]
 * \return Acceleration due to radiation pressure.                                          [m/s^2]
 * \sa computeCannonBallRadiationPressureForce().
 */
Eigen::Vector3d computeIdealRadiationPressureAcceleration(const double radiationPressure,
        const Eigen::Vector3d& vectorFromSource,
        const Eigen::Vector3d& normalToSail,
        const double area,
        const double radiationPressureCoefficient,
        const double mass );

//! Ideal Radiation pressure acceleration model class.

class IdealRadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel3d
{
private:

    //! Typedef for double-returning function.
    typedef boost::function< double( ) > DoubleReturningFunction;

    //! Typedef for Eigen::Vector3d-returning function.
    typedef boost::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

    //! Typedef for Eigen::Vector2d-returning function.
    typedef boost::function< Eigen::Vector2d( const double ) > Vector2dReturningFunction;

public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor taking function pointers for all variables.
    /*!
     * Constructor taking function pointers for all variables.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing
     *          radiation pressure acceleration.
     * \param radiationPressureFunction Function returning current radiation pressure.
     * \param radiationPressureCoefficientFunction Function returning current radiation pressure
     *          coefficient.
     * \param areaFunction Function returning current area assumed to undergo radiation pressure.
     * \param massFunction Function returning current mass of body undergoing acceleration.
     */
    IdealRadiationPressureAcceleration(
            Vector3dReturningFunction sourcePositionFunction,
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            Vector2dReturningFunction sailAnglesFunction,
            DoubleReturningFunction radiationPressureFunction,
            DoubleReturningFunction radiationPressureCoefficientFunction,
            DoubleReturningFunction areaFunction,
            DoubleReturningFunction massFunction )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          sailAnglesFunction_( sailAnglesFunction ),
          radiationPressureFunction_( radiationPressureFunction ),
          radiationPressureCoefficientFunction_( radiationPressureCoefficientFunction ),
          areaFunction_( areaFunction ),
          massFunction_( massFunction )
    {
    }

    //! Constructor taking functions pointers and constant values for parameters.
    /*!
     * Constructor taking function pointers for position vectors and radiation pressure and
     * constant values for other parameters.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing
     *              radiation pressure acceleration.
     * \param radiationPressureFunction Function returning current radiation pressure.
     * \param radiationPressureCoefficient Constant radiation pressure coefficient.
     * \param area Constant area assumed to undergo radiation pressure.
     * \param mass Constant mass of body undergoing acceleration.
     */
    IdealRadiationPressureAcceleration(
            Vector3dReturningFunction sourcePositionFunction,
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            Vector2dReturningFunction sailAnglesFunction,
            DoubleReturningFunction radiationPressureFunction,
            const double radiationPressureCoefficient,
            const double area,
            const double mass )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          sailAnglesFunction_( sailAnglesFunction ),
          radiationPressureFunction_( radiationPressureFunction ),
          radiationPressureCoefficientFunction_(
              boost::lambda::constant( radiationPressureCoefficient ) ),
          areaFunction_( boost::lambda::constant( area ) ),
          massFunction_( boost::lambda::constant( mass ) )
    {
    }

    //! Get radiation pressure acceleration.
    /*!
     * Returns the radiation pressure acceleration. No arguments are passed to this function.
     * Instead, all data required for computation is to be obtained from pointers to functions/
     * classes/structs, etc., which are to be set in a derived class and evaluated by the
     * updateMembers() function below. This function is essentially a wrapper for the free
     * function that computes the radiation pressure acceleration.
     * \return Radiation pressure acceleration.
     * \sa computeCannonBallRadiationPressureAcceleration().
     */
    Eigen::Vector3d getAcceleration( )
    {
        return computeIdealRadiationPressureAcceleration(
                    currentRadiationPressure_, currentVectorFromSource_, currentSailNormal_,
                    currentArea_, currentRadiationPressureCoefficient_, currentMass_ );
    }

    //! Update member variables used by the radiation pressure acceleration model.
    /*!
     * Updates member variables used by the acceleration model. This function evaluates all
     * dependent variables to the 'current' values of these parameters. Only these current values,
     * not the function-pointers are then used by the getAcceleration( ) function.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentVectorFromSource_ = ( acceleratedBodyPositionFunction_( ) -
                                         sourcePositionFunction_( )).normalized( );

           std::cout << " \n accelerated body position \n" << acceleratedBodyPositionFunction_( );
           std::cout << " \n source position \n" << sourcePositionFunction_( );


            currentRadiationPressure_ = radiationPressureFunction_( );
            currentSailAngles_ = sailAnglesFunction_( currentTime );

            setSailNormalVector( );

            currentRadiationPressureCoefficient_ = radiationPressureCoefficientFunction_( );
            currentArea_ = areaFunction_( );
            currentMass_ = massFunction_( );

            currentTime_ = currentTime;
        }
    }

    //! Function to retrieve the function pointer returning mass of accelerated body.
    /*!
     * Function to retrieve the function pointer returning mass of accelerated body.
     * \return Function pointer returning mass of accelerated body.
     */
    boost::function< double( ) > getMassFunction( )
    {
        return massFunction_;
    }

private:

    void setSailNormalVector( )
    {

        double gamma,beta;
        Eigen::MatrixXd R_12(3,3);
        Eigen::Vector3d intermediateNormal;

        std::cout << "\n Current vector from source" << currentVectorFromSource_;

        gamma = std::asin(currentVectorFromSource_(2)/currentVectorFromSource_.norm());
        beta = std::atan2(currentVectorFromSource_(1),currentVectorFromSource_(0));

        std::cout << "\n gamma" << gamma;
        std::cout << "\n beta" << beta;

        /*R_12(0,0) = cos(beta);
        R_12(0,1) = -sin(beta);
        R_12(0,2) = -cos(beta)*sin(gamma);
        R_12(1,0) = sin(beta)*cos(gamma);
        R_12(1,1) = cos(beta);
        R_12(1,2) = -sin(beta)*sin(gamma);
        R_12(2,0) = sin(gamma);
        R_12(2,1) = 0;
        R_12(2,2) = cos(gamma);*/

        R_12(0,0) = cos(beta)*cos(gamma);
        R_12(0,1) = -sin(beta);
        R_12(0,2) = -sin(gamma)*cos(beta);
        R_12(1,0) = cos(gamma)*sin(beta);
        R_12(1,1) = cos(beta);
        R_12(1,2) = -sin(beta)*sin(gamma);
        R_12(2,0) = sin(gamma);
        R_12(2,1) = 0;
        R_12(2,2) = cos(gamma);

        intermediateNormal(0) = cos(currentSailAngles_(0));
        intermediateNormal(1) = sin(currentSailAngles_(0))*sin(currentSailAngles_(1));
        intermediateNormal(2) = sin(currentSailAngles_(0))*cos(currentSailAngles_(1));

        currentSailNormal_ = R_12*intermediateNormal;
    }

    //! Function pointer returning position of source.
    /*!
     * Function pointer returning position of source (3D vector).
     */
    const Vector3dReturningFunction sourcePositionFunction_;

    //! Function pointer returning position of accelerated body.
    /*!
     * Function pointer returning position of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyPositionFunction_;

    //! Function pointer returning attitude angles of the sail.
    const Vector2dReturningFunction sailAnglesFunction_;

    //! Function pointer returning radiation pressure.
    /*!
     * Function pointer returning radiation pressure [N/m^{2}].
     */
    const DoubleReturningFunction radiationPressureFunction_;

    //! Function pointer returning radiation pressure coefficient.
    /*!
     * Function pointer returning radiation pressure coefficient [-].
     */
    const DoubleReturningFunction radiationPressureCoefficientFunction_;

    //! Function pointer returning area on which radiation pressure is acting.
    /*!
     * Function pointer returning area on which radiation pressure is acting [m^{2}].
     */
    const DoubleReturningFunction areaFunction_;

    //! Function pointer returning mass of accelerated body.
    /*!
     * Function pointer returning mass of accelerated body [kg].
     */
    const DoubleReturningFunction massFunction_;

    //! Current vector from accelerated body to source.
    /*!
     * Current vector from accelerated body to source (3D vector).
     */
    Eigen::Vector3d currentVectorFromSource_;

    //! Current sail attitude angles
    Eigen::Vector2d currentSailAngles_;

    // Current vector normal to the sail
    Eigen::Vector3d currentSailNormal_;

    //! Current radiation pressure.
    /*!
     * Current radiation pressure [N/m^{2}].
     */
    double currentRadiationPressure_;

    //! Current radiation pressure coefficient.
    /*!
     *  Current radiation pressure coefficient [-].
     */
    double currentRadiationPressureCoefficient_;

    //! Current area on which radiation pressure is acting.
    /*!
     * Current area on which radiation pressure is acting [m^{2}].
     */
    double currentArea_;

    //! Current mass of accelerated body.
    /*!
     * Current mass of accelerated body [kg].
     */
    double currentMass_;
};

//! Typedef for shared-pointer to IdealRadiationPressureAcceleration.
typedef boost::shared_ptr< IdealRadiationPressureAcceleration > IdealRadiationPressurePointer;

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_IDEAL_RADIATION_PRESSURE_ACCELERATION_H
