/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::waveModels::irregular

Description
    First-order wave model.

    Reference:
    \verbatim
        "On the Theory of Oscillatory Waves"
        G G Stokes
        Transactions of the Cambridge Philosophical Society (1847), Volume 8,
        Pages 441-455, Equations 18 and 19 (leading terms only)
    \endverbatim

SourceFiles
    irregular.C

\*---------------------------------------------------------------------------*/

#ifndef irregular_H
#define irregular_H

#include "waveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                           Class irregular Declaration
\*---------------------------------------------------------------------------*/

class irregular
:
    public waveModel
{
    //- Private data

        //- Peak-to-peak length [m]
        const scalar length_;

        //- Phase offset [rad]
        const scalar phase_;

        //- Depth [m]
        const scalar depth_;
        
        List<scalar> L;
        List<scalar> A;
        List<scalar> Phase;


protected:

    // Protected Member Functions

        //- The angular wavenumber [rad/m]
        scalar k() const;
        
        //- The angular wavenumber [rad/m]
        scalar kindex(int index) const;

        //- The wave celerity [m/s]
        scalar celerity() const;
        
        //- The wave celerity [m/s]
        scalar celerity_i(int index) const;

        //- Angle of the oscillation [rad]
        tmp<scalarField> angle
        (
            const scalar t,
            const scalar u,
            const scalarField& x
        ) const;
        
        //- Angle of the oscillation [rad]
        tmp<scalarField> angle_i
        (
            const scalar t,
            const scalar u,
            const scalarField& x,
            int index
        ) const;

        //- Return whether shallow and intermdiate effects are to be omitted
        bool deep() const;
        
        //- Return whether shallow and intermdiate effects are to be omitted
        bool deep_i(int index) const;

        //- Return the non-dimensionalised i-th harmonic of the velocity
        tmp<vector2DField> vi
        (
            const label i,
            const scalar t,
            const scalar u,
            const vector2DField& xz
        ) const;
        
        //- Return the non-dimensionalised i-th harmonic of the velocity
        tmp<vector2DField> vi_i
        (
            const label i,
            const scalar t,
            const scalar u,
            const vector2DField& xz,
            int index
        ) const;


public:

    //- Runtime type information
    TypeName("irregular");


    // Constructors

        //- Construct a copy
        irregular(const irregular& wave);

        //- Construct from a database and a dictionary
        irregular(const objectRegistry& db, const dictionary& dict);

        //- Construct a clone
        virtual autoPtr<waveModel> clone() const
        {
            return autoPtr<waveModel>(new irregular(*this));
        }


    //- Destructor
    virtual ~irregular();


    // Member Functions

        // Access

            //- Get the length
            scalar length() const
            {
                return length_;
            }

            //- Get the phase
            scalar phase() const
            {
                return phase_;
            }

            //- Get the depth
            scalar depth() const
            {
                return depth_;
            }

        //- Get the wave elevation at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity.
        virtual tmp<scalarField> elevation
        (
            const scalar t,
            const scalar u,
            const scalarField& x
        ) const;
        
        //- Get the wave elevation at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity.
        virtual tmp<scalarField> elevation_i
        (
            const scalar t,
            const scalar u,
            const scalarField& x,
            int index
        ) const;

        //- Get the wave velocity at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity, and z with
        //  negative gravity.
        virtual tmp<vector2DField> velocity
        (
            const scalar t,
            const scalar u,
            const vector2DField& xz
        ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
