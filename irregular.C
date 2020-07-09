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

\*---------------------------------------------------------------------------*/

#include "irregular.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(irregular, 0);
    addToRunTimeSelectionTable(waveModel, irregular, objectRegistry);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::irregular::k() const
{
    return 2*Foam::constant::mathematical::pi/length_;
}

Foam::scalar Foam::waveModels::irregular::celerity() const
{
    return sqrt(g()/k()*tanh(k()*depth()));
}

Foam::tmp<Foam::scalarField> Foam::waveModels::irregular::angle
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    return phase_ + k()*(x - (u + celerity())*t);
}

bool Foam::waveModels::irregular::deep() const
{
    return k()*depth() > log(GREAT);
}

Foam::tmp<Foam::vector2DField> Foam::waveModels::irregular::vi
(
    const label i,
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    const scalarField x(xz.component(0));
    const scalarField z(xz.component(1));

    const scalarField phi(angle(t, u, x));
    const scalarField kz(- k()*z);

    if (deep())
    {
        return i*exp(- mag(kz))*zip(cos(i*phi), sin(i*phi));
    }
    else
    {
        const scalar kd = k()*depth();
        const scalarField kdz(max(scalar(0), kd - mag(kz)));
        return i*zip(cosh(i*kdz)*cos(i*phi), sinh(i*kdz)*sin(i*phi))/sinh(kd);
    }
}

Foam::scalar Foam::waveModels::irregular::kindex(int index) const
{
    return 2*Foam::constant::mathematical::pi/L[index];
}


Foam::scalar Foam::waveModels::irregular::celerity_i(int index) const
{
    return sqrt(g()/kindex(index)*tanh(kindex(index)*depth()));
}


Foam::tmp<Foam::scalarField> Foam::waveModels::irregular::angle_i
(
    const scalar t,
    const scalar u,
    const scalarField& x,
    int index
) const
{
    //scalar omega = celerity_i(index) + kindex(index)*u;
    //return kindex(index)*x - omega*t+Phase[index];
    return Phase[index] + kindex(index)*(x - (u + celerity_i(index))*t);
}


bool Foam::waveModels::irregular::deep_i(int index) const
{
    if(fabs(depth()-log(2*GREAT)/k())<1e-5)
       return true;
    else 
       return kindex(index)*depth() > log(GREAT);
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::irregular::vi_i
(
    const label i,
    const scalar t,
    const scalar u,
    const vector2DField& xz,
    int index
) const
{
    const scalarField x(xz.component(0));
    const scalarField z(xz.component(1));

    const scalarField phi(angle_i(t, u, x,index));
    
    //scalarField z1
    //scalarField z1(z+elevation(t,u,x)); //(z = z1-elevation)
    scalarField z2(z-elevation_i(t,u,x,index));
    const scalarField kz(- kindex(index)*mag(z2));
    //const scalarField kz(- kindex(index)*z);

    if (deep_i(index))
    {
        return i*exp(- mag(kz))*zip(cos(i*phi), sin(i*phi));
    }
    else
    {
        const scalar kd = kindex(index)*depth();
        const scalarField kdz(max(scalar(0), kd - mag(kz)));
        return i*zip(cosh(i*kdz)*cos(i*phi), sinh(i*kdz)*sin(i*phi))/sinh(kd);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::irregular::irregular(const irregular& wave)
:
    waveModel(wave),
    length_(wave.length_),
    phase_(wave.phase_),
    depth_(wave.depth_)
{

}


Foam::waveModels::irregular::irregular
(
    const objectRegistry& db,
    const dictionary& dict
)
:
    waveModel(db, dict),
    length_(readScalar(dict.lookup("length"))),
    phase_(readScalar(dict.lookup("phase"))),
    depth_(dict.lookupOrDefault<scalar>("depth", log(2*GREAT)/k()))
{
     Info<<"read the wave properties"<<endl;
     dictionary dict2(IFstream("constant/irregularWaveProperties")());
    
     List<scalar> lengthW_(dict2.subDict("length").lookup("value"));
     L.swap(lengthW_);
     
     List<scalar> aW_(dict2.subDict("amplitude").lookup("value"));
     A.swap(aW_);
     
     List<scalar> p_(dict2.subDict("phase").lookup("value"));
     Phase.swap(p_);
     
     // test the length wether the same
     if( L.size() == A.size() )
     {
         if(L.size() ==Phase.size())
            Info<<"the wave component: "<<L.size()<<endl;
         else
         {
            Info<<"the Phase component is not compatible, it is: "
            <<Phase.size()<<endl;
         }
     }
     else
     {
        Info<<"the amplitude component is not compatible, it is: "
            <<A.size()<<endl;
     }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::irregular::~irregular()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveModels::irregular::elevation
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    
    scalarField eleAll(x.size(), 0);
    int waveCount = L.size();
    for(int i=0; i<waveCount; i++)
    {
        eleAll += elevation_i(t,u,x,i);
    }
    //return eleAll;
    return eleAll;
}

Foam::tmp<Foam::scalarField> Foam::waveModels::irregular::elevation_i
(
    const scalar t,
    const scalar u,
    const scalarField& x,
    int index
) const
{
    return A[index]*cos(angle_i(t, u, x,index));
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::irregular::velocity
(
    const scalar t,
    const scalar u,
    const vector2DField& xz
) const
{
    //const scalar ka = k()*amplitude(t);

    //return celerity()*ka*vi(1, t, u, xz);
    vector2DField uv(xz*0);
    scalar wa = 0;
    for(int i=0; i<L.size(); i++)
    {
        //wa = kindex(i)*(celerity_i(i)+u)*A[i];
        wa = kindex(i)*A[i];
        uv += wa*celerity_i(i)*vi_i(1, t, u, xz,i);
    }
    return uv;
}


void Foam::waveModels::irregular::write(Ostream& os) const
{
    waveModel::write(os);

    os.writeKeyword("length") << length_ << token::END_STATEMENT << nl;
    os.writeKeyword("phase") << phase_ << token::END_STATEMENT << nl;
    if (!deep())
    {
        os.writeKeyword("depth") << depth_ << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
