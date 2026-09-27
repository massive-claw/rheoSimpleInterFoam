/*---------------------------------------------------------------------------*\
  calcTauElastic.H
\*---------------------------------------------------------------------------*/

{
    // 1. Velocity gradient L and rate of strain tensor D
    volTensorField L(fvc::grad(U));
    volSymmTensorField D(symm(L));

    // 2. Shear rate: gammaDot = sqrt(2 * D:D)
    volScalarField gammaDot(Foam::sqrt(2.0) * mag(D));

    // 3. Carreau-type shear-attenuation for Psi1
    volScalarField Psi1Eff
    (
        Psi1 / pow(scalar(1.0) + sqr(lambda * gammaDot), (scalar(1.0) - powerN.value()) / 2.0)
    );

    // 4. Tensor terms: D^2 and algebraical A2 approximation
    volSymmTensorField D2(symm(D & D));
    volSymmTensorField A2_approx(4.0 * symm(D & L));

    // 5. Algebraic elastic stress tensor
    tauElastic = - 0.5 * Psi1Eff * A2_approx + 4.0 * (Psi1Eff + Psi2) * D2;

    // 6. Stress limiter against singularity
    volScalarField tauLim
    (
        maxStressRatio * (2.0 * mixture.mu() * gammaDot + dimensionedScalar("eps", dimPressure, 1e-4))
    );
    volScalarField magTau(mag(tauElastic));
    volScalarField limiter(min(scalar(1.0), tauLim / (magTau + dimensionedScalar("eps", dimPressure, 1e-6))));
    tauElastic *= limiter;

    // 7. VoF phase masking and boundary condition update
    tauElastic *= alpha1;
    tauElastic.correctBoundaryConditions();
}
