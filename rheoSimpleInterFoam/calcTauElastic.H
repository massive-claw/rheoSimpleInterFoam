/*---------------------------------------------------------------------------*\
  calcTauElastic.H
\*---------------------------------------------------------------------------*/

{
    // 1. Velocity gradient L and rate of strain tensor D
    volTensorField LRaw(fvc::grad(U));
    volSymmTensorField DRaw(symm(LRaw));

    // 2. Shear rate with upper bound clamping against wall/corner singularity
    volScalarField gammaDotRaw(Foam::sqrt(2.0) * mag(DRaw));
    volScalarField gammaScale
    (
        min
        (
            scalar(1.0),
            maxGammaDot / (gammaDotRaw + dimensionedScalar("eps", dimless/dimTime, 1e-6))
        )
    );

    volTensorField L(LRaw * gammaScale);
    volSymmTensorField D(DRaw * gammaScale);
    volScalarField gammaDot(gammaDotRaw * gammaScale);

    // 3. Carreau-type shear-attenuation for Psi1 and Psi2
    volScalarField attenuation
    (
        pow(scalar(1.0) + sqr(lambda * gammaDot), (scalar(1.0) - powerN.value()) / 2.0)
    );
    volScalarField Psi1Eff(Psi1 / attenuation);
    volScalarField Psi2Eff(Psi2 / attenuation);

    // 4. Tensor terms: D^2 and algebraical A2 approximation
    volSymmTensorField D2(symm(D & D));
    volSymmTensorField A2_approx(4.0 * symm(D & L.T()));

    // 5. Algebraic elastic stress tensor
    volSymmTensorField tauNew
    (
        - 0.5 * Psi1Eff * A2_approx + 4.0 * (Psi1Eff + Psi2Eff) * D2
    );

    // 6. Stress limiter against singularity (based on phase-1 viscosity)
    volScalarField mu1(rho1 * mixture.nuModel1().nu());
    volScalarField tauLim
    (
        maxStressRatio * (2.0 * mu1 * gammaDot + dimensionedScalar("eps", dimPressure, 1e-4))
    );
    volScalarField magTau(mag(tauNew));
    volScalarField limiter
    (
        min(scalar(1.0), tauLim / (magTau + dimensionedScalar("eps", dimPressure, 1e-6)))
    );
    tauNew *= limiter;

    // 7. VoF phase masking, under-relaxation, and boundary condition update
    tauNew *= alpha1;
    tauElastic = (scalar(1.0) - relaxTau) * tauElastic + relaxTau * tauNew;
    tauElastic.correctBoundaryConditions();
}
