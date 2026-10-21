    MRF.correctBoundaryVelocity(U);

    // 1. 最新の U から弾性応力テンソル tauElastic を計算
    #include "calcTauElastic.H"

    volScalarField muNumEff("muNumEff", alpha1 * muNum);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      - fvm::laplacian(muNumEff, U)
     ==
        fvOptions(rho, U)
      + fvc::div(tauElastic)
      - fvc::laplacian(muNumEff, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
