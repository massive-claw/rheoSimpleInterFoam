    MRF.correctBoundaryVelocity(U);

    // 1. 最新の U から弾性応力テンソル tauElastic を計算
    #include "calcTauElastic.H"

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U) // 溶媒粘性項（陰的解法）はそのまま活かす
     ==
        fvOptions(rho, U)
      + fvc::div(tauElastic)              // 2. 弾性付加応力の発散項を右辺に追加
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
