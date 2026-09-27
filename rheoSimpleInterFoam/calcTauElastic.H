/*---------------------------------------------------------------------------*\
  calcTauElastic.H
  速度勾配テンソルから第2リブリン・エリクセンテンソル（A2）を代数近似し、
  弾性応力テンソル（tauElastic）を計算する。
\*---------------------------------------------------------------------------*/

{
    // 1. 速度勾配テンソル L = grad(U) [1/s]
    volTensorField L(fvc::grad(U));

    // 2. 変形速度テンソル D = 0.5 * (L + L^T) = symm(L) [1/s]
    volSymmTensorField D(symm(L));

    // 3. テンソル自乗項 D^2 = symm(D & D) [1/s^2]
    volSymmTensorField D2(symm(D & D));

    // 4. Rivlin-Ericksen A2 テンソルの代数近似 [1/s^2]
    //    A2 ≈ 2 * symm(L^T . D + D . L)
    //    ※ OpenFOAM の GeometricField の転置には T(L) を使用します
    volSymmTensorField A2_approx(2.0 * symm((T(L) & D) + (D & L)));

    // 5. 2次流体（CEF）代数モデルによる弾性応力の算出 [Pa = kg/(m s^2)]
    //    tauElastic = - 0.5 * Psi1 * A2 + 4 * (Psi1 + Psi2) * D^2
    tauElastic = - 0.5 * Psi1 * A2_approx + 4.0 * (Psi1 + Psi2) * D2;

    // 6. 二相流（VoF）マスキング処理
    //    空気相（alpha1 ≈ 0）で非物理的な弾性応力が発生するのを防ぐ
    tauElastic *= alpha1;

    // 7. 境界パッチの応力値を更新
    tauElastic.correctBoundaryConditions();
}

