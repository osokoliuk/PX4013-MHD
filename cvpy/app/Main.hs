#!/usr/bin/env cabal
{- cabal:
build-depends:
  base,
  parallel
ghc-options: -main-is Fourier -O2 -threaded
-}
{-# LANGUAGE ExplicitNamespaces #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

module Main where

{-
Module      : Main.hs
Description : Discrete Fourier methods module
Copyright   : (c) Oleksii Sokoliuk, 2026
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

A module that defines several useful functions computed in parallel,
namely Discrete Fourier Transform, Power Spectrum and higher order
polyspectra
-}

-- Module imports

import Control.Parallel.Strategies
import Data.Complex (Complex (..), conjugate, imagPart, realPart)
import Data.List (transpose, zipWith4)
import qualified Data.Map as MP
import Data.Maybe (fromMaybe)

-- Some type definitions
type WindowFunction = Complex Double -> Complex Double
type SignalFunction = Complex Double -> Complex Double
type Index_k = Complex Double
type Index_l = Complex Double
type Index_p = Complex Double
type N = Complex Double
type T = Complex Double
type M = Complex Double
type T0 = Complex Double
type Ak = Complex Double
type Al = Complex Double
type Ap = Complex Double
type Akl = Complex Double
type Aklp = Complex Double
type OmegaT = Complex Double
type C2 = Complex Double
type C3 = Complex Double
type C4 = Complex Double
type S2 = Complex Double
type S3 = Complex Double
type S4 = Complex Double

-- | Define a shortcut for a number with Im(x) = 0
real ::
    Double ->
    Complex Double
real x = x :+ 0

-- | Average value from an array of numbers
average ::
    [Complex Double] ->
    Complex Double
average x = sum x / (real . fromIntegral . length $ x)

twoAverage ::
    [Complex Double] ->
    [Complex Double] ->
    Complex Double
twoAverage x y = average (zipWith (\x y -> x * y) x y)

threeAverage ::
    [Complex Double] ->
    [Complex Double] ->
    [Complex Double] ->
    Complex Double
threeAverage x y z = average (zipWith3 (\x y z -> x * y * z) x y z)

fourAverage ::
    [Complex Double] ->
    [Complex Double] ->
    [Complex Double] ->
    [Complex Double] ->
    Complex Double
fourAverage x y z w = average (zipWith4 (\x y z w -> x * y * z * w) x y z w)

-- | Gaussian function with an amplitude 1 and standard deviation t * omegaT
gaussian ::
    N ->
    T ->
    OmegaT ->
    Complex Double ->
    Complex Double
gaussian n t omegaT x = exp (-((x - n / 2) / (2 * t * omegaT)) ** 2)

-- | Approximate confined Gaussian window function
windowFunction ::
    N ->
    T ->
    OmegaT ->
    Complex Double ->
    Complex Double
windowFunction n t omegaT x =
    gaussian t n omegaT x
        - ( gaussian t n omegaT (real (-1 / 2))
                * (gaussian t n omegaT (x + t) + gaussian t n omegaT (x - t))
          )
            / (gaussian t n omegaT (real (-1 / 2) + t) + gaussian t n omegaT (real (-1 / 2) - t))

{- | Calculate Fourier coefficients for a given function and a given window function
where the size of the window is T, containing N samples and starting from t = t0
with a wavenumber k
-}
fourierCoefficient ::
    T0 ->
    N ->
    T ->
    Index_k ->
    WindowFunction ->
    SignalFunction ->
    Ak
fourierCoefficient t0 n t k g f =
    let i = (0 :+ 1)
     in (t / n)
            * sum
                [ g (j * t / n)
                    * f (j * t / n - t0)
                    * exp (2 * pi * i * k * j / n)
                    * exp (-2 * pi * i * k * t0 / t)
                | j <- real <$> [0, 1 .. realPart n - 1]
                ]

{- | Calculate the values for the unbiased cumulants C2, C3, C4 and the
corresponding polyspectra using the formulas
given in 10.1016/j.dsp.2026.105893
-}
calculatePolyspectra ::
    M ->
    N ->
    T ->
    OmegaT ->
    [[Ak]] ->
    [Index_k] ->
    [Index_l] ->
    [Index_p] ->
    ([S2], [[S3]], [[S4]])
calculatePolyspectra m n t omegaT ak k l p =
    let
        -- Create map for a fast lookup of Fourier coefficients
        akMap :: MP.Map Double [Ak]
        akMap = MP.fromList $ zip ([0, 1 .. realPart n - 1]) (transpose ak)

        -- Shortcut for a window function
        g = windowFunction n t omegaT

        -- Second order cumulant
        c2 :: M -> [Ak] -> [Ak] -> C2
        c2 m x y =
            m / (m - 1) * (twoAverage x y - average x * average y)

        -- Third order cumulant
        c3 :: M -> [Ak] -> [Ak] -> [Ak] -> C3
        c3 m x y z =
            (m ** 2 / ((m - 1) * (m - 2)))
                * ( threeAverage x y z
                        - twoAverage x y * average z
                        - twoAverage x z * average y
                        - twoAverage y z * average x
                        + 2 * average x * average y * average z
                  )

        -- Fourth order cumulant
        c4 :: M -> [Ak] -> [Ak] -> [Ak] -> [Ak] -> C4
        c4 m x y z w =
            ( m ** 2
                / ((m - 1) * (m - 2) * (m - 3))
                * ((m + 1) * fourAverage x y z w)
                - (m + 1) * (threeAverage x y z * average w + threeAverage x y w * average z + threeAverage x z w * average y + threeAverage y z w * average x)
                - (m - 1) * (twoAverage x y * twoAverage z w + twoAverage x z * twoAverage y w + twoAverage x w * twoAverage y z)
                + 2 * m * (twoAverage x y * average z * average w + twoAverage x z * average y * average w + twoAverage x w * average y * average z + twoAverage y z * average x * average w + twoAverage y w * average x * average z + twoAverage z w * average x * average y)
                - 6 * m * average x * average y * average z * average w
            )

        -- Second order polyspectrum (i.e. power spectrum)
        s2 :: M -> N -> T -> [Ak] -> S2
        s2 m n t ak =
            (n * (c2 m ak (conjugate <$> ak)))
                / (t * sum [g (j * t / n) * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Third order polyspectrum (i.e. bispectrum)
        s3 :: M -> N -> T -> MP.Map Double [Ak] -> Index_k -> Index_l -> S3
        s3 m n t akMap k l =
            let ak = fromMaybe [real 0] $ MP.lookup (realPart k) akMap
                al = fromMaybe [real 0] $ MP.lookup (realPart l) akMap
                akl = fromMaybe [real 0] $ MP.lookup (realPart (k + l)) akMap
             in (n * (c3 m ak al (conjugate <$> akl)))
                    / (t * sum [g (j * t / n) ** 2 * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Fourth order polyspectrum (i.e. bispectrum)
        s4 :: M -> N -> T -> MP.Map Double [Ak] -> Index_k -> Index_l -> Index_p -> S4
        s4 m n t akMap k l p =
            let ak = fromMaybe [real 0] $ MP.lookup (realPart k) akMap
                al = fromMaybe [real 0] $ MP.lookup (realPart l) akMap
                ap = fromMaybe [real 0] $ MP.lookup (realPart p) akMap
                aklp = fromMaybe [real 0] $ MP.lookup (realPart (k + l + p)) akMap
             in (n * (c4 m ak al ap (conjugate <$> aklp)))
                    / (t * sum [g (j * t / n) ** 3 * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Slice of a trispectrum along the diagonal S4(w1,-w1,w2), retains most of the information
        -- while being two dimensional and much cheaper to calculate
        compactS4 :: M -> N -> T -> MP.Map Double [Ak] -> Index_k -> Index_l -> S4
        compactS4 m n t akMap k l =
            let ak = fromMaybe [real 0] $ MP.lookup (realPart k) akMap
                al = fromMaybe [real 0] $ MP.lookup (realPart l) akMap
             in (n * (c4 m ak (conjugate <$> ak) al (conjugate <$> al)))
                    / (t * sum [g (j * t / n) ** 2 * (conjugate (g (j * t / n))) ** 2 | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Evaluate all of the polyspectra in parallel
        evaluatedS2 = (parMap rdeepseq (\x -> s2 m n t x) ak)
        evaluatedS3 = parMap rdeepseq (\y -> parMap rdeepseq (\x -> s3 m n t akMap x y) k) l
        evaluatedCompactS4 = parMap rdeepseq (\y -> parMap rdeepseq (\x -> compactS4 m n t akMap x y) k) l
     in
        (evaluatedS2, evaluatedS3, evaluatedCompactS4)

main :: IO ()
main = do
    let m = real 100
        n = real 100
        t = real 0.1
        omegaT = 0.14
        k = real <$> [0, 1 .. realPart n - 1]
        l = k
        p = k
        ak = parMap rdeepseq (\t0 -> (\k -> fourierCoefficient t0 n t k (windowFunction n t omegaT) (cos)) <$> (real <$> [0, 1 .. realPart n - 1])) (real <$> [0.0, realPart t .. realPart $ m * t])
        (s2, s3, s4) = (calculatePolyspectra m n t omegaT (transpose ak) k l p)
    print $ (fmap . fmap) realPart s4
    -- print $ realPart <$> s2
    print ""
